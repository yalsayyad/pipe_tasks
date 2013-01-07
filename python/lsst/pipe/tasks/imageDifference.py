#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2012 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#
import numpy as num

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.daf.base as dafBase
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
from lsst.meas.algorithms import SourceDetectionTask, SourceMeasurementTask, SourceDeblendTask, starSelectorRegistry
from lsst.ip.diffim import ImagePsfMatchTask
import lsst.afw.display.ds9 as ds9                                                                                                                                                                          

class ImageDifferenceConfig(pexConfig.Config):
    """Config for ImageDifferenceTask
    """
    doSelectSources = pexConfig.Field(dtype=bool, default=True, doc = "Select stars to use for kernel fitting")
    doSubtract = pexConfig.Field(dtype=bool, default=True, doc = "Compute subtracted exposure?")
    doDetection = pexConfig.Field(dtype=bool, default=True, doc = "Detect sources?")
    doDeblend = pexConfig.Field(dtype=bool, default=False,
        doc = "Deblend sources? Off by default because it may not be useful")
    doMeasurement = pexConfig.Field(dtype=bool, default=True, doc = "Measure sources?")
    doWriteSubtractedExp = pexConfig.Field(dtype=bool, default=True, doc = "Write difference exposure?")
    doWriteMatchedExp = pexConfig.Field(dtype=bool, default=False,
        doc = "Write warped and PSF-matched template coadd exposure?")
    doWriteSources = pexConfig.Field(dtype=bool, default=True, doc = "Write sources?")
    doWriteHeavyFootprintsInSources = pexConfig.Field(dtype=bool, default=False,
        doc = "Include HeavyFootprint data in source table?")
                                                      
    coaddName = pexConfig.Field(
        doc = "coadd name: typically one of deep or goodSeeing",
        dtype = str,
        default = "deep",
    )
    templateFwhm = pexConfig.Field(
        doc = "FWHM (arcsec) of Psf in template",
        dtype = float,
        default = 1.0
    )
    swapImageToConvolve = pexConfig.Field(
        doc = "Swap the order of which image gets convolved (default = template)",
        dtype = bool,
        default = False
    )

    starSelector = starSelectorRegistry.makeField("Star selection algorithm", default="secondMoment")
    selectDetection = pexConfig.ConfigurableField(
        target = SourceDetectionTask,
        doc = "Initial detections used to feed stars to kernel fitting",
    )
    selectMeasurement = pexConfig.ConfigurableField(
        target = SourceMeasurementTask,
        doc = "Initial measurements used to feed stars to kernel fitting",
    )

    subtract = pexConfig.ConfigurableField(
        target = ImagePsfMatchTask,
        doc = "Warp and PSF match template to exposure, then subtract",
    )
    detection = pexConfig.ConfigurableField(
        target = SourceDetectionTask,
        doc = "Low-threshold detection for final measurement",
    )
    deblend = pexConfig.ConfigurableField(
        target = SourceDeblendTask,
        doc = "Split blended sources into their components",
    )
    measurement = pexConfig.ConfigurableField(
        target = SourceMeasurementTask,
        doc = "Final source measurement on low-threshold detections",
    )
    
    def setDefaults(self):
        # TODO HERE: 
        # Make the minimal set of algorithsm to select stars, for speed optimization
        self.selectDetection.reEstimateBackground = False
        self.selectMeasurement.prefix = "select."
        self.selectMeasurement.doApplyApCorr = False
        self.starSelector["secondMoment"].clumpNSigma  = 2.0
        self.starSelector['secondMoment'].badFlags = [self.selectMeasurement.prefix+x for x in self.starSelector['secondMoment'].badFlags]

        self.detection.thresholdPolarity = "both"
        self.detection.reEstimateBackground = False

    def validate(self):
        pexConfig.Config.validate(self)
        if self.doMeasurement and not self.doDetection:
            raise ValueError("Cannot run source measurement without source detection.")
        if self.doDeblend and not self.doDetection:
            raise ValueError("Cannot run source deblending without source detection.")
        if self.doWriteHeavyFootprintsInSources and not self.doWriteSources:
            raise ValueError("Cannot write HeavyFootprints (doWriteHeavyFootprintsInSources) without doWriteSources")


class ImageDifferenceTask(pipeBase.CmdLineTask):
    """Subtract an image from a template coadd and measure the result
    """
    ConfigClass = ImageDifferenceConfig
    _DefaultName = "imageDifference"

    def __init__(self, **kwargs):
        pipeBase.CmdLineTask.__init__(self, **kwargs)
        self.makeSubtask("subtract")

        if self.config.doSelectSources:
            self.selectSchema = afwTable.SourceTable.makeMinimalSchema()
            self.selectAlgMetadata = dafBase.PropertyList()
            self.starSelector = self.config.starSelector.apply()
            self.makeSubtask("selectDetection", schema=self.selectSchema)
            self.makeSubtask("selectMeasurement", schema=self.selectSchema, algMetadata=self.selectAlgMetadata)

        self.schema = afwTable.SourceTable.makeMinimalSchema()
        self.algMetadata = dafBase.PropertyList()
        if self.config.doDetection:
            self.makeSubtask("detection", schema=self.schema)
        if self.config.doDeblend:
            self.makeSubtask("deblend", schema=self.schema)
        if self.config.doMeasurement:
            self.makeSubtask("measurement", schema=self.schema, algMetadata=self.algMetadata)

    @pipeBase.timeMethod
    def run(self, sensorRef):
        """Subtract an image from a template coadd and measure the result
    
        Steps include:
        - warp template coadd to match WCS of image
        - PSF match image to warped template
        - subtract image from PSF-matched, warped template
        - persist difference image
        - detect sources
        - deblend sources (disabled by default)
        - measure sources
        
        @param sensorRef: sensor-level butler data reference, used for the following data products:
        Input only:
        - calexp
        - psf
        - ccdExposureId
        - ccdExposureId_bits
        - apCorr
        - self.config.coaddName + "Coadd_skyMap"
        - self.config.coaddName + "Coadd"
        Input or output, depending on config:
        - self.config.coaddName + "Diff_subtractedExp"
        Output, depending on config:
        - self.config.coaddName + "Diff_matchedExp"
        - self.config.coaddName + "Diff_src"
            
        @return pipe_base Struct containing these fields:
        - subtractedExposure: exposure after subtracting template;
            the unpersisted version if subtraction not run but detection run
            None if neither subtraction nor detection run (i.e. nothing useful done)
        - subtractRes: results of subtraction task; None if subtraction not run
        - sources: detected and possibly measured and deblended sources; None if detection not run
        """
        self.log.log(self.log.INFO, "Processing %s" % (sensorRef.dataId))

        # initialize outputs
        subtractedExposure = None
        subtractRes = None
        sources = None

        # We make one IdFactory that will be used by both icSrc and src datasets;
        # I don't know if this is the way we ultimately want to do things, but at least
        # this ensures the source IDs are fully unique.
        expBits = sensorRef.get("ccdExposureId_bits")
        expId = long(sensorRef.get("ccdExposureId"))
        idFactory = afwTable.IdFactory.makeSource(expId, 64 - expBits)
        
        exposure = sensorRef.get("calexp")
        psf = sensorRef.get("psf")
        exposure.setPsf(psf)

        subtractedExposureName = self.config.coaddName + "Diff_subtractedExp"
        
        if self.config.doSubtract:
            templateExposure = self.getTemplate(exposure, sensorRef)
            if self.config.templateFwhm:
                wcs = templateExposure.getWcs()
                fwhmPixels = self.config.templateFwhm / wcs.pixelScale().asArcseconds()
            else:
                fwhmPixels = None

            if self.config.doSelectSources:
                # Detect on exposure since that ensures maximal astrometric coverage
                table = afwTable.SourceTable.make(self.selectSchema, idFactory)
                table.setMetadata(self.selectAlgMetadata) 
                detRet = self.selectDetection.makeSourceCatalog(table, exposure)
                selectSources = detRet.sources
                self.selectMeasurement.measure(exposure, selectSources)
                kernelCandidateList = self.starSelector.selectStars(exposure, selectSources)
                kernelSources = [x.getSource() for x in kernelCandidateList]
            else:
                selectSources = None
                kernelSources = None

            # warp template exposure to match exposure,
            # PSF match template exposure to exposure,
            # then return the difference
            subtractRes = self.subtract.subtractExposures(
                exposureToConvolve = templateExposure,
                exposureToNotConvolve = exposure,
                psfFwhmPixTc = fwhmPixels,
                candidateList = kernelSources,
                swapImageToConvolve = self.config.swapImageToConvolve
            )
            subtractedExposure = subtractRes.subtractedExposure

            if self.config.doWriteSubtractedExp:
                sensorRef.put(subtractedExposure, subtractedExposureName)
            if self.config.doWriteMatchedExp:
                sensorRef.put(subtractRes.matchedExposure, self.config.coaddName + "Diff_matchedExp")

        if self.config.doDetection:
            if subtractedExposure is None:
                subtractedExposure = sensorRef.get(subtractedExposureName)
            
            # Get from the calexp!
            if not subtractedExposure.hasPsf():
                psf = sensorRef.get("psf")
                subtractedExposure.setPsf(psf)

            # Erase existing detection mask planes
            mask  = subtractedExposure.getMaskedImage().getMask()
            mask &= ~(mask.getPlaneBitMask("DETECTED") | mask.getPlaneBitMask("DETECTED_NEGATIVE"))

            table = afwTable.SourceTable.make(self.schema, idFactory)
            table.setMetadata(self.algMetadata)
            sources = self.detection.makeSourceCatalog(table, subtractedExposure).sources

            if self.config.doDeblend:
               self.deblend.run(subtractedExposure, sources, psf)
    
            if self.config.doMeasurement:
                apCorr = sensorRef.get("apCorr")
                self.measurement.run(subtractedExposure, sources, apCorr)
    
            if sources is not None and self.config.doWriteSources:
                if self.config.doWriteHeavyFootprintsInSources:
                    sources.setWriteHeavyFootprints(True)
                sensorRef.put(sources, self.config.coaddName + "Diff_src")
 
        self.runDebug(exposure, subtractRes, selectSources, kernelSources, sources)
        return pipeBase.Struct(
            subtractedExposure = subtractedExposure,
            subtractRes = subtractRes,
            sources = sources,
        )

    def runDebug(self, exposure, subtractRes, selectSources, kernelSources, sources):
        import lsstDebug
        display = lsstDebug.Info(__name__).display 
        showSubtracted = lsstDebug.Info(__name__).showSubtracted
        showPixelResiduals = lsstDebug.Info(__name__).showPixelResiduals
        showDiaSources = lsstDebug.Info(__name__).showDiaSources
        maskTransparency = lsstDebug.Info(__name__).maskTransparency   
        if not maskTransparency:
            maskTransparency = 0
        ds9.setMaskTransparency(maskTransparency)

        if display and showSubtracted:
            ds9.mtv(subtractRes.subtractedExposure, frame=lsstDebug.frame, title="Subtracted image")
            mi = subtractRes.subtractedExposure.getMaskedImage()
            x0, y0 = mi.getX0(), mi.getY0()
            with ds9.Buffering():
                for s in sources:
                    x, y = s.getX() - x0, s.getY() - y0
                    ctype = "red" if s.get("flags.negative") else "yellow"
                    if (s.get("flags.pixel.interpolated.center") or s.get("flags.pixel.saturated.center") or
                        s.get("flags.pixel.cr.center")):
                        ptype = "x"
                    elif (s.get("flags.pixel.interpolated.any") or s.get("flags.pixel.saturated.any") or
                          s.get("flags.pixel.cr.any")):
                        ptype = "+"
                    else:
                        ptype = "o"
                    ds9.dot(ptype, x, y, size=4, frame=lsstDebug.frame, ctype=ctype)
            lsstDebug.frame += 1

        if display and showPixelResiduals and selectSources:
            import lsst.ip.diffim.utils as diUtils
            nonKernelSources = []
            for source in selectSources:
                if not source in kernelSources:
                    nonKernelSources.append(source)

            diUtils.plotPixelResiduals(exposure,
                                       subtractRes.warpedExposure,
                                       subtractRes.subtractedExposure,
                                       subtractRes.kernelCellSet,
                                       subtractRes.psfMatchingKernel,
                                       subtractRes.backgroundModel,
                                       nonKernelSources,
                                       self.subtract.config.kernel.active.detectionConfig,
                                       origVariance = False)
            diUtils.plotPixelResiduals(exposure,
                                       subtractRes.warpedExposure,
                                       subtractRes.subtractedExposure,
                                       subtractRes.kernelCellSet,
                                       subtractRes.psfMatchingKernel,
                                       subtractRes.backgroundModel,
                                       nonKernelSources,
                                       self.subtract.config.kernel.active.detectionConfig, 
                                       origVariance = True)
        if display and showDiaSources:
            import lsst.ip.diffim.diffimTools as diffimTools
            flagChecker   = diffimTools.SourceFlagChecker(sources)
            dipoleChecker = diffimTools.DipoleChecker(sources)
            isFlagged     = [flagChecker(x) for x in sources]
            isDipole      = [dipoleChecker(x) for x in sources]
            diUtils.showDiaSources(sources, subtractRes.subtractedExposure, isFlagged, isDipole, frame=lsstDebug.frame)
            lsstDebug.frame += 1
        
    def getTemplate(self, exposure, sensorRef):
        """Return a template coadd exposure that overlaps the exposure
        
        @param[in] exposure: exposure
        @param[in] sensorRef: a Butler data reference that can be used to obtain coadd data

        @return coaddExposure: a template coadd exposure assembled out of patches
        
        @note: the coadd consists of whole patches stitched together, so it may be larger than necessary
        """
        skyMap = sensorRef.get(datasetType=self.config.coaddName + "Coadd_skyMap")
        expWcs = exposure.getWcs()
        expBoxD = afwGeom.Box2D(exposure.getBBox(afwImage.PARENT))
        ctrSkyPos = expWcs.pixelToSky(expBoxD.getCenter())
        tractInfo = skyMap.findTract(ctrSkyPos)
        self.log.info("Using skyMap tract %s" % (tractInfo.getId(),))
        skyCorners = [expWcs.pixelToSky(pixPos) for pixPos in expBoxD.getCorners()]
        patchList = tractInfo.findPatchList(skyCorners)
        if not patchList:
            raise RuntimeError("No suitable tract found")
        self.log.info("Assembling %s coadd patches" % (len(patchList),))
        # compute inclusive bounding box
        coaddBBox = afwGeom.Box2I()
        for patchInfo in patchList:
            outerBBox = patchInfo.getOuterBBox()
            for corner in outerBBox.getCorners():
                coaddBBox.include(corner)
        self.log.info("exposure dimensions=%s; coadd dimensions=%s" % \
            (exposure.getDimensions(), coaddBBox.getDimensions()))
        
        coaddExposure = afwImage.ExposureF(coaddBBox, tractInfo.getWcs())
        edgeMask = afwImage.MaskU.getPlaneBitMask("EDGE")
        coaddExposure.getMaskedImage().set(num.nan, edgeMask, num.nan)
        nPatchesFound = 0
        for patchInfo in patchList:
            patchArgDict = dict(
                datasetType = self.config.coaddName + "Coadd",
                tract = tractInfo.getId(),
                patch = "%s,%s" % (patchInfo.getIndex()[0], patchInfo.getIndex()[1]),
            )
            if not sensorRef.datasetExists(**patchArgDict):
                self.log.warn("%(datasetType)s, tract=%(tract)s, patch=%(patch)s does not exist; skipping" % patchArgDict)
                continue

            nPatchesFound += 1
            self.log.info("Reading patch %s" % patchArgDict)
            coaddPatch = sensorRef.get(**patchArgDict)
            coaddView = afwImage.MaskedImageF(coaddExposure.getMaskedImage(),
                patchInfo.getOuterBBox(), afwImage.PARENT)
            coaddView <<= coaddPatch.getMaskedImage()
        
        if nPatchesFound == 0:
            raise RuntimeError("No patches found!")

        return coaddExposure

    def _getConfigName(self):
        """Return the name of the config dataset
        """
        return "%sDiff_config" % (self.config.coaddName,)
    
    def _getMetadataName(self):
        """Return the name of the metadata dataset
        """
        return "%sDiff_metadata" % (self.config.coaddName,)

    @classmethod
    def _makeArgumentParser(cls):
        """Create an argument parser
        """
        return pipeBase.ArgumentParser(name=cls._DefaultName, datasetType="calexp")
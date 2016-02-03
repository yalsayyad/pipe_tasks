#
# LSST Data Management System
# Copyright 2012-2016 LSST Corporation.
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

from .processImage import ProcessImageTask

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.afw.table as afwTable
import lsst.daf.base as dafBase
import lsst.meas.astrom as measAstrom
from lsst.meas.base import SingleFrameMeasurementTask
from lsst.meas.algorithms import SourceDetectionTask
from lsst.ip.diffim import  SourceFlagChecker, KernelCandidateF,\
     cast_KernelCandidateF, DipoleAnalysis


class ProcessDiffimConfig(ProcessImageTask.ConfigClass):
    """Config for ProcessDiffim"""
    doMerge = pexConfig.Field(dtype=bool, default=True,
        doc="Merge positive and negative diaSources with grow radius set by growFootprint")
    growFootprint = pexConfig.Field(dtype=int, default=2,
        doc="Grow positive and negative footprints by this amount before merging")

    doMatchSources = pexConfig.Field(dtype=bool, default=True,
        doc="Match diaSources with input calexp sources and ref catalog sources")
    doAddMetrics = pexConfig.Field(dtype=bool, default=True,
        doc="Add columns to the source table to hold analysis metrics?")
    doWriteSources = pexConfig.Field(dtype=bool, default=True, doc="Write sources?")
    detection = pexConfig.ConfigurableField(
        target=SourceDetectionTask,
        doc="Low-threshold detection for final measurement",
    )
    measurement = pexConfig.ConfigurableField(
        target=SingleFrameMeasurementTask,
        doc="Final source measurement on low-threshold detections; dipole fitting enabled.",
    )
    #REPEATED!
    doPreConvolve = pexConfig.Field(dtype=bool, default=False,
        doc="Convolve science image by its PSF before PSF-matching?")
    diaSourceMatchRadius = pexConfig.Field(dtype=float, default=0.5,
        doc="Match radius (in arcseconds) for DiaSource to Source association")

    #REPEATED
    coaddName = pexConfig.Field(
        doc="coadd name: typically one of deep or goodSeeing",
        dtype=str,
        default="deep",
    )

    def setDefaults(self):
        ProcessImageTask.ConfigClass.setDefaults(self)
        self.doCalibrate = False
        self.detection.thresholdPolarity = "both"
        self.detection.thresholdValue = 5.5
        self.detection.reEstimateBackground = False
        self.detection.thresholdType = "pixel_stdev"

        # Only run algorithms/plugins that make sense on an image difference
        self.measurement.plugins = ["base_PsfFlux",
                                    "base_CircularApertureFlux",
                                    "ip_diffim_NaiveDipoleCentroid",
                                    "ip_diffim_NaiveDipoleFlux",
                                    "ip_diffim_PsfDipoleFlux",
                                    "ip_diffim_ClassificationDipole",
                                    "base_SkyCoord",
                                    ]

        self.measurement.slots.calibFlux = None
        self.measurement.slots.modelFlux = None
        self.measurement.slots.instFlux = None
        self.measurement.slots.shape = None
        self.measurement.slots.centroid = "ip_diffim_NaiveDipoleCentroid"
        self.measurement.doReplaceWithNoise = False

        # Add filtered flux measurement, the correct measurement for pre-convolved images.
        # Enable all measurements, regardless of doPreConvolved, as it makes data harvesting easier.
        # To change that you must modify algorithms.names in the task's applyOverrides method,
        # after the user has set doPreConvolved.
        self.measurement.algorithms.names.add('base_PeakLikelihoodFlux')

        #self.detection.isotropicGrow = True
        #self.detection.returnOriginalFootprints = False
        #self.doWriteSourceMatches = True
        #self.measurement.doReplaceWithNoise = True
        # diffims do not yet have ap corr data; once they do, delete the following to enable the
        # ProcessImageConfig default of yes; meanwhile the warning may help remind us
        self.measurement.doApplyApCorr = "noButWarn"
        #self.deblend.maxNumberOfPeaks = 20
        #self.astrometry.forceKnownWcs = True

        #self.measurement.plugins.names |= ['base_InputCount']
        # The following line must be set if clipped pixel flags are to be added to the output table
        # The clipped mask plane is added by running SafeClipAssembleCoaddTask
        #self.measurement.plugins['base_PixelFlags'].masksFpAnywhere = ['CLIPPED']

        self.doPreConvolve = False
        self.doMatchSources = False
        self.doAddMetrics = False

class ProcessDiffimTask(ProcessImageTask):
    """Process a Image Difference
    """
    ConfigClass = ProcessDiffimConfig
    _DefaultName = "processDiffim"

    def __init__(self, schema=None, **kwargs):
        # Configuration isn't parsed until after initializing the command line task
        # Therefore we cannot pass schema to processImageTask constructor
        pipeBase.CmdLineTask.__init__(self, **kwargs)
        self.algMetadata = dafBase.PropertyList()
        if schema is None:
            schema = afwTable.SourceTable.makeMinimalSchema()
        self.schema = schema
        if self.config.doMatchSources:
            self.schema.addField("refMatchId", "L", "unique id of reference catalog match")
            self.schema.addField("srcMatchId", "L", "unique id of source match")

        if self.config.doDetection:
            self.makeSubtask("detection", schema=self.schema)
        if self.config.doDeblend:
            self.makeSubtask("deblend", schema=self.schema)
        if self.config.doMeasurement:
            self.makeSubtask("measurement", schema=self.schema, algMetadata=self.algMetadata)

        self.dataPrefix = self.config.coaddName

    def makeIdFactory(self, dataRef):
        # We make one IdFactory that will be used by both icSrc and src datasets;
        # I don't know if this is the way we ultimately want to do things, but at least
        # this ensures the source IDs are fully unique.
        expBits = dataRef.get("ccdExposureId_bits")
        expId = long(dataRef.get("ccdExposureId"))
        return afwTable.IdFactory.makeSource(expId, 64 - expBits)

    def getExposureId(self, dataRef):
        return long(dataRef.get(self.config.coaddName + "Diff_differenceExp"))

    def getAstrometer(self):
        return self.astrometry

    @pipeBase.timeMethod
    def run(self, dataRef):
        """Process Image Difference

        @param dataRef: butler data reference corresponding to coadd patch
        @return pipe_base Struct containing these fields:
        - exposure: input exposure, as modified in the course of processing
        - sources: detected source if config.doDetection, else None
        """

        if self.config.doDetection:
            self.log.info("Running diaSource detection")
            subtractedExposure = dataRef.get(self.config.coaddName + "Diff_differenceExp")

            # Get Psf from the appropriate input image if it doesn't exist
            if not subtractedExposure.hasPsf():
                raise RuntimeError("Image Difference does not have a PSF")
                """
                if self.config.convolveTemplate:# You can do without this!!!
                    subtractedExposure.setPsf(exposure.getPsf())
                else:
                    if templateExposure is None:
                        template = self.getTemplate.run(exposure, dataRef, templateIdList=templateIdList)
                    subtractedExposure.setPsf(template.exposure.getPsf())
                """

            idFactory = self.makeIdFactory(dataRef)

            # Erase existing detection mask planes
            mask  = subtractedExposure.getMaskedImage().getMask()
            mask &= ~(mask.getPlaneBitMask("DETECTED") | mask.getPlaneBitMask("DETECTED_NEGATIVE"))

            table = afwTable.SourceTable.make(self.schema, idFactory)
            table.setMetadata(self.algMetadata)
            results = self.detection.makeSourceCatalog(
                table=table,
                exposure=subtractedExposure,
                doSmooth=not self.config.doPreConvolve  #  can do without this too. 
                )

            if self.config.doMerge:
                fpSet = results.fpSets.positive
                fpSet.merge(results.fpSets.negative, self.config.growFootprint,
                            self.config.growFootprint, False)
                diaSources = afwTable.SourceCatalog(table)
                fpSet.makeSources(diaSources)
                self.log.info("Merging detections into %d sources" % (len(diaSources)))
            else:
                diaSources = results.sources

            if self.config.doMeasurement:
                self.log.info("Running diaSource measurement")
                self.measurement.run(diaSources, subtractedExposure)

            # Match with the calexp sources if possible
            if self.config.doMatchSources:
                if dataRef.datasetExists("src"):
                    # Create key,val pair where key=diaSourceId and val=sourceId
                    matchRadAsec = self.config.diaSourceMatchRadius
                    exposure = dataRef.get('calexp')
                    matchRadPixel = matchRadAsec / exposure.getWcs().pixelScale().asArcseconds()
                    srcMatches = afwTable.matchXy(dataRef.get("src"), diaSources, matchRadPixel, True)
                    srcMatchDict = dict([(srcMatch.second.getId(), srcMatch.first.getId()) for 
                                         srcMatch in srcMatches])
                    self.log.info("Matched %d / %d diaSources to sources" % (len(srcMatchDict),
                                                                             len(diaSources)))
                else:
                    self.log.warn("Src product does not exist; cannot match with diaSources")
                    srcMatchDict = {}

                # Create key,val pair where key=diaSourceId and val=refId
                refAstromConfig = measAstrom.AstrometryConfig()
                refAstromConfig.matcher.maxMatchDistArcSec = matchRadAsec
                refAstrometer = measAstrom.AstrometryTask(refAstromConfig)
                astromRet = refAstrometer.run(exposure=exposure, sourceCat=diaSources)
                refMatches = astromRet.matches
                if refMatches is None:
                    self.log.warn("No diaSource matches with reference catalog")
                    refMatchDict = {}
                else:
                    self.log.info("Matched %d / %d diaSources to reference catalog" % (len(refMatches),
                                                                                       len(diaSources)))
                    refMatchDict = dict([(refMatch.second.getId(), refMatch.first.getId()) for \
                                             refMatch in refMatches])

                # Assign source Ids
                for diaSource in diaSources:
                    sid = diaSource.getId()
                    if srcMatchDict.has_key(sid):
                        diaSource.set("srcMatchId", srcMatchDict[sid])
                    if refMatchDict.has_key(sid):
                        diaSource.set("refMatchId", refMatchDict[sid])

            if diaSources is not None and self.config.doWriteSources:
                dataRef.put(diaSources, self.config.coaddName + "Diff_diaSrc")

    @classmethod
    def _makeArgumentParser(cls):
        parser = pipeBase.ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "calexp", help="argh!")
        return parser

    def _getConfigName(self):
        """Return the name of the config dataset
        """
        return "%s_processDiffim_config" % (self.config.coaddName,)

    def _getMetadataName(self):
        """Return the name of the metadata dataset
        """
        return "%s_processDiffim_metadata" % (self.config.coaddName,)

    def getSchemaCatalogs(self):
        """Return a dict of empty catalogs for each catalog dataset produced by this task."""
        diaSrc = afwTable.SourceCatalog(self.schema)
        diaSrc.getTable().setMetadata(self.algMetadata)
        return {self.config.coaddName + "Diff_diaSrc": diaSrc}

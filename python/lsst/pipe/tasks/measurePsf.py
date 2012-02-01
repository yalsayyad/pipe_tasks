# 
# LSST Data Management System
# Copyright 2008, 2009, 2010, 2011 LSST Corporation.
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the LSST License Statement and 
# the GNU General Public License along with this program.  If not, 
# see <http://www.lsstcorp.org/LegalNotices/>.
#
import lsst.afw.detection as afwDet
import lsst.meas.algorithms as measAlg
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

def propagateFlag(flag, old, new):
    """Propagate a flag from one source to another"""
    if old.getFlagForDetection() & flag:
        new.setFlagForDetection(new.getFlagForDetection() | flag)

class MeasurePsfConfig(pexConfig.Config):
    starSelector = measAlg.starSelectorRegistry.makeField("Star selection algorithm")
    psfDeterminer = measAlg.psfDeterminerRegistry.makeField("PSF Determination algorithm")

class MeasurePsfTask(pipeBase.Task):
    """Conversion notes:
    
    Split out of Calibrate since it seemed a good self-contained task
    
    @warning
    - I'm not sure I'm using metadata correctly (to replace old sdqa code)
    - The star selector and psf determiner registries will have to be modified to return a class,
      which has a ConfigClass attribute and can be instantiated with a config. Until then, there's no
      obvious way to get a registry algorithm's Config from another Config.
    """
    ConfigClass = MeasurePsfConfig

    def __init__(self, *args, **kwargs):
        pipeBase.Task.__init__(self, *args, **kwargs)
        self.starSelector = self.config.starSelector.apply()
        self.psfDeterminer = self.config.psfDeterminer.apply()
        
    @pipeBase.timeMethod
    def run(self, exposure, sources):
        """Measure the PSF

        @param exposure Exposure to process
        @param sources Measured sources on exposure
        """
        assert exposure, "No exposure provided"
        assert sources, "No sources provided"
        self.log.log(self.log.INFO, "Measuring PSF")

        #
        # Run an extra detection step to mask out faint stars
        #
        if False:
            print "RHL is cleaning faint sources"

            import lsst.afw.math as afwMath

            sigma = 1.0
            gaussFunc = afwMath.GaussianFunction1D(sigma)
            gaussKernel = afwMath.SeparableKernel(15, 15, gaussFunc, gaussFunc)

            im = exposure.getMaskedImage().getImage()
            convolvedImage = im.Factory(im.getDimensions())
            afwMath.convolve(convolvedImage, im, gaussKernel)
            del im

            fs = afwDet.makeFootprintSet(convolvedImage, afwDet.createThreshold(4, "stdev"))
            fs = afwDet.makeFootprintSet(fs, 3, True)
            fs.setMask(exposure.getMaskedImage().getMask(), "DETECTED")

        psfCandidateList = self.starSelector.selectStars(exposure, sources)

        psf, cellSet = self.psfDeterminer.determinePsf(exposure, psfCandidateList, self.metadata)
        self.log.log(self.log.INFO, "PSF determination using %d/%d stars." % 
                     (self.metadata.get("numGoodStars"), self.metadata.get("numAvailStars")))

        # The PSF candidates contain a copy of the source, and so we need to explicitly propagate new flags
        for cand in psfCandidateList:
            cand = measAlg.cast_PsfCandidateF(cand)
            src = cand.getSource()
            if src.getFlagForDetection() & measAlg.Flags.PSFSTAR:
                ident = src.getId()
                src = sources[ident]
                assert src.getId() == ident
                src.setFlagForDetection(src.getFlagForDetection() | measAlg.Flags.PSFSTAR)

        exposure.setPsf(psf)
        return pipeBase.Struct(
            psf = psf,
            cellSet = cellSet,
        )

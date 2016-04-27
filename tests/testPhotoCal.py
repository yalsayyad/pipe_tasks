#!/usr/bin/env python

#
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
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

import os
import unittest

import numpy as np

import lsst.meas.astrom            as measAstrom
import lsst.afw.geom               as afwGeom
import lsst.afw.table              as afwTable
import lsst.afw.image              as afwImage
import lsst.utils.tests            as utilsTests
from lsst.log import Log
from lsst.pipe.tasks.photoCal import PhotoCalTask, PhotoCalConfig

import testFindAstrometryNetDataDir as helper

# Quiet down meas_astrom logging, so we can see PhotoCal logs better
Log.setLevel("meas.astrom.astrometry_net", Log.WARN)
Log.setLevel("meas.astrom.sip", Log.WARN)
Log.setLevel("astrometricSolver", Log.WARN)

class PhotoCalTest(unittest.TestCase):

    def setUp(self):
        self.conf = measAstrom.AstrometryConfig()

        # Load sample input from disk
        testDir=os.path.dirname(__file__)
        self.srcCat = afwTable.SourceCatalog.readFits(os.path.join(testDir, "data", "v695833-e0-c000.xy.fits"))
        self.srcCat["slot_ApFlux_fluxSigma"] = 1
        self.srcCat["slot_PsfFlux_fluxSigma"] = 1

        # The .xy.fits file has sources in the range ~ [0,2000],[0,4500]
        # which is bigger than the exposure
        self.bbox = afwGeom.Box2I(afwGeom.Point2I(0,0), afwGeom.Extent2I(2048, 4612))
        smallExposure = afwImage.ExposureF(os.path.join(testDir, "data", "v695833-e0-c000-a00.sci.fits"))
        self.exposure = afwImage.ExposureF(self.bbox)
        self.exposure.setWcs(smallExposure.getWcs())
        self.exposure.setFilter(smallExposure.getFilter())
        self.exposure.setCalib(smallExposure.getCalib())

        # Set up local astrometry_net_data
        helper.setupAstrometryNetDataDir('photocal', rootDir=testDir)

    def tearDown(self):
        del self.srcCat
        del self.conf
        del self.exposure

    def getAstrometrySolution(self, loglvl = Log.INFO):
        astromConfig = measAstrom.AstrometryTask.ConfigClass()
        astrom = measAstrom.AstrometryTask(config=astromConfig)
        # use solve instead of run because the exposure has the wrong bbox
        return astrom.solve(
            sourceCat = self.srcCat,
            exposure = self.exposure,
        )

    def testGetSolution(self):
        res = self.getAstrometrySolution(loglvl=Log.DEBUG)
        self.assertTrue(res is not None)
        self.assertGreater(len(res.matches), 50)

    def test1(self):
        res = self.getAstrometrySolution()
        matches = res.matches

        print 'Test1'

        logLevel = Log.DEBUG
        log = Log.getLogger('testPhotoCal')
        Log.setLevel(log, logLevel)

        schema = matches[0].second.schema

        config = PhotoCalConfig()

        # The test and associated data have been prepared on the basis that we
        # use the PsfFlux to perform photometry.
        config.fluxField = "base_PsfFlux_flux"

        config.doWriteOutput = False    # schema is fixed because we already loaded the data
        task = PhotoCalTask(config=config, schema=schema)
        pCal = task.run(exposure=self.exposure, matches=matches)
        print "Ref flux fields list =", pCal.arrays.refFluxFieldList
        refFluxField = pCal.arrays.refFluxFieldList[0]

        # These are *all* the matches; we don't really expect to do that well.
        diff=[]
        for m in matches:
            refFlux = m[0].get(refFluxField) # reference catalog flux
            if refFlux <= 0:
                continue
            refMag = afwImage.abMagFromFlux(refFlux) # reference catalog mag
            instFlux = m[1].getPsfFlux()    #Instrumental Flux
            if instFlux <= 0:
                continue
            instMag = pCal.calib.getMagnitude(instFlux)     #Instrumental mag
            diff.append(instMag - refMag)
        diff = np.array(diff)

        self.assertGreater(len(diff), 50)
        log.info('%i magnitude differences; mean difference %g; mean abs diff %g' %
                 (len(diff), np.mean(diff), np.mean(np.abs(diff))))
        self.assertLess(np.mean(diff), 0.6)

        # Differences of matched objects that were used in the fit.
        zp = pCal.calib.getMagnitude(1.)
        log.debug('zeropoint: %g' % zp)
        fitdiff = pCal.arrays.srcMag + zp - pCal.arrays.refMag
        log.debug('number of sources used in fit: %i' % len(fitdiff))
        log.debug('rms diff: %g' % np.mean(fitdiff**2)**0.5)
        log.debug('median abs(diff): %g' % np.median(np.abs(fitdiff)))

        # zeropoint: 31.3145
        # number of sources used in fit: 65
        # median diff: -0.009681
        # mean diff: 0.00331871
        # median abs(diff): 0.0368904
        # mean abs(diff): 0.0516589

        self.assertLess(abs(zp - 31.3145), 0.05)

        self.assertGreater(len(fitdiff), 50)
        # Tolerances are somewhat arbitrary; they're set simply to avoid regressions, and
        # are not based on we'd expect to get given the data quality.
        self.assertLess(np.mean(fitdiff**2)**0.5, 0.07)    # rms difference
        self.assertLess(np.median(np.abs(fitdiff)), 0.06)  # median absolution difference

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(PhotoCalTest)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)

if __name__ == "__main__":
    run(True)

#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2015 AURA/LSST.
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
# see <https://www.lsstcorp.org/LegalNotices/>.
#

"""
Test the basic mechanics of coaddition, coadd processing, and forced photometry.

In this test, we build a mock calexps using perfectly knowns WCSs, with the only sources
being stars created from a perfectly known PSF, then coadd them, process the coadd (using
the new measurement framework in meas_base), and then run forced photometry (again, using
the new forced measurement tasks in meas_base).

We do not check that the results of this processing is exactly what we'd expect, except in
some cases where it's easy and/or particularly important to do so (e.g. CoaddPsf); we mostly
just check that everything runs, and that the results make enough sense to let us proceed
to the next step.

NOTE: if this test fails with what looks like a failure to load a FITS file, try changing
the REUSE_DATAREPO variable below to False, as sometimes this error message indicates a
different problem that's revealed when we're not trying to cache the mock data between
tests (but set REUSE_DATAREPO back to True when done debugging, or this test will be very
slow).
"""

import unittest
import shutil
import os
import sys
import numbers

from lsst.utils import getPackageDir
import lsst.utils.tests
import lsst.afw.math
import lsst.afw.geom
import lsst.afw.image
import lsst.afw.table.io
import lsst.meas.algorithms
import lsst.pipe.tasks.mocks
import lsst.daf.persistence

try:
    import lsst.meas.base
except ImportError:
    print "meas_base could not be imported; skipping this test"
    sys.exit(0)

from lsst.pipe.tasks.assembleCoadd import AssembleCoaddConfig, SafeClipAssembleCoaddConfig
from lsst.pipe.tasks.multiBand import (DetectCoaddSourcesTask, MergeDetectionsTask,
                                       MeasureMergedCoaddSourcesTask, MergeMeasurementsTask)

REUSE_DATAREPO = True  # Retain the data repo for each test? This greatly speeds up the tests, but may
# conflate or perhaps even hide some errors. If you get suspicious results, try setting it to False.
SAVE_DATAREPO = False  # Retain data repo after the tests succeed? Only set it True for debugging.
# Warning: even if True, the repository is deleted every time this test is run, so if you want
# to keep a copy safe, be sure to move it somewhere else.

PIPE_TASKS_DIR = getPackageDir("pipe_tasks")
DATAREPO_ROOT = os.path.join(PIPE_TASKS_DIR, "tests", ".tests", "testCoadds-data")

if os.path.exists(DATAREPO_ROOT):
    print("Deleting existing repo: %r" % (DATAREPO_ROOT,))
    shutil.rmtree(DATAREPO_ROOT)


class CoaddsTestCase(lsst.utils.tests.TestCase):

    def setupClass(self):
        """Delete the data repo if it exists"""

    def setUp(self):

        # Create a task that creates simulated images and builds a coadd from them
        self.mocksTask = lsst.pipe.tasks.mocks.MockCoaddTask()

        # Create an instance of DetectCoaddSourcesTask to measure on the coadd.
        # There's no noise in these images, so we set a direct-value threshold,
        # and the background weighting (when using Approximate) to False

        detectConfig = DetectCoaddSourcesTask.ConfigClass()
        detectConfig.detection.thresholdType = "value"
        detectConfig.detection.thresholdValue = 0.01
        detectConfig.detection.background.weighting = False
        self.detectTask = DetectCoaddSourcesTask(config=detectConfig)

        if REUSE_DATAREPO and os.path.exists(os.path.join(DATAREPO_ROOT, "_mapper")):
            self.butler = lsst.daf.persistence.Butler(DATAREPO_ROOT)
        else:
            self.butler = lsst.pipe.tasks.mocks.makeDataRepo(DATAREPO_ROOT)

            self.mocksTask.buildAllInputs(self.butler)
            self.addMaskPlanes()
            self.checkMaskPlanesExist()
            self.mocksTask.buildCoadd(self.butler)
            self.mocksTask.buildMockCoadd(self.butler)
            self.detectTask.writeSchemas(self.butler)
            # Now run the seperate multiband tasks on the Coadd to make the reference
            # catalog for the forced photometry tests.
            self.runTaskOnPatches(self.detectTask)

            mergeDetConfig = MergeDetectionsTask.ConfigClass()
            mergeDetConfig.priorityList = ['r', ]
            mergeDetTask = MergeDetectionsTask(config=mergeDetConfig, butler=self.butler)
            mergeDetTask.writeSchemas(self.butler)
            self.runTaskOnPatchList(mergeDetTask)

            measMergedConfig = MeasureMergedCoaddSourcesTask.ConfigClass()
            measMergedConfig.measurement.slots.shape = "base_SdssShape"
            measMergedConfig.measurement.doApplyApCorr = "no"
            measMergedConfig.measurement.plugins['base_PixelFlags'].masksFpAnywhere = []
            measMergedConfig.propagateFlags.flags = {}  # Disable flag propagation: no flags to propagate
            measMergedConfig.doMatchSources = False  # We don't have a reference catalog available
            measMergedTask = MeasureMergedCoaddSourcesTask(config=measMergedConfig, butler=self.butler)
            measMergedTask.writeSchemas(self.butler)
            self.runTaskOnPatches(measMergedTask)

            mergeMeasConfig = MergeMeasurementsTask.ConfigClass()
            mergeMeasConfig.priorityList = ['r', ]
            mergeMeasTask = MergeMeasurementsTask(config=mergeMeasConfig, butler=self.butler)
            mergeMeasTask.writeSchemas(self.butler)
            self.runTaskOnPatchList(mergeMeasTask)

    def tearDown(self):
        if not REUSE_DATAREPO:
            shutil.rmtree(DATAREPO_ROOT)
        del self.mocksTask
        del self.detectTask
        del self.butler

    def getCalexpIds(self, butler):
        if len(butler._repos.inputs()) is not 1:
            raise RuntimeError("This test assumes one input repository, can not continue with %s"
                               % butler.repository.inputs)

        return butler._repos.inputs()[0].repo._mapper.index['calexp']['visit'][None]

    def addMaskPlanes(self):
        butler = lsst.daf.persistence.Butler(DATAREPO_ROOT)
        # Get the dataId for each calexp in the repository
        calexpDataIds = self.getCalexpIds(butler)
        # Loop over each of the calexp and add the CROSSTALK and NOT_DEBLENDED mask planes
        for Id in calexpDataIds:
            image = butler.get('calexp', Id)
            mask = image.getMaskedImage().getMask()
            mask.addMaskPlane("CROSSTALK")
            mask.addMaskPlane("NOT_DEBLENDED")
            butler.put(image, 'calexp', dataId=Id)

    def checkMaskPlanesExist(self):
        butler = lsst.daf.persistence.Butler(DATAREPO_ROOT)
        # Get the dataId for each calexp in the repository
        calexpDataIds = self.getCalexpIds(butler)
        # Loop over each Id and verify the mask planes were added
        for ID in calexpDataIds:
            image = butler.get('calexp', ID)
            mask = image.getMaskedImage().getMask()
            self.assert_('CROSSTALK' in mask.getMaskPlaneDict().keys())
            self.assert_('NOT_DEBLENDED' in mask.getMaskPlaneDict().keys())

    def runTaskOnPatches(self, task, tract=0):
        skyMap = self.butler.get(self.mocksTask.config.coaddName + "Coadd_skyMap", immediate=True)
        tractInfo = skyMap[tract]
        for dataRef in self.mocksTask.iterPatchRefs(self.butler, tractInfo):
            task.run(dataRef)

    def runTaskOnPatchList(self, task, tract=0):
        skyMap = self.butler.get(self.mocksTask.config.coaddName + "Coadd_skyMap", immediate=True)
        tractInfo = skyMap[tract]
        for dataRef in self.mocksTask.iterPatchRefs(self.butler, tractInfo):
            task.run([dataRef])

    def runTaskOnCcds(self, task, tract=0):
        catalog = self.butler.get("observations", tract=tract, immediate=True)
        visitKey = catalog.getSchema().find("visit").key
        ccdKey = catalog.getSchema().find("ccd").key
        for record in catalog:
            dataRef = self.butler.dataRef("forced_src", tract=tract, visit=record.getI(visitKey),
                                          ccd=record.getI(ccdKey))
            task.run(dataRef)

    def getObsDict(self, tract=0):
        catalog = self.butler.get("observations", tract=tract, immediate=True)
        visitKey = catalog.getSchema().find("visit").key
        ccdKey = catalog.getSchema().find("ccd").key
        obsDict = {}
        for record in catalog:
            visit = record.getI(visitKey)
            ccd = record.getI(ccdKey)
            obsDict.setdefault(visit, {})[ccd] = record
        return obsDict

    def comparePsfs(self, a, b):
        if a is None and b is None:
            return
        ak = lsst.meas.algorithms.KernelPsf.swigConvert(a).getKernel()
        bk = lsst.meas.algorithms.KernelPsf.swigConvert(b).getKernel()
        self.assertEqual(type(ak), type(bk))
        self.assertEqual(ak.getDimensions(), bk.getDimensions())
        self.assertEqual(ak.getNKernelParameters(), ak.getNKernelParameters())
        self.assertEqual(ak.getNSpatialParameters(), ak.getNSpatialParameters())
        for aFuncParams, bFuncParams in zip(ak.getSpatialParameters(), bk.getSpatialParameters()):
            for aParam, bParam in zip(aFuncParams, bFuncParams):
                self.assertEqual(aParam, bParam)

    @unittest.skip("Remove test until DM-5174 is complete")
    def testMasksRemoved(self):
        image = self.butler.get(self.mocksTask.config.coaddName + 'Coadd_mock',
                                {'filter': 'r', 'tract': 0, 'patch': '0,0'})
        keys = image.getMaskedImage().getMask().getMaskPlaneDict().keys()
        self.assert_('CROSSTALK' not in keys)
        self.assert_('NOT_DEBLENDED' not in keys)

    def testTempExpInputs(self, tract=0):
        skyMap = self.butler.get(self.mocksTask.config.coaddName + "Coadd_skyMap", immediate=True)
        tractInfo = skyMap[tract]
        for visit, obsVisitDict in self.getObsDict(tract).iteritems():
            foundOneTempExp = False
            for patchRef in self.mocksTask.iterPatchRefs(self.butler, tractInfo):
                try:
                    tempExp = patchRef.get(self.mocksTask.config.coaddName + "Coadd_tempExp", visit=visit,
                                           immediate=True)
                    foundOneTempExp = True
                except:
                    continue
                self.assertEqual(tractInfo.getWcs(), tempExp.getWcs())
                coaddInputs = tempExp.getInfo().getCoaddInputs()
                self.assertEqual(len(coaddInputs.visits), 1)
                visitRecord = coaddInputs.visits[0]
                self.assertEqual(visitRecord.getWcs(), tempExp.getWcs())
                self.assertEqual(visitRecord.getBBox(), tempExp.getBBox())
                self.assert_(len(coaddInputs.ccds) > 0)
                ccdKey = coaddInputs.ccds.getSchema().find("ccd").key
                for ccdRecord in coaddInputs.ccds:
                    ccd = ccdRecord.getI(ccdKey)
                    obsRecord = obsVisitDict[ccd]
                    self.assertEqual(obsRecord.getId(), ccdRecord.getId())
                    self.assertEqual(obsRecord.getWcs(), ccdRecord.getWcs())
                    self.assertEqual(obsRecord.getBBox(), ccdRecord.getBBox())
                    self.comparePsfs(obsRecord.getPsf(), ccdRecord.getPsf())
            self.assert_(foundOneTempExp)

    def testCoaddInputs(self, tract=0):
        skyMap = self.butler.get(self.mocksTask.config.coaddName + "Coadd_skyMap", immediate=True)
        tractInfo = skyMap[tract]
        obsCatalog = self.butler.get("observations", tract=tract, immediate=True)
        for patchRef in self.mocksTask.iterPatchRefs(self.butler, tractInfo):
            coaddExp = patchRef.get(self.mocksTask.config.coaddName + "Coadd", immediate=True)
            self.assertEqual(tractInfo.getWcs(), coaddExp.getWcs())
            coaddInputs = coaddExp.getInfo().getCoaddInputs()
            try:
                ccdVisitKey = coaddInputs.ccds.getSchema().find("visit").key
            except:
                print patchRef.dataId
                print coaddInputs.ccds.getSchema()
                raise
            for ccdRecord in coaddInputs.ccds:
                obsRecord = obsCatalog.find(ccdRecord.getId())
                self.assertEqual(obsRecord.getId(), ccdRecord.getId())
                self.assertEqual(obsRecord.getWcs(), ccdRecord.getWcs())
                self.assertEqual(obsRecord.getBBox(), ccdRecord.getBBox())
                self.comparePsfs(obsRecord.getPsf(), ccdRecord.getPsf())
                self.assert_(coaddInputs.visits.find(ccdRecord.getL(ccdVisitKey)) is not None)
            for visitRecord in coaddInputs.visits:
                nCcds = len([ccdRecord for ccdRecord in coaddInputs.ccds
                             if ccdRecord.getL(ccdVisitKey) == visitRecord.getId()])
                self.assert_(nCcds >= 1)
                self.assert_(nCcds <= 2)

    def testPsfInstallation(self, tract=0):
        skyMap = self.butler.get(self.mocksTask.config.coaddName + "Coadd_skyMap", immediate=True)
        tractInfo = skyMap[tract]
        for patchRef in self.mocksTask.iterPatchRefs(self.butler, tractInfo):
            coaddExp = patchRef.get(self.mocksTask.config.coaddName + "Coadd", immediate=True)
            ccdCat = coaddExp.getInfo().getCoaddInputs().ccds
            savedPsf = lsst.meas.algorithms.CoaddPsf.swigConvert(coaddExp.getPsf())
            newPsf = lsst.meas.algorithms.CoaddPsf(ccdCat, coaddExp.getWcs())
            self.assertEqual(savedPsf.getComponentCount(), len(ccdCat))
            self.assertEqual(newPsf.getComponentCount(), len(ccdCat))
            for n, record in enumerate(ccdCat):
                self.assert_(lsst.afw.table.io.comparePersistablePtrs(savedPsf.getPsf(n), record.getPsf()))
                self.assert_(lsst.afw.table.io.comparePersistablePtrs(newPsf.getPsf(n), record.getPsf()))
                self.assert_(lsst.afw.table.io.comparePersistablePtrs(savedPsf.getWcs(n), record.getWcs()))
                self.assert_(lsst.afw.table.io.comparePersistablePtrs(newPsf.getWcs(n), record.getWcs()))
                self.assertEqual(savedPsf.getBBox(n), record.getBBox())
                self.assertEqual(newPsf.getBBox(n), record.getBBox())

    def testCoaddPsf(self, tract=0):
        skyMap = self.butler.get(self.mocksTask.config.coaddName + "Coadd_skyMap", immediate=True)
        tractInfo = skyMap[tract]
        # Start by finding objects that never appeared on the edge of an image
        simSrcCat = self.butler.get("simsrc", tract=tract, immediate=True)
        simSrcSchema = simSrcCat.getSchema()
        objectIdKey = simSrcSchema.find("objectId").key
        centroidInBBoxKey = simSrcSchema.find("centroidInBBox").key
        partialOverlapKey = simSrcSchema.find("partialOverlap").key
        simSrcByObject = {}
        for simSrcRecord in simSrcCat:
            simSrcByObject.setdefault(simSrcRecord.getL(objectIdKey), []).append(simSrcRecord)
        pureObjectIds = set()  # set will contain objects that never appear on edges
        for objectId, simSrcRecords in simSrcByObject.iteritems():
            inAnyImages = False
            for simSrcRecord in simSrcRecords:
                if simSrcRecord.getFlag(centroidInBBoxKey):
                    if simSrcRecord.getFlag(partialOverlapKey):
                        break
                    inAnyImages = True
            else:  # only get here if we didn't break
                if inAnyImages:
                    pureObjectIds.add(objectId)

        truthCatalog = self.butler.get("truth", tract=tract, immediate=True)
        truthCatalog.sort()
        nTested = 0
        for patchRef in self.mocksTask.iterPatchRefs(self.butler, tractInfo):
            coaddExp = patchRef.get(self.mocksTask.config.coaddName + "Coadd", immediate=True)
            coaddWcs = coaddExp.getWcs()
            coaddPsf = coaddExp.getPsf()
            coaddBBox = lsst.afw.geom.Box2D(coaddExp.getBBox())
            for objectId in pureObjectIds:
                truthRecord = truthCatalog.find(objectId)
                position = coaddWcs.skyToPixel(truthRecord.getCoord())
                if not coaddBBox.contains(position):
                    continue
                try:
                    psfImage = coaddPsf.computeImage(position)
                except:
                    continue
                psfImageBBox = psfImage.getBBox()
                if not coaddExp.getBBox().contains(psfImageBBox):
                    continue
                starImage = lsst.afw.image.ImageF(coaddExp.getMaskedImage().getImage(),
                                                  psfImageBBox).convertD()
                starImage /= starImage.getArray().sum()
                psfImage /= psfImage.getArray().sum()
                residuals = lsst.afw.image.ImageD(starImage, True)
                residuals -= psfImage

                self.assertClose(starImage.getArray(), psfImage.getArray(), rtol=1E-3, atol=1E-2)
                nTested += 1
        if nTested == 0:
            print("WARNING: CoaddPsf test inconclusive (this can occur randomly, but very rarely; "
                  "first try running the test again)")

    def testForcedPhotCoaddTask(self):
        config = lsst.meas.base.ForcedPhotCoaddConfig()
        config.references.filter = 'r'
        config.measurement.doApplyApCorr = "no"
        task = lsst.meas.base.ForcedPhotCoaddTask(config=config, butler=self.butler)
        task.writeSchemas(self.butler)
        self.runTaskOnPatches(task)

    def testForcedPhotCcdTask(self):
        config = lsst.meas.base.ForcedPhotCcdConfig()
        config.references.filter = 'r'
        config.measurement.doApplyApCorr = "no"
        task = lsst.meas.base.ForcedPhotCcdTask(config=config, butler=self.butler)
        self.runTaskOnCcds(task)

    def testAlgMetadataOutput(self):
        """Test to see if algMetadata is persisted correctly from MeasureMergedCoaddSourcesTask.

        This test fails with a NotFoundError if the algorithm metadata is not persisted"""
        patchList = ['0,0', '0,1', '1,0', '1,1']
        for patch in patchList:
            cat = self.butler.get("deepCoadd_meas", filter='r', tract=0, patch=patch)
            meta = cat.getTable().getMetadata()
            for circApertureFluxRadius in meta.get('base_CircularApertureFlux_radii'):
                self.assertIsInstance(circApertureFluxRadius, numbers.Number)
            for nOffset in meta.get('NOISE_OFFSET'):
                self.assertIsInstance(nOffset, numbers.Number)
            for noiseSrc in meta.get('NOISE_SOURCE'):
                self.assertEqual(noiseSrc, 'measure')
            for noiseExpID in meta.get('NOISE_EXPOSURE_ID'):
                self.assertIsInstance(noiseExpID, numbers.Number)
            for noiseSeedMul in meta.get('NOISE_SEED_MULTIPLIER'):
                self.assertIsInstance(noiseSeedMul, numbers.Number)


class AssembleCoaddTestCase(lsst.utils.tests.TestCase):

    def testSafeClipConfig(self):
        # Test for DM-4797: ensure that AssembleCoaddConfig.setDefaults() is
        # run when SafeClipAssembleCoaddConfig.setDefaults() is run. This
        # simply sets the default value for badMaskPlanes.
        self.assertEqual(AssembleCoaddConfig().badMaskPlanes, SafeClipAssembleCoaddConfig().badMaskPlanes)


def suite():
    lsst.utils.tests.init()
    suites = []
    suites += unittest.makeSuite(CoaddsTestCase)
    suites += unittest.makeSuite(AssembleCoaddTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)


def run(shouldExit=False):
    status = lsst.utils.tests.run(suite(), False)
    if os.path.exists(DATAREPO_ROOT):
        if SAVE_DATAREPO:
            print >> sys.stderr, "SAVE_DATAREPO is True: saving %r" % (os.path.abspath(DATAREPO_ROOT),)
        elif status == 0:
            shutil.rmtree(DATAREPO_ROOT)
        else:
            # Do not delete the DATAREPO_ROOT if the test failed to allow for forensics
            print >> sys.stderr, "Tests failed: saving %r" % (os.path.abspath(DATAREPO_ROOT),)

    if shouldExit:
        sys.exit(status)

    return status

if __name__ == "__main__":
    run(True)

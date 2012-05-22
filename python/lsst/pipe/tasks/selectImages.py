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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#
"""Note: this code uses MySQLdb primarily because daf_persistence cannot call scisql.scisql_s2CPolyRegion
"""
import MySQLdb
from lsst.daf.persistence import DbAuth
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

__all__ = ["BaseSelectImagesTask", "BaseExposureInfo"]

class SelectImagesConfig(pexConfig.Config):
    """Config for BaseSelectImagesTask
    """
    host = pexConfig.Field(
        doc = "Database server host name",
        dtype = str,
    )
    port = pexConfig.Field(
        doc = "Database server port",
        dtype = int,
    )
    database = pexConfig.Field(
        doc = "Name of database",
        dtype = str,
    )
    maxExposures = pexConfig.Field(
        doc = "maximum exposures to select; intended for debugging; ignored if None",
        dtype = int,
        optional = True,
    )


class BaseExposureInfo(object):
    """Data about a selected exposure
    """
    def __init__(self, result):
        """Create exposure information from a query result from a db connection
        
        The object has the following fields:
        - dataId: data ID of exposure (a dict)
        - coordList: a list of corner coordinates of the exposure (list of afwCoord.IcrsCoord)
        plus any others items that are desired
        
        Subclasses must override _setData and getColumnNames.
        """
        self._ind = -1
        self._setData(result)

    @property
    def _nextInd(self):
        self._ind += 1
        return self._ind
    
    def _setData(self, result):
        """Set exposure information based on a query result from a db connection
        
        Must set at least the following fields:
        - dataId: data ID of exposure (a dict)
        - coordList: coordinates of the corner of the exposure (list of afwCoord.IcrsCoord)
        """
        raise NotImplementedError()

    @staticmethod
    def getColumnNames():
        """Set database query columns to be consistent with constructor
        
        For example:
        return "raftName, visit, ccdName, filterName, ra1, dec1, ra2, dec2, ra3, dec3, ra4, dec4"
        """
        raise NotImplementedError()


class BaseSelectImagesTask(pipeBase.Task):
    """Base task for selecting images suitable for coaddition
    """
    ConfigClass = SelectImagesConfig
    _DefaultName = "selectImages"
    
    @pipeBase.timeMethod
    def run(self, coordList):
        """Select images suitable for coaddition in a particular region
        
        @param[in] coordList: list of coordinates defining region of interest; if None then select all images
        subclasses may add additional keyword arguments, as required
        
        @return a pipeBase Struct containing:
        - exposureInfoList: a list of exposure information objects (subclasses of BaseExposureInfo),
            which have at least the following fields:
            - dataId: data ID dictionary
            - coordList: coordinates of the corner of the exposure (list of afwCoord.IcrsCoord)
        """
        raise NotImplementedError()
    
    def _runArgDictFromDataId(self, dataId):
        """Extract keyword arguments for run (other than coordList) from a data ID
        
        @return keyword arguments for run (other than coordList), as a dict
        """
        raise NotImplementedError()
    
    def runDataRef(self, dataRef, coordList):
        """Run based on a data reference
        
        @param[in] dataRef: data reference; must contain any extra keys needed by coordList
        @param[in] coordList: list of coordinates defining region of interest
        @return a pipeBase Struct containing:
        - dataRefList: a list of data references
        - exposureInfoList: a list of ccdInfo objects
        """
        butler = dataRef.butlerSubset.butler
        runArgDict = self._runArgDictFromDataId(dataRef.dataId)
        exposureInfoList = self.run(coordList, **runArgDict).exposureInfoList
        dataRefList = [butler.dataRef(
            datasetType = "calexp",
            dataId = ccdInfo.dataId,
        ) for ccdInfo in exposureInfoList]
        return pipeBase.Struct(
            dataRefList = dataRefList,
            exposureInfoList = exposureInfoList,
        )
    
    def searchWholeSky(self, dataRef):
        """Search the whole sky using a data reference
        @param[in] dataRef: data reference; must contain key "filter"
        @return a pipeBase Struct containing:
        - exposureInfoList: a list of ccdInfo objects
        """
        runArgDict = self._runArgDictFromDataId(dataRef.dataId)
        return self.run(coordList=None, **runArgDict)
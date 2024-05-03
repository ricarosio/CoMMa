# -*- coding: utf-8 -*-
#Authors: Riccardo Arosio and Joana Gafeira
#Institutions: University College Cork, British Geological Survey
#CoMMa version 1.2 - May 2024

import arcpy
import os
import numpy as np

from arcpy.sa import *


class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "DEM pre-processing and useful derivatives for the CoMMa tool"
        self.alias = "CoMMa_Auxiliaries"
        self.description = ""

        # List of tool classes associated with this toolbox
        self.tools = [Smoothing_filter, Filter_and_clip, Fencing, PI_calculator, MFil_calculator]


class helpers(object):
    
    # This function converts backslash (accepted through the ArcGIS tool) to forwardslash (needed in python script) in a path
    def convert_backslash_forwardslash(self,inText):
        # inText: input path
        
        inText = fr"{inText}"
        if inText.find('\t'):
            inText = inText.replace('\t', '\\t')
        elif inText.find('\n'):
            inText = inText.replace('\n', '\\n')
        elif inText.find('\r'):
            inText = inText.replace('\r', '\\r')

        inText = inText.replace('\\','/')
        return inText
    
    def cleaning(self,workspace):
        
        arcpy.env.workspace = workspace
        try:
            featureclasses = arcpy.ListFeatureClasses()
            for fc in featureclasses:
                arcpy.Delete_management(fc)
        except:
            arcpy.AddWarning('Could not delete temporary data ' + str(fc))
        try:
            rasters = arcpy.ListRasters()
            for ra in rasters:
                arcpy.Delete_management(ra)
        except:
            arcpy.AddWarning('Could not delete temporary data ' + str(ra))
        
        try:
            arcpy.Delete_management(workspace)
        except:
            arcpy.AddWarning('Could not delete temporary folder ' + str(workspace[:-1]))
            
        return


class PI_calculator(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Mean LTPs"
        self.description = "The Bathymetric Position Index (BPI) is a measure of \
            where a referenced location is relative to the locations surrounding it. \
            BPI is derived from a bathymetric dataset using a topographic position index (TPI) \
            modified by Andrew Weiss and applied in the Benthic Terrain Modeler in ArcGIS"
        self.canRunInBackground = False
        self.category = "Local Topographic Position (LTP) derivatives"

    def getParameterInfo(self):
        """Define parameter definitions"""
        parameters = None
        
        param0 = arcpy.Parameter(
            displayName="Input bathymetry raster",
            name="inputDEM",
            datatype="GPRasterLayer",
            parameterType="Required",
            direction="Input")
        
        param1 = arcpy.Parameter(
            displayName="Filter type",
            name="filt_type",
            datatype="GPString",
            parameterType="Required",
            direction="Input")
        param1.filter.type = "ValueList"
        param1.filter.list = ["BPI", "DEV"]
        
        param2 = arcpy.Parameter(
            displayName="Inner radius (unit: cell)",
            name="inner_radius",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")
        param2.value = 3
        
        param3 = arcpy.Parameter(
            displayName="Outer radius (unit: cell)",
            name="outer_radius",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")
        param3.value = 5
        
        param4 = arcpy.Parameter(
            displayName="Output position index",
            name="outRas",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Output")
        
        parameters = [param0, param1, param2, param3, param4]
        
        return parameters

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        
        if parameters[0].valueAsText is not None:
            folder = os.path.basename(os.path.dirname(parameters[0].valueAsText))
            #desc = arcpy.Describe(gdb)
            if folder.lower().endswith(('.gdb', '.mdb')):
                parameters[0].setErrorMessage("Geodatabases cannot be used in this version of the CoMMa Toolbox")
                
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        
        arcpy.env.overwriteOutput = True
        inputDEM = parameters[0].valueAsText
        filt_type = parameters[1].valueAsText
        inner_radius = parameters[2].valueAsText
        outer_radius = parameters[3].valueAsText
        outRas = parameters[4].valueAsText
        
        # enable the helper functions
        helper = helpers()
        
        inputDEM = helper.convert_backslash_forwardslash(inputDEM)
        outRas = helper.convert_backslash_forwardslash(outRas)
        
        if inputDEM.rfind("/") < 0:
            aprx = arcpy.mp.ArcGISProject("CURRENT")
            m = aprx.activeMap            
            for lyr in m.listLayers():
                if lyr.isRasterLayer:
                    if inputDEM == lyr.name: 
                        inputDEM = helper.convert_backslash_forwardslash(lyr.dataSource)
        

        if int(inner_radius) > 0:
            neighborhood = NbrAnnulus(inner_radius, outer_radius, "CELL")
            
        else:
            neighborhood = NbrCircle(outer_radius, "CELL")
        
        if filt_type == "BPI":
            out_focal_statistics = FocalStatistics(inputDEM, neighborhood, "MEAN")
            result_raster = Minus(inputDEM, out_focal_statistics)
        
        elif filt_type == "DEV":
            out_focal_statistics = FocalStatistics(inputDEM, neighborhood, "MEAN")
            std_ = FocalStatistics(inputDEM, neighborhood, "STD")
            result_raster = (inputDEM - out_focal_statistics) / std_
        
        result_raster.save(outRas) 
        
        return



class MFil_calculator(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Median LTPs"
        self.description = "This filer lets you choose between a simple median or a modified minimum, maximum or directional median \
        (minimum median was ideated originally by Adams et al. 2005). The modified medians first calculate a minimum \
        or maximum value surface running a convolution for a neighborhood defined as 1/4 of the \
        wavelength of interest. Subsequently they calculate the median value for the minimised/maximised surface. \
        Directional median calculation: it creates 16 median focal statistics using wedge neighborhoods, \
        then finds the median of the bowtie (opposite wedges), finally, merges the resulting rasters keeping the minimum \
        value for each pixel. Based on Kim and Wessel 2008 - computer-intensive with large kernel sizes"
        self.canRunInBackground = False
        self.category = "Local Topographic Position (LTP) derivatives"

    def getParameterInfo(self):
        """Define parameter definitions"""
        parameters = None
        
        param0 = arcpy.Parameter(
            displayName="Input bathymetry raster",
            name="inputDEM",
            datatype="GPRasterLayer",
            parameterType="Required",
            direction="Input")
        
        param1 = arcpy.Parameter(
            displayName="Filter type",
            name="filt_type",
            datatype="GPString",
            parameterType="Required",
            direction="Input")
        param1.filter.type = "ValueList"
        param1.filter.list = ["MEDIAN", "MIN MEDIAN", "MAX MEDIAN", "DIRECTIONAL MEDIAN"]
        
        param2 = arcpy.Parameter(
            displayName="Inner radius (unit: cell)",
            name="inner_radius",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")
        param2.value = 3
        
        param3 = arcpy.Parameter(
            displayName="Outer radius (unit: cell)",
            name="outer_radius",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")
        param3.value = 5
        
        param4 = arcpy.Parameter(
            displayName="Output Median Residual",
            name="outRas",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Output")
        
        parameters = [param0, param1, param2, param3, param4]
        
        return parameters

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        
        if parameters[1].value == "MEDIAN":
            
            parameters[2].enabled = True
        
        else:
            
            parameters[2].enabled = False
            
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        
        if parameters[0].valueAsText is not None:
            folder = os.path.basename(os.path.dirname(parameters[0].valueAsText))
            #desc = arcpy.Describe(gdb)
            if folder.lower().endswith(('.gdb', '.mdb')):
                parameters[0].setErrorMessage("Geodatabases cannot be used in this version of the CoMMa Toolbox")
                
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        
        arcpy.env.overwriteOutput = True
        inputDEM = parameters[0].valueAsText
        filt_type = parameters[1].valueAsText
        inner_radius = parameters[2].valueAsText
        outer_radius = parameters[3].valueAsText
        outRas = parameters[4].valueAsText
        
        # enable the helper functions
        helper = helpers()
        
        inputDEM = helper.convert_backslash_forwardslash(inputDEM)
        outRas = helper.convert_backslash_forwardslash(outRas)
        
        if inputDEM.rfind("/") < 0:
            aprx = arcpy.mp.ArcGISProject("CURRENT")
            m = aprx.activeMap            
            for lyr in m.listLayers():
                if lyr.isRasterLayer:
                    if inputDEM == lyr.name: 
                        inputDEM = helper.convert_backslash_forwardslash(lyr.dataSource)
        
        if int(inner_radius) > 0:
            medianWindow = NbrAnnulus(inner_radius, outer_radius, "CELL")
        else:
            medianWindow = NbrCircle(outer_radius, "CELL")
        
        if filt_type == "MEDIAN":
            fcMed = FocalStatistics(inputDEM,medianWindow,"MEDIAN")
            fcMed = inputDEM - fcMed
            
        elif filt_type == "MIN MEDIAN":
            quartWindow = int(int(outer_radius)/4)
            minimWindow = NbrCircle(quartWindow, "CELL")
            minimRas = FocalStatistics(inputDEM, minimWindow, "MINIMUM")
            fcMed = FocalStatistics(minimRas,medianWindow,"MEDIAN")
            fcMed = inputDEM - fcMed

        elif filt_type == "MAX MEDIAN":
            quartWindow = int(int(outer_radius)/4)
            maxWindow = NbrCircle(quartWindow, "CELL")
            maximRas = FocalStatistics(inputDEM, maxWindow, "MAXIMUM")
            fcMed = FocalStatistics(maximRas,medianWindow,"MEDIAN")
            fcMed = inputDEM - fcMed
        
        elif filt_type == "DIRECTIONAL MEDIAN":
            Diz_wedges = {}
            Diz_bowtie = {}
            c = 0
            a = [11.25, 33.75, 56.25, 78.75, 101.25, 123.75, 146.25, 168.75, 
                 191.25, 213.75, 236.25, 258.75, 281.25, 303.75, 326.25, 348.75]
            b = [348.75, 11.25, 33.75, 56.25, 78.75, 101.25, 123.75, 146.25, 
                 168.75, 191.25, 213.75, 236.25, 258.75, 281.25, 303.75, 326.25]
            
            n = [1,2,3,4,5,6,7,8]
            m = [9,10,11,12,13,14,15,16]
             
            angles = zip(a, b)
            couples = zip(n, m)
            
            for start, end in angles:
                c+=1
                wedge = NbrWedge(outer_radius, start, end, "CELL") 
                wedgeMed = FocalStatistics(inputDEM, wedge, "MEDIAN")
                
                Diz_wedges["wedge_{0}".format(c)] = wedgeMed
                
            c = 0 
            for n, m in couples:
                c+=1
                bowtie = (Diz_wedges.get("wedge_{0}".format(n)) + Diz_wedges.get("wedge_{0}".format(m))) / 2
                Diz_bowtie["bowtie_{0}".format(c)] = bowtie
            
            Bowtie_list = list(Diz_bowtie.values())
            fcMed = arcpy.ia.Merge(Bowtie_list, "MIN")
            fcMed = inputDEM - fcMed

        fcMed.save(outRas) 
        
        return
    

class Smoothing_filter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Smoothing Filters"
        self.description = "A set of useful filters to clean the data"
        self.canRunInBackground = False
        self.category = "Pre-processing"

    def getParameterInfo(self):
        """Define parameter definitions"""
        parameters = None
        
        param0 = arcpy.Parameter(
            displayName="Input bathymetry raster",
            name="inputDEM",
            datatype="GPRasterLayer",
            parameterType="Required",
            direction="Input")
        
        param1 = arcpy.Parameter(
            displayName="Filter type",
            name="filt_type",
            datatype="GPString",
            parameterType="Required",
            direction="Input")
        param1.filter.type = "ValueList"
        param1.filter.list = ["MEAN", "MEDIAN", "LOW-PASS 3x3", "LOW-PASS 5x5"]
        
        param2 = arcpy.Parameter(        
            displayName="Neighbourhood type",
            name="neigh_type",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        param2.filter.type = "ValueList"
        param2.filter.list = ["CIRCLE", "RECTANGLE", "ANNULUS", "WEDGE"]
        
        param3 = arcpy.Parameter(
            displayName="Distance 1 (unit: cell)",
            name="distance_1",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")
        param3.value = 3
        
        param4 = arcpy.Parameter(
            displayName="Distance 2 (unit: cell)",
            name="distance_2",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")
        
        param5 = arcpy.Parameter(
            displayName="Angle 1 (unit: deg)",
            name="angle_1",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")
        
        param6 = arcpy.Parameter(
            displayName="Angle 2 (unit: deg)",
            name="angle_2",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")
        
        param7 = arcpy.Parameter(
            displayName="Output Filtered",
            name="outRas",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Output")
        
        parameters = [param0, param1, param2, param3, param4, param5, param6, param7]
        
        return parameters

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        
        if parameters[1].value == "MEAN" or parameters[1].value == "MEDIAN":
            
            parameters[2].enabled = True
        
        else:
            
            parameters[2].enabled = False
            
        
        if parameters[2].enabled == True and parameters[2].value == "CIRCLE": 
            
            parameters[3].enabled = True
            parameters[4].enabled = False
            parameters[5].enabled = False
            parameters[6].enabled = False
        
        elif parameters[2].enabled == True and parameters[2].value == "WEDGE": 
            
            parameters[3].enabled = True
            parameters[4].enabled = False
            parameters[5].enabled = True
            parameters[6].enabled = True
           
        elif parameters[2].enabled == True and [parameters[2].value == "RECTANGLE" 
                                                or parameters[2].value == "ANNULUS"]:
           parameters[3].enabled = True
           parameters[4].enabled = True
           parameters[5].enabled = False
           parameters[6].enabled = False
           
        else:
            
            parameters[3].enabled = False
            parameters[4].enabled = False
            parameters[5].enabled = False
            parameters[6].enabled = False
        
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        
        if parameters[0].valueAsText is not None:
            folder = os.path.basename(os.path.dirname(parameters[0].valueAsText))
            #desc = arcpy.Describe(gdb)
            if folder.lower().endswith(('.gdb', '.mdb')):
                parameters[0].setErrorMessage("Geodatabases cannot be used in this version of the CoMMa Toolbox")
                
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        
        arcpy.env.overwriteOutput = True
        inputDEM = parameters[0].valueAsText
        filt_type = parameters[1].valueAsText
        neigh_type = parameters[2].valueAsText
        distance_1 = parameters[3].valueAsText
        distance_2 = parameters[4].valueAsText
        angle_1 = parameters[5].valueAsText
        angle_2 = parameters[6].valueAsText
        outRas = parameters[7].valueAsText
        
        # enable the helper functions
        helper = helpers()
        
        inputDEM = helper.convert_backslash_forwardslash(inputDEM)
        outRas = helper.convert_backslash_forwardslash(outRas)
        
        if inputDEM.rfind("/") < 0:
            aprx = arcpy.mp.ArcGISProject("CURRENT")
            m = aprx.activeMap            
            for lyr in m.listLayers():
                if lyr.isRasterLayer:
                    if inputDEM == lyr.name: 
                        inputDEM = helper.convert_backslash_forwardslash(lyr.dataSource)
        
        if filt_type == "MEAN" or filt_type == "MEDIAN":
            if neigh_type == "ANNULUS":
                Window = NbrAnnulus(distance_1, distance_2, "CELL")
            elif neigh_type =="CIRCLE":
                Window = NbrCircle(distance_1, "CELL")
            elif neigh_type =="RECTANGLE":
                Window = NbrRectangle(distance_1, distance_2, "CELL")
            else:
                Window = NbrWedge(distance_1, angle_1, angle_2, "CELL")
        else:
            pass
        
        if filt_type == "MEAN":
            fc = FocalStatistics(inputDEM, Window, "MEAN")
            
        elif filt_type == "MEDIAN":
            fc = FocalStatistics(inputDEM, Window, "MEDIAN")
            
        elif filt_type == "LOW-PASS 3x3":
            fc = Convolution(inputDEM, 11)
            
        elif filt_type == "LOW-PASS 5x5":
            fc = Convolution(inputDEM, 12)
            
        fc.save(outRas) 
        
        return
    

class Filter_and_clip(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Filter and Clip tool"
        self.description = "In areas with complex seabed morphology or features connected by regional morphology \
            (e.g. pockmarks in basins), it may be needed to prepare the data to allow the best results when running \
            an automated mapping tool. Both a High and a Low Pass filter are applied to the original DEM. \
            These filters allow to highlight certain features in the data and to exclude areas of the DEM with gentle local variations."
            
        self.canRunInBackground = False
        self.category = "Pre-processing"

    def getParameterInfo(self):
        """Define parameter definitions"""
        parameters = None
        
        param0 = arcpy.Parameter(
            displayName="Input raster",
            name="inputDEM",
            datatype="GPRasterLayer",
            parameterType="Required",
            direction="Input")
        
        param1 = arcpy.Parameter(
            displayName="Filter threshold for reclassification",
            name="Filter_t",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        param1.value = 0.2
        
        param2 = arcpy.Parameter(
            displayName="Buffer size (unit: m)",
            name="Buffer_s",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        
        param3 = arcpy.Parameter(
            displayName="Workspace",
            name="workspace",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")  
        
        param4 = arcpy.Parameter(
            displayName="Output filtered raster name",
            name="outRas",
            datatype="GPString",
            parameterType="Required",
            direction="Output")
        
        param5 = arcpy.Parameter(
            displayName="Do you want to delete the temporary files?",
            name="delTemp",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")
        
        parameters = [param0, param1, param2, param3, param4, param5]
        
        return parameters

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        
        if parameters[0].valueAsText is not None:
            folder = os.path.basename(os.path.dirname(parameters[0].valueAsText))
            #desc = arcpy.Describe(gdb)
            if folder.lower().endswith(('.gdb', '.mdb')):
                parameters[0].setErrorMessage("Geodatabases cannot be used in this version of the CoMMa Toolbox")
        
        if parameters[3].valueAsText is not None:
            folder = os.path.basename(parameters[3].valueAsText)
            #desc = arcpy.Describe(gdb)
            if folder.lower().endswith(('.gdb', '.mdb')):
                parameters[3].setErrorMessage("Geodatabases cannot be used in this version of the CoMMa Toolbox")
                
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        
        arcpy.env.overwriteOutput = True
        inputDEM = parameters[0].valueAsText
        Filter_t = parameters[1].valueAsText
        Buffer_s = parameters[2].valueAsText
        workspace = parameters[3].valueAsText
        outRas = parameters[4].valueAsText
        delTemp = parameters[5].valueAsText
        
        
        if outRas[-4:] != '.tif':
            outRas = outRas + '.tif'
        # enable the helper functions
        helper = helpers()
        
        inputDEM = helper.convert_backslash_forwardslash(inputDEM)
        workspace = helper.convert_backslash_forwardslash(workspace)
        
        # if the input bathyRas is selected from a drop-down list, the bathyRas does not contain the full path
        # In this case, the full path needs to be obtained from the map layer
        if inputDEM.rfind("/") < 0:
            aprx = arcpy.mp.ArcGISProject("CURRENT")
            m = aprx.activeMap            
            for lyr in m.listLayers():
                if lyr.isRasterLayer:
                    if inputDEM == lyr.name: 
                        inputDEM = helper.convert_backslash_forwardslash(lyr.dataSource)
        
        Cellsize_Xresult = arcpy.GetRasterProperties_management(inputDEM, "CELLSIZEX")
        Cellsize_X = float(Cellsize_Xresult.getOutput(0))
        cutoffArea = (Cellsize_X*Cellsize_X)*3.0
        
        tempWS = workspace + "\\temp_fclip\\"
        arcpy.env.overwriteOutput = True
        if not os.path.exists(tempWS):
            os.mkdir(tempWS)
              
        arcpy.env.qualifiedFieldNames = False
        arcpy.env.workspace = workspace
        
        # Process: Filter; High Pass Filter applied to input DEM
        arcpy.AddMessage('Running the high pass filter ...')
        Filter1out = abs(Filter(inputDEM, "HIGH", "NODATA"))
        Filter1out.save(tempWS + "1H")
        
        # Process: Filter (2); Low Pass Filter applied to DEM created by High Pass Filter
        arcpy.AddMessage('Running the low pass filter ...')
        Filter2out = Filter(Filter1out, "LOW", "NODATA")
        Filter2out.save(tempWS + "2HL")
        
        # Process: Reclassify (user define Threshold reclassify range) in raster; -1 and lower classed as 1,
        # values greater than -1 classed as no data
        reclassifyRange = Filter_t + " 999 1"
        reclassed = Reclassify(tempWS + "2HL", "Value", reclassifyRange, "NODATA")
               
        reclassed.save(tempWS + "3Rc")
        # Process: Raster to Polygon
        arcpy.RasterToPolygon_conversion(tempWS + "3Rc", tempWS + "4P.shp",
                                         "NO_SIMPLIFY", "")
        
        # Add field with AREA values
        arcpy.AddField_management(tempWS + "4P.shp", "Area", "DOUBLE", "8", "1")
        areaExpression = "float(!SHAPE.AREA!)"
        arcpy.CalculateField_management(tempWS + "4P.shp", "Area", areaExpression,
                                        "PYTHON3")
        
        # Get the cell size of the input raster to remove small flecks
        with arcpy.da.UpdateCursor(tempWS + "4P.shp", "Area") as cursor:
            for row in cursor:
                if row[0] < cutoffArea:
                    cursor.deleteRow()
        
        # Process: Buffer (user defined buffer size, default 30m)
        
        arcpy.conversion.PolygonToRaster(tempWS + "4P.shp", 'gridcode', tempWS + "4Pras.tif", 'CELL_CENTER', 
                                         '', Cellsize_X)
        
        #transform buffer size in metres to buffer in pixel units as accepted by the "Expand" tool
        Buffer_s_p = float(Buffer_s) / Cellsize_X
        Buff_ras = Expand(tempWS + "4Pras.tif", Buffer_s_p, 1)
        Buff_ras.save(tempWS + "5B.tif")
        
        # Process: Clip
        # Clips the original DEM based on filtered DEM and keeping only the areas with higher seabed variations
        
        filtclipped = ExtractByMask(inputDEM, tempWS + "5B.tif")
        filtclipped.save(workspace + "/" + outRas)
        
        #Delete temporary files if requested
        if str(delTemp) == 'true':
            arcpy.AddMessage('Deleting temporary files ...')
            helper.cleaning(tempWS[:-1])
            
        return


class Fencing(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Fencing tool"
        self.description = "This script assists the Seabed-feature delineation creating a buffer of artificial pixels \
            around the boundaries of a raster. The buffer should prevent the misclassification of positive or negative features \
            that are cut by boundaries (i.e. nodata, which creates breaches in the Fill tool (see Seabed-feature delineation function below)"
            
        self.canRunInBackground = False
        self.category = "Pre-processing"
        
    def getParameterInfo(self):
        """Define parameter definitions"""
        parameters = None
        
        param0 = arcpy.Parameter(
            displayName="Input raster",
            name="inputDEM",
            datatype="GPRasterLayer",
            parameterType="Required",
            direction="Input")
        
        param1 = arcpy.Parameter(
            displayName="Are you mapping positive or negative features?",
            name="in_fill_direc",
            datatype="GPString",
            parameterType="Required",
            direction="Input")
        
        # Set a value
        param1.filter.type = "ValueList"
        param1.filter.list = ["Positive", "Negative"]
        
        param2 = arcpy.Parameter(
            displayName="Workspace",
            name="workspace",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")  
        
        param3 = arcpy.Parameter(
            displayName="Output fenced raster name",
            name="outRas",
            datatype="GPString",
            parameterType="Required",
            direction="Output")
        
        
        parameters = [param0, param1, param2, param3]
    
        return parameters
        
    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True
        
    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return
        
    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        
        if parameters[0].valueAsText is not None:
            folder = os.path.basename(os.path.dirname(parameters[0].valueAsText))
            #desc = arcpy.Describe(gdb)
            if folder.lower().endswith(('.gdb', '.mdb')):
                parameters[0].setErrorMessage("Geodatabases cannot be used in this version of the CoMMa Toolbox")
        
        if parameters[2].valueAsText is not None:
            folder = os.path.basename(parameters[2].valueAsText)
            #desc = arcpy.Describe(gdb)
            if folder.lower().endswith(('.gdb', '.mdb')):
                parameters[2].setErrorMessage("Geodatabases cannot be used in this version of the CoMMa Toolbox")
                
        return
        
    def execute(self, parameters, messages):
        """The source code of the tool."""
        
        arcpy.env.overwriteOutput = True
        inputDEM = parameters[0].valueAsText
        fillDirec = parameters[1].valueAsText
        workspace = parameters[2].valueAsText
        outRas = parameters[3].valueAsText
        
        if outRas[-4:] != '.tif':
            outRas = outRas + '.tif'
        # enable the helper functions
        helper = helpers()
        
        inputDEM = helper.convert_backslash_forwardslash(inputDEM)
        workspace = helper.convert_backslash_forwardslash(workspace)

        # if the input bathyRas is selected from a drop-down list, the bathyRas does not contain the full path
        # In this case, the full path needs to be obtained from the map layer
        if inputDEM.rfind("/") < 0:
            aprx = arcpy.mp.ArcGISProject("CURRENT")
            m = aprx.activeMap            
            for lyr in m.listLayers():
                if lyr.isRasterLayer:
                    if inputDEM == lyr.name: 
                        inputDEM = helper.convert_backslash_forwardslash(lyr.dataSource)
        
        Cellsize_Xresult = arcpy.GetRasterProperties_management(inputDEM, "CELLSIZEX")
        Cellsize_X = float(Cellsize_Xresult.getOutput(0))
        cutoffArea = (Cellsize_X*Cellsize_X)*3
        
        if fillDirec == "Positive":
            Value_Zresult = arcpy.GetRasterProperties_management(inputDEM, "MINIMUM")
            Value_Z = float(Value_Zresult.getOutput(0))
        
        else:
            Value_Zresult = arcpy.GetRasterProperties_management(inputDEM, "MAXIMUM")
            Value_Z = float(Value_Zresult.getOutput(0))
        
        tempWS = workspace + "\\temp_fencing\\"
        arcpy.env.overwriteOutput = True
        
        if not os.path.exists(tempWS):
            os.mkdir(tempWS)
            
        arcpy.env.qualifiedFieldNames = False
        arcpy.env.workspace = workspace
        
        #Calculation for buffer extent to be applied (twice the pixel dimension)
        Buffed = str((float(Cellsize_X))*4)+" Meters"
        
        #Calculation: simple integer raster
        IntegRas = Int(Raster(inputDEM)/Raster(inputDEM))
        IntegRas.save(tempWS + "int.tif")
                
        #Process: transform integer raster to polygon
        arcpy.RasterToPolygon_conversion(IntegRas, tempWS + "IntP.shp", "NO_SIMPLIFY", "")
        arcpy.management.RepairGeometry(tempWS + "IntP.shp")
        
        arcpy.AddMessage("Creating the buffer ...")        
        #Creates the buffer using a predetermined dimension
        arcpy.Buffer_analysis (tempWS + "IntP.shp", tempWS + "BuffP.shp", 
                               Buffed, "", "", "ALL")
                    
        # Add and calculate the value of the buffer pixels
        arcpy.AddField_management(tempWS + "BuffP.shp", "Wall_h", "DOUBLE", "8", "2")
                
        arcpy.CalculateField_management(tempWS + "BuffP.shp", "Wall_h", Value_Z, "PYTHON3")
                
        #Convert the polygon to a homogeneous raster with the same value (the wall height) assigned to its cells 
        bufraster = arcpy.PolygonToRaster_conversion(tempWS + "BuffP.shp", "Wall_h", tempWS + "_wall.tif", 
                                                     "CELL_CENTER", "", Cellsize_X)
        
        arcpy.AddMessage("Building the wall ...")         
        #Create a mosaic of the original raster and the buffer, so that the buffer is applied to the raster boundaries
        closed_ra = arcpy.MosaicToNewRaster_management([bufraster, inputDEM], workspace, outRas, "", 
                                                       "32_BIT_FLOAT", Cellsize_X, "1", "LAST", "LAST")
        
        return
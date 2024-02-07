# -*- coding: utf-8 -*-
#Authors: Riccardo Arosio and Joana Gafeira
#Institutions: University College Cork, British Geological Survey
#CoMMa version 1.0 - 2023

import arcpy
import os
import re
import math

from arcpy.sa import *


class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "Confined morphologies descriptor"
        self.alias = "CoMMa_Description"

        # List of tool classes associated with this toolbox
        self.tools = [Basic_descriptor, Texture_descriptor, Volume_descriptor]


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
    
    
    def BPI_calc(self,inputDEM,inner_radius,outer_radius):
        
        if int(inner_radius) > 0:
            neighborhood = NbrAnnulus(inner_radius, outer_radius, "CELL")
            
        else:
            neighborhood = NbrCircle(outer_radius, "CELL")
            
        out_focal_statistics = FocalStatistics(inputDEM, neighborhood, "MEAN")
        result_raster = Plus(Minus(inputDEM, out_focal_statistics), 0.5)
        
        return result_raster
    
    
    def LocDevfGlob_calc(self,inputDEM,analysisType,outer_radius):
                   
        zoneRaster = Int(Divide(inputDEM,inputDEM))
        neighborhood = NbrCircle(outer_radius, "CELL")
        
        if analysisType == "Mean":
            globalRas = ZonalStatistics(zoneRaster,"value",inputDEM,"MEAN")
            fcMean = FocalStatistics(inputDEM,neighborhood,"MEAN")
            outRaster = Minus(globalRas,fcMean)
   
        
        elif analysisType == "Median":

            globalRas = ZonalStatistics(zoneRaster,"value",inputDEM,"MEDIAN")
            fcMed = FocalStatistics(inputDEM,neighborhood,"MEDIAN")
            outRaster = Minus(globalRas,fcMed)
        
        return outRaster
    

    def Percentile_calc(self,inputFeatures):
        
        import numpy as np
        
        arr = arcpy.da.FeatureClassToNumPyArray(inputFeatures, 'LDfG_MEAN')
        arr = np.array(arr,np.float64)

        p1 = np.percentile(arr, 20)  # rank = 0
        p2 = np.percentile(arr, 40)  # rank = 1
        p3 = np.percentile(arr, 60)  # rank = 2
        p4 = np.percentile(arr, 80)  # rank = 3
        p5 = np.percentile(arr, 100)  # rank = 4
    
        #use cursor to update the new rank field
        with arcpy.da.UpdateCursor(inputFeatures, ['LDfG_MEAN', 'PeLDfG_Ran']) as cursor:
            for row in cursor:
                if row[0] < p1:
                    row[1] = 0  #rank 0
                elif p1 <= row[0] and row[0] < p2:
                    row[1] = 1
                elif p2 <= row[0] and row[0] < p3:
                    row[1] = 2
                elif p3 <= row[0] and row[0] < p4:
                    row[1] = 3
                else:
                    row[1] = 4
                    
                cursor.updateRow(row)    
        return


    def Ruggd_calc(self,neighborhood_size,out_aspect,out_slope):
        """
        Compute terrain ruggedness, using the vector ruggedness measure (VRM),
        as described in:
            Sappington et al., 2007. Quantifying Landscape Ruggedness for
            Animal Habitat Analysis: A Case Study Using Bighorn Sheep in the
            Mojave Desert. Journal of Wildlife Management. 71(5): 1419 -1426.
        """
  
        hood_size = int(neighborhood_size)
    
        try:
            # Convert Slope and Aspect rasters to radians
            #print("Converting slope and aspect to radians...")
            slope_rad = out_slope * (math.pi / 180)
            aspect_rad = out_aspect * (math.pi / 180)
    
            # Calculate x, y, and z rasters
            #print("Calculating x, y, and z rasters...")
            xy_raster_calc = Sin(slope_rad)
            z_raster_calc = Cos(slope_rad)
            x_raster_calc = Con(out_aspect == -1, 0, Sin(aspect_rad)) * xy_raster_calc
            y_raster_calc = Con(out_aspect == -1, 0, Cos(aspect_rad)) * xy_raster_calc
    
            # Calculate sums of x, y, and z rasters for selected neighborhood size
            #print("Calculating sums of x, y, and z rasters in neighborhood...")
            hood = NbrRectangle(hood_size, hood_size, "CELL")
            x_sum_calc = FocalStatistics(x_raster_calc, hood, "SUM", "NODATA")
            y_sum_calc = FocalStatistics(y_raster_calc, hood, "SUM", "NODATA")
            z_sum_calc = FocalStatistics(z_raster_calc, hood, "SUM", "NODATA")
    
            # Calculate the resultant vector
            #print("Calculating the resultant vector...")
            result_vect = (x_sum_calc**2 + y_sum_calc**2 + z_sum_calc**2)**0.5
    
            arcpy.env.rasterStatistics = "STATISTICS"
            # Calculate the Ruggedness raster
            #print("Calculating the final ruggedness raster...")
            ruggedness = 1 - (result_vect / hood_size**2)
        
        except Exception as e:
            print(e)    
        
        return ruggedness
  

    def AVI(self,out_slope,out_aspect):
        
        """
        AVI (Aspect Variability Index) is a measure of topographic complexity 
        used in Neilson et al. that measures variability in aspects.
        Nielsen SE, Herrero S, Boyce MS, Mace RD, Benn B, Gibeau ML, and 
        Jevons S: Modelling the spatial distribution of human-caused grizzly 
        bear mortalities in the Central Rockies ecosystem of Canada. 
        Biological Conservation 2004, 120:101-113.
        """
        
        arcpy.env.rasterStatistics = "STATISTICS"
        
        try:
            neighborhood = NbrRectangle(3, 3, "CELL")
                
            out_Var_aspect = FocalStatistics(out_aspect, neighborhood, "VAR")
            out_Mean_slope = FocalStatistics(out_slope, neighborhood, "MEAN")
            AVI_raster = (out_Var_aspect * out_Mean_slope)/(out_Var_aspect + out_Mean_slope)
    
        except Exception as e:
            print(e)
        
        return AVI_raster

    
    def field_removal(self,inputFeatures,keepFields):
        
        allFields = [field.name for field in arcpy.ListFields(inputFeatures)]
        
        for f in allFields:
            if f in keepFields:
                pass
            else:
                try:
                    arcpy.DeleteField_management(inputFeatures, f)
                except arcpy.ExecuteError:
                    pass
        
        return inputFeatures
    
    
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
    
    
    def makeNullZero(self,i):
            
        if i is None:
            return 0
        else:
            return i


class Basic_descriptor(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Basic descriptors"
        self.description = "Basic and geometrical description tool"
        self.canRunInBackground = False

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
            displayName="Input delineated morphologies",
            name="in_features",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")
        
        param2 = arcpy.Parameter(
            displayName="Workspace",
            name="workspace",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input") 

        param3 = arcpy.Parameter(
            displayName="Output Feature name",
            name="outFeat",
            datatype="GPString",
            parameterType="Required",
            direction="Output")
        
        param4 = arcpy.Parameter(
            displayName="Are you mapping positive or negative features?",
            name="in_fill_direc",
            datatype="GPString",
            parameterType="Required",
            direction="Input")
        
        # Set a value
        param4.filter.type = "ValueList"
        param4.filter.list = ["Positive", "Negative"]

        param5 = arcpy.Parameter(
            displayName="Geomorphons raster file",
            name="inputGeo",
            datatype="GPRasterLayer",
            parameterType="Optional",
            direction="Input")
        
        param6 = arcpy.Parameter(
            displayName="Do you want to delete the temporary files?",
            name="delTemp",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")
        
        parameters = [param0, param1, param2, param3, param4, param5, param6]
        
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

        if parameters[1].valueAsText is not None:
            folder = os.path.basename(os.path.dirname(parameters[1].valueAsText))
            #desc = arcpy.Describe(gdb)
            if folder.lower().endswith(('.gdb', '.mdb')):
                parameters[1].setErrorMessage("Geodatabases cannot be used in this version of the CoMMa Toolbox")
        
        if parameters[2].valueAsText is not None:
            folder = os.path.basename(parameters[2].valueAsText)
            #desc = arcpy.Describe(gdb)
            if folder.lower().endswith(('.gdb', '.mdb')):
                parameters[2].setErrorMessage("Geodatabases cannot be used in this version of the CoMMa Toolbox")
        
        if parameters[5].valueAsText is not None:
            folder = os.path.basename(os.path.dirname(parameters[5].valueAsText))
            #desc = arcpy.Describe(gdb)
            if folder.lower().endswith(('.gdb', '.mdb')):
                parameters[5].setErrorMessage("Geodatabases cannot be used in this version of the CoMMa Toolbox")
                
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        
        arcpy.env.overwriteOutput = True
        inputDEM = parameters[0].valueAsText
        in_features = parameters[1].valueAsText
        workspace = parameters[2].valueAsText
        outFeat = parameters[3].valueAsText
        fillDirec = parameters[4].valueAsText
        inputGeo = parameters[5].valueAsText
        delTemp = parameters[6].valueAsText
        
        # enable the helper functions
        helper = helpers()
        
        inputDEM = helper.convert_backslash_forwardslash(inputDEM)
        in_features = helper.convert_backslash_forwardslash(in_features)
        workspace = helper.convert_backslash_forwardslash(workspace)
        inputGeo = helper.convert_backslash_forwardslash(inputGeo)
        
        if outFeat[-4:] != '.shp':
            outFeat = outFeat + '.shp'
            
        if inputDEM.rfind("/") < 0:
            aprx = arcpy.mp.ArcGISProject("CURRENT")
            m = aprx.activeMap            
            for lyr in m.listLayers():
                if lyr.isRasterLayer:
                    if inputDEM == lyr.name: 
                        inputDEM = helper.convert_backslash_forwardslash(lyr.dataSource)
                        
        if in_features.rfind("/") < 0:
            aprx = arcpy.mp.ArcGISProject("CURRENT")
            m = aprx.activeMap            
            for lyr in m.listLayers():
                if lyr.isFeatureLayer:
                    if in_features == lyr.name: 
                        in_features = helper.convert_backslash_forwardslash(lyr.dataSource)
        
        if parameters[5].value is not None:
            if inputGeo.rfind("/") < 0:
                aprx = arcpy.mp.ArcGISProject("CURRENT")
                m = aprx.activeMap            
                for lyr in m.listLayers():
                    if lyr.isRasterLayer:
                        if inputGeo == lyr.name: 
                            inputGeo = helper.convert_backslash_forwardslash(lyr.dataSource)
        
        
        Cellsize_Xresult = arcpy.GetRasterProperties_management(inputDEM, "CELLSIZEX")
        Cellsize_X2 = (float(Cellsize_Xresult.getOutput(0)))**2
        
        tempWS = workspace + "\\temp\\"
        arcpy.env.overwriteOutput = True
        if not os.path.exists(tempWS):
            os.mkdir(tempWS)
        
        arcpy.Copy_management(in_features, tempWS + "deli.shp")
              
        arcpy.env.qualifiedFieldNames = False
        arcpy.env.workspace = workspace
        
        #list of fields that will be saved - check if fields external to the CoMMa toolbox are present and saves them
        keepFields = []
        CM_delim_fields = ["Area","Perimeter","MBG_Width","MBG_Length","MBG_W_L","MBG_Orient"]
        fie = [field.name for field in arcpy.ListFields(tempWS + "deli.shp")]
        for f in fie:
            if f not in CM_delim_fields:
                keepFields.append(f)
        
        arcpy.AddMessage("Calculating derivatives ...")
        
        #Vertical relief
        arcpy.AddMessage("... Vertical relief")
        if fillDirec == "Positive":
            featuresDEM = Negate(inputDEM)  
            clip_rel = ExtractByMask(featuresDEM, tempWS + "deli.shp")
            filledDEM = Fill(clip_rel)
            Vert_rel = Minus(filledDEM, clip_rel)
        else:
            clip_rel = ExtractByMask(inputDEM, tempWS + "deli.shp")
            filledDEM = Fill(clip_rel)
            Vert_rel = Minus(filledDEM, clip_rel)
        
        # Slope
        arcpy.AddMessage("... Slope")
        Slope_ = Slope(inputDEM, "DEGREE", 1)

        #Local deviation from global median
        arcpy.AddMessage("... Deviation from Global Median")
        LDfG_ = helper.LocDevfGlob_calc(inputDEM,"Median",3)
        
        List_dev = [Vert_rel, Slope_, LDfG_, inputDEM]
        List_nam = ["Vert_rel", "Slope_", "LDfG_", "Depth"]
        c=0
        for ra in List_dev:
            ra_name = List_nam[c]
            ZonalStatisticsAsTable(tempWS + "deli.shp", "FID", ra, tempWS + ra_name, 
                                   "DATA", "MIN_MAX_MEAN")
            c+=1
        
        arcpy.MakeFeatureLayer_management(tempWS + "deli.shp", tempWS + "deli_lyr")
        arcpy.AddJoin_management(tempWS + "deli_lyr", "FID", 
                                 tempWS + "Vert_rel", "FID")
        arcpy.AddJoin_management(tempWS + "deli_lyr", "deli.FID", 
                                 tempWS + "Slope_", "FID")
        arcpy.AddJoin_management(tempWS + "deli_lyr", "deli.FID", 
                                 tempWS + "LDfG_", "FID")
        arcpy.AddJoin_management(tempWS + "deli_lyr", "deli.FID", 
                                 tempWS + "Depth", "FID")
        
        arcpy.AddMessage("Changing names to original fields ...")
        
        arcpy.CreateFileGDB_management(tempWS[:-1], "field_alter.gdb")
        arcpy.FeatureClassToGeodatabase_conversion(tempWS + "deli_lyr",  tempWS + "field_alter.gdb")
        base_lyr = tempWS + "field_alter.gdb/deli_lyr"
        
        arcpy.AlterField_management(base_lyr, "MIN", "VeRe_MIN", "VeRe_MIN")
        arcpy.AlterField_management(base_lyr, "MAX", "Conf_VR", "Conf_VR")
        arcpy.AlterField_management(base_lyr, "MEAN", "VeRe_MEAN", "VeRe_MEAN")
        arcpy.AlterField_management(base_lyr, "MIN_1", "Slope_MIN", "Slope_MIN")
        arcpy.AlterField_management(base_lyr, "MAX_1", "Slope_MAX", "Slope_MAX")
        arcpy.AlterField_management(base_lyr, "MEAN_1", "Slope_MEAN", "Slope_MEAN")
        arcpy.AlterField_management(base_lyr, "MIN_12", "LDfG_MIN", "LDfG_MIN")
        arcpy.AlterField_management(base_lyr, "MAX_12", "LDfG_MAX", "LDfG_MAX")
        arcpy.AlterField_management(base_lyr, "MEAN_12", "LDfG_MEAN", "LDfG_MEAN")
        arcpy.AlterField_management(base_lyr, "MIN_12_13", "Depth_MIN", "Depth_MIN")
        arcpy.AlterField_management(base_lyr, "MAX_12_13", "Depth_MAX", "Depth_MAX")
        arcpy.AlterField_management(base_lyr, "MEAN_12_13", "Depth_MEAN", "Depth_MEAN")
        
        arcpy.conversion.FeatureClassToShapefile(base_lyr, tempWS[:-1])
        arcpy.Copy_management(tempWS + "deli_lyr.shp", tempWS + "sh.shp")
        
        #dropFields = ["MBG_Width", "MBG_Length", "MBG_Orient", "MBG_W_L"]
        #arcpy.DeleteField_management(tempWS + "sh.shp", CM_delim_fields)
        
        
        if parameters[5].value is not None:
            
            arcpy.AddMessage("Processing geomorphons ...")
            
            keepFields.extend(["Conf_VR","Slope_MIN","Slope_MAX","Slope_MEAN",
                               "LDfG_MAX","LDfG_MEAN","Depth_MIN","Depth_MAX",
                               "Depth_MEAN"])
            helper.field_removal(tempWS + "sh.shp", keepFields)
            
            #inputGeo_MjFi = MajorityFilter(inputGeo, "EIGHT", "HALF")
            
            if fillDirec == "Positive":
                remapPEAK = "1 NODATA;2 1; 3 10 NODATA"
                remapRID = "1 2 NODATA;3 1; 4 10 NODATA"
                remapSHOU = "1 3 NODATA;4 1; 5 10 NODATA"
                remapSPUR = "1 4 NODATA;5 1; 6 10 NODATA"
                remapSLO = "1 5 NODATA;6 1; 7 10 NODATA"
            
                OutRecla_Pea = Reclassify(inputGeo, "Value", remapPEAK, "NODATA")
                OutRecla_Rid = Reclassify(inputGeo, "Value", remapRID, "NODATA")
                OutRecla_Shou = Reclassify(inputGeo, "Value", remapSHOU, "NODATA")
                OutRecla_Spur = Reclassify(inputGeo, "Value", remapSPUR, "NODATA")
                OutRecla_Slo = Reclassify(inputGeo, "Value", remapSLO, "NODATA")
            
                arcpy.AddField_management(tempWS + "sh.shp", "gm_peak", "DOUBLE", "8", "2")
                arcpy.AddField_management(tempWS + "sh.shp", "gm_ridg", "DOUBLE", "8", "2")
                arcpy.AddField_management(tempWS + "sh.shp", "gm_shou", "DOUBLE", "8", "2")
                arcpy.AddField_management(tempWS + "sh.shp", "gm_spur", "DOUBLE", "8", "2")
                arcpy.AddField_management(tempWS + "sh.shp", "gm_slop", "DOUBLE", "8", "2")
            
                #extract the sum of each polygon feature
                ZonalStatisticsAsTable(tempWS + "sh.shp", "FID", 
                                       OutRecla_Pea, tempWS + "peak", 
                                       "DATA", "SUM")
                ZonalStatisticsAsTable(tempWS + "sh.shp", "FID", 
                                       OutRecla_Rid, tempWS + "ridg", 
                                       "DATA", "SUM")
                ZonalStatisticsAsTable(tempWS + "sh.shp", "FID", 
                                       OutRecla_Shou, tempWS + "shou", 
                                       "DATA", "SUM")
                ZonalStatisticsAsTable(tempWS + "sh.shp", "FID", 
                                       OutRecla_Spur, tempWS + "spur", 
                                       "DATA", "SUM")
                ZonalStatisticsAsTable(tempWS + "sh.shp", "FID", 
                                       OutRecla_Slo, tempWS + "slop", 
                                       "DATA", "SUM")
            
                if not arcpy.Exists(tempWS + "sh_lyr"):
                    arcpy.MakeFeatureLayer_management(tempWS + "sh.shp", tempWS + "sh_lyr")
            
                arcpy.AddJoin_management(tempWS + "sh_lyr", "FID", 
                                         tempWS + "peak", "FID", "KEEP_ALL")
                arcpy.AddJoin_management(tempWS + "sh_lyr", "sh.FID", 
                                         tempWS + "ridg", "FID", "KEEP_ALL")
                arcpy.AddJoin_management(tempWS + "sh_lyr", "sh.FID", 
                                         tempWS + "shou", "FID", "KEEP_ALL")
                arcpy.AddJoin_management(tempWS + "sh_lyr", "sh.FID", 
                                         tempWS + "spur", "FID", "KEEP_ALL")
                arcpy.AddJoin_management(tempWS + "sh_lyr", "sh.FID", 
                                         tempWS + "slop", "FID", "KEEP_ALL")
                
                arcpy.CreateFileGDB_management(tempWS[:-1], "field_alter_geom.gdb")
                arcpy.FeatureClassToGeodatabase_conversion(tempWS + "sh_lyr",  tempWS + "field_alter_geom.gdb")
                base_lyr = tempWS + "field_alter_geom.gdb/sh_lyr"
                
                arcpy.AlterField_management(base_lyr, "SUM", "sum_pea", "sum_pea")
                arcpy.AlterField_management(base_lyr, "SUM_1", "sum_rid", "sum_rid")
                arcpy.AlterField_management(base_lyr, "SUM_12", "sum_sho", "sum_sho")
                arcpy.AlterField_management(base_lyr, "SUM_12_13", "sum_spu", "sum_spu")
                arcpy.AlterField_management(base_lyr, "SUM_12_13_14", "sum_slo", "sum_slo")
                
            
                relExpression = "(float(!sum_pea!)*{} / float(!SHAPE.AREA!))*100".format(Cellsize_X2)
                arcpy.CalculateField_management(base_lyr, "sum_pea", "'' if !sum_pea! is None else !sum_pea!", 
                                                "PYTHON3")
                arcpy.CalculateField_management(base_lyr, "gm_peak", relExpression, "PYTHON3")
                relExpression = "(float(!sum_rid!)*{} / float(!SHAPE.AREA!))*100".format(Cellsize_X2)
                arcpy.CalculateField_management(base_lyr, "sum_rid", "'' if !sum_rid! is None else !sum_rid!", 
                                                "PYTHON3")
                arcpy.CalculateField_management(base_lyr, "gm_ridg", relExpression, "PYTHON3")
                relExpression = "(float(!sum_sho!)*{} / float(!SHAPE.AREA!))*100".format(Cellsize_X2)
                arcpy.CalculateField_management(base_lyr, "sum_sho", "'' if !sum_sho! is None else !sum_sho!", 
                                                "PYTHON3")
                arcpy.CalculateField_management(base_lyr, "gm_shou", relExpression, "PYTHON3")
                relExpression = "(float(!sum_spu!)*{} / float(!SHAPE.AREA!))*100".format(Cellsize_X2)
                arcpy.CalculateField_management(base_lyr, "sum_spu", "'' if !sum_spu! is None else !sum_spu!", 
                                "PYTHON3")
                arcpy.CalculateField_management(base_lyr, "gm_spur", relExpression, "PYTHON3")
                relExpression = "(float(!sum_slo!)*{} / float(!SHAPE.AREA!))*100".format(Cellsize_X2)
                arcpy.CalculateField_management(base_lyr, "sum_slo", "'' if !sum_slo! is None else !sum_slo!", 
                                "PYTHON3")
                arcpy.CalculateField_management(base_lyr, "gm_slop", relExpression, "PYTHON3")
            
                arcpy.conversion.FeatureClassToShapefile(base_lyr, tempWS[:-1])
                arcpy.Copy_management(tempWS + "sh_lyr.shp", tempWS + "shb.shp")
                
                
                #Count the number of PEAKS in each feature
                arcpy.AddMessage("Counting peaks ...")
                arcpy.RasterToPolygon_conversion(OutRecla_Pea, tempWS + "peaks.shp", "NO_SIMPLIFY")
                arcpy.AddField_management(tempWS + "peaks.shp", "Area", "DOUBLE", "8", "1")
                areaExpression = "float(!SHAPE.AREA!)"
                arcpy.CalculateField_management(tempWS + "peaks.shp", "Area", areaExpression, "PYTHON3")
                
                if not arcpy.Exists(tempWS + "peaks_lyr"):
                    arcpy.MakeFeatureLayer_management(tempWS + "peaks.shp", tempWS + "peaks_lyr")
                
                Areatresh = str(Cellsize_X2*4)
                arcpy.SelectLayerByAttribute_management(tempWS + "peaks_lyr", "NEW_SELECTION", '("Area" < ' + Areatresh + ')')
                
                #Delete polygons with an area inferior to MinArea or with an elongation inferior to MinRatio
                arcpy.DeleteRows_management(tempWS + "peaks_lyr")
                  
                arcpy.SpatialJoin_analysis(tempWS + "shb.shp", tempWS + "peaks.shp", tempWS + "shbp.shp",
                                           "JOIN_ONE_TO_ONE","KEEP_ALL", "","COMPLETELY_CONTAINS")
                
                arcpy.AddField_management(tempWS + "shbp.shp", "peak_no", "DOUBLE", "8", "2")
                 
                relExpression = "int(!Join_Count!)"
                arcpy.CalculateField_management(tempWS + "shbp.shp", "peak_no", relExpression, "PYTHON3")
                
                arcpy.Delete_management(tempWS + "peaks_lyr")
                
                poly1 = tempWS + "shbp.shp"
            
            else:
                remapHOLL = "1 6 NODATA;7 1; 8 10 NODATA"
                remapVALL = "1 8 NODATA;9 1; 10 NODATA"
                remapPIT = "1 9 NODATA;10 1"
            
                OutRecla_Hol = Reclassify(inputGeo, "Value", remapHOLL, "NODATA")
                OutRecla_Val = Reclassify(inputGeo, "Value", remapVALL, "NODATA")
                OutRecla_Pit = Reclassify(inputGeo, "Value", remapPIT, "NODATA")
            
                arcpy.AddField_management(tempWS + "sh.shp", "gm_holl", "DOUBLE", "8", "2")
                arcpy.AddField_management(tempWS + "sh.shp", "gm_vall", "DOUBLE", "8", "2")
                arcpy.AddField_management(tempWS + "sh.shp", "gm_pit", "DOUBLE", "8", "2")
            
                #extract the sum of each polygon feature
                ZonalStatisticsAsTable(tempWS + "sh.shp", "FID", 
                                       OutRecla_Hol, tempWS + "hol", 
                                       "DATA", "SUM")
                ZonalStatisticsAsTable(tempWS + "sh.shp", "FID", 
                                       OutRecla_Val, tempWS + "val", 
                                       "DATA", "SUM")
                ZonalStatisticsAsTable(tempWS + "sh.shp", "FID", 
                                       OutRecla_Pit, tempWS + "pi", 
                                       "DATA", "SUM")
                
                if not arcpy.Exists(tempWS + "sh_lyr"):
                    arcpy.MakeFeatureLayer_management(tempWS + "sh.shp", tempWS + "sh_lyr")
            
                arcpy.AddJoin_management(tempWS + "sh_lyr", "FID", 
                                         tempWS + "hol", "FID", "KEEP_ALL")
                arcpy.AddJoin_management(tempWS + "sh_lyr", "sh.FID", 
                                         tempWS + "val", "FID", "KEEP_ALL")
                arcpy.AddJoin_management(tempWS + "sh_lyr", "sh.FID", 
                                         tempWS + "pi", "FID", "KEEP_ALL")
                
                arcpy.CreateFileGDB_management(tempWS[:-1], "field_alter_geom.gdb")
                arcpy.FeatureClassToGeodatabase_conversion(tempWS + "sh_lyr",  tempWS + "field_alter_geom.gdb")
                base_lyr = tempWS + "field_alter_geom.gdb/sh_lyr"
                
                arcpy.AlterField_management(base_lyr, "SUM", "sum_hol", "sum_hol")
                arcpy.AlterField_management(base_lyr, "SUM_1", "sum_val", "sum_val")
                arcpy.AlterField_management(base_lyr, "SUM_12", "sum_pit", "sum_pit")
            
                relExpression = "(float(!sum_hol!)*{} / float(!SHAPE.AREA!))*100".format(Cellsize_X2)
                arcpy.CalculateField_management(base_lyr, "sum_hol", "'' if !sum_hol! is None else !sum_hol!", 
                                                "PYTHON3")
                arcpy.CalculateField_management(base_lyr, "gm_holl", relExpression, "PYTHON3")
                relExpression = "(float(!sum_val!)*{} / float(!SHAPE.AREA!))*100".format(Cellsize_X2)
                arcpy.CalculateField_management(base_lyr, "sum_val", "'' if !sum_val! is None else !sum_val!", 
                                                "PYTHON3")
                arcpy.CalculateField_management(base_lyr, "gm_vall", relExpression, "PYTHON3")
                relExpression = "(float(!sum_pit!)*{} / float(!SHAPE.AREA!))*100".format(Cellsize_X2)
                arcpy.CalculateField_management(base_lyr, "sum_pit", "'' if !sum_pit! is None else !sum_pit!", 
                                                "PYTHON3")
                arcpy.CalculateField_management(base_lyr, "gm_pit", relExpression, "PYTHON3")
                
                arcpy.conversion.FeatureClassToShapefile(base_lyr, tempWS[:-1])
                arcpy.Copy_management(tempWS + "sh_lyr.shp", tempWS + "shb.shp")

                #Count the number of PITS in each feature
                arcpy.RasterToPolygon_conversion(OutRecla_Pit, tempWS + "pits.shp", "NO_SIMPLIFY")
                arcpy.AddField_management(tempWS + "pits.shp", "Area", "DOUBLE", "8", "1")
                areaExpression = "float(!SHAPE.AREA!)"
                arcpy.CalculateField_management(tempWS + "pits.shp", "Area", areaExpression, "PYTHON3")
                
                if not arcpy.Exists(tempWS + "pits_lyr"):
                    arcpy.MakeFeatureLayer_management(tempWS + "pits.shp", tempWS + "pits_lyr")
                
                Areatresh = str(Cellsize_X2*4)
                arcpy.AddMessage("selection ...")
                arcpy.SelectLayerByAttribute_management(tempWS + "pits_lyr", "NEW_SELECTION", '("Area" < ' + Areatresh + ')')
                
                #Delete polygons with an area inferior to MinArea or with an elongation inferior to MinRatio
                arcpy.DeleteRows_management(tempWS + "pits_lyr")
                arcpy.Delete_management(tempWS + "pits_lyr")
                    
                arcpy.SpatialJoin_analysis(tempWS + "shb.shp", tempWS + "pits.shp", tempWS + "shbp.shp",
                                           "JOIN_ONE_TO_ONE","KEEP_ALL", "","COMPLETELY_CONTAINS")
                
                arcpy.AddField_management(tempWS + "shbp.shp", "pit_no", "DOUBLE", "8", "2")
                relExpression = "int(!Join_Count!)"
                arcpy.CalculateField_management(tempWS + "shbp.shp", "pit_no", relExpression, "PYTHON3")
                
                poly1 = tempWS + "shbp.shp"
                
        else:
            poly1 = tempWS + "sh.shp"
        
        keepFields.extend(["MBG_Width","MBG_Length","MBG_W_L","MBG_Orient",
                           "Geo_rat","Conf_VR","Slope_MIN","Slope_MAX",
                           "Slope_MEAN","LDfG_MAX","LDfG_MEAN","Depth_MIN","Depth_MAX",
                           "Depth_MEAN","gm_peak","gm_ridg","gm_shou","gm_spur","gm_slop",
                           "gm_holl","gm_vall","gm_pit","peak_no","pit_no"])
        helper.field_removal(poly1, keepFields)
        
          
        arcpy.AddMessage("Calculating area, perimeter and W/L ratio and orientation ...")
        arcpy.DeleteField_management(poly1, "AREA")
        arcpy.AddField_management(poly1, "Area", "DOUBLE", "8", "1")
        areaExpression = "float(!SHAPE.AREA!)"
        arcpy.CalculateField_management(poly1, "Area", areaExpression, "PYTHON3")
    
        arcpy.DeleteField_management(poly1, "PERIMETER")
        arcpy.AddField_management(poly1, "Perimeter", "DOUBLE", "8", "1") 
        perimeterExpression = "float(!SHAPE.LENGTH!)"
        arcpy.CalculateField_management(poly1, "Perimeter", perimeterExpression, "PYTHON3")
        
        
        # Calculate the MBG for each polygon and adding the following field: MBG_Width, MBG_Length, MBG_Orientation
        arcpy.MinimumBoundingGeometry_management(poly1, tempWS + "mbg.shp",
                                                 "RECTANGLE_BY_WIDTH", "NONE", "", "MBG_FIELDS")
        # Add and calculate the MBG Width/Length field
        #arcpy.AddField_management(tempWS + "mbg.shp", "MBG_W_L", "DOUBLE", "8", "2")
        arcpy.CalculateField_management(tempWS + "mbg.shp", "MBG_W_L", '!MBG_Width!/!MBG_Length!', "PYTHON3")
          
        # Join the MBG fields to the feature outline polygon shapefile
        arcpy.SpatialJoin_analysis(poly1, tempWS + "mbg.shp", tempWS + "joi.shp",
                                   "JOIN_ONE_TO_MANY","KEEP_ALL", "","WITHIN")
        
        # Delete polygons that got the attributes from another MBG
        if not arcpy.Exists(tempWS + "joi_lyr"):
            arcpy.MakeFeatureLayer_management(tempWS + "joi.shp", tempWS + "joi_lyr")
            
        arcpy.SelectLayerByAttribute_management(tempWS + "joi_lyr", "NEW_SELECTION", '("TARGET_FID" <> "JOIN_FID")')
        arcpy.DeleteRows_management(tempWS + "joi_lyr")
        
        arcpy.AddMessage("Calculating Convex Hull area ...")
        # Add and calculate the MBG convex hull field
        arcpy.MinimumBoundingGeometry_management(tempWS + "joi.shp", tempWS + "ch.shp", "CONVEX_HULL", "NONE")
        
        arcpy.AddField_management(tempWS + "ch.shp", "Area_CH", "DOUBLE", "8", "2")
        areaExpression = "float(!SHAPE.AREA!)"
        arcpy.CalculateField_management(tempWS + "ch.shp", "Area_CH", areaExpression, "PYTHON3")
        
        arcpy.SpatialJoin_analysis(tempWS + "joi.shp", tempWS + "ch.shp", tempWS + "chm.shp", "JOIN_ONE_TO_MANY","KEEP_ALL", "","WITHIN")    

        # Delete polygons that got the attributes from another MBG 
        if not arcpy.Exists(workspace + "chm_lyr"):
            arcpy.MakeFeatureLayer_management(tempWS + "chm.shp", tempWS + "chm_lyr")
            
        arcpy.SelectLayerByAttribute_management(tempWS + "chm_lyr", "NEW_SELECTION", '("TARGET_FID" <> "JOIN_FID")')
        arcpy.DeleteRows_management(tempWS + "chm_lyr")
        arcpy.CopyFeatures_management(tempWS + "chm.shp", tempWS + "b.shp")

        
        arcpy.AddMessage("Calculating Dissection, Depth range, CH and PP scores, BPI variance ...")
        
        #area of footprint's convex hull, and ranges between 0 (lack of boundary concavities) and ~1.
        arcpy.AddField_management(tempWS + "b.shp", "CH_Score", "DOUBLE", "8", "2")
        arcpy.CalculateField_management(tempWS + "b.shp", "CH_Score", '!SHAPE.AREA!/!Area_CH!', "PYTHON3") 
        
        # Create a field to store Prof_Ind values
        arcpy.AddField_management(tempWS + "b.shp", "Dissect", "DOUBLE", "4", "2")
        arcpy.CalculateField_management(tempWS + "b.shp", "Dissect", '(!Depth_MIN! - !Depth_MEAN!) / ((!Depth_MIN! - !Depth_MAX!) + 0.1)', "PYTHON3")
        
        arcpy.AddField_management(tempWS + "b.shp", "Depth_RAN", "DOUBLE", "4", "2")
        arcpy.CalculateField_management(tempWS + "b.shp", "Depth_RAN", '(abs(!Depth_MIN!) - abs(!Depth_MAX!))', "PYTHON3")
        
        # Calculate the compactness of an object using the Polsby-Popper (PP) score. 
        arcpy.AddField_management(tempWS + "b.shp", "PP_Score", "DOUBLE", "8", "2") 
        arcpy.CalculateField_management(tempWS + "b.shp", "PP_Score", '4*3.1492*!Area!/(!Perimeter!*!Perimeter!)', "PYTHON3")
        
        #Calculate the variance from LDfG mean for each feature
        arcpy.AddField_management(tempWS + "b.shp", "LDfG_VAR", "DOUBLE", "8", "2") 
        
        arcpy.analysis.Statistics(tempWS + "b.shp", tempWS + "std_LDfG", [["LDfG_MEAN", "MEAN"]])
        SC = arcpy.SearchCursor(tempWS + "std_LDfG")
        for row in SC:
             mean = row.getValue("MEAN_LDfG_MEAN")
        
        arcpy.CalculateField_management(tempWS + "b.shp", "LDfG_VAR", '((!LDfG_MEAN! - %d)**2)'%mean, "PYTHON3")
        
        #Calculate the mean BPI percentiles
        arcpy.AddField_management(tempWS + "b.shp", "PeLDfG_Ran", "LONG")
        helper.Percentile_calc(tempWS + "b.shp")
    
        #Save file and final cleaning
        
        #Delete useless fields
        arcpy.AddMessage('Final cleaning ...') 
        keepFields.extend(["Perimeter","Area","Area_CH","CH_Score","Dissect","Depth_RAN",
                           "PP_Score","LDfG_VAR","PeLDfG_Ran"])
        helper.field_removal(tempWS + "b.shp", keepFields)
        
        arcpy.AddMessage('Saving official shapefile ...')
        
        arcpy.Copy_management(tempWS + "/" + "b.shp", workspace + "/" + outFeat)
        #arcpy.DeleteField_management(workspace + "/" + outFeat, "Shape_Length;Shape_Area;InPoly_F_1;Shape_Leng;gm_scor;gm_pnum")
        
        #Delete temporary files if requested
        if str(delTemp) == 'true':
            arcpy.AddMessage('Deleting temporary files ...')
            helper.cleaning(tempWS[:-1])
            
        #Print file report
        txtFile = open(workspace + "/" + str(outFeat[:-4]) + "_Info_Basic_Desc.txt", "w")
        txtFile.write("Script: CoMMa Delineation ToolBox v1.0" "\n")
        txtFile.write("\n")
        txtFile.write("File Name: " + str(outFeat) + "\n")
        txtFile.write("Input DEM: " + os.path.basename(inputDEM) + "\n")
        txtFile.write("Input features: " + os.path.basename(in_features) + "\n")
                    
        txtFile.write("\n")
        txtFile.write("Relief type: " + fillDirec + "\n")
        txtFile.write("\n")
        #txtFile.write("BPI inner radius: " + inner_radius + "\n")
        #txtFile.write("BPI outer radius: " + outer_radius + "\n")
        
        txtFile.write("\n")
        if parameters[5].value is not None:
            txtFile.write("Input geomorphons file: " + os.path.basename(inputGeo) + "\n")
        else:
            txtFile.write("The geomorphons add-on was not implemented")

        txtFile.close()
        
        
        return


class Texture_descriptor(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Texture descriptors"
        self.description = "Texture description tool"
        self.canRunInBackground = False

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
            displayName="Input backscatter raster (if available)",
            name="in_backsc",
            datatype="GPRasterLayer",
            parameterType="Optional",
            direction="Input")
        
        param2 = arcpy.Parameter(
            displayName="Input delineated morphologies",
            name="in_features",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")
        
        param3 = arcpy.Parameter(
            displayName="Workspace",
            name="workspace",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input") 

        param4 = arcpy.Parameter(
            displayName="Output Feature name",
            name="outFeat",
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

        if parameters[1].valueAsText is not None:
            folder = os.path.basename(os.path.dirname(parameters[1].valueAsText))
            #desc = arcpy.Describe(gdb)
            if folder.lower().endswith(('.gdb', '.mdb')):
                parameters[1].setErrorMessage("Geodatabases cannot be used in this version of the CoMMa Toolbox")
        
        if parameters[2].valueAsText is not None:
            folder = os.path.basename(os.path.dirname(parameters[2].valueAsText))
            #desc = arcpy.Describe(gdb)
            if folder.lower().endswith(('.gdb', '.mdb')):
                parameters[2].setErrorMessage("Geodatabases cannot be used in this version of the CoMMa Toolbox")
        
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
        in_backsc = parameters[1].valueAsText
        in_features = parameters[2].valueAsText
        workspace = parameters[3].valueAsText
        outFeat = parameters[4].valueAsText
        delTemp = parameters[5].valueAsText
        
        # enable the helper functions
        helper = helpers()
        
        inputDEM = helper.convert_backslash_forwardslash(inputDEM)
        in_backsc = helper.convert_backslash_forwardslash(in_backsc)
        in_features = helper.convert_backslash_forwardslash(in_features)
        workspace = helper.convert_backslash_forwardslash(workspace)
        
        if outFeat[-4:] != '.shp':
            outFeat = outFeat + '.shp'
            
        if inputDEM.rfind("/") < 0:
            aprx = arcpy.mp.ArcGISProject("CURRENT")
            m = aprx.activeMap            
            for lyr in m.listLayers():
                if lyr.isRasterLayer:
                    if inputDEM == lyr.name: 
                        inputDEM = helper.convert_backslash_forwardslash(lyr.dataSource)
        
        if parameters[1].value is not None:
            if in_backsc.rfind("/") < 0:
                aprx = arcpy.mp.ArcGISProject("CURRENT")
                m = aprx.activeMap            
                for lyr in m.listLayers():
                    if lyr.isRasterLayer:
                        if in_backsc == lyr.name: 
                            in_backsc = helper.convert_backslash_forwardslash(lyr.dataSource)
        else:
            pass
        
        if in_features.rfind("/") < 0:
            aprx = arcpy.mp.ArcGISProject("CURRENT")
            m = aprx.activeMap            
            for lyr in m.listLayers():
                if lyr.isFeatureLayer:
                    if in_features == lyr.name: 
                        in_features = helper.convert_backslash_forwardslash(lyr.dataSource)
    
            
        tempWS = workspace + "\\temp_tex\\"
        arcpy.env.overwriteOutput = True
        if not os.path.exists(tempWS):
            os.mkdir(tempWS)
        
        arcpy.Copy_management(in_features, tempWS + "deli.shp")
              
        arcpy.env.qualifiedFieldNames = False
        arcpy.env.workspace = workspace
        
        #list of fields that will be saved - check if fields external to the CoMMa toolbox are present and saves them
        keepFields = [field.name for field in arcpy.ListFields(tempWS + "deli.shp")]
        
        #Calculate different texture parameters: ruggedness
        arcpy.AddMessage("Calculating texture parameters ...")
        
        Aspect_ = Aspect(inputDEM)
        
        Slope_ = Slope(inputDEM, "DEGREE", 1)
        
        Rugg_ = helper.Ruggd_calc(3,Aspect_,Slope_)
        Avi_ = helper.AVI(Slope_,Aspect_)
        
        if parameters[1].value is not None:
            List_dev = [Rugg_, Avi_, in_backsc]
            List_nam = ["Rugg", "Avi", "Bsc"]
        
        else:
            List_dev = [Rugg_, Avi_]
            List_nam = ["Rugg", "Avi"]
            
        c=0
        for ra in List_dev:
            ra_name = List_nam[c]
            ZonalStatisticsAsTable(tempWS + "deli.shp", "FID", ra, tempWS + ra_name, 
                                   "DATA", "MIN_MAX_MEAN")
            c+=1
        
        arcpy.MakeFeatureLayer_management(tempWS + "deli.shp", tempWS + "deli_lyr")
        arcpy.AddJoin_management(tempWS + "deli_lyr", "FID", 
                                 tempWS + "Rugg", "FID")
        arcpy.AddJoin_management(tempWS + "deli_lyr", "deli.FID", 
                                 tempWS + "Avi", "FID")
        
        if parameters[1].value is not None:
            arcpy.AddJoin_management(tempWS + "deli_lyr", "deli.FID", 
                                     tempWS + "Bsc", "FID")
        
        arcpy.AddMessage("Changing names to original fields ...")
        
        arcpy.CreateFileGDB_management(tempWS[:-1], "field_alter.gdb")
        arcpy.FeatureClassToGeodatabase_conversion(tempWS + "deli_lyr",  tempWS + "field_alter.gdb")
        base_lyr = tempWS + "field_alter.gdb/deli_lyr"
        
        arcpy.AlterField_management(base_lyr, "MIN", "Rugg_MIN", "Rugg_MIN")
        arcpy.AlterField_management(base_lyr, "MAX", "Rugg_MAX", "Rugg_MAX")
        arcpy.AlterField_management(base_lyr, "MEAN", "Rugg_MEAN", "Rugg_MEAN")
        arcpy.AlterField_management(base_lyr, "MIN_1", "Avi_MIN", "Avi_MIN")
        arcpy.AlterField_management(base_lyr, "MAX_1", "Avi_MAX", "Avi_MAX")
        arcpy.AlterField_management(base_lyr, "MEAN_1", "Avi_MEAN", "Avi_MEAN")
        
        if parameters[1].value is not None:
            arcpy.AlterField_management(base_lyr, "MIN_1_2", "Bsc_MIN", "Bsc_MIN")
            arcpy.AlterField_management(base_lyr, "MAX_1_2", "Bsc_MAX", "Bsc_MAX")
            arcpy.AlterField_management(base_lyr, "MEAN_1_2", "Bsc_MEAN", "Bsc_MEAN")
        
        arcpy.conversion.FeatureClassToShapefile(base_lyr, tempWS[:-1])
        arcpy.Copy_management(tempWS + "deli_lyr.shp", tempWS + "sh.shp")

        
        #Save file and final cleaning
        arcpy.AddMessage('Final cleaning ...') 
        keepFields.extend(["Rugg_MIN","Rugg_MAX","Rugg_MEAN","Avi_MIN","Avi_MAX",
                           "Avi_MEAN", "Bsc_MIN", "Bsc_MEAN", "Bsc_MAX"])
        helper.field_removal(tempWS + "sh.shp", keepFields)
        
        arcpy.AddMessage('Saving official shapefile ...')
        arcpy.Copy_management(tempWS + "sh.shp", workspace + "/" + outFeat)
        
        #Delete temporary files if requested
        if str(delTemp) == 'true':
            arcpy.AddMessage('Deleting temporary files ...')
            helper.cleaning(tempWS[:-1])
            
        #Print file report
        txtFile = open(workspace + "/" + str(outFeat[:-4]) + "_Info_Textu_Desc.txt", "w")
        txtFile.write("Script: CoMMa Delineation ToolBox v1.0" "\n")
        txtFile.write("\n")
        txtFile.write("File Name: " + str(outFeat) + "\n")
        txtFile.write("Input DEM: " + os.path.basename(inputDEM) + "\n")
        txtFile.write("Input features: " + os.path.basename(in_features) + "\n")
        
        txtFile.write("\n")
        if parameters[1].value is not None:
            txtFile.write("Input backscatter file: " + os.path.basename(in_backsc) + "\n")
        else:
            txtFile.write("Backscatter was not used")

        txtFile.close()
        
        
        return


class Volume_descriptor(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Volume descriptor"
        self.description = "Volumetric description tool"
        self.canRunInBackground = False

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
            displayName="Input delineated morphologies",
            name="in_features",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")
        
        param2 = arcpy.Parameter(
            displayName="Are you mapping positive or negative features?",
            name="in_fill_direc",
            datatype="GPString",
            parameterType="Required",
            direction="Input")
        
        # Set a value list of 1, 10 and 100
        param2.filter.type = "ValueList"
        param2.filter.list = ["Positive", "Negative"]
        
        param3 = arcpy.Parameter(
            displayName="Workspace",
            name="workspace",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input") 

        param4 = arcpy.Parameter(
            displayName="Output Feature name",
            name="outFeat",
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

        if parameters[1].valueAsText is not None:
            folder = os.path.basename(os.path.dirname(parameters[1].valueAsText))
            #desc = arcpy.Describe(gdb)
            if folder.lower().endswith(('.gdb', '.mdb')):
                parameters[1].setErrorMessage("Geodatabases cannot be used in this version of the CoMMa Toolbox")
        
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
        in_features = parameters[1].valueAsText
        fillDirec = parameters[2].valueAsText
        workspace = parameters[3].valueAsText
        outFeat = parameters[4].valueAsText
        delTemp = parameters[5].valueAsText
        
        # enable the helper functions
        helper = helpers()
        
        inputDEM = helper.convert_backslash_forwardslash(inputDEM)
        in_features = helper.convert_backslash_forwardslash(in_features)
        workspace = helper.convert_backslash_forwardslash(workspace)
        
        if outFeat[-4:] != '.shp':
            outFeat = outFeat + '.shp'
            
        if inputDEM.rfind("/") < 0:
            aprx = arcpy.mp.ArcGISProject("CURRENT")
            m = aprx.activeMap            
            for lyr in m.listLayers():
                if lyr.isRasterLayer:
                    if inputDEM == lyr.name: 
                        inputDEM = helper.convert_backslash_forwardslash(lyr.dataSource)
                        
        if in_features.rfind("/") < 0:
            aprx = arcpy.mp.ArcGISProject("CURRENT")
            m = aprx.activeMap            
            for lyr in m.listLayers():
                if lyr.isFeatureLayer:
                    if in_features == lyr.name: 
                        in_features = helper.convert_backslash_forwardslash(lyr.dataSource)
    
        tempWS = workspace + "\\temp_vol\\"
        arcpy.env.overwriteOutput = True
        if not os.path.exists(tempWS):
            os.mkdir(tempWS)
        
        Cellsize_Xresult = arcpy.GetRasterProperties_management(inputDEM, "CELLSIZEX")
        Cellsize_X = float((Cellsize_Xresult.getOutput(0)))
        
        arcpy.CreateFileGDB_management(tempWS[:-1], "volum_calc.gdb")
        
        arcpy.FeatureClassToGeodatabase_conversion(in_features, tempWS + "volum_calc.gdb")
        arcpy.RasterToGeodatabase_conversion(inputDEM, tempWS + "volum_calc.gdb")
        
        arcpy.env.qualifiedFieldNames = False
        #arcpy.env.workspace = workspace
        arcpy.env.workspace = tempWS + "volum_calc.gdb"
        
        for fc in arcpy.ListFeatureClasses():
            poly1 = tempWS + "volum_calc.gdb/{}".format(fc)
            arcpy.AddMessage(fc)
        for raster in arcpy.ListRasters():
            arcpy.AddMessage(raster)
            ras1 = tempWS + "volum_calc.gdb/{}".format(raster)
        
        #list of fields that will be saved - check if fields external to the 
        #CoMMa toolbox are present and saves them
        keepFields = [field.name for field in arcpy.ListFields(poly1)]
        
        #Check if maximum depth is present or if it needs to be calculated
        fie = [field.name for field in arcpy.ListFields(poly1)]
        
        arcpy.AddMessage("Calculating volumes and optimal height measurement...")
        
        arcpy.AddField_management(poly1, "Optim_VR", "DOUBLE", "8", "2")
        arcpy.AddField_management(poly1, "Volume", "DOUBLE", "8", "6") 
        
        with arcpy.da.UpdateCursor(poly1, ['SHAPE@', 'OBJECTID', 'Optim_VR', 'Volume']) as cursor:
            
            for row in cursor:
                
                arcpy.AddMessage('processing object '+str(row[1]))
                
                #clip the mound
                clip = arcpy.Clip_management(ras1, "", 'clip', row[0], '-9999', 'ClippingGeometry')
                
                arcpy.env.extent = inputDEM
                
                #expand periphery with boundary values
                neighborhood = NbrCircle(2, "CELL")
                out_focal_statistics = Con(IsNull('clip'), 
                                           FocalStatistics('clip', neighborhood, "MEDIAN"), 'clip')
                for n in range(1,3):
                    out_focal_statistics = Con(IsNull(out_focal_statistics), 
                                               FocalStatistics(out_focal_statistics, neighborhood, "MEDIAN"), out_focal_statistics)
                               
                #extract the created periphery to interpolate the basal surface
                arcpy.RasterToPolygon_conversion(Int(out_focal_statistics), 'expanded', "NO_SIMPLIFY")
                arcpy.analysis.PairwiseErase('expanded', row[0], 'donut')
                
                clip2 = arcpy.Clip_management(out_focal_statistics, "", 'donutr', 'donut', '-9999', 'ClippingGeometry')
                
                #convert the periphery to points
                arcpy.RasterToPoint_conversion('donutr', 'donut_p', "VALUE")
                
                arcpy.env.extent = 'donutr'
                
                #interpolate using a spline
                outSpline = Spline('donut_p', "grid_code", Cellsize_X, "TENSION", 10)
                clip3 = arcpy.Clip_management(outSpline, "", 'splined', row[0], '-9999', 'ClippingGeometry')
                
                subtracted = Minus(clip, clip3)
                
                if fillDirec == "Positive":
                    volume = arcpy.SurfaceVolume_3d(subtracted, "", "ABOVE", 0)
                    maxRSTresult = arcpy.GetRasterProperties_management(subtracted, "MAXIMUM")
                    maxRST = float(maxRSTresult.getOutput(0))
                else:
                    volume = arcpy.SurfaceVolume_3d(subtracted, "", "BELOW", 0)
                    maxRSTresult = arcpy.GetRasterProperties_management(subtracted, "MINIMUM")
                    maxRST = abs(float(maxRSTresult.getOutput(0)))
                    
                message = volume.getMessages()
                                    
                m = re.search('Volume=(\d+\.\d+)', message)
                if m:
                    found = float((m.group(1)))
                else:
                    found = 0
                
                row[2] = maxRST
                row[3] = found
                
                cursor.updateRow(row)
        
        del cursor, row
            
        #Save file and final cleaning
        arcpy.AddMessage('Final cleaning ...') 
        keepFields.extend(["Volume", "Optim_VR"])
        helper.field_removal(poly1, keepFields)
        
        arcpy.env.extent = inputDEM
        
        arcpy.AddMessage('Saving official shapefile ...')
        arcpy.conversion.ExportFeatures(poly1, workspace + "/" + outFeat, "", "NOT_USE_ALIAS")
        
        #Delete temporary files if requested
        if str(delTemp) == 'true':
            arcpy.AddMessage('Deleting temporary files ...')
            helper.cleaning(tempWS[:-1])
        
                
        return

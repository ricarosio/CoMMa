# -*- coding: utf-8 -*-
#Authors: Riccardo Arosio and Joana Gafeira
#Institutions: University College Cork, British Geological Survey
#CoMMa version 1.2 - May 2024

import arcpy
import os

from arcpy.sa import *


class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "Confined morphologies delineators"
        self.alias = "CoMMa_Delineation"
        self.description = ""

        # List of tool classes associated with this toolbox
        self.tools = [Delineate, Delineate_geomorphons]


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
    
    
    def create_remap_string(self,selected_numbers):
        
        remap_list = []
        
        for num in range(1, 11):
            if num in selected_numbers:
                remap_list.append([num, 1])
            else:
                remap_list.append([num, 0])
                
        return remap_list

    
    
    def convert_to_classnum(self,selected_classes):
        
        classes = ["Flat", "Peak", "Ridge", "Shoulder", "Spur", "Slope", "Hollow", "Footslope", "Valley", "Pit"]
        numbers = list(range(1, 11))
        conversion_dict = dict(zip(classes, numbers))
        selected_numbers = [conversion_dict[class_name] for class_name in selected_classes if class_name in conversion_dict]
        
        return selected_numbers

        

class Delineate(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Boundary-based delineation"
        self.description = "This script applies a fill algorithm to delineate confined features on bathymetry data. \
            Three thresholds values have to be defined to run the tool; these are the Minimum Relief, Minimum Width and Minimum Size Ratio. \
            The user will have to also define the Buffer Distance and if the holes inside a delineated feature are removed."
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        parameters = None
        
        param0 = arcpy.Parameter(
            displayName="Input raster DEM",
            name="inputDEM",
            datatype="GPRasterLayer",
            parameterType="Required",
            direction="Input")
        
        param1 = arcpy.Parameter(
            displayName="If you want to delineate using a LPT derivative, insert the raster here",
            name="inputDer",
            datatype="GPRasterLayer",
            parameterType="Optional",
            direction="Input")
        
        param2 = arcpy.Parameter(
            displayName="Select whether you want to map positive or negative features",
            name="in_fill_direc",
            datatype="GPString",
            parameterType="Required",
            direction="Input")
        
        # Set a value list of 1, 10 and 100
        param2.filter.type = "ValueList"
        param2.filter.list = ["Positive", "Negative"]

        param3 = arcpy.Parameter(
            displayName="Vertical Cutoff (unit: m or derivative unit if the latter is used)",
            name="cutoffVR",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        
        param4 = arcpy.Parameter(
            displayName="Minimum Vertical Threshold (unit: m)",
            name="minVR",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        
        param5 = arcpy.Parameter(
            displayName="Minimum Width (unit: m)",
            name="minWidth",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        
        param6 = arcpy.Parameter(
            displayName="Minimum W/L Ratio",
            name="minRatio",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        param6.value = 0
        
        param7 = arcpy.Parameter(
            displayName="Buffer to apply to the delineation (unit: m)",
            name="delBuffer",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        
        param8 = arcpy.Parameter(
            displayName="Select whether you want to simplify and smooth the delineations",
            name="in_smoosimp",
            datatype="GPString",
            parameterType="Required",
            direction="Input")
        
        # Set a value list of 1, 10 and 100
        param8.filter.type = "ValueList"
        param8.filter.list = ["No simplification or smoothing", 
                              "Light simplification and smoothing",
                              "Heavy simplification and smoothing"]
        
        param9 = arcpy.Parameter(
            displayName="Workspace",
            name="workspace",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")  
        
        param10 = arcpy.Parameter(
            displayName="Output Feature name",
            name="outFeat",
            datatype="GPString",
            parameterType="Required",
            direction="Output")
        
        param11 = arcpy.Parameter(
            displayName="Do you want to apply the geomorphons filter?",
            name="WBTyn",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")
        
        param12 = arcpy.Parameter(
            displayName="Insert the geomorphons raster file",
            name="inputGeo",
            datatype="GPRasterLayer",
            parameterType="Optional",
            direction="Input")

        param13 = arcpy.Parameter(
            displayName="Minimum Geomorphons area over total area ratio",
            name="minGeoRatio",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")   
        
        param14 = arcpy.Parameter(
            displayName="Do you want to delete internal holes in the polygons?",
            name="delHoles",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")
        
        param15 = arcpy.Parameter(
            displayName="Do you want to delete the temporary files?",
            name="delTemp",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")
        
        
        parameters = [param0, param1, param2, param3, param4, param5, param6, 
                      param7, param8, param9, param10, param11, param12, param13, 
                      param14, param15]
    
        return parameters
    

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        
        if parameters[11].value == True: #checking if the checkmark box is checked
           #if the box is checked (true), enable parameters 11 and 12
           parameters[12].enabled = True
           parameters[13].enabled = True
           
        else:
           parameters[12].enabled = False
           parameters[13].enabled = False
           
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
        
        if parameters[9].valueAsText is not None:
            folder = os.path.basename(parameters[9].valueAsText)
            #desc = arcpy.Describe(gdb)
            if folder.lower().endswith(('.gdb', '.mdb')):
                parameters[9].setErrorMessage("Geodatabases cannot be used in this version of the CoMMa Toolbox")

        if parameters[12].valueAsText is not None:
            folder = os.path.basename(os.path.dirname(parameters[11].valueAsText))
            #desc = arcpy.Describe(gdb)
            if folder.lower().endswith(('.gdb', '.mdb')):
                parameters[12].setErrorMessage("Geodatabases cannot be used in this version of the CoMMa Toolbox")
            
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        
        arcpy.env.overwriteOutput = True
        inputDEM = parameters[0].valueAsText
        inputDer = parameters[1].valueAsText
        fillDirec = parameters[2].valueAsText
        cutoffVR = parameters[3].valueAsText
        minVR = parameters[4].valueAsText
        minWidth = parameters[5].valueAsText
        minRatio = parameters[6].valueAsText
        delBuffer = parameters[7].valueAsText
        in_smoosimp = parameters[8].valueAsText
        workspace = parameters[9].valueAsText
        outFeat = parameters[10].valueAsText
        WBTyn = parameters[11].valueAsText
        inputGeo = parameters[12].valueAsText
        minGeoRatio = parameters[13].valueAsText
        delHoles = parameters[14].valueAsText
        delTemp = parameters[15].valueAsText
        
        # enable the helper functions
        helper = helpers()
        
        inputDEM = helper.convert_backslash_forwardslash(inputDEM)
        inputDer = helper.convert_backslash_forwardslash(inputDer)
        inputGeo = helper.convert_backslash_forwardslash(inputGeo)
        workspace = helper.convert_backslash_forwardslash(workspace)
        
        if outFeat[-4:] != '.shp':
            outFeat = outFeat + '.shp'
        # if the input bathyRas is selected from a drop-down list, the bathyRas does not contain the full path
        # In this case, the full path needs to be obtained from the map layer
        if inputDEM.rfind("/") < 0:
            aprx = arcpy.mp.ArcGISProject("CURRENT")
            m = aprx.activeMap            
            for lyr in m.listLayers():
                if lyr.isRasterLayer:
                    if inputDEM == lyr.name: 
                        inputDEM = helper.convert_backslash_forwardslash(lyr.dataSource)
        
        if parameters[1].value is not None:
            if inputDer.rfind("/") < 0:
                aprx = arcpy.mp.ArcGISProject("CURRENT")
                m = aprx.activeMap            
                for lyr in m.listLayers():
                    if lyr.isRasterLayer:
                        if inputDer == lyr.name: 
                            inputDer = helper.convert_backslash_forwardslash(lyr.dataSource)
        else:
            pass
        
        if parameters[12].value is not None:
            if inputGeo.rfind("/") < 0:
                aprx = arcpy.mp.ArcGISProject("CURRENT")
                m = aprx.activeMap            
                for lyr in m.listLayers():
                    if lyr.isRasterLayer:
                        if inputGeo == lyr.name: 
                            inputGeo = helper.convert_backslash_forwardslash(lyr.dataSource)
        else:
            pass
        
        Cellsize_Xresult = arcpy.GetRasterProperties_management(inputDEM, "CELLSIZEX")
        Cellsize_X = float(Cellsize_Xresult.getOutput(0))
        
        tempWS = workspace + "\\temp\\"
        arcpy.env.overwriteOutput = True
        if not os.path.exists(tempWS):
            os.mkdir(tempWS)
              
        arcpy.env.qualifiedFieldNames = False
        arcpy.env.workspace = workspace
        
        #list of fields that will be saved
        keepFields = []

        # Process: Invert DEM and fill
        arcpy.AddMessage("Filling bathymetry ...")
        # Invert the input DEM to map positive topographic features
        if fillDirec == "Positive":           
            # Execute Negate
            invertedInputDEM = Negate(inputDEM)
            featuresDEM = invertedInputDEM
            # Process: Fill
            filledDEM = Fill(featuresDEM)
            filledDEM.save(tempWS + "a1f.tif")
            
        else: 
        # Process: Fill
            featuresDEM = inputDEM
            filledDEM = Fill(featuresDEM)
            filledDEM.save(tempWS + "a1f.tif")
        
        # Process: Minus 
        arcpy.AddMessage("Subtracting ...")
        minusDEM = Minus(filledDEM, featuresDEM)
        minusDEM.save(tempWS + "a1m.tif")
        
        arcpy.AddMessage("Delineating features ...")
        #Delineating features
        #RECLASSIFY the raster using the Minimum Vertical Threshold.
        #Areas with values < cutoffVR will be classed as 1 and areas with values > Cutoff will be classed as 2
        reclassifyRanges = "" + cutoffVR + " 9000 1"
        
        if parameters[1].value is not None:
            arcpy.AddMessage("Preparing derivative ...")
            if fillDirec == "Negative":  
                # Execute Negate
                featuresDer = Negate(inputDer)
                saveDer = featuresDer
                 
            else:
                featuresDer = inputDer
                saveDer = Raster(featuresDer)

            saveDer.save(tempWS + "a1md.tif")
            
            OutReclassed = Reclassify(featuresDer, "Value", reclassifyRanges, "NODATA")
            
        else:
            
            OutReclassed = Reclassify(tempWS + "a1m.tif", "Value", reclassifyRanges, "NODATA")
        
        
        #Convert to feature class
        arcpy.RasterToPolygon_conversion(OutReclassed, tempWS + "b3p.shp", "NO_SIMPLIFY")
        
        if arcpy.management.GetCount(tempWS + "b3p.shp")[0] == "0":
            
            arcpy.AddMessage("No features detected, please change the parameters ...")
        
        else:
            
            if float(delBuffer) > 0:
                arcpy.AddMessage("Buffering the delineated features ...")
                arcpy.analysis.Buffer(tempWS + "b3p.shp", tempWS + "b3pbu.shp", delBuffer + " Meters", 
                                      "FULL", "ROUND", "NONE")
                arcpy.management.Dissolve(tempWS + "b3pbu.shp", tempWS + "b3pf.shp", "", "", 
                                          "SINGLE_PART")
                
            else:
                arcpy.management.Rename(tempWS + "b3p.shp", tempWS + "b3pf.shp")
                   
            
            #Characterise the features' geometry _ Part 2
            #Calculate the minimum bounding geometry (MBG) for each polygon and adding the 
            #following field: MBG_Width, MBG_Length, MBG_Orientation
            
            arcpy.MinimumBoundingGeometry_management(tempWS + "b3pf.shp", tempWS + "c4b.shp",
                                                     "RECTANGLE_BY_WIDTH", "NONE", "", "MBG_FIELDS")
    
            #Add and calculate the MBG Width/Length field
            arcpy.AddField_management(tempWS + "c4b.shp", "MBG_W_L", "DOUBLE", "8", "2")
            
            arcpy.CalculateField_management(tempWS + "c4b.shp", "MBG_W_L",
                                            '!MBG_Width!/!MBG_Length!', "PYTHON3")
            
            #Join the MBG fields to the feature outline polygon shapefile
            arcpy.SpatialJoin_analysis(tempWS + "b3pf.shp", tempWS + "c4b.shp", tempWS + "d5g.shp",
                                       "JOIN_ONE_TO_MANY","KEEP_ALL", "","WITHIN")
            
            #Delete polygons that got the attributes from another MBG
            if not arcpy.Exists(tempWS + "d5g_lyr"):
                arcpy.MakeFeatureLayer_management(tempWS + "d5g.shp", tempWS + "d5g_lyr")
            
            arcpy.SelectLayerByAttribute_management(tempWS + "d5g_lyr", "NEW_SELECTION",
                                                    '("TARGET_FID" <> "JOIN_FID")')
            
            arcpy.DeleteRows_management(tempWS + "d5g_lyr")
            
            #Delete useless field
            
            arcpy.Delete_management(tempWS + "d5g_lyr")
            
            # Process: Selection on area and elongation
            # Select features with an area less than Minimum Area threshold and with an atypical shape ratio
            arcpy.AddMessage("Deleting features with values below the given area and elongation thresholds ...")
               
            arcpy.CopyFeatures_management(tempWS + "d5g.shp", tempWS + "e6d.shp")
            
            if not arcpy.Exists(tempWS + "e6d_lyr"):
                arcpy.MakeFeatureLayer_management(tempWS + "e6d.shp", tempWS + "e6d_lyr")
            
            arcpy.SelectLayerByAttribute_management(tempWS + "e6d_lyr", "NEW_SELECTION",
                                                    '("MBG_Width" < ' + minWidth + ') OR ("MBG_W_L" < ' + minRatio + ')')
            
            #Delete polygons with an area inferior to MinArea or with an elongation inferior to MinRatio
            arcpy.DeleteRows_management(tempWS + "e6d_lyr")
            arcpy.Delete_management(tempWS + "e6d_lyr")
            
            if not arcpy.Exists(tempWS + "e6d_lyr"):
                arcpy.MakeFeatureLayer_management(tempWS + "e6d.shp", tempWS + "e6d_lyr")
            
            
            #Reshaping polygon containing holes
            if delHoles:
                arcpy.AddMessage("Deleting holes ...")
                arcpy.EliminatePolygonPart_management(tempWS + "e6d.shp", tempWS + "e6h.shp","PERCENT",
                                                      "","55","CONTAINED_ONLY")
            else:
                arcpy.CopyFeatures_management(tempWS + "e6d.shp", tempWS + "e6h.shp")
            
            arcpy.Delete_management(tempWS + "e6d_lyr")
            
            #Getting Maximum Relief values from Mosaic raster
            #Use the created shapefile as a mask for extraction of the Maximum Relief 
            #value from the Mosaic raster
            arcpy.AddMessage("Extracting Relief ...")
            
            #Add field for the Maximum Relief values
            arcpy.AddField_management(tempWS + "e6h.shp", "Relief", "DOUBLE", "8", "2")
            
            #extract the maximum height of each polygon feature
            ZonalStatisticsAsTable(tempWS + "e6h.shp", "FID", 
                                   tempWS + "a1m.tif", tempWS + "f7vt", 
                                   "NODATA", "MAXIMUM")
            
            if not arcpy.Exists(tempWS + "e6h_lyr"):
                arcpy.MakeFeatureLayer_management(tempWS + "e6h.shp", tempWS + "e6h_lyr")
            
            arcpy.AddJoin_management(tempWS + "e6h_lyr", "FID", 
                                     tempWS + "f7vt", "FID", "KEEP_ALL")
            
            arcpy.CopyFeatures_management(tempWS + "e6h_lyr", tempWS + "f7jo.shp")
            
            relExpression = "float(!MAX!)"
            
            arcpy.CalculateField_management(tempWS + "f7jo.shp", "Relief", relExpression, "PYTHON3")
    
            
            #Delete polygons with a relief smaller than the Minimum Vertical Threshold (MinVR)
            if not arcpy.Exists(tempWS + "f7jo_lyr"):
                arcpy.MakeFeatureLayer_management(tempWS + "f7jo.shp", tempWS + "f7jo_lyr")
            
            arcpy.SelectLayerByAttribute_management(tempWS + "f7jo_lyr", 
                                                    "NEW_SELECTION", '"Relief" <  %d' %float(minVR))
            arcpy.DeleteRows_management(tempWS + "f7jo_lyr")
            
            arcpy.CopyFeatures_management(tempWS + "f7jo_lyr", tempWS + "f7se.shp")
            
            keepFields.extend(["Relief","MBG_Width","MBG_Length",
                               "MBG_W_L","MBG_Orient"])
            helper.field_removal(tempWS + "f7se.shp", keepFields)
            
            #Add field with AREA values
            arcpy.AddMessage("Characterising features geometry ...")
    
            arcpy.AddField_management(tempWS + "f7se.shp", "Area", "DOUBLE", "8", "1")
            areaExpression = "float(!SHAPE.AREA!)"
            arcpy.CalculateField_management(tempWS + "f7se.shp", "Area", areaExpression,
                                            "PYTHON3")
            
            #Add field with the PERIMETER values
            arcpy.AddField_management(tempWS + "f7se.shp", "Perimeter", "DOUBLE", "8", "1")
            perimeterExpression = "float(!SHAPE.LENGTH!)"
            arcpy.CalculateField_management(tempWS + "f7se.shp", "Perimeter", perimeterExpression, "PYTHON3")
            
            ######################################################################
            ######################################################################
            # START GEOMORPHONS###################################################
            
            if str(WBTyn) == 'true':
                
                arcpy.AddMessage("Processing geomorphons ...")
                
                if fillDirec == "Positive":
                    remapString  = "2 5 1; 1 0; 6 10 0"
                else:
                    remapString  = "1 6 0; 7 1; 8 0; 9 10 1"
                
                OutRecla_Rid = Reclassify(inputGeo, "Value", remapString, "NODATA")
                
                outMajFilt_Rid = MajorityFilter(OutRecla_Rid, "EIGHT", "HALF")
                
                arcpy.AddField_management(tempWS + "f7se.shp", "Geo_rat", "DOUBLE", "8", "2")
                
                ZonalStatisticsAsTable(tempWS + "f7se.shp", "FID", 
                                       outMajFilt_Rid, tempWS + "f7ge", 
                                       "NODATA", "SUM")
                
                if not arcpy.Exists(tempWS + "rcalc5_lyr"):
                    arcpy.MakeFeatureLayer_management(tempWS + "f7se.shp", tempWS + "rcalc5_lyr")
                
                arcpy.AddJoin_management(tempWS + "rcalc5_lyr", "FID", 
                                         tempWS + "f7ge", "FID", "KEEP_ALL")
                
                arcpy.CopyFeatures_management(tempWS + "rcalc5_lyr", tempWS + "f7geo.shp")
                
                Pixel_area = Cellsize_X*Cellsize_X
                relExpression = "(float(!SUM!) * %d) / !Area!" %Pixel_area
                
                arcpy.CalculateField_management(tempWS + "f7geo.shp", "Geo_rat", 
                                                relExpression, "PYTHON3")
                
                # Delete polygons with a ratio geomorphons positive classes/ area smaller than 
                #the Minimum Geomorphons ratio threshold (minGeoRatio)
                if not arcpy.Exists(tempWS + "f7geo_lyr"):
                    arcpy.MakeFeatureLayer_management(tempWS + "f7geo.shp", tempWS + "f7geo_lyr")
                
                arcpy.SelectLayerByAttribute_management(tempWS + "f7geo_lyr", "NEW_SELECTION", 
                                                        '"Geo_rat" <  %d' %float(minGeoRatio))
                arcpy.DeleteRows_management(tempWS + "f7geo_lyr")
                
                arcpy.CopyFeatures_management(tempWS + "f7geo_lyr", tempWS + "f7seb.shp")
                poly1 = tempWS + "f7seb.shp"
            
            else:
                poly1 = tempWS + "f7se.shp"
    
            ######################################################################
            ######################################################################
            # END GEOMORPHONS#####################################################                  

            #Reshaping the delineated polygons
            arcpy.AddMessage("Reshaping the delineated polygons ...")
            
            if not arcpy.Exists(tempWS + "g8b_lyr"):
                arcpy.MakeFeatureLayer_management(poly1, tempWS + "g8b_lyr")
            
            fields = ['MBG_Width']

            all_rows = [i[0] for i in arcpy.da.SearchCursor(poly1,fields[0])]
            min_width = min(all_rows)
            arcpy.AddMessage(min_width)
            
            if in_smoosimp == "No simplification or smoothing":
                arcpy.AddMessage("No simplification or smoothing")
                arcpy.CopyFeatures_management(tempWS + "g8b_lyr", tempWS + "h9z.shp")
            
            elif in_smoosimp == "Light simplification and smoothing":
                
                arcpy.AddMessage("Light simplification and smoothing...")
        
                # Simplifying
                arcpy.SimplifyPolygon_cartography(tempWS + "g8b_lyr", tempWS + "h9s.shp", "WEIGHTED_AREA", "%d" %float(min_width/10))
                arcpy.Delete_management(tempWS + "g8b_lyr")
                
                # Smoothing
                if not arcpy.Exists(tempWS + "h9s_lyr"):
                    arcpy.MakeFeatureLayer_management(tempWS + "h9s.shp", tempWS + "h9s_lyr")
                
                arcpy.SmoothPolygon_cartography(tempWS + "h9s_lyr", tempWS + "h9z.shp", 
                                                "PAEK", "%d Meters" %float(min_width/5), "NO_FIXED", "NO_CHECK")
        
                arcpy.Delete_management(tempWS + "h9s_lyr")
                
            elif in_smoosimp == "Heavy simplification and smoothing":
                
                arcpy.AddMessage("Heavy simplification and smoothing...")
        
                # Simplifying
                arcpy.SimplifyPolygon_cartography(tempWS + "g8b_lyr", tempWS + "h9s.shp", "WEIGHTED_AREA", "%d" %float(min_width/5))
                arcpy.Delete_management(tempWS + "g8b_lyr")
                
                # Smoothing
                if not arcpy.Exists(tempWS + "h9s_lyr"):
                    arcpy.MakeFeatureLayer_management(tempWS + "h9s.shp", tempWS + "h9s_lyr")
                
                arcpy.SmoothPolygon_cartography(tempWS + "h9s_lyr", tempWS + "h9z.shp", 
                                                "PAEK", "%d Meters" %float(min_width/2), "NO_FIXED", "NO_CHECK")
        
                arcpy.Delete_management(tempWS + "h9s_lyr")
           
            #Final operations
            arcpy.AddMessage('Final cleaning ...')    
            
            #Delete useless fields
            keepFields.extend(["Area","Perimeter","Relief","MBG_Width","MBG_Length",
                               "MBG_W_L","MBG_Orient","Geo_rat"])
            helper.field_removal(tempWS + "h9z.shp", keepFields)
            
            #Save file
            arcpy.AddMessage('Saving official shapefile ...')
            arcpy.Copy_management(tempWS + "/" + "h9z.shp", workspace + "/" + outFeat)
            
            #Delete temporary files if requested
            if str(delTemp) == 'true':
                arcpy.AddMessage('Deleting temporary files ...')
                helper.cleaning(tempWS[:-1])
                
            #Print file report
            txtFile = open(workspace + "/" + str(outFeat[:-4]) + "_Info.txt", "w")
            txtFile.write("Script: CoMMa Delineation ToolBox v1.0" "\n")
            txtFile.write("\n")
            txtFile.write("File Name: " + str(outFeat) + "\n")
            txtFile.write("Input DEM: " + os.path.basename(inputDEM) + "\n")
            
            if parameters[1].value is not None:
                txtFile.write("Delineation type: Fill algorithm on derivative" + "\n")
                txtFile.write("Input derivative: " + os.path.basename(inputDer) + "\n")
            else:
                txtFile.write("Delineation type: Fill algorithm on bathymetry" + "\n")
                txtFile.write("Input derivative: none" + "\n")
                
            txtFile.write("\n")
            txtFile.write("Relief type: " + fillDirec + "\n")
            txtFile.write("\n")
            txtFile.write("Input cell size: " + str(Cellsize_X) + "\n")
            txtFile.write("Cutoff value: " + cutoffVR + "\n")
            txtFile.write("Minimum relief value: " + minVR + "\n")
            txtFile.write("Minimum Width (m): " + minWidth + "\n")
            txtFile.write("Minimum W/L Ratio: " + minRatio + "\n")
            txtFile.write("Buffer extent (m): " + delBuffer + "\n")
            txtFile.write("Simplification or smoothing?: " + in_smoosimp + "\n")
            txtFile.write("\n")
            if WBTyn:
                txtFile.write("Input geomorphons file: " + os.path.basename(inputGeo) + "\n")
                txtFile.write("Geomorphons positives vs area ratio: " + minGeoRatio + "\n")
            else:
                txtFile.write("The geomorphons add-on was not implemented")
    
            txtFile.close()
        
        '''
        pw = workspace + "/" + outFeat
        aprx = arcpy.mp.ArcGISProject("CURRENT")
        aprxMap = aprx.listMaps("MainMap")[0] 
        aprxMap.addDataFromPath(pw)
        '''  
        return


class Delineate_geomorphons(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Elements-based delineation"
        self.description = "This script allows the semi-automated delineation of confined features directly from the bathymetric data. \
            The geomorphons tool is used to create 10 land surface facets that are grouped together to create the positive or negative relief \
            The user will have to also define the Buffer Distance and if the holes inside a delineated feature are removed."
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        parameters = None
        
        param0 = arcpy.Parameter(
            displayName="Input raster DEM",
            name="inputDEM",
            datatype="GPRasterLayer",
            parameterType="Required",
            direction="Input")

        param1 = arcpy.Parameter(
            displayName="Input Geomorphons raster",
            name="inputGeom",
            datatype="GPRasterLayer",
            parameterType="Required",
            direction="Input")
        
        param2 = arcpy.Parameter(
            displayName="Select the core geomorphons you wish to map",
            name="co_gmf_class",
            datatype="GPString",
            parameterType="Required",
            direction="Input",
            multiValue="True")
        
        param2.filter.type = "ValueList"
        param2.filter.list = ["Flat", "Peak", "Ridge", "Shoulder", "Spur",
                              "Slope", "Hollow", "Footslope", "Valley", "Pit"]
        
        param3 = arcpy.Parameter(
            displayName="Select the subordinate geomorphons you wish to map",
            name="su_gmf_class",
            datatype="GPString",
            parameterType="Optional",
            direction="Input",
            multiValue="True")
        
        param3.filter.type = "ValueList"
        param3.filter.list = ["Flat", "Peak", "Ridge", "Shoulder", "Spur",
                              "Slope", "Hollow", "Footslope", "Valley", "Pit"]
        
        param4 = arcpy.Parameter(
            displayName="Minimum Vertical Threshold (unit: m)",
            name="minVR",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        
        param5 = arcpy.Parameter(
            displayName="Minimum Width (unit: m)",
            name="minWidth",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        
        param6 = arcpy.Parameter(
            displayName="Minimum W/L Ratio",
            name="minRatio",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        
        param7 = arcpy.Parameter(
            displayName="Buffer to apply to the delineation (unit: m)",
            name="delBuffer",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        
        param8 = arcpy.Parameter(
            displayName="Select whether you want to simplify and smooth the delineations",
            name="in_smoosimp",
            datatype="GPString",
            parameterType="Required",
            direction="Input")
        
        # Set a value list of 1, 10 and 100
        param8.filter.type = "ValueList"
        param8.filter.list = ["No simplification or smoothing", 
                              "Light simplification and smoothing",
                              "Heavy simplification and smoothing"]
        
        param9 = arcpy.Parameter(
            displayName="Workspace",
            name="workspace",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")  
        
        param10 = arcpy.Parameter(
            displayName="Output Feature name",
            name="outFeat",
            datatype="GPString",
            parameterType="Required",
            direction="Output")
                
        param11 = arcpy.Parameter(
            displayName="Do you want to delete internal holes in the polygons?",
            name="delHoles",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")
        
        param12 = arcpy.Parameter(
            displayName="Do you want to delete the temporary files?",
            name="delTemp",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")
        
        
        parameters = [param0, param1, param2, param3, param4, param5, param6, 
                      param7, param8, param9, param10, param11, param12]
    
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
        
        if parameters[9].valueAsText is not None:
            folder = os.path.basename(parameters[9].valueAsText)
            #desc = arcpy.Describe(gdb)
            if folder.lower().endswith(('.gdb', '.mdb')):
                parameters[9].setErrorMessage("Geodatabases cannot be used in this version of the CoMMa Toolbox")
                
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        
        arcpy.env.overwriteOutput = True
        inputDEM = parameters[0].valueAsText
        inputGeom = parameters[1].valueAsText
        co_gmf_class = parameters[2].valueAsText
        su_gmf_class = parameters[3].valueAsText
        minVR = parameters[4].valueAsText
        minWidth = parameters[5].valueAsText
        minRatio = parameters[6].valueAsText
        delBuffer = parameters[7].valueAsText
        in_smoosimp = parameters[8].valueAsText
        workspace = parameters[9].valueAsText
        outFeat = parameters[10].valueAsText
        delHoles = parameters[11].valueAsText
        delTemp = parameters[12].valueAsText
        
        # enable the helper functions
        helper = helpers()
        
        inputDEM = helper.convert_backslash_forwardslash(inputDEM)
        inputGeom = helper.convert_backslash_forwardslash(inputGeom)
        workspace = helper.convert_backslash_forwardslash(workspace)
        
        if outFeat[-4:] != '.shp':
            outFeat = outFeat + '.shp'
        # if the input bathyRas is selected from a drop-down list, the bathyRas does not contain the full path
        # In this case, the full path needs to be obtained from the map layer
        if inputDEM.rfind("/") < 0:
            aprx = arcpy.mp.ArcGISProject("CURRENT")
            m = aprx.activeMap            
            for lyr in m.listLayers():
                if lyr.isRasterLayer:
                    if inputDEM == lyr.name: 
                        inputDEM = helper.convert_backslash_forwardslash(lyr.dataSource)
                        
        if parameters[1].value is not None:
            if inputGeom.rfind("/") < 0:
                aprx = arcpy.mp.ArcGISProject("CURRENT")
                m = aprx.activeMap            
                for lyr in m.listLayers():
                    if lyr.isRasterLayer:
                        if inputGeom == lyr.name: 
                            inputGeom = helper.convert_backslash_forwardslash(lyr.dataSource)

        else:
            pass
        
        Cellsize_Xresult = arcpy.GetRasterProperties_management(inputDEM, "CELLSIZEX")
        Cellsize_X = float(Cellsize_Xresult.getOutput(0))
        
        tempWS = workspace + "\\temp\\"
        arcpy.env.overwriteOutput = True
        if not os.path.exists(tempWS):
            os.mkdir(tempWS)
              
        arcpy.env.qualifiedFieldNames = False
        arcpy.env.workspace = workspace
        
        #list of fields that will be saved
        keepFields = []
        
        #Delineating features
        arcpy.AddMessage("Delineating features ...")
        
        co_gmf_class = co_gmf_class.split(";")
        su_gmf_class = su_gmf_class.split(";")

        arcpy.AddMessage(co_gmf_class)
        arcpy.AddMessage(su_gmf_class)
        selected_core = helper.convert_to_classnum(co_gmf_class)
        selected_subord = helper.convert_to_classnum(su_gmf_class)
        
        size_1 = len(selected_core)
        size_2 = len(selected_subord)
        res = []
        i, j = 0, 0
         
        while i < size_1 and j < size_2:
            if selected_core[i] < selected_subord[j]:
                res.append(selected_core[i])
                i += 1
         
            else:
                res.append(selected_subord[j])
                j += 1
         
        selected_total = res + selected_core[i:] + selected_subord[j:]
        
        arcpy.AddMessage(selected_total)
        remap_string_to = helper.create_remap_string(selected_total)
        arcpy.AddMessage(remap_string_to)
        remap_string_co = helper.create_remap_string(selected_core)
        arcpy.AddMessage(remap_string_co)

        #remapString  = "2 5 1; 1 0; 6 10 0"
        
        core_reclass = Reclassify(inputGeom, "Value", RemapValue(remap_string_co), "NODATA")
        total_reclass = Reclassify(inputGeom, "Value", RemapValue(remap_string_to), "NODATA")
        
        total_reclass_filt = MajorityFilter(total_reclass, "EIGHT", "HALF")
        total_reclass_filt.save(tempWS + "b3p.tif")
        
        remapNODATA = "0 NODATA;1 1"
        total_reclass_filt2 = Reclassify(total_reclass_filt, "Value", remapNODATA, "NODATA")
        core_reclass2 = Reclassify(core_reclass, "Value", remapNODATA, "NODATA")
      
        #Convert to feature class
        arcpy.RasterToPolygon_conversion(total_reclass_filt2, tempWS + "b3p.shp", "NO_SIMPLIFY")
        arcpy.RasterToPolygon_conversion(core_reclass2, tempWS + "b3c.shp", "NO_SIMPLIFY")
        
        if arcpy.management.GetCount(tempWS + "b3p.shp")[0] == "0":
            
            arcpy.AddMessage("No features detected, please change the parameters ...")
        
        else:
            
            if not arcpy.Exists(tempWS + "b3p_lyr"):
                arcpy.MakeFeatureLayer_management(tempWS + "b3p.shp", tempWS + "b3p_lyr")
            
            cored_subordinates = arcpy.management.SelectLayerByLocation(tempWS + "b3p_lyr", 
                                                                        'INTERSECT', tempWS + "b3c.shp")
            
            arcpy.CopyFeatures_management(cored_subordinates, tempWS + "b3f.shp")
            
            if float(delBuffer) > 0:
                arcpy.AddMessage("Buffering the delineated features ...")
                arcpy.analysis.Buffer(tempWS + "b3f.shp", tempWS + "b3fbu.shp", delBuffer + " Meters", 
                                      "FULL", "ROUND", "NONE")
                arcpy.management.Dissolve(tempWS + "b3fbu.shp", tempWS + "b3ff.shp", "", "", 
                                          "SINGLE_PART")
            else:
                arcpy.management.Rename(tempWS + "b3f.shp", tempWS + "b3ff.shp")
            
            #Characterise the features' geometry _ Part 2
            #Calculate the minimum bounding geometry (MBG) for each polygon and adding the 
            #following field: MBG_Width, MBG_Length, MBG_Orientation
            
            arcpy.MinimumBoundingGeometry_management(tempWS + "b3ff.shp", tempWS + "c4b.shp",
                                                     "RECTANGLE_BY_WIDTH", "NONE", "", "MBG_FIELDS")
    
            #Add and calculate the MBG Width/Length field
            arcpy.AddField_management(tempWS + "c4b.shp", "MBG_W_L", "DOUBLE", "8", "2")
            
            arcpy.CalculateField_management(tempWS + "c4b.shp", "MBG_W_L",
                                            '!MBG_Width!/!MBG_Length!', "PYTHON3")
            
            #Join the MBG fields to the feature outline polygon shapefile
            arcpy.SpatialJoin_analysis(tempWS + "b3ff.shp", tempWS + "c4b.shp", tempWS + "d5g.shp",
                                       "JOIN_ONE_TO_MANY","KEEP_ALL", "","WITHIN")
            
            #Delete polygons that got the attributes from another MBG
            if not arcpy.Exists(tempWS + "d5g_lyr"):
                arcpy.MakeFeatureLayer_management(tempWS + "d5g.shp", tempWS + "d5g_lyr")
            
            arcpy.SelectLayerByAttribute_management(tempWS + "d5g_lyr", "NEW_SELECTION",
                                                    '("TARGET_FID" <> "JOIN_FID")')
            
            arcpy.DeleteRows_management(tempWS + "d5g_lyr")
            
            #Delete useless field
            
            arcpy.Delete_management(tempWS + "d5g_lyr")
            
            # Process: Selection on area and elongation
            # Select features with an area less than Minimum Area threshold and with an atypical shape ratio
            arcpy.AddMessage("Deleting features with values below the given area and elongation thresholds ...")
               
            arcpy.CopyFeatures_management(tempWS + "d5g.shp", tempWS + "e6d.shp")
            
            if not arcpy.Exists(tempWS + "e6d_lyr"):
                arcpy.MakeFeatureLayer_management(tempWS + "e6d.shp", tempWS + "e6d_lyr")
            
            arcpy.SelectLayerByAttribute_management(tempWS + "e6d_lyr", "NEW_SELECTION",
                                                    '("MBG_Width" < ' + minWidth + ') OR ("MBG_W_L" < ' + minRatio + ')')
            
            #Delete polygons with an area inferior to MinArea or with an elongation inferior to MinRatio
            arcpy.DeleteRows_management(tempWS + "e6d_lyr")
            arcpy.Delete_management(tempWS + "e6d_lyr")
            
            if not arcpy.Exists(tempWS + "e6d_lyr"):
                arcpy.MakeFeatureLayer_management(tempWS + "e6d.shp", tempWS + "e6d_lyr")
            
            
            #Reshaping polygon containing holes
            if delHoles:
                arcpy.AddMessage("Deleting holes ...")
                arcpy.EliminatePolygonPart_management(tempWS + "e6d.shp", tempWS + "e6h.shp","PERCENT",
                                                      "","55","CONTAINED_ONLY")
            else:
                arcpy.CopyFeatures_management(tempWS + "e6d.shp", tempWS + "e6h.shp")
            
            arcpy.Delete_management(tempWS + "e6d_lyr")
            
            #Getting Maximum Relief values from Mosaic raster
            #Use the created shapefile as a mask for extraction of the Maximum Relief 
            #value from the Mosaic raster
            arcpy.AddMessage("Extracting Relief ...")
            
            #Add field for the Maximum Relief values
            arcpy.AddField_management(tempWS + "e6h.shp", "Relief", "DOUBLE", "8", "2")
            
            #extract the minimum and maximum depth of each polygon feature
            ZonalStatisticsAsTable(tempWS + "e6h.shp", "FID", inputDEM, 
                                   tempWS + "f7vt", "DATA", "MIN_MAX")
            
            if not arcpy.Exists(tempWS + "e6h_lyr"):
                arcpy.MakeFeatureLayer_management(tempWS + "e6h.shp", tempWS + "e6h_lyr")
            
            arcpy.AddJoin_management(tempWS + "e6h_lyr", "FID", 
                                     tempWS + "f7vt", "FID", "KEEP_ALL")
            
            arcpy.CopyFeatures_management(tempWS + "e6h_lyr", tempWS + "f7jo.shp")
            
            relExpression = "float(!MAX!) - float(!MIN!)"
            
            arcpy.CalculateField_management(tempWS + "f7jo.shp", "Relief", relExpression, "PYTHON3")
    
            
            #Delete polygons with a relief smaller than the Minimum Relief Threshold (MinVR)
            if not arcpy.Exists(tempWS + "f7jo_lyr"):
                arcpy.MakeFeatureLayer_management(tempWS + "f7jo.shp", tempWS + "f7jo_lyr")
            
            arcpy.SelectLayerByAttribute_management(tempWS + "f7jo_lyr", 
                                                    "NEW_SELECTION", '"Relief" <  %d' %float(minVR))
            arcpy.DeleteRows_management(tempWS + "f7jo_lyr")
            
            arcpy.CopyFeatures_management(tempWS + "f7jo_lyr", tempWS + "f7se.shp")
            
            keepFields.extend(["Relief","MBG_Width","MBG_Length",
                               "MBG_W_L","MBG_Orient"])
            helper.field_removal(tempWS + "f7se.shp", keepFields)
            
            #Add field with AREA values
            arcpy.AddMessage("Characterising features geometry ...")
    
            arcpy.AddField_management(tempWS + "f7se.shp", "Area", "DOUBLE", "8", "1")
            areaExpression = "float(!SHAPE.AREA!)"
            arcpy.CalculateField_management(tempWS + "f7se.shp", "Area", areaExpression,
                                            "PYTHON3")
            
            #Add field with the PERIMETER values
            arcpy.AddField_management(tempWS + "f7se.shp", "Perimeter", "DOUBLE", "8", "1")
            perimeterExpression = "float(!SHAPE.LENGTH!)"
            arcpy.CalculateField_management(tempWS + "f7se.shp", "Perimeter", perimeterExpression, "PYTHON3")
            
            
            #Reshaping the delineated polygons
            arcpy.AddMessage("Reshaping the delineated polygons ...")
            
            if not arcpy.Exists(tempWS + "g8b_lyr"):
                arcpy.MakeFeatureLayer_management(tempWS + "f7se.shp", tempWS + "g8b_lyr")
            
            fields = ['MBG_Width']

            all_rows = [i[0] for i in arcpy.da.SearchCursor(tempWS + "g8b_lyr",fields[0])]
            min_width = min(all_rows)
            arcpy.AddMessage(min_width)
            
            if in_smoosimp == "No simplification or smoothing":
                arcpy.AddMessage("No simplification or smoothing")
                arcpy.CopyFeatures_management(tempWS + "g8b_lyr", tempWS + "h9z.shp")
            
            elif in_smoosimp == "Light simplification and smoothing":
                
                arcpy.AddMessage("Light simplification and smoothing...")
        
                # Simplifying
                arcpy.SimplifyPolygon_cartography(tempWS + "g8b_lyr", tempWS + "h9s.shp", "WEIGHTED_AREA", "%d" %float(min_width/10))
                arcpy.Delete_management(tempWS + "g8b_lyr")
                
                # Smoothing
                if not arcpy.Exists(tempWS + "h9s_lyr"):
                    arcpy.MakeFeatureLayer_management(tempWS + "h9s.shp", tempWS + "h9s_lyr")
                
                arcpy.SmoothPolygon_cartography(tempWS + "h9s_lyr", tempWS + "h9z.shp", 
                                                "PAEK", "%d Meters" %float(min_width/5), "NO_FIXED", "NO_CHECK")
        
                arcpy.Delete_management(tempWS + "h9s_lyr")
                
            elif in_smoosimp == "Heavy simplification and smoothing":
                
                arcpy.AddMessage("Heavy simplification and smoothing...")
        
                # Simplifying
                arcpy.SimplifyPolygon_cartography(tempWS + "g8b_lyr", tempWS + "h9s.shp", "WEIGHTED_AREA", "%d" %float(min_width/5))
                arcpy.Delete_management(tempWS + "g8b_lyr")
                
                # Smoothing
                if not arcpy.Exists(tempWS + "h9s_lyr"):
                    arcpy.MakeFeatureLayer_management(tempWS + "h9s.shp", tempWS + "h9s_lyr")
                
                arcpy.SmoothPolygon_cartography(tempWS + "h9s_lyr", tempWS + "h9z.shp", 
                                                "PAEK", "%d Meters" %float(min_width/2), "NO_FIXED", "NO_CHECK")
        
                arcpy.Delete_management(tempWS + "h9s_lyr")
            
           
            #Final operations
            arcpy.AddMessage('Final cleaning ...')    
            
            #Delete useless fields
            keepFields.extend(["Area","Perimeter","Relief","MBG_Width",
                               "MBG_Length","MBG_W_L","MBG_Orient"])
            helper.field_removal(tempWS + "h9z.shp", keepFields)
            
            #Save file
            arcpy.AddMessage('Saving official shapefile ...')
            arcpy.Copy_management(tempWS + "/" + "h9z.shp", workspace + "/" + outFeat)
            
            #Delete temporary files if requested
            if str(delTemp) == 'true':
                arcpy.AddMessage('Deleting temporary files ...')
                helper.cleaning(tempWS[:-1])
                
            #Print file report
            txtFile = open(workspace + "/" + str(outFeat[:-4]) + "_Info.txt", "w")
            txtFile.write("Script: CoMMa Delineation ToolBox v1.0" "\n")
            txtFile.write("\n")
            txtFile.write("File Name: " + str(outFeat) + "\n")
            txtFile.write("Input DEM: " + os.path.basename(inputDEM) + "\n")
            txtFile.write("Input Geomorphons layer: " + os.path.basename(inputGeom) + "\n")
                            
            txtFile.write("\n")
            txtFile.write("Core classes: " + str(co_gmf_class) + "\n")
            txtFile.write("Subordinate classes: " + str(su_gmf_class) + "\n")
            txtFile.write("\n")
            txtFile.write("Input cell size: " + str(Cellsize_X) + "\n")
            txtFile.write("Minimum relief value (m): " + minVR + "\n")
            txtFile.write("Minimum Width (m): " + minWidth + "\n")
            txtFile.write("Minimum W/L Ratio: " + minRatio + "\n")
            txtFile.write("Buffer extent (m): " + delBuffer + "\n")
            txtFile.write("Simplification or smoothing?: " + in_smoosimp + "\n")
            txtFile.write("\n")
    
            txtFile.close()
        
        return
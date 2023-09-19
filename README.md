# CoMMa Toolbox

## Overview 
Confined Morphologies Mapping (**CoMMa**) Toolbox is an **ArcGIS Pro** toolbox created for semi-automated morphology mapping.

It includes a selection of tools for the **delineation** and **description** of any type of enclosed features on a DEM, either negative or positive.
The CoMMa Toolbox is made up of individual Python scripts that use a sequence of pre-existing ArcGIS geoprocessing tools and do not require the installation of any new Python package.

Riccardo Arosio (University College Cork) and Joana Gafeira (British Geological Survey) conceived the original idea of the new ArcGIS Pro based on a previous toolbox created by Joana Gafeira, the BGS Seabed Mapping Toolbox (Gafeira, J., 2017). 
The tools development was mainly funded by INFOMAR through the Irish Marine Institute’s research grant PDOC 19/08/03. 
The British Geological Survey and iAtlantic have also supported the creation of the toolbox.

##  Documentation
The documentation for CoMMa is available at….. , including a user guide, example code, and gallery.

## Auxiliaries Tools
The CoMMa Toolbox works on DEMs, obtained from multibeam echosounder data or other geophysical and optical instruments (e.g., Lidar, 3D seismic etc.). 
Datasets may be affected by artefacts that can hinder a correct delineation of the features of interest, for example, vessel motion-related artefacts. 
A degree of data preparation, such as cleaning the initial data to remove artefacts, could be advised to enhance the performance of both the delineation and characterisation tools. 
The CoMMa Toolbox includes five tools devoted to data preparation, that can be found within the sub-toolbox CoMMa_Auxiliaries.

###Local Topographic Position (LTP) derivatives

|    Tool     |                   Description|
|-------------|:----------------------------------------------------------------------------------------------------------:|
| Mean LTPs   |   Local topographic position index metrics based on the absolute and relative mean of the neighbourhood.   |
| Median LTPs |   Local topographic position index metrics based on the absolute and relative median of the neighbourhood. |

## Delineation Tools
There are two available delineation tools in CoMMa Toolbox, the “Fill-based” and the “Geomorphons-based” delineation tools. 
The delineation tools can use directly the DEM or a derived raster (such as BPI or Geomorphons) and will be based on a set of user-defined thresholds to best delineate the target features.   

## Description Tools
The Description tools in CoMMa Toolbox calculate a series of basic geometrical and statistical attributes, additional metrics (such as zonal vector ruggedness and aspect variability index), backscatter statistics and the volume for each shape contained in the delineation shapefile. 
These tools can be used to characterise both features that were mapped automatically or manually mapped.

## Installation
The easiest way to ……:

## Copyright and licence
CoMMa Toolbox may be freely distributed, modified and used commercially under the terms of its GNU LGPLv3 license.

## Acknowledgements 

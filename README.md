
![image](https://github.com/ricarosio/CoMMa/assets/145455310/87c96df5-d4bd-4dcc-8eac-d837a9fb275e)

# Table of contents

- [Overview ](#overview)
- [Documentation](#documentation)
- [Data Preparation Tools](#data-preparation-tools)
  * [Local Topographic Position (LTP) derivatives](#local-topographic-position-(ltp)-derivatives)
- [Pre-processing](#pre-processing)
  * [Fencing Tool](#fencingtool)
- [Delineation Tools](#delineation-tools)
- [Description Tools](#description-tools)
- [Installation](#installation)
- [Copyright and licence](#Copyright-and-licence)
- [Acknowledgements](#acknowledgements)
    

# Overview 
Confined Morphologies Mapping (**CoMMa**) Toolbox is an **ArcGIS Pro** toolbox created for semi-automated morphology mapping.

It includes a selection of tools for the **delineation** and **description** of any type of enclosed features on a DEM, either negative or positive.
The CoMMa Toolbox is made up of individual Python scripts that use a sequence of pre-existing ArcGIS geoprocessing tools and do not require the installation of any new Python package.

Riccardo Arosio (University College Cork) and Joana Gafeira (British Geological Survey) conceived the original idea of the new ArcGIS Pro based on a previous toolbox created by Joana Gafeira, the BGS Seabed Mapping Toolbox (Gafeira, J., 2017). 
The tools development was mainly funded by INFOMAR through the Irish Marine Institute’s research grant PDOC 19/08/03. 
The British Geological Survey and iAtlantic have also supported the creation of the toolbox.

#  Documentation
The documentation for CoMMa is available at….. , including a user guide, example code, and gallery.

# Data Preparation Tools
The CoMMa Toolbox works on DEMs, obtained from multibeam echosounder data or other geophysical and optical instruments (e.g., Lidar, 3D seismic etc.). 
Datasets may be affected by artefacts that can hinder a correct delineation of the features of interest, for example, vessel motion-related artefacts. 
A degree of data preparation, such as cleaning the initial data to remove artefacts, could be advised to enhance the performance of both the delineation and characterisation tools. 
The CoMMa Toolbox includes five tools devoted to data preparation, that can be found within the sub-toolbox CoMMa_Auxiliaries.

## Local Topographic Position (LTP) derivatives

|    Tool     |                   Description|
|-------------|:----------------------------------------------------------------------------------------------------------:|
| Mean LTPs   |   Local topographic position index metrics based on the absolute and relative mean of the neighbourhood.   |
| Median LTPs |   Local topographic position index metrics based on the absolute and relative median of the neighbourhood. |


## Pre-processing

### Fencing Tool	
This tool creates an artificial containing “fence” around the perimeter of the input DEM, allowing the delineation of landforms that are at the boundary of the dataset and that otherwise would be considered, by the “Boundary-based Delineation” tool, as unconfined morphological features. 
This script assists the feature delineation tools by creating a buffer around the Input Raster. For target features with positive relief, the minimum value of the input DEM will be used for the artificial fence set by the buffer, whereas for negative target features, it will use the maximum value. This artificial fence should allow the delineation of landforms only partially captured within the dataset and that otherwise would be considered, by the “Boundary-based Delineation” tool, as unconfined. 

![image](https://github.com/ricarosio/CoMMa/assets/145455310/f726129e-fb85-42a5-804a-9b6241458c9e) 

# Delineation Tools
There are two available delineation tools in CoMMa Toolbox, the “Fill-based” and the “Geomorphons-based” delineation tools. 
The delineation tools can use directly the DEM or a derived raster (such as BPI or Geomorphons) and will be based on a set of user-defined thresholds to best delineate the target features.   

# Description Tools
The Description tools in CoMMa Toolbox calculate a series of basic geometrical and statistical attributes, additional metrics (such as zonal vector ruggedness and aspect variability index), backscatter statistics and the volume for each shape contained in the delineation shapefile. 
These tools can be used to characterise both features that were mapped automatically or manually mapped.

# Installation
The easiest way to ……:

# Copyright and licence
CoMMa Toolbox may be freely distributed, modified and used commercially under the terms of its GNU LGPLv3 license.

# Acknowledgements 

Riccardo Arosio (University College Cork) and Joana Gafeira (British Geological Survey) conceived the original idea of the new ArcGIS Pro based on a previous toolbox created by Joana Gafeira, the BGS Seabed Mapping Toolbox (Gafeira, J., 2017). Riccardo Arosio wrote the Python scripts while Joana Gafeira and Laurence De Clippele performed extensive testing.

The tools development was mainly funded by INFOMAR through the Irish Marine Institute’s research grant PDOC 19/08/03. The British Geological Survey and iAtlantic have also supported the creation of the toolbox.


![image](https://github.com/ricarosio/CoMMa/assets/145455310/5fdd2e27-a0cd-4895-82c1-f5eebc690f1c)       ![image](https://github.com/ricarosio/CoMMa/assets/145455310/84f767ec-1421-4513-870d-3b502776e568)        ![image](https://github.com/ricarosio/CoMMa/assets/145455310/ecdc36f2-b59e-409f-92f9-09439b80b6d2)  
![image](https://github.com/ricarosio/CoMMa/assets/145455310/f8b95095-86d1-40b5-9c4b-8e2714579ef6)         ![image](https://github.com/ricarosio/CoMMa/assets/145455310/a30e5c28-1591-47ac-9224-32018afd9076)  




[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11109024.svg)](https://doi.org/10.5281/zenodo.11109024)

<p align="center" width="100%">
    <img width="50%" src="https://github.com/ricarosio/CoMMa/blob/main/CoMMa_Toolbox/CoMMa_Version_1_1/Documentation/Images/CoMMaToolbox.png">   
</p>
<p align="center" width="100%">
**Version 1.2**
</p>

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
- [Copyright and licence](#copyright-and-licence)
- [Acknowledgements](#acknowledgements)
    

# Overview 
**Co**nfined **M**orphologies **Ma**pping (**CoMMa**) Toolbox is an **ArcGIS Pro** toolbox created for semi-automated morphology mapping.

It includes a selection of tools for the **delineation** and **description** of any type of enclosed features on a DEM, either negative or positive. The CoMMA toolbox includes three toolsets that allow the user to pre-process DEM data and calculate local topographic parameters, delineate potential features using the Delineation tools and describe the morphological characteristics of confined features using the basic, texture and volume descriptors tools. 
The CoMMa Toolbox is made up of individual Python scripts that use a sequence of pre-existing ArcGIS geoprocessing tools and do not require the installation of any new Python package.

Riccardo Arosio (University College Cork) and Joana Gafeira (Kelpie Geoscience) conceived the original idea of the new ArcGIS Pro based on a previous toolbox created by Joana Gafeira, the BGS Seabed Mapping Toolbox (Gafeira, J., 2017). 
The tools development was mainly funded by **INFOMAR** through the **Irish Marine Institute**’s research grant PDOC 19/08/03. 
The **British Geological Survey** and **iAtlantic Project** have also supported the creation of the toolbox.

#  Documentation
Several support material for the CoMMa Toolbox is available within the documentation folder, including CoMMa Toolbox summary flyer, tools metadata, and [**User Guide**](CoMMa_Toolbox/CoMMa_Version_1_2/Documentation_v1_2/Arosioetal_CoMMa_Supp_CoMMa_User_guide_V1_2.pdf). 
The paper [CoMMa: A GIS geomorphometry toolbox to map and measure confined landforms](https://www.sciencedirect.com/science/article/pii/S0169555X24001776?via%3Dihub) also provides an assessment of the toolbox performance and additional information. 

# Data Preparation Tools
The CoMMa Toolbox works on DEMs, obtained from multibeam echosounder data or other geophysical and optical instruments (e.g., Lidar, 3D seismic etc.). 
Datasets may be affected by artefacts that can hinder a correct delineation of the features of interest, for example, vessel motion-related artefacts. 
A degree of data preparation, such as cleaning the initial data to remove artefacts, could be advised to enhance the performance of both the delineation and characterisation tools. 
The CoMMa Toolbox includes five tools devoted to data preparation, that can be found within the sub-toolbox CoMMa’s Data Preparation.

## Local Topographic Position (LTP) derivatives

-  **Mean LTPs:**   Local topographic position index metrics based on the absolute and relative mean of the neighbourhood.
    - The Bathymetry Position Index (BPI)
    - The Deviation from mean elevation (DEV)
    
-  **Median LTPs:**     Local topographic position index metrics based on the absolute and relative median of the neighbourhood.
    - Median Bathymetry Position Index 
    - Minimum Median Bathymetry Position Index
    - Maximum Median Bathymetry Position Index
    - Directional Median Bathymetry Position Index  

## Pre-processing tools

### Fencing Tool	
This tool creates an artificial containing “fence” around the perimeter of the input DEM, allowing the delineation of landforms that are at the boundary of the dataset and that otherwise would be considered, by the “Boundary-based Delineation” tool, as unconfined morphological features. 
This script assists the feature delineation tools by creating a buffer around the Input Raster. For target features with positive relief, the minimum value of the input DEM will be used for the artificial fence set by the buffer, whereas for negative target features, it will use the maximum value. This artificial fence should allow the delineation of landforms only partially captured within the dataset and that otherwise would be considered, by the “Boundary-based Delineation” tool, as unconfined. 
<p align="center" width="100%">
    <img width="40%" src="https://github.com/ricarosio/CoMMa/blob/main/CoMMa_Toolbox/CoMMa_Version_1_1/Documentation/Images/FencingTool.jpg">
</p>

### Filter and clip
This tool removes the flat or featureless areas in the DEM and preserves areas of the seabed where the features are more likely to occur. The application of this tool is particularly useful to remove the effects of broad-scale topography on Fill-based delineations. 

### Smoothing filters
The Smoothing filters tool provides a selection of standard filtering algorithms that can be used smooth the DEM and remove, at least partially, noise and artefacts. Although they will also subdue the real signal proportionally to their respective aggressiveness. These filters are: 
- MEAN — Calculates the mean (average value) of the cells in the neighbourhood set by the user.
- MEDIAN — Calculates the median of the cells in the neighbourhood set by the user.
- LOW-PASS 3x3 — It calculates the average value for each 3 x 3 neighbourhood.
- LOW-PASS 5x5 — It calculates the average value for each 5 x 5 neighbourhood.


# Delineation Tools
There are two available delineation tools in CoMMa, the **“Boundary-based delineation”** and the **“Elements-based delineation”** tools. 
The **“Boundary-based delineation”** tool focuses on recognising the landform boundary and can be applied directly to bathymetric DEM and, when the general seascape is otherwise flat, it might be sufficient to isolate and correctly delineate the targeted features. However, this is often not the case, and the interference of sloping topography or other underlying large-scale landforms can distort the signal of the targets. In this case, a LTP derivative can be use to isolate a specific wavelength thought to best delineate the feature of interest. 
The **“Elements-based delineation”** relies on land surface units created by the geomorphon landforms algorithm. This tool delineates confined landforms by aggregating positive or negative land surface elements (geomorphons).

<p align="center" width="100%">
    <img width="65%" src="https://github.com/ricarosio/CoMMa/blob/main/CoMMa_Toolbox/CoMMa_Version_1_1/Documentation/Images/DelineationExamples.jpg"> 
</p>    

# Description Tools
The Description tools in CoMMa Toolbox calculate a series of basic geometrical and statistical attributes, additional metrics (such as zonal vector ruggedness and aspect variability index), backscatter statistics and the volume for each shape contained in the delineation shapefile. 
These tools can be used to characterise both features that were mapped automatically or manually mapped.

-	**Basic Descriptors**: Calculates a series of geometrical and statistical attributes for each shape contained in the delineation shapefile.
-	**Texture Descriptors**: Calculates a few additional metrics, such as zonal vector ruggedness and aspect variability index, and optionally backscatter statistics.
-	**Volume Descriptor**: Calculates the volume and more accurately the height for each shape contained in the delineation shapefile.


# Installation
The CoMMa toolbox is comprised of three Python toolboxes, that can be loaded to an ArcGIS Pro project as any standard toolbox. 

To add the three CoMMa python toolboxes to a project, in the Catalog pane, right-click on Toolboxes and select “Add Toolbox”
Then, navigate to the folder where CoMMa Toolbox was saved and select the following .pyt files: 
<p align="center" width="100%">
- CoMMa Data Preparation.pyt,  CoMMa Delineation.pyt  and  CoMMa Description.pyt.
</p>
<p align="center" width="100%">
    <img width="70%" src="https://github.com/ricarosio/CoMMa/blob/main/CoMMa_Toolbox/CoMMa_Version_1_1/Documentation/Images/AddToolbox_L.gif">
</p>
A reference to the toolsets is saved within the project and they will be in the Toolbox node of the Catalog pane the next time the project is open.

# Copyright and licence
CoMMa Toolbox may be freely distributed, modified and used commercially under the terms of its GNU LGPLv3 license.

# Acknowledgements 

Riccardo Arosio (University College Cork) and Joana Gafeira (Kelpie Geoscience) conceived the original idea of the new ArcGIS Pro based on a previous toolbox created by Joana Gafeira, the BGS Seabed Mapping Toolbox (Gafeira, J., 2017). Riccardo Arosio wrote the Python scripts while Joana Gafeira and Laurence De Clippele performed extensive testing.

The tools development was mainly funded by INFOMAR through the Irish Marine Institute’s research grant PDOC 19/08/03. The British Geological Survey and EU H2020/iAtlantic have also supported the creation of the toolbox.


![image](https://github.com/ricarosio/CoMMa/assets/145455310/5fdd2e27-a0cd-4895-82c1-f5eebc690f1c)       ![image](https://github.com/ricarosio/CoMMa/assets/145455310/84f767ec-1421-4513-870d-3b502776e568)        ![image](https://github.com/ricarosio/CoMMa/assets/145455310/ecdc36f2-b59e-409f-92f9-09439b80b6d2)  
![image](https://github.com/ricarosio/CoMMa/assets/145455310/f8b95095-86d1-40b5-9c4b-8e2714579ef6)         ![image](https://github.com/ricarosio/CoMMa/assets/145455310/a30e5c28-1591-47ac-9224-32018afd9076)  




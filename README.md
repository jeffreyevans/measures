# measures
Code in support of the MACP indigenous measures work

**Measuring Our Impacts – Workflow**

This repository is in support of The Nature Conservancy’s Measuring Our Impact effort from the Chief Strategy Office. There are limited ways established to track metrics of ecological condition through time and those that are used are limited geographically or to specific time frames. This work is an effort to establish a workflow for creating a hierarchical modeling framework for understanding ecosystem characteristics. At coarse scales (e.g., ~500 m) some albeit limited inference can be drawn, but finer spatial-scale data and ancillary data related to local threats and actual desired conditions can be used to draw more robust inference.

**Workflow for Ecological Condition**

1.  **Data acquisition **

    1.  Download MODIS/Copernicus current 500m multispectral

        1.  Only a 2 bands are available at 250 m; good summary of availability is found [here](https://ladsweb.modaps.eosdis.nasa.gov/missions-and-measurements/modis/)

    2.  Download MODIS/Copernicus timeseries indices (LAI, FPAR, & Maybe VCI)

    3.  Download MODIS Landcover

    4.  Download climate metrics (30yr norms), [Rehfeldt data](http://charcoal.cnre.vt.edu/climate/) for North America; else WorldClim. 

    5.  Download 90m SRTM Digital Elevation Model (DEM)

        1.  Note – elevatr package might do this; Jeff Hollister from USEPA has added functionality to download SRTM specifically - devtools::install\_github("jhollist/elevatr", "opentopo-gl3")

* Note - For US, data reprojected to USGS Albers Equal Area ( conus.aea = "+proj=aea +lat\_1=29.5 +lat\_2=45.5 +lat\_0=37.5 +lon\_0=-96 +x\_0=0 +y\_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no\_defs"); For now, in other areas we’ll use the a global Albers Equal Area projection (global.aea = "+proj=aea +lat\_1=29.5 +lat\_2=42.5")

2.  **Data processing **

    1.  QC Band Masking

    2.  NoData (NA) imputation and timeseries smoothing of spectral temporal indices

3.  **Creation of *x* parameter(s)**

    1.  Calculate topographic metrics (e.g., CTI, slope/aspect interaction) using RSAGA

    2.  If necessary, derive climate norms and other metrics (NA for North America)

    3.  Run PCA’s on topographic and climate variables, respectively, to create rasters representing collapsed variation of climate and abiotic variables.

4.  **Stratification (training sample)**

    1.  Create classification(s) rasters of stratifying variables (e.g., using quartiles) representing a temperature/moisture and biophysical condition gradient.

    2.  Create a single stratification raster representing combinatorics of classified stratifying variables.

    3.  Draw a spatial stratified random sample representing *p* of the population (ie., number of pixels).

5.  **Creation of *y* parameter and model matrix**

    1.  Calculate weighted Kappa on landcover timeseries (w/in n x n window) 

    2.  Calculate Kendall’s Tau on timeseries of spectral index (standardize slope ranges to observed data through whole raster, \[-1- 1\]) 

    3.  Calculate max of spectral index for current growing season or use a long-term statistic such as Vegetation Condition Index (VCI). (TBD)

    4.  Standardize the above dependent response variables to common range and expected response 

    5.  Create function (e.g., linear) to combine above dependent variables into a single response (*y*). We can also try a multiple imputation, using separate dependent variables.

    6.  Assign values from the y and x variables to the spatial stratified random sample to create the model matrix used in the imputation model.

6.  **Imputation Model**

    1.  Specify imputation model using random forest proximities for the multivariate distances

        1.  Apply back-prediction evaluation and validation (eg., RMSE) of model

        2.  Compare to other distance measures (eg., Most Similar Neighbor, Gradient Nearest Neighbor, Mahalanobis).

    2.  Predict/impute model to landscape-level independent variables (rasters) creating a matrix of *k*=1 and of *k*=(n-1) specifically pulling distances of each *k*.

    3.  Run additional models focused on specific ranges of our *y* variable indicating ecological condition, this follows the idea of using “reference conditions”.

7.  **Results summary and visualization**

    1.  Create maps of distances based on reference condition maps

        1.  Calculate proportion of pixels, within threshold distance, of given reference condition(s) across landscape.

        2.  Calculate proportion of pixels, within threshold distance, of given reference condition(s) within TNC interventions and control locations
        
**Workflow for Intervention Analysis**

1.  **Sampling**

    1.  Utilizing Dynamic Time Warping (DTW) to draw random sample of temporal evaluation metric

    2.  Evaluate land tenure status of stratified random samples

    3.  Specify a Causal Impact Analysis (Brodersen et al., 2015) model for evaluating interventions 


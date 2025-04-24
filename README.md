# measures
Code in support of the retrospective PFP analysiis

Contact:
Jeffrey Evans, PhD., 
The Nature Conservancy
Global Protect Science
jeffrey_evans@tnc.org

**Measuring Our Impacts – Modeling Pipeline**

This repository is in support of The Nature Conservancy’s retrospecive evaluation of 6 PFP engagment. There are limited ways to track metrics of ecological condition through time and those that are used are limited geographically or to specific time frames. This work is an effort to establish a workflow for creating a hierarchical modeling framework for understanding ecosystem characteristics. At coarse scales (e.g., ~250 m) some albeit limited inference can be drawn, but finer spatial-scale data and ancillary data related to local threats and actual desired conditions can be used to draw more robust inference.

**Workflow for evaluation of PFP outcomes using Ecological Condition metrics**

1.  **Data acquisition and processing**

    1.  Download and process Lear Area Index (LAI) and Fractional Cover (fCOV) data from the Univerity of Hong Kong [GLASS repostory](https://glass.hku.hk/). The "download_GLASS.R" script facilitates download, in the nested HDF archive format, and processing of this data. 
	
    2.  Download [WorldClime data](http://charcoal.cnre.vt.edu/climate/) ()

    3.  Download IUCN Protected Area polygon vector data, merge with provided PFP polygon data (see build_knn_covariates.R). Note; there is an option to update these boundaries in the BuildTreatments.R script. 

    4.  Download Open Street Map polygon, line and point data for settletments, roads, rivers, and lakes and calculate Eculidian distance rasters (see build_knn_covariates.R) 
   
    5. 	Download and process the [World Terrestrial Ecosystems](https://www.usgs.gov/centers/geosciences-and-environmental-change-science-center/science/global-ecosystems-global-data) raster, subset and reproject to each study area

    6.  Perform pre-processing on temporal data including; NoData (NA) imputation, outlier removal using LOESS smooting, and detrending periodicity/seasionality using BEAST model (detrend_timeseries_BEAST.R)   

	7.  Create a forest mask, using LAI timeseries, by identifying pixels exceeding defined forest threshold(s) that have a 10% (n=103) contigious run during any point in the timeseries (see create_forest_mask.R) 

2.  **Select kNN controls**

Note; these processing steps are implementing in the BuildTreatments.R script. It is currently written with multi-threading and database output. 

    1.  Assign an "in or out" attribute to all forest mask pixels based on intervention units. Split pixels/point centroids that are within interventions to a interventions dataset and outside of interventions to control dataset

    2.	Assign ecoosystems (World Terrestrial Ecosystems), climate PCA and ecological indicators raster values to both control and interventions datasets  	 

    3.  For control selection, start by stratifying based on the ecosystem type 

	4.  Using a Mahalanobis kNN, calculate multivarite distances between all interventions/controls observsations and select n=100 smallest distances in the control data, for each intervntion observsations

	5.  For intervention observsation, calculate the pre-intervention temporal distance (using Dynamic Time Warpping) for each n=100 candiate controls. Select one control observsation with the smallest temporal distance.  	

3.  **Creation of *y* parameters and model design matrix**

Note; these processing steps are implementing in the BuildTreatments.R script.

    1.  For "current" paramter, calcualte the expected median LAI and fCOV for the growing seasions (or non-raniy season) representing the last year of data in the timeseries. 

    2.  For "change" paramter, calcualte the expected median LAI and fCOV for the growing seasions (or non-raniy season) representing the first year and last year of data in the timeseries. The change is the delta of current and baseline.  

    3.  For "volume" paramter, calcualte the area under the curve for the entire timeseries. 

    4.  Feature engineering of the timeseries for the model design matrix (covariates)	

4.  **Trend Analysis**

Note; these processing steps are implementing in the tau_trend_analysis.R script are is written to be multi-threaded

    1.  Subset timeseries data to pre and post interventions and growing seasons within each year for LAI and fCOV. 

    2.  Calculate the Kendall Tau, for each pixel timeseries, for the pre/post intervention timeseries 

5.  **Causal Model**

Note; the causal model is implemented in the CausalForestModel.R script

    1.  Stratify models by intervention unit and ecosystem type so, each model represents the ecosystem type within a given intervention unit. The esablished minimum sample size for a model is n=100 

    2.  Using Causal Forest (ref) we estimate treatment effect size for each intervention observsation, standardized to an effect size. Estimates are written to corresponding pixel to create an effect size  raster output. 


6.  **Results summary and visualization**

Note; these processing steps are implementing in the Compile_Classify_Results.R script

    1.  Tabulate results for effect sizes
	
		1. Classify effect size results into 0.1 bins 
		
		2. Output tables of classified effect size aggregating by;
		
			1. Interventions
			
			2. IUCN protection for 1a and 1b

	2.	Tabulate results for Tau trends

		1. Classify pre and post intervention Tau trends results into 0.1 bins 
		
		2.  Classify direction of change in pre and post intervention trends

		3. Output tables of classified trends aggregating by;
		
			1. Interventions
			
			2. IUCN protection for 1a and 1b
			
		4. Output tables of trend direction aggregating by;
		
			1. Interventions
			
			2. IUCN protection for 1a and 1b			

	3. Create rasters (in single tif) representing response variables for treatment and control pixels.

	4. Evaluate representation of ecosystems (WTE's) of interventions compared to entire country, output table (see process_elu.R & representation.R). 

	5. Evaluate trend and effect size relationships with climate. The CHELSA_climate_analysis.R and climate_correlations.R scripts support a thorough analysis of climate change metrics, including download, processing and derivation of climate metrics. The relationships are characterized using a nonlinear regression. The ttopographic_analysis.R script implements a more general univariate assesment with an option for elevation, slope*cos(aspect) (slope aspect interaction transformation), trasp (aspect transformation) or the aridity index.  
  	
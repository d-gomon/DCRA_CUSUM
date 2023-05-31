# R_code
Code used to analyze the [DCRA](https://dica.nl/dcra/home) data set.


## The six steps of the application.

- **Step 1. Data preparation** 

Removing patients from the analysis, converting to the correct format. This process is described in the file `1_LoadingData.R`.

- **Step 2. Calculating control limits for the CUSUM charts**

To determine whether a CUSUM chart has produced a signal, we first need to calculate control limits for the considered hospitals. This process is described in `2_DCRA_control_limit_calculation.R`.

- **Step 3. Constructing the CUSUM charts + funnel plot**

The CUSUM control charts and funnel plot are constructed on the data in `3_DCRA_charts_calculation.Rmd`.

- **Step 4. Anonimizing hospital numbers**

As the data used for this analysis contains sensitive information, we anonimize hospital identifiers in this file. Some of the code has been removed from the script `4_Anonimize_Charts_Public.Rmd`

- **Step 5. Determining hospital detections and graphical display**

In this step, we determine which hospitals have been detected using which methods and summarize the information graphically in `5_DCRA_Detection_anonimised.Rmd`

- **Step 6. Producing additional materials**

In the file `6_AdditionalFiles.Rmd`, we produce some additional materials for readers of the article. These can be found in the subdirectory `Additional Materials`.


## Additional Materials

The four resulting files from Step 6 are attached in the subdirectory `Additional Materials`.


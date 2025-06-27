# Data and code accompanying "Should you use data integration for your distribution model?"


## Simulation

The simulation code is organized as follows. 

The folder "`main_code/`" contains files to execute each of the three 
simulations (described in Appendix A) called `run_sim1.R`, `run_sim2.R`, and 
`run_sim3.R`. On execution, each of these will produce outputs and save them
to the folder `sim_output/`. The code files `vis_sim1.R`, `vis_sim2.R`, and
`vis_sim3.R` are used to produce figures visualizing the outcomes of the
simulation. There are two other relevant scripts, `run_apriori.R` and
`vis_aposteriori.R`, which produce *a priori* and *a posteriori* results
described in Appendixes 2 and 3.

The folder `support_fn/` contains files that are called my scripts in `main_code`.
These files store the functions used to simulate and estimate the different
models.


## Worked example

The code to reproduce the Worked Example is organized in the folder `worked_example`. 

The important scripts are numbered: `1_data_summaries.R`, 
`2_fit_model_theta1fixed.R`, and `3_fit_model_NDVI_in_det.R`. These scripts 
create *a priori* data summaries, execute the first model, and execute the modified
model with NDVI affecting detection, respectively. These scripts depend on the 
helper script `preprocess_data.R` which loads inputs from four data files,
`WE1_iNat.csv`,
`WE1_sequences.csv`,
`WE1_deployments.csv`, and
`Avg_NDVI_Landsat_NYC.tif`, which provide the iNaturalist observations, SNAPSHOT
sequences, SNAPSHOT deployment metadata, and NDVI raster used in the analysis.





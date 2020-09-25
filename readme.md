#  Code for Tracing the influence of land-use change on water quality and coral reefs using a Bayesian model

2020-09-25

Citation:

[Brown CJ, Jupiter SD, Albert S, Klein CJ, Mangubhai S, Maina JM, Mumby P, Olley J, Stewart-Koster B, Tulloch V, Wenger A. Tracing the influence of land-use change on water quality and coral reefs using a Bayesian model. Scientific reports. 2017 Jul 6;7(1):1-0.](https://www.nature.com/articles/s41598-017-05031-7)

Feel free to email me if you have queries about running these models (chris.brown@griffith.edu.au). The model was created in 2015, so R has progressed considerably since then. If I was to do this again, I'd use updated spatial packages (updates to raster and pkg sf, instead of sp), and probably write the model with STAN which is much faster at fitting models than JAGS.

## Brief description of files and folders

/wq_bayes_testsim has scripts for running the simulation tests described in the paper

The core of the Bayesian model (the JAGS code) is in BayesianTracingModel.txt

B3_Turbidity_model_power_hierarchical_2alpha.R is the script I used to run the main model for Fiji presented in the paper (includes model pre and post processing). The three datafiles for this script are also here, as .csv files.

B4_model_predictions_power_hierarchical_2alpha.R demonstrates some of the post-processing and how to do model predictions (again this could be updated )

There are other steps before script B3 that relate to downloading the MERIS turbidity images and processing those to average over time and also to resample them to a lower resolution (JAGS slows down drastically the larger your dataset).
I haven't shared that code here, as it will be specific to each project, I haven't documented it very well and the R spatial methods used are not up to date with the latest spatial packages.

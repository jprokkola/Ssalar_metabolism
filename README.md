#README
The folder contains R scripts for data analysis associated with Prokkola et al. 2021 (Genetic coupling of life-history and aerobic performance in juvenile Atlantic salmon). No guarantees that these can be used for any other purpose, but feel free to try. 

File names shown like `this`.

To run the analysis, download the Data folder from Zenodo (see data availability in the manuscript), and place it in this folder.

## Folder: SMR analysis
Analysing standard metabolic rates. Data has been corrected for background respiration and flush phases excluded prior to this analysis.


> `Functions_for_filtering_SMR.R`
> 
> Functions to identify linear slopes using a combination of quadratic slopes, linear regression residuals, and R^2 -value filtering.
> 
>  
> `Plot_functions_SMR_analysis.R`
> 
> Functions for plots to use in multi-batch respirometry analysis (here multiple batches of 16 chambers). Use together with the filtering functions.
> 
> 
> `Highfood_SMR_analysis.R` and `Lowfood_SMR_analysis.R`
> 
> Calculating SMR from background-corrected oxygen saturation data, high food and low food treatments in separate scripts. These use functions specified in the two files above. Output is the SMR for each individual combined with other sample information from All families info -file (see data availability).
> 
> `Combine_SMR.R`
> 
>Combining High food and Low food SMR data. Compare SMR from quantile and MLND approaches, and get N per group.
>
>`SMR_models_revised.R`
>
> Statistical analysis of SMR with linear mixed models and model averaging. Calculating mass-scaling exponents and predicted means for plotting.


## Folder: MMR analysis

Analysing maximum metabolic rate (MMR) with two approaches: the "spline method", and 1min sliding windows using the R package respR. Data has been corrected for background respiration prior to this analysis.

Sub-folders High food and Low food for fish different treatments. Each has scripts 1-3:

> `1.` `MMR_dataPrep_TAmethod_highfood.R` or `MMR_dataPrep_TAmethod_lowfood.R`
> 
> Data preparation for calculating MMR (identifying the point of beginning of slope and excluding poor slopes).
> 
> `2.` `MMR_spline_method_Fishresp_highfood.R` or `MMR_spline_method_Fishresp_lowfood.R`
> 
> Calculating MMR using the spline method from data prepared in TAmethod script.

>  `3.` `MMR_respr1min_highfood.R` or `MMR_respr1min_lowfood.R`
> 
> Calculating MMR from 1min sliding windows using respR and FishResp from data prepared in TAmethod script. 


> `4.` `Combine_MMR.R`
>
>Combining High food and Low food MMR data. 
>
>`5.` `MMR_models_revised.R.`
>
>Statistical analysis with linear mixed models and model averaging. Calculating predicted means for plotting.


## Folder: Aerobic scope analysis

> `AS_dataprep.R`
> 
> Combining SMR and MMR data, and calculating absolute aerobic scope and mass- and family-residuals for each trait.
> 
> `AbsAS_models_revised.R`
> 
> Statistical analysis with linear mixed models and model averaging. Calculating predicted means for plotting.


## Folder: Plot all traits
> `Fig2_script.R`
> 
Compilation of manuscript Fig. 2 from tables of predicted means.




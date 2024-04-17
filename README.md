# Bayesian Heirarchical stock recruitment model for Scottish salmon stocks

This model pools salmon egg-to-egg stock recruitment data for 11 Scottish stocks and 1 English stock
to allow the transportation of biological reference points (e.g. spawning stock at maximum sustainable yield)
to non-monitored stocks. Informative explanatory covariates are identified through cross-validation.

## Publications

Extending Bayesian hierarchical stock recruitment models for setting conservation limits for Atlantic salmon stocks.

James Ounsley, Nora N. Hanson, Gordon W. Smith, Jonathan P. Gillson and Stuart J. Middlemas

In prep.

## Dependencies

R libraries:

  - optparse
  - tidyverse
  - rstan
  - bayesplot

## Setup

1. TODO: Data - how are we providing this to the user?
2. Run through scripts in sequence. 

    + 01_wrangle_data.R - wrangles the raw data from use by the models.
    + 02_run_model.R - runs a model with specified covariates. Easiest to
     run from the command line using Rscript.exe passing in arguments through 
      optparse. E.g Linux:
    
          cd path/to/repo
          Rscript R/02_run_model.R -c "lat pc1"
      
        E.g. Windows, from cmd prompt:
    
          cd path\to\repo
          path\to\R\Rscript.exe R\02_run_model.R -c "lat pc1"
    
        Run all required models here before moving to next script.
      
    + 03_cross_validate_model.R - take a fitted model from 02 and perform leave one
      river out cross validation. Again, best run through Rscript.exe.
      Run all required models before moving to next script.
    + 04_create_cv_summary.R - run through all cv fits and derive summary values
      and plots.
    + 05_create_prediction_summary.R - Create posterior predictions from a 
      specified model, hard coded to use the best fit model in cv analysis.
      Modify if a different model is required, not all plots may be relevant
      however.
    + 06_create_waic_summary.R - calculate PSIS loo for specified models. 
      Defaults to those with only spatial smoother as used in the ms.
      
      

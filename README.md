# SurvTMLE
Targeted maximum likelihood estimation of the effect of an exposure on a time-to-event outcome

This repository contains supplementary material for the article "Guidelines and best practices for the use of targeted maximum likelihood and machine learning when estimating causal effects of exposures on time-to-event outcomes" by Denis Talbot, Awa Diop, Miceline Mésidor, Yohann Chiu, Caroline Sirois, Andrew J. Spieker, Antoine Pariente, Pernelle Noize, Marc Simard, Miguel Angel Luque-Fernandez, Michael Schomaker, Kenji Fujita, Danijela Gnjidic and Mireille E Schnitzer. 

The files included are as follows:

- The survTMLE_v*.R file contains an R function implementing targeted maximum likelihood estimation for estimating the effect of a single exposure on a time-to-event outcome as counterfactual survival curves and a working marginal structural model.  
- The Examples.R file contains a few examples of usage of the survTMLE function based on simulated data.
- The TutorialBoxes.R file contains the complete R code for a "handmade" implementation of TMLE, as presented in the boxes of our associated tutorial. 

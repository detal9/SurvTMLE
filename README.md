# SurvTMLE
Targeted maximum likelihood estimation of the effect of an exposure on a time-to-event outcome

This repository contains supplementary material for the article "Guidelines and best practices for the use of targeted maximum likelihood and machine learning when estimating causal effects of exposures on time-to-event outcomes" by Denis Talbot, Awa Diop, Miceline Mésidor, Yohann Chiu, Caroline Sirois, Andrew J. Spieker, Antoine Pariente, Pernelle Noize, Marc Simard, Miguel Angel Luque-Fernandez, Michael Schomaker, Kenji Fujita, Danijela Gnjidic and Mireille E Schnitzer. 

The files included are:

- survTMLE.R contains an R function implementing targeted maximum likelihood estimation for estimating the effect of a single exposure on a time-to-event outcome as counterfactual survival curves, average treatment effect at different time points, and a working marginal structural model.  
- Examples.R contains a few examples of usage of the survTMLE function based on simulated data.
- TutorialBoxes.R contains the complete R code for a "handmade" implementation of TMLE, as presented in the boxes of our associated tutorial.
- ComparisonLTMLE.R uses both the ltmle package and survTMLE in Example 1, illustrating some similarities and differences

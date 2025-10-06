############################################################
#                                                          #
#  Survival LTMLE with a time-fixed exposure and           #
#  time-varying censoring                                  #
#  V1.0.1                                                  #
#  date: 2025-10-06                                        #
#  Maintainer: Denis Talbot                                #
#              denis.talbot@fmed.ulaval.ca                 #
#                                                          #
############################################################

require(mlogit);
require(polspline);
require(SuperLearner);


#### Description
# Estimates counterfactual survival curves, average treatment effects and parameters
# of marginal structural hazard model using a targeted maximum likelihood estimator
# for time-to-event data subject to right censoring. 
#
# The function currently supports binary or multilevel exposures. The assumed data
# structure is (L1, A, C1, Y1, L2, C2, Y2, ..., LK, CK, YK), where L are covariates or
# potential confounders, A is the time-fixed exposure, C are indicators of censoring or
# loss to follow-up and Y are event indicators. As such, continuous follow-up time need
# to be discretized before using the function.
#
# The function does not currently support missing data other than those due to censoring.
# Please see our tutorial and examples for more information on how to use this function.


#### Arguments
# dat:           A dataframe containing the data to be used. 
#                The data must be supplied in a wide format (one row per individual).
#
# Yvar:          Either a character vector of the names of the Y variables or
#                a numeric vector of the indices of the Y variables in dat.
#                The Y variable must be coded 0/1, where 1 indicates that the event
#                has occured. Once Y_t = 1 for a given subject, the data for times
#                k > t are not used for that subject. The variables must be temporally
#                ordered in Yvar.
#
# Cvar:          Either a character vector of the names of the C variables or
#                a numeric vector of the indices of the C variables in dat.
#                The C variable must be coded 0/1, where 1 indicates that the follow-up
#                of the subject has been censored. Once C_t = 1 for a given subject,
#                the data for times k > t are not used for that subject. The variables
#                must be temporally ordered in Cvar.
#
# Avar:          Either the name of the exposure variable or its index in dat.
#                The exposure can either be a binary or a categorical
#                variable. Continuous exposures are not supported. If the exposure has more
#                than two levels, it must be of type factor.
#
# Lvar:          Either a list of character vectors indicating the names of the time-varying
#                L covariates measured at each time-point or a list of numeric vectors of the
#                indices of the L covariates. The sets of variables (the elements of the list)
#                must be temporally ordered in Lvar. For example, if there are two time points
#                and (L11, L12) are the covariates at time 1 and (L21, L22, L23) are the covariates
#                at time 2, then Lvar = list(c(L11, L12), c(L21, L22, L23)). If there are only
#                baseline covariates, it is possible to repeat the same covariates at all times
#                with lookback = 1.
#
# L0var:         A character vector of time-fixed L variables measured at baseline or
#                a numeric vector of the indices of the time-fixed L variables in dat (optional).
#
# lookback:      A numeric value indicating how far back in time should time-varying covariates
#                be considered when modeling the outcome indicators (Y) and the censoring
#                indicators (C). For example, if lookback = 2, Y_t and C_t will be modeled as
#                a function of L_t and L_{t-1}. The default is NULL, indicating that all previous
#                time-varying covariates should be used.
#
# lookbackY:     See description for lookback. LookbackY indicates the outcome-specific lookback.
#                The default is to use the value provided in lookback.
#
# lookbackC:     See description for lookback. LookbackC indicates the censoring-specific lookback.
#                The default is to use the value provided in lookback.
#
# Ymod:          The approach to be used to model the outcome indicators, either "parametric" or
#                "SL" (Super Learner). The default is "parametric". 
#                The choice "parametric" results in using logistic regressions.
#
# Cmod:          The approach to be used to model the censoring indicators, either "parametric" or
#                "SL" (Super Learner). The default is "parametric".
#                The choice "parametric" results in using logistic regressions.
#
# Amod:          The approach to be used to model the exposure, either "parametric" or
#                "SL" (Super Learner). The default is "parametric".
#                The choice "parametric" results in using a logistic regression if the exposure
#                is binary and a multinomial regression otherwise. 
#                If the exposure is categorical (>2 levels), Amod = "SL" results in using a
#                polychotomous regression and classification.
#
# SL.library:    A character vector or a list of character vectors indicating the learners and screeners
#                to be used within the Super Learner.
#                The default is c("SL.glm", "SL.glm.interaction"); The list of all default
#                learners can be viewed with listWrappers(); 
#
# Y.SL.library:  A character vector indicating the learners to be used within the Super Learner
#                for modeling the outcome. The default is the same as SL.library.
#
# C.SL.library:  A character vector indicating the learners to be used within the Super Learner
#                for modeling the outcome. The default is the same as SL.library.
#
# A.SL.library:  A character vector indicating the learners to be used within the Super Learner
#                for modeling the exposure. The default is the same as SL.library if A is binary.
#                If A is multilevel (>2 levels) then this option is ignored and
#                polychotomous regression and classification is used.
#
# gbound:        Bound for the g estimates (A and C). 
#                For exposure A, all predicted probabilities < gbound or > 1 - gbound
#                are truncated. For censoring C, only predicted probabilities 
#                P(C = 0|A, L) < gbound are truncated. Defaults is 0.005.
#
# V:             Number of folds for the cross-validation when using the Super Learner. 
#                Default is adapted as a function of effective sample size in each model.
#                See 'Phillips, R. V., van der Laan, M. J., Lee, H., & Gruber, S. (2023).
#                Practical considerations for specifying a super learner.
#                International Journal of Epidemiology, 52(4), 1276-1285.' for details
#
# MSM.form:      The right hand side of a formula for the MSM relating the hazards 
#                to the exposure and time (optional).
#                The formula can only involve terms for the exposure (with the same name
#                as in dat) and "time", which is coded as time = 1, ..., K. 
#                If the name of Avar in dat is "trt" some examples are ~ trt + time,
#                ~ trt + as.factor(time) or trt + time + trt*time.
#
# Print:         TRUE or FALSE, whether the function should print the main results (default = TRUE)
#
# gA:            Optional argument for supplying external treatment probabilities. 
#                If A is binary gA must be a vector of P(A = 1|L0, Lvar[[1]]).
#                If A is categorical gA must be a matrix with column j being P(A = j|L0, Lvar[[1]]).
#                Default is NULL, meaning that gA is estimated within the function


#### Value
# The function returns a list with the following objects
# St:            A matrix of the estimated survival probabilities, their standard error, and 95%
#                confidence intervals at the different time points according to possible exposure levels
# MSM:           A matrix of the estimated parameters of the working marginal structural model, their
#                standard error, and 95% confidence intervals. Only returned if MSM.form was supplied.
# vcov:          The estimated variance covariance matrix of the working marginal structural model.
#                Only returned if MSM.form was supplied.
# ATE:		       Estimated average treatment effects at the different time points.
# SE.ATE:	       Standard error of the average treatment effects.
# LL.ATE:        Lower limit of the 95% confidence intervals for the average treatment effects.
# UL.ATE:        Upper limit of the 95% confidence intervals for the average treatment effects.
# CV.details:    Number of cross-validation folds that were effectively used if SuperLearner was used.
#                VA for modeling the treatment probability, VC for modeling the censoring probabilities
#                VY for modeling the outcome probabilities.
# gA:            The predicted treatment probabilities. If A is binary gA is a vector of
#                P(A = 1|L0, Lvar[[1]]). If A is categorical gA is a matrix with column j being
#                P(A = j|L0, Lvar[[1]]);
# gC:            The predicted censoring probabilities. gC is a list of length K, where each element
#                is a vector of length = n with P(C_t = 0 | C_{t-1} = 0, A, L_t, Y_{t-1} = 0).
#                NAs are inserted where C_{t-1} = 1 or Y_{t-1} = 1.



#### Details
# * On the use of glmnet and screen.glmnet for modeling outcome probabilities:
# The TMLE algorithm coded here uses iterated conditional expectations and
# needs to model some outcome probabilities. Those outcome probabilities are
# continuous and bounded in [0,1]. The default version of glmnet and screen.glmnet
# do not allow this. Here, we use a modified version where the outcome probabilities
# are treated as continuous and predictions outside the bounds are brought 
# back at the bounds.


#### Function

surv.TMLE = function(dat, Yvar, Cvar, Avar, Lvar, L0var = NULL,
                     lookback = NULL, lookbackY = lookback, lookbackC = lookback,
                     Ymod = "parametric", Cmod = "parametric", Amod = "parametric",
                     SL.library = c("SL.glm", "SL.glm.interaction"),
                     Y.SL.library = SL.library, C.SL.library = SL.library,
                     A.SL.library = SL.library, gbound = 0.025, V = NULL,
                     MSM.form = NULL, Print = TRUE, gA = NULL){

  ### Sample size (n) and number of time points (K)
  n = nrow(dat);
  K = length(Yvar);
  
  
  #### Error checks
  ## Verifications for dat
  if(!is.data.frame(dat)) stop("dat must be a data frame");

  ## Verifications for Yvar
  if(!is.character(Yvar) & !is.numeric(Yvar)) stop("Yvar must either be a character vector or a numeric vector");
  if(is.numeric(Yvar)){
    if(min(Yvar) <= 0) stop("At least one of the indices of Yvar is <= 0...");
    if(max(Yvar%%1) > 0) stop("At least one of the indices of Yvar is not an integer.");
  }else{ # Yvar is a character
    if(min(Yvar %in% names(dat)) == 0) stop("At least one of the names supplied for Yvar is not a variable name in dat.");
  }
  if(max(dat[,Yvar], na.rm = TRUE) > 1 |
     min(dat[,Yvar], na.rm = TRUE) < 0 |
     max(dat[,Yvar] %% 1, na.rm = TRUE) > 0) stop("Yvar is not coded 0/1"); 

  ## Verifications for Cvar  
  if(!is.character(Cvar) & !is.numeric(Cvar)) stop("Cvar must either be a character vector or a numeric vector");
  if(is.numeric(Cvar)){
    if(min(Cvar) <= 0) stop("At least one of the indices of Cvar is <= 0...");
    if(max(Cvar%%1) > 0) stop("At least one of the indices of Cvar is not an integer.");
  }else{ # Cvar is a character
    if(min(Cvar %in% names(dat)) == 0) stop("At least one of the names supplied for Cvar is not a variable name in dat.");
  }
  if(max(dat[,Cvar], na.rm = TRUE) > 1 |
     min(dat[,Cvar], na.rm = TRUE) < 0 |
     max(dat[,Cvar] %% 1, na.rm = TRUE) > 0) stop("Cvar is not coded 0/1"); 
  if(length(Cvar) != length(Yvar)) stop("Yvar and Cvar should have the same length");


  ## Verifications for Avar
  if(!is.numeric(Avar) & !is.character(Avar)) stop("Avar must a numeric value or a character");
  if(length(Avar) > 1) stop("Avar must a numeric value or a character of length 1");
  if(is.numeric(Avar)){
    if(min(Avar) <= 0) stop("Avar is <= 0...");
    if(max(Avar%%1) > 0) stop("Avar is not an integer.");
  }else{ # Avar is a character
    if(!(Avar %in% names(dat))) stop("The name supplied for Avar is not a variable name in dat.");
  }
  if(min(table(dat[,Avar])) == 1){ 
    stop("One level of A has a cell count of 1. Either A is continuous (not supported) or data is too sparse...");
  }
  if(length(table(dat[,Avar])) > 2){
    if(!is.factor(dat[,Avar])) stop("When the exposure has more than 2 levels, it must be of type factor");
  }else{ # Recode Avar as 0/1
    if(max(dat[,Avar]) != 1 | min(dat[,Avar]) != 0){
      cat(paste0("Avar has been recoded such that Avar = ", max(dat[,Avar]),
                 " => Avar = 1, and Avar = ", min(dat[,Avar]), " => Avar = 0.\n")); 
    }
    dat[,Avar] = 1*(dat[,Avar] == max(dat[,Avar]));
  }

  ## Verifications for Lvar  
  if(!is.list(Lvar)) stop("Lvar must be a list");
  if(length(Lvar) != length(Yvar)) stop("Lvar and Yvar must have the same length");
  for(i in 1:length(Lvar)){
    if(!is.character(Lvar[[i]]) & !is.numeric(Lvar[[i]])) stop("Lvar must either be a list of character vectors or a list of numeric vectors");
    if(is.numeric(Lvar[[i]])){
      if(min(Lvar[[i]]) <= 0) stop("At least one of the indices of Lvar is <= 0...");
      if(max(Lvar[[i]]%%1) > 0) stop("At least one of the indices of Lvar is not an integer.");
    }else{ # Lvar is a character
      if(min(Lvar[[i]] %in% names(dat)) == 0) stop("At least one of the names supplied for Lvar is not a variable name in dat.");
    }
  }

  ## Verifications for L0var
  if(!is.null(L0var)){
    if(!is.character(L0var) & !is.numeric(L0var)) stop("When supplied, L0var must either be a character vector or a numeric vector");
    if(is.numeric(L0var)){
      if(min(L0var) <= 0) stop("At least one of the indices of L0var is <= 0...");
      if(max(L0var%%1) > 0) stop("At least one of the indices of L0var is not an integer.");
    }else{ # L0var is a character
      if(min(L0var %in% names(dat)) == 0) stop("At least one of the names supplied for L0var is not a variable name in dat.");
    }
  }

  ## Verifications for lookback
  if(!is.null(lookback)){
    if(!is.numeric(lookback)) stop("When supplied, lookback must be a numeric value");
    if(length(lookback) > 1) stop("lookback should be of length 1");
    if(lookback < 1) stop("lookback should be >= 1");
    if(lookback %% 1 != 0){
      warning("lookback was not an integer and has been rounded up");
      lookback = ceiling(lookback);
    }
  }

  ## Verifications for lookbackY
  if(!is.null(lookbackY)){
    if(!is.numeric(lookbackY)) stop("When supplied, lookbackY must be a numeric value");
    if(length(lookbackY) > 1) stop("lookbackY should be of length 1");
    if(lookbackY < 1) stop("lookbackY should be >= 1");
    if(lookbackY %% 1 != 0){
      warning("lookbackY was not an integer and has been rounded up");
      lookbackY = ceiling(lookbackY);
    }
  }

  ## Verifications for lookbackC
  if(!is.null(lookbackC)){
    if(!is.numeric(lookbackC)) stop("When supplied, lookbackC must be a numeric value");
    if(length(lookbackC) > 1) stop("lookbackC should be of length 1");
    if(lookbackC < 1) stop("lookbackC should be >= 1");
    if(lookbackC %% 1 != 0){
      warning("lookbackC was not an integer and has been rounded up");
      lookbackC = ceiling(lookbackC);
    }
  }

  ## Verifications for Ymod
  if(!is.character(Ymod)) stop("Ymod must either be 'parametric' or 'SL'");
  if(length(Ymod) > 1) stop("Ymod must be of length 1");
  if(Ymod != "parametric" & Ymod != "SL") stop("Ymod must either be 'parametric' or 'SL'");

  ## Verifications for Cmod
  if(!is.character(Cmod)) stop("Cmod must either be 'parametric' or 'SL'");
  if(length(Cmod) > 1) stop("Cmod must be of length 1");
  if(Cmod != "parametric" & Cmod != "SL") stop("Cmod must either be 'parametric' or 'SL'");

  ## Verifications for Amod
  if(!is.character(Amod)) stop("Amod must either be 'parametric' or 'SL'");
  if(length(Amod) > 1) stop("Amod must be of length 1");
  if(Amod != "parametric" & Amod != "SL") stop("Amod must either be 'parametric' or 'SL'");

  ## Verifications for SL.library
  if(Ymod == "SL" | Cmod == "SL" | Amod == "SL"){
    if(is.list(SL.library)){
      learners = unlist(SL.library);
    } else {learners = SL.library};
    if(!is.character(learners)) stop("SL.library must be a character vector or a list of character vectors");

    # The following code was extracted/adapted from listWrappers
    everything = sort(getNamespaceExports("SuperLearner"));
    SL.algos = c(everything[grepl(pattern = "^[S]L", everything)], everything[grepl(pattern = "screen", everything)]);

    if(min(learners %in% SL.algos) == 0) warning("One of the learner supplied in SL.library is not in listWrappers(). \nThis may cause an error if the learner is not appropriate for SuperLearner.");
  }

  ## Verifications for Y.SL.library
  if(Ymod == "SL"){
    if(is.list(Y.SL.library)){
      ylearners = unlist(Y.SL.library);
    } else {ylearners = Y.SL.library};
    if(!is.character(ylearners)) stop("Y.SL.library must be a character vector or a list of character vectors");

    # The following code was extracted/adapted from listWrapper
    everything = sort(getNamespaceExports("SuperLearner"));
    SL.algos = c(everything[grepl(pattern = "^[S]L", everything)], everything[grepl(pattern = "screen", everything)]);

    if(min(ylearners %in% SL.algos) == 0) warning("One of the learner supplied in Y.SL.library is not in listWrappers(). \nThis may cause an error if the learner is not appropriate for SuperLearner.");
  }

  ## Verifications for C.SL.library
  if(Cmod == "SL"){
    if(is.list(C.SL.library)){
      clearners = unlist(C.SL.library);
    } else {clearners = C.SL.library};
    if(!is.character(clearners)) stop("C.SL.library must be a character vector or a list of character vectors");

    # The following code was extracted/adapted from listWrapper
    everything = sort(getNamespaceExports("SuperLearner"));
    SL.algos = c(everything[grepl(pattern = "^[S]L", everything)], everything[grepl(pattern = "screen", everything)]);

    if(min(clearners %in% SL.algos) == 0) warning("One of the learner supplied in C.SL.library is not in listWrappers(). \nThis may cause an error if the learner is not appropriate for SuperLearner.");
  }

  ## Verifications for A.SL.library
  if(Amod == "SL" & nlevels(as.factor(dat[,Avar])) == 2){
    if(is.list(A.SL.library)){
      alearners = unlist(A.SL.library);
    } else {alearners = A.SL.library};
    if(!is.character(alearners)) stop("A.SL.library must be a character vector or a list of character vectors");

    # The following code was extracted/adapted from listWrapper
    everything = sort(getNamespaceExports("SuperLearner"));
    SL.algos = c(everything[grepl(pattern = "^[S]L", everything)], everything[grepl(pattern = "screen", everything)]);

    if(min(alearners %in% SL.algos) == 0) warning("One of the learner supplied in A.SL.library is not in listWrappers(). \nThis may cause an error if the learner is not appropriate for SuperLearner.");
  }

  ## Verifications for gbound
  if(!is.numeric(gbound)) stop("gbound must be a numeric value between 0 and 1");
  if(length(gbound) > 1) stop("gbound must be a numeric value between 0 and 1");
  if(gbound < 0 | gbound > 1) stop("gbound must be a numeric value between 0 and 1");


  ## Verifications for V
  if(is.null(V)){ # Default data-adaptive behavior
    VA = NULL; VY = NULL; VC = NULL;
  } else{
    if(!is.numeric(V)) stop("V must be a numeric value >= 2");
    if(length(V) > 1) stop("V must be a numeric value >= 2");
    if(V < 2) stop("V must be a numeric value >= 2");    
    VA = V; VY = V; VC = V;
  }


  ## Verification for MSM.form
  if(!is.null(MSM.form)){
    if(!inherits(MSM.form, "formula")) stop("When supplied, MSM.form must be a formula");
    if(any(!(all.vars(MSM.form) %in% c("time", Avar)))){
      stop("The formula for MSM.form includes variable names other than the name of the exposure
            and time");
    }
  }

  ## Verification of Print
  if(!(Print == TRUE | Print == FALSE)) stop("Print must either be TRUE or FALSE");


  ### Compute lower and upper bounds for g estimates
  gbounds = c(min(gbound, 1 - gbound), max(gbound, 1 - gbound)); # Bounds for g estimate

  
  ### Initialize an object to output number of cross-validation folds
  CV.details = list();

  
  ### If Y_{t-1} = 1, replace Y_t by 1 and C_t by 0
  ### If C_{t-1} = 1, replace Y_t by 0 and C_t by 1
  ### If C_t = 1, replace Y_t by 0
  for(j in 2:K){
    dat[,Yvar[j]][dat[,Yvar[j-1]] == 1] = 1;
    dat[,Cvar[j]][dat[,Yvar[j-1]] == 1] = 0;
    dat[,Yvar[j]][dat[,Cvar[j-1]] == 1] = 0;
    dat[,Cvar[j]][dat[,Cvar[j-1]] == 1] = 1;
  }
  for(j in 1:K){
    dat[,Yvar[j]][dat[,Cvar[j]] == 1] = 0;
  }
  

  ### Verify if there is any missing data
  if(anyNA(dat[, Avar]) | anyNA(dat[, Yvar]) |
     anyNA(dat[, Cvar]) | anyNA(dat[, L0var]) |
     anyNA(dat[, unlist(Lvar[[1]])])) stop("dat contains missing data for at least one of the variables to be used.\nMissing data are not currently allowed. You may consider performing multiple imputations.");
  for(j in 2:K){
    if(any(dat[,Yvar[j-1]] == 0 & dat[,Cvar[j-1]] == 0 &
       is.na(dat[, unlist(Lvar[[j]])]))) stop("dat contains missing data for at least one of the variables to be used.\nMissing data are not currently allowed. You may consider performing multiple imputations.");
  }
  
  ### Verification for gA
  if(!is.null(gA)){
    if(!is.numeric(gA)) stop("gA must contain numeric values only");
    if(anyNA(gA)) stop("gA must not contain missing values");
    if(length(table(dat[,Avar])) == 2){ # A is binary
      if(length(gA) != n) stop("When Avar is binary, gA must be of length n");
    } else{ # A is categorical
      if(!is.matrix(gA)) stop("When Avar is categorical, gA must be a matrix");
      if(ncol(gA) != length(table(dat[,Avar]))) stop("Then number of columns of gA does not 
                                                     correspond with the number of levels of Avar");
      if(nrow(gA) != n) stop("The number of rows of gA is not equal to the sample size n");
    }
    if(max(gA) >= 1 | min(gA) <= 0) stop("gA values must be greater than 0 and lower than 1");
  }
  
  
  #### Modeling the exposure
  if(is.null(gA)){
    if(nlevels(as.factor(dat[,Avar])) == 2){ # A is binary
      if(is.null(L0var)){ # No time-fixed covariates
        Aform = paste(names(dat[,Avar, drop = FALSE]), " ~ ",
                      paste(names(dat[, Lvar[[1]], drop = FALSE]), collapse = " + "), sep = "");
      }else{ # Some time-fixed covariates
        Aform = paste(names(dat[,Avar, drop = FALSE]), " ~ ",
                      paste(names(dat[, L0var, drop = FALSE]),
                            names(dat[, Lvar[[1]], drop = FALSE]), collapse = " + ", sep = " + "), sep = "");
      }
      if(Amod == "parametric"){
        gA = glm(Aform, data = dat, family = "binomial", maxit = 500)$fitted;
      }else{ #Amod == "SL"
        if(is.null(V)){ # Data-adaptive choice of VA
          neff = min(5*sum(dat[, Avar]),5*sum(1 - dat[, Avar]), n);
          VA = neff;
          if(neff >= 30) VA = 20;
          if(neff >= 500) VA = 10;
          if(neff >= 5000) VA = 5;
          if(neff >= 10000) VA = 2;
          CV.details$VA = VA;
        }
        X = as.data.frame(model.matrix(as.formula(Aform), data = dat)[,-1]); # Obtain a matrix to accommodate factors
        names(X) = paste0("X", 1:ncol(X));
        mod.A = SuperLearner(Y = dat[, Avar], X = X, family = "binomial",
                             SL.library = A.SL.library, cvControl = list(V = VA));
        gA = predict(mod.A, OnlySL = TRUE)$pred; 
      } 
      gA = pmax(pmin(gA, gbounds[2]), gbounds[1]); # Bounding gA
    }else{ # A is categorical
      nlevel = nlevels(dat[, Avar]); # Number of levels of A
      levelA = levels(dat[, Avar]); # Levels of A
      if(Amod == "parametric"){
        dat$id = 1:n;
        Avarname = names(dat[,Avar, drop = FALSE]);
        ds = mlogit.data(data = dat, shape = "wide", choice = Avarname, varying = NULL, idvar = id);
        if(is.null(L0var)){ # No time-fixed covariates
          Aform = paste(names(dat[,Avar, drop = FALSE]), " ~ 1 | ",
                        paste(names(dat[, Lvar[[1]], drop = FALSE]), collapse = " + "), sep = "");
        }else{ # Some time-fixed covariates
          Aform = paste(names(dat[,Avar, drop = FALSE]), " ~ 1 | ",
                        paste(names(dat[, L0var, drop = FALSE]),
                              names(dat[, Lvar[[1]], drop = FALSE]), sep = " + ", collapse = " + "), sep = "");
        }
        Aform = as.formula(Aform);      
        mod.A = mlogit(Aform, data = ds);
        gA = predict(mod.A, type = "probs", newdata = ds[order(ds$id),]);
        Indicator = matrix(NA, nrow = n, ncol = nlevel); # Indicator matrix I(A = a)
        for(i in 1:nlevel){
          Indicator[,i] = (dat[, Avar] == levelA[i]);
        }
        gA = c(as.matrix(gA)); # Vector form of gA;
        gA = pmin(gA, gbounds[2]); gA = pmax(gA, gbounds[1]); # Bounding gA
        gA = matrix(gA, ncol = nlevel, nrow = n);
      }else{ #Amod == "SL"
        if(is.null(L0var)){ # No time-fixed covariates
          Aform = paste(names(dat[,Avar, drop = FALSE]), " ~ ",
                        paste(names(dat[, Lvar[[1]], drop = FALSE]), collapse = " + "), sep = "");
        }else{ # Some time-fixed covariates
          Aform = paste(names(dat[,Avar, drop = FALSE]), " ~ ",
                        paste(names(dat[, L0var, drop = FALSE]),
                              names(dat[, Lvar[[1]], drop = FALSE]), sep = " + ", collapse = " + "), sep = "");
        }
        X = model.matrix(as.formula(Aform), data = dat)[,-1]; # Obtain a matrix to accomodate factors
        names(X) = paste0("X", 1:ncol(X));
        mod.A = polyclass(dat[,Avar], X);
        gA = matrix(NA, nrow = n, ncol = nlevel); 
        for(i in 1:nlevel){ # Classes are ordered the same as the levels of A
          gA[,i] = ppolyclass(i, X, mod.A);
        } 
        gA = c(as.matrix(gA)); # Vector form of gA;
        gA = pmin(gA, gbounds[2]); gA = pmax(gA, gbounds[1]); # Bounding gA
        gA = matrix(gA, ncol = nlevel, nrow = n);
      }
    }
  }
  ## Note :
  # If A is binary gA is a vector of P(A = 1|L0, Lvar[[1]]).
  # If A is categorical gA is a matrix with column j being P(A = j|L0, Lvar[[1]]);


  #### Modeling the censoring
  gC = list();
  if(Cmod == "SL" & is.null(V)) CV.details$VC = rep(NA, K);
  if(is.null(lookbackC)) lookbackC == Inf;
  if(is.null(L0var)){ # No time-fixed covariates
    if(Cmod == "parametric"){
      for(j in 1:K){
        j_lookback = max(1, j - lookbackC + 1); # First index of the Lvar to use
        Lvarj = unlist(Lvar[j_lookback:j]); # The variables to be used

        # The formula for Cj:
        Cform = paste(names(dat[,Cvar[j], drop = FALSE]), " ~ ",
                      names(dat[,Avar, drop = FALSE]), " + ",
                      paste(names(dat[, Lvarj, drop = FALSE]), collapse = " + "), sep = "");

        # P(Cj = 0|A,L, C_j-1 = 0, Y_j-1 = 0):      
        if(j == 1){ # First time point
          if(max(dat[,Cvar[j]]) == 0){ # No censoring
            gC[[j]] = rep(1, n);
          }else{ # At least some censoring
            gC[[j]] = pmax(1 - glm(Cform, data = dat, family = "binomial", maxit = 500)$fitted, gbounds[1]); 
          }
        }else{ # Not the first time point
          Cvarj_1 = names(dat[,Cvar[j-1], drop = FALSE]); # previous Cvar
          Yvarj_1 = names(dat[,Yvar[j-1], drop = FALSE]); # previous Yvar
          if(max(dat[,Cvar[j]][dat[,Cvar[j-1]] == 0]) == 0){ # No censoring
            gC[[j]] = rep(NA, n);
            gC[[j]][dat[, Cvarj_1] == 0 & dat[,Yvarj_1] == 0] = 1;
          }else{  # At least some censoring
            gC[[j]] = rep(NA, n);
            used = which(dat[,Cvarj_1] == 0 & dat[,Yvarj_1] == 0);
            gC[[j]][used] = pmax(1 - glm(Cform, data = dat[used,], family = "binomial", maxit = 500)$fitted, gbounds[1]);
          }
        }
      } # End of loop on time points 
    }else{ # Cmod == "SL"
      for(j in 1:K){
        j_lookback = max(1, j - lookbackC + 1); # First index of the Lvar to use
        Lvarj = unlist(Lvar[j_lookback:j]); # The variables to be used

        # The formula for Cj:
        Cform = paste(names(dat[,Cvar[j], drop = FALSE]), " ~ ",
                      names(dat[,Avar, drop = FALSE]), " + ",
                      paste(names(dat[, Lvarj, drop = FALSE]), collapse = " + "), sep = "");

        # P(Cj = 0|A,L, C_j-1 = 0, Y_j-1 = 0):      
        if(j == 1){ # First time point
          if(max(dat[,Cvar[j]]) == 0){ # No censoring
            gC[[j]] = rep(1, n);
          }else{ # At least some censoring
            if(is.null(V)){ # Data-adaptive choice of VC
              neff = min(5*sum(dat[, Cvar[[j]]], na.rm = TRUE),
                         5*sum(1 - dat[, Cvar[[j]]], na.rm = TRUE), n);
              VC = neff;
              if(neff >= 30) VC = 20;
              if(neff >= 500) VC = 10;
              if(neff >= 5000) VC = 5;
              if(neff >= 10000) VC = 2;
              CV.details$VC[j] = VC;
            }
            # The design matrix
            X = as.data.frame(model.matrix(as.formula(Cform), data = dat)[,-1]);
            names(X) = paste0("X", 1:ncol(X));
            suppressWarnings({
            mod.C = SuperLearner(Y = dat[, Cvar[[j]]], X = X, family = "binomial",
                                 SL.library = C.SL.library, cvControl = list(V = VC));});
            gC[[j]] = pmax(1 - predict(mod.C, OnlySL = TRUE)$pred, gbounds[1]); 
          }
        }else{ # not the first time point
          Cvarj_1 = names(dat[,Cvar[j-1], drop = FALSE]); # previous Cvar
          Yvarj_1 = names(dat[,Yvar[j-1], drop = FALSE]); # previous Yvar
          if(max(dat[,Cvar[j]][dat[,Cvar[j-1]] == 0]) == 0){ # No censoring
            gC[[j]] = rep(NA, n);
            gC[[j]][dat[, Cvarj_1] == 0 & dat[,Yvarj_1] == 0] = 1;
          }else{  # At least some censoring
            gC[[j]] = rep(NA, n);
            used = which(dat[, Cvarj_1] == 0 & dat[, Yvarj_1] == 0);
            X = as.data.frame(model.matrix(as.formula(Cform), data = dat[used,])[,-1]);
            names(X) = paste0("X", 1:ncol(X));
            if(is.null(V)){ # Data-adaptive choice of VC
              neff = min(5*sum(dat[, Cvar[[j]]], na.rm = TRUE),
                         5*sum(1 - dat[, Cvar[[j]]], na.rm = TRUE), sum(used));
              VC = neff;
              if(neff >= 30) VC = 20;
              if(neff >= 500) VC = 10;
              if(neff >= 5000) VC = 5;
              if(neff >= 10000) VC = 2;
              CV.details$VC[j] = VC;
            }
            suppressWarnings({
            mod.C = SuperLearner(Y = dat[used, Cvar[[j]]], X = X, family = "binomial",
                                 SL.library = C.SL.library, cvControl = list(V = VC));});
            gC[[j]][used] = pmax(1 - predict(mod.C, OnlySL = TRUE)$pred, gbounds[1]); 
          }
        }
      } # End of loop on time points 
    } # End of Cmod == "SL"
  }else{ # There are time-fixed covariates
    if(Cmod == "parametric"){
      for(j in 1:K){
        j_lookback = max(1, j - lookbackC + 1); # First index of the Lvar to use
        Lvarj = unlist(Lvar[j_lookback:j]); # The variables to be used

        # The formula for Cj:
        Cform = paste(names(dat[,Cvar[j], drop = FALSE]), " ~ ",
                      names(dat[,Avar, drop = FALSE]), " + ",
                      paste(names(dat[, L0var, drop = FALSE]), collapse = " + "), " + ",
                      paste(names(dat[, Lvarj, drop = FALSE]), collapse = " + "), sep = "");

        # P(Cj = 0|A,L, C_j-1 = 0, Y_j-1 = 0):      
        if(j == 1){ # First time point
          if(max(dat[,Cvar[j]]) == 0){ # No censoring
            gC[[j]] = rep(1, n);
          }else{ # At least some censoring
            gC[[j]] = pmax(1 - glm(Cform, data = dat, family = "binomial", maxit = 500)$fitted, gbounds[1]); 
          }
        }else{ # Not the first time point
          Cvarj_1 = names(dat[,Cvar[j-1], drop = FALSE]); # previous Cvar
          Yvarj_1 = names(dat[,Yvar[j-1], drop = FALSE]); # previous Yvar
          if(max(dat[,Cvar[j]][dat[,Cvar[j-1]] == 0]) == 0){ # No censoring
            gC[[j]] = rep(NA, n);
            gC[[j]][dat[, Cvarj_1] == 0 & dat[,Yvarj_1] == 0] = 1;
          }else{  # At least some censoring
            gC[[j]] = rep(NA, n);
            used = which(dat[, Cvarj_1] == 0 & dat[, Yvarj_1] == 0);
            gC[[j]][used] = pmax(1 - glm(Cform, data = dat[used,], family = "binomial", maxit = 500)$fitted, gbounds[1]);
          }
        }
      } # End of loop on time points 
    }else{ # Cmod == "SL"
      for(j in 1:K){
        j_lookback = max(1, j - lookbackC + 1); # First index of the Lvar to use
        Lvarj = unlist(Lvar[j_lookback:j]); # The variables to be used

        # The formula for Cj:
        Cform = paste(names(dat[,Cvar[j], drop = FALSE]), " ~ ",
                      names(dat[,Avar, drop = FALSE]), " + ",
                      paste(names(dat[, L0var, drop = FALSE]), collapse = " + "), " + ",
                      paste(names(dat[, Lvarj, drop = FALSE]), collapse = " + "), sep = "");

        # P(Cj = 0|A,L, C_j-1 = 0, Y_j-1 = 0):      
        if(j == 1){ # First time point
          if(max(dat[,Cvar[j]]) == 0){ # No censoring
            gC[[j]] = rep(1, n);
          } else{ # At least some censoring
            if(is.null(V)){ # Data-adaptive choice of VC
              neff = min(5*sum(dat[, Cvar[[j]]], na.rm = TRUE),
                         5*sum(1 - dat[, Cvar[[j]]], na.rm = TRUE), n);
              VC = neff;
              if(neff >= 30) VC = 20;
              if(neff >= 500) VC = 10;
              if(neff >= 5000) VC = 5;
              if(neff >= 10000) VC = 2;
              CV.details$VC[j] = VC;
            }
            # The design matrix
            X = as.data.frame(model.matrix(as.formula(Cform), data = dat)[,-1]);
            names(X) = paste0("X", 1:ncol(X));
            suppressWarnings({
            mod.C = SuperLearner(Y = dat[, Cvar[[j]]], X = X, family = "binomial",
                                 SL.library = C.SL.library, cvControl = list(V = VC));});
            gC[[j]] = pmax(1 - predict(mod.C, OnlySL = TRUE)$pred, gbounds[1]); 
          }
        }else{ # Not the first time point
          Cvarj_1 = names(dat[,Cvar[j-1], drop = FALSE]); # previous Cvar
          Yvarj_1 = names(dat[,Yvar[j-1], drop = FALSE]); # previous Yvar
          if(max(dat[,Cvar[j]][dat[,Cvar[j-1]] == 0]) == 0){ # No censoring
            gC[[j]] = rep(NA, n);
            gC[[j]][dat[, Cvarj_1] == 0 & dat[,Yvarj_1] == 0] = 1;
          }else{  # At least some censoring
            gC[[j]] = rep(NA, n);
            used = which(dat[, Cvarj_1] == 0 & dat[, Yvarj_1] == 0);
            X = as.data.frame(model.matrix(as.formula(Cform), data = dat[used,])[,-1]);
            names(X) = paste0("X", 1:ncol(X));
            if(is.null(V)){ # Data-adaptive choice of VC
              neff = min(5*sum(dat[, Cvar[[j]]], na.rm = TRUE),
                         5*sum(1 - dat[, Cvar[[j]]], na.rm = TRUE), sum(used));
              VC = neff;
              if(neff >= 30) VC = 20;
              if(neff >= 500) VC = 10;
              if(neff >= 5000) VC = 5;
              if(neff >= 10000) VC = 2;
              CV.details$VC[j] = VC;
            }
            suppressWarnings({
            mod.C = SuperLearner(Y = dat[used, Cvar[[j]]], X = X, family = "binomial",
                                 SL.library = C.SL.library, cvControl = list(V = VC));});
            gC[[j]][used] = pmax(1 - predict(mod.C, OnlySL = TRUE)$pred, gbounds[1]); 
          }
        }
      } # End of loop on time points 
    } # End of Cmod == "SL"
  } 
  ## Note:
  # gC is a list of length K, 
  # each element is a vector of length = n
  # P(C_t = 0 | C_{t-1} = 0, A, L_t, Y_{t-1} = 0)
  # NAs are inserted where C_{t-1} = 1 or Y_{t-1} = 1


  #### Compute the cumulative product of gC
  gC.cumul = gC;
  for(j in 2:K){
    gC.cumul[[j]] = gC[[j]]*gC.cumul[[j-1]];
  }


  #### Replace glmnet by glmnet by SL.glmnet_1_100_10 to accomodate continuous [0,1] data in ICE if needed

  if("SL.glmnet" %in% unlist(Y.SL.library) | "screen.glmnet" %in% unlist(Y.SL.library)){
    ## Replace "SL.glmnet" with "SL.glmnet_1_100_10"  
    replace_string <- function(input_vector) {
      gsub("SL.glmnet", "SL.glmnet_1_100_10", input_vector)
    }
    
    if(is.character(Y.SL.library)){
      Y.SL.library = replace_string(Y.SL.library); 
    } else{ # Y.SL.library is a list of character vectors
      Y.SL.library = lapply(Y.SL.library, replace_string);
    }
    
    ## Replace "screen.glmnet" with "screen.glmnet_1_10_200"
    replace_string <- function(input_vector) {
      gsub("screen.glmnet", "screen.glmnet_1_10_200", input_vector)
    }
    
    if(is.character(Y.SL.library)){
      Y.SL.library = replace_string(Y.SL.library); 
    } else{ # Y.SL.library is a list of character vectors
      Y.SL.library = lapply(Y.SL.library, replace_string);
    }
    
    ## The code to create SL.glmnet_1_100_10, supplied by Michael Shomaker
    SL.glmnet_base <- function (Y, X, newX, family, obsWeights, id, alpha = 1, nfolds = 10, 
                                nlambda = 100, useMin = TRUE, loss = "deviance", verbose=T, ...) 
    {
      #
      if(verbose==T){cat("SL.glmnet started with alpha=", alpha, ", ", nfolds, "-fold CV, and ", nlambda, " candidate lambdas. ", sep="")}
      start_time <- Sys.time()
      SuperLearner:::.SL.require("glmnet")
      # for ltmle 
      fam.init <- family$family
      Y <- as.vector(as.matrix(Y))
      if (all(Y == 0 | Y == 1)) {
        family$family <- "binomial"
      } else {
        family$family <- "gaussian"
      }
      fam.end <- family$family
      #
      if (!is.matrix(X)) {
        X <- model.matrix(~-1 + ., X)
        newX <- model.matrix(~-1 + ., newX)
      }
      fitCV <- glmnet::cv.glmnet(x = X, y = Y, weights = obsWeights, 
                                 lambda = NULL, type.measure = loss, nfolds = nfolds, 
                                 family = family$family, alpha = alpha, nlambda = nlambda, 
                                 ...)
      pred <- predict(fitCV, newx = newX, type = "response", 
                      s = ifelse(useMin, "lambda.min", "lambda.1se"))
      #
      if(fam.init=="binomial" & fam.end=="gaussian"){if(any(pred<0)){pred[pred<0]<-0};if(any(pred>1)){pred[pred>1]<-1};if(verbose==T){cat("Note: predictions falling outside [0,1] have been set as 0/1")}}
      #
      fit <- list(object = fitCV, useMin = useMin)
      class(fit) <- "SL.glmnet"
      out <- list(pred = pred, fit = fit)
      #
      end_time <- Sys.time()
      if(verbose==T){cat("SL.glmnet finished. Time:", round(difftime(end_time, start_time, units="mins"), digits=4), "mins \n\n")}
      #
      return(out)
    }
    assign("SL.glmnet_base", SL.glmnet_base, envir = .GlobalEnv)
    
    #
    make.SL.glmnet <- function(alpha=1, nlambda=100, nfolds=10, verbose=T){
      
      tuneGrid <- expand.grid(alpha=alpha, nlambda=nlambda, nfolds=nfolds)
      
      for (q in seq(nrow(tuneGrid))) {
        eval(parse(text = paste0("SL.glmnet_", paste(tuneGrid[q, ],collapse="_"),
                                 "<- function(..., alpha=", tuneGrid[q, 1], ", nlambda=",tuneGrid[q, 2], ", nfolds=", tuneGrid[q, 3],
                                 " , verbose = ", verbose,  ")
                         {SL.glmnet_base(..., alpha=alpha, nlambda=nlambda, nfolds=nfolds, verbose=verbose)}"
        )), envir = .GlobalEnv)
        
      }  
    }
    make.SL.glmnet(verbose = FALSE);
    
    ## The code to create screen.glmnet, also supplied by Michael Shomaker
    screen.cramersv_base <- function(Y, X, nscreen = 4, num_cat = 10, verbose = T, ...) {
      if(verbose==T){cat("screen.cramersv for", nscreen, "variables; ")}
      start_time <- Sys.time()
      SuperLearner:::.SL.require("vcd")
      #
      if (ncol(X) > nscreen) {
        dat <- cbind(Y, X)
        contin_var <- apply(dat, 2, function(var) length(unique(var)) > num_cat)
        make_categ <- function(var) {
          num_qu <- length(unique(quantile(var, prob = seq(0, 1, 0.2))))
          if(num_qu >2){ret<-cut(var, unique(quantile(var, prob = seq(0, 1, 0.2))), include.lowest = T)}
          if(num_qu<=2){ret<-cut(var, c(-Inf,unique(quantile(var, prob = seq(0, 1, 0.2))),Inf), include.lowest = T)}
          ret
        }
        if (any(contin_var)) {
          dat[, contin_var] <- apply(dat[, contin_var, drop = FALSE], 2, make_categ)
        }
        calc_cram_v <- function(x_var, y_var) vcd::assocstats(table(y_var, x_var))$cramer
        cramers_v <- apply(dat[, !colnames(dat) %in% "Y"], 2, calc_cram_v, y_var = dat[, "Y"])
        if(verbose==T){cat("screened:",colnames(X)[unname(rank(-cramers_v) <= nscreen)],"\n",sep=" ");
          cat("sample size:",dim(X)[1],"; "); 
        }
        whichVariable <- unname(rank(-cramers_v) <= nscreen)
      }else{
        if(verbose==T){cat("screened all", ncol(X), "variables \n")}
        whichVariable <- rep(TRUE, ncol(X))}
      end_time <- Sys.time()
      if(verbose==T){cat("time:", round(difftime(end_time, start_time, units="mins"), digits=4), "mins \n")}
      return(whichVariable)
    }
    assign("screen.cramersv_base", screen.cramersv_base, envir = .GlobalEnv)
    
    screen.glmnet_base <- function(Y, X, family, alpha = 1, verbose=T,
                                   nfolds = 10, nlambda = 200, nscreen=2, ...){
      if(verbose==T){cat("screen.glmnet with alpha=", alpha, " and ", nfolds, " fold CV \n", sep="")}
      start_time <- Sys.time()
      SuperLearner:::.SL.require("glmnet")
      # relevant for column names but shouldn't be a matrix anyways
      X <- as.data.frame(X)
      # chose family dependent upon response
      Y <- as.vector(as.matrix(Y))
      if (all(Y == 0 | Y == 1)) {
        family$family <- "binomial"
      } else {
        family$family <- "gaussian"
      }
      saveY<-Y;saveX<-X
      # needed for var names to select from levels of factors later on
      if (ncol(X) > 26 * 27) stop("Find further column names for X!")
      let <- c(letters, sort(do.call("paste0", expand.grid(letters, letters[1:26]))))
      names(X) <- let[1:ncol(X)]
      # factors are coded as dummies which are standardized in cv.glmnet()
      # intercept is not in model.matrix() because its already in cv.glmnet()
      is_fact_var <- sapply(X, is.factor)
      X <- try(model.matrix(~ -1 + ., data = X), silent = FALSE)
      successfulfit <- FALSE
      fitCV <- try(glmnet::cv.glmnet(
        x = X, y = Y, lambda = NULL, type.measure = "deviance",
        nfolds = nfolds, family = family$family, alpha = alpha,
        nlambda = nlambda, keep = T
      ), silent = TRUE)
      # if no variable was selected, penalization might have been too strong, try log(lambda)
      if (all(fitCV$nzero == 0) | all(is.na(fitCV$nzero))) {
        fitCV <- try(glmnet::cv.glmnet(
          x = X, y = Y, lambda = log(fitCV$glmnet.fit$lambda + 1), type.measure = "deviance",
          nfolds = nfolds, family = family$family, alpha = alpha, keep = T
        ), silent = TRUE)
      }
      if(class(fitCV)=="try-error"){successfulfit <- FALSE}else{successfulfit <- TRUE}
      whichVariable <- NULL
      if(successfulfit==TRUE){
        coefs <- coef(fitCV$glmnet.fit, s = fitCV$lambda.min)
        if(all(coefs[-1]==0)){whichVariable<-screen.cramersv_base(Y=saveY,X=saveX,nscreen=nscreen)
        if(verbose==T){cat("Lasso screened away all variables and screening was thus based on Cramer's V \n")}}else{  
          var_nms <- coefs@Dimnames[[1]]
          # Instead of Group Lasso:
          # If any level of a dummy coded factor is selected, the whole factor is selected
          if (any(is_fact_var)) {
            nms_fac <- names(which(is_fact_var))
            is_selected <- coefs[-1] != 0 # drop intercept
            # model.matrix adds numbers to dummy coded factors which we need to get rid of
            var_nms_sel <- gsub("[^::a-z::]", "", var_nms[-1][is_selected])
            sel_fac <- nms_fac[nms_fac %in% var_nms_sel]
            sel_numer <- var_nms_sel[!var_nms_sel %in% sel_fac]
            all_sel_vars <- c(sel_fac, sel_numer)
            whichVariable <- names(is_fact_var) %in% all_sel_vars
          } else {
            # metric variables only
            whichVariable <- coefs[-1] != 0
          }}
      }  
      if(is.null(whichVariable)){
        whichVariable<-screen.cramersv_base(Y,X)
        if(verbose==T){cat("Lasso failed and screening was based on Cramer's V\n")}}
      if(verbose==T){cat("screened ",sum(whichVariable)," variables: ",paste(colnames(saveX)[whichVariable],sep=" "),"\n");
        cat("sample size:",dim(X)[1],"; ")
      }
      end_time <- Sys.time()
      if(verbose==T){cat("time:", round(difftime(end_time, start_time, units="mins"), digits=4), "mins \n")}
      return(whichVariable)
    }
    assign("screen.glmnet_base", screen.glmnet_base, envir = .GlobalEnv)
    
    make.glmnet <- function(alpha=1, verbose=T, nfolds = 10, nlambda = 200){
      
      tuneGrid <- expand.grid(alpha=alpha, nfolds=nfolds, nlambda=nlambda)
      
      for (q in seq(nrow(tuneGrid))) {
        eval(parse(text = paste0("screen.glmnet_", paste(tuneGrid[q, ],collapse="_"),
                                 "<- function(..., alpha = ", tuneGrid[q, 1], ", verbose = ", verbose, ", nfolds =", tuneGrid[q, 2],
                                 ", nlambda=", tuneGrid[q, 3], ")
                         {screen.glmnet_base(..., alpha=alpha, verbose=verbose, nfolds=nfolds, nlambda=nlambda)}"
        )), envir = .GlobalEnv)
      }  
    }
    make.glmnet(verbose = FALSE);
    
    message("Note: SL.glmnet and screen.glmnet have been replaced by modified versions that accomodate non-binary data in the iterated conditional expectation step.\nSee the details section of the function's documentation for more on this.");
  }  


  #### Modeling the outcome
  if(Ymod == "SL" & is.null(V)) CV.details$VY = numeric(nlevels(as.factor(dat[,Avar]))*K*(K+1)/2);
  u = 1;
  if(is.null(lookbackY)) lookbackY == Inf;
  St = matrix(NA, nrow = K, ncol = nlevels(as.factor(dat[, Avar]))); # Object that will contain the estimated survival curves
  ICt = array(0, dim = c(n, ncol = nlevels(as.factor(dat[, Avar])), K)); # Object that will contain the empirical efficient influence curve
  for(k in K:1){
    Q = Qs = d = list(); # Initializing objects for the Qs, Q-stars and ds
    Qs[[k+1]] = matrix(dat[, Yvar[k]], 
                       nrow = n,
                       ncol = nlevels(as.factor(dat[,Avar])),
                       byrow = FALSE); 
    for(j in k:1){
      Q[[j]] = matrix(NA, nrow = n, ncol = nlevels(as.factor(dat[,Avar])));
      Qs[[j]] = matrix(NA, nrow = n, ncol = nlevels(as.factor(dat[,Avar])));
      d[[j]] = matrix(NA, nrow = n, ncol = nlevels(as.factor(dat[,Avar])));
      j_lookback = max(1, j - lookbackY + 1); # First index of the Lvar to use
      Lvarj = unlist(Lvar[j_lookback:j]); # The variables to be used
      ak = 1;
      for(a in sort(unique(dat[,Avar]))){
        if(Ymod == "parametric"){
          if(is.null(L0var)){ # No time-fixed covariates
            ## The formula for Qj:
            Qform = paste("YSL ~ ",
                           names(dat[,Avar, drop = FALSE]), " + ",
                           paste(names(dat[, Lvarj, drop = FALSE]), collapse = " + "), sep = "");
          }else{# There are time-fixed covariates
            Qform = paste("YSL ~ ",
                           names(dat[,Avar, drop = FALSE]), " + ",
                           paste(names(dat[, L0var, drop = FALSE]), collapse = " + "), " + ",
                           paste(names(dat[, Lvarj, drop = FALSE]), collapse = " + "), sep = "");
          }

          ## The Q-model
          if(j == 1){ # First time point
            dat2 = dat;
            dat2$YSL = Qs[[j+1]][,ak];
            dat2 = dat2[dat[, Cvar[[j]]] == 0,];
            modQ = suppressWarnings(glm(Qform, data = dat2, family = "binomial", maxit = 500));
          }else{ # Note the first time point
            dat2 = dat;
            dat2$YSL = Qs[[j+1]][,ak];
            dat2 = dat2[dat[, Cvar[[j]]] == 0 & dat[,Yvar[[j-1]]] == 0,];
            modQ = suppressWarnings(glm(Qform, data = dat2, family = "binomial", maxit = 500));
          }
          ## Computing Qj
          if(j == 1){ # First time point
            newdat = dat;
          }else{
            newdat = dat[dat[,Cvar[[j-1]]] == 0,];
          }
          newdat[, Avar] = a;
          if(j == 1){
            Q[[j]][,ak] = predict(modQ, newdata = newdat, type = "res");
          }else{
            Q[[j]][dat[,Cvar[[j-1]]] == 0,ak] = predict(modQ, newdata = newdat, type = "res");
          }
          if(j != 1) Q[[j]][dat[,Yvar[[j-1]]] == 1 & dat[,Cvar[[j-1]]] == 0, ak] = 1; # If Y_{j-1} = 1 then Q_j = 1
        }else{#SL
          if(j == 1){
            YSL = Qs[[j+1]][dat[Cvar[[j]]] == 0,ak];
          }else{
            YSL = Qs[[j+1]][dat[Cvar[[j]]] == 0 & dat[,Yvar[[j-1]]] == 0,ak];
          }
          if(is.null(L0var)){ # No time-fixed covariates
            ## The formula for Qj:
            Qform = paste("YSL ~ ",
                           names(dat[,Avar, drop = FALSE]), " + ",
                           paste(names(dat[, Lvarj, drop = FALSE]), collapse = " + "), sep = "");
          }else{# There are time-fixed covariates
            Qform = paste("YSL ~ ",
                           names(dat[,Avar, drop = FALSE]), " + ",
                           paste(names(dat[, L0var, drop = FALSE]), collapse = " + "), " + ",
                           paste(names(dat[, Lvarj, drop = FALSE]), collapse = " + "), sep = "");
          }

          ## The Q-model
          # The design matrix
          if(j == 1){
            X = as.data.frame(model.matrix(as.formula(Qform), data = dat[dat[Cvar[[j]]] == 0,])[,-1]);
          }else{
            X = as.data.frame(model.matrix(as.formula(Qform), data = dat[dat[Cvar[[j]]] == 0 & dat[,Yvar[[j-1]]] == 0,])[,-1]);
          }
          if(is.null(V)){ # Data-adaptive choice of VY
            neff = min(5*sum(YSL), 5*sum(1 - YSL), length(YSL));
            if(is.na(neff)) neff = 0;
            VY = neff;
            if(neff >= 30) VY = 20;
            if(neff >= 500) VY = 10;
            if(neff >= 5000) VY = 5;
            if(neff >= 10000) VY = 2;
            CV.details$VY[u] = VY;
            u = u + 1;
          }
          if(min(5*sum(YSL), 5*sum(1 - YSL), length(YSL)) > 0){ # There were both events and survivals
            modQ = suppressWarnings(SuperLearner(Y = YSL, X = X, family = "binomial",
                                                 SL.library = Y.SL.library, cvControl = list(V = VY)));
          } else{ # Either all events or all survival
            modQ = NULL;
          }

          ## Computing Qj
          if(j == 1){ # First time point
            newdat = dat;
          }else{
            newdat = dat[dat[,Cvar[[j-1]]] == 0,];
          }
          newdat[, Avar] = a;
          if(is.null(modQ)){ # either all events or all survival; prediction = only observed value
            if(j == 1){
              Q[[j]][,ak] = YSL[1];
            }else{
              Q[[j]][dat[,Cvar[[j-1]]] == 0,ak] = YSL[1];
            }
          } else{ # There were both events and survivals
            if(j == 1){
              suppressMessages({
              attach(newdat); # Attaching and detaching seem to prevents some errors that happen with screeners
              Q[[j]][,ak] = predict(modQ, OnlySL = TRUE, newdata = newdat)$pred;
              detach(newdat);});
            }else{
              suppressMessages({
              attach(newdat);
              Q[[j]][dat[,Cvar[[j-1]]] == 0,ak] = predict(modQ, OnlySL = TRUE, newdata = newdat)$pred;
              detach(newdat);});
            }
          }
        }#End else-SL

        ## Computing Hj
        if(nlevels(as.factor(dat[,Avar])) == 2){
          ga = (dat[,Avar] == 1)*gA + (dat[,Avar] == 0)*(1 - gA);
        }else{
          ga = gA[, ak];
        }
        Hj = (dat[, Avar] == a)*(dat[, Cvar[[j]]] == 0)/(ga*gC.cumul[[j]]);
        Hj[is.na(Hj)] = 0; # Hj = 0 if Cj = 0 or Y_{j-1} = 0

        ## Computing Qj-star
        if(!is.null(modQ)){
          if(j==1){
            epsilon = suppressWarnings(coef(glm(Qs[[j+1]][,ak] ~ 1 + offset(qlogis(Q[[j]][,ak])),
                                                weights = Hj,
                                                family = "binomial",
                                                subset = dat[,Cvar[[j]]] == 0)));
          }else{
            epsilon = suppressWarnings(coef(glm(Qs[[j+1]][,ak] ~ 1 + offset(qlogis(Q[[j]][,ak])),
                                                weights = Hj,
                                                family = "binomial",
                                                subset = dat[,Cvar[[j]]] == 0 & dat[,Yvar[[j-1]]] == 0)));
  
          }
        }else{
          epsilon = 0;
        }
        Qs[[j]][,ak] = plogis(qlogis(Q[[j]][,ak]) + epsilon);
        if(j != 1) Qs[[j]][dat[,Yvar[[j-1]]] == 1 & dat[,Cvar[[j-1]]] == 0, ak] = 1;
     
        ## Computing dj
        d[[j]][,ak] = Hj*(Qs[[j]][,ak] - Qs[[j+1]][,ak]);
        d[[j]][Hj == 0, ak] = 0;
        ICt[,ak,k] = ICt[,ak,k] + d[[j]][,ak];
        ak = ak + 1;
      } # End of loop on a
    } # End of loop on j
    St[k,] = 1 - colMeans(Qs[[1]]);
    ICt[,,k] = ICt[,,k] + matrix(St[k,], nrow = n, ncol = 2, byrow = TRUE) - Qs[[1]];
  } # End of loop on k
  SE.St = t(apply(ICt, c(2,3), function(x){var(x)/n}))**0.5;
  colnames(St) = colnames(SE.St) = levels(as.factor(dat[,Avar]));
  St = rbind(1, St);
  for(k in 1:ncol(St)){
    for(j in 2:(K+1)){
      if(St[j, k] > St[j-1, k]) St[j,k] = St[j-1, k];
    }
  }
  SE.St = rbind(0, SE.St);

  nlevel = nlevels(as.factor(dat[, Avar]));
  ATE = SE.ATE = LL.ATE = UL.ATE =
    matrix(NA, nrow = K + 1, ncol = nlevel*(nlevel-1)/2);
  column = 1;
  for(i in 1:(nlevel-1)){
    for(j in (i+1):nlevel){
      colnames(ATE)[column] = paste0(levels(as.factor(dat[,Avar]))[i],"-",levels(as.factor(dat[,Avar]))[j]);
      column = column + 1;
    }
  }

  for(r in 1:(K+1)){
    column = 1;
    for(i in 1:(nlevel-1)){
      for(j in (i+1):nlevel){
        ATE[r, column] = St[r, i] - St[r, j];
        if(r!=1) SE.ATE[r, column] = sqrt(var(ICt[, i, r-1] - ICt[, j, r-1])/n);
        column = column + 1;
      }
    }
  }

  LL.ATE = ATE - 1.96*SE.ATE;
  UL.ATE = ATE + 1.96*SE.ATE;


  #### Modeling the hazard (MSM)

  if(!is.null(MSM.form)){
    ## Initializing some objects
    max.k = K*nlevels(as.factor(dat[,Avar]));
    lambda.t = w.t = numeric(max.k);
    X.t = data.frame(matrix(NA, nrow = max.k, ncol = 2));
    k = 1;
    
    ## Building lambda, the weights, and the exposure and time matrix
    for(j in 2:(K+1)){
      ak = 1;
      for(a in sort(unique(dat[,Avar]))){
        lambda.t[k] = (St[j-1, ak] - St[j, ak])/St[j-1, ak];
        w.t[k] = St[j-1, ak];
        X.t[k,] = data.frame(a, j-1);
        ak = ak + 1;
        k = k + 1;
      } # End of loop on a
    } # End of loop on j
  
    
    ## Estimating the MSM parameters
    names(X.t) = c(names(dat[,Avar, drop = FALSE]), "time");
    dat.lambda = data.frame(lambda.t, w.t, X.t);
    X.t = model.matrix(MSM.form, data = dat.lambda);
    mod.lambda = suppressWarnings(glm(lambda.t ~ X.t[,-1],
                                      family = "binomial",
                                      weights = w.t,
                                      data = dat.lambda));
    lambda = coef(mod.lambda);
    names(lambda)[-1] = substr(names(lambda)[-1], 10, 150); # Renames the column of lambda
  
    ## Computing the variance of the MSM parameters
    t1 = t2 = 0;
    nlevelA = nlevels(as.factor(dat[,Avar]));
    k = 1;
    for(j in 1:K){
      for(a in 1:nlevelA){
        t1 = t1 + as.numeric(w.t[k]*exp(t(X.t[k,])%*%lambda)/(1 + exp(t(X.t[k,])%*%lambda))**2)*
                  X.t[k,]%*%t(X.t[k,]);
        if(k + nlevelA <= max.k){
          t2 = t2 + (-X.t[k,] + X.t[k+nlevelA,]%*%solve(1 + exp(t(X.t[k+nlevelA,])%*%lambda)))%*%t(ICt[,a,j]);
        }else{
          t2 = t2 + -X.t[k,]%*%t(ICt[,a,j]);
        }
        k = k + 1;
      }
    }
    IC.lambda = solve(t1)%*%t2;
    Var.lambda = var(t(IC.lambda))/n;
  }
 
  #### Printing results
  results.St = data.frame(S = St, SE.S = SE.St);
  if(!is.null(MSM.form)){
    SE.lambda = sqrt(diag(Var.lambda));
    results.lambda = data.frame(Coef = lambda, 
                                SE = SE.lambda,
                                lower95 = lambda - qnorm(0.975)*SE.lambda,
                                upper95 = lambda + qnorm(0.975)*SE.lambda);
  }

  if(Print){
    cat("\n Estimated survival probabilities:\n ---------------------------------\n");
    print(results.St, digits = 3);
    cat("\n\n");

    if(!is.null(MSM.form)){
      cat("\n Estimated MSM parameters:\n ---------------------------------\n");
      print(results.lambda, digits = 3);
    }
  }
  if(!is.null(MSM.form)){
    invisible(list(St = results.St, MSM = results.lambda, vcov = Var.lambda,
                   ATE = ATE, SE.ATE = SE.ATE, LL.ATE = LL.ATE, UL.ATE = UL.ATE,
                   CV.details = CV.details));
  }else{
    invisible(list(St = results.St,
                   ATE = ATE, SE.ATE = SE.ATE, LL.ATE = LL.ATE, UL.ATE = UL.ATE,
                   CV.details = CV.details, gA = gA, gC = gC));
  }
} # End of function

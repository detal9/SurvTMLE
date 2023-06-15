############################################################
#                                                          #
#  Survival LTMLE with a time-fixed exposure and           #
#  time-varying censoring                                  #
#  V0.6                                                    #
#  date: 2023-06-14                                        #
#  Maintainer: Denis Talbot                                #
#              denis.talbot@fmed.ulaval.ca                 #
#                                                          #
############################################################

require(mlogit);
require(polspline);
require(SuperLearner);


#### Things to be added in future versions:
# How to handle missing data
# Eliminate non-integer #successes in binomial glm warnings
# Add a description and value section 


#### Description
# ...


#### Arguments
# dat:           A dataframe containing the data to be used. 
#                The data must be supplied in a wide format (one row per individual).
# Yvar:          Either a character vector of the names of the Y variables or
#                a numeric vector of the indices of the Y variables in dat.
#                The Y variable must be coded 0/1, where 1 indicates that the event
#                has occured. Once Y_t = 1 for a given subject, the data for times
#                k > t are not used for that subject. The variables must be temporally
#                ordered in Yvar.
# Cvar:          Either a character vector of the names of the C variables or
#                a numeric vector of the indices of the C variables in dat.
#                The C variable must be coded 0/1, where 1 indicates that the follow-up
#                of the subject has been censored. Once C_t = 1 for a given subject,
#                the data for times k > t are not used for that subject. The variables
#                must be temporally ordered in Cvar.
# Avar:          Either the name of the exposure variable or its index in dat.
#                The exposure can either be a binary or a categorical
#                variable. Continuous exposures are not supported. If the exposure has more
#                than two levels, it must be of type factor.
# Lvar:          Either a list of character vectors indicating the names of the time-varying
#                L covariates measured at each time-point or a list of numeric vectors of the
#                indices of the L covariates. The sets of variables (the elements of the list)
#                must be temporally ordered in Lvar. For example, if there are two time points
#                and (L11, L12) are the covariates at time 1 and (L21, L22, L23) are the covariates
#                at time 3, then Lvar = list(c(L11, L12), c(L21, L22, L23)). If there are only
#                baseline covariates, it is possible to repeat the same covariates at all times
#                with lookback = 1.
# L0var:         A character vector of time-fixed L variables measured at baseline or
#                a numeric vector of the indices of the time-fixed L variables in dat (optional).
# lookback:      A numeric value indicating how far back in time should time-varying covariates
#                be considered when modeling the outcome indicators (Y) and the censoring
#                indicators (C). For example, if lookback = 2, Y_t and C_t will be modeled as
#                a function of L_t and L_{t-1}. The default is NULL, indicating that all previous
#                time-varying covariates should be used.
# lookbackY:     See description for lookback. LookbackY indicates the outcome-specific lookback.
#                The default is to use the value provided in lookback.
# lookbackC:     See description for lookback. LookbackC indicates the censoring-specific lookback.
#                The default is to use the value provided in lookback.
# Ymod:          The approach to be used to model the outcome indicators, either "parametric" or
#                "SL" (Super Learner). The default is "parametric". 
#                The choice "parametric" results in using logistic regressions.
# Cmod:          The approach to be used to model the censoring indicators, either "parametric" or
#                "SL" (Super Learner). The default is "parametric".
#                The choice "parametric" results in using logistic regressions.
# Amod:          The approach to be used to model the exposure, either "parametric" or
#                "SL" (Super Learner). The default is "parametric".
#                The choice "parametric" results in using a logistic regression if the exposure
#                is binary and a multinomial regression otherwise. 
#                If the exposure is categorical (>2 levels), Amod = "SL" results in using a
#                polychotomous regression and classification.
# SL.library:    A character vector indicating the learners to be used within the Super Learner.
#                The default is c("SL.glm", "SL.glm.interaction"); The list of all possible
#                learners can be viewed with listWrappers();
# Y.SL.library:  A character vector indicating the learners to be used within the Super Learner
#                for modeling the outcome. The default is the same as SL.library.
# C.SL.library:  A character vector indicating the learners to be used within the Super Learner
#                for modeling the outcome. The default is the same as SL.library.
# A.SL.library:  A character vector indicating the learners to be used within the Super Learner
#                for modeling the exposure. The default is the same as SL.library if A is binary.
#                If A is multilevel (>2 levels) then this option is ignored and
#                polychotomous regression and classification is used.
# gbound:        Bound for the g estimates (A and C). 
#                For exposure A, all predicted probabilities < gbound or > 1 - gbound
#                are truncated. For censoring C, only predicted probabilities 
#                P(C = 0|A, L) < gbound are truncated. Defaults is 0.005.
# V:             Number of folds for the cross-validation when using the Super Learner. 
#                Default is 5.
# MSM.form:      The right hand side of a formula for the MSM relating the hazards 
#                to the exposure and time (optional).
#                The formula can only involve terms for the exposure (with the same name
#                as in dat) and "time", which is coded as time = 1, ..., K. 
#                If the name of Avar in dat is "trt" some examples are ~ trt + time,
#                ~ trt + as.factor(time) or trt + time + trt*time.
# Print:         TRUE or FALSE, whether the function should print the main results (default = TRUE)





#### Value
# ...



#### Details
# Main terms only, but possible to include product terms or polynomial terms or even spline terms
#
# Missing data are not supported. Users may consider using multiple imputation. 


#### Function

surv.TMLE = function(dat, Yvar, Cvar, Avar, Lvar, L0var = NULL,
                     lookback = NULL, lookbackY = lookback, lookbackC = lookback,
                     Ymod = "parametric", Cmod = "parametric", Amod = "parametric",
                     SL.library = c("SL.glm", "SL.glm.interaction"),
                     Y.SL.library = SL.library, C.SL.library = SL.library,
                     A.SL.library = SL.library, gbound = 0.025, V = 5,
                     MSM.form = NULL, Print = TRUE){

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
    if(!is.character(SL.library)) stop("SL.library must a character vector");

    # The following code was extracted/adapted from listWrapper
    everything = sort(getNamespaceExports("SuperLearner"));
    SL.algos = everything[grepl(pattern = "^[S]L", everything)];

    if(min(SL.library %in% SL.algos) == 0) stop("One of the learner supplied in SL.library is invalid");
  }

  ## Verifications for Y.SL.library
  if(Ymod == "SL"){
    if(!is.character(Y.SL.library)) stop("Y.SL.library must a character vector");

    # The following code was extracted/adapted from listWrapper
    everything = sort(getNamespaceExports("SuperLearner"));
    SL.algos = everything[grepl(pattern = "^[S]L", everything)];

    if(min(Y.SL.library %in% SL.algos) == 0) stop("One of the learner supplied in Y.SL.library is invalid");
  }

  ## Verifications for C.SL.library
  if(Cmod == "SL"){
    if(!is.character(C.SL.library)) stop("C.SL.library must a character vector");

    # The following code was extracted/adapted from listWrapper
    everything = sort(getNamespaceExports("SuperLearner"));
    SL.algos = everything[grepl(pattern = "^[S]L", everything)];

    if(min(C.SL.library %in% SL.algos) == 0) stop("One of the learner supplied in C.SL.library is invalid");
  }

  ## Verifications for A.SL.library
  if(Amod == "SL" & nlevels(as.factor(dat[,Avar])) == 2){
    if(!is.character(A.SL.library)) stop("A.SL.library must a character vector");

    # The following code was extracted/adapted from listWrapper
    everything = sort(getNamespaceExports("SuperLearner"));
    SL.algos = everything[grepl(pattern = "^[S]L", everything)];

    if(min(A.SL.library %in% SL.algos) == 0) stop("One of the learner supplied in A.SL.library is invalid");
  }

  ## Verifications for gbound
  if(!is.numeric(gbound)) stop("gbound must be a numeric value between 0 and 1");
  if(length(gbound) > 1) stop("gbound must be a numeric value between 0 and 1");
  if(gbound < 0 | gbound > 1) stop("gbound must be a numeric value between 0 and 1");


  ## Verifications for V
  if(!is.numeric(V)) stop("V must be a numeric value >= 2");
  if(length(V) > 1) stop("V must be a numeric value >= 2");
  if(V < 2) stop("V must be a numeric value >= 2");


  ## Verification for MSM.form
  if(!is.null(MSM.form)){
    if(!inherits(MSM.form, "formula")) stop("When supplied, MSM.form must be a formula");
  }

  ## Verification of Print
  if(!(Print == TRUE | Print == FALSE)) stop("Print must either be TRUE or FALSE");


  ### Compute lower and upper bounds for g estimates
  gbounds = c(min(gbound, 1 - gbound), max(gbound, 1 - gbound)); # Bounds for g estimate


  #### Modeling the exposure
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
      X = as.data.frame(model.matrix(as.formula(Aform), data = dat)[,-1]); # Obtain a matrix to accomodate factors
      names(X) = paste0("X", 1:ncol(X));
      mod.A = SuperLearner(Y = dat[, Avar], X = X, family = "binomial",
                           SL.library = A.SL.library, cvControl = list(V = V));
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
  ## Note :
  # If A is binary gA is a vector of P(A = 1|L0, Lvar[[1]]).
  # If A is categorical gA is a matrix with column j being P(A = j|L0, Lvar[[1]]);


  #### Modeling the censoring
  gC = list();
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
          if(max(dat[,Cvar[j]], na.rm = TRUE) == 0){ # No censoring
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
            # The design matrix
            X = as.data.frame(model.matrix(as.formula(Cform), data = dat)[,-1]);
            names(X) = paste0("X", 1:ncol(X));
            mod.C = SuperLearner(Y = dat[, Cvar[[j]]], X = X, family = "binomial",
                                 SL.library = C.SL.library, cvControl = list(V = V));
            gC[[j]] = pmax(1 - predict(mod.C, OnlySL = TRUE)$pred, gbounds[1]); 
          }
        }else{ # not the first time point
          Cvarj_1 = names(dat[,Cvar[j-1], drop = FALSE]); # previous Cvar
          Yvarj_1 = names(dat[,Yvar[j-1], drop = FALSE]); # previous Yvar
          if(max(dat[,Cvar[j]], na.rm = TRUE) == 0){ # No censoring
            gC[[j]] = rep(NA, n);
            gC[[j]][dat[, Cvarj_1] == 0 & dat[,Yvarj_1] == 0] = 1;
          }else{  # At least some censoring
            gC[[j]] = rep(NA, n);
            used = which(dat[, Cvarj_1] == 0 & dat[, Yvarj_1] == 0);
            X = as.data.frame(model.matrix(as.formula(Cform), data = dat[used,])[,-1]);
            names(X) = paste0("X", 1:ncol(X));
            mod.C = SuperLearner(Y = dat[used, Cvar[[j]]], X = X, family = "binomial",
                                 SL.library = C.SL.library, cvControl = list(V = V));
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
          if(max(dat[,Cvar[j]], na.rm = TRUE) == 0){ # No censoring
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
            # The design matrix
            X = as.data.frame(model.matrix(as.formula(Cform), data = dat)[,-1]);
            names(X) = paste0("X", 1:ncol(X));
            mod.C = SuperLearner(Y = dat[, Cvar[[j]]], X = X, family = "binomial",
                                 SL.library = C.SL.library, cvControl = list(V = V));
            gC[[j]] = pmax(1 - predict(mod.C, OnlySL = TRUE)$pred, gbounds[1]); 
          }
        }else{ # Not the first time point
          Cvarj_1 = names(dat[,Cvar[j-1], drop = FALSE]); # previous Cvar
          Yvarj_1 = names(dat[,Yvar[j-1], drop = FALSE]); # previous Yvar
          if(max(dat[,Cvar[j]], na.rm = TRUE) == 0){ # No censoring
            gC[[j]] = rep(NA, n);
            gC[[j]][dat[, Cvarj_1] == 0 & dat[,Yvarj_1] == 0] = 1;
          }else{  # At least some censoring
            gC[[j]] = rep(NA, n);
            used = which(dat[, Cvarj_1] == 0 & dat[, Yvarj_1] == 0);
            X = as.data.frame(model.matrix(as.formula(Cform), data = dat[used,])[,-1]);
            names(X) = paste0("X", 1:ncol(X));
            mod.C = SuperLearner(Y = dat[used, Cvar[[j]]], X = X, family = "binomial",
                                 SL.library = C.SL.library, cvControl = list(V = V));
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


  #### Modeling the outcome
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
          modQ = SuperLearner(Y = YSL, X = X, family = "binomial",
                              SL.library = Y.SL.library, cvControl = list(V = V));

    
          ## Computing Qj
          if(j == 1){ # First time point
            newdat = dat;
          }else{
            newdat = dat[dat[,Cvar[[j-1]]] == 0,];
          }
          newdat[, Avar] = a;
          if(j == 1){
            Q[[j]][,ak] = predict(modQ, OnlySL = TRUE, newdata = newdat)$pred;
          }else{
            Q[[j]][dat[,Cvar[[j-1]]] == 0,ak] = predict(modQ, OnlySL = TRUE, newdata = newdat)$pred;
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

  ## Initializing some objects
  max.k = K*nlevels(as.factor(dat[,Avar]));
  lambda.t = w.t = numeric(max.k);
  X.t = data.frame(matrix(, nrow = max.k, ncol = 2));
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

 
  #### Printing results
  results.St = data.frame(S = St, SE.S = SE.St);
  SE.lambda = sqrt(diag(Var.lambda));
  results.lambda = data.frame(Coef = lambda, 
                              SE = SE.lambda,
                              lower95 = lambda - qnorm(0.975)*SE.lambda,
                              upper95 = lambda + qnorm(0.975)*SE.lambda);                              

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
                   ATE = ATE, SE.ATE = SE.ATE, LL.ATE = LL.ATE, UL.ATE = UL.ATE));
  }else{
    invisible(list(St = results.St,
                   ATE = ATE, SE.ATE = SE.ATE, LL.ATE = LL.ATE, UL.ATE = UL.ATE));
  }
} # End of function

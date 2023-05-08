# install.packages("mlogit");
# install.packages("SuperLearner");
# install.packages("polspline");
source("...\survTMLE_v0_4.R");


######## Example 1 - Single baseline binary covariate 

#### Function to generate data

generateData <- function(n){
  expit <- plogis;

  ## Generate baseline data
  L <- rbinom(n, size = 1, prob = 0.5);
  A <- rbinom(n, size = 1, prob = plogis(-3 + 0.6*L));
 
  ## Generate counterfactual outcome
  #time 1
  py1.1 <- expit(-2 - 1 + 0.25*L);
  py1.0 <- expit(-2 + 0 + 0.25*L);
  Y1.1 <- rbinom(n, 1, py1.1);
  Y1.0 <- rbinom(n, 1, py1.0);

  #time 2
  Y2.1 <- Y2.0 <- rep(1, n);
  py2.1 <- expit(-2 - 1 + 0.25*L)[Y1.1 == 0];
  py2.0 <- expit(-2 + 0 + 0.25*L)[Y1.0 == 0];
  Y2.1[Y1.1 == 0] <- rbinom(n = length(py2.1), 1, py2.1);
  Y2.0[Y1.0 == 0] <- rbinom(n = length(py2.0), 1, py2.0);

  #time 3
  Y3.1 <- Y3.0 <- rep(1, n);
  py3.1 <- expit(-2 - 1 + 0.25*L)[Y2.1 == 0];
  py3.0 <- expit(-2 + 0 + 0.25*L)[Y2.0 == 0];
  Y3.1[Y2.1 == 0] <- rbinom(n = length(py3.1), 1, py3.1);
  Y3.0[Y2.0 == 0] <- rbinom(n = length(py3.0), 1, py3.0);

  #time 4
  Y4.1 <- Y4.0 <- rep(1, n);
  py4.1 <- expit(-2 - 1 + 0.25*L)[Y3.1 == 0];
  py4.0 <- expit(-2 + 0 + 0.25*L)[Y3.0 == 0];
  Y4.1[Y3.1 == 0] <- rbinom(n = length(py4.1), 1, py4.1);
  Y4.0[Y3.0 == 0] <- rbinom(n = length(py4.0), 1, py4.0);


  ## Generate censoring and observed outcome
  #time 1
  Y1 <- Y1.1*A + Y1.0*(1 - A);
  censor.prob1 <- expit(-5 + 0.2*A + 0.2*L);
  Censor1 <- rbinom(n, 1, censor.prob1);

  #time 2
  Censor2 <- rep(1, n);
  censor.prob2 <- expit(-5 + 0.2*A + 0.2*L)[Censor1 == 0 & Y1 == 0];
  Censor2[Censor1 == 0 & Y1 == 0] <- rbinom(n = length(censor.prob2), 1, censor.prob2);
  Y2 <- Y2.1*A + Y2.0*(1 - A);

  #time 3
  Censor3 <- rep(1, n);
  censor.prob3 <- expit(-5 + 0.2*A + 0.2*L)[Censor2 == 0 & Y2 == 0];
  Censor3[Censor2 == 0 & Y2 == 0] <- rbinom(n = length(censor.prob3), 1,
                                            censor.prob3);
  Y3 <- Y3.1*A + Y3.0*(1 - A);

  #time 4
  Censor4 <- rep(1, n);
  censor.prob4 <- expit(-5 + 0.2*A + 0.2*L)[Censor3==0 & Y3 == 0];
  Censor4[Censor3 == 0 & Y3 == 0] <- rbinom(n = length(censor.prob4), 1,
                                            censor.prob4);
  Y4 <- Y4.1*A + Y4.0*(1 - A);


  ## The observed outcome is missing if censored
  Y1[Censor1==1] <- NA; Y2[Censor2==1] <- NA; Y3[Censor3==1] <- NA; Y4[Censor4==1] <- NA;

  ## Once an event occurred, future Y = 1
  Y2[Y1 == 1] = 1;  Y3[Y2 == 1] = 1;  Y4[Y3 == 1] = 1; 

  ## return counterfactual and observed data
  data.frame(id = 1:n, L, A, Censor1, Censor2, Censor3, Censor4, Y1, Y2, Y3, Y4, Y1.1, Y2.1, Y3.1,
              Y4.1, Y1.0, Y2.0, Y3.0, Y4.0);
}

#### Determine true values 

set.seed(7777);
n = 5000000;
FullData <- generateData(n);

## True survival curves
with(FullData, (1-colMeans(cbind(Y1.1, Y2.1, Y3.1, Y4.1))));
with(FullData, (1-colMeans(cbind(Y1.0, Y2.0, Y3.0, Y4.0))));

# > with(FullData,(1-colMeans(cbind(Y1.1, Y2.1, Y3.1, Y4.1))));
#      Y1.1      Y2.1      Y3.1      Y4.1 
# 0.9459612 0.8951800 0.8471490 0.8017134 
# > with(FullData,(1-colMeans(cbind(Y1.0, Y2.0, Y3.0, Y4.0))));
#      Y1.0      Y2.0      Y3.0      Y4.0 
# 0.8662326 0.7509436 0.6510594 0.5645078 


with(FullData,
     {plot(1-colMeans(cbind(Y1.1, Y2.1, Y3.1, Y4.1)), type = "b", ylim = c(0,1),
           ylab = "S(t)", xlab = "t", main = "True survival curves");
      lines(1-colMeans(cbind(Y1.0, Y2.0, Y3.0, Y4.0)), type = "b", lty = 2, ylim = c(0,1));
      legend(x = 3,y = 1, lty = c(1,2), legend = c("A = 1", "A = 0"), bty = "n")});

trueS <- matrix(c(with(FullData,(1-colMeans(cbind(Y1.0, Y2.0, Y3.0, Y4.0)))),
                  with(FullData,(1-colMeans(cbind(Y1.1, Y2.1, Y3.1, Y4.1))))),
                nrow = 1, byrow = TRUE); 

## True hazard odds ratio
# Create lambda vector
lambdas <- c(1 - trueS[1],                    # lambda_1(1)
             (trueS[1] - trueS[2])/trueS[1],  # lambda_1(2)
             (trueS[2] - trueS[3])/trueS[2],  # lambda_1(3)
             (trueS[3] - trueS[4])/trueS[3],  # lambda_1(4)
             1 - trueS[5],                    # lambda_0(1)
             (trueS[5] - trueS[6])/trueS[5],  # lambda_0(2)
             (trueS[6] - trueS[7])/trueS[6],  # lambda_0(3)
             (trueS[7] - trueS[8])/trueS[7]); # lambda_0(4)

# Create a data frame with lambdas, A, t and weights
data_St_True <- data.frame(lambdas, A = rep(c(0,1), each = 4), Time = rep(1:4,2),
                           weights = c(1, trueS[1:3], 1, trueS[5:7]));

# Fit the hazard model
fit.mod.True = glm(lambdas ~ A + factor(Time),
                   family = "binomial",
                   data = data_St_True,
                   weights = weights); 
trueHR <- coef(fit.mod.True)


##### Estimation
n = 5000;

set.seed(1234);

dat = generateData(n)[,1:11];
results = surv.TMLE(dat = dat, Yvar = c("Y1", "Y2", "Y3", "Y4"),
                    Cvar = c("Censor1", "Censor2", "Censor3", "Censor4"),
                    Avar = "A", Lvar = list("L", "L", "L", "L"),
                    L0var = NULL, lookback = 1,
                    Ymod = "parametric", Cmod = "parametric", Amod = "parametric",
                    gbound = 0.005, V = 5, MSM.form = ~A+as.factor(time), Print = FALSE);

results;



######## Example 2 - Single baseline binary covariate + single time-varying binary covariate

generateData2 <- function(n, ntime = 10){
  expit <- plogis;

  ## Generate baseline data
  L <- rbinom(n, size = 1, prob = 0.5);
  L1 <- rbinom(n, size = 1, prob = 0.3 + 0.4*L);
  A <- rbinom(n, size = 1, prob = plogis(-3 + 0.6*L + 0.6*L1));
  py1.1 <- expit(-2 - 1 + 0.25*L + 0.25*L1);
  py1.0 <- expit(-2 + 0 + 0.25*L + 0.25*L1);
  Y1.1 <- rbinom(n, 1, py1.1);
  Y1.0 <- rbinom(n, 1, py1.0);
  Y1 <- Y1.1*A + Y1.0*(1 - A);
  censor.prob1 <- expit(-5 + 0.2*A + 0.2*L + 0.2*L1);
  Censor1 <- rbinom(n, 1, censor.prob1);

  for(j in 2:ntime){
    Lj_1 = get(paste0("L", j-1));
    Cj_1 = get(paste0("Censor", j-1));
    Yj_1 = get(paste0("Y", j-1));
    nameL = paste0("L", j);
    namepy1 = paste0("py", j, ".1");
    namepy0 = paste0("py", j, ".0");
    nameY1 = paste0("Y", j, ".1");
    nameY0 = paste0("Y", j, ".0");
    nameY = paste0("Y", j);
    nameCprob = paste0("censor.prob", j);
    nameCensor = paste0("Censor", j); 

    assign(nameL, rbinom(n, size = 1, prob = 0.2 + 0.3*L + 0.3*Lj_1));
    Lj = get(nameL);
    assign(namepy1, ifelse(Yj_1 == 0, expit(-2 - 1 + 0.25*L + 0.25*Lj), 1));
    assign(namepy0, ifelse(Yj_1 == 0, expit(-2 + 0 + 0.25*L + 0.25*Lj), 1));
    pyj.1 = get(namepy1);
    pyj.0 = get(namepy0);
    assign(nameY1, rbinom(n, 1, pyj.1));
    assign(nameY0, rbinom(n, 1, pyj.0));
    assign(nameCprob, ifelse(Cj_1 == 0 & Yj_1 == 0, expit(-5 + 0.2*A + 0.2*L + 0.2*Lj), 1));
    pCj = get(nameCprob);
    assign(nameCensor, rbinom(n, 1, pCj));
    assign(nameY, get(nameY1)*A + get(nameY0)*(1 - A));
  }
  id = 1:n;
  var.names = c("id", "L", "A", paste0("Censor", 1:ntime), paste0("Y", 1:ntime),
                paste0("L", 1:ntime), paste0("Y", 1:ntime, ".1"), paste0("Y", 1:ntime, ".0")); 
  ds = matrix(NA, nrow = n, ncol = length(var.names));
  for(j in 1:length(var.names)){
    ds[,j] = get(var.names[j]);
  }
  ds = data.frame(ds);
  names(ds) = var.names;
  return(ds);
}

dat = generateData2(n, ntime = 20);
results = surv.TMLE(dat = dat,
                    Yvar = c("Y1", "Y2", "Y3", "Y4", "Y5", "Y6", "Y7", "Y8", "Y9", "Y10", "Y11", "Y12", "Y13", "Y14", "Y15", "Y16", "Y17", "Y18", "Y19", "Y20"),
                    Cvar = c("Censor1", "Censor2", "Censor3", "Censor4", "Censor5", "Censor6", "Censor7", "Censor8", "Censor9", "Censor10",
                             "Censor11", "Censor12", "Censor13", "Censor14", "Censor15", "Censor16", "Censor17", "Censor18", "Censor19", "Censor20"),
                    Avar = "A",
                    Lvar = list("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10", "L11", "L12", "L13", "L14", "L15", "L16", "L17", "L18", "L19", "L20"),
                    L0var = "L",
                    lookback = 1,
                    Ymod = "parametric", Cmod = "parametric", Amod = "parametric",
                    gbound = 0.005, V = 5, MSM.form = ~A+time, Print = FALSE);



######## Example 3 - Single baseline continuous covariate 
########             + single time-varying continuous covariate
########             with lagged effect and non-linear terms

n = 10000;
generateData3 <- function(n, ntime = 10){
  expit <- plogis;

  ## Generate baseline data
  L <- rnorm(n);
  L1 <- rnorm(n, mean = 0.4*L, sd = sqrt(1 - 0.4**2));
  A <- rbinom(n, size = 1, prob = plogis(-3 + 0.6*L + 0.6*L1));
  py1.1 <- expit(-2 - 1 + 0.25*L + 0.25*L1);
  py1.0 <- expit(-2 + 0 + 0.25*L + 0.25*L1);
  Y1.1 <- rbinom(n, 1, py1.1);
  Y1.0 <- rbinom(n, 1, py1.0);
  Y1 <- Y1.1*A + Y1.0*(1 - A);
  censor.prob1 <- expit(-5 + 0.2*A + 0.2*L + 0.2*L1);
  Censor1 <- rbinom(n, 1, censor.prob1);

  for(j in 2:ntime){
    Lj_1 = get(paste0("L", j-1));
    Cj_1 = get(paste0("Censor", j-1));
    Yj_1 = get(paste0("Y", j-1));
    nameL = paste0("L", j);
    namepy1 = paste0("py", j, ".1");
    namepy0 = paste0("py", j, ".0");
    nameY1 = paste0("Y", j, ".1");
    nameY0 = paste0("Y", j, ".0");
    nameY = paste0("Y", j);
    nameCprob = paste0("censor.prob", j);
    nameCensor = paste0("Censor", j); 

    assign(nameL, rnorm(n, mean = 0.2 + 0.3*L + 0.3*Lj_1, sd = sqrt(1 - 0.3**2 - 0.3**2)));
    Lj = get(nameL);
    assign(namepy1, ifelse(Yj_1 == 0,
                           expit(-2 - 1 + 0.25*L + 0.25*Lj + 0.1*L*Lj - 0.2*L*1 + 0.2*Lj_1 + 0.1*Lj**2),
                           1));
    assign(namepy0, ifelse(Yj_1 == 0, 
                           expit(-2 + 0 + 0.25*L + 0.25*Lj + 0.1*L*Lj - 0.2*L*0 + 0.2*Lj_1 + 0.1*Lj**2),
                           1));
    pyj.1 = get(namepy1);
    pyj.0 = get(namepy0);
    assign(nameY1, rbinom(n, 1, pyj.1));
    assign(nameY0, rbinom(n, 1, pyj.0));
    assign(nameCprob, ifelse(Cj_1 == 0 & Yj_1 == 0,
                             expit(-5 + 0.2*A + 0.2*L + 0.2*Lj + 0.1*Lj_1 + 0.1*Lj**2 - 0.2*A*Lj),
                             1));
    pCj = get(nameCprob);
    assign(nameCensor, rbinom(n, 1, pCj));
    assign(nameY, get(nameY1)*A + get(nameY0)*(1 - A));
  }
  id = 1:n;
  var.names = c("id", "L", "A", paste0("Censor", 1:ntime), paste0("Y", 1:ntime),
                paste0("L", 1:ntime), paste0("Y", 1:ntime, ".1"), paste0("Y", 1:ntime, ".0")); 
  ds = matrix(NA, nrow = n, ncol = length(var.names));
  for(j in 1:length(var.names)){
    ds[,j] = get(var.names[j]);
  }
  ds = data.frame(ds);
  names(ds) = var.names;
  return(ds);
}

dat = generateData3(n, ntime = 12);

# With Super Learner
results = surv.TMLE(dat = dat,
                    Yvar = c("Y1", "Y2", "Y3", "Y4", "Y5", "Y6", "Y7", "Y8", "Y9", "Y10", "Y11", "Y12", "Y13", "Y14", "Y15", "Y16", "Y17", "Y18", "Y19", "Y20"),
                    Cvar = c("Censor1", "Censor2", "Censor3", "Censor4", "Censor5", "Censor6", "Censor7", "Censor8", "Censor9", "Censor10",
                             "Censor11", "Censor12", "Censor13", "Censor14", "Censor15", "Censor16", "Censor17", "Censor18", "Censor19", "Censor20"),
                    Avar = "A",
                    Lvar = list("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10", "L11", "L12", "L13", "L14", "L15", "L16", "L17", "L18", "L19", "L20"),
                    L0var = "L",
                    lookback = 2,
                    Ymod = "SL", Cmod = "SL", Amod = "SL",
                    SL.library = c("SL.glm", "SL.gam", "SL.earth"),
                    gbound = 0.005, V = 5, MSM.form = ~A+time, Print = FALSE);

# With parametric models
results = surv.TMLE(dat = dat,
                    Yvar = c("Y1", "Y2", "Y3", "Y4", "Y5", "Y6", "Y7", "Y8", "Y9", "Y10", "Y11", "Y12", "Y13", "Y14", "Y15", "Y16", "Y17", "Y18", "Y19", "Y20"),
                    Cvar = c("Censor1", "Censor2", "Censor3", "Censor4", "Censor5", "Censor6", "Censor7", "Censor8", "Censor9", "Censor10",
                             "Censor11", "Censor12", "Censor13", "Censor14", "Censor15", "Censor16", "Censor17", "Censor18", "Censor19", "Censor20"),
                    Avar = "A",
                    Lvar = list("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10", "L11", "L12", "L13", "L14", "L15", "L16", "L17", "L18", "L19", "L20"),
                    L0var = "L",
                    lookback = 2,
                    Ymod = "parametric", Cmod = "parametric", Amod = "parametric",
                    gbound = 0.005, V = 5, MSM.form = ~A+time, Print = FALSE);




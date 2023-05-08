# Note: The warnings 'non-integer #successes in a binomial glm!' can be ignored

###### Box 1 -- Function to generate the data

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

  # time 2
  Y2.1 <- Y2.0 <- rep(1, n);
  py2.1 <- expit(-2 - 1 + 0.25*L)[Y1.1 == 0];
  py2.0 <- expit(-2 + 0 + 0.25*L)[Y1.0 == 0];
  Y2.1[Y1.1 == 0] <- rbinom(n = length(py2.1), 1, py2.1);
  Y2.0[Y1.0 == 0] <- rbinom(n = length(py2.0), 1, py2.0);

  # time 3
  Y3.1 <- Y3.0 <- rep(1, n);
  py3.1 <- expit(-2 - 1 + 0.25*L)[Y2.1 == 0];
  py3.0 <- expit(-2 + 0 + 0.25*L)[Y2.0 == 0];
  Y3.1[Y2.1 == 0] <- rbinom(n = length(py3.1), 1, py3.1);
  Y3.0[Y2.0 == 0] <- rbinom(n = length(py3.0), 1, py3.0);

  # time 4
  Y4.1 <- Y4.0 <- rep(1, n);
  py4.1 <- expit(-2 - 1 + 0.25*L)[Y3.1 == 0];
  py4.0 <- expit(-2 + 0 + 0.25*L)[Y3.0 == 0];
  Y4.1[Y3.1 == 0] <- rbinom(n = length(py4.1), 1, py4.1);
  Y4.0[Y3.0 == 0] <- rbinom(n = length(py4.0), 1, py4.0);


  ## Generate censoring and observed outcome
  # time 1
  Y1 <- Y1.1*A + Y1.0*(1 - A);
  censor.prob1 <- expit(-5 + 0.2*A + 0.2*L);
  Censor1 <- rbinom(n, 1, censor.prob1);

  # time 2
  Censor2 <- rep(1, n);
  censor.prob2 <- expit(-5 + 0.2*A + 0.2*L)[Censor1 == 0 & Y1 == 0];
  Censor2[Censor1 == 0 & Y1 == 0] <- rbinom(n = length(censor.prob2), 1, censor.prob2);
  Y2 <- Y2.1*A + Y2.0*(1 - A);

  # time 3
  Censor3 <- rep(1, n);
  censor.prob3 <- expit(-5 + 0.2*A + 0.2*L)[Censor2 == 0 & Y2 == 0];
  Censor3[Censor2 == 0 & Y2 == 0] <- rbinom(n = length(censor.prob3), 1,
                                            censor.prob3);
  Y3 <- Y3.1*A + Y3.0*(1 - A);

  # time 4
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



###### Estimate true values by generating large counterfactual sample
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

trueS <- matrix(c(with(FullData,(1-colMeans(cbind(Y1.1, Y2.1, Y3.1, Y4.1)))),
                  with(FullData,(1-colMeans(cbind(Y1.0, Y2.0, Y3.0, Y4.0))))),
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
data_St_True <- data.frame(lambdas, A = rep(c(1,0), each = 4), Time = rep(1:4,2),
                           weights = c(1, trueS[1:3], 1, trueS[5:7]));

# Fit the hazard model
fit.mod.True = glm(lambdas ~ A + factor(Time),
                   family = "binomial",
                   data = data_St_True,
                   weights = weights); 
trueHR <- coef(fit.mod.True)



###### Generate sample data
set.seed(1234);
n = 5000;
ObsData <- generateData(n)[,1:11];



###### Box 2 -- Estimating $S^a(1)$

## Step 1 -- Outcome model
modQ1 <- glm(Y1 ~ A + L,
             data = subset(ObsData, Censor1==0),
             family = "binomial");
newdat1 <- newdat0 <- ObsData;
newdat1$A <- 1;
newdat0$A <- 0;
Q1_1 <- predict(modQ1, newdata = newdat1, type = "response");
Q1_0 <- predict(modQ1, newdata = newdat0, type = "response");

  
## Step 2 -- Exposure and censoring models  
gA <- glm(A ~ L, family = "binomial", data = ObsData)$fitted;
gC1 <- glm(Censor1 ~ A + L, family = "binomial", data = ObsData)$fitted;


## Step 3 -- Updating the initial estimate
# Computing the clever covariates
ObsData$H1_1 <- (ObsData$A == 1 & ObsData$Censor1 == 0)*(1/(gA*(1-gC1)));
ObsData$H1_0 <- (ObsData$A == 0 & ObsData$Censor1 == 0)*(1/((1 - gA)*(1-gC1)));

# Estimating the update coefficient
mod.update1_1 <- glm(Y1 ~ 1 + offset(qlogis(Q1_1)),
                     weight = H1_1, family = "binomial",
                     data = ObsData);
mod.update1_0 <- glm(Y1 ~ 1 + offset(qlogis(Q1_0)),
                     weight = H1_0, family = "binomial",
                     data = ObsData);
  
# Performing the update
Q1_1.star <- plogis(qlogis(Q1_1) + coef(mod.update1_1));
Q1_0.star <- plogis(qlogis(Q1_0) + coef(mod.update1_0));
S1_1.star <- 1 - mean(Q1_1.star); S1_0.star <- 1 - mean(Q1_0.star);

  
## Step 4 -- Estimate the variance
d1_1 <- with(ObsData, H1_1*(Q1_1.star - Y1));
d1_0 <- with(ObsData, H1_0*(Q1_0.star - Y1));
d1_1[ObsData$H1_1 == 0] = 0; # To deal with NAs in Y1
d1_0[ObsData$H1_0 == 0] = 0; # To deal with NAs in Y1
IC1_1 <- d1_1 + S1_1.star - Q1_1.star;
IC1_0 <- d1_0 + S1_0.star - Q1_0.star;
VarS1_1 <- var(IC1_1)/n;
VarS1_0 <- var(IC1_0)/n;
S1_1.star + c(-1.96, 1.96)*sqrt(VarS1_1);
S1_0.star + c(-1.96, 1.96)*sqrt(VarS1_0);


###### Estimating $S^a(2)$

#### j = t = 2
## Step 1 -- Outcome model
modQ2 <- glm(Y2 ~ A + L,
             data = ObsData, subset = (Censor2 == 0 & Y1 == 0),
             family = "binomial");
newdat1 <- newdat0 <- ObsData[ObsData$Censor1 == 0,];
Q2_1 <- Q2_0 <- rep(NA, n);
newdat1$A <- 1;
newdat0$A <- 0;
Q2_1[ObsData$Censor1 == 0] <- predict(modQ2, newdata = newdat1, type = "response");
Q2_1[ObsData$Y1 == 1] <- 1;
Q2_0[ObsData$Censor1 == 0] <- predict(modQ2, newdata = newdat0, type = "response");
Q2_0[ObsData$Y1 == 1] <- 1;

  
## Step 2 -- Exposure and censoring models  
# The exposure model has already been fitted to estimate S^a(1)
gC2 <- rep(NA,n);
gC2[ObsData$Censor1 == 0 & ObsData$Y1 == 0] <- glm(Censor2 ~ A + L, family = "binomial",
                                                   data = subset(ObsData, Censor1 == 0 & Y1 == 0))$fitted;
  
  
## Step 3 -- Updating the initial estimate
# Computing the clever covariates
ObsData$H2_1 <- (ObsData$A == 1 & ObsData$Censor2 == 0)*(1/(gA*(1 - gC1)*(1 - gC2)));
ObsData$H2_0 <- (ObsData$A == 0 & ObsData$Censor2 == 0)*(1/((1 - gA)*(1 - gC1)*(1 - gC2)));
ObsData$H2_1[ObsData$Censor2 == 1] <- 0; # To deal with NAs in gC2 
ObsData$H2_0[ObsData$Censor2 == 1] <- 0; # To deal with NAs in gC2



# Estimating the update coefficient
mod.update2_1 <- glm(Y2 ~ 1 + offset(qlogis(Q2_1)), weights = H2_1,
                     family = "binomial",
                     data = ObsData, subset = (Y1 == 0));
mod.update2_0 <- glm(Y2 ~ 1 + offset(qlogis(Q2_0)), weights = H2_0,
                     family = "binomial",
                     data = ObsData, subset = (Y1 == 0));

# Performing the update 
Q2_1.star2 <- plogis(qlogis(Q2_1) + coef(mod.update2_1));
Q2_1.star2[ObsData$Y1==1] <- 1;

Q2_0.star2 <- plogis(qlogis(Q2_0) + coef(mod.update2_0));
Q2_0.star2[ObsData$Y1==1] <- 1;

    
## Step 4 -- Estimate the variance component
d2_1_2 <- with(ObsData, H2_1*(Q2_1.star2 - Y2));
d2_0_2 <- with(ObsData, H2_0*(Q2_0.star2 - Y2));
d2_1_2[ObsData$H2_1 == 0] <- 0;
d2_0_2[ObsData$H2_0 == 0] <- 0;
  
  
#### j = 1
## Step 1 -- Outcome model
newdat1 <- newdat0 <- ObsData;
newdat1$A <- 1;
newdat0$A <- 0;
modQ2 <- glm(Q2_1.star2 ~ A + L,
             data = ObsData,
             family = "binomial", subset = (Censor1 == 0));
Q2_1 <- predict(modQ2, newdata = newdat1, type = "response");
  
modQ2 <- glm(Q2_0.star2 ~ A + L,
             data = ObsData,
             family = "binomial", subset = (Censor1 == 0));
Q2_0 <- predict(modQ2, newdata = newdat0, type = "response");
  
  
## Step 2 -- Exposure and censoring models  
# An exposure model was already fitted when estimating S^a(1)
# A model for C_1 was already fitted when estimating S^a(1)
  
  
## Step 3 -- Updating the initial estimate
# Computing the clever covariates was done when estimating S^a(1)
  
# Estimating the update coefficient
mod.update2_1 <- glm(Q2_1.star2 ~ 1 + offset(qlogis(Q2_1)), weights = H1_1,
                     family = "binomial",
                     data = ObsData);
mod.update2_0 <- glm(Q2_0.star2 ~ 1 + offset(qlogis(Q2_0)), weights = H1_0,
                     family = "binomial",
                      data = ObsData);

# Performing the update 
Q2_1.star <- plogis(qlogis(Q2_1) + coef(mod.update2_1));
Q2_0.star <- plogis(qlogis(Q2_0) + coef(mod.update2_0));


## Step 4 -- Estimate the variance component
d2_1_1 <- with(ObsData, H1_1*(Q2_1.star - Q2_1.star2));
d2_0_1 <- with(ObsData, H1_0*(Q2_0.star - Q2_0.star2));
d2_1_1[ObsData$H1_1 == 0] <- 0;
d2_0_1[ObsData$H1_0 == 0] <- 0;
  
  
#### Final estimate and variance
S2_1.star <- 1 - mean(Q2_1.star); S2_0.star <- 1 - mean(Q2_0.star);
IC2_1 <- d2_1_2 + d2_1_1 + S2_1.star - Q2_1.star;
IC2_0 <- d2_0_2 + d2_0_1 + S2_0.star - Q2_0.star;
VarS2_1 <- var(IC2_1)/n;
VarS2_0 <- var(IC2_0)/n;
S2_1.star + c(-1.96, 1.96)*sqrt(VarS2_1);
S2_0.star + c(-1.96, 1.96)*sqrt(VarS2_0); 





###### Estimating $S^a(3)$

#### j = t = 3
## Step 1 -- Outcome model
modQ3 <- glm(Y3 ~ A + L,
             data = ObsData, subset = (Censor3 == 0 & Y2 == 0),
             family = "binomial");
newdat1 <- newdat0 <- ObsData[ObsData$Censor2 == 0,];
Q3_1 <- Q3_0 <- rep(NA, n);
newdat1$A <- 1;
newdat0$A <- 0;
Q3_1[ObsData$Censor2 == 0] <- predict(modQ3, newdata = newdat1, type = "response");
Q3_1[ObsData$Y2 == 1] <- 1;
Q3_0[ObsData$Censor2 == 0] <- predict(modQ3, newdata = newdat0, type = "response");
Q3_0[ObsData$Y2 == 1] <- 1;

  
## Step 2 -- Exposure and censoring models  
# The exposure model has already been fitted
gC3 <- rep(NA,n);
gC3[ObsData$Censor2 == 0 & ObsData$Y2 == 0] <- glm(Censor3 ~ A + L, family = "binomial",
                                                   data = subset(ObsData, Censor2 == 0 & Y2 == 0))$fitted;
  
  
## Step 3 -- Updating the initial estimate
# Computing the clever covariates
ObsData$H3_1 <- (ObsData$A == 1 & ObsData$Censor3 == 0)*(1/(gA*(1 - gC1)*(1 - gC2)*(1 - gC3)));
ObsData$H3_0 <- (ObsData$A == 0 & ObsData$Censor3 == 0)*(1/((1 - gA)*(1 - gC1)*(1 - gC2)*(1 - gC3)));
ObsData$H3_1[ObsData$Censor3 == 1] <- 0; # To deal with NAs in gC3 
ObsData$H3_0[ObsData$Censor3 == 1] <- 0; # To deal with NAs in gC3



# Estimating the update coefficient
mod.update3_1 <- glm(Y3 ~ 1 + offset(qlogis(Q3_1)), weights = H3_1,
                     family = "binomial",
                     data = ObsData, subset = (Y2 == 0));
mod.update3_0 <- glm(Y3 ~ 1 + offset(qlogis(Q3_0)), weights = H3_0,
                     family = "binomial",
                     data = ObsData, subset = (Y2 == 0));

# Performing the update 
Q3_1.star3 <- plogis(qlogis(Q3_1) + coef(mod.update3_1));
Q3_1.star3[ObsData$Y2 == 1] <- 1;

Q3_0.star3 <- plogis(qlogis(Q3_0) + coef(mod.update3_0));
Q3_0.star3[ObsData$Y2 == 1] <- 1;

    
## Step 4 -- Estimate the variance component
d3_1_3 <- with(ObsData, H3_1*(Q3_1.star3 - Y3));
d3_0_3 <- with(ObsData, H3_0*(Q3_0.star3 - Y3));
d3_1_3[ObsData$H3_1 == 0] <- 0;
d3_0_3[ObsData$H3_0 == 0] <- 0;
  


#### j = 2
## Step 1 -- Outcome model
newdat1 <- newdat0 <- ObsData[ObsData$Censor1 == 0,];
newdat1$A <- 1;
newdat0$A <- 0;
modQ3 <- glm(Q3_1.star3 ~ A + L,
             data = ObsData,
             family = "binomial", subset = (Censor2 == 0 & Y1 == 0));
Q3_1[ObsData$Censor1 == 0] <- predict(modQ3, newdata = newdat1, type = "response");
  
modQ3 <- glm(Q3_0.star3 ~ A + L,
             data = ObsData,
             family = "binomial", subset = (Censor2 == 0  & Y1 == 0));
Q3_0[ObsData$Censor1 == 0] <- predict(modQ3, newdata = newdat0, type = "response");
  
  
## Step 2 -- Exposure and censoring models  
# An exposure model was already fitted
# A model for C_2 was already fitted
  
  
## Step 3 -- Updating the initial estimate
# Computing the clever covariates was already done
  
# Estimating the update coefficient
mod.update3_1 <- glm(Q3_1.star3 ~ 1 + offset(qlogis(Q3_1)), weights = H2_1,
                     family = "binomial",
                     data = ObsData);
mod.update3_0 <- glm(Q3_0.star3 ~ 1 + offset(qlogis(Q3_0)), weights = H2_0,
                     family = "binomial",
                     data = ObsData);

# Performing the update 
Q3_1.star2 <- plogis(qlogis(Q3_1) + coef(mod.update3_1));
Q3_1.star2[ObsData$Y1 == 1] <- 1;
Q3_0.star2 <- plogis(qlogis(Q3_0) + coef(mod.update3_0));
Q3_0.star2[ObsData$Y1 == 1] <- 1;



## Step 4 -- Estimate the variance component
d3_1_2 <- with(ObsData, H2_1*(Q3_1.star2 - Q3_1.star3));
d3_0_2 <- with(ObsData, H2_0*(Q3_0.star2 - Q3_0.star3));
d3_1_2[ObsData$H2_1 == 0] <- 0;
d3_0_2[ObsData$H2_0 == 0] <- 0;



#### j = 1
## Step 1 -- Outcome model
newdat1 <- newdat0 <- ObsData;
newdat1$A <- 1;
newdat0$A <- 0;
modQ3 <- glm(Q3_1.star2 ~ A + L,
             data = ObsData,
             family = "binomial", subset = (Censor1 == 0));
Q3_1 <- predict(modQ3, newdata = newdat1, type = "response");
  
modQ3 <- glm(Q3_0.star2 ~ A + L,
             data = ObsData,
             family = "binomial", subset = (Censor1 == 0));
Q3_0 <- predict(modQ3, newdata = newdat0, type = "response");
   
  
## Step 2 -- Exposure and censoring models  
# An exposure model was already fitted
# A model for C_1 was already fitted
  
  
## Step 3 -- Updating the initial estimate
# Computing the clever covariates was already done
  
# Estimating the update coefficient
mod.update3_1 <- glm(Q3_1.star2 ~ 1 + offset(qlogis(Q3_1)), weights = H1_1,
                     family = "binomial",
                     data = ObsData);
mod.update3_0 <- glm(Q3_0.star2 ~ 1 + offset(qlogis(Q3_0)), weights = H1_0,
                     family = "binomial",
                      data = ObsData);

# Performing the update 
Q3_1.star <- plogis(qlogis(Q3_1) + coef(mod.update3_1));
Q3_0.star <- plogis(qlogis(Q3_0) + coef(mod.update3_0));


## Step 4 -- Estimate the variance component
d3_1_1 <- with(ObsData, H1_1*(Q3_1.star - Q3_1.star2));
d3_0_1 <- with(ObsData, H1_0*(Q3_0.star - Q3_0.star2));
d3_1_1[ObsData$H1_1 == 0] <- 0;
d3_0_1[ObsData$H1_0 == 0] <- 0;

  
#### Final estimate and variance
S3_1.star <- 1 - mean(Q3_1.star); S3_0.star <- 1 - mean(Q3_0.star);
IC3_1 <- d3_1_3 + d3_1_2 + d3_1_1 + S3_1.star - Q3_1.star;
IC3_0 <- d3_0_3 + d3_0_2 + d3_0_1 + S3_0.star - Q3_0.star;
VarS3_1 <- var(IC3_1)/n;
VarS3_0 <- var(IC3_0)/n;
S3_1.star + c(-1.96, 1.96)*sqrt(VarS3_1);
S3_0.star + c(-1.96, 1.96)*sqrt(VarS3_0); 




###### Estimating $S^a(4)$

#### j = t = 4
## Step 1 -- Outcome model
Q4_1 <- Q4_0 <- rep(NA, n);
modQ4 <- glm(Y4 ~ A + L,
             data = ObsData, subset = (Censor4 == 0 & Y3 == 0),
             family = "binomial");
newdat1 <- newdat0 <- ObsData[ObsData$Censor3 == 0,];
newdat1$A <- 1;
newdat0$A <- 0;
Q4_1[ObsData$Censor3 == 0] <- predict(modQ4, newdata = newdat1, type = "response");
Q4_1[ObsData$Y3 == 1] <- 1;
Q4_0[ObsData$Censor3 == 0] <- predict(modQ3, newdata = newdat0, type = "response");
Q4_0[ObsData$Y3 == 1] <- 1;

  
## Step 2 -- Exposure and censoring models  
# The exposure model has already been fitted
gC4 <- rep(NA,n);
gC4[ObsData$Censor3 == 0 & ObsData$Y3 == 0] <- glm(Censor4 ~ A + L, family = "binomial",
                                                   data = subset(ObsData, Censor3 == 0 & Y3 == 0))$fitted;
  
  
## Step 3 -- Updating the initial estimate
# Computing the clever covariates
ObsData$H4_1 <- (ObsData$A == 1 & ObsData$Censor4 == 0)*(1/(gA*(1 - gC1)*(1 - gC2)*(1 - gC3)*(1 - gC4)));
ObsData$H4_0 <- (ObsData$A == 0 & ObsData$Censor4 == 0)*(1/((1 - gA)*(1 - gC1)*(1 - gC2)*(1 - gC3)*(1 - gC4)));
ObsData$H4_1[ObsData$Censor4 == 1] <- 0; # To deal with NAs in gC4 
ObsData$H4_0[ObsData$Censor4 == 1] <- 0; # To deal with NAs in gC4



# Estimating the update coefficient
mod.update4_1 <- glm(Y4 ~ 1 + offset(qlogis(Q4_1)), weights = H4_1,
                     family = "binomial",
                     data = ObsData, subset = (Y3 == 0));
mod.update4_0 <- glm(Y4 ~ 1 + offset(qlogis(Q4_0)), weights = H4_0,
                     family = "binomial",
                     data = ObsData, subset = (Y3 == 0));

# Performing the update 
Q4_1.star4 <- plogis(qlogis(Q4_1) + coef(mod.update4_1));
Q4_1.star4[ObsData$Y3 == 1] <- 1;

Q4_0.star4 <- plogis(qlogis(Q4_0) + coef(mod.update4_0));
Q4_0.star4[ObsData$Y3 == 1] <- 1;

    
## Step 4 -- Estimate the variance component
d4_1_4 <- with(ObsData, H4_1*(Q4_1.star4 - Y4));
d4_0_4 <- with(ObsData, H4_0*(Q4_0.star4 - Y4));
d4_1_4[ObsData$H4_1 == 0] <- 0;
d4_0_4[ObsData$H4_0 == 0] <- 0;
  

#### j = 3
## Step 1 -- Outcome model
newdat1 <- newdat0 <- ObsData[ObsData$Censor2 == 0,];
newdat1$A <- 1;
newdat0$A <- 0;
modQ4 <- glm(Q4_1.star4 ~ A + L,
             data = ObsData, subset = (Censor3 == 0 & Y2 == 0),
             family = "binomial");
Q4_1[ObsData$Censor2 == 0] <- predict(modQ4, newdata = newdat1, type = "response");

modQ4 <- glm(Q4_0.star4 ~ A + L,
             data = ObsData, subset = (Censor3 == 0 & Y2 == 0),
             family = "binomial");
Q4_0[ObsData$Censor2 == 0] <- predict(modQ4, newdata = newdat1, type = "response");



  
## Step 2 -- Exposure and censoring models  
# The exposure model has already been fitted
# The censoring model has already been fitted
  
  
## Step 3 -- Updating the initial estimate
# The clever covariates have already been computed


# Estimating the update coefficient
mod.update4_1 <- glm(Q4_1.star4 ~ 1 + offset(qlogis(Q4_1)), weights = H3_1,
                     family = "binomial",
                     data = ObsData, subset = (Y2 == 0));
mod.update4_0 <- glm(Q4_0.star4 ~ 1 + offset(qlogis(Q4_0)), weights = H3_0,
                     family = "binomial",
                     data = ObsData, subset = (Y2 == 0));

# Performing the update 
Q4_1.star3 <- plogis(qlogis(Q4_1) + coef(mod.update4_1));
Q4_1.star3[ObsData$Y2 == 1] <- 1;

Q4_0.star3 <- plogis(qlogis(Q4_0) + coef(mod.update4_0));
Q4_0.star3[ObsData$Y2 == 1] <- 1;

    
## Step 4 -- Estimate the variance component
d4_1_3 <- with(ObsData, H3_1*(Q4_1.star3 - Q4_1.star4));
d4_0_3 <- with(ObsData, H3_0*(Q4_0.star3 - Q4_0.star4));
d4_1_3[ObsData$H3_1 == 0] <- 0;
d4_0_3[ObsData$H3_0 == 0] <- 0;
  


#### j = 2
## Step 1 -- Outcome model
newdat1 <- newdat0 <- ObsData[ObsData$Censor1 == 0,];
newdat1$A <- 1;
newdat0$A <- 0;
modQ4 <- glm(Q4_1.star3 ~ A + L,
             data = ObsData,
             family = "binomial", subset = (Censor2 == 0 & Y1 == 0));
Q4_1[ObsData$Censor1 == 0] <- predict(modQ4, newdata = newdat1, type = "response");
  
modQ4 <- glm(Q4_0.star3 ~ A + L,
             data = ObsData,
             family = "binomial", subset = (Censor2 == 0  & Y1 == 0));
Q4_0[ObsData$Censor1 == 0] <- predict(modQ4, newdata = newdat0, type = "response");

  
  
## Step 2 -- Exposure and censoring models  
# An exposure model was already fitted
# A model for C_2 was already fitted
  
  
## Step 3 -- Updating the initial estimate
# Computing the clever covariates was already done
  
# Estimating the update coefficient
mod.update4_1 <- glm(Q4_1.star3 ~ 1 + offset(qlogis(Q4_1)), weights = H2_1,
                     family = "binomial",
                     data = ObsData);
mod.update4_0 <- glm(Q4_0.star3 ~ 1 + offset(qlogis(Q4_0)), weights = H2_0,
                     family = "binomial",
                     data = ObsData);

# Performing the update 
Q4_1.star2 <- plogis(qlogis(Q4_1) + coef(mod.update4_1));
Q4_1.star2[ObsData$Y1 == 1] <- 1;
Q4_0.star2 <- plogis(qlogis(Q4_0) + coef(mod.update4_0));
Q4_0.star2[ObsData$Y1 == 1] <- 1;



## Step 4 -- Estimate the variance component
d4_1_2 <- with(ObsData, H2_1*(Q4_1.star2 - Q4_1.star3));
d4_0_2 <- with(ObsData, H2_0*(Q4_0.star2 - Q4_0.star3));
d4_1_2[ObsData$H2_1 == 0] <- 0;
d4_0_2[ObsData$H2_0 == 0] <- 0;



#### j = 1
## Step 1 -- Outcome model
newdat1 <- newdat0 <- ObsData;
newdat1$A <- 1;
newdat0$A <- 0;
modQ4 <- glm(Q4_1.star2 ~ A + L,
             data = ObsData,
             family = "binomial", subset = (Censor1 == 0));
Q4_1 <- predict(modQ4, newdata = newdat1, type = "response");
  
modQ4 <- glm(Q4_0.star2 ~ A + L,
             data = ObsData,
             family = "binomial", subset = (Censor1 == 0));
Q4_0 <- predict(modQ4, newdata = newdat0, type = "response");

  
  
## Step 2 -- Exposure and censoring models  
# An exposure model was already fitted
# A model for C_1 was already fitted
  
  
## Step 3 -- Updating the initial estimate
# Computing the clever covariates was already done
  
# Estimating the update coefficient
mod.update4_1 <- glm(Q4_1.star2 ~ 1 + offset(qlogis(Q4_1)), weights = H1_1,
                     family = "binomial",
                     data = ObsData);
mod.update4_0 <- glm(Q4_0.star2 ~ 1 + offset(qlogis(Q4_0)), weights = H1_0,
                     family = "binomial",
                      data = ObsData);

# Performing the update 
Q4_1.star <- plogis(qlogis(Q4_1) + coef(mod.update4_1));
Q4_0.star <- plogis(qlogis(Q4_0) + coef(mod.update4_0));


## Step 4 -- Estimate the variance component
d4_1_1 <- with(ObsData, H1_1*(Q4_1.star - Q4_1.star2));
d4_0_1 <- with(ObsData, H1_0*(Q4_0.star - Q4_0.star2));
d4_1_1[ObsData$H1_1 == 0] <- 0;
d4_0_1[ObsData$H1_0 == 0] <- 0;

  
#### Final estimate and variance
S4_1.star <- 1 - mean(Q4_1.star); S4_0.star <- 1 - mean(Q4_0.star);
IC4_1 <- d4_1_4 + d4_1_3 + d4_1_2 + d4_1_1 + S4_1.star - Q4_1.star;
IC4_0 <- d4_0_4 + d4_0_3 + d4_0_2 + d4_0_1 + S4_0.star - Q4_0.star;
VarS4_1 <- var(IC4_1)/n;
VarS4_0 <- var(IC4_0)/n;
S4_1.star + c(-1.96, 1.96)*sqrt(VarS4_1);
S4_0.star + c(-1.96, 1.96)*sqrt(VarS4_0); 





###### Box 4 -- Estimating the hazard model 
lambdas <- c(1 - S1_1.star,
             (S1_1.star - S2_1.star)/S1_1.star,
             (S2_1.star - S3_1.star)/S2_1.star,
             (S3_1.star - S4_1.star)/S3_1.star,
             1 - S1_0.star,
             (S1_0.star - S2_0.star)/S1_0.star,
             (S2_0.star - S3_0.star)/S2_0.star,
             (S3_0.star - S4_0.star)/S3_0.star);
lambdas[lambdas < 0] <- 0;

data_St <- data.frame(lambdas = lambdas,
                      A = rep(c(1,0),each = 4),
                      Time = rep(1:4,2),
                      weights = c(1, S1_1.star, S2_1.star, S3_1.star,
                                  1, S1_0.star, S2_0.star, S3_0.star));


## Logistic regression
fit.mod <- glm(lambdas ~ A + factor(Time), family = "binomial", weights = weights, data = data_St);
B <- matrix(coef(fit.mod),nrow=5)
X <- matrix(model.matrix(fit.mod), nrow = 8, byrow = FALSE);

## Efficient influence curve
# Objects for the rows of X
X1.1 <- t(X[1,, drop = FALSE]);
X1.0 <- t(X[5,, drop = FALSE]);
X2.1 <- t(X[2,, drop = FALSE]);
X2.0 <- t(X[6,, drop = FALSE]);
X3.1 <- t(X[3,, drop = FALSE]);
X3.0 <- t(X[7,, drop = FALSE]);
X4.1 <- t(X[4,, drop = FALSE]);
X4.0 <- t(X[8,, drop = FALSE]);

#Time 1
term1_1_t1 <- as.numeric(exp(t(X1.1)%*%B)/(1 + exp(t(X1.1)%*%B))**2)*X1.1%*%t(X1.1);
term1_0_t1 <- as.numeric(exp(t(X1.0)%*%B)/(1 + exp(t(X1.0)%*%B))**2)*X1.0%*%t(X1.0);
term1_t1 <-  term1_0_t1 + term1_1_t1
term2_1_t1 <- (-X1.1 + X2.1%*%(1 + exp(t(X2.1)%*%B))**-1)%*%t(IC1_1);
term2_0_t1 <- (-X1.0 + X2.0%*%(1 + exp(t(X2.0)%*%B))**-1)%*%t(IC1_0);
term2_t1 <- term2_1_t1 + term2_0_t1;

#Time 2
term1_1_t2 <- (1 - mean(Q1_1.star))*as.numeric(exp(t(X2.1)%*%B)/(1 + exp(t(X2.1)%*%B))**2)*
               X2.1%*%t(X2.1);
term1_0_t2 <- (1 - mean(Q1_0.star))*as.numeric(exp(t(X2.0)%*%B)/(1 + exp(t(X2.0)%*%B))**2)*
              X2.0%*%t(X2.0);
term1_t2 <- term1_0_t2 + term1_1_t2
term2_1_t2 <- (-X2.1 + X3.1%*%(1 + exp(t(X3.1)%*%B))**-1)%*%t(IC2_1);
term2_0_t2 <- (-X2.0 + X3.0%*%(1 + exp(t(X3.0)%*%B))**-1)%*%t(IC2_0);
term2_t2 <- term2_1_t2 + term2_0_t2;

#Time 3
term1_1_t3 <- (1 - mean(Q2_1.star))*as.numeric(exp(t(X3.1)%*%B)/(1 + exp(t(X3.1)%*%B))**2)*
               X3.1%*%t(X3.1);
term1_0_t3 <- (1 - mean(Q2_0.star))*as.numeric(exp(t(X3.0)%*%B)/(1 + exp(t(X3.0)%*%B))**2)*
               X3.0%*%t(X3.0);
term1_t3 <- term1_0_t3 + term1_1_t3
term2_1_t3 <- (-X3.1 + X4.1%*%(1 + exp(t(X4.1)%*%B))**-1)%*%t(IC3_1);
term2_0_t3 <- (-X3.0 + X4.0%*%(1 + exp(t(X4.0)%*%B))**-1)%*%t(IC3_0);
term2_t3 <- term2_1_t3 + term2_0_t3;

#Time 4
term1_1_t4 <- (1 - mean(Q3_1.star))*as.numeric(exp(t(X4.1)%*%B)/(1 + exp(t(X4.1)%*%B))**2)*
               X4.1%*%t(X4.1);
term1_0_t4 <- (1 - mean(Q3_0.star))*as.numeric(exp(t(X4.0)%*%B)/(1 + exp(t(X4.0)%*%B))**2)*
               X4.0%*%t(X4.0);
term1_t4 <- term1_0_t4 + term1_1_t4
term2_1_t4 <- (-X4.1)%*%t(IC4_1);
term2_0_t4 <- (-X4.0)%*%t(IC4_0);
term2_t4 <- term2_1_t4 + term2_0_t4;

#Variance
var.hr <- var(t(solve(term1_t1 + term1_t2 + term1_t3 + term1_t4)%*%
              (term2_t1 + term2_t2 + term2_t3 + term2_t4)))/n; 
cbind(B, B - 1.96*sqrt(diag(var.hr)), B + 1.96*sqrt(diag(var.hr)));







library(Rlab)
library(dplyr)

# set seed
set.seed(123)

lbb <- function(iter, n, p, mean, sd, con, 
                bL, bY, or1, or2, or3, or4, or5, nU)
  
{
  ## iter   number of iterations
  ## n      sample size for each of the generated datasets
  ## p      prevalence of U 
  ## mean   mean of normally distributed exposure
  ## sd     standard deviation of normally distributed exposure
  ## con    exposure contrast (i.e., increment of exposure for associations)
  ## bL     baseline odds of loss
  ## bY     baseline odds of outcome
  ## or1    odds ratio for exposure-loss effect
  ## or2    odds ratio for U-loss effect
  ## or3    odds ratio for exposure*U loss interaction effect
  ## or4    odds ratio for U-outcome effect
  ## or5    odds ratio for exposure-outcome effect
  ## nU    number of Us
  
  # create matrix to store simulation results
  results.cs <- rep(NA, iter)
  results.dos <- rep(NA, iter)
  results.both <- rep(NA, iter)
  
  for (i in 1:iter){
    set.seed(i)
    
    # monitor simulation
    if((i %% 10) == 0) print(c(iter, i))
    
    # specify contrast
    mean_std <- mean/con
    sd_std <- sd/con
    
    ############################################
    ### Mechanism 1: collider-stratification ###
    ############################################
    
    # data generating mechanism 
    cs <- data.frame("id" = 1:n) %>%
      mutate(no2 = rnorm(n, mean_std, sd_std), # exposure is gaussian
             U = rbern(n, p), 
             prob_loss = plogis(log(bL) + log(or1)*no2 + nU*log(or2)*U),
             loss = rbern(n, prob_loss),
             pY = plogis(log(bY) + nU*log(or4)*U + log(or5)*no2),
             Y = rbern(n, pY))
    
    # fit a logistic model among live births
    model1 <- cs %>%
      glm(formula = Y ~ no2,
          family = binomial(link = "logit"),
          data = .,
          subset = loss==0)
    
    results.cs[i] <- exp(model1$coef[2])
    
    ##############################################
    ### Mechanism 2: depletion of susceptibles ###
    ##############################################
    
    # data generating mechanism
    dos <- data.frame("id" = 1:n) %>%
      mutate(no2 = rnorm(n, mean_std, sd_std),
             U = rbern(n, p),
             prob_loss = plogis(log(bL) + nU*log(or3)*no2*U),
             loss = rbern(n, prob_loss),
             pY = plogis(log(bY) + nU*log(or4)*U + log(or5)*no2),
             Y = rbern(n, pY))
    
    # fit a logistic model among live births
    model2 <- dos %>%
      glm(formula = Y ~ no2,
          family = binomial(link = "logit"),
          data = .,
          subset = loss==0)
    
    results.dos[i] <- exp(model2$coef[2])
    
    #######################
    ### Both Mechanisms ###
    #######################
    
    # data generating mechanism
    both <- data.frame("id" = 1:n) %>%
      mutate(no2 = rnorm(n, mean_std, sd_std),
             U = rbern(n, p),
             prob_loss = plogis(log(bL) + log(or1)*no2 + nU*log(or2)*U + nU*log(or3)*no2*U),
             loss = rbern(n, prob_loss),
             pY = plogis(log(bY) + nU*log(or4)*U + log(or5)*no2),
             Y = rbern(n, pY))
    
    # fit a logistic model among live births
    model3 <- both %>%
      glm(formula = Y ~ no2,
          family = binomial(link = "logit"),
          data = .,
          subset = loss==0)
    
    results.both[i] <- exp(model3$coef[2])
    
  }
  
  # turn results into data frame
  simResults <- data.frame(results.cs, results.dos, results.both) 
  
  # get mean ORs for ASD-no2 association
  csOR <- mean(simResults$results.cs)
  dosOR <- mean(simResults$results.dos)
  bothOR <- mean(simResults$results.both)
  
  # get the 95% simulation intervals
  csOR.si <- quantile(simResults$results.cs, c(0.025,0.975))
  dosOR.si <- quantile(simResults$results.dos, c(0.025,0.975))
  bothOR.si <- quantile(simResults$results.both, c(0.025,0.975))
  
  # put results together
  lbb.cs <- paste(formatC(csOR,digits=2,format="f"),
                  " (",formatC(csOR.si[1],digits=2,format="f"),
                  ",",formatC(csOR.si[2],digits=2,format="f"),")",sep="")
  lbb.dos <- paste(formatC(dosOR,digits=2,format="f"),
                   " (",formatC(dosOR.si[1],digits=2,format="f"),
                   ",",formatC(dosOR.si[2],digits=2,format="f"),")",sep="")
  lbb.both <- paste(formatC(bothOR,digits=2,format="f"),
                    " (",formatC(bothOR.si[1],digits=2,format="f"),
                    ",",formatC(bothOR.si[2],digits=2,format="f"),")",sep="")
  lbb.results <- cbind(lbb.cs,lbb.dos,lbb.both)
  lbb.results
    
}

lbb(iter=1000, n=100000, p=0.75, mean=16.7, sd=4.3, con=5.85,
    bL=0.05, bY=0.015, or1=3, or2=3, or3=3, or4=3, or5=1, nU=1)



##############################################################
# Purity function Q(alpha) in Eq 17
# Input: alpha
# Need to change name: 
#   cv_train: the data set that used to calculate purity
#             the groups of the data are named as 'pbo', 'drg'
##############################################################

alphaKL = function(alpha){
  
  # make sure alpha has norm 1
  alpha = alpha / sqrt(sum(alpha^2))
  alpha = matrix(alpha, p, 1)
  
  # get the data sets of the two groups: 'pbo' and 'drg'
  datatemp = cv_train
  datatemp$W = covarlist%*%alpha
  dat_pbo_est = datatemp[datatemp$group == 'pbo', ]
  dat_drg_est = datatemp[datatemp$group == 'drg', ]
  
  # fit mixed effect model with slope and concavity 
  # fixed and random effects 
  fit_pbo_est = myTryCatch(lmer(Y ~ week + I(week^2) +
                                  W + W * week +
                                  W * I(week^2) + 
                                  (week+I(week^2)|subj),
                                data = dat_pbo_est, REML = FALSE))
  fit_drg_est = myTryCatch(lmer(Y ~ week + I(week^2) + 
                                  W + W * week +
                                  W * I(week^2) + 
                                  (week+I(week^2)|subj),
                                data = dat_drg_est, REML = FALSE))
  
  # pull out the fixed effect coefficients
  # pull out the covariance matrix from random effects
  beta1 = as.matrix(fixef(fit_drg_est$value))[2:3]
  gamma1 = as.matrix(fixef(fit_drg_est$value))[5:6] 
  D1 = as.matrix(VarCorr(fit_drg_est$value)$subj)[2:3, 2:3] 
  beta2 = as.matrix(fixef(fit_pbo_est$value))[2:3] 
  gamma2 = as.matrix(fixef(fit_pbo_est$value))[5:6] 
  D2 = as.matrix(VarCorr(fit_pbo_est$value)$subj)[2:3, 2:3] 
  
  # calcualte a1, a2, a3
  a1 = - dim(D1)[1] + 0.5*sum(diag(ginv(D1)%*%D2)) + 
    0.5*sum(diag(ginv(D2)%*%D1)) +
    0.5*t(beta1 - beta2)%*%(ginv(D1) + ginv(D2))%*%(beta1 - beta2)
  a2 = t(beta1 - beta2)%*%(ginv(D1) + ginv(D2))%*%(gamma1 - gamma2)
  a3 = 0.5*t(gamma1 - gamma2)%*%(ginv(D1) + ginv(D2))%*%(gamma1 - gamma2)
  
  # calculate the purity
  res = a1 + a2 * c(t(estmux) %*% alpha) + a3 * c(
    t(alpha) %*% (estmux %*% t(estmux)+ estsigmax) %*% alpha)
  return(-res)
}

##############################################################
# myTryCatch 
##############################################################

myTryCatch <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error=function(e) {
      err <<- e
      NULL
    }), warning=function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    })
  list(value=value, warning=warn, error=err)
}

##############################################################
# Data generation function 
# n: number of subject in one group
# p: dimension of baseline covariates 
# beta.drg, beta.pbo: fixed effects coefficient
# gamma.drg, gamma.pbo: fixed effects coefficient
# D.drg, D.pbo: covariance matrix of random effects 
# e.drg, e.pbo: random error
# mux, sigmax: the 
##############################################################

datageneration = function(
  n, p, alpha, 
  beta.drg, beta.pbo, 
  gamma.drg, gamma.pbo, 
  D.drg, D.pbo,
  e.drg, e.pbo, mux, sigmax
){
  
  ### simulation of baseline covariates
  n0 = n 
  n = 2 * n
  X = mvrnorm(n, mux, sigmax)
  w0 = X %*% (alpha)
  b.drg = mvrnorm(n, rep(0,3), D.drg)
  b.pbo = mvrnorm(n, rep(0,3), D.pbo)
  
  ### simulation of outcomes with only fixed effects 
  Y.drg.fix = XX %*% t(matrix(beta.drg, n, 3, byrow = TRUE) +
                         w0 %x% t(gamma.drg)) 
  Y.drg.fix = matrix(Y.drg.fix, n*nt, 1) 
  Y.pbo.fix = XX %*% t(matrix(beta.pbo, n, 3, byrow = TRUE) +
                         w0 %x% t(gamma.pbo)) 
  Y.pbo.fix = matrix(Y.pbo.fix, n*nt, 1) 
  ### simulation of outcomes with fixed and random effects
  Y.drg.ran = XX %*% t(matrix(beta.drg, n, 3, byrow = TRUE) + b.drg + 
                         w0 %x% t(gamma.drg)) 
  Y.drg.ran = matrix(Y.drg.ran, n*nt, 1)
  Y.pbo.ran = XX %*% t(matrix(beta.pbo, n, 3, byrow = TRUE) + b.pbo + 
                         w0 %x% t(gamma.pbo)) 
  Y.pbo.ran = matrix(Y.pbo.ran, n*nt, 1)
  ### simulation of outcomes with fixed, random effects and random error.
  Y.drg = Y.drg.ran + rnorm(length(Y.drg.ran), 0, e.drg)
  Y.pbo = Y.pbo.ran + rnorm(length(Y.pbo.ran), 0, e.pbo)
  
  # data set
  data = data.frame( Y.drg = Y.drg, Y.pbo = Y.pbo, 
                     Y.drg.fix = Y.drg.fix, Y.pbo.fix = Y.pbo.fix,
                     Y.drg.ran = Y.drg.ran, Y.pbo.ran = Y.pbo.ran,
                     week = rep(t, n), 
                     subj = rep(1:n, each = nt),
                     w0 = rep(w0, each = nt),
                     matrix(rep(X, each = nt), nt*n, p), 
                     b.drg = matrix(rep(b.drg, each = nt), nt*n, 3),
                     b.pbo = matrix(rep(b.pbo, each = nt), nt*n, 3))
  
  data$group = c(rep('drg', nt*n0), rep('pbo', nt*n0))
  
  ### The truegroup assignment 
  data = data %>% group_by(subj) %>% 
    mutate(last.drg = Y.drg[length(Y.drg)],
           first.drg = Y.drg[1],
           last.drg.fix = Y.drg.fix[length(Y.drg.fix)],
           first.drg.fix = Y.drg.fix[1],
           last.drg.ran = Y.drg.ran[length(Y.drg)],
           first.drg.ran = Y.drg.ran[1],
           last.pbo = Y.pbo[length(Y.drg)],
           first.pbo = Y.pbo[1],
           last.pbo.fix = Y.pbo.fix[length(Y.drg)],
           first.pbo.fix = Y.pbo.fix[1],
           last.pbo.ran = Y.pbo.ran[length(Y.drg)],
           first.pbo.ran = Y.pbo.ran[1]) %>% 
    mutate(cs.drg = last.drg - first.drg,
           cs.pbo = last.pbo - first.pbo, 
           cs.drg.fix = last.drg.fix - first.drg.fix,
           cs.pbo.fix = last.pbo.fix - first.pbo.fix,
           cs.drg.ran = last.drg.ran - first.drg.ran,
           cs.pbo.ran = last.pbo.ran - first.pbo.ran)
  data = data.frame(data)
  data$truegroup.cs = ifelse(data$cs.drg > data$cs.pbo, 'drg', 'pbo')
  data$truegroup.fix = ifelse(data$cs.drg.fix > data$cs.pbo.fix, 'drg', 'pbo')
  data$truegroup.ran = ifelse(data$cs.drg.ran > data$cs.pbo.ran, 'drg', 'pbo')
  
  ### The observed outcome
  data$Y = ifelse(data$group == 'drg', data$Y.drg, data$Y.pbo)
  
  ### dataset with baseline characteristics and change score
  uni.data = unique(data[,c('subj','w0',
                            paste("X", 1:p, sep=''),
                            'cs.drg','cs.pbo',
                            'cs.drg.fix', 'cs.pbo.fix',
                            'cs.drg.ran', 'cs.pbo.ran',
                            'group',
                            'truegroup.cs','truegroup.ran','truegroup.fix')])
  rownames(uni.data)= NULL
  
  return(list(data = data, uni.data = uni.data))
}

##############################################################
# Proportion of correct decision calculation
# Input: alpha, the training dataset called data, and the 
#        testing dataset called data.test
# The groups are named as 'pbo', 'drg'
# The matrix of baseline covariates called: covarlist
##############################################################

# pcd calculation
pcdfun = function(alpha, data, data.test){
  
  # calculate the biosignatures for each subject
  datatemp = data
  datatemp$W = covarlist%*%alpha
  dat_pbo_est = datatemp[datatemp$group == 'pbo', ]
  dat_drg_est = datatemp[datatemp$group == 'drg', ]
  
  # fit mixed effect model with linear and quadaratic terms of time
  fit_pbo_est = myTryCatch(lmer(Y ~ week + I(week^2) + W + W * week +
                                  W * I(week^2) + (week+I(week^2)|subj),
                                data = dat_pbo_est, REML = FALSE))
  fit_drg_est = myTryCatch(lmer(Y ~ week + I(week^2) + W + W * week +
                                  W * I(week^2) + (week+I(week^2)|subj),
                                data = dat_drg_est, REML = FALSE))
  
  # pull out the fixed effects 
  beta1 = as.matrix(fixef(fit_drg_est$value))[1:3]
  gamma1 = as.matrix(fixef(fit_drg_est$value))[4:6] 
  beta2 = as.matrix(fixef(fit_pbo_est$value))[1:3] 
  gamma2 = as.matrix(fixef(fit_pbo_est$value))[4:6] 
  gdiff = (cbind(1, t, t^2)[length(t),] - cbind(1, t, t^2)[1,])/ (max(t) - min(t))
  gdiff = matrix(gdiff, 3, 1)
  
  # calculate the average tangent slope
  estgroup = c()
  for(i in unique(data.test$subj)){
    temp = c(matrix(unlist(data.test[i, paste('X', 1:p, sep = '')]),1,p) %*% alpha)
    # print(temp)
    if (t(gdiff) %*% (beta1 + temp * gamma1) > t(gdiff) %*% (beta2 + temp * gamma2)){
      estgroup = c(estgroup,'drg')
    }else{
      estgroup = c(estgroup, 'pbo')
    }
  }
  
  # calculate the proportion of correct decision
  result = list(cs = sum(estgroup == data.test$truegroup.cs) / length(estgroup), 
                fix = sum(estgroup == data.test$truegroup.fix) / length(estgroup),
                ran = sum(estgroup == data.test$truegroup.ran) / length(estgroup))
  
  return(result)
}

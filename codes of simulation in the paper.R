
##################################
# load in functions
##################################

# change the file path
# setwd('')
source('OWL_codes_DTRLearn.R')
source('functions in the paper.R')

##################################
# parameters 
##################################

# p = (select one value from 2, 10, 20, 30)
# theta_angle = (select one value from 0,1,2,5)

n = 100
mux = p:1; mux[1:floor(p/2)] = -mux[1:floor(p/2)]
sigmax = matrix(0.5,p,p); diag(sigmax) = 1
t = 0:7; nt = length(t)
XX = cbind(1, t, t^2)
alpha = matrix(1:p, p, 1)
alpha = alpha / sqrt(sum(alpha^2))
D.drg = matrix(-0.1, 3, 3); diag(D.drg) = c(0.5,0.5,0.01)
D.drg[2,3] = D.drg[3,2] = -0.01
D.drg[1,3] = D.drg[3,1] = -0.01
D.pbo = matrix(-0.12, 3, 3); diag(D.pbo) = c(0.4,0.5,0.01)
D.pbo[2,3] = D.pbo[3,2] = -0.01
D.pbo[1,3] = D.pbo[3,1] = -0.01
e.drg = e.pbo = 1
beta.drg = matrix(c(0,3,-0.5),3,1)
beta.pbo = matrix(c(0,2.3,-0.4),3,1)
gamma.pbo = matrix(c(0, cos(theta_angle/ 180 * pi), 
                     sin(theta_angle/ 180 * pi)),3,1)
gamma.drg = matrix(c(0, cos(theta_angle/ 180 * pi), 
                     -sin(theta_angle/ 180 * pi)),3,1)
covar_names = paste('X', 1:p, sep='')


##################################
# repetition 200 times
# when there is no missing data
##################################

names = paste('simulationrdata.',p,
              '.a.',theta_angle,
              '.RData', sep='')

alphalingem = lingem = c()
alphasimml = simml = c()
alphakl = resultalpha = c()
owlpcd = c()

for(iters in 1:200){
  
  set.seed(iters)
  
  # generate training data
  tmp = datageneration(n, p, alpha, beta.drg, beta.pbo, 
                       gamma.drg, gamma.pbo, D.drg, D.pbo,
                       e.drg, e.pbo, mux, sigmax)
  data = tmp$data
  uni.data = tmp$uni.data
  uni.data$cs = ifelse(uni.data$group == 'drg', 
                       uni.data$cs.drg, uni.data$cs.pbo)
  cv_train = data
  covarlist = matrix(
    unlist(data[,paste('X',1:p,sep='')]), dim(data)[1], p)
  estmux = matrix(apply(unique(covarlist), 2, mean), p, 1)
  estsigmax = cov(unique(covarlist))
  table(uni.data$group, uni.data$truegroup.cs)
  
  # generate test data
  tmp = datageneration(1000, p, alpha, beta.drg, beta.pbo, gamma.drg, gamma.pbo, D.drg, D.pbo,
                       e.drg, e.pbo, mux, sigmax)
  
  data.test = tmp$data
  data.test = tmp$uni.data
  
  #######################################
  # Different methods for PCD calculation
  #######################################
  
  ########################
  # LS-KLD
  ########################
  
  # use optim function to find the 
  # solution of alpha that maximize the purity function 
  reskl = optim(rep(1,p), alphaKL, 
                control = list(reltol = 1e-2))
  
  # save the alpha vector
  alphakl = rbind(alphakl, 
                  c(reskl$par / sqrt(sum(reskl$par^2))))
  colnames(alphakl) = paste('a',1:p,sep='')
  
  # save the PCD
  tmptrue = pcdfun(alpha, data, data.test)
  tmpkl = pcdfun(reskl$par, data, data.test)
  resultalpha = rbind(resultalpha, 
                      c(tmptrue$cs, tmpkl$cs, 
                        tmptrue$fix, tmpkl$fix,
                        tmptrue$ran, tmpkl$ran,
                        c(cosine(c(alpha), c(reskl$par)))))
  colnames(resultalpha) = c(
    'truecs','klcs',
    'truefix','klfix',
    'trueran','klran',
    'cosinekl'
  )
  
  ########################
  # Linear GEM
  ########################
  
  # linear model 
  lmform = 'cs ~ -1 + X1'
  for(i in paste('X', 2:p, sep = '')){
    lmform = paste(lmform,'+',i)
  }
  lmform = formula(lmform)
  
  # fit linear model
  l_est_drg = lm(lmform, data = uni.data[uni.data$group == 'drg', ])
  l_est_pbo = lm(lmform, data = uni.data[uni.data$group == 'pbo', ])
  # test results
  test_drg = predict(l_est_drg, newdat = data.test[,paste('X', 1:p, sep='')])
  test_pbo = predict(l_est_pbo, newdat = data.test[,paste('X', 1:p, sep='')])
  lin.group = ifelse(test_drg > test_pbo, 'drg', 'pbo')
  # save the coefficent values
  beta1 = l_est_drg$coefficients
  beta2 = l_est_pbo$coefficients
  alphatemp = (beta1 - beta2) / as.numeric(sqrt(t(beta1 - beta2) %*% estsigmax %*% (beta1 - beta2)))
  alphatemp = alphatemp / sqrt(sum(alphatemp^2))
  alphatemp = matrix(alphatemp, length(alphatemp), 1)
  alphalingem = rbind(alphalingem, c(alphatemp))
  colnames(alphalingem) = paste('a',1:p,sep='')
  
  # save the results 
  lingem = rbind(lingem, c(
    sum(lin.group == data.test$truegroup.cs) / dim(data.test)[1],
    sum(lin.group == data.test$truegroup.fix)/ dim(data.test)[1],
    sum(lin.group == data.test$truegroup.ran)/ dim(data.test)[1],
    cosine(c(alpha), c(alphatemp))))
  colnames(lingem) = c('cs','fix','ran','cosine')
  
  ########################
  # SIMML 
  # used the code in simml library
  ########################
  
  yy1 = uni.data$cs
  trt1 = ifelse(uni.data$group == 'drg', 1, 2)
  xx1 = matrix(unlist(uni.data[, paste('X',1:p, sep ='')]), 
               dim(uni.data)[1], p)
  xx1.test = matrix(unlist(data.test[,paste('X',1:p, sep ='')]),
                    dim(data.test )[1], p)
  
  # fit simml
  glm.fit = glm(yy1 ~ xx1, family = "gaussian")
  mu.hat = as.vector(predict(glm.fit, newX = xx1, type="link"))
  simml.obj = simml(yy1, trt1, xx1, Xm = mu.hat, family = "gaussian")
  simml.trt.rule = pred.simml(simml.obj, 
                               newX = xx1.test)$trt.rule
  simml.trt.rule = ifelse(simml.trt.rule == 1, 'drg', 'pbo')
 
  # save the alpha coefficient
  temp = simml.obj$beta.coef
  temp = temp/sqrt(sum(temp^2))
  
  # save the PCD
  simml = rbind(simml, c(
    sum(simml.trt.rule == data.test$truegroup.cs) / dim(data.test)[1],
    sum(simml.trt.rule == data.test$truegroup.fix)/ dim(data.test)[1],
    sum(simml.trt.rule == data.test$truegroup.ran)/ dim(data.test)[1],
    cosine(c(alpha), c(temp))
  ))
  
  alphasimml = rbind(alphasimml, c(temp))
  colnames(alphasimml) = paste('a',1:p,sep='')
  colnames(simml) = c('cs','fix','ran','cosine')
  
  ########################
  # Outcome Weighted Learning (OWL)
  ########################
  
  tmp1 = data.frame(data %>% group_by(subj) %>% slice(1))
  tmp2 = data.frame(data %>% group_by(subj) %>% slice(n()))
  tmp = tmp1[,c(paste('X',1:p,sep=''),'truegroup.cs','truegroup.fix','truegroup.ran' )]
  tmp$cs = tmp2$Y - tmp1$Y
  head(tmp)
  tmptr = ifelse(tmp$truegroup.ran == 'drg', 1, -1)
  tmpy = tmp$cs
  tmpx = tmp[, paste('X',1:p,sep='')]
  tmpx = as.matrix(tmpx)
  
  testlabel = ifelse(data.test$truegroup.ran == 'drg', 1, -1)
  xtest = data.test[,paste('X', 1:p, sep='')]
  xtest = as.matrix(xtest)
  
  # fit the owl learning with gaussian kernal
  owl.gaussian.fit = Olearning_Single(tmpx, tmptr, tmpy, 
                                      kernel ="rbf", m=5,
                                      clinear=2.^(-2:2),
                                       sigma = c(0.01, 0.04, 0.08, 0.32, 1.28),
                                       ridge.reg=TRUE)
  trt.rule.owl.gaussian =  predict(owl.gaussian.fit, xtest)
  
  owlpcd = c(owlpcd, sum(diag(table(testlabel, trt.rule.owl.gaussian))) / 
            sum(table(testlabel, trt.rule.owl.gaussian)))
  
  result = list(lingem = lingem, simml = simml, 
                owlpcd = owlpcd, 
                alphalingem = alphalingem, 
                alphasimml = alphasimml, 
                alphakl = alphakl, 
                resultalpha = resultalpha)
 
  save(result, file = names)
}




##################################
# repetition 200 times
# when there is missing data
# MCAR 30%
##################################

names = paste('simulationrdata.mcar.',p,
              '.a.',theta_angle,
              '.RData', sep='')

alphalingem = lingem = c()
alphasimml = simml = c()
alphakl = resultalpha = c()
owlpcd = c()

for(iters in 1:200){
  
  set.seed(iters)
  
  # generate training data
  tmp = datageneration(n, p, alpha, beta.drg, beta.pbo, 
                       gamma.drg, gamma.pbo, D.drg, D.pbo,
                       e.drg, e.pbo, mux, sigmax)
  
  data = tmp$data
  
  mis = c()
  for(misssss in 1:(2*n)){
    mis = c(mis, 1, 1, sample(c(0,1),6, replace = TRUE))
  }
  mis = ifelse(mis == 0, NA, mis)
  data$mis = mis
  data = na.omit(data)
  
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
  data$Y = ifelse(data$group == 'drg', data$Y.drg, data$Y.pbo)
  
  uni.data = unique(data[,c('subj','w0',
                            paste("X", 1:p, sep=''),
                            'cs.drg','cs.pbo',
                            'cs.drg.fix', 'cs.pbo.fix',
                            'cs.drg.ran', 'cs.pbo.ran',
                            'group',
                            'truegroup.cs','truegroup.ran','truegroup.fix')])
  rownames(uni.data)= NULL
  
  uni.data$cs = ifelse(uni.data$group == 'drg', 
                       uni.data$cs.drg, uni.data$cs.pbo)
  
  cv_train = data
  covarlist = matrix(unlist(data[,paste('X',1:p,sep='')]), dim(data)[1], p)
  estmux = matrix(apply(unique(covarlist), 2, mean), p, 1)
  estsigmax = cov(unique(covarlist))
  
  # generate test data
  tmp = datageneration(1000, p, alpha, beta.drg, beta.pbo, gamma.drg, gamma.pbo, D.drg, D.pbo,
                       e.drg, e.pbo, mux, sigmax)
  
  data.test = tmp$data
  data.test = tmp$uni.data
  
  #######################################
  # Different methods for PCD calculation
  #######################################
  
  ########################
  # LS-KLD
  ########################
  
  # use optim function to find the 
  # solution of alpha that maximize the purity function 
  reskl = optim(rep(1,p), alphaKL, 
                control = list(reltol = 1e-2))
  
  # save the alpha vector
  alphakl = rbind(alphakl, 
                  c(reskl$par / sqrt(sum(reskl$par^2))))
  colnames(alphakl) = paste('a',1:p,sep='')
  
  # save the PCD
  tmptrue = pcdfun(alpha, data, data.test)
  tmpkl = pcdfun(reskl$par, data, data.test)
  resultalpha = rbind(resultalpha, 
                      c(tmptrue$cs, tmpkl$cs, 
                        tmptrue$fix, tmpkl$fix,
                        tmptrue$ran, tmpkl$ran,
                        c(cosine(c(alpha), c(reskl$par)))))
  colnames(resultalpha) = c(
    'truecs','klcs',
    'truefix','klfix',
    'trueran','klran',
    'cosinekl'
  )
  
  ########################
  # Linear GEM
  ########################
  
  # linear model 
  lmform = 'cs ~ -1 + X1'
  for(i in paste('X', 2:p, sep = '')){
    lmform = paste(lmform,'+',i)
  }
  lmform = formula(lmform)
  
  # fit linear model
  l_est_drg = lm(lmform, data = uni.data[uni.data$group == 'drg', ])
  l_est_pbo = lm(lmform, data = uni.data[uni.data$group == 'pbo', ])
  # test results
  test_drg = predict(l_est_drg, newdat = data.test[,paste('X', 1:p, sep='')])
  test_pbo = predict(l_est_pbo, newdat = data.test[,paste('X', 1:p, sep='')])
  lin.group = ifelse(test_drg > test_pbo, 'drg', 'pbo')
  # save the coefficent values
  beta1 = l_est_drg$coefficients
  beta2 = l_est_pbo$coefficients
  alphatemp = (beta1 - beta2) / as.numeric(sqrt(t(beta1 - beta2) %*% estsigmax %*% (beta1 - beta2)))
  alphatemp = alphatemp / sqrt(sum(alphatemp^2))
  alphatemp = matrix(alphatemp, length(alphatemp), 1)
  alphalingem = rbind(alphalingem, c(alphatemp))
  colnames(alphalingem) = paste('a',1:p,sep='')
  
  # save the results 
  lingem = rbind(lingem, c(
    sum(lin.group == data.test$truegroup.cs) / dim(data.test)[1],
    sum(lin.group == data.test$truegroup.fix)/ dim(data.test)[1],
    sum(lin.group == data.test$truegroup.ran)/ dim(data.test)[1],
    cosine(c(alpha), c(alphatemp))))
  colnames(lingem) = c('cs','fix','ran','cosine')
  
  ########################
  # SIMML 
  # used the code in simml library
  ########################
  
  yy1 = uni.data$cs
  trt1 = ifelse(uni.data$group == 'drg', 1, 2)
  xx1 = matrix(unlist(uni.data[, paste('X',1:p, sep ='')]), 
               dim(uni.data)[1], p)
  xx1.test = matrix(unlist(data.test[,paste('X',1:p, sep ='')]),
                    dim(data.test )[1], p)
  
  # fit simml
  glm.fit = glm(yy1 ~ xx1, family = "gaussian")
  mu.hat = as.vector(predict(glm.fit, newX = xx1, type="link"))
  simml.obj = simml(yy1, trt1, xx1, Xm = mu.hat, family = "gaussian")
  simml.trt.rule = pred.simml(simml.obj, 
                              newX = xx1.test)$trt.rule
  simml.trt.rule = ifelse(simml.trt.rule == 1, 'drg', 'pbo')
  
  # save the alpha coefficient
  temp = simml.obj$beta.coef
  temp = temp/sqrt(sum(temp^2))
  
  # save the PCD
  simml = rbind(simml, c(
    sum(simml.trt.rule == data.test$truegroup.cs) / dim(data.test)[1],
    sum(simml.trt.rule == data.test$truegroup.fix)/ dim(data.test)[1],
    sum(simml.trt.rule == data.test$truegroup.ran)/ dim(data.test)[1],
    cosine(c(alpha), c(temp))
  ))
  
  alphasimml = rbind(alphasimml, c(temp))
  colnames(alphasimml) = paste('a',1:p,sep='')
  colnames(simml) = c('cs','fix','ran','cosine')
  
  ########################
  # Outcome Weighted Learning (OWL)
  ########################
  
  tmp1 = data.frame(data %>% group_by(subj) %>% slice(1))
  tmp2 = data.frame(data %>% group_by(subj) %>% slice(n()))
  tmp = tmp1[,c(paste('X',1:p,sep=''),'truegroup.cs','truegroup.fix','truegroup.ran' )]
  tmp$cs = tmp2$Y - tmp1$Y
  head(tmp)
  tmptr = ifelse(tmp$truegroup.ran == 'drg', 1, -1)
  tmpy = tmp$cs
  tmpx = tmp[, paste('X',1:p,sep='')]
  tmpx = as.matrix(tmpx)
  
  testlabel = ifelse(data.test$truegroup.ran == 'drg', 1, -1)
  xtest = data.test[,paste('X', 1:p, sep='')]
  xtest = as.matrix(xtest)
  
  # fit the owl learning with gaussian kernal
  owl.gaussian.fit = Olearning_Single(tmpx, tmptr, tmpy, 
                                      kernel ="rbf", m=5,
                                      clinear=2.^(-2:2),
                                      sigma = c(0.01, 0.04, 0.08, 0.32, 1.28),
                                      ridge.reg=TRUE)
  trt.rule.owl.gaussian =  predict(owl.gaussian.fit, xtest)
  
  owlpcd = c(owlpcd, sum(diag(table(testlabel, trt.rule.owl.gaussian))) / 
               sum(table(testlabel, trt.rule.owl.gaussian)))
  
  result = list(lingem = lingem, simml = simml, 
                owlpcd = owlpcd, 
                alphalingem = alphalingem, 
                alphasimml = alphasimml, 
                alphakl = alphakl, 
                resultalpha = resultalpha)
  
  save(result, file = names)
}


##################################
# repetition 200 times
# when there is missing data
# Dropout
##################################

names = paste('simulationrdata.drop.',p,
              '.a.',theta_angle,
              '.RData', sep='')

alphalingem = lingem = c()
alphasimml = simml = c()
alphakl = resultalpha = c()
owlpcd = c()

for(iters in 1:200){
  
  set.seed(iters)
  
  # generate training data
  tmp = datageneration(n, p, alpha, beta.drg, beta.pbo, 
                       gamma.drg, gamma.pbo, D.drg, D.pbo,
                       e.drg, e.pbo, mux, sigmax)
  mis = c()
  for(misssss in 1:(2*n)){
    if(sample(0:1)[1] == 0){
      mis = c(mis, rep(1,2), rep(NA, 6))
    }else{
      mis = c(mis, rep(1,8))
    }
  }
  # mis = ifelse(mis == 0, NA, mis)
  data$mis = mis
  data = na.omit(data)
  
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
  data$Y = ifelse(data$group == 'drg', data$Y.drg, data$Y.pbo)
  
  uni.data = unique(data[,c('subj','w0',
                            paste("X", 1:p, sep=''),
                            'cs.drg','cs.pbo',
                            'cs.drg.fix', 'cs.pbo.fix',
                            'cs.drg.ran', 'cs.pbo.ran',
                            'group',
                            'truegroup.cs','truegroup.ran','truegroup.fix')])
  rownames(uni.data)= NULL
  
  uni.data$cs = ifelse(uni.data$group == 'drg', 
                       uni.data$cs.drg, uni.data$cs.pbo)
  
  cv_train = data
  covarlist = matrix(unlist(data[,paste('X',1:p,sep='')]), dim(data)[1], p)
  estmux = matrix(apply(unique(covarlist), 2, mean), p, 1)
  estsigmax = cov(unique(covarlist))
  table(uni.data$group, uni.data$truegroup.cs)
  
  # generate test data
  tmp = datageneration(1000, p, alpha, beta.drg, beta.pbo, gamma.drg, gamma.pbo, D.drg, D.pbo,
                       e.drg, e.pbo, mux, sigmax)
  
  data.test = tmp$data
  data.test = tmp$uni.data
  
  #######################################
  # Different methods for PCD calculation
  #######################################
  
  ########################
  # LS-KLD
  ########################
  
  # use optim function to find the 
  # solution of alpha that maximize the purity function 
  reskl = optim(rep(1,p), alphaKL, 
                control = list(reltol = 1e-2))
  
  # save the alpha vector
  alphakl = rbind(alphakl, 
                  c(reskl$par / sqrt(sum(reskl$par^2))))
  colnames(alphakl) = paste('a',1:p,sep='')
  
  # save the PCD
  tmptrue = pcdfun(alpha, data, data.test)
  tmpkl = pcdfun(reskl$par, data, data.test)
  resultalpha = rbind(resultalpha, 
                      c(tmptrue$cs, tmpkl$cs, 
                        tmptrue$fix, tmpkl$fix,
                        tmptrue$ran, tmpkl$ran,
                        c(cosine(c(alpha), c(reskl$par)))))
  colnames(resultalpha) = c(
    'truecs','klcs',
    'truefix','klfix',
    'trueran','klran',
    'cosinekl'
  )
  
  ########################
  # Linear GEM
  ########################
  
  # linear model 
  lmform = 'cs ~ -1 + X1'
  for(i in paste('X', 2:p, sep = '')){
    lmform = paste(lmform,'+',i)
  }
  lmform = formula(lmform)
  
  # fit linear model
  l_est_drg = lm(lmform, data = uni.data[uni.data$group == 'drg', ])
  l_est_pbo = lm(lmform, data = uni.data[uni.data$group == 'pbo', ])
  # test results
  test_drg = predict(l_est_drg, newdat = data.test[,paste('X', 1:p, sep='')])
  test_pbo = predict(l_est_pbo, newdat = data.test[,paste('X', 1:p, sep='')])
  lin.group = ifelse(test_drg > test_pbo, 'drg', 'pbo')
  # save the coefficent values
  beta1 = l_est_drg$coefficients
  beta2 = l_est_pbo$coefficients
  alphatemp = (beta1 - beta2) / as.numeric(sqrt(t(beta1 - beta2) %*% estsigmax %*% (beta1 - beta2)))
  alphatemp = alphatemp / sqrt(sum(alphatemp^2))
  alphatemp = matrix(alphatemp, length(alphatemp), 1)
  alphalingem = rbind(alphalingem, c(alphatemp))
  colnames(alphalingem) = paste('a',1:p,sep='')
  
  # save the results 
  lingem = rbind(lingem, c(
    sum(lin.group == data.test$truegroup.cs) / dim(data.test)[1],
    sum(lin.group == data.test$truegroup.fix)/ dim(data.test)[1],
    sum(lin.group == data.test$truegroup.ran)/ dim(data.test)[1],
    cosine(c(alpha), c(alphatemp))))
  colnames(lingem) = c('cs','fix','ran','cosine')
  
  ########################
  # SIMML 
  # used the code in simml library
  ########################
  
  yy1 = uni.data$cs
  trt1 = ifelse(uni.data$group == 'drg', 1, 2)
  xx1 = matrix(unlist(uni.data[, paste('X',1:p, sep ='')]), 
               dim(uni.data)[1], p)
  xx1.test = matrix(unlist(data.test[,paste('X',1:p, sep ='')]),
                    dim(data.test )[1], p)
  
  # fit simml
  glm.fit = glm(yy1 ~ xx1, family = "gaussian")
  mu.hat = as.vector(predict(glm.fit, newX = xx1, type="link"))
  simml.obj = simml(yy1, trt1, xx1, Xm = mu.hat, family = "gaussian")
  simml.trt.rule = pred.simml(simml.obj, 
                              newX = xx1.test)$trt.rule
  simml.trt.rule = ifelse(simml.trt.rule == 1, 'drg', 'pbo')
  
  # save the alpha coefficient
  temp = simml.obj$beta.coef
  temp = temp/sqrt(sum(temp^2))
  
  # save the PCD
  simml = rbind(simml, c(
    sum(simml.trt.rule == data.test$truegroup.cs) / dim(data.test)[1],
    sum(simml.trt.rule == data.test$truegroup.fix)/ dim(data.test)[1],
    sum(simml.trt.rule == data.test$truegroup.ran)/ dim(data.test)[1],
    cosine(c(alpha), c(temp))
  ))
  
  alphasimml = rbind(alphasimml, c(temp))
  colnames(alphasimml) = paste('a',1:p,sep='')
  colnames(simml) = c('cs','fix','ran','cosine')
  
  ########################
  # Outcome Weighted Learning (OWL)
  ########################
  
  tmp1 = data.frame(data %>% group_by(subj) %>% slice(1))
  tmp2 = data.frame(data %>% group_by(subj) %>% slice(n()))
  tmp = tmp1[,c(paste('X',1:p,sep=''),'truegroup.cs','truegroup.fix','truegroup.ran' )]
  tmp$cs = tmp2$Y - tmp1$Y
  head(tmp)
  tmptr = ifelse(tmp$truegroup.ran == 'drg', 1, -1)
  tmpy = tmp$cs
  tmpx = tmp[, paste('X',1:p,sep='')]
  tmpx = as.matrix(tmpx)
  
  testlabel = ifelse(data.test$truegroup.ran == 'drg', 1, -1)
  xtest = data.test[,paste('X', 1:p, sep='')]
  xtest = as.matrix(xtest)
  
  # fit the owl learning with gaussian kernal
  owl.gaussian.fit = Olearning_Single(tmpx, tmptr, tmpy, 
                                      kernel ="rbf", m=5,
                                      clinear=2.^(-2:2),
                                      sigma = c(0.01, 0.04, 0.08, 0.32, 1.28),
                                      ridge.reg=TRUE)
  trt.rule.owl.gaussian =  predict(owl.gaussian.fit, xtest)
  
  owlpcd = c(owlpcd, sum(diag(table(testlabel, trt.rule.owl.gaussian))) / 
               sum(table(testlabel, trt.rule.owl.gaussian)))
  
  result = list(lingem = lingem, simml = simml, 
                owlpcd = owlpcd, 
                alphalingem = alphalingem, 
                alphasimml = alphasimml, 
                alphakl = alphakl, 
                resultalpha = resultalpha)
  
  save(result, file = names)
}


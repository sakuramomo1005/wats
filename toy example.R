#################################
# A toy example to illustrate
# the prediction results from 
# LS-KLD method
#################################

# change the file path
# setwd('')
source('OWL_codes_DTRLearn.R')
source('functions in the paper.R')

library(dplyr)
library(MASS)
library(lme4)
library(lsa)
library(simml)

# parameters
p = 2
theta_angle = 5
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

# data generation
# generate training data
set.seed(123)
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

# generate test data
tmp = datageneration(1000, p, alpha, beta.drg, beta.pbo, gamma.drg, gamma.pbo, D.drg, D.pbo,
                     e.drg, e.pbo, mux, sigmax)
data.test = tmp$data
data.test = tmp$uni.data

reskl = optim(rep(1,p), alphaKL, 
              control = list(reltol = 1e-2))

# estimate the proportion of correct decision
# true alpha
pcdfun(alpha, data, data.test)
# estimate the proportion of correct decision
# estimated alpha using LS-KLD
pcdfun(reskl$par, data, data.test)
# cosine similarity between true alpha and the 
# alpha estimated by LS-KLD
cosine(c(reskl$par), c(alpha))











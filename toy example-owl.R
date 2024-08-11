setwd('/Users/lanqiuyao/Desktop/NYU/Research/FDA.2021/47. code in paper/supply')
source("OWL_codes_DTRLearn.R")
source("csim_codes.R")
library(glmnet)

set.seed(123)
n = 250
p = 10 
delta = 1
xi = 0
w = 1
correlationX = 0

# a training dataset
data.train <- generate.data(n= n, p=p, delta = delta, xi = xi, sigma=0.5, correlationX=correlationX, w=w)
data.train$SNR  # the ratio of the interactions ("signal") vs. the main effects ("noise") in the canonical parameter
Tr.train <- data.train$Tr
y.train  <- data.train$y
X.train <- data.train$X

table(Tr.train) # 1 2
apply(X.train, 2, mean)
apply(X.train, 2, sd)
range(y.train)
# a (large) testing dataset
data.test <- generate.data(n=10^3, p=p, delta = delta, xi = xi, sigma=0.5, correlationX=correlationX, w=w)
Tr.test <- data.test$Tr
y.test  <- data.test$y
X.test <- data.test$X
data.test$value.opt  # the optimal "Value"
data.test$SNR
potental.outcome <- data.test$potental.outcome 

print("the outcome weighted learning tuning parameter selection; takes few minutes")
system.time(
  owl.fit.tuning <- Olearning_Single(X.train, (-1)^Tr.train, y.train,
                                     clinear=2.^(-2:2),
                                     sigma = c(0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28),
                                     kernel ="rbf", m=5)
)
bestC <- owl.fit.tuning$bestC
bestSig <- owl.fit.tuning$bestSig

bestC = 4
bestSig = 0.08
# 4. Outcome weighted learning (residualized)
cvfit=cv.glmnet(X.train, y.train, nfolds=5)
co=as.matrix(predict(cvfit,s="lambda.min",type="coeff"))
y.train.r <- y.train - cbind(1, X.train)%*%co

# a) OWL-linear
owl.linear.fit   <-  wsvm(X.train, (-1)^Tr.train, y.train.r, C=bestC, e=1e-5, ridge.reg=TRUE)
# apply the fitted classifiers to the testset
trt.rule.owl.linear  <- predict(owl.linear.fit, X.test)
# estimate the value (computed from the testset)
trt.rule <- trt.rule.owl.linear/2 + 1.5
# Value estimation
tmp <- rep(0,length(Tr.test))
for(i in 1:length(Tr.test)) tmp[i] <- potental.outcome[i, trt.rule[i]]
owl.linear.value <- mean(tmp)
owl.linear.value
#  0.7616705
sum(trt.rule == Tr.test) / length(Tr.test) 
# 0.519

# b) OWL-gaussian
owl.gaussian.fit <- wsvm(X.train,
                         (-1)^Tr.train, 
                         y.train.r, 'rbf', 
                         bestSig, C= 1,
                         e=1e-5, ridge.reg=TRUE)
# apply the fitted classifiers to the testset
trt.rule.owl.gaussian <-  predict(owl.gaussian.fit, X.test)
# estimate the value (computed from the testset)
trt.rule <- trt.rule.owl.gaussian/2 + 1.5
# Value estimation
tmp <- rep(0,length(Tr.test))
for(i in 1:length(Tr.test)) tmp[i] <- potental.outcome[i, trt.rule[i]]
owl.gaussian.value <- mean(tmp)
owl.gaussian.value
# 0.7391567
sum(trt.rule == Tr.test) / length(Tr.test) # ] 0.512


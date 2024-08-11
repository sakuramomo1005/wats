library(mgcv)

#' Single-index models with multiple-links (the main function)
#'
#' \code{simml} is the wrapper function for Single-index models with multiple-links (SIMML).
#' The function estimates a linear combination (a single-index) of covariates X, and models the treatment-specific outcome y, via treatment-specific nonparametrically-defined link functions.
#'
#'
#' SIMML captures the effect of covariates via a single-index and their interaction with the treatment via nonparametric link functions.
#' Interaction effects are determined by distinct shapes of the link functions.
#' The estimated single-index is useful for comparing differential treatment efficacy.
#' The resulting \code{simml} object can be used to estimate an optimal treatment decision rule
#' for a new patient with pretreatment clinical information.
#'
#' @param y   a n-by-1 vector of treatment outcomes; y is assumed to follow an exponential family distribution; any distribution supported by \code{mgcv::gam}.
#' @param Tr  a n-by-1 vector of treatment indicators; each element represents one of the L(>1) treatment conditions; e.g., c(1,2,1,1,1...); can be a factor-valued.
#' @param X   a n-by-p matrix of pre-treatment covarates.
#' @param ortho.constr  separates the interaction effects from the main effect (without this, the interaction effect can be confounded by the main effect; the default is \code{TRUE}.
#' @param family  specifies the distribution of y; e.g., "gaussian", "binomial", "poisson"; the defult is "gaussian"; can be any family supported by \code{mgcv::gam}.
#' @param mu.hat   a n-by-1 vector for efficinecy augmentation provided by the user; the defult is \code{NULL}; the optimal choice for this vector is h(E(y|X)), where h is the canonical link function.
#' @param alpha.ini  an initial solution of \code{alpha.coef}; a p-by-1 vector; the defult is \code{NULL}.
#' @param ind.to.be.positive  for identifiability of the solution \code{alpha.coef}, we restrict the jth component of \code{alpha.coef} to be positive; by default \code{j=1}.
#' @param bs type of basis for representing the treatment-specific smooths; the defult is "ps" (p-splines); any basis supported by \code{mgcv::gam} can be used, e.g., "cr" (cubic regression splines)
#' @param k  basis dimension; the same number (k) is used for all treatment groups, however, the smooths of different treatments have different roughness parameters.
#' @param pen.order 0 indicates the ridge penalty; 1 indicates the 1st difference penalty; 2 indicates the 2nd difference penalty, used in a penalized least squares (LS) estimation of \code{alpha.coef}.
#' @param lambda  a regularziation parameter associated with the penalized LS of \code{alpha.coef}.
#' @param max.iter  an integer specifying the maximum number of iterations for \code{alpha.coef} update.
#' @param eps.iter a value specifying the convergence criterion of algorithm.
#' @param bootstrap if \code{TRUE}, compute bootstrap confidence intervals for the single-index coefficients, \code{alpha.coef}; the default is \code{FALSE}.
#' @param boot.alpha  specifies bootstrap CI percentiles; e.g., 0.05 gives 95\% CIs; 0.1 gives 90\% CIs.
#' @param nboot  when \code{bootstrap=TRUE}, a value specifying the number of bootstrap replications.
#' @param seed  when  \code{bootstrap=TRUE}, randomization seed used in bootstrap resampling.
#' @param trace.iter if \code{TRUE}, trace the estimation process and print the differences in \code{alpha.coef}.
#' @param obs d
#'
#' @return a list of information of the fitted SIMML including
#'  \item{alpha.coef}{ the estimated single-index coefficients.} \item{g.fit}{a \code{mgcv:gam} object containing information about the estimated treatment-specific link functions.} \item{alpha.ini}{the initial value used in the estimation of \code{alpha.coef}} \item{alpha.path}{solution path of \code{alpha.coef} over the iterations} \item{d.alpha}{records the change in \code{alpha.coef} over the solution path, \code{alpha.path}} \item{scale.X}{sd of pretreatment covariates X} \item{center.X}{mean of pretreatment covariates X} \item{L}{number of different treatment options} \item{p}{number of pretreatment covariates X} \item{n}{number of subjects} \item{boot.ci}{(1-boot.alpha/2) percentile bootstrap CIs (LB, UB) associated with \code{alpha.coef}}
#'
#' @author Park, Petkova, Tarpey, Ogden
#' @import mgcv plyr
#' @seealso \code{pred.simml},  \code{fit.simml}
#' @export
#'
#' @examples
#'
#' ## application of SIMML (on a simulated dataset) (see help(generate.data) for data generation).
#'
#' family <- "gaussian"   #"poisson"
#' delta = 1              # moderate main effect
#' w=2                    # if w=2 (w=1), a nonlinear (linear) contrast function
#' n=500                  # number of subjects
#' p=10                   # number of pretreatment covariates
#'
#' # generate a training dataset
#' data <- generate.data(n= n, p=p, delta = delta, w= w, family = family)
#' data$SNR  # the ratio of interactions("signal") vs. main effects("noise") in the canonical param.
#' Tr <- data$Tr
#' y  <- data$y
#' X <- data$X
#'
#' # generate a (large, 10^5) testing dataset
#' data.test <- generate.data(n=10^5, p=p, delta = delta,  w= w, family = family)
#' Tr.test <- data.test$Tr
#' y.test  <- data.test$y
#' X.test <- data.test$X
#' data.test$value.opt     # the optimal "value"
#'
#'
#' ##  estimate SIMML
#' #1) SIMML without efficiency augmenation
#' simml.obj1 <- simml(y, Tr, X, family = family)
#'
#' #2) SIMML with efficiency augmenation
#' # we can improove efficinecy by using the efficiency augmentation term, mu.hat.
#' # mu.hat is estimated by a main effect only model (y~ X).
#' glm.fit <- glm(y ~ X, family=family)  # could also use cv.glmnet to obtain a mu.hat
#' mu.hat <- as.vector(predict(glm.fit, newx =X, type="link"))
#' simml.obj2 <- simml(y, Tr, X, mu.hat = mu.hat, family = family)
#'
#'
#' ## apply the estimated SIMMLs to the testing set and obtain treatment assignment rules.
#' simml.trt.rule1 <- pred.simml(simml.obj1, newx= X.test)$trt.rule
#' # "value" estimation (estimated by IPWE)
#' simml.value1 <-  mean(y.test[simml.trt.rule1 == Tr.test])
#' simml.value1
#'
#' simml.trt.rule2 <- pred.simml(simml.obj2, newx= X.test)$trt.rule
#' # "value" estimation (estimated by IPWE)
#' simml.value2 <-  mean(y.test[simml.trt.rule2 == Tr.test])
#' simml.value2
#'
#' # compare these to the optimal "value"
#' data.test$value.opt
#'
#'
#' ## estimate the MC (modified covariates) model of Tien et al 2014
#'
#' n.t <- summary(as.factor(Tr)); pi.t <- n.t/sum(n.t)
#' mc <-  (as.numeric(Tr) + pi.t[1] -2) *cbind(1, X)   # 0.5*(-1)^as.numeric(Tr) *cbind(1, X)
#' mc.coef  <-  coef(glm(y ~ mc, family =  family))
#' mc.trt.rule <- (cbind(1, X.test) %*% mc.coef[-1] > 0) +1
#' # "value" estimation (estimated by IPWE)
#' mc.value  <-  mean(y.test[mc.trt.rule == Tr.test])
#' mc.value
#'
#'
#' ## visualization of the estimated treatment-specific link functions of SIMML
#' simml.obj1$alpha.coef      # estimated single-index coefficients
#' g.fit <- simml.obj1$g.fit   # estimated trt-specific link functions; "g.fit" is a mgcv::gam object.
#' #plot(g.fit)
#'
#'
#' ## by using the package "mgcViz", we can improve the visualization.
#' # install.packages("mgcViz")
#' # mgcViz depends on "rgl". "rgl" depends on XQuartz, which you can download from xquartz.org
#' # library(mgcViz)
#'
#' ## transform the "mgcv::gam" object to a "mgcViz" object (to improve visualization)
#' #g.fit <- getViz(g.fit)
#'
#' #plot1  <- plot( sm(g.fit,1) )  # for treatment group 1
#' #plot1 + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
#' #  l_ciLine(mul = 5, colour = "blue", linetype = 2) + l_points(shape = 19, size = 1, alpha = 0.1) +
#' #  xlab(expression(paste("z = ", alpha*minute, "x")))  +  ylab("y") +
#' #  ggtitle("Treatment group 1 (Tr =1)") +  theme_classic()
#'
#' #plot2 <- plot( sm(g.fit,2) )   # for treatment group 2
#' #plot2 + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
#' #  l_ciLine(mul = 5, colour = "blue", linetype = 2) + l_points(shape = 19, size = 1, alpha = 0.1) +
#' #  xlab(expression(paste("z = ", alpha*minute, "x"))) +ylab("y") +
#' #  ggtitle("Treatment group 2 (Tr =2)") + theme_classic()
#'
#'
#' #trans = function(x) x + g.fit$coefficients[2]
#'
#' #plotDiff(s1 = sm(g.fit, 2), s2 = sm(g.fit, 1), trans=trans) +  l_ciPoly() +
#' #  l_fitLine() + geom_hline(yintercept = 0, linetype = 2) +
#' #  xlab(expression(paste("z = ", alpha*minute, "x")) ) +
#' #  ylab("(Treatment 2 effect) - (Treatment 1 effect)") +
#' #  ggtitle("Contrast between two treatment effects") +
#' #  #geom_vline(xintercept=-0.45, linetype="dotted", color = "red", size=0.8) +
#' #  theme_classic()
#'
#'
#' # another way of visualization, using ggplot2
#' #library(ggplot2)
#' #dat  <- data.frame(y= simml.obj1$g.fit$model$y,
#' #                   x= simml.obj1$g.fit$model$single.index,
#' #                   Treatment= simml.obj1$g.fit$model$Tr)
#' #g.plot <- ggplot(dat, aes(x=x, y=y, color=Treatment, shape=Treatment, linetype=Treatment)) +
#' #   geom_point(aes(color=Treatment, shape=Treatment), size=1, fill="white") +
#' #   scale_colour_brewer(palette="Set1", direction=-1) + theme_classic() +
#' #   xlab(expression(paste(alpha*minute,"x"))) + ylab("y")
#' #g.plot + geom_smooth(method=gam, formula= y~ s(x, bs=simml.obj1$bs, k=simml.obj1$k),
#' #                     se=TRUE, fullrange=TRUE, alpha = 0.35)
#'
#'
#' #### can obtain bootstrap CIs associated with the single-index coefficients (alpha.coef).
#' glm.fit <- glm(y ~ X, family=family)  # could also use cv.glmnet.
#' mu.hat <- as.vector(predict(glm.fit, newx= X, type= "link"))    # efficiency augmentation vector
#' simml.obj <- simml(y,Tr,X, mu.hat=mu.hat, family=family, bootstrap =TRUE, nboot=15, max.iter=7)
#' # the default is to use 200 bootstrap replications.
#' simml.obj$alpha.coef
#' simml.obj$boot.ci    # displays a set of (1-boot.alpha/2) percentile bootstrap CIs (LB, UB).
#'
#' # compare the estimates to the true alpha.coef.
#' data$true.alpha
simml <- function(y, Tr, X,
                  mu.hat = NULL,
                  family = "gaussian",
                  ortho.constr = TRUE,
                  bs ="ps",
                  k = 8,
                  sp = NULL,
                  alpha.ini = NULL,
                  ind.to.be.positive = 1,
                  pen.order = 0,
                  lambda = 0,
                  max.iter = 30,
                  eps.iter = 0.01,
                  trace.iter = TRUE,
                  bootstrap = FALSE,
                  nboot= 200,
                  boot.alpha= 0.05,
                  seed= 1357,
                  obs = FALSE)
{
  
  simml.obj <- fit.simml(y=y, Tr=Tr, X=X, mu.hat=mu.hat, family=family,
                         ortho.constr =ortho.constr, bs =bs, k = k, sp=sp,
                         alpha.ini = alpha.ini, ind.to.be.positive=ind.to.be.positive,
                         pen.order = pen.order, lambda = lambda, max.iter = max.iter, trace.iter=trace.iter, obs=obs)
  
  boot.mat = boot.ci <- NULL
  if(bootstrap){
    set.seed(seed)
    indices <- 1:simml.obj$n
    boot.mat <- matrix(0, nboot, simml.obj$p)
    for(i in 1:nboot){
      boot.indices <- sample(indices, simml.obj$n, replace = TRUE)
      boot.mat[i,] <- fit.simml(y=y[boot.indices], Tr = Tr[boot.indices], X = X[boot.indices,], 
                                mu.hat = mu.hat[boot.indices],
                                family=family, ortho.constr =ortho.constr, bs =bs, k = k, sp=sp,
                                alpha.ini = alpha.ini, ind.to.be.positive=ind.to.be.positive,
                                pen.order = pen.order, lambda = lambda, max.iter = max.iter, trace.iter=trace.iter, obs=obs)$alpha.coef
    }
    ## ordinary percentile boot.ci
    #simml.obj$boot.ci <- apply(boot.mat, 2, function(x) quantile(x, probs = c(0.05, 0.95)))
    ## bias-corrected boot.ci
    prop.ind <-  matrix(0, nboot, simml.obj$p)
    for(i in 1:nboot)  prop.ind[i, ] <- (boot.mat[i,] < simml.obj$alpha.coef)
    boot.prop <- apply(prop.ind, 2, mean)
    boot.ci <- matrix(0, simml.obj$p, 2)
    for(j in 1:simml.obj$p){
      boot.ci[j, ] <- stats::quantile(boot.mat[,j], c(stats::pnorm( 2*stats::qnorm(boot.prop[j]) + stats::qnorm(boot.alpha/2)),
                                                      stats::pnorm( 2*stats::qnorm(boot.prop[j]) + stats::qnorm(1-boot.alpha/2)) ) )
    }
    boot.ci <- cbind(simml.obj$alpha.coef, boot.ci, (boot.ci[,1] > 0 | boot.ci[,2] < 0) )
    colnames(boot.ci) <- c("coef", "LB", "UB", " ***")
    rownames(boot.ci) <- colnames(X)
  }
  simml.obj$boot.mat <- boot.mat
  simml.obj$boot.ci <- boot.ci
  
  return(simml.obj)
}





#' Single-index models with multiple-links (the workhorse function)
#'
#' \code{fit.simml} is the workhorse function for Single-index models with multiple-links (SIMML).
#' The function estimates a linear combination (a single-index) of covariates X, and models the treatment-specific outcome y, via treatment-specific nonparametrically-defined link functions.
#'
#' SIMML captures the effect of covariates via a single-index and their interaction with the treatment via nonparametric link functions.
#' Interaction effects are determined by distinct shapes of the link functions.
#' The estimated single-index is useful for comparing differential treatment efficacy.
#' The resulting \code{simml} object can be used to estimate an optimal treatment decision rule
#' for a new patient with pretreatment clinical information.
#'
#' @param y   a n-by-1 vector of treatment outcomes; y is assumed to follow an exponential family distribution; any distribution supported by \code{mgcv::gam}.
#' @param Tr  a n-by-1 vector of treatment indicators; each element represents one of the L(>1) treatment conditions; e.g., c(1,2,1,1,1...); can be a factor-valued.
#' @param X   a n-by-p matrix of pre-treatment covarates.
#' @param ortho.constr  separates the interaction effects from the main effect (without this, the interaction effect can be confounded by the main effect; the default is \code{TRUE}.
#' @param family  specifies the distribution of y; e.g., "gaussian", "binomial", "poisson"; the defult is "gaussian"; can be any family supported by \code{mgcv::gam}.
#' @param mu.hat   a n-by-1 vector for efficinecy augmentation provided by the user; the defult is \code{NULL}; the optimal choice for this vector is h(E(y|X)), where h is the canonical link function.
#' @param alpha.ini  an initial solution of \code{alpha.coef}; a p-by-1 vector; the defult is \code{NULL}.
#' @param ind.to.be.positive  for identifiability of the solution \code{alpha.coef}, we restrict the jth component of \code{alpha.coef} to be positive; by default \code{j=1}.
#' @param bs type of basis for representing the treatment-specific smooths; the defult is "ps" (p-splines); any basis supported by \code{mgcv::gam} can be used, e.g., "cr" (cubic regression splines)
#' @param k  basis dimension; the same number (k) is used for all treatment groups, however, the smooths of different treatments have different roughness parameters.
#' @param pen.order 0 indicates the ridge penalty; 1 indicates the 1st difference penalty; 2 indicates the 2nd difference penalty, used in a penalized least squares (LS) estimation of \code{alpha.coef}.
#' @param lambda  a regularziation parameter associated with the penalized LS of \code{alpha.coef}.
#' @param max.iter  an integer specifying the maximum number of iterations for \code{alpha.coef} update.
#' @param eps.iter a value specifying the convergence criterion of algorithm.
#' @param trace.iter if \code{TRUE}, trace the estimation process and print the differences in \code{alpha.coef}.
#' @param obs d
#'
#' @return a list of information of the fitted SIMML including
#'  \item{alpha.coef}{ the estimated single-index coefficients.} \item{g.fit}{a \code{mgcv:gam} object containing information about the estimated treatment-specific link functions.} \item{alpha.ini}{the initial value used in the estimation of \code{alpha.coef}} \item{alpha.path}{solution path of \code{alpha.coef} over the iterations} \item{d.alpha}{records the magnitude of change in \code{alpha.coef} over the solution path, \code{alpha.path}} \item{scale.X}{sd of pretreatment covariates X} \item{center.X}{mean of pretreatment covariates X} \item{L}{number of different treatment options} \item{p}{number of pretreatment covariates X} \item{n}{number of subjects}
#'
#' @author Park, Petkova, Tarpey, Ogden
#' @import mgcv plyr
#' @seealso \code{pred.simml},  \code{fit.simml}
#' @export
#'
fit.simml  <- function(y, Tr, X,
                       mu.hat = NULL,
                       family = "gaussian",
                       ortho.constr = TRUE,
                       bs ="ps",
                       k = 8,
                       sp = NULL,
                       alpha.ini = NULL,
                       ind.to.be.positive =1,
                       pen.order = 0,
                       lambda = 0,
                       max.iter = 30,
                       eps.iter = 0.01,
                       trace.iter= TRUE,
                       obs=FALSE)
{
  
  Tr <- as.factor(Tr)
  L <- length(levels(Tr))
  n <- length(y)
  p <- ncol(X)
  n.t <- summary(Tr)
  pi.t <- n.t/sum(n.t)
  
  ## Center and scale X
  Xc <- scale(X, center = TRUE, scale = TRUE)
  scale.X <-  attr(Xc, "scaled:scale")
  center.X <- attr(Xc, "scaled:center")
  
  ## If not provided the user, the efficiency augmentation vector is set to be a zero vector.
  if(is.null(mu.hat)) mu.hat <- rep(0, n)
  
  ## Specify a penalty matrix associated with the penalized LS for estimating alpha.coef.
  D <- diag(p);  if(pen.order != 0)  for(j in 1:pen.order) D <- diff(D);
  Pen <- sqrt(lambda)*D
  
  ## Initialization
  if(is.null(alpha.ini)){
    if(L==2){
      # initialize alpha.coef by the modified covariate method of Tien et al 2014.
      mc <-  (as.numeric(Tr) + pi.t[1] -2) *cbind(1, Xc)  # 0.5*(-1)^as.numeric(Tr) *cbind(1, Xc)
      alpha.ini <-  stats::coef(stats::glm(y ~ mc, family =  family))[-c(1,2)]
    }else{
      # initialize alpha.coef by the linear GEM method of Petkova et al 2017 (using the "numerator" criterion).
      dat <- data.frame(Tr=Tr, y=y, Xc = Xc)
      dat.list <- dlply(dat, .(dat$Tr), function(x) {as.matrix(x[,-1])})
      B <- matrix(0, p, L)
      for(t in 1:L) B[,t] <- stats::coef(stats::glm(dat.list[[t]][,1]  ~ dat.list[[t]][,-1], family =  family))[-1];
      Bc <- B - apply(B, 1, function(x) stats::weighted.mean(x, w=pi.t))
      alpha.ini <- eigen(Bc %*% t(Bc))$vectors[,1]
    }
    names(alpha.ini) <- colnames(X)
  }
  alpha.coef <- alpha.ini/sqrt(sum(alpha.ini^2))  # enforce unit L2 norm
  if(alpha.coef[ind.to.be.positive] < 0) alpha.coef <- -1*alpha.coef       # for the (sign) identifiability
  # initialize the single index
  single.index <- as.vector(Xc %*% alpha.coef)
  
  
  
  ## Given a single.index, estimate treatment-specific smooths (subject to the interaction effect identifiability constraint).
  g.fit <- gam(y ~ Tr + s(single.index, by = Tr,  bs=bs, k=k), family =  family, gamma=1.4, sp = sp)
  smooths.coef <- stats::coef(g.fit)[-c(1:L)]
  B <- matrix(smooths.coef, ncol=L)
  # the following enforces the identifiability constraint for the interaction component
  if(ortho.constr){
    if(obs){
      g0.fit <- gam(y ~ s(single.index, bs=bs, k=k), family =  family, gamma=1.4,  sp = sp)
      B <- B - stats::coef(g0.fit)[-1]
    }else{
      B <- B - apply(B, 1, function(x) stats::weighted.mean(x, w=pi.t))
    }
  }
  g.fit$coefficients[-c(1:L)] <- as.vector(B)
  
  ## Compute the adjusted responses (adjusted, for the nonlinearity associated with the GLM canonical link)
  # obtain the 1st derivative of the inverse canonical link, w.r.t. the "linear" predictor
  h.prime <- g.fit$family$mu.eta(predict.gam(g.fit, type="link"))
  adjusted.responses <- (y - predict.gam(g.fit, type = "response")) /h.prime + predict.gam(g.fit, type="link")
  # take the 1st deriavative of the treatment-specific smooths, w.r.t. the single.index.
  g.der <- der.link(g.fit)
  
  alpha.path <- alpha.coef
  d.alpha <- NULL
  
  ## Start iteration
  for (it in 2:max.iter)
  {
    ## Update alpha.coef and intercept through lsfit
    # adjusted responses, adjusted for the nonlinearity associated with the treatment-specific smooths
    y.star    <- adjusted.responses - predict.gam(g.fit, type="link")  + g.der*single.index - mu.hat
    # adjusetd covariates, adjusted for the nonlinearity of the the treatment smooths
    X.tilda   <- diag(g.der) %*% Xc
    nix       <- rep(0, nrow(D))
    X.p       <- rbind(X.tilda, Pen)
    y.p       <- c(y.star, nix)
    # perform a (penalized) WLS
    alpha.fit   <- stats::lsfit(X.p, y.p, wt =c(g.fit$weights, (nix + 1)))
    # for the identifiability
    alpha.new <- alpha.fit$coef[-1]/sqrt(sum(alpha.fit$coef[-1]^2))
    if(alpha.new[ind.to.be.positive] < 0) alpha.new <- -1*alpha.new
    alpha.path <- rbind(alpha.path, alpha.new)
    
    ## Check the convergence of alpha
    d.alpha   <- c(d.alpha, sum((alpha.new-alpha.coef)^2)/ sum(alpha.new^2))
    if(trace.iter){
      cat("iter:", it, " "); cat(" difference in alpha: ", d.alpha[(it-1)], "\n")
    }
    if (d.alpha[(it-1)] < eps.iter)
      break
    alpha.coef <- alpha.new
    single.index <- as.vector(Xc %*% alpha.coef)
    
    ## Given a single.index, estimate treatment-specific smooths (subject to the interaction effect identifiability constraint).
    g.fit <- gam(y ~ Tr + s(single.index, by = Tr, bs=bs, k=k), family =  family, gamma=1.4,  sp = sp)
    smooths.coef <- stats::coef(g.fit)[-c(1:L)]
    B <- matrix(smooths.coef, ncol=L)
    # the following enforces the identifiability constraint for the interaction component
    if(ortho.constr){
      if(obs){
        g0.fit <- gam(y ~ s(single.index, bs=bs, k=k), family =  family, gamma=1.4,  sp = sp)
        B <- B - stats::coef(g0.fit)[-1]
      }else{
        B <- B - apply(B, 1, function(x) stats::weighted.mean(x, w=pi.t))
      }
    }
    g.fit$coefficients[-c(1:L)] <- as.vector(B)
    
    ## compute the adjusted responses (adjusted, for the nonlinearity associated with the GLM canonical link)
    # obtain the 1st derivative of the inverse canonical link, w.r.t. the "linear" predictor
    h.prime <- g.fit$family$mu.eta(predict.gam(g.fit, type="link"))
    adjusted.responses <- (y - predict.gam(g.fit, type = "response")) /h.prime + predict.gam(g.fit, type="link")
    # take the 1st deriavative of the treatment-specific smooths, w.r.t. the single.index.
    g.der <- der.link(g.fit)
  }
  
  results <- list(alpha.coef = round(alpha.coef,3),
                  alpha.ini = alpha.ini, d.alpha=d.alpha, alpha.path=alpha.path,
                  g.fit= g.fit,
                  scale.X=scale.X, center.X = center.X,
                  L=L, p=p, n=n, bs=bs, k=k)
  class(results) <- c("simml", "list")
  
  return(results)
}


#' A subfunction used in estimation
#'
#' This function computes the 1st derivative of the treatment-specific "smooth" w.r.t. the single index, using finite difference.
#'
#' @param g.fit  a \code{mgcv::gam} object
#' @param arg.number the argument of \code{g.fit} that we take derivative w.r.t.; the default is \code{arg.number=2} (i.e., take deriviative w.r.t. the single-index.)
#' @param eps a small finite difference used in numerical differentiation.
#' @seealso \code{fit.simml}, \code{simml}
#'
der.link <- function(g.fit, arg.number =2, eps=10^(-6))
{
  m.terms <- attr(stats::terms(g.fit), "term.labels")
  newD <- stats::model.frame(g.fit)[, m.terms, drop = FALSE]
  newDF <- data.frame(newD)  # needs to be a data frame for predict
  X0 <- predict.gam(g.fit, newDF, type = "lpmatrix")
  newDF[,arg.number] <- newDF[,arg.number] + eps
  X1 <- predict.gam(g.fit, newDF, type = "lpmatrix")
  Xp <- (X1 - X0) / eps
  Xi <- Xp * 0
  want <- grep(m.terms[arg.number], colnames(X1))
  Xi[, want] <- Xp[, want]
  g.der  <- as.vector(Xi %*% stats::coef(g.fit))  # the first derivative of eta
  return(g.der)
}


#' SIMML prediction function
#'
#' This function makes predictions from an estimated SIMML, given a (new) set of pretreatment covariates.
#' The function returns a set of predicted outcomes for each treatment condition and a set of recommended treatment assignments (assuming a larger value of the outcome is better).
#'
#' @param simml.obj  a \code{simml} object
#' @param newx a (n-by-p) matrix of new values for the pretreatment covariates X at which predictions are to be made.
#' @param type the type of prediction required; the default "response" is on the scale of the response variable;  the alternative "link" is on the scale of the linear predictors.
#' @param maximize the default is \code{TRUE}, assuming a larger value of the outcome is better; if \code{FALSE}, a smaller value is assumed to be prefered.
#'
#' @return
#' \item{pred.new}{a (n-by-L) matrix of predicted values; each column represents a treatment option.}
#' \item{trt.rule}{a (n-by-1) vector of suggested treatment assignments}
#'
#'
#' @author Park, Petkova, Tarpey, Ogden
#' @seealso \code{simml},\code{fit.simml}
#' @export
#'
pred.simml  <-  function(simml.obj, newx, type = "response", maximize=TRUE)
{
  #if(!inherits(simml.obj, "simml"))   # checks input
  #  stop("obj must be of class `simml'")
  
  if(ncol(newx) != simml.obj$p)
    stop("newx needs to be of p columns ")
  
  alpha.coef  <- simml.obj$alpha.coef
  g.fit <- simml.obj$g.fit
  newx.scaled <- scale(newx, center = simml.obj$center.X, scale = simml.obj$scale.X)
  single.index  <- newx.scaled %*% alpha.coef
  
  L <- simml.obj$L  # the number of treatment options (levels)
  
  # compute treatment-specific predicted outcomes
  pred.new <- matrix(0, nrow(newx), L)
  for(t in 1:L){
    newd <- data.frame(Tr= t, single.index=single.index)
    pred.new[ ,t] <- predict.gam(g.fit, newd, type =type)
    rm(newd)
  }
  
  # compute treatment assignment
  if(maximize){
    trt.rule <- apply(pred.new, 1, which.max)  # assume a larger value of y is prefered.
  }else{
    trt.rule <- apply(pred.new, 1, which.min)  # assume a smaller value of y is prefered.
  }
  
  if(L==2)  colnames(pred.new) <- c("Tr1", "Tr2")
  
  return(list(trt.rule = trt.rule, pred.new = pred.new))
}



#' A dataset generation function
#'
#' \code{generate.data} generates an example dataset from a mean model that has a "main" effect component and a treatment-by-covariates interaction effect component (and a random component for noise).
#'
#' @param n  sample size.
#' @param p  dimension of covariates.
#' @param family specifies the distribution of the outcome y;  "gaussian", "binomial", "poisson"; the defult is "gaussian"
#' @param sigma  standard deviation of the random noise term (for gaussian response).
#' @param sigmaX  standard deviation of the covariates.
#' @param correlationX  correlation among the covariates.
#' @param pi.1  probability of being assigned to the treatment 1
#' @param w  controls the nonliarity of the treatment-specific link functions that define the interaction effect component.
#' \describe{
#' \item{\code{w=1}}{linear}
#' \item{\code{w=2}}{nonlinear}
#' }
#' @param delta  controls the intensity of the main effect; can take any intermediate value, e.g., \code{delta= 1.4}.
#' \describe{
#' \item{\code{delta=1}}{moderate main effect}
#' \item{\code{delta=2}}{big main effect}
#' }
#' @param xi     if nonzero, the interaction-effect part of the model deviates from a GEM (single-index) model.
#' @param true.alpha  a p-by-1 vector of the true single-index coefficients (associated with the interaction effect component); if \code{NULL}, \code{true.alpha} is set to be \code{(1, 0.5, 0.25, 0.125, 0,...0)}' (only the first 4 elements are nonzero).
#' @param true.eta   a p-by-1 vector of the true main effect coefficients; if \code{NULL}, \code{true.eta} is set to be \code{(0,..., 0.125, 0.25, 0.25, 1)}' (only the last 4 elements are nonzero).
#'
#'
#' @return
#' \item{y}{a n-by-1 vector of treatment outcomes.}
#' \item{Tr}{a n-by-1 vector of treatment indicators.}
#' \item{X}{a n-by-p matrix of pretreatment covariates.}
#' \item{SNR}{the "signal" (interaction effect) to "nuisance" (main effect) variance ratio (SNR) in the canonical parameter function.}
#' \item{true.alpha}{the true single-index coefficient vector.}
#' \item{true.eta}{the true main effect coefficient vector.}
#' \item{optTr}{a n-by-1 vector of treatments, indicating the optimal treatment selections.}
#' \item{value.opt}{the "value" implied by the optimal treatment decision rule, \code{optTr}.}
#' @export
#'
generate.data <- function(n = 200, # number of observations
                          p = 10,  # number of covariates
                          family = "gaussian",  # the distribution of the outcome y
                          correlationX= 0, # correlation among pretreatment covariates X
                          sigmaX = 1, # pretreatment covariate sd
                          sigma = 0.5, # error sd (for gaussian response)
                          w = 1, # shape of the interaction effect curves (1 bell-shaped; 2 cosine; 3: sine)
                          delta = 1,  # magnitude of the main effect
                          xi = 0,   # controls the deviation from a single-index model 
                          pi.1 = 0.5,  # probability of being assigned to the treatment 1
                          true.alpha= NULL,
                          true.eta= NULL,
                          obs= FALSE)
{
  
  if(is.null(true.alpha)){
    true.alpha <- c(c(1, 0.5, 0.25, 0.125), rep(0, p-4))  # only the first 4 components are nonzero.
    true.alpha <- true.alpha/sqrt(sum(true.alpha^2))    # the true single index coefficients
  }
  if(length(true.alpha)!= p)   stop("true.alpha must be of length p")
  true.alpha2 <- c(1,1,1,1,1,1, rep(0, p-6))/sqrt(6)
  
  if(is.null(true.eta)){
    if(correlationX!=0){
      eta.hold <- runif(6, -1,1)
    }else{
      eta.hold <- c(-1, 1, -1, 1, -1, 1) 
    }
    eta.hold  <- eta.hold /sqrt(sum(eta.hold^2))
    true.eta <- c(eta.hold, rep(0, p-6))
  }
  if(length(true.eta)!= p)   stop("true.eta must be of length p")
  
  
  # Pre-treatment covariates (the inputs)
  if(correlationX!=0){
    Psix <- sigmaX*(diag(1 - correlationX, nrow = p, ncol = p) + matrix(correlationX, nrow = p, ncol = p) )   # X covariance matrix.
    ePsix <- eigen(Psix)
    tmp <- sapply(1:p, function(x) stats::rnorm(n)) %*% diag(sqrt(ePsix$values)) %*% t(ePsix$vectors)
    tmp2 <- apply(tmp, 2, pnorm)
    X <- tmp2 *pi - pi/2
  }else{
    X <- sapply(1:p, function(x) runif(n, -pi/2, pi/2))
  }
  
  #sd(drop(X %*% true.eta))
  
  # Treatment variable (the controling variable)
  if(obs){
    mu.Tr  <- -1 + 0.2*X[,1]^2 + 0.2*X[,2]^2 + 0.2*X[,3]^2 + 0.2*X[,4]^2
    p.Tr <- 1/(1+ exp(-mu.Tr) )
    Tr <-  rbinom(n, 1, p.Tr) +1
  }else{
    Tr <- drop(stats::rbinom(n, 1, pi.1) + 1)  # generate treatment variables
  }
  
  # delta is the intensity parameter (delta = 1: moderate main effect; delta=2: big main effect)
  # xi controls the deviation from a single-index Tr-by-X interaction model 
  if(w==1){  # simulation set A
    contrast <- exp(-(drop(X %*% true.alpha)-0.5)^2 - (drop(X %*% true.alpha2)-0.5)^2 *xi) - 0.5
    main.effect <- 0.8*delta*cos(drop(X %*% true.eta))
  }
  if(w==2){  # simulation set B
    contrast <- cos( drop(X %*% true.alpha) + (Tr==2)*drop(X %*% true.alpha2)*xi - pi/8) - 0.5
    main.effect <- delta*cos(drop(X %*% true.eta))
  }
  if(w==3){  # an additional simulation setting 
    contrast <- sin( drop(X %*% true.alpha) + (Tr==2)*drop(X %*% true.alpha2)*xi - pi/8) - 0.5
    main.effect <- delta*cos(drop(X %*% true.eta))
  }
  interaction.effect <- 2*(as.numeric(Tr) + pi.1 -2) * contrast   #  (-1)^Tr* contrast 
  
  
  # the hypothetical outcomes, attributable to the interaction
  if(family == "gaussian"){
    mu.inter1 <-  2*(pi.1 - 1) *contrast   #-contrast         # if Tr =1
    mu.inter2 <-  2*pi.1*contrast   # contrast        # if Tr =2
  }
  if(family == "binomial"){
    mu.inter1 <-   1/(1+ exp(-(2*(pi.1 - 1) *contrast )))    # 1/(1+ exp(-(-contrast)))  # if Tr =1
    mu.inter2 <-   1/(1+ exp(- 2*pi.1*contrast ))    # 1/(1+ exp(- contrast))   # if Tr =2
  }
  if(family == "poisson"){
    mu.inter1 <-  exp(2*(pi.1 - 1)*contrast)  # exp(-contrast)   # if Tr =1
    mu.inter2 <-  exp(2*pi.1*contrast)  # exp( contrast)   # if Tr =2
  }
  
  # the canonical parameter
  theta <- main.effect + interaction.effect
  
  # the "signal" to "noise" ratio
  noise <- sigma * stats::rnorm(n)
  SNR <- stats::var(interaction.effect)/stats::var(main.effect)
  
  if(family == "gaussian"){
    mu <- theta
    y <-  mu  + noise
  }
  if(family == "binomial"){
    mu <- 1/(1+ exp(-theta))
    y <- stats::rbinom(n, size=1, prob= mu)
  }
  if(family == "poisson"){
    mu <-  exp(theta)
    y <- stats::rpois(n, lambda= mu)
  }
  
  optTr <- as.numeric(mu.inter2 > mu.inter1) + 1  # this takes 1 or 2
  
  value.opt.IPWE <- mean(mu[Tr == optTr])
  value.opt.IPWE
  
  potental.outcome <- main.effect  + cbind(mu.inter1, mu.inter2)
  value.opt <- mean(apply(potental.outcome, 1, max))
  
  
  return(list(y=y, Tr=Tr, X=X, SNR=SNR, true.alpha=true.alpha, true.eta=true.eta, delta=delta, w=w, xi=xi, 
              optTr=optTr, value.opt=value.opt, potental.outcome=potental.outcome, value.opt.IPWE=value.opt.IPWE))
}






######################################################################
## END OF THE FILE
######################################################################

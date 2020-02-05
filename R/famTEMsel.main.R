#--------------------
# Package: famTEMsel
#--------------------

#library(samTEMsel)
#library(mgcv)
#library(splines)
#library(SAM)


#' Functional Additive Models for Treatment Effect-Modifier Selection (cross-validation function)
#'
#' Does k-fold cross-validation for \code{\link{famTEMsel}}, selects an optimal regularization parameter index, \code{lambda.opt.index},
#' and returns the estimated constrained functional additive model given the optimal regularization parameter index \code{lambda.opt.index}.
#'
#'
#'
#' @param y   a n-by-1 vector of responses
#' @param A   a n-by-1 vector of treatment variable; each element represents one of the L(>1) treatment conditions; e.g., c(1,2,1,1,3...); can be a factor-valued
#' @param X   a length-p list of functional-valued covariates, with its jth element corresponding to a n-by-n.eval[j] matrix of the observed jth functional covariates; n.eval[j] represents the number of evaluation points of the jth functional covariates
#' @param Z   a n-by-q matrix of scalar-valued covaraites
#' @param mu.hat   a n-by-1 vector of the fitted (X,Z)-main effect term of the model provided by the user; defult is \code{NULL}, in which case \code{mu.hat} is taken to be a vector of zeros; the optimal choice for this vector is E(y|X,Z)
#' @param n.folds  number of folds for cross-validation; the default is 5.
#' @param d     number of basis spline functions to be used for each component function; the default value is 3; d=1 corresponds to the linear model
#' @param k     dimension of the basis for representing each single-index coefficient function; see \code{mgcv::gam} for detail; the default value is 6.
#' @param bs    type of basis for representing the single-index coefficient functions; the defult is "ps" (p-splines); any basis supported by mgcv::gam can be used, e.g., "cr" (cubic regression splines).
#' @param sp    smoothing parameter associated with the single-index coefficient function; the default is \code{NULL}, in which case the smoothing parameter is estimated based on generalized cross-validation.
#' @param lambda.opt.index  a user-supplied optimal regularization parameter index to be used; the default is \code{NULL}, in which case n.folds cross-validation is performed to select an optimal index.
#' @param lambda   a user-supplied regularization parameter sequence; typical usage is to have the program compute its own lambda sequence based on nlambda and lambda.min.ratio.
#' @param nlambda  total number of lambda values; the default value is 30.
#' @param lambda.min.ratio  the smallest value for lambda, as a fraction of lambda.max, the (data derived) entry value (i.e. the smallest value for which all coefficients are zero); the default is 1e-2.
#' @param lambda.index.grid  a set of indices of \code{lambda}, in which the search for an optimal regularization parameter is to be conducted.
#' @param cv1sd  if \code{TRUE}, an optimal regularization parameter is chosen based on: the mean cross-validated error + 1 SD of the mean cross-validated error, which typically results in an increase in regularization; the defualt is \code{FALSE}.
#' @param thol  stopping precision for the coordinate-descent algorithm; the default value is 1e-5.
#' @param max.ite  number of maximum iterations for the coordinate-descent procedure used in estimating the component functions; the default value is 1e+5.
#' @param regfunc  type of the regularizer for variable selection; the default is "L1"; can also be "MCP" or "SCAD".
#' @param max.iter  number of maximum iterations for the iterative procedure used in estimating the single-index coefficient functions; the default value is 1e+1.
#' @param eps.iter  	a value specifying the convergence criterion for the iterative procedure used in estimating the single-index coefficient functions; the defult is 1e-2.
#' @param eps.num.deriv  a small value that defines a finite difference used in computing the numerical (1st) derivative of the estimated component function; the default is 1e-4.
#' @param plots  if \code{TRUE}, produce a cross-validation plot of the estimated mean squared error versus the regulariation parameter index.
#' @param trace.iter if \code{TRUE}, trace the estimation process by printing the difference in the estimated single-index basis coefficients (as compared to the previous iteration), and the functional norms of the estimated component functions, for each iteration.
#'
#'
#' @return a list of information of the fitted constrained functional additive model including
#'  \item{famTEMsel.obj}{an object of class \code{famTEMsel}, which contains the sequence of the set of fitted component functions \code{samTEMsel.obj} implied by the sequence of the regularization parameters \code{lambda} and the corresponding set of fitted single-index coefficient functions \code{si.fit}; see \code{\link{famTEMsel}} for detail.}
#'  \item{lambda.opt.index}{an index number, indicating the index of the estimated optimal regularization parameter in \code{lambda}.}
#'  \item{nonzero.index}{a set of numbers, indicating the indices of estimated nonzero component functions, evalated at the regularization parameter index \code{lambda.opt.index}.}
#'  \item{nonzero.X.index}{a set of numbers, indicating the indices of estimated nonzero component functions associated with the p functional covariates, evalated at the regularization parameter index \code{lambda.opt.index}.}
#'  \item{func_norm.opt}{a p-by-1 vector, indicating the norms of the estimated component functions evaluatd at the regularization parameter index \code{lambda.opt.index}, with each element corresponding to the norm of each estimated component function.}
#'  \item{cv.storage}{a n.folds-by-nlambda matrix of the estimated mean squared errors, with each column corresponding to each of the regularization parameters in \code{lambda} and each row corresponding to each of the n.folds folds.}
#'  \item{mean.cv}{a nlambda-by-1 vector of the estimated mean squared errors, with each element corresponding to each of the regularization parameters in \code{lambda}.}
#'  \item{sd.cv}{a nlambda-by-1 vector of the standard deviation of the estimated mean squared errors, with each element corresponding to each of the regularization parameters in \code{lambda}.}
#'
#'
#' @author Park, Petkova, Tarpey, Ogden
#' @importFrom splines ns
#' @import SAM mgcv graphics stats
#' @seealso \code{famTEMsel}, \code{predict_famTEMsel}, \code{plot_famTEMsel}
#' @export
#' @examples
#'
#' p = q = 10 # p and q are the numbers of functional and scalar pretreatment covariates, respectively.
#' n.train = 300  # training set sample size
#' n.test = 1000  # testing set sample size
#' n = n.train + n.test
#'
#' # generate p pretreatment functional covariates X by first seting up functional basis:
#' n.eval = 50; s = seq(0, 1, length.out = n.eval)  # a grid of support points
#' b1 = sqrt(2)*sin(2*pi*s)
#' b2 = sqrt(2)*cos(2*pi*s)
#' b3 = sqrt(2)*sin(4*pi*s)
#' b4 = sqrt(2)*cos(4*pi*s)
#' B  = cbind(b1, b2, b3, b4)   # a (n.eval-by-4) basis matrix
#' # randomly generate basis coefficients, and then add measurement noise
#' X = vector("list", length= p);  for(j in 1:p){
#'   X[[j]] = matrix(rnorm(n*4, 0, 1), n, 4) %*% t(B) +
#'     matrix(rnorm(n*n.eval, 0, 0.25), n, n.eval)    # measurement noise
#' }
#' Z  = matrix(rnorm(n*q, 0, 1), n, q)   # q scalar covariates
#' A  = rbinom(n, 1, 0.5) + 1  # treatment variable taking a value in {1,2} with equal prob.
#'
#' # X main effect on y; depends on the first 5 covariates
#' # the effect is generated randomly; randomly generated basis coefficients, scaled to unit L2 norm.
#' tmp = apply(matrix(rnorm(4*5), 4,5), 2, function(s) s/sqrt(sum(s^2) ))
#' main.effect = rep(0, n); for(j in 1:5){
#'   main.effect = main.effect + cos(X[[j]]%*% B %*%tmp[,j]/n.eval) # nonlinear effect (cosine)
#' }; rm(tmp)
#' # Z main effect on y; also depends on first 5 covariates
#' for(k in 1:5){
#'   main.effect = main.effect + cos(Z[,k])
#' }
#'
#' # define (interaction effect) coefficient functions associted with X[[1]] and X[[2]]
#' beta1 =  B %*% c(0.5,0.5,0.5,0.5)
#' beta2 =  B %*% c(0.5,-0.5,0.5,-0.5)
#' # A-by-X ineraction effect on y; depends only on X[[1]] and X[[2]].
#' interaction.effect = (A-1.5)*( 2*sin(X[[1]]%*%beta1/n.eval) + 2*sin(X[[2]]%*%beta2/n.eval))
#' # A-by-Z ineraction effect on y; depends only on Z[,1] and Z[,2].
#' interaction.effect =  interaction.effect +  (A-1.5)*(Z[,1] + 2*sin(Z[,2]))
#'
#' # generate outcome y
#' noise = rnorm(n, 0, 0.5)
#' y = main.effect  + interaction.effect + noise
#'
#' var.main <- var(main.effect)
#' var.interaction <- var(interaction.effect)
#' var.noise <- var(noise)
#' SNR <- var.interaction/ (var.main + var.noise)
#' SNR  # "signal-to-noise" ratio
#'
#'
#' # train/test set splitting
#' train.index = 1:n.train
#' y.train = y[train.index]
#' X.train= X.test = vector("list", p);  for(j in 1:p){
#' X.train[[j]] = X[[j]][train.index,]
#' X.test[[j]]  = X[[j]][-train.index,]
#' }
#' A.train = A[train.index]
#' A.test  = A[-train.index]
#' y.train = y[train.index]
#' y.test  = y[-train.index]
#' Z.train = Z[train.index,]
#' Z.test  = Z[-train.index,]
#'
#'
#' # obtain an optimal regularization parameter and the corresponding model by running cv.famTEMsel().
#' \donttest{cv.obj = cv.famTEMsel(y.train, A.train, X.train, Z.train)
#' lambda.opt.index = cv.obj$lambda.opt.index   # optimal regularization parameter index
#' cv.obj$func_norm.opt  # L2 norm of the component functions, associated with lambda.opt.index.
#' famTEMsel.obj = cv.obj$famTEMsel.obj  # extract the fitted model associted with lambda.opt.index.
#' # see also, famTEMsel() for the detail of famTEMsel.obj.
#'
#' famTEMsel.obj$nonzero.index  # set of indices for the component functions estimated as nonzero
#' # plot the component functions estimated as nonzero
#' plot_famTEMsel(famTEMsel.obj, which.index = famTEMsel.obj$nonzero.index)
#'
#' # make ITRs for subjects with pretreatment characteristics, X.test and Z.test
#' trt.rule = make_ITR_famTEMsel(famTEMsel.obj, newX = X.test, newZ = Z.test)$trt.rule
#' head(trt.rule)
#'
#' # an (IPWE) estimate of the "value" of this particualr treatment rule, trt.rule:
#' mean(y.test[A.test==trt.rule])
#'
#' # compare the above value to the following estimated "values" of "naive" treatment rules:
#' mean(y.test[A.test==1])   # a rule that assigns everyone to A=1
#' mean(y.test[A.test==2])   # a rule that assigns everyone to A=2
#'}
#'
# cross-validation function for the sparsity tuning parameter selection for famTEMsel()
# the optimal lambda index is searched through the set of indices, lambda.index.grid.
cv.famTEMsel <- function(y, A, X, Z= NULL, mu.hat = NULL,
                         n.folds = 5, d = 3, k = 6, bs="ps", sp = NULL,
                         lambda.opt.index = NULL, lambda = NULL,
                         nlambda = 30, lambda.min.ratio = 1e-02,
                         lambda.index.grid = 1:floor(nlambda/3), cv1sd = FALSE,
                         thol = 1e-05, max.ite = 1e+05, regfunc = "L1",
                         max.iter=1e+1, eps.iter=1e-02, eps.num.deriv = 1e-04,
                         trace.iter = TRUE, plots = TRUE)
{

  if(is.null(lambda.opt.index)){
    folds = sample(cut(seq(1, length(y)), breaks= n.folds, labels=FALSE))
    cv.storage  = matrix(NA, n.folds, length(lambda.index.grid))
    if(is.null(mu.hat))  mu.hat = rep(0, length(y))
    for(kk in 1:n.folds){
      testIndexes = which(folds==kk, arr.ind=TRUE)
      X.train = X.test = vector("list", length(X)); for(j in 1:length(X)){
        X.train[[j]] = X[[j]][-testIndexes,]
        X.test[[j]]  = X[[j]][ testIndexes,]
      }
      Z.train = Z[-testIndexes,]
      Z.test  = Z[ testIndexes,]
      y.train = y[-testIndexes]
      y.test  = y[ testIndexes]
      A.test  = A[testIndexes]
      A.train = A[-testIndexes]
      mu.hat.train = mu.hat[-testIndexes]
      mu.hat.test  = mu.hat[ testIndexes]
      for(ll in seq_along(lambda.index.grid)){
        train.fit = famTEMsel(y.train, A.train, X.train, Z.train, mu.hat.train, d=d, k=k, bs=bs, sp=sp,
                              lambda.index = lambda.index.grid[ll],
                              lambda = lambda, nlambda=nlambda, lambda.min.ratio=lambda.min.ratio,
                              thol=thol, max.ite=max.ite, regfunc=regfunc,
                              max.iter = max.iter, eps.iter=eps.iter, eps.num.deriv=eps.num.deriv,
                              trace.iter=trace.iter)
        test.pred <- predict_famTEMsel(train.fit, newX= X.test, newZ= Z.test, newA= A.test,
                                       lambda.index= lambda.index.grid[ll])$predicted
        cv.storage[kk, ll] <- mean((y.test - mu.hat.test - test.pred)^2)
      }
    }
    mean.cv  = apply(cv.storage, 2, mean)
    sd.cv    = apply(cv.storage, 2, sd)/sqrt(n.folds)
    if(plots){
      plot(mean.cv, xlab=expression(paste("Regularization parameter ", lambda, " index")), ylab = "Mean squared error",cex.lab=1,lwd=2)
      lines(mean.cv+sd.cv, type = "l",  col ="red", lwd= 0.2)
      lines(mean.cv-sd.cv, type = "l",  col ="red", lwd= 0.2)
    }
    if(cv1sd){
      cv1sd = mean.cv[which.min(mean.cv)] + sd.cv[which.min(mean.cv)]  # move the index number in the direction of increasing regularization
      opt.ll = which(cv1sd > mean.cv)[1]
    }else{
      opt.ll = which.min(mean.cv)
    }
    lambda.opt.index = lambda.index.grid[opt.ll]
  }else{
    cv.storage= mean.cv = NULL
  }

  famTEMsel.obj = famTEMsel(y, A, X, Z, mu.hat, d=d, k = k, bs=bs, sp=sp,
                            lambda.index = lambda.opt.index, lambda = lambda,
                            nlambda=nlambda, lambda.min.ratio=lambda.min.ratio,
                            thol=thol, max.ite=max.ite, regfunc=regfunc,
                            max.iter = max.iter, eps.iter = eps.iter, eps.num.deriv=eps.num.deriv,
                            trace.iter=trace.iter)

  res = list(famTEMsel.obj = famTEMsel.obj, lambda.opt.index= lambda.opt.index, lambda = famTEMsel.obj$lambda,
             nonzero.index = famTEMsel.obj$nonzero.index, nonzero.X.index = famTEMsel.obj$nonzero.X.index,
             func_norm.opt = famTEMsel.obj$func_norm.final, mean.cv=mean.cv, sd.cv=sd.cv, cv.storage=cv.storage)
  return(res)
}





#' Functional Additive Models for Treatment Effect-Modifier Selection (main function)
#'
#' The function \code{famTEMsel} implements estimation of a constrained functional additve model.
#'
#' A constrained functional model represents the joint effects of treatment, pretreatment p functional covariates and q scalar covariates on an outcome variable via a sum of treatment-specific additive flexible component functions defined over the (p + q) covariates, subject to the constraint that the expected value of the outcome given the covariates equals zero, while leaving the main effects of the covariates unspecified.
#' The p pretreatment functional covariates appear in the model as 1-dimensional projections, via inner products with corresponding single-index coefficient functions.
#' Under this model, the treatment-by-covariates interaction effects are determined by distinct shapes (across treatment levels) of the treatment-specific flexible component functions.
#' Optimized under a penalized least square criterion with a L1 (or SCAD/MCP) penalty, the constrained functional additive model can effectively identify/select treatment effect-modifiers (from the p functional and q scalar covariates) that exhibit possibly nonlinear interactions with the treatment variable; this is achieved by producing a sparse set of estimated component functions of the model.
#' The estimated nonzero component functions and single-index coefficient functions (available from the returned \code{famTEMsel} object) can be used to make individualized treatment recommendations (ITRs) for future subjects; see also \code{make_ITR_famTEMsel} for such ITRs.
#'
#'
#' The regularization path for the component functions is computed at a grid of values for the regularization parameter \code{lambda}.
#'
#'
#' @param y   a n-by-1 vector of responses
#' @param A   a n-by-1 vector of treatment variable; each element represents one of the L(>1) treatment conditions; e.g., c(1,2,1,1,3...); can be a factor-valued
#' @param X   a length-p list of functional-valued covariates, with its jth element corresponding to a n-by-n.eval[j] matrix of the observed jth functional covariates; n.eval[j] represents the number of evaluation points of the jth functional covariates
#' @param Z   a n-by-q matrix of scalar-valued covaraites
#' @param mu.hat   a n-by-1 vector of the fitted (X,Z)-main effect term of the model provided by the user; defult is \code{NULL}, in which case \code{mu.hat} is taken to be a vector of zeros; the optimal choice for this vector is E(y|X,Z)
#' @param d     number of basis spline functions to be used for each component function; the default value is 3; d=1 corresponds to the linear model
#' @param k     number of basis spline functions to be used for each single-index coefficient function associated with each functional covariate;
#' @param bs    type of basis for representing the single-index coefficient functions; the defult is "ps" (p-splines); any basis supported by mgcv::gam can be used, e.g., "cr" (cubic regression splines)
#' @param sp    smoothing parameter associated with the single-index coefficient function; the default is NULL, in which case the smoothing parameter is estimated based on generalized cross-validation
#' @param lambda.index  a user-supplied regularization parameter index to be used; the default is floor(nlambda/3).
#' @param lambda   a user-supplied regularization parameter sequence; typical usage is to have the program compute its own lambda sequence based on nlambda and lambda.min.ratio.
#' @param nlambda  total number of lambda values; the default value is 30.
#' @param lambda.min.ratio  the smallest value for lambda, as a fraction of lambda.max, the (data derived) entry value (i.e. the smallest value for which all coefficients are zero); the default is 0.01.
#' @param thol      stopping precision for the coordinate-descent algorithm; the default value is 1e-5.
#' @param max.ite   number of maximum iterations for the coordinate-descent procedure in fitting the component functions; the default value is 1e+5.
#' @param regfunc   type of the regularizer; the default is "L1"; can also be "MCP" or "SCAD".
#' @param max.iter  number of maximum iterations for the iterative procedure in fitting the single-index coefficient functions; the default value is 1e+1.
#' @param eps.iter  	   a value specifying the convergence criterion for the iterative procedure in fitting the single-index coefficient functions; the defult is 1e-2.
#' @param eps.num.deriv  a small value used in the finite difference method for computing the numerical (1st) derivatives of the estimated component functions; the default is 1e-4.
#' @param trace.iter if \code{TRUE}, trace the estimation process by printing, for each iteration, the difference from the previous iteration in the estimated single-index basis coefficients and the functional norms of the estimated component functions.
#'
#'
#'
#' @return a list of information of the fitted constrained functional additive models including
#'  \item{samTEMsel.obj}{an object of class \code{samTEMsel}, which contains the sequence of the set of fitted component functions implied by the sequence of the regularization parameters \code{lambda}; the sparse additive models are fitted over the set of the p functional covariates projected onto the estimated single-index coefficient functions (stored in \code{si.fit}) and the set of q scalar covariates; the object \code{samTEMsel.obj} includes the residuals of the fitted models and the fitted values for the response variable; see \code{samTEMsel::samTEMsel} for detail of the \code{samTEMsel} object.}
#'  \item{si.fit}{the length-p list of the single-index coefficient function estimate objects; each element is a \code{mgcv::gam} object; the jth element corresponds to the estimated single-index coefficient function associated with the jth functional covariate.}
#'  \item{si.coef.path}{the length-p list, where the jth element is a (iter-by-k) matrix, with the lth row corresponding to the basis coefficient vector estimate associated with the jth single-index coefficient function at the lth iteration of the fitting procedure.}
#'  \item{mean.fn}{the length-p list of mean functions (averaged across n observations), where the jth element is a n.eval[j]-by-1 vector of the evaluation of the estimated mean of the jth functional covariate.}
#'  \item{n.eval}{a length-p vector, where its jth element represents the number of evaluation points of the jth functional covariate.}
#'  \item{func_norm.record}{the iter-by-(p+q) matrix, with its lth row corresponding to the vector of the estimated (p+q) component functions' L2 norms at the lth iteration.}
#'  \item{func_norm}{a length (p+q) vector of the estimated (p+q) component functions' L2 norms, at the final iteration.}
#'  \item{lambda}{the sequence of regularization parameters used in the object \code{samTEMsel.obj}.}
#'  \item{lambda.index}{an index number, indicating the index of the regularization parameter in \code{lambda} used in obtaining the fitted model (including the single-index coefficient functions).}
#'  \item{nonzero.index}{a set of numbers, indicating the indices of estimated nonzero component functions of this particular fit under the regularization parameter index \code{lambda.index}.}
#'  \item{nonzero.X.index}{a set of numbers, indicating the indices of estimated nonzero component functions associated with the p functional covariates, based on this particular fit under the regularization parameter index \code{lambda.index}.}
#'
#'
#' @author Park, Petkova, Tarpey, Ogden
#' @import SAM samTEMsel
#' @seealso \code{cv.famTEMsel}, \code{predict_famTEMsel}, \code{plot_famTEMsel}, \code{make_ITR_famTEMsel}
#' @export
#'
#' @examples
#'
#' p = q = 10 # p and q are the numbers of functional and scalar pretreatment covariates, respectively.
#' n.train = 300  # training set sample size
#' n.test = 1000  # testing set sample size
#' n = n.train + n.test
#'
#' # generate p pretreatment functional covariates X by first seting up functional basis:
#' n.eval = 50; s = seq(0, 1, length.out = n.eval)  # a grid of support points
#' b1 = sqrt(2)*sin(2*pi*s)
#' b2 = sqrt(2)*cos(2*pi*s)
#' b3 = sqrt(2)*sin(4*pi*s)
#' b4 = sqrt(2)*cos(4*pi*s)
#' B  = cbind(b1, b2, b3, b4)   # a (n.eval-by-4) basis matrix
#' # randomly generate basis coefficients, and then add measurement noise
#' X = vector("list", length= p);  for(j in 1:p){
#'   X[[j]] = matrix(rnorm(n*4, 0, 1), n, 4) %*% t(B) +
#'     matrix(rnorm(n*n.eval, 0, 0.25), n, n.eval)    # measurement noise
#' }
#' Z  = matrix(rnorm(n*q, 0, 1), n, q)   # q scalar covariates
#' A  = rbinom(n, 1, 0.5) + 1  # treatment variable taking a value in {1,2} with equal prob.
#'
#' # X main effect on y; depends on the first 5 covariates
#' # the effect is generated randomly; randomly generated basis coefficients, scaled to unit L2 norm.
#' tmp = apply(matrix(rnorm(4*5), 4,5), 2, function(s) s/sqrt(sum(s^2) ))
#' main.effect = rep(0, n); for(j in 1:5){
#'   main.effect = main.effect + cos(X[[j]]%*% B%*%tmp[,j]/n.eval) # nonlinear effect (cosine)
#' }; rm(tmp)
#' # Z main effect on y; also depends on first 5 covariates
#' for(k in 1:5){
#'   main.effect = main.effect + cos(Z[,k])
#' }
#'
#' # define (interaction effect) coefficient functions associted with X[[1]] and X[[2]]
#' beta1 =  B %*% c(0.5,0.5,0.5,0.5)
#' beta2 =  B %*% c(0.5,-0.5,0.5,-0.5)
#' # A-by-X ineraction effect on y; depends only on X[[1]] and X[[2]].
#' interaction.effect = (A-1.5)*( 2*sin(X[[1]]%*%beta1/n.eval) + 2*sin(X[[2]]%*%beta2/n.eval))
#' # A-by-Z ineraction effect on y; depends only on Z[,1] and Z[,2].
#' interaction.effect =  interaction.effect + (A-1.5)*(Z[,1] + 2*sin(Z[,2]))
#'
#' # generate outcome y
#' noise = rnorm(n, 0, 0.5)
#' y = main.effect  + interaction.effect + noise
#'
#' var.main <- var(main.effect)
#' var.interaction <- var(interaction.effect)
#' var.noise <- var(noise)
#' SNR <- var.interaction/ (var.main + var.noise)
#' SNR  # "signal-to-noise" ratio
#'
#' # train/test set splitting
#' train.index = 1:n.train
#' y.train = y[train.index]
#' X.train= X.test = vector("list", p);  for(j in 1:p){
#'   X.train[[j]] = X[[j]][train.index,]
#'   X.test[[j]]  = X[[j]][-train.index,]
#' }
#' A.train = A[train.index]
#' A.test  = A[-train.index]
#' y.train = y[train.index]
#' y.test  = y[-train.index]
#' Z.train = Z[train.index,]
#' Z.test  = Z[-train.index,]
#'
#'
#' # fit a model with some regularization parameter index, say, lambda.index = 10.
#' # (an optimal regularization parameter can be estimated by running cv.famTEMsel().)
#' famTEMsel.obj = famTEMsel(y.train, A.train, X.train, Z.train, lambda.index=10)
#' famTEMsel.obj$func_norm  # L2 norm of the estimated component functions of the model
#' famTEMsel.obj$nonzero.index  # set of indices for the component functions estimated as nonzero
#' # plot the component functions estimated as nonzero and the single-index functions
#' plot_famTEMsel(famTEMsel.obj, which.index = famTEMsel.obj$nonzero.index)
#'
#' # make ITRs for subjects with pretreatment characteristics, X.test and Z.test
#' trt.rule = make_ITR_famTEMsel(famTEMsel.obj, newX = X.test, newZ = Z.test)$trt.rule
#' head(trt.rule)
#'
#' # an (IPWE) estimate of the "value" of this particualr treatment rule, trt.rule:
#' mean(y.test[A.test==trt.rule])
#'
#' # compare the above value to the following estimated "values" of "naive" treatment rules:
#' mean(y.test[A.test==1])   # a rule that assigns everyone to A=1
#' mean(y.test[A.test==2])   # a rule that assigns everyone to A=2
#'
# this function fits a constrained functional additive model for treatment effect-modifier selection
# given a fixed lambda.index (the sparsity parameter), fit a constrained funcitonal additive model
famTEMsel <- function(y, A, X, Z= NULL, mu.hat= NULL,
                      d = 3, k = 6, bs="ps", sp = NULL,
                      lambda = NULL, nlambda = 30, lambda.min.ratio = 1e-02,
                      lambda.index = floor(nlambda/3),
                      thol = 1e-05, max.ite = 1e+05, regfunc = "L1",
                      eps.iter=1e-02, max.iter=1e+1,
                      eps.num.deriv = 1e-04,
                      trace.iter = TRUE)
{

  n = length(y)
  p = length(X)
  q = ncol(Z); if(is.null(q)) q = 0
  pr.A = summary(as.factor(A))/n
  A = as.numeric(as.factor(A))  # so that the variable A takes a value in 1, 2, 3,..
  A.unique = unique(A)
  L = length(A.unique)

  n.eval = sapply(X, function(x) ncol(x))
  mean.fn = Xc = X.grid = vector("list", p); for(j in 1:p){
    mean.fn[[j]] = apply(X[[j]], 2, mean)
    Xc[[j]] = X[[j]] - matrix(mean.fn[[j]], n, n.eval[j], byrow=TRUE)   # center X
    Xc[[j]] = as.matrix(Xc[[j]])
    X.grid[[j]] = matrix(seq(0,1, length=n.eval[j]), n, n.eval[j], byrow=TRUE)
  }

  ## center Y (per each treatment level)
  #if(is.null(mu.hat))  mu.hat = rep(0, length(y))
  #if(is.null(mu.hat)){y.tilde = y} else y.tilde = y - mu.hat
  #y.mean = vector("list", L);  for(a in 1:L){
  #  y.mean[[a]] = mean(y.tilde[A==A.unique[a]])
  #  y.tilde = y.tilde - y.mean[[a]]*(A==A.unique[a])  # now, y.tilde does not have the A main effects in it
  #}

  # initialize the single-indices (si) (j=1,...,p) by taking scalar-averages
  si = vector("list", p); for(j in 1:p){
    si[[j]] = apply(Xc[[j]], 1, mean)  # simply take average of each function Xc[[j]] over its domain
  }

  # the initial single-index variables
  U = cbind(do.call(cbind, si), Z)

  # iterate until convergence (i.e., the single-index coefficients change less than eps) or reaching max.iter
  sdiff = 10^10   # set an arbitrary large number for the "starting difference" (sdiff)
  iter = 0
  diff = sdiff
  func_norm.record = NULL  # will record the L2 norms of the estimated component functions (across iter)
  si.fit = si.coef = si.coef.path = vector("list", p); for(j in 1:p){
    si.coef[[j]] = rep(1,k)/sqrt(sum(k^2))
  }

  while(diff > eps.iter){
    iter = iter + 1  # for each iteration, will optimize the constrained additive model by a coordinate descent, and then optimize U
    si.coef.old = si.coef

    # 1) Update samTEMsel.obj
    samTEMsel.obj = samTEMsel(y, A, U, mu.hat = mu.hat, d = d, lambda=lambda, nlambda = nlambda, thol = thol, max.ite = max.ite, regfunc = regfunc)
    tmp.func_norm = samTEMsel.obj$func_norm[, lambda.index]
    func_norm.record = rbind(func_norm.record, tmp.func_norm)  # record the L2 norms of the fitted component functions at lambda.index
    nonzero.X.index  = which(tmp.func_norm[1:p] > 1e-1 * max(tmp.func_norm))  # identify nonzero functional components (for this particular iter)
    rm(tmp.func_norm)

    # 2) compute/update si.coef
    for(j in nonzero.X.index){
      index = (j - 1) * (d*L) +  1:(d*L)
      g.der = matrix(0, n, 1)  # will store the 1st derivative of the jth component function here
      tmp = ns(samTEMsel.obj$X.tilde[,j] + eps.num.deriv,  # compute basis splines associated with the perturbed variable: samTEMsel.obj$X.tilde[,j] + eps.num.deriv
               df= samTEMsel.obj$d, knots= samTEMsel.obj$knots[,j],
               Boundary.knots= samTEMsel.obj$Boundary.knots[,j])
      # basis.eps is the basis spline model matrix associated with the perturbed variable: samTEMsel.obj$X.tilde[,j] + eps.num.deriv
      basis.eps = NULL; for(a in 1:L){
        basis.eps = cbind(basis.eps, tmp*samTEMsel.obj$A.indicator[,a])
      }
      # compute the 1st derivative of the basis splines (using finite-difference)
      basis.partial = (basis.eps - samTEMsel.obj$basis[,index])/eps.num.deriv
      g.der   = as.vector(basis.partial %*% samTEMsel.obj$w[index, lambda.index]) /samTEMsel.obj$X.ran[j]
      # adjusted responses, adjusted for the nonlinearity associated with the jth component function gj
      y.star  = samTEMsel.obj$residuals[,lambda.index] + g.der * si[[j]] #- mu.hat
      # adjusetd covariates, adjusted for the nonlinearity of the jth component function gj
      X.star  = diag(g.der) %*% Xc[[j]]
      X.grid.tmp   = X.grid[[j]]
      si.fit[[j]]  = gam(y.star ~ s(X.grid.tmp, by=X.star, k=k, bs= bs, sp=sp))
      si.coef.tmp  = si.fit[[j]]$coefficients[-1]
      si.coef.tmp  = si.coef.tmp/sqrt(sum(si.coef.tmp^2))
      if(si.coef.tmp[1] < 0)  si.coef.tmp = -si.coef.tmp   # for identifiability
      si.coef[[j]] = si.fit[[j]]$coefficients[-1] = si.coef.tmp
      si[[j]] =   predict(si.fit[[j]], type = "terms") /g.der  # update the jth single-index
      si.coef.path[[j]] = rbind(si.coef.path[[j]], si.coef[[j]])
      rm(tmp, basis.eps, basis.partial, g.der, y.star, X.star, X.grid.tmp, si.coef.tmp)
    }

    # Update the matrix U, based on the updated single-index si.
    U = cbind(do.call(cbind, si), Z)

    # Keep track of the iterative procedure of fitting the single-index basis coefficients
    ediff = 0;  for(j in nonzero.X.index){
      ediff = ediff + max(abs(si.coef[[j]] - si.coef.old[[j]]))  # the "end difference"
    }
    diff  = abs(sdiff - ediff)
    sdiff = ediff
    if(trace.iter){
      cat("iter", iter, "diff", diff, "\n")
      cat("functional norm: ", round(samTEMsel.obj$func_norm[, lambda.index], 2), "\n")
    }
    if(iter >= max.iter | length(nonzero.X.index) < 1)  break
  }

  func_norm = samTEMsel.obj$func_norm[, lambda.index]
  nonzero.index = which(func_norm > 1e-1 * max(func_norm))  # identify nonzero components

  res =  list(samTEMsel.obj = samTEMsel.obj,
              lambda= samTEMsel.obj$lambda,
              lambda.index  = lambda.index,
              nonzero.index = nonzero.index,
              nonzero.X.index = nonzero.X.index,
              func_norm = func_norm, func_norm.record = func_norm.record,
              si.fit = si.fit, si.coef.path= si.coef.path, mean.fn = mean.fn,
              n.eval=n.eval, n=n, p=p, q=q, L=L, X=X, Z=Z, A=A, y=y, iter = iter)
  class(res) = "famTEMsel"

  return(res)
}





#' \code{famTEMsel} prediction function
#'
#' \code{predict_famTEMsel} makes predictions given a (new) set of functional covariates \code{newX}, a (new) set of scalar covariates \code{newZ} and a (new) vector of treatment indicators \code{newA} based on a constrained functional additive model \code{famTEMsel.obj}.
#' Specifically, \code{predict_famTEMsel} predicts the responses y based on the (X,Z)-by-A interaction effect (plus the A main effect) portion of the full model that includes the unspecified X main effect term.
#'
#' @param famTEMsel.obj  a \code{famTEMsel} object
#' @param newX a (n by p) list of new values for the functional covariates X at which predictions are to be made; the jth element of the list corresponds to a n-by-n.eval[j] matrix of the observed jth functional covariates; n.eval[j] represents the number of evaluation points of the jth functional covariates; if \code{NULL}, X from the training set is used.
#' @param newZ a (n by q) matrix of new values for the scalar covariates Z at which predictions are to be made; if \code{NULL}, Z from the training set is used.
#' @param newA a (n by 1) vector of new values for the treatment A at which predictions are to be made; if \code{NULL}, A from the training set is used.
#' @param type the type of prediction required; the default "response" gives the predicted responses y based on the whole model; the alternative "terms" gives the component-wise predicted responses from each of the p components (and plus the treatment-specific intercepts) of the model.
#' @param lambda.index an index of the tuning parameter \code{lambda} at which predictions are to be made; one can supply \code{lambda.opt.index} obtained from the function \code{cv.samTEMsel}; the default is \code{NULL}, in which case the predictions based on the most non-sparse model is returned.
#'
#'
#' @return
#' \item{predicted}{a (n-by-\code{length(lambda.index)}) matrix of predicted values; if type =  "terms", then a (n-by-\code{length(lambda.index)}*(p+q+1)) matrix of predicted values, where the last column corresponds to the (treatment-specific) intercept.}
#' \item{U}{a n-by-(p+q) matrix of the index variables; the first p columns correspond to the 1-D projections of the p functional covariates and the last q columns correspond to the q scalar covariates.}
#'
#' @author Park, Petkova, Tarpey, Ogden
#' @seealso \code{famTEMsel},\code{cv.famTEMsel}, \code{plot_famTEMsel}
#' @export
#'
predict_famTEMsel = function(famTEMsel.obj, newX=NULL, newZ= NULL, newA=NULL,
                             type = "response",
                             lambda.index = NULL)
{
  if(!inherits(famTEMsel.obj, "famTEMsel"))   # checks input
    stop("object must be of class `famTEMsel'")

  lambda  = famTEMsel.obj$samTEMsel.obj$lambda
  nlambda = length(lambda)
  if (is.null(lambda.index)) {
    lambda.index = seq(nlambda)
  }else {
    lambda.index = intersect(lambda.index, seq(nlambda))
  }
  gcinfo(FALSE)
  if(is.null(newX)) newX = famTEMsel.obj$X
  if(is.null(newZ)) newZ = famTEMsel.obj$Z
  if(is.null(newA)) newA = famTEMsel.obj$A

  L  = famTEMsel.obj$samTEMsel.obj$L
  p  = famTEMsel.obj$p
  n  = nrow(newX[[1]])
  si = vector("list", p); for(j in 1:p){
    Xc.tmp     = newX[[j]] - matrix(famTEMsel.obj$mean.fn[[j]], n, famTEMsel.obj$n.eval[j], byrow=TRUE)
    if(is.null(famTEMsel.obj$si.fit[[j]])){
      si[[j]]  = apply(Xc.tmp, 1, mean)  # simply take average of each function Xc.tmp over its domain
    }else{
      X.grid.tmp = matrix(seq(0,1, length=famTEMsel.obj$n.eval[j]), n, famTEMsel.obj$n.eval[j], byrow=TRUE)
      dat        = list(X.grid.tmp = X.grid.tmp, X.star = as.matrix(Xc.tmp))
      si[[j]]    = predict(famTEMsel.obj$si.fit[[j]], newdata= dat, type = "terms")[,1]
      rm(dat, X.grid.tmp)
    }
    rm(Xc.tmp)
  }

  U = cbind(do.call(cbind, si), newZ)
  U.min.rep = matrix(rep(famTEMsel.obj$samTEMsel.obj$X.min, n), nrow=n, byrow=TRUE)
  U.ran.rep = matrix(rep(famTEMsel.obj$samTEMsel.obj$X.ran, n), nrow=n, byrow=TRUE)
  U.tilde = (U-U.min.rep)/U.ran.rep
  U.tilde[U.tilde < 0] = 0.01
  U.tilde[U.tilde > 1] = 0.99

  A.unique = unique(famTEMsel.obj$samTEMsel.obj$A)
  A.indicator = matrix(0, n, L); for(a in 1:L){
    A.indicator[,a] = as.numeric(newA == A.unique[a])
  }

  d = famTEMsel.obj$samTEMsel.obj$d
  v = famTEMsel.obj$samTEMsel.obj$p
  basis <- matrix(0,n, d*v*L)
  for(j in 1:v){  # construct the jth model matrix
    index = (j - 1) * d*L + c(1:(d*L))
    tmp0 = ns(U.tilde[,j], df= d,
              knots= famTEMsel.obj$samTEMsel.obj$knots[,j],
              Boundary.knots= famTEMsel.obj$samTEMsel.obj$Boundary.knots[,j])
    tmp1 = NULL
    for(a in 1:L){
      tmp1 = cbind(tmp1, tmp0*A.indicator[,a])
    }
    basis[,index] = tmp1
  }

  # compute the treatment-specific intercepts
  intercepts = NULL; for(a in 1:L){
    intercepts = rbind(intercepts, famTEMsel.obj$samTEMsel.obj$intercept[[a]][, lambda.index])
  }

  # make predictions
  if(type=="terms"){
      y.hat.list = vector("list", length = v);  for (j in 1:v){
      index = (j - 1) *d*L + c(1:(d*L))
      y.hat.list[[j]] = basis[,index] %*% as.matrix(famTEMsel.obj$samTEMsel.obj$w[index, lambda.index])
    }
    predicted = cbind(do.call("cbind", y.hat.list), A.indicator %*% intercepts)
  }else{
    predicted = cbind(basis, A.indicator) %*% rbind(as.matrix(famTEMsel.obj$samTEMsel.obj$w[, lambda.index]), intercepts)
  }

  res = list(predicted = predicted, U = U)
  return(res)
}




#' make individualized treatment recommendations (ITRs) based on a \code{famTEMsel} object
#'
#' The function \code{make_ITR_famTEMsel} returns individualized treatment decision recommendations for subjects with pretreatment characteristics \code{newX} and \code{newZ}, given a \code{famTEMsel} object, \code{famTEMsel.obj}, and an (optimal) regularization parameter index, \code{lambda.index}.
#'
#'
#' @param famTEMsel.obj  a \code{famTEMsel} object, containing the fitted constrained functional additive models.
#' @param newX a length-p list of new values for the functional-valued covariates X, where the jth element is a (n-by-n.eval[j]) matrix of the observed jth function, at which predictions are to be made; if \code{NULL}, X from the training set is used.
#' @param newZ a (n-by-q) matrix of new values for the scalar-valued covariates Z at which predictions are to be made; if \code{NULL}, Z from the training set is used.
#' @param lambda.index an index of the regularization parameter \code{lambda} at which predictions are to be made; one can supply \code{lambda.opt.index} obtained from the function \code{cv.famTEMsel()}; the default is \code{NULL}, in which case the predictions are made at the \code{lambda.index} used in obtaining \code{famTEMsel.obj}.
#' @param maximize default is \code{TRUE}, assuming a larger value of the outcome is better; if \code{FALSE}, a smaller value is assumed to be prefered.
#'
#' @return
#' \item{pred.new}{a (n-by-L) matrix of predicted values, with each column representing one of the L treatment options.}
#' \item{trt.rule}{a (n-by-1) vector of the individualized treatment recommendations}
#'
#' @author Park, Petkova, Tarpey, Ogden
#' @seealso \code{famTEMsel},\code{cv.famTEMsel}, \code{predict_famTEMsel}
#' @export
#'
# make treatment decisions given newX and nexZ, based on the estimated sequence of models, famTEMsel.obj, and a given regularization parameter index, lambda.index.
make_ITR_famTEMsel  = function(famTEMsel.obj, newX = NULL, newZ= NULL, lambda.index = NULL, maximize=TRUE){

  if(!inherits(famTEMsel.obj, "famTEMsel"))   # checks input
    stop("famTEMsel.obj must be of class `famTEMsel'")

  if(is.null(lambda.index))  lambda.index = famTEMsel.obj$lambda.index

  if(is.null(newX) & is.null(newZ)){
    newX = famTEMsel.obj$X
    newZ = famTEMsel.obj$Z
  }
  L = famTEMsel.obj$L
  n = nrow(newX[[1]])

  # compute potential outcome for each of the L treatment conditions, given newX and newZ.
  pred.new= matrix(0, nrow=n, ncol = L);  for(a in 1:L){
    pred.new[,a] = predict_famTEMsel(famTEMsel.obj, newX= newX, newZ=newZ, newA=rep(a, n),
                                    lambda.index = lambda.index)$predicted
  }
  # compute treatment assignment
  if(maximize){
    trt.rule = apply(pred.new, 1, which.max)
  }else{
    trt.rule = apply(pred.new, 1, which.min)
  }
  if(L==2)  colnames(pred.new) <- c("Tr1", "Tr2")

  return(list(trt.rule = trt.rule, pred.new = pred.new))
}






#' plot component functions from a \code{famTEMsel} object
#'
#' Produces plots of the component functions and the single-index coefficient functions from a \code{famTEMsel} object.
#'
#' @param famTEMsel.obj  a \code{famTEMsel} object
#' @param newX a (n by p) list of new values for the functional covariates X at which predictions are to be made; the jth element of the list corresponds to a n-by-n.eval[j] matrix of the observed jth functional covariates; n.eval[j] represents the number of evaluation points of the jth functional covariates; if \code{NULL}, X from the training set is used.
#' @param newZ a (n by q) matrix of new values for the scalar covariates Z at which predictions are to be made; if \code{NULL}, Z from the training set is used.
#' @param newA a (n-by-1) vector of new values for the treatment A at which plots are to be made; the default is \code{NULL}, in which case A is taken from the training set.
#' @param lambda.index an index of the tuning parameter \code{lambda} at which plots are to be made; one can supply \code{lambda.opt.index} obtained from the function \code{cv.samTEMsel}; the default is \code{NULL}, in which case \code{plot_samTEMsel} utilizes the most non-sparse model.
#' @param which.index  this specifies which component functions are to be plotted; the default is only nonzero component functions.
#' @param scatter.plot if \code{TRUE}, draw scatter plots of partial residuals versus the covariates; these scatter plots are made based on the training observations;  the default is \code{TRUE}.
#' @param ylims this specifies the vertical range of the plots, e.g., c(-10, 10).
#' @param single.index.plot if \code{TRUE}, draw the plots of the estimated single-index coefficient functions; the default is \code{FALSE}.
#' @param solution.path  if \code{TRUE}, draw the functional norms of the fitted component functions (based on the training set) versus the regularization parameter; the default is \code{FALSE}.
#'
#'
#' @author Park, Petkova, Tarpey, Ogden
#' @seealso \code{famTEMsel},\code{predict_famTEMsel}, \code{cv.famTEMsel}
#' @export
#'
plot_famTEMsel <- function(famTEMsel.obj,
                           newX = NULL, newZ = NULL, newA=NULL,
                           scatter.plot = TRUE,
                           lambda.index =famTEMsel.obj$lambda.index,
                           which.index  =famTEMsel.obj$nonzero.index,
                           ylims,
                           single.index.plot=FALSE, solution.path = FALSE)
{

  v = famTEMsel.obj$samTEMsel.obj$p
  L = famTEMsel.obj$L

  # compute component-wise predictions, given newX and newA
  pred.new  = predict_famTEMsel(famTEMsel.obj, newX=newX, newZ = newZ, newA=newA, type = "terms", lambda.index = lambda.index)
  predicted.new = pred.new$predicted  # the matrix of component-wise predicted responses
  newU  = pred.new$U   # the matrix of single-indices
  if(is.null(newA))  newA = famTEMsel.obj$A

  if(scatter.plot){
    pred.train = predict_famTEMsel(famTEMsel.obj, type = "terms", lambda.index = lambda.index)
    predicted.train = pred.train$predicted
    U.train = pred.train$U
    resid.list =  vector("list", length = v);  for(j in 1:v){  # compute the jth partial residuals
      resid.list[[j]] = famTEMsel.obj$samTEMsel.obj$residuals[, lambda.index] + predicted.train[,j]
    }
    A.train = famTEMsel.obj$A
  }

  if(missing(ylims)){
    if(scatter.plot){
      ylims =  range(do.call("cbind", resid.list))
    }else{
      ylims = range(predicted.new[, which.index])
    }
  }

  colvals    = 1:L
  colvals[1] = "royalblue3"
  colvals[L] = "tomato3"
  if(L > 2) colvals[2] = "darkseagreen"

  for(j in which.index){
    o = order(newU[newA==1, j])
    plot(matrix(newU[newA==1,],ncol=v)[o, j], predicted.new[newA==1, j][o], type = "l",
         ylab = paste("f(X", j, ")", sep = ""), xlab = paste("X", j, sep = ""),
         ylim = ylims, c(range(newU[,j])[1], range(newU[,j])[2]),
         col = colvals[1], lwd = 2)
    if(scatter.plot)  points(U.train[A.train==1,j],
                             resid.list[[j]][A.train==1],
                             col = colvals[1])
    if(L > 1){
      for(a in 2:L){
        o = order(newU[newA==a, j])
        lines(matrix(newU[newA==a,],ncol=v)[o, j], predicted.new[newA==a, j][o], col = colvals[a], lwd = 2)
        if(scatter.plot)  points(U.train[A.train==a,j],
                                 resid.list[[j]][A.train==a],
                                 col = colvals[a], pch = a, lwd = 1)
      }
    }
  }

  if(single.index.plot){
    for (j in intersect(famTEMsel.obj$nonzero.X.index, which.index)) {
      plot(famTEMsel.obj$si.fit[[j]], shade=TRUE,
           ylab = paste("beta", j, "(s)", sep = ""), xlab = "s")
    }
  }

  if(solution.path){
    par = par(omi = c(0.0, 0.0, 0, 0), mai = c(1, 1, 0.1, 0.1))
    matplot(famTEMsel.obj$samTEMsel.obj$lambda[length(famTEMsel.obj$samTEMsel.obj$lambda):1],
            t(famTEMsel.obj$samTEMsel.obj$func_norm), type="l",
            xlab=expression(paste("Regularization parameter ", lambda)),
            ylab = "Functional norms",cex.lab=1,log="x",lwd=2)
  }

}


######################################################################
## END OF THE FILE
######################################################################

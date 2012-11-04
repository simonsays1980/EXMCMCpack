"MCMCFinMixmetrop1R" <- function(fun, theta.init, weights.start, burnin = 500, mcmc = 20000,
                                 thin = 1, tune = 1, verbose = 0, seed = NA, developer = FALSE, 
                                 logfun = TRUE, force.samp = FALSE, V = NULL, 
                                 optim.control = list(fnscale = -1,
                                 trace = 0, REPORT = 10, maxit = 500), 
                                 hin = NULL, hin.jac = NULL, control.outer = list(optim.method = "BFGS", 
                                 optim.lower = -Inf, optim.upper = Inf, fnscale = -1), ...) {
  ## packages
  require(alabama)
  
  ## error checking here
  MCMCpack:::check.offset(list(...))
  MCMCpack:::check.mcmc.parameters(burnin, mcmc, tune)
  
  ## form the tuning vector
  tune <- MCMCpack:::vector.tune(tune, length(theta.init))
  
  ## form seed 
  seeds <- MCMCpack:::form.seeds(seed)
  lecuyer <- seeds[[1]]
  seed.array <- seeds[[2]]
  lecuyer.stream <- seeds[[3]]
  
  ## setup the environment so that fun can see the things passed as ...
  userfun <- function(ttt) fun(ttt, ...)
  my.env <- environment(fun = userfun)
  
  ## setup function for maximization based on value of logfun
  if ( logfun ) {
    maxfun <- userfun
  }
  else if ( logfun == FALSE ) {
    maxfun <- function(ttt, ...) log(fun(ttt, ...))
  }
  else {
    cat("logfun not a logical value.\n")
    stop("Respecify and call MCMCFinMixmetrop1R() again. \n",
         call. = FALSE)
  }
  
  ## hessian matrix
  if( is.null(V) ) {
    ## find approx mode and Hessian using auglag()
    opt.out <- auglag(theta.init, maxfun, my.env, hin = hin, hin.jac = hin.jac, 
                      control.optim = optim.control)

    if( opt.out$convergence != 0) { ## no convergence
      warning("Mode and Hessian were not found with call to auglag().\nSampling proceeded anyway. \n")
    }
    
    CC <- NULL
    try(CC <- chol(-1 * opt.out$hessian), silent = TRUE)
    hess.new <- opt.out$hessian
    hess.flag <- 0
   
    if ( force.samp == TRUE ) {
      ## enforce negative definiteness 
      if ( max(diag(opt.out$hessian) == 0) ) {
        for (i in 1:nrow(hess.new) ) {
          if ( hess.new[i, i] == 0 ) {
            hess.new[i, i] <- -1e-6
          }
        }
      }
      
      while (is.null(CC)) {
        hess.flag <- 1
        hess.new <- hess.new - diag(diag(0.01 * abs(opt.out$hessian)))
        try(CC <- chol(-1 * hess.new), silent = TRUE)
      }
    }
    
    else { ## force.samp == FALSE
      if ( is.null(CC) ) {
        hess.flag <- 2
      }
      
    }
    
    if ( hess.flag == 1) {
      warning("Hessian from call to auglag() not negative definite.\nSampling proceeded after enforcing
              negative definiteness. \n")
    }
    
    if ( hess.flag == 2 ) {
      cat("Hessian from call to auglag() not negative definite.\n")
      cat("Sampling (as.specified) cannot proceed.\n")
      stop("Check data and fun() and call MCMCFinMixmetrop1R() again. \n", 
           call. = FALSE)
    }
    
    V <- tune %*% solve(-1 * hess.new) %*% tune
  }
  
  else { ## V defined by user
    if( nrow(V) != ncol(V) || nrow(V) != (length(theta.init) + 1) ) {
      cat("V not of appropriate dimension.\n")
      stop("Check V and theta.init and call MCMCFinMixmetrop1R() again. \n",
           call. = FALSE)
    }
    
    CC <- NULL
    try(CC <- chol(V), silent = TRUE) 
    if ( is.null(CC) ) {
      cat("V not positive definite.\n")
      stop("Check V and call MCMCFinMixmetrop1R() again. \n",
           call. = FALSE)
    } 
    
    V <- tune %*% V %*% tune
  }
  
  ## Call the C++ function to do the MCMC sampling
  sample <- .Call("MCMCFinMixmetrop1R_cc", userfun, as.double(theta.init),
                  as.integer(weights.start), my.env, as.integer(burnin), as.integer(mcmc), 
                  as.integer(thin), as.integer(verbose),
                  lecuyer = as.integer(lecuyer),
                  seedarray = as.integer(seed.array),
                  lecuyerstream = as.integer(lecuyer.stream),
                  as.logical(logfun), 
                  as.matrix(V),
                  PACKAGE = "MCMCpack")
  
  ## turn sample into an mcmc object
  sample <- mcmc(data = sample, start = burnin + 1, end = burnin + mcmc, thin = thin)
  return(sample)
}
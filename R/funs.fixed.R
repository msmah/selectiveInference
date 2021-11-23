# Lasso inference function (for fixed lambda). Note: here we are providing inference
# for the solution of
# min 1/2 || y - \beta_0 - X \beta ||_2^2 + \lambda || \beta ||_1

fixedLassoInf <- function(x, y, beta, 
                          lambda, 
                          family=c("gaussian","binomial","cox"),
                          intercept=TRUE, 
                          add.targets=NULL, 
                          status=NULL,
                          sigma=NULL, 
                          alpha=0.1,
                          type=c("partial", "full"), 
                          tol.beta=1e-5, 
                          tol.kkt=0.1,
                          gridrange=c(-100,100), 
			  bits=NULL, 
			  verbose=FALSE, 
                          linesearch.try=10) {

  family = match.arg(family)
  this.call = match.call()
  type = match.arg(type)
  if(family=="binomial")  {
    if(type!="partial") stop("Only type= partial allowed with binomial family")
    out=fixedLogitLassoInf(x,y,beta,lambda,alpha=alpha, type=type, tol.beta=tol.beta, tol.kkt=tol.kkt,
                           gridrange=gridrange, bits=bits, verbose=verbose,
                           linesearch.try=linesearch.try,
                           this.call=this.call)
    return(out)
  }
  else if(family=="cox")  {
    if(type!="partial") stop("Only type= partial allowed with Cox family")
    out=fixedCoxLassoInf(x,y,status,beta,lambda,alpha=alpha, type="partial",tol.beta=tol.beta,
                         tol.kkt=tol.kkt, gridrange=gridrange, bits=bits, verbose=verbose,this.call=this.call)
    return(out)
  }
  
  else{
    
    checkargs.xy(x,y)
    if (missing(beta) || is.null(beta)) stop("Must supply the solution beta")
    if (missing(lambda) || is.null(lambda)) stop("Must supply the tuning parameter value lambda")
    
    n = nrow(x)
    p = ncol(x)
    beta = as.numeric(beta)
    if (type == "full") {
      if (p > n) {
        # need intercept (if there is one) for debiased lasso
        hbeta = beta
        if (intercept == T) {
          if (length(beta) != p + 1) {
            stop("Since type='full', p > n, and intercept=TRUE, beta must have length equal to ncol(x)+1")
          }
          # remove intercept if included
          beta = beta[-1]
        } else if (length(beta) != p) {
          stop("Since family='gaussian', type='full' and intercept=FALSE, beta must have length equal to ncol(x)")
        }
      }
    } else if (length(beta) != p) {
      stop("Since family='gaussian' and type='partial', beta must have length equal to ncol(x)")
    }
    
    checkargs.misc(beta=beta,lambda=lambda,sigma=sigma,alpha=alpha,
                   gridrange=gridrange,tol.beta=tol.beta,tol.kkt=tol.kkt)
    if (!is.null(bits) && !requireNamespace("Rmpfr",quietly=TRUE)) {
      warning("Package Rmpfr is not installed, reverting to standard precision")
      bits = NULL
    }
    
    if (!is.null(add.targets) && (!is.vector(add.targets)
                                  || !all(is.numeric(add.targets)) || !all(add.targets==floor(add.targets))
                                  || !all(add.targets >= 1 && add.targets <= p))) {
      stop("'add.targets' must be a vector of integers between 1 and p")
    }
    
    # If glmnet was run with an intercept term, center x and y
    if (intercept==TRUE) {
      obj = standardize(x,y,TRUE,FALSE)
      x = obj$x
      y = obj$y
    }
    
    # Check the KKT conditions
    g = t(x)%*%(y-x%*%beta) / lambda
    if (any(abs(g) > 1+tol.kkt * sqrt(sum(y^2))))
      warning(paste("Solution beta does not satisfy the KKT conditions",
                    "(to within specified tolerances)"))
    
    tol.coef = tol.beta * sqrt(n / colSums(x^2))
    vars = which(abs(beta) > tol.coef)
    sign_vars = sign(beta[vars])

    if(sum(vars)==0){
      cat("Empty model",fill=T)
      return()
    }

    if (any(sign(g[vars]) != sign_vars)) {
      warning(paste("Solution beta does not satisfy the KKT conditions",
                    "(to within specified tolerances). You might try rerunning",
                    "glmnet with a lower setting of the",
                    "'thresh' parameter, for a more accurate convergence."))
    }

    # Get lasso polyhedral region, of form Gy >= u

    logical.vars=rep(FALSE,p)
    logical.vars[vars]=TRUE
    
    if (type == 'full') {
       out = fixedLassoPoly(x, y, lambda, beta, logical.vars, inactive=TRUE)
    } 
    else {
       out = fixedLassoPoly(x, y, lambda, beta, logical.vars)
    }
    
    A = out$A
    b = out$b
    
    # Check polyhedral region
    tol.poly = 0.01
    if (max(A %*% y - b) > tol.poly * sqrt(sum(y^2)))
      stop(paste("Polyhedral constraints not satisfied; you must recompute beta",
                 "more accurately. With glmnet, make sure to use exact=TRUE in coef(),",
                 "and check whether the specified value of lambda is too small",
                 "(beyond the grid of values visited by glmnet).",
                 "You might also try rerunning glmnet with a lower setting of the",
                 "'thresh' parameter, for a more accurate convergence."))
    
    # Estimate sigma
    if (is.null(sigma)) {
      if (n >= 2*p) {
        oo = intercept
        sigma = sqrt(sum(lsfit2(x,y,intercept=oo)$res^2)/(n-p-oo))
      }
      else {
        sigma = sd(y)
        warning(paste(sprintf("p > n/2, and sd(y) = %0.3f used as an estimate of sigma;",sigma),
                      "you may want to use the estimateSigma function"))
      }
    }
    
    # add additional targets for inference if provided
    if (!is.null(add.targets)) {
       # vars is boolean...
       old_vars = vars & TRUE
       vars[add.targets] = TRUE
       sign_vars = sign(beta[vars]) 
       sign_vars[!old_vars] = NA
       stop("`add.targets` not fully implemented yet")
    }

    k = length(vars)
    pv = vlo = vup = numeric(k)
    vmat = matrix(0,k,n)
    ci = tailarea = matrix(0,k,2)
      
    if (type=="full" & p > n) {
      if (intercept == TRUE) {
        pp=p+1
        Xint <- cbind(rep(1,n),x)
        # indices of selected predictors
        S = c(1,vars + 1)
      } else {
        pp=p
        Xint <- x
        # indices of selected predictors
        S = vars
        # notS = which(abs(beta) <= tol.coef)
      }
      
      notS = setdiff(1:pp,S)
      
      XS = Xint[,S]
      hbetaS = hbeta[S]
      
      # Reorder so that active set S is first
      Xordered = Xint[,c(S,notS,recursive=T)]
      hsigmaS = 1/n*(t(XS)%*%XS) # hsigma[S,S]
      hsigmaSinv = solve(hsigmaS) # pinv(hsigmaS)
      
      FS = rbind(diag(length(S)),matrix(0,pp-length(S),length(S)))
      GS = cbind(diag(length(S)),matrix(0,length(S),pp-length(S)))

      is_wide = n < (2 * p) # somewhat arbitrary decision -- it is really for when we don't want to form with pxp matrices

      # Approximate inverse covariance matrix for when (n < p) from lasso_Inference.R
      if (!is_wide) {
           hsigma = 1/n*(t(Xordered)%*%Xordered)
           htheta = debiasingMatrix(hsigma, 
                                    is_wide, 
                                    n, 
                                    1:length(S), 
                                    verbose=FALSE, 
                                    max_try=linesearch.try, 
                                    warn_kkt=TRUE)
           ithetasigma = (GS-(htheta%*%hsigma))
      } else {
           htheta = debiasingMatrix(Xordered, 
                                    is_wide, 
                                    n, 
                                    1:length(S), 
                                    verbose=FALSE, 
                                    max_try=linesearch.try, 
                                    warn_kkt=TRUE)
           ithetasigma = (GS-((htheta%*%t(Xordered)) %*% Xordered)/n)
      }

      M <- (((htheta%*%t(Xordered))+ithetasigma%*%FS%*%hsigmaSinv%*%t(XS))/n)

      # vector which is offset for testing debiased beta's
      null_value <- (((ithetasigma%*%FS%*%hsigmaSinv)%*%sign(hbetaS))*lambda/n)

      if (intercept == T) {
        M = M[-1,] # remove intercept row
        null_value = null_value[-1] # remove intercept element
      }
    } else if (type=="partial" || p > n) {
      xa = x[,vars,drop=F]
      M = pinv(crossprod(xa)) %*% t(xa)
      null_value = rep(0,k)
    } else {
      M = pinv(crossprod(x)) %*% t(x)
      M = M[vars,,drop=F]
      null_value = rep(0,k)
    }

  for (j in 1:k) {
    if (verbose) cat(sprintf("Inference for variable %i ...\n",vars[j]))

    vj = M[j,]
    mj = sqrt(sum(vj^2))
    vj = vj / mj        # Standardize (divide by norm of vj)

    if (!is.na(sign_vars[j])) {
        vj = sign_vars[j] * vj
    }

    limits.info = TG.limits(y, A, b, vj, Sigma=diag(rep(sigma^2, n)))
    a = TG.pvalue.base(limits.info, null_value=null_value[j], bits=bits)
    pv[j] = a$pv
    
    if (is.na(sign_vars[j])) { # for variables not in the active set, report 2-sided pvalue
       pv[j] = 2 * min(pv[j], 1 - pv[j])
    }
    
    vlo[j] = a$vlo * mj # Unstandardize (mult by norm of vj)
    vup[j] = a$vup * mj # Unstandardize (mult by norm of vj)
    if (!is.na(sign_vars[j])) { 
        vmat[j,] = vj * mj * sign_vars[j]  # Unstandardize (mult by norm of vj) and fix sign
    } else {
        vmat[j,] = vj * mj # Unstandardize (mult by norm of vj)
    }
    a = TG.interval.base(limits.info, 
                         alpha=alpha,
                         gridrange=gridrange,
			 flip=(sign_vars[j]==-1),
                         bits=bits)
    ci[j,] = (a$int-null_value[j]) * mj # Unstandardize (mult by norm of vj)
    tailarea[j,] = a$tailarea
  }

  out = list(type=type,
             lambda=lambda,
             pv=pv,
             ci=ci,
             tailarea=tailarea,
             vlo=vlo,
             vup=vup,
             vmat=vmat,
             y=y,
             vars=vars,
             sign=sign_vars,
             sigma=sigma,
             alpha=alpha,
             sd=sigma*sqrt(rowSums(vmat^2)),
             coef0=vmat%*%y,
             call=this.call)

  class(out) = "fixedLassoInf"
  return(out)
}
}

#############################
#  File src/library/stats/R/lsfit.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 1995-2016 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  https://www.R-project.org/Licenses/

lsfit2 <- function(x, y, wt = NULL, intercept = TRUE, tolerance = 1e-07,
                  yname = NULL)
{
    ## find names of x variables (design matrix)

    x <- as.matrix(x)
    y <- as.matrix(y)
    xnames <- colnames(x)
    if( is.null(xnames) ) {
	if(ncol(x) == 1L) xnames <- "X"
	else xnames <- paste0("X", 1L:ncol(x))
    }
    if( intercept ) {
	x <- cbind(1, x)
	xnames <- c("Intercept", xnames)
    }

    ## find names of y variables (responses)

    if(is.null(yname) && ncol(y) > 1) yname <- paste0("Y", 1L:ncol(y))

    ## remove missing values

    #good <- complete.cases(x, y, wt)
    #dimy <- dim(as.matrix(y))
    #if( any(!good) ) {
    #    warning(sprintf(ngettext(sum(!good),
    #                             "%d missing value deleted",
    #                             "%d missing values deleted"),
    #                    sum(!good)), domain = NA)
	#x <- as.matrix(x)[good, , drop=FALSE]
	#y <- as.matrix(y)[good, , drop=FALSE]
	#wt <- wt[good]
    #}

    ## check for compatible lengths

    nrx <- NROW(x)
    ncx <- NCOL(x)
    nry <- NROW(y)
    ncy <- NCOL(y)
    nwts <- length(wt)
    if(nry != nrx)
        stop(sprintf(paste0(ngettext(nrx,
                       "'X' matrix has %d case (row)",
                       "'X' matrix has %d cases (rows)"),
              ", ",
              ngettext(nry,
                       "'Y' has %d case (row)",
                       "'Y' has %d cases (rows)")),
                       nrx, nry),
                       domain = NA)
    if(nry < ncx)
        stop(sprintf(paste0(ngettext(nry,
                              "only %d case",
                              "only %d cases"),
                     ", ",
                     ngettext(ncx,
                              "but %d variable",
                              "but %d variables")),
                     nry, ncx),
             domain = NA)
    ## check weights if necessary
    if( !is.null(wt) ) {
	if(any(wt < 0)) stop("negative weights not allowed")
	if(nwts != nry)
            stop(gettextf("number of weights = %d should equal %d (number of responses)", nwts, nry), domain = NA)
	wtmult <- sqrt(wt)
	if(any(wt == 0)) {
	    xzero <- as.matrix(x)[wt == 0, ]
	    yzero <- as.matrix(y)[wt == 0, ]
	}
	x <- x*wtmult
	y <- y*wtmult
	invmult <- 1/ifelse(wt == 0, 1, wtmult)
    }

    # Here y is a matrix, so z$residuals and z$effects will be
    z <- .Call(C_Cdqrls, x, y, tolerance, FALSE)

    resids <- array(NA, dim = dimy)
    dim(z$residuals) <- c(nry, ncy)
    if(!is.null(wt)) {
	if(any(wt == 0)) {
	    if(ncx == 1L) fitted.zeros <- xzero * z$coefficients
	    else fitted.zeros <- xzero %*% z$coefficients
	    z$residuals[wt == 0, ] <- yzero - fitted.zeros
	}
	z$residuals <- z$residuals*invmult
    }
    resids[good, ] <- z$residuals
    if(dimy[2L] == 1 && is.null(yname)) {
	resids <- drop(resids)
	names(z$coefficients) <- xnames
    } else {
	colnames(resids) <- yname
	colnames(z$effects) <- yname
	dim(z$coefficients) <- c(ncx, ncy)
	dimnames(z$coefficients) <- list(xnames, yname)
    }
    z$qr <- as.matrix(z$qr)
    colnames(z$qr) <- xnames
    output <- list(coefficients = z$coefficients, residuals = resids)

    ## if X matrix was collinear, then the columns may have been
    ## pivoted hence xnames may need to be corrected

    if( z$rank != ncx ) {
	xnames <- xnames[z$pivot]
	dimnames(z$qr) <- list(NULL, xnames)
	warning("'X' matrix was collinear")
    }

    ## return weights if necessary

    if (!is.null(wt) ) {
	weights <- rep.int(NA, dimy[1L])
	weights[good] <- wt
	output <- c(output, list(wt=weights))
    }

    ## return rest of output

    ## Neither qt nor tol are documented to be there.
    rqr <- list(qt = drop(z$effects), qr = z$qr, qraux = z$qraux, rank = z$rank,
		pivot = z$pivot, tol = z$tol)
    oldClass(rqr) <- "qr"
    output <- c(output, list(intercept = intercept, qr = rqr))
    return(output)
}

ls.diag <- function(ls.out)
{
    resids <- as.matrix(ls.out$residuals)
    d0 <- dim(resids)
    xnames <- colnames(ls.out$qr$qr)
    yname <- colnames(resids)

    ## remove any missing values

    good <- complete.cases(resids, ls.out$wt)
    if( any(!good) ) {
	warning("missing observations deleted")
	resids <- resids[good, , drop = FALSE]
    }

    ## adjust residuals if needed

    if( !is.null(ls.out$wt) ) {
	if( any(ls.out$wt[good] == 0) )
	    warning("observations with 0 weight not used in calculating standard deviation")
	resids <- resids * sqrt(ls.out$wt[good])
    }

    ## initialize

    p <- ls.out$qr$rank
    n <- nrow(resids)
    hatdiag <- rep.int(NA, n)
    stats <- array(NA, dim = d0)
    colnames(stats) <- yname
    stdres <- studres <- dfits <- Cooks <- stats

    ## calculate hat matrix diagonals

    q <- qr.qy(ls.out$qr, rbind(diag(p), matrix(0, nrow=n-p, ncol=p)))
    hatdiag[good] <- rowSums(as.matrix(q^2))

    ## calculate diagnostics

    stddev <- sqrt(colSums(as.matrix(resids^2))/(n - p))
    stddevmat <- matrix(stddev, nrow=sum(good), ncol=ncol(resids), byrow=TRUE)
    stdres[good, ] <- resids/(sqrt(1-hatdiag[good]) * stddevmat)
    studres[good, ] <- (stdres[good, ]*stddevmat) /
        sqrt(((n-p)*stddevmat^2 - resids^2/(1-hatdiag[good]))/(n-p-1))
    dfits[good, ] <- sqrt(hatdiag[good]/(1-hatdiag[good])) * studres[good, ]
    Cooks[good, ] <- ((stdres[good, ]^2 * hatdiag[good])/p)/(1-hatdiag[good])
    if(ncol(resids)==1 && is.null(yname)) {
	stdres <- as.vector(stdres)
	Cooks <- as.vector(Cooks)
	studres <- as.vector(studres)
	dfits <- as.vector(dfits)
    }

    ## calculate unscaled covariance matrix

    qr <- as.matrix(ls.out$qr$qr[1L:p, 1L:p])
    qr[row(qr)>col(qr)] <- 0
    qrinv <- solve(qr)
    covmat.unscaled <- qrinv%*%t(qrinv)
    dimnames(covmat.unscaled) <- list(xnames, xnames)

    ## calculate scaled covariance matrix

    covmat.scaled <- sum(stddev^2) * covmat.unscaled

    ## calculate correlation matrix

    cormat <- covmat.scaled /
	sqrt(outer(diag(covmat.scaled), diag(covmat.scaled)))

    ## calculate standard error

    stderr <- outer(diag(covmat.unscaled)^0.5, stddev)
    dimnames(stderr) <- list(xnames, yname)

    return(list(std.dev=stddev, hat=hatdiag, std.res=stdres,
		stud.res=studres, cooks=Cooks, dfits=dfits,
		correlation=cormat, std.err=stderr,
		cov.scaled=covmat.scaled, cov.unscaled=covmat.unscaled))
}

ls.print <- function(ls.out, digits = 4L, print.it = TRUE)
{
    ## calculate residuals to be used

    resids <- as.matrix(ls.out$residuals)
    if( !is.null(ls.out$wt) ) {
	if(any(ls.out$wt == 0))
	    warning("observations with 0 weights not used")
	resids <- resids * sqrt(ls.out$wt)
    }
    n <- apply(resids, 2L, length) - colSums(is.na(resids))
    lsqr <- ls.out$qr
    p <- lsqr$rank

    ## calculate total sum sq and df

    if(ls.out$intercept) {
	if(is.matrix(lsqr$qt))
	    totss <- colSums(lsqr$qt[-1L, ]^2)
	else totss <- sum(lsqr$qt[-1L]^2)
	degfree <- p - 1
    } else {
	totss <- colSums(as.matrix(lsqr$qt^2))
	degfree <- p
    }

    ## calculate residual sum sq and regression sum sq

    resss <- colSums(resids^2, na.rm=TRUE)
    resse <- (resss/(n-p))^.5
    regss <- totss - resss
    rsquared <- regss/totss
    fstat <- (regss/degfree)/(resss/(n-p))
    pvalue <- pf(fstat, degfree, (n-p), lower.tail = FALSE)

    ## construct summary

    Ynames <- colnames(resids)
    summary <- cbind(format(round(resse, digits)),
		     format(round(rsquared, digits)),
		     format(round(fstat, digits)),
		     format(degfree),
		     format(n-p),
		     format(round(pvalue, digits)))
    dimnames(summary) <- list(Ynames,
			      c("Mean Sum Sq", "R Squared",
				"F-value", "Df 1", "Df 2", "Pr(>F)"))
    mat <- as.matrix(lsqr$qr[1L:p, 1L:p])
    mat[row(mat)>col(mat)] <- 0
    qrinv <- solve(mat)

    ## construct coef table

    m.y <- ncol(resids)
    coef.table <- as.list(1L:m.y)
    if(m.y==1) coef <- matrix(ls.out$coefficients, ncol=1)
    else coef <- ls.out$coefficients
    for(i in 1L:m.y) {
	covmat <- (resss[i]/(n[i]-p)) * (qrinv%*%t(qrinv))
	se <- diag(covmat)^.5
	coef.table[[i]] <- cbind(coef[, i], se, coef[, i]/se,
				 2*pt(abs(coef[, i]/se), n[i]-p,
                                      lower.tail = FALSE))
	dimnames(coef.table[[i]]) <-
	    list(colnames(lsqr$qr),
		 c("Estimate", "Std.Err", "t-value", "Pr(>|t|)"))

	##-- print results --

	if(print.it) {
	    if(m.y>1)
		cat("Response:", Ynames[i], "\n\n")
	    cat(paste("Residual Standard Error=",
                      format(round(resse[i], digits)), "\nR-Square=",
                      format(round(rsquared[i], digits)), "\nF-statistic (df=",
		      format(degfree), ", ", format(n[i]-p), ")=",
		      format(round(fstat[i], digits)), "\np-value=",
		      format(round(pvalue[i], digits)), "\n\n", sep=""))
	    print(round(coef.table[[i]], digits))
	    cat("\n\n")
	}
    }
    names(coef.table) <- Ynames

    invisible(list(summary = summary, coef.table = coef.table))
}

#############################


fixedLassoPoly =
  function(X, y, lambda, beta, active, inactive = FALSE) {

    XA = X[, active, drop=FALSE]
    XI = X[, !active, drop=FALSE]
    XAi = pinv(crossprod(XA))
    XAp = XAi %*% t(XA)
    Ir = t(XI) %*% t(XAp)  # matrix in the "irrepresentable" condition

    if(length(lambda)>1) {
       lambdaA= lambda[active]
       lambdaI = lambda[!active]
    } else {
       lambdaA = rep(lambda, sum(active))
       lambdaI = rep(lambda, sum(!active))
    }

    penalized = lambdaA != 0
    signA = sign(beta[active])
    active_subgrad = signA * lambdaA
    if (length(signA)>1) sign_diag = diag(signA)
    if (length(signA)==1) sign_diag = matrix(signA, 1, 1)
    
    if (inactive) { # should we include the inactive constraints?
      RA = diag(rep(1, nrow(XA))) - XA %*% XAp # RA is residual forming matrix of selected model
      
      A = rbind(
        t(XI) %*% RA,
        -t(XI) %*% RA,
        -(sign_diag %*% XAp)[penalized,] # no constraints for unpenalized
      )

      b = c(
        lambdaI - Ir %*% active_subgrad,
        lambdaI + Ir %*% active_subgrad,
        -(sign_diag %*% XAi %*% active_subgrad)[penalized])
    } else {
      A = -(sign_diag %*% XAp)[penalized,]  # no constraints for unpenalized
      b = -(sign_diag %*% XAi %*% active_subgrad)[penalized]
    }
    
    return(list(A=A, b=b))
  }

##############################

## Approximates inverse covariance matrix theta
## using coordinate descent 

debiasingMatrix = function(Xinfo,               # could be X or t(X) %*% X / n 
                                                # depending on is_wide
                           is_wide,
                           nsample, 
                           rows, 
                           verbose=FALSE, 
                           bound=NULL,          # starting value of bound
                           linesearch=TRUE,     # do a linesearch?
                           scaling_factor=1.5,  # multiplicative factor for linesearch
                           max_active=NULL,     # how big can active set get?
                           max_try=10,          # how many steps in linesearch?
                           warn_kkt=FALSE,      # warn if KKT does not seem to be satisfied?
                           max_iter=50,         # how many iterations for each optimization problem
                           kkt_stop=TRUE,       # stop based on KKT conditions?
                           parameter_stop=TRUE, # stop based on relative convergence of parameter?
                           objective_stop=TRUE, # stop based on relative decrease in objective?
                           kkt_tol=1.e-4,       # tolerance for the KKT conditions
                           parameter_tol=1.e-4, # tolerance for relative convergence of parameter
                           objective_tol=1.e-4  # tolerance for relative decrease in objective
                           ) {


  if (is.null(max_active)) {
     max_active = max(50, 0.3 * nsample)
  } 

  p = ncol(Xinfo);
  M = matrix(0, length(rows), p);

  if (is.null(bound)) {
      bound = (1/sqrt(nsample)) * qnorm(1-(0.1/(p^2)))
  }
 
  xperc = 0;
  xp = round(p/10);
  idx = 1;
  for (row in rows) {
    if ((idx %% xp)==0){
      xperc = xperc+10;
      if (verbose) {
        print(paste(xperc,"% done",sep="")); }
    }

    output = debiasingRow(Xinfo,               # could be X or t(X) %*% X / n 
                                               # depending on is_wide
                          is_wide,
                          row,
                          bound,
                          linesearch=linesearch,
                          scaling_factor=scaling_factor,
                          max_active=max_active,
			  max_try=max_try,
			  warn_kkt=FALSE,
			  max_iter=max_iter,
			  kkt_stop=kkt_stop,
			  parameter_stop=parameter_stop,
			  objective_stop=objective_stop,
			  kkt_tol=kkt_tol,
			  parameter_tol=parameter_tol,
			  objective_tol=objective_tol)

    if (warn_kkt && (!output$kkt_check)) {
       warning("Solution for row of M does not seem to be feasible")
    } 
  
    if (!is.null(output$soln)) {
        M[idx,] = output$soln;
    } else {
        stop(paste("Unable to approximate inverse row ", row));
    }

    idx = idx + 1;
  }
  return(M)
}

# Find one row of the debiasing matrix -- assuming X^TX/n is not too large -- i.e. X is tall

debiasingRow = function (Xinfo,               # could be X or t(X) %*% X / n depending on is_wide
                         is_wide, 
                         row, 
                         bound, 
	                 linesearch=TRUE,     # do a linesearch?
		         scaling_factor=1.5,  # multiplicative factor for linesearch
		         max_active=NULL,     # how big can active set get?
			 max_try=10,          # how many steps in linesearch?
			 warn_kkt=FALSE,      # warn if KKT does not seem to be satisfied?
			 max_iter=50,         # how many iterations for each optimization problem
                         kkt_stop=TRUE,       # stop based on KKT conditions?
                         parameter_stop=TRUE, # stop based on relative convergence of parameter?
                         objective_stop=TRUE, # stop based on relative decrease in objective?
                         kkt_tol=1.e-4,       # tolerance for the KKT conditions
			 parameter_tol=1.e-4, # tolerance for relative convergence of parameter
			 objective_tol=1.e-4  # tolerance for relative decrease in objective
                         ) {

  p = ncol(Xinfo)

  if (is.null(max_active)) {
      max_active = min(nrow(Xinfo), ncol(Xinfo))
  }

   
  # Initialize variables 

  soln = rep(0, p)
  soln = as.numeric(soln)
  ever_active = rep(0, p)
  ever_active[1] = row      # 1-based
  ever_active = as.integer(ever_active)
  nactive = as.integer(1)

  linear_func = rep(0, p)
  linear_func[row] = -1
  linear_func = as.numeric(linear_func)
  gradient = 1. * linear_func 

  counter_idx = 1;
  incr = 0;

  last_output = NULL

  if (is_wide) {
     Xsoln = as.numeric(rep(0, nrow(Xinfo)))
  }

  while (counter_idx < max_try) {

      if (!is_wide) {
          result = solve_QP(Xinfo, # this is non-neg-def matrix
                            bound, 
                            max_iter, 
                            soln, 
                            linear_func, 
                            gradient, 
                            ever_active, 
                            nactive, 
                            kkt_tol, 
                            objective_tol, 
			    parameter_tol,
                            max_active,
			    kkt_stop,
			    objective_stop,
			    parameter_stop)
      } else {
          result = solve_QP_wide(Xinfo,                      # this is a design matrix
                                 as.numeric(rep(bound, p)),  # vector of Lagrange multipliers
				 0,                          # ridge_term 
                                 max_iter, 
                                 soln, 
                                 linear_func, 
                                 gradient, 
                                 Xsoln,
                                 ever_active, 
                                 nactive, 
                                 kkt_tol, 
                                 objective_tol, 
				 parameter_tol,
                                 max_active,
				 kkt_stop,
				 objective_stop,	
				 parameter_stop)

      }

      iter = result$iter

      # Logic for whether we should continue the line search

      if (!linesearch) {
        break
      }

      if (counter_idx == 1){
        if (iter == (max_iter+1)){
           incr = 1; # was the original problem feasible? 1 if not
         } else {
           incr = 0; # original problem was feasible
         }
      } 

      if (incr == 1) { # trying to find a feasible point
         if ((iter < (max_iter+1)) && (counter_idx > 1)) { 
           break;      # we've found a feasible point and solved the problem            
         }
         bound = bound * scaling_factor;
      } else {         # trying to drop the bound parameter further
         if ((iter == (max_iter + 1)) && (counter_idx > 1)) {
            result = last_output; # problem seems infeasible because we didn't solve it
   	    break;                # so we revert to previously found solution
         }
         bound = bound / scaling_factor;
      }

      # If the active set has grown to a certain size
      # then we stop, presuming problem has become
      # infeasible.

      # We revert to the previous solution
	
      if (result$max_active_check) {
	  result = last_output;
	  break;
      }
      
      counter_idx = counter_idx + 1
      last_output = list(soln=result$soln,
                         kkt_check=result$kkt_check)
    }


  # Check feasibility

  if (warn_kkt && (!result$kkt_check)) {
     warning("Solution for row of M does not seem to be feasible")
  } 

  return(list(soln=result$soln,
              kkt_check=result$kkt_check,
	      gradient=result$gradient))

}


##############################

print.fixedLassoInf <- function(x, tailarea=TRUE, ...) {
  cat("\nCall:\n")
  dput(x$call)

  cat(sprintf("\nStandard deviation of noise (specified or estimated) sigma = %0.3f\n",
              x$sigma))

  cat(sprintf("\nTesting results at lambda = %0.3f, with alpha = %0.3f\n",x$lambda,x$alpha))
  cat("",fill=T)
  tab = cbind(x$vars,
    round(x$coef0,3),
    round(x$coef0 / x$sd,3),
    round(x$pv,3),round(x$ci,3))
  colnames(tab) = c("Var", "Coef", "Z-score", "P-value", "LowConfPt", "UpConfPt")
  if (tailarea) {
    tab = cbind(tab,round(x$tailarea,3))
    colnames(tab)[(ncol(tab)-1):ncol(tab)] = c("LowTailArea","UpTailArea")
  }
  rownames(tab) = rep("",nrow(tab))
  print(tab)

  cat(sprintf("\nNote: coefficients shown are %s regression coefficients\n",
              ifelse(x$type=="partial","partial","full")))
  invisible()
}

#estimateLambda <- function(x, sigma, nsamp=1000){
#  checkargs.xy(x,rep(0,nrow(x)))
#  if(nsamp < 10) stop("More Monte Carlo samples required for estimation")
#  if (length(sigma)!=1) stop("sigma should be a number > 0")
 # if (sigma<=0) stop("sigma should be a number > 0")

 # n = nrow(x)
 # eps = sigma*matrix(rnorm(nsamp*n),n,nsamp)
 # lambda = 2*mean(apply(t(x)%*%eps,2,max))
 # return(lambda)
#}

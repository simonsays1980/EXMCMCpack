//////////////////////////////////////////////////////////////////////////
// MCMCmnlMH.cc samples from the posterior distribution of a multinomial
// logit model using a random walk Metropolis algorithm.
//
// The initial version of this file was generated by the
// auto.Scythe.call() function in the MCMCpack R package
// written by:
//
// Andrew D. Martin
// Dept. of Political Science
// Washington University in St. Louis
// admartin@wustl.edu
//
// Kevin M. Quinn
// Dept. of Government
// Harvard University
// kevin_quinn@harvard.edu
// 
// This software is distributed under the terms of the GNU GENERAL
// PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
// file for more information.
//
// This file was initially generated on Wed Dec 29 15:27:08 2004
// 12/31/2004 filled out template and got it initial version working (KQ)
// 7/27/2007 DBP ported to scythe 1.0
//
// Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
// Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
//    and Jong Hee Park
//////////////////////////////////////////////////////////////////////////


#ifndef MCMCMNLMH_CC
#define MCMCMNLMH_CC

#include "include/scythestat/matrix.h"
#include "include/scythestat/distributions.h"
#include "include/scythestat/stat.h"
#include "include/scythestat/la.h"
#include "include/scythestat/ide.h"
#include "include/scythestat/smath.h"
#include "MCMCrng.h"
#include "MCMCfcds.h"
#include "MCMCmnl.h"

#include <R.h>           // needed to use Rprintf()
#include <R_ext/Utils.h> // needed to allow user interrupts

using namespace std;
using namespace scythe;

// natural log of the multivariate-t density (up to a constant of 
//   proportionality)
// theta: eval point
// mu:    mode
// C:     Cholesky factor of inverse scale matrix
// df:    degrees of freedom
static inline double lnmulttdens(const Matrix<>& theta, 
			  const Matrix<>& mu,
			  const Matrix<>& C,
			  const double& df){
  
  const int d = theta.size();
  //const Matrix<> z = t(theta - mu) * C; 
  // C is now C' if VC mat is C C'
  const Matrix<> z = C * (theta - mu);
  double zsumsq = 0;
  for (int i=0; i<d; ++i){
    zsumsq += std::pow(z[i], 2);
  }
  
  return ( (-(df + d)/2) * std::log(1 + (1/df) * zsumsq)  );
}



template <typename RNGTYPE>
void MCMCmnlMH_impl(rng<RNGTYPE>& stream, const Matrix<>& Y, 
		    const Matrix<>& X, const Matrix<>& b0,
		    const Matrix<>& B0, const Matrix<>& V,
		    Matrix<>& beta, const Matrix<>& beta_hat, 
		    const Matrix<>& tune,
		    const unsigned int burnin, const unsigned int mcmc,
		    const unsigned int thin, const unsigned int verbose,
		    const unsigned int RW, const double tdf,
		    Matrix<>& storemat)
{
  // define constants
  const unsigned int tot_iter = burnin + mcmc;  // total iterations
  const unsigned int nstore = mcmc / thin;      // # of draws to store
  const unsigned int k = X.cols();
  
  // Initialize storage matrix
  storemat = Matrix<>(nstore, k, false);
  
  // proposal parameters
  const Matrix<> propV = tune * V * tune;
  const Matrix<> propC = cholesky(propV);    
  const Matrix<> propCinvT = t(cholesky(invpd(propV)));

  
  double logpost_cur = mnl_logpost(Y, X, beta, b0, B0);
  double logjump_cur =  lnmulttdens(beta, beta_hat, propCinvT, tdf);

  int count = 0;
  int accepts = 0;

  ///// MCMC SAMPLING OCCURS IN THIS FOR LOOP
  for (unsigned int iter = 0; iter < tot_iter; ++iter) {
    
    // sample beta
    if (RW == 0){ // Independent Metropolis-Hastings
    
      const double u = stream();
 
      if (u < 0.75){
       	const Matrix<> beta_can = beta_hat + stream.rmvt(propV, tdf);
	const double logpost_can = mnl_logpost(Y, X, beta_can, b0, B0);
	const double logjump_can = lnmulttdens(beta_can, beta_hat, 
					       propCinvT, tdf);
	
	const double ratio = std::exp( logpost_can - logjump_can - 
				       logpost_cur + logjump_cur );
	
	if (stream() < ratio) {
	  beta = beta_can;
	  logpost_cur = logpost_can;
	  logjump_cur = logjump_can;
	  ++accepts;
	}
      }
      else{
       	const Matrix<> beta_can = beta_hat + beta_hat - beta;
	const double logpost_can = mnl_logpost(Y, X, beta_can, b0, B0);
	const double logjump_can = lnmulttdens(beta_can, beta_hat, 
					       propCinvT, tdf);
	
	const double ratio = std::exp( logpost_can  - 
				       logpost_cur );
	
	if (stream() < ratio) {
	  beta = beta_can;
	  logpost_cur = logpost_can;
	  logjump_cur = logjump_can;
	  ++accepts;
	}
	
      }
      

    }
    else{ // Random Walk Metropolis
      const 
	Matrix<> beta_can = gaxpy(propC, stream.rnorm(k,1,0,1), beta);      
      const double logpost_can = mnl_logpost(Y, X, beta_can, b0, B0);
      const double ratio = std::exp(logpost_can - logpost_cur); 
      
      if (stream() < ratio) {
	beta = beta_can;
	logpost_cur = logpost_can;
	++accepts;
      }
      
    }
    
    // store values in matrices
    if (iter >= burnin && ((iter % thin) == 0)) { 
      for (unsigned int j = 0; j < k; j++)
	storemat(count, j) = beta[j];
      
      ++count;
    }
    
    // print output to stdout
    if (verbose > 0 && iter % verbose == 0) {
      Rprintf("\n\nMCMCmnl Metropolis iteration %i of %i \n", 
	      (iter+1), tot_iter);
      Rprintf("beta = \n");
      
      for (unsigned int j=0; j<k; ++j)
	Rprintf("%10.5f\n", beta[j]);
      
      Rprintf("Metropolis acceptance rate for beta = %3.5f\n\n", 
	      static_cast<double>(accepts) / static_cast<double>(iter+1));	
    }
    
    R_CheckUserInterrupt(); // allow user interrupts       
  } // end MCMC loop
}


extern "C" {
  
  // MCMC sampling for multinomial logit via Metropolis-Hastings
  void MCMCmnlMH(double *sampledata, const int *samplerow, 
		 const int *samplecol, const double *Ydata, 
		 const int *Yrow, const int *Ycol, 
		 const double *Xdata, const int *Xrow, const int *Xcol, 
		 const int *burnin, const int *mcmc, const int *thin, 
		 const double *tunedata, const int *tunerow, 
		 const int *tunecol, const int *uselecuyer, 
		 const int *seedarray, const int *lecuyerstream, 
		 const int *verbose, const double *betastartdata,
		 const int *betastartrow, const int *betastartcol,
		 const double *betamodedata,
		 const int *betamoderow, const int *betamodecol,
		 const double *b0data, const int *b0row, 
		 const int *b0col, const double *B0data, 
		 const int *B0row, const int *B0col, 
		 const double *Vdata, const int *Vrow, 
		 const int *Vcol, const int* RW, const double* tdf) 
  {
    
    // pull together Matrix objects
    // REMEMBER TO ACCESS PASSED ints AND doubles PROPERLY
    const Matrix<> Y(*Yrow, *Ycol, Ydata);
    const Matrix<> X(*Xrow, *Xcol, Xdata);
    const Matrix<> tune(*tunerow, *tunecol, tunedata);
    Matrix<> beta(*betastartrow, *betastartcol, betastartdata);     
    const Matrix<> betamode(*betamoderow, *betamodecol, betamodedata);     
    const Matrix<> b0(*b0row, *b0col, b0data);
    const Matrix<> B0(*B0row, *B0col, B0data);
    const Matrix<> V(*Vrow, *Vcol, Vdata);
    
    // storage matrix or matrices
    Matrix<> storemat;
    MCMCPACK_PASSRNG2MODEL(MCMCmnlMH_impl, Y, X, b0, B0, V, beta, betamode,
			   tune, *burnin, *mcmc, *thin, *verbose, *RW,  
			   *tdf, storemat);
    
    // load draws into sample array
    for(unsigned int i = 0; i < storemat.size(); ++i)
      sampledata[i] = storemat(i);
    
  } // end MCMCmnlMH 
} // end extern "C"

#endif

/* Copyright (C) 2005-2006  Oscar Rueda Palacio and Ramón Díaz-Uriarte */

/* This program is free software; you can redistribute it and/or */
/* modify it under the terms of the GNU General Public License */
/* as published by the Free Software Foundation; either version 2 */
/* of the License, or (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with this program; if not, write to the Free Software */
/* Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, */
/* USA. */




#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Utils.h>

#include<limits.h>


/* Turn all C++ new to R_alloc, using sed as: */

/* sed 's/double \*\([a-zA-Z0-9_-]*\) = new double\[\([ *+()a-zA-Z0-9_-]*\)\]/double *\1; \1 = (double *) R_alloc(\2, sizeof(double))/' nnhl.proto.c > nnhl-a.c */
/* sed 's/int \*\([a-zA-Z0-9_-]*\) = new int\[\([ *+()a-zA-Z0-9_-]*\)\]/int *\1; \1 = (int *) R_alloc(\2, sizeof(int))/' nnhl-a.c > nnhl-b.c */

/* Later turned to Calloc (and corresponding Free) */

/* Even if now legal, just in case turn all C++ style comments to C with
sed 's/\/\/\(.*\)/\/\* \1 \*\//' nnhl.c > nnhl2.c
*/


/*  Modified from McConnells 'Complete Code' */
#define ASSERT(condition, message) {       \
  if( !(condition)) {                      \
    printf("\n !! Assertion failed :");    \
    printf( #condition );                  \
    printf( message);                      \
    abort();                               \
  }                                        \
}                                          \



/* This was from the C++ version */
/* #ifdef DEBUG */
/*  #define PR(x) std::cout << std::endl << #x " = " << x << std::endl; */
/* #else */
/*  #define PR(x); */
/* #endif */


#ifdef DEBUG
 #define PR(x) {printf("\n %s = %f \n", #x,  (float) x); fflush(stdout);}
#else
 #define PR(x);
#endif


#ifdef DEBUG
 #define PR2(x) {printf("\n %s = %p \n", #x,  x); fflush(stdout);}
#else
 #define PR2(x);
#endif


void dinvgamma(double *x, double *shape, double *rate, int give_log, 
	       double *res) {
  *res = *shape * log(*rate) - lgamma(*shape) - (*shape+1)*log(*x) -
    (*rate / *x);
  if (!give_log) *res = exp(*res);
}

void Free_2(double **x, int nrow) { /* Free 2-D array */
  int i; 
  for (i = 0; i < nrow; ++i) { 
    Free(x[i]);
  } 
  Free(x); 
}

void normalNHHMMlikelihood(double *y, int *k, double *x, 
			   int *n, double *q, double *beta, 
			   double *stat, double *mu, double *sigma2, 
			   double *loglik)
{

#ifdef DEBUG
      printf("\n Likelihood evaluation \n\n"); fflush(stdout);
#endif



  /* TOL to assure that denominator c does not get zero */
  double res=0;
  
  int num_states = *k;
  int hidden_state = 0;
  int hidden_state_b = 0;
/*   double (**Q) = new double *[num_states]; */
/*   for (int hidden_state = 0; hidden_state< num_states; ++hidden_state) */
/*     Q[hidden_state] = new double[num_states]; */


  double **Q;
  Q = Calloc(num_states, double*);
  
  for (hidden_state = 0; hidden_state< num_states; ++hidden_state)
    Q[hidden_state] = Calloc(num_states, double);


  
  if (*k==1) { /*  with one state is just log lik. */
    int ii = 0;
    for (ii = 0; ii<*n; ++ii) {
      res += dnorm(y[ii], *mu, sqrt(*sigma2), 1);
    }
  }
  else { /*  >1 states */
    double c_i, y_i, x_i, mu_j, sd_j, rowSumsQ;
    double *filter; filter = Calloc(*k, double); 
    double *filtercond; filtercond = Calloc(*k, double);
    
    for (hidden_state=0; hidden_state < *k; ++hidden_state) {
      /*  next is filtercond for i==0;  */
      /*  later we change it for i > 0 */
      filtercond[hidden_state] = stat[hidden_state]; 
    }
	
    /*  Loop over all positions */
    int obs_index = 0;
    for (obs_index = 0; obs_index < *n; ++obs_index) {
      c_i = 0;
      y_i = y[obs_index];
      x_i = x[obs_index];
	  

/* #ifdef DEBUG */
/*       printf("\n ni:  %d", *n); */
/*       if (obs_index == 0) {printf("\n x_0:  %f", x_i); fflush(stdout);} */
/*       if (obs_index == 1) {printf("\n x_1:  %f", x_i); fflush(stdout);} */
/*       if (obs_index == 2) {printf("\n x_2:  %f", x_i); fflush(stdout);} */
/*       if (obs_index == 3) {printf("\n x_3:  %f", x_i); fflush(stdout);} */
/*       if (obs_index == 4) {printf("\n x_4:  %f", x_i); fflush(stdout);} */
/*       if (obs_index == 5) {printf("\n x_5:  %f", x_i); fflush(stdout);} */
/*       if (obs_index == 6) {printf("\n x_6:  %f", x_i); fflush(stdout);} */
/*       if (obs_index == 7) {printf("\n x_7:  %f", x_i); fflush(stdout);} */
/*       if (obs_index == 8) {printf("\n x_8:  %f", x_i); fflush(stdout);} */
/*       if (obs_index == 9) {printf("\n x_9:  %f", x_i); fflush(stdout);} */
/*       if (obs_index == 10) {printf("\n x_10:  %f", x_i); fflush(stdout);} */
/*       if (obs_index == 11) {printf("\n x_11:  %f", x_i);  fflush(stdout);} */
/*       if (obs_index == 12) {printf("\n x_12:  %f", x_i); fflush(stdout);} */
/*       if (obs_index == 13) {printf("\n x_13:  %f", x_i); fflush(stdout);} */
/*       if (obs_index == 14) {printf("\n x_14:  %f", x_i); fflush(stdout);} */
/*       if (obs_index == 15) {printf("\n x_15:  %f", x_i); fflush(stdout);} */
/*       if (obs_index == 16) {printf("\n x_16:  %f", x_i); fflush(stdout);} */
/*       if (obs_index == 17) {printf("\n x_17:  %f", x_i); fflush(stdout);} */
/*       if (obs_index == 18) {printf("\n x_18:  %f", x_i); fflush(stdout);} */
/*       if (obs_index == 19) {printf("\n x_19:  %f", x_i); fflush(stdout);} */
/*       if (obs_index == 20) {printf("\n x_20:  %f", x_i); fflush(stdout);} */
/* #endif */

		
      for (hidden_state = 0; hidden_state < *k; ++hidden_state) {
	mu_j = mu[hidden_state];
	sd_j = sqrt(sigma2[hidden_state]);
	/*  will deal with this later zz: dnorm */
	filter[hidden_state] = filtercond[hidden_state] *
	  dnorm(y_i, mu_j, sd_j, 0);
	c_i += filter[hidden_state];
      }
      /*  Check c[i] to avoid division by zero */
      /*  The next should NEVER happen. And we don't want an "==" */
      /*  with a float. Check for errors a+nyway */
      /*  zz: recall to check for boundaies in other quantities */
      /*  (filter, filtercond, etc) */
      if(c_i <= 0) {
	c_i = 10E-100;
      }

      
      res += log(c_i);
      if (obs_index == ((*n) - 1)) break;
	
     
      /*  for debugging */
      /* double Q1; */
      /* double Q2; */
      /* double Q3; */
      /* double Q4; */
      for (hidden_state = 0; hidden_state <*k; ++hidden_state) {
	rowSumsQ = 0;
	for(hidden_state_b = 0; hidden_state_b < *k; ++hidden_state_b) {
	  /*  			Q1 = exp(q[hidden_state_b * *k + hidden_state]); */
	  /*  			Q2 = exp(beta[hidden_state_b * *k + hidden_state]); */
	  /*  			Q2 = pow(Q2, x_i); problematic */
	  /*  			Q3 = exp(x_i); problematic */
	  /*  			Q4 = pow(Q3, beta[hidden_state_b * *k + hidden_state]); */
			
	  /*  			Q[hidden_state][hidden_state_b] = Q1 * Q2; */
			
	  Q[hidden_state][hidden_state_b] = 
	    exp(q[hidden_state_b * *k + hidden_state] *
		(1- x_i));

	  rowSumsQ += Q[hidden_state][hidden_state_b];
	}
	for(hidden_state_b=0; hidden_state_b < *k; ++hidden_state_b) {
	  /* prevent overflow */ 
	  Q[hidden_state][hidden_state_b] = 
	    Q[hidden_state][hidden_state_b] / rowSumsQ;
	  /* printf("%f\n", Q[j][l]); */
	}
      }
		
      for (hidden_state = 0; hidden_state < *k; ++hidden_state) {
	filtercond[hidden_state] = 0;
	filter[hidden_state] = filter[hidden_state] / c_i;
      }
      for(hidden_state=0; hidden_state <*k; ++hidden_state) {
	/*  But we can actually pull the division outside and do at end: zzz */
	/*  filtercond[j] = filtercond[j] / c_i; */
	for(hidden_state_b = 0; hidden_state_b < *k; ++hidden_state_b) {
	  filtercond[hidden_state] += filter[hidden_state_b] * Q[hidden_state_b][hidden_state];
	}
      }
    } /*  loop over index_cond */
    /* 	  delete [] filter; */
/*     NULL_DELETE_ARRAY(filter); */
/*     delete [] filtercond; */
    Free(filter);
    Free(filtercond);
  
  } /*  else {...: when k > 1, so we have more than one hidden state */
    /*     printf("desde la funcion: res=%f\n", res); */
/*   delete_2(Q, num_states); */
  Free_2(Q, num_states);

  *loglik = res;
}



void Birth(double *y, int *varEqual, int *genome, int *index, 
	   double *mu, double *sigma2, 
	   double *beta, double *stat,
	   double *statBirth, 
	   int *r, double *loglikLast, double *probK, double *pb, double *muAlfa,
	   double *muBeta, double *x, int *n, double *candidatoMu,
	   double *candidatoSigma2,  
	   double *candidatoBeta, double *loglikBirth,
	   int *accepted, double *maxVar) {

    double *candidatoQ; candidatoQ = Calloc((*r+1)*(*r+1), double);
    double probBirth = 0;
    int newState;
    double loglikCandidate = 0;
    double loglikPartial = 0;
    int i, j, k;
    int reachedMaxVar = 0;
    int nn;
    for (i=0; i< *r; ++i) {
      candidatoMu[i] = mu[i];
      /*  we can try to adapt sigma2 */
      candidatoSigma2[i] = sigma2[i];
    }
    candidatoMu[*r] = rnorm(*muAlfa, *muBeta);
    if (*varEqual) {
      candidatoSigma2[*r] = candidatoSigma2[0];
    }
    else {
      candidatoSigma2[*r] = runif(0, *maxVar);
      candidatoSigma2[*r] = pow(candidatoSigma2[*r], 2);
    }

    /*  Check that the new variance does not exceed the maximum variance     */

    if (!reachedMaxVar) {

    for (i=0; i<*r; ++i) {
      for (j=0; j<*r; ++j) {
	candidatoBeta[i*(*r+1)+j] = beta[i* *r +j];
	candidatoQ[i*(*r+1)+j] = -candidatoBeta[i*(*r+1)+j];
      }
      candidatoBeta[i* (*r+1) + *r] = rgamma(1,1);
      candidatoQ[i* (*r+1) + *r]  = - candidatoBeta[i * (*r+1) + *r];
    }
    for (i=0; i<*r; ++i) {
      candidatoBeta[*r * (*r+1) + i] = rgamma(1,1);
      candidatoQ[*r * (*r+1)  + i] = - candidatoBeta[*r * (*r+1) + i];
    }
    candidatoBeta[(*r+1) * (*r+1) -1] = 0;
    candidatoQ[(*r+1) * (*r+1) -1] = 0;

    newState = *r + 1;

    probBirth = log(probK[*r]) - log(probK[*r-1]);
    for (k=0; k < *genome; ++k) {
      nn = index[k+1] - index[k];
      
      double *yy; yy = Calloc(nn, double);
      double *xx; xx = Calloc(nn, double);
      
     
      for (j=0; j < nn-1; ++j) {
	yy[j] = y[j + index[k]];
	xx[j] = x[j + index[k]];
      }
      xx[nn-1] = 0;
      yy[nn-1] = y[nn-1 + index[k]];
      normalNHHMMlikelihood(yy, &newState, xx, &nn, candidatoQ, candidatoBeta, 
			  statBirth, candidatoMu, candidatoSigma2, 
			  &loglikPartial);
      loglikCandidate += loglikPartial;
/*       delete [] xx; */
/*       delete [] yy; */
      Free(xx);
      Free(yy);

    }

    probBirth = probBirth + loglikCandidate - *loglikLast;
    probBirth = exp(probBirth) * (1-pb[*r]) / pb[*r-1];
    if (probBirth >1) probBirth = 1;
    if (runif(0,1) <= probBirth) {
      *accepted =1;
      *loglikBirth = loglikCandidate;
          }
    }
/*    delete [] candidatoQ; */
    Free(candidatoQ);
  }

void Death(double *y, int *genome, int *index, double *mu, double *sigma2, 
	   double *beta, double *stat,
	   double *statDeath, 
	   int *r, double *loglikLast, double *probK, double *pb, 
	   double *x, int *n, double *candidatoMu,
	   double *candidatoSigma2, 
	   double *candidatoBeta, double *loglikDeath,
	   int *accepted) {

  double *candidatoQ; candidatoQ = Calloc((*r-1)*(*r-1), double);
  double probDeath = 0;
  int death, indexBeta=0, indexMu=0;
  int newState;
  double loglikCandidate = 0;
  double loglikPartial = 0;
  int nn;
  int i, j, k;
  death = (int)rint(runif(1, *r));
  if (death >1) {
    for (i=0; i< (death-1); ++i) {
      candidatoMu[indexMu] = mu[i];
      candidatoSigma2[indexMu] = sigma2[i];
      indexMu ++;
      for (j=0; j< (death-1); ++j) {
	candidatoBeta[indexBeta] = beta[i* *r + j];
	candidatoQ[indexBeta] = -candidatoBeta[indexBeta];
	indexBeta ++;

      }
      if (death < *r) {
	for (j=death; j< *r; ++j) {
	  candidatoBeta[indexBeta] = beta[i* *r + j];
	  candidatoQ[indexBeta] = -candidatoBeta[indexBeta];
	  indexBeta ++;

	}
      }
    }
  }
  if (death < *r) {
    for (i=death; i< *r; ++i) {
      candidatoMu[indexMu] = mu[i];
      candidatoSigma2[indexMu] = sigma2[i];
      indexMu ++;

      if (death > 1) {
	for (j=0; j< (death-1); ++j) {
	  candidatoBeta[indexBeta] = beta[i* *r + j];
	  candidatoQ[indexBeta] = -candidatoBeta[indexBeta];
	  indexBeta ++;
	}
      }
      for (j=death; j< *r; ++j) {
	candidatoBeta[indexBeta] = beta[i* *r + j];
	candidatoQ[indexBeta] = -candidatoBeta[indexBeta];
	indexBeta++;
      }
    }
  }
  newState = *r -1;
  probDeath = log(probK[*r-1]) - log(probK[*r-2]);
  for (k=0; k < *genome; ++k) {
    nn = index[k+1] - index[k];
    double *yy; yy = Calloc(nn, double);
    double *xx; xx = Calloc(nn, double);
    for (j=0; j < nn-1; ++j) {
      yy[j] = y[j + index[k]];
      xx[j] = x[j + index[k]];
    }
    xx[nn-1] = 0;
    yy[nn-1] = y[nn-1 + index[k]];
    normalNHHMMlikelihood(yy, &newState, xx, &nn, candidatoQ, candidatoBeta, 
			  statDeath, candidatoMu, candidatoSigma2, 
			  &loglikPartial);
    loglikCandidate += loglikPartial;
/*     delete [] xx; */
/*     delete [] yy; */

    Free(xx);
    Free(yy);
  }
  probDeath = probDeath - loglikCandidate + *loglikLast;
  probDeath = exp(probDeath) * (1-pb[*r-1]) / pb[*r-2];
  probDeath = 1 / probDeath;
  if (probDeath >1) probDeath = 1;
  if (runif(0,1) <= probDeath) {
    *accepted =1;
    *loglikDeath = loglikCandidate;
  }
/*   delete [] candidatoQ; */
  Free(candidatoQ);

}
/*  this model splits somewhat different the beta matrix */
void Split(double *y, int *varEqual, int *genome, int *index, 
	   double *mu, double *sigma2, double *beta, 
	   double *stat, double *statSplit, 
	   int *r, double *loglikLast, double *probK, double *ps, 
	   double *x, int *n, double *candidatoMu,
	   double *candidatoSigma2, 
	   double *candidatoBeta, double *loglikSplit,
	   double *muAlfa, double *muBeta,
	   double *tauSplitMu,
	   double *tauSplitBeta, int *accepted) {

#ifdef DEBUG
  printf("\n Entering split \n"); fflush(stdout);
#endif

  double *ui; ui = Calloc(*r-1, double);
  double *epj; epj = Calloc(*r-1, double);
  double *gi0; gi0 = Calloc(2, double);
  double *candidatoQ; candidatoQ = Calloc((*r+1)*(*r+1), double);

  double priorSigma2, priorCandidatoSigma2;
  double probSplit = 0;
  int split, indexBeta=0, indexMu=0, indexEpj=0, indexUi=0;
  int newState;
  double epMu, epSigma2;
  double jacobian;
  int i, j, k;
  double loglikCandidate = 0;
  double loglikPartial = 0;
  int nn;
  /*  auxiliary variables */
  epMu = rnorm(0, *tauSplitMu);
  /*  In this version epSigma2 is uniform, to reduce the variance */
  epSigma2 = rbeta(2,2);
  for(i=0; i < (*r-1); ++i) {
    ui[i] =  rbeta(2,2);
    epj[i] =  rlnorm(0, *tauSplitBeta);
  }
  gi0[0] = rgamma(1,1);
  gi0[1] = rgamma(1,1);
  
  split = (int)rint(runif(1, *r));
  if (split >1) {
    for (i=0; i< (split-1); ++i) {
      candidatoMu[indexMu] = mu[i];
      candidatoSigma2[indexMu++] = sigma2[i];
      for (j=0; j< (split-1); ++j) {
	candidatoBeta[indexBeta] = beta[i* *r + j];
	candidatoQ[indexBeta] = -candidatoBeta[indexBeta];
	indexBeta++;
      }
      candidatoBeta[indexBeta] = beta[i* *r + split-1] * ui[indexUi];
      candidatoQ[indexBeta] = - candidatoBeta[indexBeta]; indexBeta++;
      candidatoBeta[indexBeta] = beta[i* *r + split-1]  * (1-ui[indexUi++]);
      candidatoQ[indexBeta] = - candidatoBeta[indexBeta]; indexBeta++;
      if (split < *r) {
	for (j=split; j< *r; ++j) {
	  candidatoBeta[indexBeta] = beta[i* *r + j];
	  candidatoQ[indexBeta] = -candidatoBeta[indexBeta]; indexBeta++;
	}
      }
    }
  }
  candidatoMu[indexMu] = mu[split-1] - epMu;
  if (*varEqual) {
    candidatoSigma2[indexMu++] = sigma2[split-1];
  }
  else {
    candidatoSigma2[indexMu++] = sigma2[split-1] * epSigma2;
  }
  candidatoMu[indexMu] = mu[split-1] + epMu;
  if (*varEqual) {
    candidatoSigma2[indexMu++] = sigma2[split-1];
  }
  else {
    candidatoSigma2[indexMu++] = sigma2[split-1] * (1 - epSigma2);
  }

  for (j=0; j< split -1; ++j) {
    candidatoBeta[indexBeta] = beta[(split-1) * *r + j] * epj[indexEpj++];
    candidatoQ[indexBeta] = - candidatoBeta[indexBeta]; indexBeta++;
  }
  candidatoBeta[indexBeta] = 0;
  candidatoQ[indexBeta++] = 0;
  candidatoBeta[indexBeta] = gi0[0];
  candidatoQ[indexBeta] = - candidatoBeta[indexBeta]; indexBeta++;
  if (split < *r) {
    for (j=split; j< *r; ++j) {
      candidatoBeta[indexBeta] = beta[(split-1) * *r + j] * epj[indexEpj++];
      candidatoQ[indexBeta] = -candidatoBeta[indexBeta]; indexBeta++;
    }
  }
  indexUi = 0;
  indexEpj =0;
  for (j=0; j< split -1; ++j) {
    candidatoBeta[indexBeta] = beta[(split-1) * *r + j] / epj[indexEpj++];
    candidatoQ[indexBeta] = - candidatoBeta[indexBeta]; indexBeta++;
  }
  candidatoBeta[indexBeta] = gi0[1];
  candidatoQ[indexBeta] = - candidatoBeta[indexBeta]; indexBeta++;
  candidatoBeta[indexBeta] = 0;
  candidatoQ[indexBeta++] = 0;
  if (split < *r) {
    for (j=split; j< *r; ++j) {
      candidatoBeta[indexBeta] = beta[(split-1) * *r + j] / epj[indexEpj++];
      candidatoQ[indexBeta] = -candidatoBeta[indexBeta]; indexBeta++;
    }
  }
  if (split < *r) {
    for (i=split; i< *r; ++i) {
      candidatoMu[indexMu] = mu[i];
      candidatoSigma2[indexMu++] = sigma2[i];

      if (split > 1) {
	for (j=0; j< (split-1); ++j) {
	  candidatoBeta[indexBeta] = beta[i* *r + j];
	  candidatoQ[indexBeta] = -candidatoBeta[indexBeta]; indexBeta++;
	}
      }
      candidatoBeta[indexBeta] = beta[i* *r + split-1] * ui[indexUi];
      candidatoQ[indexBeta] = - candidatoBeta[indexBeta]; indexBeta++;
      candidatoBeta[indexBeta] = beta[i* *r + split-1]  * (1-ui[indexUi++]);
      candidatoQ[indexBeta] = - candidatoBeta[indexBeta]; indexBeta++;

      for (j=split; j< *r; ++j) {
	candidatoBeta[indexBeta] = beta[i* *r + j];
	candidatoQ[indexBeta] = -candidatoBeta[indexBeta]; indexBeta++;
      }
    }
  }

  ASSERT(indexBeta == ((*r + 1) * (*r + 1)), 
	 "Size of indexBeta not equal to max limit in Split" );


  /*  Check adjacency condition (Richardson and Green 1997) */

  double distSplit= fabs(candidatoMu[split] - candidatoMu[split-1]);
  int valid = 0;
  for (i=0; i<*r+1; ++i) {
    if ((i != split-1) && (i != split)) {
      if ((distSplit < fabs(candidatoMu[i] - candidatoMu[split-1])) |
	  (distSplit < fabs(candidatoMu[i] - candidatoMu[split])))
	valid = 1;
    }
  }
    
  /*  if r=1 it always does */
  if (*r==1) valid = 1;
  if (valid) {
    newState = *r + 1;
    probSplit = log(probK[*r]) - log(probK[*r-1]);
    /*  sum of priors */
    for (i=0; i<(*r+1); ++i) {
      probSplit = probSplit + dnorm(candidatoMu[i], *muAlfa, *muBeta, 1);
      probSplit = probSplit + log(1 / sqrt(candidatoSigma2[i]));
    }
    for (i=0; i<*r; ++i) {
      probSplit = probSplit - dnorm(mu[i], *muAlfa, *muBeta, 1);
      probSplit = probSplit - log(1 / sqrt(sigma2[i]));
    }
    /*  take out zeros of the diagonal */
    for (i=0; i<((*r+1) * (*r+1)); ++i) {
      if (candidatoBeta[i] > 0)
	probSplit = probSplit + dgamma(candidatoBeta[i], 1, 1, 1);
    }
    for (i=0; i<(*r * *r); ++i) {
      if (beta[i] > 0)
	probSplit = probSplit - dgamma(beta[i], 1, 1, 1);
    }
    for (k=0; k < *genome; ++k) {
      nn = index[k+1] - index[k];
      double *yy; yy = Calloc(nn, double);
      double *xx; xx = Calloc(nn, double);
      for (j=0; j < nn-1; ++j) {
	yy[j] = y[j + index[k]];
	xx[j] = x[j + index[k]];
      }
      xx[nn-1] = 0;
      yy[nn-1] = y[nn-1 + index[k]];

      normalNHHMMlikelihood(yy, &newState, xx, &nn, candidatoQ, candidatoBeta, 
			    statSplit, candidatoMu, candidatoSigma2, 
			    &loglikPartial);
#ifdef DEBUG
      printf("\n       After likelihood\n"); fflush(stdout);
#endif

      loglikCandidate += loglikPartial;
/*       delete [] xx; */
/*       delete [] yy; */
      Free(xx);
      Free(yy);

    }
    probSplit = probSplit + loglikCandidate - *loglikLast;
    probSplit = probSplit + log(*r + 1);
    probSplit = probSplit - log(2);
    probSplit = probSplit - dnorm(epMu, 0, *tauSplitMu, 1);
    for (i=0; i<(*r-1); ++i) {
      probSplit = probSplit - dlnorm(epj[i], 0, *tauSplitBeta, 1);
      probSplit = probSplit - dbeta(ui[i], 2, 2, 1);
    }
    probSplit = probSplit - dbeta(epSigma2, 2, 2, 1);
    probSplit = probSplit - dgamma(gi0[0], 1, 1, 1);
    probSplit = probSplit - dgamma(gi0[1], 1, 1, 1);

    /*  jacobian of the transformation */
    jacobian =*r * log(2) + 2 * log(sqrt(sigma2[split-1]));
    if (split > 1) {
      for (i=0; i< split-1; ++i) {
	jacobian = jacobian + log(beta[i* *r + split-1]) - log(epj[i]);
	jacobian = jacobian + log(beta[(split-1) * *r + i]);

      }
    }

    if (split < *r) {
      for (i=split; i< *r; ++i) {
	jacobian = jacobian + log(beta[i* *r + split-1]) - log(epj[i-1]);
	jacobian = jacobian + log(beta[(split-1) * *r + i]);
      }
    }

    probSplit = exp(probSplit + jacobian) * ((1-ps[*r]) / ps[*r-1]) * (double)(*r)/(double)(*r+1);
    if (probSplit >1) probSplit = 1;
    if (runif(0,1) <= probSplit) {
      *accepted =1;
      *loglikSplit = loglikCandidate;

    }
  }
 
  Free(ui);
  Free(epj);
  Free(gi0);

#ifdef DEBUG
  int dd1 = 0;
  for  (dd1 = 0; dd1 < ((*r + 1) * (*r + 1)); ++dd1) {
      PR(dd1);
      PR(candidatoQ[dd1]);
      PR(candidatoBeta[dd1]);
  }
#endif

  Free(candidatoQ);

#ifdef DEBUG
  printf("\n Exiting split \n"); fflush(stdout);
#endif
}

void Combine(double *y, int *varEqual, int *genome, int *index, 
	     double *mu, double *sigma2, 
	     double *beta, double *stat, double *statCombine, 
	     int *r, double *loglikLast, double *probK, double *ps, 
	     double *x, int *n, double *candidatoMu,
	     double *candidatoSigma2, 
	     double *candidatoBeta, double *loglikCombine,
	     double *muAlfa, double *muBeta,
	     double *tauSplitMu,
	     double *tauSplitBeta, int *accepted, double *maxVar) {

#ifdef DEBUG
  printf("\n Entering combine \n"); fflush(stdout);
#endif

  
  double *candidatoQ; candidatoQ = Calloc((*r-1)*(*r-1), double);
  double probCombine = 0;
  double priorSigma2, priorCandidatoSigma2;
  int combine, indexBeta=0, indexMu=0, indexEpj=0, indexUi=0;
  int newState;
  double *epj; epj = Calloc(*r-1, double);
  double *ui; ui = Calloc(*r-1, double);
  double epMu;
  double epSigma2;
  double *gi0; gi0 = Calloc(2, double);
  double jacobian;
  int i, j, k;
  double loglikCandidate = 0;
  double loglikPartial = 0;
  int reachedMaxVar = 0;
  int nn;

  combine = (int)rint(runif(1, *r));

  /*  choose the most closest state */
  /*  and keep combine, combine +1 */
  if (combine == *r) combine = combine - 1;
  if (combine > 1) {
    if (fabs(mu[combine-1] - mu[combine-2]) < fabs(mu[combine-1] - mu[combine]))
      combine = combine -1;
  }
    
  /*   printf("combine=%d , %d\n", combine, combine+1); */
  if (combine >1) {
    for (i=0; i< (combine-1); ++i) {
      candidatoMu[indexMu] = mu[i];
      candidatoSigma2[indexMu++] = sigma2[i];
      for (j=0; j< (combine-1); ++j) {
	candidatoBeta[indexBeta] = beta[i* *r + j];
	candidatoQ[indexBeta] = -candidatoBeta[indexBeta]; indexBeta++;
      }
      candidatoBeta[indexBeta] = beta[i* *r + combine-1] + beta[i* *r + combine];
      ui[indexUi++] = beta[i* *r + combine-1] / (beta[i* *r + combine -1] + 
						 beta[i* *r + combine]);
      candidatoQ[indexBeta] = - candidatoBeta[indexBeta]; indexBeta++;
      if (combine < *r) {
	for (j=combine +1; j< *r; ++j) {
	  candidatoBeta[indexBeta] = beta[i* *r + j];
	  candidatoQ[indexBeta] = -candidatoBeta[indexBeta]; indexBeta++;
	}
      }
    }
  }
  candidatoMu[indexMu] = (mu[combine-1] + mu[combine]) / 2;
  if (*varEqual) {
    candidatoSigma2[indexMu++] = sigma2[combine-1];
  }
  else {
    candidatoSigma2[indexMu++] = sigma2[combine-1] + sigma2[combine];
    if (candidatoSigma2[indexMu -1] > *maxVar) reachedMaxVar = 1;
  }
  /*  check the new variance does not reached maximum */
  if (!reachedMaxVar) {
    for (j=0; j< (combine -1); ++j) {
      candidatoBeta[indexBeta] = sqrt(beta[(combine-1) * *r + j] * beta[combine * *r +j]);
      epj[indexEpj++] = sqrt(beta[(combine-1) * *r + j] / beta[combine * *r +j]);
      candidatoQ[indexBeta] = - candidatoBeta[indexBeta]; indexBeta++;
    }
    candidatoBeta[indexBeta] = 0;
    candidatoQ[indexBeta++] = 0;
    if (combine < *r) {
      for (j=combine +1; j< *r; ++j) {
	candidatoBeta[indexBeta] = sqrt(beta[(combine-1) * *r + j] * beta[combine * *r +j]);
	epj[indexEpj++] = sqrt(beta[(combine-1) * *r + j] / beta[combine * *r +j]);
	candidatoQ[indexBeta] = -candidatoBeta[indexBeta]; indexBeta++;
      }
    }
    if (combine < (*r-1)) {
      for (i=combine+1; i< *r; ++i) {
	candidatoMu[indexMu] = mu[i];
	candidatoSigma2[indexMu++] = sigma2[i];
	  
	if (combine > 1) {
	  for (j=0; j< (combine-1); ++j) {
	    candidatoBeta[indexBeta] = beta[i* *r + j];
	    candidatoQ[indexBeta] = -candidatoBeta[indexBeta]; indexBeta++;
	  }
	}
	candidatoBeta[indexBeta] = beta[i* *r + combine-1] + beta[i* *r + combine];
	ui[indexUi++] = beta[i* *r + combine-1] / (beta[i* *r + combine -1] +
						   beta[i* *r + combine]);
	candidatoQ[indexBeta] = - candidatoBeta[indexBeta]; indexBeta++;

	for (j=combine+1; j< *r; ++j) {
	  candidatoBeta[indexBeta] = beta[i* *r + j];
	  candidatoQ[indexBeta] = -candidatoBeta[indexBeta]; indexBeta++;
	}
      }
    }

    ASSERT(indexBeta == ((*r - 1) * (*r - 1)), 
	   "Size of indexBeta not equal to max limit in Combine" );
      
    /*  recover the other auxiliary variables */

    epMu = (mu[combine] - mu[combine-1]) / (2*sqrt(sigma2[combine -1]));
    epSigma2 = sigma2[combine-1] / (sigma2[combine-1] + sigma2[combine]);
    gi0[0] = beta[*r * (combine-1) + combine];
    gi0[1] = beta[*r * combine + combine -1];
    newState = *r - 1;
    probCombine = log(probK[*r-1]) - log(probK[*r-2]);
    /*  sum of priors */
    for (i=0; i<*r; ++i) {
      probCombine = probCombine + dnorm(mu[i], *muAlfa, *muBeta, 1);
      probCombine = probCombine + log(1/sqrt(sigma2[i]));
    }
    for (i=0; i<(*r-1); ++i) {
      probCombine = probCombine - dnorm(candidatoMu[i], *muAlfa, *muBeta, 1);
      probCombine = probCombine - log(1/sqrt(candidatoSigma2[i]));
    }
    /*  take out zeros of the diagonal */
    for (i=0; i<(*r * *r); ++i) {
      if (beta[i] > 0)
	probCombine = probCombine + dgamma(beta[i], 1, 1, 1);
    }
      
    for (i=0; i<((*r-1) * (*r-1)); ++i) {
      if (candidatoBeta[i] > 0)
	probCombine = probCombine - dgamma(candidatoBeta[i], 1, 1, 1);
    }

    for (k=0; k < *genome; ++k) {
      nn = index[k+1] - index[k];
      double *yy; yy = Calloc(nn, double);
      double *xx; xx = Calloc(nn, double);
      for (j=0; j < nn-1; ++j) {
	yy[j] = y[j + index[k]];
	xx[j] = x[j + index[k]];
      }
      xx[nn-1] = 0;
      yy[nn-1] = y[nn-1 + index[k]];
      normalNHHMMlikelihood(yy, &newState, xx, &nn, candidatoQ, candidatoBeta, 
			    statCombine, candidatoMu, candidatoSigma2, 
			    &loglikPartial);
      loglikCandidate += loglikPartial;
/*       delete [] xx; */
/*       delete [] yy; */
      Free(xx);
      Free(yy);
    }
    probCombine = probCombine - loglikCandidate + *loglikLast;
    /*     printf("candida=%f last=%f\n", loglikCandidate, *loglikLast); */
    probCombine = probCombine + log(*r);
    probCombine = probCombine - log(2);
    probCombine = probCombine - dnorm(epMu, 0, *tauSplitMu, 1);
    probCombine = probCombine - dbeta(epSigma2, 2, 2, 1);
    if (*r >2) {
      for (i=0; i<*r-2; ++i) {
	probCombine = probCombine - dlnorm(epj[i], 0, *tauSplitBeta, 1);
	probCombine = probCombine - dbeta(ui[i], 2, 2, 1);
      }
    }
    probCombine = probCombine - dgamma(gi0[0], 1, 1, 1);
    probCombine = probCombine - dgamma(gi0[1], 1, 1, 1);
      
    /*  jacobian of the transformation */
    /*  r-1 because of a 2 in the denominator and the split equivalent is r-1 */
    jacobian = (*r-1) * log(2) + 2* log(sqrt(sigma2[combine-1]));
    if (combine > 1) {
      for (i=0; i< combine-1; ++i) {
	jacobian = jacobian + log(beta[i* *r + combine-1])- log(epj[i]);
	jacobian = jacobian + log(beta[(combine-1) * *r + i]) ;
      }
    }
    if (combine < *r-1) {
      for (i=combine; i< *r-1; ++i) {
	jacobian = jacobian + log(beta[(i+1)* *r + combine-1])- log(epj[i-1]);
	jacobian = jacobian + log(beta[(combine-1) * *r + i+1]);
      }
    }
      
    probCombine = exp(probCombine + jacobian) * ((1-ps[*r-1]) / ps[*r-2]) * 
      (double)(*r-1) / (double)(*r);
    probCombine = 1 /probCombine;
    if (probCombine >1) probCombine = 1;
    if (runif(0,1) <= probCombine) {
      *accepted =1;
	
      *loglikCombine = loglikCandidate;
    }
  }

/*   delete [] epj; */
/*   delete [] ui; */
/*   delete [] gi0; */
/*   delete [] candidatoQ; */
  Free(epj);
  Free(ui);
  Free(gi0);
  Free(candidatoQ);

#ifdef DEBUG
  printf("\n Exiting combine \n"); fflush(stdout);
#endif

}

void MetropolisUpdate(double *y, double *x, int *varEqual, 
		      int *genome, int *index, double *mu, 
		      double *sigma2, double *beta, double *stat, int *r, int *n, 
		      double *muAlfa, double *muBeta, 
		      double *sigmaTauMu,
		      double *sigmaTauSigma2, double *sigmaTauBeta,
		      double *loglikLast, double *maxVar) {
#ifdef DEBUG
  printf("\n Entering MetropolisUpdate \n"); fflush(stdout);
#endif


  int i, j, k;
  double acepProb;
  double loglikCandidate = 0;
  double loglikPartial = 0;
  double *q; q = Calloc(*r * *r, double);
  double priorSigma2, priorCandidatoSigma2;
  double *candidatoMu; candidatoMu = Calloc(*r, double);
  double *candidatoSigma2; candidatoSigma2 = Calloc(*r, double);
  double *candidatoBeta; candidatoBeta = Calloc(*r * *r, double);
  int reachedMaxVar=0;
  int nn;
  /* update mu */
  
  acepProb = 0;
  for(i=0; i<*r; ++i) {
    candidatoMu[i] = mu[i] + rnorm(0, *sigmaTauMu);
    acepProb = acepProb + dnorm(candidatoMu[i], *muAlfa, *muBeta, 1) -
      dnorm(mu[i], *muAlfa, *muBeta, 1);
  }
  
  for (i=0; i<(*r * *r); ++i) {
    q[i] = -beta[i];
  }

  for (k=0; k < *genome; ++k) {
    nn = index[k+1] - index[k];
    double *yy; yy = Calloc(nn, double);
    double *xx; xx = Calloc(nn, double);
    for (j=0; j < nn-1; ++j) {
      yy[j] = y[j + index[k]];
      xx[j] = x[j + index[k]];
    }
    xx[nn-1] = 0;
    yy[nn-1] = y[nn-1 + index[k]];
    normalNHHMMlikelihood(yy, r, xx, &nn, q, beta, stat, candidatoMu, sigma2, 
			  &loglikPartial);
    loglikCandidate += loglikPartial;
/*     delete [] xx; */
/*     delete [] yy; */
    Free(xx);
    Free(yy);
  }
      
  acepProb = acepProb + loglikCandidate - *loglikLast;
  acepProb = exp(acepProb);
  if (acepProb > 1) acepProb = 1;
  /*     printf("acepProb=%f\n", acepProb); */
  if (runif(0,1) < acepProb) {
    for (i=0; i<*r; ++i) mu[i] = candidatoMu[i];
    *loglikLast = loglikCandidate;
  }
  
  /*  Update sigma2 */
  loglikCandidate = 0;
  acepProb = 0;
  candidatoSigma2[0] = exp(log(sigma2[0]) + rnorm(0, *sigmaTauSigma2));
  if (candidatoSigma2[0] > *maxVar) reachedMaxVar = 1;
  for (i=1; i<*r; ++i) {
    if (*varEqual) {
      candidatoSigma2[i] = candidatoSigma2[0];
    }
    else {
      candidatoSigma2[i] = exp(log(sigma2[i]) + rnorm(0, *sigmaTauSigma2));
      if (candidatoSigma2[i] > *maxVar) reachedMaxVar = 1;
    }
  }
  
  for (i=0; i<(*r * *r); ++i) {
    q[i] = -beta[i];
  }
  /*  Check candidates doesn't reach Maximum variance */
  if (!reachedMaxVar) {
    for(i=0; i<*r; ++i) {
      acepProb = acepProb + log(1/sqrt(candidatoSigma2[i])) - log(1/sqrt(sigma2[i]));
    }
    for (k=0; k < *genome; ++k) {
      nn = index[k+1] - index[k];
      double *yy; yy = Calloc(nn, double);
      double *xx; xx = Calloc(nn, double);
      for (j=0; j < nn-1; ++j) {
	yy[j] = y[j + index[k]];
	xx[j] = x[j + index[k]];
      }
      xx[nn-1] = 0;
      yy[nn-1] = y[nn-1 + index[k]];      
      normalNHHMMlikelihood(yy, r, xx, &nn, q, beta, stat, mu, candidatoSigma2, 
			    &loglikPartial);
      loglikCandidate += loglikPartial;
/*       delete [] xx; */
/*       delete [] yy; */
      Free(xx);
      Free(yy);
    }
    acepProb = acepProb + loglikCandidate - *loglikLast;
    for (i=0; i<*r; ++i) {
      acepProb = acepProb + log(candidatoSigma2[i]) - log(sigma2[i]);
    }
    acepProb = exp(acepProb);
    if (acepProb > 1) acepProb = 1;
    if (runif(0,1) < acepProb) {
      for (i=0; i<*r; ++i) sigma2[i] = candidatoSigma2[i];
      *loglikLast = loglikCandidate;
    }
  }

  /* Update beta */
  loglikCandidate = 0; 	
  acepProb =0;
  if (*r>1) {
    for (i=0; i < (*r * *r); ++i) {
      if (beta[i] > 0) {
	/* log(beta[i] gets -Inf */
	candidatoBeta[i] = exp(log(beta[i]) + rnorm(0, *sigmaTauBeta));
	q[i] = -candidatoBeta[i];
	acepProb = acepProb + dgamma(candidatoBeta[i], 1, 1, 1) -
	  dgamma(beta[i], 1, 1, 1);
      }
      else {
	candidatoBeta[i] = beta[i];
	q[i] = -candidatoBeta[i];
      }
    }
    for (k=0; k < *genome; ++k) {
      nn = index[k+1] - index[k];
      double *yy; yy = Calloc(nn, double);
      double *xx; xx = Calloc(nn, double);
      for (j=0; j < nn-1; ++j) {
	yy[j] = y[j + index[k]];
	xx[j] = x[j + index[k]];
      }
      xx[nn-1] = 0;
      yy[nn-1] = y[nn-1 + index[k]];
      normalNHHMMlikelihood(yy, r, xx, &nn, q, candidatoBeta, stat, mu, sigma2, 
			    &loglikPartial);
      loglikCandidate += loglikPartial;
/*       delete [] xx; */
/*       delete [] yy; */
      Free(xx);
      Free(yy);
    }
    acepProb = acepProb + loglikCandidate - *loglikLast;
    /* correct jacobian for multiplicative random walk */
    for (i=0; i< *r * *r; ++i) {
      if (i % (*r+1)) acepProb = acepProb + log(candidatoBeta[i]) - log(beta[i]);
    }
    acepProb = exp(acepProb);
    if (acepProb > 1) acepProb = 1;
    if (runif(0,1) < acepProb) {
      for (i=0; i< (*r * *r); ++i) beta[i] = candidatoBeta[i];
      *loglikLast = loglikCandidate;
    }
  }
    
/*   delete [] q; */
/*   delete [] candidatoMu; */
/*   delete [] candidatoSigma2; */
/*   delete [] candidatoBeta; */
  Free(q);
  Free(candidatoMu);
  Free(candidatoSigma2);
  Free(candidatoBeta);

#ifdef DEBUG
  printf("\n Exiting MetropolisUpdate \n"); fflush(stdout);
#endif


}
void viterbi(double *y, double *x, int *genome, int *index, 
	     int *k, int *n, double *mu, double *sigma2,
	     double *beta, double *stat, int *states) {
  int i, j, l, g, nn;
  double *Q; Q = Calloc(*k * *k, double);
  double rowSumsQ = 0;
  for (g=0; g < *genome; ++g) {
    nn = index[g+1] - index[g];
    double *yy; yy = Calloc(nn, double);
    double *xx; xx = Calloc(nn, double);
    double *m; m = Calloc(nn* *k, double);
    int *b; b = Calloc(nn* *k, int);
    int *bAux; bAux = Calloc(*k, int);
    double *mAux; mAux = Calloc(*k, double);
    for (i=0; i < nn-1; ++i) {
      yy[i] = y[i + index[g]];
      xx[i] = x[i + index[g]];
    }
    xx[nn-1] = 0;
    yy[nn-1] = y[nn-1 + index[g]];

    /*  Forward recursion */

    for(j=0; j<*k; ++j) {
      m[j] = log(stat[j]) + dnorm(yy[0], mu[j], sqrt(sigma2[j]), 1);
    }
    for(i=1; i<nn; ++i) {
      for (j=0; j <*k; ++j) {
	rowSumsQ = 0;
	for(l=0; l < *k; ++l) {
	  Q[j* *k + l] = exp(beta[l * *k + j] * (xx[i-1] - 1));
	  rowSumsQ += Q[j * *k + l];
	}
	for(l=0; l < *k; ++l) {
	  Q[j * *k + l] = Q[j * *k + l] / rowSumsQ;
	}
      }
      for(j=0; j < *k; ++j) {
	for (l=0; l < *k; ++l) {
	  mAux[l] = m[(*k * (i-1)) + l] + log(Q[l * *k + j]);
	  bAux[l] = l + 1;
	}
	revsort(mAux, bAux, *k);
	b[i * *k + j] = bAux[0];
	m[i * *k + j] = mAux[0] + 
	  dnorm(yy[i], mu[j], sqrt(sigma2[j]), 1); 
      }
    }

    /*  Backward recursion */

    for (j=0; j < *k; ++j) {
      mAux[j] = m[((nn * *k)  - *k)  +j];
      bAux[j] = j + 1;
    }
    revsort(mAux, bAux, *k);
    states[index[g] + nn-1] = bAux[0];
    for (i=nn-2; i>=0; --i) {
      states[index[g] + i] = b[(i+1) * *k + states[index[g] + i + 1]-1];
    }
/*     delete [] xx; */
/*     delete [] yy; */
/*     delete [] m; */
/*     delete [] b; */
/*     delete [] bAux; */
/*     delete [] mAux; */

     Free(xx);
     Free(yy);
     Free(m);
     Free(b);
     Free(bAux);
     Free(mAux);



  }
  Free(Q);
/*   delete [] Q; */
}

/*preliminary method. It should be much faster*/
void probseqC(double *y, double *x, int *genome, int *index,
	      int *k, int *n, int *N, double *mu, double *sigma2, 
	      double *beta, double *stat, int *start, 
	      int *end, int *alteration, int *nAlteration, double *prob) {
  
  int i, j, r, cond;
  int *states;
  double *AuxMu, *AuxSigma2, *AuxBeta;
  double freq=0.0;

  states = (int *) R_alloc(*n, sizeof(int));
  AuxMu = (double *) R_alloc(*k, sizeof(double));
  AuxSigma2 = (double *) R_alloc(*k, sizeof(double)); 
  AuxBeta = (double *) R_alloc(*k * *k, sizeof(double));
  
  for(i=0; i < *N; i++) {
    for(j=0; j < *k; j++) {
      AuxMu[j] = mu[i* *k + j];
      AuxSigma2[j] = sigma2[i* *k + j];
    }
    for(j=0; j < *k * *k; j++) {
      AuxBeta[j] = beta[i * *k * *k + j];
    }
    viterbi(y, x, genome, index, k, n, AuxMu, AuxSigma2, AuxBeta, stat,
		    states);
    cond = 0;
    for (j= *start-1; j< *end; j++) {
      for (r=0; r< *nAlteration; r++) {
	if (states[j]==alteration[r]) cond ++;
      }
    }
    if (cond==(*end - *start + 1)) freq++;
  }
  *prob = freq / *N;
}

void wholeViterbi(double *y, double *x, int *genome, int *index,
		  int *k, int *n, int *N, double *mu, double *sigma2,
		  double *beta, double *stat, char **filename) {
  
/*   Format of the file: */
/*   every line show the state and the last spot with that state */
/*   this is repeated for every breakpoint */
/*   the last number of every line is a -1 followed by */
/*   the times that sequence is repeated */

  int i, j, r, cond, countrep, repeated;
  int *states, *lastSeq, lastState;
  double *AuxMu, *AuxSigma2, *AuxBeta;
  FILE *seqFile;
  states = (int *) R_alloc(*n, sizeof(int));
  lastSeq = (int *) R_alloc(*n, sizeof(int));
  AuxMu = (double *) R_alloc(*k, sizeof(double));
  AuxSigma2 = (double *) R_alloc(*k, sizeof(double));
  AuxBeta = (double *) R_alloc(*k * *k, sizeof(double));
  seqFile = fopen(*filename, "w");
  if (seqFile == NULL) {
    Rprintf("Can't open file for writing\n");
    exit(0);
  }
  else {
    i = 0;
    /*    First sequence */
    for(j=0; j < *k; j++) {
      AuxMu[j] = mu[i* *k + j];
      AuxSigma2[j] = sigma2[i* *k + j];
    }
    for(j=0; j < *k * *k; j++) {
      AuxBeta[j] = beta[i * *k * *k + j];
    }
    viterbi(y, x, genome, index, k, n, AuxMu, AuxSigma2, AuxBeta, stat,
	    states);
    lastState = states[0];
    lastSeq[0] = states[0];
    fprintf(seqFile, "%d\t", lastState);
    for (j=1; j < *n; j++) {
      lastSeq[j] = states[j];
      if (states[j] != lastState) {
	fprintf(seqFile, "%d\t%d\t", j, states[j]);
	lastState = states[j];
      }
    }
    lastSeq[*n-1] = states[*n-1];
    fprintf(seqFile, "%d\t", *n);
    countrep = 1;
    /* next sequences */
    for(i=1; i < *N; i++) {
      for(j=0; j < *k; j++) {
	AuxMu[j] = mu[i* *k + j];
	AuxSigma2[j] = sigma2[i* *k + j];
      }
      for(j=0; j < *k * *k; j++) {
	AuxBeta[j] = beta[i * *k * *k + j];
      }
      viterbi(y, x, genome, index, k, n, AuxMu, AuxSigma2, AuxBeta, stat,
	      states);
      repeated = 1;
      for (j=0; j<*n; j++) {
	if (lastSeq[j] != states[j]) {
	  repeated = 0;
	  break;
	}
      }
      if (repeated) {
	countrep ++;
      }
      else {
	fprintf(seqFile, "%d\t%d\n", -1, countrep);
	countrep = 1;
	lastState = states[0];
	lastSeq[0] = states[0];
	fprintf(seqFile, "%d\t", lastState);
	for (j=1; j < *n; j++) {
	  lastSeq[j] = states[j];
	  if (states[j] != lastState) {
	    fprintf(seqFile, "%d\t%d\t", j, states[j]);
	    lastState = states[j];
	  }
	}
	fprintf(seqFile, "%d\t", *n);
	lastSeq[*n-1] = states[*n-1];
      }
    }
  }
  fprintf(seqFile, "%d\t%d\n", -1, countrep);
  fclose(seqFile);
}




void edges(double *y, double *x, int *genome, int *index,
	   int *k, int *n, int *N, double *mu, double *sigma2,
	   double *beta, double *stat, int *count_edge) {
  
  int i, j, r, cond, countrep, repeated;
  int *states, *lastSeq, lastState;
  double *AuxMu, *AuxSigma2, *AuxBeta;
  
  states = (int *) R_alloc(*n, sizeof(int));
  
  lastSeq = (int *) R_alloc(*n, sizeof(int));
  AuxMu = (double *) R_alloc(*k, sizeof(double));
  AuxSigma2 = (double *) R_alloc(*k, sizeof(double));
  AuxBeta = (double *) R_alloc(*k * *k, sizeof(double));
  
  for (j=0; j < *n; j++) count_edge[j] = 0; /*initialize */
  
  for(i=0; i < *N; i++) {

    for(j=0; j < *k; j++) {
      AuxMu[j] = mu[i* *k + j];
      AuxSigma2[j] = sigma2[i* *k + j];
    }
    for(j=0; j < *k * *k; j++) {
      AuxBeta[j] = beta[i * *k * *k + j];
    }
    viterbi(y, x, genome, index, k, n, AuxMu, AuxSigma2, AuxBeta, stat,
	    states);
    for (j=1; j < *n; j++) {
      if(states[j - 1] != states[j]) {
	count_edge[j] += 1;
/* 	printf("\n Inside inner loop  i: %d", i); */
/* 	printf("    j: %d", j);  */
/* 	printf("\n                  count_edge[j]: %d", count_edge[j]); */
      }
    }
  }
}











/* void wholeViterbi(double *y, double *x, int *genome, int *index, */
/* 		  int *k, int *n, int *N, double *mu, double *sigma2, */
/* 		  double *beta, double *stat, char **filename) { */
  
/*   /\* Format of the file: *\/ */
/*   /\* every line show the state and the last spot with that state *\/ */
/*   /\* this is repeated for every breakpoint *\/ */
/*   /\* the last number of every line is a -1 followed by *\/ */
/*   /\* the times that sequence is repeated *\/ */

/*     int i, j, r, cond; */
/*     int *states, lastState; */
/*   double *AuxMu, *AuxSigma2, *AuxBeta; */
/*   FILE *seqFile; */
/*   states = (int *) R_alloc(*n, sizeof(int)); */
/*   AuxMu = (double *) R_alloc(*k, sizeof(double)); */
/*   AuxSigma2 = (double *) R_alloc(*k, sizeof(double)); */
/*   AuxBeta = (double *) R_alloc(*k * *k, sizeof(double)); */
/*   seqFile = fopen(*filename, "w"); */
/*   if (seqFile == NULL) { */
/*     Rprintf("Can't open file for writing\n"); */
/*     exit(0); */
/*   } */
/*   else { */
/*     for(i=0; i < *N; i++) { */
/*       for(j=0; j < *k; j++) { */
/* 	AuxMu[j] = mu[i* *k + j]; */
/* 	AuxSigma2[j] = sigma2[i* *k + j]; */
/*       } */
/*       for(j=0; j < *k * *k; j++) { */
/* 	AuxBeta[j] = beta[i * *k * *k + j]; */
/*       } */
/*       viterbi(y, x, genome, index, k, n, AuxMu, AuxSigma2, AuxBeta, stat, */
/* 	      states); */
/*       lastState = states[0]; */
/*       fprintf(seqFile, "%d\t", lastState); */
/*       for (j=1; j < *n-1; j++) { */
/* 	if (states[j] != lastState) { */
/* 	  fprintf(seqFile, "%d\t%d\t", j+1, states[j]); */
/* 	  lastState = states[j]; */
/* 	} */
/*       } */
/*       fprintf(seqFile, "%d\n", *n); */
/*     } */
/*   } */
/*   fclose(seqFile); */
/* } */
  


void MetropolisSweep(double *y, double *x, int *varEqual, int *genome, 
		     int *index, int *kMax, int *n, int *burnin,
		     int *TOT, int *times, int *burninTimes, int *probB,
		     int *probD, int *probS, int *probC, double *probK,
		     double *pb, double *ps,
		     double *muAlfa, double *muBeta,
		     double *sigmaTauMu, double *sigmaTauSigma2,
		     double *sigmaTauBeta, double *tauSplitMu,
		     double *tauSplitBeta,
		     int *k, double *mu, double *sigma2,
		     double *beta, double *stat, 
		     int *startK, int *RJ, double *maxVar, 
		     double *probStates, double *loglik)  {

  /*  Verify we can run */
/*   We assume at least 32-bit integers, which is good enough. FIXME: do the check*/
/*   unsigned long int max_int = (unsigned long int)>(std::numeric_limits<int>::max()); */
/*   if (max_int < 200000000) { */
/*     std::cout << "************   WARNING !!!!!!!!!! **************" << */
/*       "With this size of ints you can run into trouble. " << */
/*       "Modify the C++ code to use long ints (or change to another computer)." << std::endl; */
/*   } */
      
  /*  Initializa RNG */
  GetRNGstate();
#ifdef DEBUG
  double dummy_random_number;
  dummy_random_number = runif(0, 1);
  PR(dummy_random_number);
#endif
  
  int i,j,m;
  int t;
  int r, nn;
  int *states; states = Calloc(*n, int);
  int *accepted; accepted = Calloc(1, int);
  double *loglikLast; loglikLast = Calloc(*kMax, double);
  /* index to permutations */
  int *indexPerm; indexPerm = Calloc(*kMax, int);
  /* index to start every object */
  int *indexMu; indexMu = Calloc(*kMax, int); /* same for mu and sigma2 */
  int *indexBeta; indexBeta = Calloc(*kMax, int);
  int *indexStat; indexStat = Calloc(*kMax, int);
  /* old parameters (max limit) */
  double *OldStat; OldStat = Calloc(*kMax, double);
  double *OldMu; OldMu = Calloc(*kMax, double);
  double *OldSigma2; OldSigma2 = Calloc(*kMax, double);
  double *OldBeta; OldBeta = Calloc(*kMax * *kMax, double);
  /* new parameters (max limit) */

  double *NewStat; NewStat = Calloc(*kMax, double);
  double *NewMu; NewMu = Calloc(*kMax, double);
  double *NewSigma2; NewSigma2 = Calloc(*kMax, double);
  double *NewBeta; NewBeta = Calloc(*kMax * *kMax, double);
    
  /* auxiliaries for viterbi */
  double *AuxSigma2; AuxSigma2 = Calloc(*kMax, double);
  double *AuxBeta; AuxBeta = Calloc(*kMax * *kMax, double);

  double *q; q = Calloc(*kMax * *kMax, double);

  /*          Initializations */
  *accepted = 0;
  for (i=0; i<*kMax; ++i) indexPerm[i] = i+1;
  for (i=0; i<*kMax; ++i) loglikLast[i] = 0;
  for (i=0; i<*kMax; ++i) states[i] = -99;

  /* FIXME  do this appropriately: fill with true junk zz*/
  int junkindex = 0;
  for(junkindex = 0; junkindex < ((*kMax) * (*kMax)); ++junkindex) {
    q[junkindex] = -99;
    OldBeta[junkindex] = -99;
    NewBeta[junkindex] = -99;
    AuxBeta[junkindex] = -99;
  }
  for(junkindex = 0; junkindex < (*kMax) ; ++junkindex) {
    indexMu[junkindex] = -99;
    indexBeta[junkindex] = -99;
    indexStat[junkindex] = -99;
    OldStat[junkindex] = -99;
    OldMu[junkindex] = -99;
    OldSigma2[junkindex] = -99;
    NewStat[junkindex] = -99;
    NewMu[junkindex] = -99;
    NewSigma2[junkindex] = -99;
    AuxSigma2[junkindex] = -99;
  }


  indexMu[0] = 0;
  indexBeta[0] = 0;
  indexStat[0] = 0;
  for(i=1; i<*kMax; ++i) {
    indexMu[i] = (*TOT)*i + indexMu[i-1];
    indexBeta[i] = (*TOT)*i*i + indexBeta[i-1];
    indexStat[i] = i + indexStat[i-1];
  }
  /*  Overdispersed init */
  for (i=0; i<*kMax; ++i) {
    for (j=0; j<=i; ++j) {
      /*  this could be automatic */
      OldMu[j] = runif(-2, 2);
    }
    sigma2[indexMu[i]] = runif(0, *maxVar);
    for (j=1; j<=i; ++j) {
      if (varEqual) {
	sigma2[indexMu[i] + j] = sigma2[indexMu[j]];
      }
      else {
	sigma2[indexMu[i] + j] = runif(0, *maxVar);
      }
    }
    /*  relabelling */
    R_rsort(OldMu, i+1);
    for (j=0; j<=i; ++j) {
      mu[indexMu[i] + j] = OldMu[j];
    }
      
    for (j=0; j<((i+1)*(i+1)); ++j) {
      if (!(j % (i+2))) beta[indexBeta[i] +j] = 0;
      else beta[indexBeta[i] + j] = runif(0,100);
    }
  }

  /*  loglik of start values */
  for (i=1; i<=*kMax;++i) {
    for (j=0; j<i; ++j) {
      OldStat[j] = stat[indexStat[i-1] + j];
      OldMu[j] = mu[indexMu[i-1] + j];
      OldSigma2[j] = sigma2[indexMu[i-1] + j];
    }
    for (j=0; j <i*i; ++j) {
      OldBeta[j] = beta[indexBeta[i-1] + j];
      q[j] = -OldBeta[j];
    }
    loglikLast[i-1] = 0;
    for (m=0; m < *genome; ++m) {
      nn = index[m+1] - index[m];
      double *yy; yy = Calloc(nn, double);
      double *xx; xx = Calloc(nn, double);

      /*  Initialization */
      for (j=0; j < nn-1; ++j) {
	yy[j] = -9999;
	xx[j] = -9999;
      }


      double loglikPartial = 0;
      for (j=0; j < nn-1; ++j) {
#ifdef DEBUG
	printf("\n Inside MSweep; y[j + index[m]]: %f", y[j + index[m]]);
	printf("\n Inside MSweep; x[j + index[m]]: %f", x[j + index[m]]); /* esto está pocho */
#endif
	yy[j] = y[j + index[m]];
	xx[j] = x[j + index[m]];
      }
      xx[nn-1] = 0;
      yy[nn-1] = y[index[m + 1] - 1];

#ifdef DEBUG
      printf("\n    Calling normal NHHMMlikelihhod from MetropolisSweep\n");
      printf("\n nn eb MetropolisSweep:  %d\n", nn); /*  zz: cómo que 20!!! */
#endif
      /*       #ifdef DEBUG */
      /*  		if (i == 0) printf("\n x_0:  %f", x_i); */
      /*  		if (i == 1) printf("\n x_1:  %f", x_i); */
      /*  		if (i == 2) printf("\n x_2:  %f", x_i); */
      /*  		if (i == 3) printf("\n x_3:  %f", x_i); */
      /*  		if (i == 4) printf("\n x_4:  %f", x_i); */
      /*  		if (i == 5) printf("\n x_5:  %f", x_i); */
      /*  		if (i == 6) printf("\n x_6:  %f", x_i); */
      /*  		if (i == 7) printf("\n x_7:  %f", x_i); */
      /*  		if (i == 8) printf("\n x_8:  %f", x_i); */
      /*  		if (i == 9) printf("\n x_9:  %f", x_i); */
      /*  		if (i == 10) printf("\n x_10:  %f", x_i); */
      /*  		if (i == 11) printf("\n x_11:  %f", x_i); */
      /*  		if (i == 12) printf("\n x_12:  %f", x_i); */
      /*  		if (i == 13) printf("\n x_13:  %f", x_i); */
      /*  		if (i == 14) printf("\n x_14:  %f", x_i); */
      /*  		if (i == 15) printf("\n x_15:  %f", x_i); */
      /*  		if (i == 16) printf("\n x_16:  %f", x_i); */
      /*  		if (i == 17) printf("\n x_17:  %f", x_i); */
      /*  		if (i == 18) printf("\n x_18:  %f", x_i); */
      /*  		if (i == 19) printf("\n x_19:  %f", x_i); */
      /*  		if (i == 20) printf("\n x_20:  %f", x_i); */
      /*               #endif */



      normalNHHMMlikelihood(yy, &i, xx, &nn, q, OldBeta, OldStat, OldMu, OldSigma2, 
			    &loglikPartial);
      loglikLast[i-1] += loglikPartial;
/*       delete [] xx; */
/*       delete [] yy; */
      Free(xx);
      Free(yy);
      
    }

    loglik[(i-1)* *TOT] = loglikLast[i-1];
    /*     viterbi*/
    if (i >1) {
      viterbi(y, x, genome, index, &i, n, OldMu, OldSigma2, OldBeta, OldStat,
	      states);
      for (m=0; m < i-1; ++m) {
	for (j=0; j < *n; ++j) {
	  /* Marginal probabilities */
	  if(states[j]==m+1) {
	    probStates[(*n * (m + (i-1) * (i-2) / 2)) + j] =
	      (1.0 / (double)(times[i-1]+1)) +
	      probStates[(*n * (m + (i-1) * (i-2) / 2)) + j] *
	      (double)(times[i-1]+1 -1) / (double)(times[i-1]+1);
	  }
	  else {
	    probStates[(*n * (m + (i-1) * (i-2) / 2)) + j] =
	      probStates[(*n * (m + (i-1) * (i-2) / 2)) + j] *
	      (double)(times[i-1]+1 -1) / (double)(times[i-1]+1);
	  }
	}
      }
    }
  }
  
  if (*startK==0) r = (int)rint(runif(1, *kMax));
  else r = *startK;

  k[0] = r;
  /*  Next free position */
  for (i=0; i <*kMax; ++i) times[i] = 1;

  /*  Loop MCMC iterations */
  for(t=1; t<*TOT; ++t) {

    /* Allow R interrupts; check every 100 iterations */
    if (!(t % 100))
      R_CheckUserInterrupt(); 

    /* METROPOLIS UPDATE */
    /* old parameters */
    for (i=0; i<r; ++i) {
      OldStat[i] = stat[indexStat[r-1] + i];
      OldSigma2[i] = sigma2[indexMu[r-1] + ((times[r-1]-1)*r) + i];
      OldMu[i] = mu[indexMu[r-1] + ((times[r-1]-1)*r) + i];
    }
    for (i=0; i<r * r; ++i) {
      OldBeta[i] = beta[indexBeta[r-1] + ((times[r-1]-1)*r*r) + i];
    }
    PR(r-1);
    PR(sigmaTauMu[r-1]);
    MetropolisUpdate(y, x, varEqual, genome, index, OldMu, 
		     OldSigma2, OldBeta, OldStat, &r, 
		     n, muAlfa, muBeta, &sigmaTauMu[r-1], &sigmaTauSigma2[r-1], 
		     &sigmaTauBeta[r-1], &loglikLast[r-1], maxVar);
    loglik[*TOT * (r-1) + times[r-1]] = loglikLast[r-1];
    rsort_with_index(OldMu, indexPerm, r);

    /*  Auxiliaries correctly permuted for viterbi */
    for (i=0; i< r; ++i) 
      AuxSigma2[i] = OldSigma2[indexPerm[i]-1];

    for (i=0; i<r; ++i) 
      for (j=0; j<r; ++j) 
	AuxBeta[i*r +j] = 
	  OldBeta[(indexPerm[i]-1)*r + indexPerm[j]-1];


    /*  viterbi */
    if ((t > *burnin) && (r >1)) {
      viterbi(y, x, genome, index, &r, n, OldMu, AuxSigma2, AuxBeta, OldStat,
	      states);
      for (i=0; i < r-1; ++i) {
	for (j=0; j < *n; ++j) {
	  if(states[j]==i+1) {
	    probStates[(*n * (i + (r-1) * (r-2) / 2)) + j] =
	      (1.0 / (double)(times[r-1]+1 - burninTimes[r-1])) +
	      probStates[(*n * (i + (r-1) * (r-2) / 2)) + j] *
	      (double)(times[r-1]+1 - burninTimes[r-1]-1) /
	      (double)(times[r-1]+1 - burninTimes[r-1]);
	  }
	  else {
	    probStates[(*n * (i + (r-1) * (r-2) / 2)) + j] =
	      probStates[(*n * (i + (r-1) * (r-2) / 2)) + j] *
	      (double)(times[r-1]+1 - burninTimes[r-1] -1) /
	      (double)(times[r-1]+1 - burninTimes[r-1]);
	  }
	}
      }
    }
  
    /*  save updates */
    for (i=0; i <r; ++i) {
      mu[indexMu[r-1] + (times[r-1]*r) + i] = OldMu[i];
      sigma2[indexMu[r-1] + (times[r-1]*r) + i] = OldSigma2[indexPerm[i]-1];
    }
    for (i=0; i<r; ++i) {
      for (j=0; j<r; ++j) {
	beta[indexBeta[r-1] + (times[r-1]*r*r) + i*r +j] = OldBeta[(indexPerm[i]-1)*r + indexPerm[j]-1];
	q[i*r +j] = -OldBeta[i*r + j];
      }
    }
    for(i=0; i<*kMax; ++i) indexPerm[i]= i+1;
    times[r-1]++;
    if (*RJ) {
      /*  Birth or death */
      if(runif(0,1) <= pb[r-1]) {
	for (i=0; i<(r+1); ++i) NewStat[i] = stat[indexStat[r] + i];
	Birth(y, varEqual, genome, index, OldMu, OldSigma2, 
	      OldBeta, OldStat, NewStat, 
	      &r, &loglikLast[r-1], probK,
	      pb, muAlfa, muBeta, x, n, NewMu, NewSigma2, NewBeta,
	      &loglikLast[r], accepted, maxVar);
	  
	if (*accepted) {
	  r = r+1;
	  *probB = *probB + 1;
	  *accepted = 0;
	  /* Relabelling */
	  rsort_with_index(NewMu, indexPerm, r);

	  /*  Auxiliaries correctly permuted for viterbi */
	  for (i=0; i< r; ++i) 
	    AuxSigma2[i] = NewSigma2[indexPerm[i]-1];

	  for (i=0; i<r; ++i) 
	    for (j=0; j<r; ++j) 
	      AuxBeta[i*r +j] = 
		NewBeta[(indexPerm[i]-1)*r + indexPerm[j]-1];
	    
	  /* 	    viterbi  */
	  if ((t > *burnin) && (r >1)) {
	    viterbi(y, x, genome, index, &r, n, NewMu, AuxSigma2, AuxBeta, NewStat,
		    states);
	    for (i=0; i < r-1; ++i) {
	      for (j=0; j < *n; ++j) {
		if(states[j]==i+1) {
		  probStates[(*n * (i + (r-1) * (r-2) / 2)) + j] =
		    (1.0 / (double)(times[r-1]+1 - burninTimes[r-1])) +
		    probStates[(*n * (i + (r-1) * (r-2) / 2)) + j] *
		    (double)(times[r-1]+1 - burninTimes[r-1] -1) /
		    (double)(times[r-1]+1 - burninTimes[r-1]);
		}
		else {
		  probStates[(*n * (i + (r-1) * (r-2) / 2)) + j] =
		    probStates[(*n * (i + (r-1) * (r-2) / 2)) + j] *
		    (double)(times[r-1]+1 - burninTimes[r-1] -1) /
		    (double)(times[r-1]+1 - burninTimes[r-1]);
		}
	      }
	    }
	  }
	
	  for (i=0; i<r; ++i) {
	  mu[indexMu[r-1] + (times[r-1]*r) + i] = NewMu[i];
	  sigma2[indexMu[r-1] + (times[r-1]*r) +i] = NewSigma2[indexPerm[i]-1];
	}
	for (i=0; i<r; ++i) {
	    for(j=0; j<r; ++j) {	  
	      beta[indexBeta[r-1] + (times[r-1]*r*r) + i* r + j] = 
		NewBeta[(indexPerm[i]-1)* r + indexPerm[j]-1];
	    }
	  }
	  /* Recover original permutations */
	  for (i=0; i<*kMax; ++i) indexPerm[i] = i+1;
	  loglik[*TOT * (r-1) + times[r-1]] = loglikLast[r-1];
	  times[r-1]++;
	}
      }
      else {
	for (i=0; i<(r-1); ++i) 	
	  NewStat[i] = stat[indexStat[r-2] + i];
	Death(y, genome, index, OldMu, OldSigma2, OldBeta, OldStat, NewStat,
	      &r, &loglikLast[r-1], probK, pb, x, n, NewMu, NewSigma2, 
	      NewBeta, &loglikLast[r-2], accepted);
	  
	if (*accepted) {
	  r = r-1;
	  *accepted = 0;
	  *probD = *probD + 1;
/* 	   viterbi  */
	  if ((t > *burnin) && (r >1)) {
	    viterbi(y, x, genome, index, &r, n, NewMu, NewSigma2, NewBeta, NewStat,
		    states);

	    for (i=0; i < r-1; ++i) {
	      for (j=0; j < *n; ++j) {
		if(states[j]==i+1) {
		  probStates[(*n * (i + (r-1) * (r-2) / 2)) + j] =
		    (1.0 / (double)(times[r-1]+1 - burninTimes[r-1])) +
		    probStates[(*n * (i + (r-1) * (r-2) / 2)) + j] *
		    (double)(times[r-1]+1 - burninTimes[r-1] - 1) /
		    (double)(times[r-1]+1 - burninTimes[r-1]);
		}
		else {
		  probStates[(*n * (i + (r-1) * (r-2) / 2)) + j] =
		    probStates[(*n * (i + (r-1) * (r-2) / 2)) + j] *
		    (double)(times[r-1]+1 - burninTimes[r-1] -1) /
		    (double)(times[r-1]+1 - burninTimes[r-1]);
		}
	      }
	    }
	  }
	
	  for (i=0; i<r; ++i) {
	    mu[indexMu[r-1] + (times[r-1]*r) + i] = NewMu[i];
	    sigma2[indexMu[r-1] + (times[r-1]*r) + i] = NewSigma2[i];
	  }
	  for (i=0; i<r*r; ++i) {
	    beta[indexBeta[r-1] + (times[r-1]*r*r) + i] = NewBeta[i];
	  }
	  loglik[*TOT * (r-1) + times[r-1]] = loglikLast[r-1];
	  times[r-1]++;
	}
	  
      }
      k[2*t-1] = r;
	
      /*  Split or combine */
      /*  Old parameters */
      for (i=0; i<r; ++i) {
	OldStat[i] = stat[indexStat[r-1] + i];
	OldSigma2[i] = sigma2[indexMu[r-1] + ((times[r-1]-1)*r) + i];
	OldMu[i] = mu[indexMu[r-1] + ((times[r-1]-1)*r) + i];
      }
      for (i=0; i<r * r; ++i) {
	OldBeta[i] = beta[indexBeta[r-1] + ((times[r-1]-1)*r*r) + i];
	q[i] = -OldBeta[i];
      }
	
      if(runif(0,1) <= ps[r-1]) {
	  
	for (i=0; i<(r+1); ++i) NewStat[i] = stat[indexStat[r] + i];
	  
	Split(y, varEqual, genome, index, OldMu, OldSigma2, 
	      OldBeta, OldStat, NewStat, 
	      &r, &loglikLast[r-1], probK, ps, x, n, NewMu, 
	      NewSigma2, NewBeta, &loglikLast[r], 
	      muAlfa, muBeta, tauSplitMu,
	      tauSplitBeta, accepted);
	  
	if (*accepted) {
	  r = r+1;
	  *accepted = 0;
	  *probS = *probS + 1;
	  rsort_with_index(NewMu, indexPerm, r);

	  /*  Auxiliaries correctly permuted for viterbi */
	  for (i=0; i< r; ++i) 
	    AuxSigma2[i] = NewSigma2[indexPerm[i]-1];
	    
	  for (i=0; i<r; ++i) 
	    for (j=0; j<r; ++j) 
	      AuxBeta[i*r +j] = 
		NewBeta[(indexPerm[i]-1)*r + indexPerm[j]-1];

	  /* 	    viterbi  */
	  if ((t > *burnin) && (r >1)) {
	    viterbi(y, x, genome, index, &r, n, NewMu, AuxSigma2, AuxBeta, NewStat,
		    states);
	    for (i=0; i < r-1; ++i) {
	      for (j=0; j < *n; ++j) {
		if(states[j]==i+1) {
		  probStates[(*n * (i + (r-1) * (r-2) / 2)) + j] =
		    (1.0 / (double)(times[r-1]+1 - burninTimes[r-1])) +
		    probStates[(*n * (i + (r-1) * (r-2) / 2)) + j] *
		    (double)(times[r-1]+1 - burninTimes[r-1] -1) /
		    (double)(times[r-1]+1 - burninTimes[r-1]);
		}
		else {
		  probStates[(*n * (i + (r-1) * (r-2) / 2)) + j] =
		    probStates[(*n * (i + (r-1) * (r-2) / 2)) + j] *
		    (double)(times[r-1]+1 - burninTimes[r-1] -1) /
		    (double)(times[r-1]+1 - burninTimes[r-1]);
		}
	      }
	    }
	  }
	
	  for (i=0; i<r; ++i) {
	    mu[indexMu[r-1] + (times[r-1]*r) + i] = NewMu[i];
	    sigma2[indexMu[r-1] + (times[r-1]*r) + i] = NewSigma2[indexPerm[i]-1];
	  }
	  for (i=0; i<r; ++i) {
	    for (j=0; j<r; ++j) {
	      beta[indexBeta[r-1] + (times[r-1]*r*r) + i*r +j] = NewBeta[(indexPerm[i]-1)* r + indexPerm[j]-1];
	    }
	  }
	  for(i=0; i<*kMax; ++i) indexPerm[i] = i+1;
	  loglik[*TOT * (r-1) + times[r-1]] = loglikLast[r-1];
	  times[r-1]++;
	}
	  
      }
      else {
	for (i=0; i<(r-1); ++i) NewStat[i] = stat[indexStat[r-2] + i];
	  
	Combine(y, varEqual, genome, index, OldMu, OldSigma2, 
		OldBeta, OldStat, NewStat, 
		&r, &loglikLast[r-1], probK, ps, x, n, NewMu, 
		NewSigma2, NewBeta, &loglikLast[r-2], muAlfa, 
		muBeta, tauSplitMu,
		tauSplitBeta, accepted, maxVar);
	  
	if (*accepted) {
	  r = r-1;
	  *probC = *probC +1;
	  *accepted = 0;
/* 	    viterbi  */
	  if ((t > *burnin) && (r >1)) {
	    viterbi(y, x, genome, index, &r, n, NewMu, NewSigma2, NewBeta, NewStat,
		    states);
	    for (i=0; i < r-1; ++i) {
	      for (j=0; j < *n; ++j) {
		if(states[j]==i+1) {
		  probStates[(*n * (i + (r-1) * (r-2) / 2)) + j] =
		    (1.0 / (double)(times[r-1]+1 - burninTimes[r-1])) +
		    probStates[(*n * (i + (r-1) * (r-2) / 2)) + j] *
		    (double)(times[r-1]+1 - burninTimes[r-1] -1) /
		    (double)(times[r-1]+1 - burninTimes[r-1]);
		}
		else {
		  probStates[(*n * (i + (r-1) * (r-2) / 2)) + j] =
		    probStates[(*n * (i + (r-1) * (r-2) / 2)) + j] *
		    (double)(times[r-1]+1 - burninTimes[r-1] -1) /
		    (double)(times[r-1]+1 - burninTimes[r-1]);
		}
	      }
	    }
	  }
	
	  for (i=0; i<r; ++i) {
	    mu[indexMu[r-1] + (times[r-1]*r) + i] = NewMu[i];
	    sigma2[indexMu[r-1] + (times[r-1]*r) +i] = NewSigma2[i];
	  }
	  for (i=0; i<r*r; ++i) {
	    beta[indexBeta[r-1] + (times[r-1]*r*r) + i] = NewBeta[i];
	  }
	  loglik[*TOT * (r-1) + times[r-1]] = loglikLast[r-1];
	  times[r-1]++;
	}
	
      }
      k[2*t] = r;

      /* check if it's burnin' time! */
      }
    if (t==*burnin) {
      for (i=0; i<*kMax; ++i) {
	burninTimes[i] = times[i];
      }
    }
    }
/*   delete []states; */
/*   delete []loglikLast; */
/*   delete []accepted; */
/*   delete []q; */
/*   delete []indexMu; */
/*   delete []indexBeta; */
/*   delete []indexStat; */
/*   delete []OldStat; */
/*   delete []OldMu; */
/*   delete []OldSigma2; */
/*   delete []OldBeta; */
/*   delete []NewStat; */
/*   delete []NewMu; */
/*   delete []NewSigma2; */
/*   delete []NewBeta; */
/*   delete []AuxSigma2; */
/*   delete []AuxBeta; */
/*   delete []indexPerm; */



  Free(states);
  Free(loglikLast);
  Free(accepted);
  Free(q);
  Free(indexMu);
  Free(indexBeta);
  Free(indexStat);
  Free(OldStat);
  Free(OldMu);
  Free(OldSigma2);
  Free(OldBeta);
  Free(NewStat);
  Free(NewMu);
  Free(NewSigma2);
  Free(NewBeta);
  Free(AuxSigma2);
  Free(AuxBeta);
  Free(indexPerm);
  
  PutRNGstate();
  }

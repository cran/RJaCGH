/* Copyright (C) 2005-2008  Oscar Rueda Palacio and Ramon Diaz-Uriarte */

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


// FIXME:
// Quitar toas las alocaciones a yy, xx temporales. (los rodeados por P1, P2)
// prob.states: no calcularlo? hacerlo como viterbi?
// beta: solo mediana
  


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Utils.h>

#include<limits.h>

#include <zlib.h>

/*  Modified from McConnells 'Complete Code' */
#define ASSERT(condition, message) {       \
  if( !(condition)) {                      \
    Rprintf("\n !! Assertion failed :");    \
    Rprintf( #condition );                  \
    Rprintf( message);                      \
    error("\nFatal error in RJaCGH. Please let us know about it.\n");                               \
  }                                        \
}                                          \


// From include/Rmath.h0:
#define M_1_SQRT_2PI 0.398942280401432677939946059934        /* 1/sqrt(2pi) */



// FIXME: compress mu, sigma, etc, too.

#define PR4(x) {Rprintf("\n %s = %f \n", #x,  (float) x); fflush(stdout);}
#define PR5(x, y) {Rprintf("\n At %s %s = %f \n", #y, #x,  (float) x); fflush(stdout);}
#define P1(x) {Rprintf("\n        ALLOCATION TO DO  ....... %s \n", #x); fflush(stdout);}
#define P2(x) {Rprintf("\n        ALLOCATION  DONE  ....... %s \n", #x); fflush(stdout);}



#ifdef DEBUGW
 #define PR(x) {Rprintf("\n %s = %f \n", #x,  (float) x); fflush(stdout);}
#else
 #define PR(x);
#endif


#ifdef DEBUGW
 #define PR2(x) {Rprintf("\n %s = %p \n", #x,  x); fflush(stdout);}
#else
 #define PR2(x);
#endif

#ifdef DEBUGW
 #define PR3(x) {Rprintf("\n %s = %i \n", #x,  x); fflush(stdout);}
#else
 #define PR3(x);
#endif


#define CHECK_NUM(x, y) {			\
    if(!R_FINITE(x)) {				\
      Rprintf("\n ERROR in %s :", #y);		\
      Rprintf("   !R_FINITE %s \n", #x);		\
      error("\nFatal error in RJaCGH. Please let us know.\n");					\
    }						\
  }                                            

/* Blows up in both laptop and desktop CHECK_NUM(1.0/10E-310, A); */


struct Sequence {
  // For writing out the Viterbi sequences
  // collapsed or compacted sequence
  int k;
  int MCMC_num;
  int len;
  struct Sequence *next;
  unsigned int state_count[];
};


struct regionS {
  // To store pREC-S
  int Start;
  int End;
  int num_arrays; 
  int region_number; //yes, redundant but cleanest way to keep track
  int sum_num_arrays; //ditto
  struct regionS *next;
  int arrays[];
};


// Yes, a global variable. To deal with a second call
// to C from R. This way, we use a different, simpler,
// function than pREC-S.

static struct regionS *regS = NULL;



//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//        Functions for handling pREC-S                 //
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

void add_regionS(struct regionS **regionRef,
		 const int Start,
		 const int End,
		 const int num_arrays,
		 const int *arrays) {
  struct regionS *new_region_ptr;
  // We do NOT want R_alloc, because we want the structure
  // to persist after the first .C call, and not be affected
  // by garbage collection.
  /*   new_region_ptr = (struct regionS *) R_alloc(1, */
/* 					      sizeof(struct regionS) + */
/* 					      num_arrays * sizeof(int)); */

  new_region_ptr = malloc( sizeof(struct regionS) +
					      num_arrays * sizeof(int));
  
  if(new_region_ptr == NULL) {
    Rprintf("\n ERROR allocating new regionS\n");
    error("\nFatal error in RJaCGH. Please let us know.\n");
  }
  new_region_ptr->num_arrays = num_arrays;
  new_region_ptr->Start = Start;
  new_region_ptr->End = End;
  for(int i = 0; i < num_arrays; i++) {
    new_region_ptr->arrays[i] = arrays[i];
  }
  if((*regionRef) == NULL) {
    new_region_ptr->region_number = 1;
    new_region_ptr->sum_num_arrays = num_arrays;
  } else { 
    new_region_ptr->region_number = ((*regionRef)->region_number) + 1;
    new_region_ptr->sum_num_arrays = ((*regionRef)->sum_num_arrays) + num_arrays;
  }
  new_region_ptr->next = *regionRef;
  *regionRef = new_region_ptr;
}


int is_region_subset(struct regionS *region,
		     const int start,
		     const int end,
		     const int narrays,
		     const int *arrays) {
  // pREC-S: is region just found a subset of a
  // previously found one?
  if(start > (region->End)) return(0);
  if(end > (region->End)) return(0);
  if(narrays > (region->num_arrays)) return(0);
  //Arrays are always sorted or stored in same order.
  for(int i = 0; i < narrays; i++) {
    if(arrays[i] != (region->arrays[i])) return(0);
  }
  return(1);
}


void updateRegionS(struct regionS **regionRef,
		   const int Start,
		   const int End,
		   const int num_arrays,
		   const int *arrays) {
  // Add a new pREC-S region to regionRef,
  // but only if the candidate sequence to add is not
  // a subset of an existing region in regionRef.
  int contained = 0;
  struct regionS *current = *regionRef;
  while(1) {
    if(current == NULL) {
      add_regionS(regionRef, Start, End, num_arrays, arrays);
      break;
    }
    contained = is_region_subset(current, Start, End, num_arrays, arrays);
    if(contained) break;
    current = current -> next;
  }
}


void return_pREC_S(int *regionsStart, int *regionsEnd,
		   int *regionsNarrays, int *allArrays) {
  // Called from R; returns start, end, arrays.
  // And we free the memory, since we malloc'ed it.
  int i = 0;
  int j = 0;
  void *last_ptr;
  while(regS != NULL) {
    regionsStart[i] = (regS->Start);
    regionsEnd[i]   = (regS->End);
    regionsNarrays[i] = (regS->num_arrays);

    // could use sum_num_arrays instead of j
    for(int m = 0; m < (regS->num_arrays); m++) {
      allArrays[j + m] = regS->arrays[m];
    }
    j += (regS->num_arrays);
    i++;
    last_ptr = regS;
    regS = regS->next;
    free(last_ptr);
  }
}


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//        Functions for Viterbi seqs                    //
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////


void addSeq(struct Sequence **seqRef, 
	    const int k, 
	    const int len,
	    unsigned int *state_count) {
  // Adds a new compacted sequence with given state_count data. 
  // MCMC counter is initialized (to 1).

  struct Sequence  *new_seq_ptr;
  // If not using R's memory alloc. facilites, use next line.
  // (Would also need to deallocate).
  /*   new_seq_ptr = malloc( sizeof(struct Sequence) + len * sizeof(unsigned int)); */
  new_seq_ptr = (struct Sequence *) R_alloc(1, sizeof(struct Sequence) +
					    len * sizeof(unsigned int));
  if(new_seq_ptr == NULL) {
    Rprintf("\n ERROR allocating new Sequence\n");
    error("\nFatal error in RJaCGH. Please let us know.\n");
  }
  new_seq_ptr->len = len;
  new_seq_ptr->k = k;
  new_seq_ptr->MCMC_num = 1;
  for(int i = 0; i < len; i++) {
    new_seq_ptr->state_count[i] = state_count[i];
  }
  new_seq_ptr->next = *seqRef;
  *seqRef = new_seq_ptr;
}


int CompareSeq_tmp(struct Sequence *seq, const int tmp_k,
		   const int tmp_len,
		   unsigned int *tmp_state_count) {
  // Is a vector of state_counts equal to a compacted sequence?
  if((seq->k != tmp_k) || (seq->len != tmp_len)) return(0);
  for(int i = 0; i < tmp_len; i++) {
    if(seq->state_count[i] != tmp_state_count[i])
      return(0);
  }
  return(1);
}


void Compare_Add_Seq(struct Sequence **seq, const int tmp_k,
		     const int tmp_len,
		     unsigned int *tmp_state_count) {
  // Compare a tmp vector of state_counts to the already processed sequences.
  // If same as a previous one, increase that one's counter.
  // Otherwise (i.e., we reach end of linked list without finding that sequence)
  // add to the set of sequences.
  int same_seq;
  struct Sequence *current = *seq;
  
  while(1) {
    if(current == NULL) {
      addSeq(seq, tmp_k, tmp_len, tmp_state_count);
      break;
    }
    same_seq = CompareSeq_tmp(current, tmp_k, tmp_len, tmp_state_count);
    if(same_seq) {
      current->MCMC_num++;
      break;
    }
    current = current->next;
  }
}


void create_constant_Sequence(struct Sequence **seqRef, const int n) {
  // When number of hidden states = 1, create and compare this sequence
  unsigned int tmp_state_count[2] = {1, n};
  Compare_Add_Seq(seqRef, 1, 2, tmp_state_count);
}

void add_constant_Sequence(int *genome, int *index, 
			   int *n, struct Sequence **seqRef, 
			   struct Sequence **arraySeqRefs,
			   int write_seq, unsigned int *viterbi_counts) {
  if(write_seq) {
    viterbi_counts[0]++;
    if((*genome) > 1) {
      int g, nn, index_g;
      for (g=0; g < (*genome); g++) {
	index_g = index[g];
	nn = index[g+1] - index_g;
	create_constant_Sequence(arraySeqRefs + g, nn); 
      }
    } else {
      create_constant_Sequence(seqRef, *n);
    }
  }
}


void viterbi_to_Sequence(struct Sequence **seqRef,  const int k,  
			 const int n, int *states) {
  unsigned int * restrict tmp_state_count;
  tmp_state_count = (unsigned int *) R_alloc(2 * n, sizeof(unsigned int));
  int lastState;
  int i = 0;

  lastState = states[0];
  for(int j = 1; j < n; j++) {
    if(lastState != states[j]) {
      tmp_state_count[i] = lastState;
      tmp_state_count[i + 1] = j;
      lastState = states[j];
      i += 2;
    }
  }
  tmp_state_count[i] = lastState;
  tmp_state_count[i + 1] = n;

  Compare_Add_Seq(seqRef, k, i + 2, tmp_state_count);
}

void Sequence_to_gzfile(const char *path, struct Sequence *seq,
			int *num_sequences, unsigned int *viterbi_counts,
			unsigned int *k_sum, unsigned int Total_k) {
  gzFile outfile;
  outfile = gzopen(path, "wb9");
  int closing, writing;

  while(1) {
    if(seq == NULL) {
      closing = gzclose(outfile);
      if(closing < 0) Rprintf("\n ERROR writing Sequence to gz file\n");
      break;
    }
    writing = gzprintf(outfile, "%i ", seq->k);
    if(writing < 0) Rprintf("\n ERROR writing value to gz file\n");
    writing = gzprintf(outfile, "%i ", seq->MCMC_num);
    if(writing < 0) Rprintf("\n ERROR writing value to gz file\n");
    writing = gzprintf(outfile, "%i ", viterbi_counts[(seq->k) - 1]);
    if(writing < 0) Rprintf("\n ERROR writing value to gz file\n");
    writing = gzprintf(outfile, "%i ", k_sum[(seq->k) - 1]);
    if(writing < 0) Rprintf("\n ERROR writing value to gz file\n");
    writing = gzprintf(outfile, "%i ", Total_k);
    if(writing < 0) Rprintf("\n ERROR writing value to gz file\n");
    
    for(int i = 0; i < ((seq->len) - 1); i++) {
      writing = gzprintf(outfile, "%i ", seq->state_count[i]);
      if(writing < 0) Rprintf("\n ERROR writing value to gz file\n");
    }
    writing = gzprintf(outfile, "%i\n", seq->state_count[(seq->len) - 1]);
    if(writing < 0) Rprintf("\n ERROR writing value to gz file\n");

    seq = seq->next;
    (*num_sequences)++;
  }
}



//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////


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

void Free_2_int(int **x, int nrow) { /* Free 2-D array */
  int i; 
  for (i = 0; i < nrow; ++i) { 
    Free(x[i]);
  } 
  Free(x); 
}


double state_to_prob(const int state, const int k, const int alteration,
		     const double *state_probs, const double prob_seq) {
  // alteration: +1: Gain; -1: Loss
  // To get a value: add base index of the k * 3 subarray, then row, then
  // move left or right according to alteration.
  // the "base index" is obtained from sum of arithmetic progression to n terms
  // (Abramowitz and Stegun, p. 10; a = 0; d = 3).
  // state_probs[ (3 * k * (k - 1) / 2) + (state - 1 ) * 3 + (1 + alteration)]
  
  // write down some explicit cases
  // FIXME: we could "unroll", separating an initial if with alteration = 1
  // and alteration = -1, and avoid the addition.

  // FIXME: we multiply tmp by prob_seq, but we could just have
  // multiplied the state_probs initially. But prob_seq is computed
  // by row. Hummm... So we'd end up with different state_probs per row. 
  double tmp = 0.0;
  if(k == 1) tmp = state_probs[1 + alteration];
  else if(k == 2) tmp = state_probs[3 * state + 1 + alteration];
  else if(k == 3) tmp = state_probs[3 * state + 7 + alteration];
  else if(k == 4) tmp = state_probs[3 * state + 16 + alteration];
  else if(k == 5) tmp = state_probs[3 * state + 28 + alteration];
  else tmp = state_probs[ (3 * k * (k - 1) / 2) + (state - 1 ) * 3 + (1 + alteration)];
  tmp *= prob_seq;
  CHECK_NUM(tmp, state_to_prob);
  return tmp;
}


void read_convert_prob_seq(double *prob_row_seq,
			   double **stretched,
			   const char *path,
			   const double *state_probs,
			   const int alteration,
			   const int num_sequences) {
  int k;
  int file_row, column;
  int count_viterbi, count_k, sum_viterbi, sum_k;
  char SEPCHARS[] = " \t";
  int SIZE_READ = 999999;
  char inb[SIZE_READ];
  char *token;
  char *token2;
  gzFile infile;

  if(!( (alteration == 1) || (alteration == -1))) {
    Rprintf("\n ERROR in read_convert_prob_seq: alteration != {-1, +1}\n");
    error("\nFatal error in RJaCGH. Please let us know.\n");
  }
  
  infile = gzopen(path, "rb");
  file_row = 0;
  while(1) {
    gzgets(infile, inb, SIZE_READ);
    token = strtok(inb, SEPCHARS);
    // if(gzeof(infile) == 1) break; won't work: can stop one line before last.
    if(token == NULL) break;
    if(token[0] == '\n') break;

    column = 0;
    k = atoi(token);
    token = strtok(NULL, SEPCHARS);
    count_viterbi = atoi(token);
    token = strtok(NULL, SEPCHARS);
    sum_viterbi = atoi(token);
    token = strtok(NULL, SEPCHARS);
    count_k = atoi(token);
    token = strtok(NULL, SEPCHARS);
    sum_k = atoi(token);

    prob_row_seq[file_row] =
      ((float) (count_viterbi * count_k) / (sum_viterbi * sum_k)); // * prior_k; 
    
    CHECK_NUM(prob_row_seq[file_row], read_convert_prob_seq);

    while(1) {
      token = strtok(NULL, SEPCHARS);
      if(token == NULL){
	break;
      }
      token2 = strtok(NULL, SEPCHARS);
      if(token2 == NULL) {
	Rprintf("\n ERROR: there should not be a lonely token!\n");
	error("\nFatal error in RJaCGH. Please let us know.\n");
      }
      // We read a pair; first token is state, second is last position
      // of that state. Fill up to that column, and translate to probs.
      for(; column < atoi(token2); column++) {
	stretched[file_row][column] = state_to_prob(atoi(token), k,
						    alteration, state_probs,
						    prob_row_seq[file_row]);
      }
     }
    file_row++;
  }

  // FIXME: remove this from production code??
  if (file_row != num_sequences) {
    if(file_row > num_sequences) {
      Rprintf("ERROR in read_convert_prob_seq: file_row > num_sequences");
      error("\nFatal error in RJaCGH. Please let us know.\n");
    } else {
      Rprintf("ERROR in read_convert_prob_seq: file_row < num_sequences");
      error("\nFatal error in RJaCGH. Please let us know.\n");
    }
  }
  gzclose(infile);
}


void print_regions_pREC_A(const int numregions,
			  const int *regionsStart,
			  const int *regionsEnd,
			  const double *regionsProb) {
  Rprintf("\n    --- Regions found --- \n");
  for(int i = 0; i < numregions; i++) {
    Rprintf("\n Region: %i. Start = %i, End = %i, Prob = %f",
	   i, regionsStart[i] + 1, regionsEnd[i] + 1, 
	   regionsProb[i]);
  }
  Rprintf("\n\n");
  
}
		  


void print_stretched(const int num_probes,
		     const int total_num_sequences,
		     const int num_arrays,
		     const int *indices_sequences,
		     const double *prob_row_seq,
		     double **stretched) {
  // BEWARE: entries are now multiplied by prob_row_seq!!!
  int arrindex = 1;
  Rprintf("\n  ----  Stretched sequence (probabilities) --- ");
  for(int i = 0; i < total_num_sequences; i++) {
    if(i >= indices_sequences[arrindex]) arrindex++;
    if(arrindex > num_arrays) {
      Rprintf("\nprint_stretched: We've got a problem: arrindex > num_arrays.\n");
      error("\nFatal error in RJaCGH. Please let us know.\n");
    }
    Rprintf("\n Array= %i, Row= %i, Prob= %f || ", arrindex, i, 
	    prob_row_seq[i]);
    for(int j = 0; j < num_probes; j++) {
      Rprintf("%f ", stretched[i][j]);
    }
  }
  Rprintf("\n");
}
  
double prob_seq_array(const double *prob_alteration_seq,
		      const int total_num_sequences) {
  // prob_alteration_seq: vector where each row is the probability
  // of a given alteration in a sequence of probes.

  int i = 0;
  double tmp_res = 0.0;

  for(i = 0; i < total_num_sequences; i++) {
    tmp_res += prob_alteration_seq[i];
  }
  // FIXME: remove when well tested; should be caught before.
  CHECK_NUM(tmp_res, prob_seq_array);
  return tmp_res;
} 

void update_prob_alteration_seq(const double *probe_probs_vector,
				const int total_num_sequences,
				double *prob_alteration_seq) {
  //BEWARE: this only updates one probe at a time!!!
  for(int nr = 0; nr < total_num_sequences; nr++) {
    prob_alteration_seq[nr] = fmin(prob_alteration_seq[nr], 
				   probe_probs_vector[nr]);	   
  }
}


double prob_seq_all(const double *prob_alteration_seq,
		    const double *array_weights,
		    const int *indices_sequences,
		    const int total_num_sequences,
		    const int num_arrays) {
/* 		    const int printout) { */
  // prob_alteration_seq: vector where each row is the probability
  // of a given alteration in a sequence of probes.
  // prob_row_seq: the probability of that "row" (viterbi * k).

  int arrindex = 1;
  int i = 0;
  double res = 0.0;
  double tmp_res = 0.0;
  if(array_weights[0] < 0) { //No array weights
    for(i = 0; i < total_num_sequences; i++) {
      tmp_res += prob_alteration_seq[i]; // * prob_row_seq[i];
      if((i + 1) >= indices_sequences[arrindex]) {
	res += tmp_res;
	tmp_res = 0.0;
	arrindex++;
      }
    }
    res /= num_arrays;
  } else { // using array weights 
    for(i = 0; i < total_num_sequences; i++) {
      tmp_res += prob_alteration_seq[i]; // * prob_row_seq[i];
      if((i + 1) >= indices_sequences[arrindex]) {
	res += tmp_res * array_weights[(arrindex - 1)];
	tmp_res = 0.0;
	arrindex++;
      }
    }
  }

  // FIXME: remove later, when well tested?
  if(arrindex != (num_arrays + 1)) {
    if(arrindex > (num_arrays + 1)) {
      Rprintf("\n prob_seq_all: We've got a problem: arrindex > num_arrays. \n");
      error("\nFatal error in RJaCGH. Please let us know.\n");
    } else {
      Rprintf("\n prob_seq_all: We've got a problem: arrindex < num_arrays. \n");
      error("\nFatal error in RJaCGH. Please let us know.\n");
    }
  }
  // FIXME: remove later
  CHECK_NUM(res, prob_seq_all);
  return res;
} 



void pREC_A(const int num_arrays, const int num_probes,
	    const int total_num_sequences,
	    const double threshold,
	    const double *array_weights,
	    double **stretched,
	    const int *starting_indices_sequences,
	    int *numregions,
	    int *regionsStart, int *regionsEnd, double *regionsProb) {

  int start;
  int end;
  double tp_1 = 0.0;
  double tp_2 = 0.0;
  
  // To allow fast sequential access by probe (column in original stretched)
  // we transpose (so stretched_tr is probe by sequence).
  double **stretched_tr; stretched_tr = Calloc(num_probes, double*);
  for(int ii = 0; ii < num_probes; ii++) {
    stretched_tr[ii] = Calloc(total_num_sequences, double);
  }
  for(int row = 0; row < total_num_sequences; row++) {
    for(int column = 0; column < num_probes; column++) {
      stretched_tr[column][row] = stretched[row][column];
    }
  }

  *numregions = 0;
  start = 0;
  while(start < num_probes) {
    tp_1 = prob_seq_all(stretched_tr[start], array_weights,
			  starting_indices_sequences, total_num_sequences,
			  num_arrays);
    if(tp_1 >= threshold) {
      end = start + 1;
      while(end < num_probes) {
	update_prob_alteration_seq(stretched_tr[end], total_num_sequences, 
				   stretched_tr[start]);
	tp_2 = prob_seq_all(stretched_tr[start],
			    array_weights, 
			    starting_indices_sequences,
			    total_num_sequences, num_arrays);
	if(tp_2 < threshold) {
	  break;
	} else {
	  tp_1 = tp_2;
	  end++;
	}
      }
      // Add this region; we fall-through to this place either when we can no longer
      // add probes (tp_2 < threshold) or when we are at the end of the set of probes
      // (end = num_probes, so we are out of the while loop).
      regionsStart[*numregions] = start;
      regionsEnd[*numregions]   = end - 1;
      regionsProb[*numregions] = tp_1;
      (*numregions)++;
      start = end; //update HERE: we only get here if tp_2 < threshold or if
      // no more probes;
    } else {
      start++;
    }
  }
  Free_2(stretched_tr, num_probes);
}

void print_regions_pREC_S(struct regionS *seq) {
  int num_seqs = 0;
  int i;
  Rprintf("\n Initial pointer to pREC_S: %p\n", seq);
  Rprintf("\nPointer\tRegion\tStart\tEnd\tProbes\tArrays\n");
  while(seq != NULL) {
    num_seqs++;
    Rprintf("%p\t%i\t%i\t%i\t%i\t", seq, num_seqs, (seq->Start) + 1, 
	   (seq->End) + 1,
	   (seq->End) - (seq->Start) + 1);
    for(i = 0; i < (seq->num_arrays); i++) {
      Rprintf("%i;", (seq->arrays[i]) + 1);
    }
    Rprintf("\n");
    seq = seq->next;
  }
}
  
// A R devolver: Start, End, vector con numero de arrays
// en cada secuencias y 
// todo junto (colapsado) un vector con los indices de arrays.
// al leer en R, corto por los índices

void pREC_S(const int num_arrays, const int num_probes,
	    const int total_num_sequences,
	    const double p_w,
	    const int freq_arrays,
	    double **stretched,
	    const int *starting_indices_sequences) {
  int probe, arr, next_probe;
  int offset, arr_n_seqs;
  int candidate0[num_arrays];
  int candidate1[num_arrays];
  int ok_a0, ok_a1;
  double *marginal_probs; 
  marginal_probs = (double *) R_alloc(num_probes, sizeof(double));
  double *current_prob_alteration_seq;
  current_prob_alteration_seq = (double *) R_alloc(total_num_sequences,
						   sizeof(double));

  // To allow fast sequential access by probe (column in original stretched)
  // we transpose: stretched_tr is probe by sequence. 
  double **stretched_tr; stretched_tr = Calloc(num_probes, double*);
  for(int ii = 0; ii < num_probes; ii++) {
    stretched_tr[ii] = Calloc(total_num_sequences, double);
  }
  for(int row = 0; row < total_num_sequences; row++) {
    for(int column = 0; column < num_probes; column++) {
      stretched_tr[column][row] = stretched[row][column];
    }
  }

  for(probe = 0; probe < num_probes; probe++) {
    // Which arrays are to be considered?
    ok_a0 = 0;
    for(arr = 0; arr < num_arrays; arr++) {
      offset = starting_indices_sequences[arr];
      arr_n_seqs = starting_indices_sequences[arr + 1] - offset;
      if(prob_seq_array(&stretched_tr[probe][offset],
			arr_n_seqs) >= p_w) {
	candidate0[ok_a0] = arr;
	ok_a0++;
      }
    }
    // if(ok_a0 < freq_arrays) continue; // obscure
    if(ok_a0 >= freq_arrays) { // add probes
      next_probe = probe + 1;
      while(next_probe < num_probes) {
	ok_a1 = 0;
	for(int ia = 0; ia < ok_a0; ia++) {
	  arr = candidate0[ia];
	  offset = starting_indices_sequences[arr];
	  arr_n_seqs = starting_indices_sequences[arr + 1] - offset;
	  // Store updates in stretched_tr[probe]; "probe" is first probe of a 
	  // sequence, so never revisited in subsequent iterations.
	  update_prob_alteration_seq(&stretched_tr[next_probe][offset],
				     arr_n_seqs, &stretched_tr[probe][offset]);
	  if(prob_seq_array(&stretched_tr[probe][offset], arr_n_seqs) >= p_w) {
	    candidate1[ok_a1] = arr;
	    ok_a1++;
	  } 
	}
	if(ok_a1 < freq_arrays) break; // no more looping over probes; 
                                       //fall through and write ok_a0 arrays 
	else {
	  if(ok_a1 < ok_a0) {
	    // Fewer arrays meet pw in this step: acceptable region has decreased.
	    // Thus, previous step marks the end of a common region.
	    updateRegionS(&regS, probe, next_probe - 1, ok_a0, candidate0);
	    ok_a0 = ok_a1;
	    for(int aa = 0; aa < ok_a1; aa++) 
	      candidate0[aa] = candidate1[aa];
	  }
	  next_probe++;
	}
      }
      updateRegionS(&regS, probe, next_probe - 1, ok_a0,
			candidate0);
    } // ok_a0 > freq_arrays
  } // end loop over probes
  Free_2(stretched_tr, num_probes);
}

void wrap_pREC(const int *alteration,
	       const int *numarrays,
	       const int *num_sequences,
	       const int *num_probes,
	       const int *starting_indices_sequences,
	       const int *starting_indices_state_probs,
	       char **filename,
	       const double *threshold,
	       const int *freq_arrays, //only for pREC-S
	       const double *array_weights, //only for pREC-A
	       const double *state_probs,
	       int *numregions,
	       int *total_narrays, //only for pREC-S
	       int *regionsStart,  //here only for pREC-A
	       int *regionsEnd,    //ditto
	       double *regionsProb, //ditto
	       const int *verboseC, 
	       const int *method_prec) {
  
  // Common wrapper to pREC_A and pREC_S; most of dealing with
  // files and obtaining the viterbi sequences as probabilities.

  
  int total_num_sequences = starting_indices_sequences[(*numarrays)];
  int index_seq, index_state_probs;
  double *prob_row_seq;
  double **stretched;
  
  prob_row_seq = (double *) R_alloc(total_num_sequences, sizeof(double));
  stretched = (double **) Calloc(total_num_sequences, double*);
  for(int ii = 0; ii < total_num_sequences; ii++) {
    stretched[ii] = Calloc(*num_probes, double);
  }

  char *filenames[(*numarrays)];
  filenames[0] = strtok(*filename, "\n");
  for(int narray = 1; narray < (*numarrays); narray++) {
    filenames[narray] = strtok(NULL, "\n");
  }

  for(int i = 0; i < (*numarrays); i++) {
    index_seq = starting_indices_sequences[i];
    index_state_probs = starting_indices_state_probs[i];
    read_convert_prob_seq(prob_row_seq + index_seq,
			  stretched + index_seq,
			  filenames[i],
			  state_probs + index_state_probs,
			  *alteration,  num_sequences[i]);
  }

  if((*method_prec) == 0) {
    pREC_A(*numarrays, *num_probes, total_num_sequences, *threshold,
	   array_weights, stretched,
	   starting_indices_sequences, numregions,
	   regionsStart, regionsEnd, regionsProb);
    if(*verboseC) { 
      print_stretched(*num_probes, total_num_sequences, *numarrays,
		      starting_indices_sequences, prob_row_seq, stretched);
      print_regions_pREC_A(*numregions, regionsStart, regionsEnd, 
			   regionsProb);
    }
  } else if((*method_prec == 1)) {
    pREC_S(*numarrays, *num_probes, total_num_sequences, *threshold,
	   *freq_arrays, stretched,
	   starting_indices_sequences);
    if(regS != NULL) {
      //Return to R the size of the vectors it will need for second call
      *numregions = regS->region_number;
      *total_narrays = regS->sum_num_arrays;
    } else {
      *numregions = 0;
      *total_narrays = 0;
    }
    if(*verboseC) {
      print_regions_pREC_S(regS);
    }

  } else {
    error("\n method_prec has to be one of 0 --pREC_A--- or 1 ---pREC_S---\n");
  }
  Free_2(stretched, total_num_sequences);
}


double dnorm5(double x, const double mu, const double sigma)
{ /* From dnorm4 in R */
    x = (x - mu) / sigma;
    return (M_1_SQRT_2PI * exp(-0.5 * x * x)  /   sigma);
}


double dnorm6(double x, const double mu, const double sigma)
{ /* From dnorm4 in R; for logged density */
    x = (x - mu) / sigma;
    return (-(M_LN_SQRT_2PI  +  0.5 * x * x + log(sigma)));
}


void normalNHHMMlikelihood(const double *y, const int *k, const double *x, 
			   const int *n, const double *q,  
			   const double *stat, const double *mu,
			   const double *sigma2, 
			   double *loglik) {
#ifdef DEBUG
  Rprintf("\n Likelihood evaluation \n\n"); fflush(stdout);
#endif

  double res=0;
  int hidden_state = 0;
  
  if ((*k)==1) { /*  with one state is just log lik. */
    double sqrtsigma = sqrt(*sigma2);
    for (int ii = 0; ii < (*n); ++ii) {
      res += dnorm6(y[ii], *mu, sqrtsigma);
    }
  }
  else { /*  >1 states */
    int k2 = (*k) * (*k);
    int k_mult;
    int Q_index;
    double c_i, y_i;
    double one_minus_x;
    double one_over_c_i;
    double *filter; filter = Calloc((*k), double); 
    double *filtercond; filtercond = Calloc((*k), double);
    double *sd_j; sd_j = Calloc((*k), double);
    double *rowSumsQ; rowSumsQ = Calloc((*k), double);
    double *Q; Q = Calloc(k2, double);
    double *q_prime; q_prime = Calloc(k2, double);

    // transpose q just once so elements sequential and
    // same order as Q.
    Q_index = 0;
    int row;
    int col;
    for(row = 0; row < (*k); row++) {
	k_mult = 0;
      for(col = 0; col < (*k); col++) {
	q_prime[Q_index] = q[row + k_mult];
	k_mult += (*k);
	Q_index++;
      }
    }
    //initialize filtercond and (logically unrelated) precompute standard dev.
    for (hidden_state=0; hidden_state < (*k); hidden_state++) {
      filtercond[hidden_state] = stat[hidden_state];
      sd_j[hidden_state] = sqrt(sigma2[hidden_state]);
    }

    // loop over all cases except last
    for (int obs_index = 0; obs_index < ((*n) - 1); obs_index++) {
      c_i = 0.0;
      y_i = y[obs_index];

      for(int j = 0; j < (*k); j++) {
	filter[j] = filtercond[j] * dnorm5(y_i, mu[j], sd_j[j]);
 	c_i += filter[j]; 
      }
      if(c_i <= 10E-300){
	c_i = 10E-300;
      }

      //do sum of logs, not log of products, for numerical stability
      res += log(c_i);

      one_minus_x = 1. - x[obs_index];
      one_over_c_i = 1. / c_i;

/*       CHECK_NUM(c_i, B); */
/*       CHECK_NUM(one_over_c_i, B); */

      Q_index = 0;
      for(int row = 0; row < (*k); row++) {
	rowSumsQ[row] = 0.0;
	for(int column = 0; column < (*k); column++) {
	  if(row != column) Q[Q_index] = exp(q_prime[Q_index] * one_minus_x);
	  else Q[Q_index] = 1.0;
	  rowSumsQ[row] += Q[Q_index];
	  Q_index++;
	}
	filter[row] /= rowSumsQ[row];
      }

      for(int row = 0; row < (*k); row++) {
	k_mult = row;
	filtercond[row] = filter[0] * Q[k_mult];
	for(int col = 1; col < (*k); col++) {
	  k_mult += (*k);
	  filtercond[row] += filter[col] * Q[k_mult];
	}
	filtercond[row] *= one_over_c_i;
      }	  
    }

    // last case: only filter computation.
    y_i = y[(*n) - 1];
    c_i = 0.0;
    for(int j = 0; j < (*k); j++) {
      filter[j] = filtercond[j] * dnorm5(y_i, mu[j], sd_j[j]);
      c_i += filter[j];
    }
    if(c_i <= 10E-300) {
      c_i = 10E-300;
    }

    res += log(c_i);
    Free(filter);
    Free(filtercond);
    Free(sd_j); 
    Free(rowSumsQ);
    Free(q_prime); 
    Free(Q);
  }
  CHECK_NUM(res, normalNHHMlikelihood);
  *loglik = res;
}
  

void Birth(double *y, double *x, int *varEqual, int *genome, int *index, 
	   double *mu, double *sigma2, 
	   double *beta, double *stat,
	   double *statBirth, 
	   int *r, double *loglikLast, double *probK, double *pb, double *muAlfa,
	   double *muBeta, int *n, double *candidatoMu,
	   double *candidatoSigma2,  
	   double *candidatoBeta, double *loglikBirth,
	   int *accepted, double *maxVar, double *s1, double *s2, double heat) {

  /* What factors are affected by heat? */
  /* By now, only the likelihood (as suggested by Green) */

  /* parameters for different gamma distributions for betas */
  double gam1=1; double gam2=1;
  /**********************************************************/
  double *candidatoQ; candidatoQ = Calloc((*r+1)*(*r+1), double);
  double *candidatoQ2; candidatoQ2 = Calloc((*r+1)*(*r+1), double);
  double *candidatoBeta2; candidatoBeta2 = Calloc((*r+1)*(*r+1), double);
    double *candidatoMu2; candidatoMu2 = Calloc(*r+1, double);
    double *candidatoSigma22; candidatoSigma22 = Calloc(*r+1, double);
    double probBirth=0; double probBirth2=0;  double probBirth3=0;
    int newState;
    double loglikCandidate = 0;
    double loglikCandidate2 = 0;
    double loglikPartial = 0;
    int i, j, k;
    int nn;

    /* Delayed rejection */
    /* First stage x -> y */

    for (i=0; i< *r; ++i) {
      candidatoMu[i] = mu[i];
      candidatoSigma2[i] = sigma2[i];
    }
    candidatoMu[*r] = rnorm(*muAlfa, *s1);
    if (*varEqual) {
      candidatoSigma2[*r] = candidatoSigma2[0];
    }
    else {
      candidatoSigma2[*r] = runif(0, *maxVar);
      candidatoSigma2[*r] = pow(candidatoSigma2[*r], 2);
    }
    for (i=0; i<*r; ++i) {
      for (j=0; j<*r; ++j) {
	candidatoBeta[i*(*r+1)+j] = beta[i* *r +j];
	candidatoQ[i*(*r+1)+j] = -candidatoBeta[i*(*r+1)+j];
      }
      candidatoBeta[i* (*r+1) + *r] = rgamma(gam1,gam2);
      candidatoQ[i* (*r+1) + *r]  = - candidatoBeta[i * (*r+1) + *r];
      probBirth  = probBirth + dgamma(candidatoBeta[i * (*r+1) + *r], 1, 1, 1);
      probBirth  = probBirth - dgamma(candidatoBeta[i * (*r+1) + *r], gam1, gam2, 1);
    }
    for (i=0; i<*r; ++i) {
      candidatoBeta[*r * (*r+1) + i] = rgamma(gam1,gam2);
      candidatoQ[*r * (*r+1)  + i] = - candidatoBeta[*r * (*r+1) + i];
      probBirth = probBirth + dgamma(candidatoBeta[*r * (*r+1) + i], 1, 1, 1);
      probBirth = probBirth - dgamma(candidatoBeta[*r * (*r+1) + i], gam1, gam2, 1);
    }
    candidatoBeta[(*r+1) * (*r+1) -1] = 0;
    candidatoQ[(*r+1) * (*r+1) -1] = 0;

    newState = *r + 1;
    
    probBirth = probBirth + log(probK[*r]) - log(probK[*r-1]);
    probBirth = probBirth + dnorm(candidatoMu[*r], *muAlfa, *muBeta, 1); 
    probBirth = probBirth - dnorm(candidatoMu[*r], *muAlfa, *s1, 1); 
    
    for (k=0; k < *genome; ++k) {
      nn = index[k+1] - index[k];
      //P1(1098);
      double *yy; yy = Calloc(nn, double);
      double *xx; xx = Calloc(nn, double);
      //P2(1098);
     
      for (j=0; j < nn-1; ++j) {
	yy[j] = y[j + index[k]];
	xx[j] = x[j + index[k]];
      }
      xx[nn-1] = 0;
      yy[nn-1] = y[nn-1 + index[k]];
      normalNHHMMlikelihood(yy, &newState, xx, &nn, candidatoQ,  
			  statBirth, candidatoMu, candidatoSigma2, 
			  &loglikPartial);
      loglikCandidate += loglikPartial;
      Free(xx);
      Free(yy);

    }
    probBirth = probBirth + heat * (loglikCandidate - *loglikLast);
    probBirth = probBirth + log(1-pb[*r]) - log(pb[*r-1]);
    if(probBirth >= 0) probBirth = 1; else probBirth = exp(probBirth);
    CHECK_NUM(probBirth, Birth);
/*     if (probBirth >1) probBirth = 1; */
    if (runif(0,1) <= probBirth) {
      *accepted =1;
      *loglikBirth = loglikCandidate;
    }
    /* Second stage */
    else {
      for (i=0; i< *r; ++i) {
	candidatoMu2[i] = mu[i];
	candidatoSigma22[i] = sigma2[i];
      }
      candidatoMu2[*r] = rnorm(*muAlfa, *s2);
      if (*varEqual) {
	candidatoSigma22[*r] = candidatoSigma22[0];
      }
      else {
	candidatoSigma22[*r] = runif(0, *maxVar);
	candidatoSigma22[*r] = pow(candidatoSigma22[*r], 2);
      }
      for (i=0; i<*r; ++i) {
	for (j=0; j<*r; ++j) {
	  candidatoBeta2[i*(*r+1)+j] = beta[i* *r +j];
	  candidatoQ2[i*(*r+1)+j] = -candidatoBeta2[i*(*r+1)+j];
	}
	candidatoBeta2[i* (*r+1) + *r] = rgamma(gam1,gam2);
	candidatoQ2[i* (*r+1) + *r]  = - candidatoBeta2[i * (*r+1) + *r];
	probBirth2 = probBirth2 + dgamma(candidatoBeta2[i* (*r+1) + *r], gam1, gam2, 1);
	probBirth2 = probBirth2 - dgamma(candidatoBeta2[i* (*r+1) + *r], 1, 1, 1);
	probBirth3 = probBirth3 + dgamma(candidatoBeta2[i* (*r+1) + *r], 1, 1, 1);
	probBirth3 = probBirth3 - dgamma(candidatoBeta2[i* (*r+1) + *r], gam1, gam2, 1);

      }
      for (i=0; i<*r; ++i) {
	candidatoBeta2[*r * (*r+1) + i] = rgamma(1,1);
	candidatoQ2[*r * (*r+1)  + i] = - candidatoBeta2[*r * (*r+1) + i];
	probBirth2 = probBirth2 + dgamma(candidatoBeta2[*r * (*r+1) + i], gam1, gam2, 1);
	probBirth2 = probBirth2 - dgamma(candidatoBeta2[*r * (*r+1) + i], 1, 1, 1);
	probBirth3 = probBirth3 + dgamma(candidatoBeta2[*r * (*r+1) + i], 1, 1, 1);
	probBirth3 = probBirth3 - dgamma(candidatoBeta2[*r * (*r+1) + i], gam1, gam2, 1);
      }
      candidatoBeta2[(*r+1) * (*r+1) -1] = 0;
      candidatoQ2[(*r+1) * (*r+1) -1] = 0;
      
      /* Virtual move z -> y* */
      probBirth2 = probBirth2 + log(probK[*r-1]) - log(probK[*r]);
      probBirth2 = probBirth2 + dnorm(candidatoMu2[*r], *muAlfa, *s1, 1);
      probBirth2 = probBirth2 - dnorm(candidatoMu2[*r], *muAlfa, *muBeta, 1);
      for (k=0; k < *genome; ++k) {
	nn = index[k+1] - index[k];
	//P1(1170);
	double *yy; yy = Calloc(nn, double);
	double *xx; xx = Calloc(nn, double);
	//P2(1170);
	
	for (j=0; j < nn-1; ++j) {
	  yy[j] = y[j + index[k]];
	  xx[j] = x[j + index[k]];
	}
	xx[nn-1] = 0;
	yy[nn-1] = y[nn-1 + index[k]];
	normalNHHMMlikelihood(yy, &newState, xx, &nn, candidatoQ2, 
			      statBirth, candidatoMu2, candidatoSigma22, 
			      &loglikPartial);
	loglikCandidate2 += loglikPartial;
	Free(xx);
	Free(yy);
    }
      probBirth2 = probBirth2 + heat * (*loglikLast - loglikCandidate2);
      probBirth2 = probBirth2 + log(pb[*r-1]) - log(1-pb[*r]);
      if(probBirth2 >= 0) probBirth2 = 1; else probBirth2 = exp(probBirth2);
      CHECK_NUM(probBirth2, Birth);
/*       if (probBirth2 >1) probBirth2 = 1;     */
      /* x -> z */
      probBirth3 = probBirth3 + log(probK[*r]) - log(probK[*r-1]);
      probBirth3 = probBirth3 + dnorm(candidatoMu2[*r], *muAlfa, *muBeta, 1);
      probBirth3 = probBirth3 + dnorm(candidatoMu[*r], *muAlfa, *s2, 1);
      probBirth3 = probBirth3 - dnorm(candidatoMu[*r], *muAlfa, *s1, 1);
      probBirth3 = probBirth3 - dnorm(candidatoMu2[*r], *muAlfa, *s2, 1);
      probBirth3 = probBirth3 + heat * (loglikCandidate2 - *loglikLast) ;
      probBirth3 = probBirth3 - log(*r);
      probBirth3 = probBirth3 + log(1-pb[*r]) + log(1-pb[*r]);
      probBirth3 = probBirth3 - log(pb[*r-1]) - log(pb[*r-1]);
      probBirth3 = probBirth3 + log(1-probBirth2) - log(1-probBirth);
      if(probBirth3 >= 0) probBirth3 = 1; else probBirth3 = exp(probBirth3);
      CHECK_NUM(probBirth3, Birth);
/*       if (probBirth3 >1) probBirth3 = 1; */
      if (runif(0,1) <= probBirth3) {
	*accepted =2;
	*loglikBirth = loglikCandidate2;
	for (i=0; i < *r+1; i++) {
	  candidatoMu[i] = candidatoMu2[i];
	  candidatoSigma2[i] = candidatoSigma22[i];
	}
	for (i=0; i < (*r+1) * (*r+1); i++) {
	  candidatoBeta[i] = candidatoBeta2[i];
	}
      }
    }
/*     printf("logBirth=%f\n", *loglikBirth); */
    Free(candidatoSigma22);    
    Free(candidatoMu2);  
    Free(candidatoBeta2);
    Free(candidatoQ2);
    Free(candidatoQ);
  }

void Death(double *y, double *x, int *genome, int *index, double *mu, double *sigma2, 
	   double *beta, double *stat,
	   double *statDeath, double *muAlfa, double *muBeta, 
	   int *r, double *loglikLast, double *probK, double *pb, 
	   int *n, double *candidatoMu,
	   double *candidatoSigma2, 
	   double *candidatoBeta, double *loglikDeath,
	   int *accepted, double *s1, double *s2, double heat) {

  /* What factors are affected by heat? */
  /* By now, only priors and likelihoods */

  /* parameters for different betas */
  double gam1=1; double gam2=1;
  /**********************************/

  double *candidatoQ; candidatoQ = Calloc((*r-1)*(*r-1), double);
  double probDeath=0; double probDeath2=0; double probDeath3=0;
  int death, death2, indexBeta=0, indexMu=0;
  int newState;
  double loglikCandidate = 0;
  double loglikCandidate2 = 0;
  double loglikPartial = 0;
  int nn;
  int i, j, k;
  double u1bp, u2b;

  /* think carefully about these values */
/*   double s1 = *muBeta * 5; double s2 = *muBeta/ 1; */

  /* Delayed rejection */
  /* First stage z -> y* */

  death = (int)rint(runif(1, *r));
  u1bp = mu[death-1];
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
    for (j=0; j < (death-1); ++j) {
      probDeath = probDeath + dgamma(beta[j * *r + (death -1)], gam1, gam2, 1);
      probDeath = probDeath - dgamma(beta[j * *r + (death -1)], 1, 1, 1);
      probDeath = probDeath + dgamma(beta[(death-1) * *r + j], gam1, gam2, 1);
      probDeath = probDeath - dgamma(beta[(death -1) * *r +j], 1, 1, 1);
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
    for (j=death; j < *r; ++j) {
      probDeath = probDeath + dgamma(beta[j * *r + (death -1)], gam1, gam2, 1);
      probDeath = probDeath - dgamma(beta[j * *r + (death -1)], 1, 1, 1);
      probDeath = probDeath + dgamma(beta[(death -1) * *r + j], gam1, gam2, 1);
      probDeath = probDeath - dgamma(beta[(death -1) * *r + j], 1, 1, 1);
    }
  }
  newState = *r -1;
  probDeath = probDeath + log(probK[*r-2]) - log(probK[*r-1]);
  for (k=0; k < *genome; ++k) {
    nn = index[k+1] - index[k];
    //P1(1316);
    double *yy; yy = Calloc(nn, double);
    double *xx; xx = Calloc(nn, double);
    //P2(1316);
    for (j=0; j < nn-1; ++j) {
      yy[j] = y[j + index[k]];
      xx[j] = x[j + index[k]];
    }
    xx[nn-1] = 0;
    yy[nn-1] = y[nn-1 + index[k]];
    normalNHHMMlikelihood(yy, &newState, xx, &nn, candidatoQ,  
			  statDeath, candidatoMu, candidatoSigma2, 
			  &loglikPartial);
    loglikCandidate += loglikPartial;
    Free(xx);
    Free(yy);
  }
  probDeath = probDeath + dnorm(u1bp, *muAlfa, *s1, 1);
  probDeath = probDeath - dnorm(u1bp, *muAlfa, *muBeta, 1);
  probDeath = probDeath + heat * (loglikCandidate - *loglikLast);
  probDeath = probDeath + log(pb[*r-2]) - log(1-pb[*r-1]);
  if(probDeath >= 0) probDeath = 1; else probDeath = exp(probDeath);
  CHECK_NUM(probDeath, Death);
/*   if (probDeath >1) probDeath = 1; */
  if (runif(0,1) <= probDeath) {
    *accepted =1;
    *loglikDeath = loglikCandidate;
  }
  /* Second stage */
  else {
    death2 = death;
    while (death2 == death) {
      death2 = (int)rint(runif(1, *r));
    }
    indexMu=0; indexBeta=0;
    u2b = mu[death2-1];
    if (death2 >1) {
      for (i=0; i< (death2-1); ++i) {
	candidatoMu[indexMu] = mu[i];
	candidatoSigma2[indexMu] = sigma2[i];
	indexMu ++;
	for (j=0; j< (death2-1); ++j) {
	  candidatoBeta[indexBeta] = beta[i* *r + j];
	  candidatoQ[indexBeta] = -candidatoBeta[indexBeta];
	  indexBeta ++;
	}
	if (death2 < *r) {
	  for (j=death2; j< *r; ++j) {
	    candidatoBeta[indexBeta] = beta[i* *r + j];
	    candidatoQ[indexBeta] = -candidatoBeta[indexBeta];
	    indexBeta ++;
	  }
	}
      }

      for (j=0; j < (death2-1); ++j) {
	probDeath2 = probDeath2 + dgamma(beta[j * *r + (death2 -1)], 1, 1, 1);
	probDeath2 = probDeath2 - dgamma(beta[j * *r + (death2 -1)], gam1, gam2, 1);
	probDeath2 = probDeath2 + dgamma(beta[(death2 -1) * *r + j], 1, 1, 1);
	probDeath2 = probDeath2 - dgamma(beta[(death2 -1) * *r +j], gam1, gam2, 1);
	probDeath3 = probDeath3 + dgamma(beta[j * *r + (death2 -1)], gam1, gam2, 1);
	probDeath3 = probDeath3 - dgamma(beta[j * *r + (death2 -1)], 1, 1, 1);
	probDeath3 = probDeath3 + dgamma(beta[(death2 -1) * *r + j], gam1, gam2, 1);
	probDeath3 = probDeath3 - dgamma(beta[(death2 -1) * *r +j], 1, 1, 1);
      }
    }
    if (death2 < *r) {
      for (i=death2; i< *r; ++i) {
	candidatoMu[indexMu] = mu[i];
	candidatoSigma2[indexMu] = sigma2[i];
	indexMu ++;
	
	if (death2 > 1) {
	  for (j=0; j< (death2-1); ++j) {
	    candidatoBeta[indexBeta] = beta[i* *r + j];
	    candidatoQ[indexBeta] = -candidatoBeta[indexBeta];
	    indexBeta ++;
	  }
	}
	for (j=death2; j< *r; ++j) {
	  candidatoBeta[indexBeta] = beta[i* *r + j];
	  candidatoQ[indexBeta] = -candidatoBeta[indexBeta];
	  indexBeta++;
	}
      }
      for (j=death2; j < *r; ++j) {
	probDeath2 = probDeath2 + dgamma(beta[j * *r + (death2 -1)], 1, 1, 1);
	probDeath2 = probDeath2 - dgamma(beta[j * *r + (death2 -1)], gam1, gam2, 1);
	probDeath2 = probDeath2 + dgamma(beta[(death2 -1) * *r + j], 1, 1, 1);
	probDeath2 = probDeath2 - dgamma(beta[(death2 -1) * *r + j], gam1, gam2, 1);
	probDeath3 = probDeath3 + dgamma(beta[j * *r + (death2 -1)], gam1, gam2, 1);
	probDeath3 = probDeath3 - dgamma(beta[j * *r + (death2 -1)], 1, 1, 1);
	probDeath3 = probDeath3 + dgamma(beta[(death2 -1) * *r + j], gam1, gam2, 1);
	probDeath3 = probDeath3 - dgamma(beta[(death2 -1) * *r + j], 1, 1, 1);
      }
    }

    /* Virtual move x -> y */
    probDeath2 = probDeath2 + log(probK[*r-2]) - log(probK[*r-1]);
    probDeath2 = probDeath2 + dnorm(u2b, *muAlfa, *muBeta, 1);
    probDeath2 = probDeath2 - dnorm(u2b, *muAlfa, *s1, 1);
    for (k=0; k < *genome; ++k) {
      nn = index[k+1] - index[k];
      //P1(1417);
      double *yy; yy = Calloc(nn, double);
      double *xx; xx = Calloc(nn, double);
      //P2(1417);
      for (j=0; j < nn-1; ++j) {
	yy[j] = y[j + index[k]];
	xx[j] = x[j + index[k]];
      }
      xx[nn-1] = 0;
      yy[nn-1] = y[nn-1 + index[k]];
      normalNHHMMlikelihood(yy, &newState, xx, &nn, candidatoQ,  
			    statDeath, candidatoMu, candidatoSigma2, 
			    &loglikPartial);
      loglikCandidate2 += loglikPartial;
      Free(xx);
      Free(yy);
    }
    probDeath2 = probDeath2 + heat * (*loglikLast - loglikCandidate2); 
    probDeath2 = probDeath2 + log(1.0 - pb[*r-1]) - log(pb[*r-2]);
    if(probDeath2 >= 0) probDeath2 = 1; else probDeath2 = exp(probDeath2);
    CHECK_NUM(probDeath2, Death);
/*     if (probDeath2 > 1) probDeath2 = 1; */
    /* z -> x */
    probDeath3 = probDeath3 + log(probK[*r-2]) - log(probK[*r-1]);
    probDeath3 = probDeath3 + dnorm(mu[death-1], *muAlfa, *s1, 1);
    probDeath3 = probDeath3 + dnorm(mu[death2-1], *muAlfa, *s2, 1);
    probDeath3 = probDeath3 - dnorm(mu[death-1], *muAlfa, *s2, 1);
    probDeath3 = probDeath3 - dnorm(mu[death2-1], *muAlfa, *muBeta, 1);
    probDeath3 = probDeath3 + heat * (loglikCandidate2 - *loglikLast); 
    probDeath3 = probDeath3 + log(*r-1);
    probDeath3 = probDeath3 + log(pb[*r-2]) + log(pb[*r-2]);
    probDeath3 = probDeath3 - log(1-pb[*r-1]) - log(1-pb[*r-1]);
    probDeath3 = probDeath3 + log(1-probDeath2) - log(1-probDeath);
    if(probDeath3 >= 0) probDeath3 = 1; else probDeath3 = exp(probDeath3);
    CHECK_NUM(probDeath3, Death);
/*     if (probDeath3 > 1) probDeath3 = 1; */
    if (runif(0,1) <= probDeath3) {
      *accepted = 2;
      *loglikDeath = loglikCandidate2;
    }
  }
  Free(candidatoQ);
}
/*  this model splits somewhat different the beta matrix */
/* tauSplitBeta is now the constant 1/(1-x) */
void Split(double *y, double *x, int *varEqual, int *genome, int *index, 
	   double *mu, double *sigma2, double *beta, 
	   double *stat, double *statSplit, 
	   int *r, double *loglikLast, double *probK, double *ps, 
	   int *n, double *candidatoMu,
	   double *candidatoSigma2, 
	   double *candidatoBeta, double *loglikSplit,
	   double *muAlfa, double *muBeta,
	   double *tauSplitMu,
	   double *tauSplitBeta, int *accepted, 
	   double *maxVar, double heat) {


#ifdef DEBUG
  Rprintf("\n Entering split \n"); fflush(stdout);
#endif

  double *ui; ui = Calloc(*r-1, double);
  double *epj; epj = Calloc(*r-1, double);
  double *gi0; gi0 = Calloc(2, double);
  double *candidatoQ; candidatoQ = Calloc((*r+1)*(*r+1), double);

  double probSplit;
  int split;
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
  epSigma2 = runif(0,1);
  for(i=0; i < (*r-1); ++i) {
    ui[i] =  runif(0,1);
    epj[i] =  runif(0, 1);
  }
  gi0[0] = rgamma(1,1);
  gi0[1] = rgamma(1,1);
  
  split = (int)rint(runif(1, *r));
  /* common elements */

  for(i=0; i < (split -1); ++i) {
    candidatoMu[i] = mu[i];
    candidatoSigma2[i] = sigma2[i];
    for (j=0; j< (split-1); ++j) {
      candidatoBeta[i* (*r+1) + j] = beta[i* *r + j];
      candidatoQ[i* (*r+1) + j] = -candidatoBeta[i* (*r+1) + j];
    }
    for (j=split; j< *r; ++j) {
      candidatoBeta[i* (*r+1) + j+1] = beta[i* *r + j];
      candidatoQ[i* (*r+1) + j+1] = -candidatoBeta[i* (*r+1) + j+1];
    }
  }
  for(i=split; i < *r; ++i) {
    candidatoMu[i+1] = mu[i];
    candidatoSigma2[i+1] = sigma2[i];
    for (j=0; j< (split-1); ++j) {
      candidatoBeta[(i+1)* (*r+1) + j] = beta[i* *r + j];
      candidatoQ[(i+1)* (*r+1) + j] = -candidatoBeta[(i+1)* (*r+1) + j];
    }
    for (j=split; j< *r; ++j) {
      candidatoBeta[(i+1)* (*r+1) + (j+1)] = beta[i* *r + j];
      candidatoQ[(i+1)* (*r+1) + (j+1)] = -candidatoBeta[(i+1)* (*r+1) + (j+1)];
    }
  }

  /* rows i,i0 */

  for(i=0; i < (split-1); i++) {
    candidatoBeta[(split-1) * (*r+1) + i] = 
      beta[(split-1) * *r + i] - *tauSplitBeta * log(ui[i]);
    candidatoQ[(split-1) * (*r+1) + i] = - candidatoBeta[(split-1) * (*r+1) + i];
    candidatoBeta[split * (*r+1) + i] = 
      beta[(split-1) * *r + i] - *tauSplitBeta * log(1-ui[i]);
    candidatoQ[split * (*r+1) + i] = - candidatoBeta[split * (*r+1) + i];
  }

  for(i=split; i < *r; i++) {
    candidatoBeta[(split-1) * (*r+1) + i+1] = 
      beta[(split-1) * *r + i] - *tauSplitBeta * log(ui[i-1]);
    candidatoQ[(split-1) * (*r+1) + i+1] = - candidatoBeta[(split-1) * (*r+1)+ i+1];
    candidatoBeta[split * (*r+1) + i+1] = 
      beta[(split-1) * *r + i] - *tauSplitBeta * log(1-ui[i-1]);
    candidatoQ[split * (*r+1) + i+1] = - candidatoBeta[split * (*r+1) + i+1];
  }

  /* cols i0,j */

  for(i=0; i < (split-1); i++) {
    candidatoBeta[i* (*r+1) + (split-1)] = beta[i* *r + (split-1)] * epj[i];
    candidatoQ[i* (*r+1) + (split-1)] = - candidatoBeta[i* (*r+1) + (split-1)];
    candidatoBeta[i* (*r+1) + split] = beta[i* *r + (split-1)] * (1-epj[i]);
    candidatoQ[i* (*r+1) + split] = - candidatoBeta[i* (*r+1) + split];
  }

  for(i=split; i < *r; i++) {
    candidatoBeta[(i+1)* (*r+1) + (split-1)] = beta[i* *r + (split-1)] * epj[i-1];
    candidatoQ[(i+1)* (*r+1) + (split-1)] = - candidatoBeta[(i+1)* (*r+1) + (split-1)];
    candidatoBeta[(i+1)* (*r+1) + split] = beta[i* *r + (split-1)] * (1-epj[i-1]);
    candidatoQ[(i+1)* (*r+1) + split] = - candidatoBeta[(i+1)* (*r+1) + split];
  }

  /* center values */

  candidatoMu[split-1] = mu[split-1] - epMu;
  candidatoMu[split] = mu[split-1] + epMu;
  if (*varEqual) {
    candidatoSigma2[split-1] = sigma2[split-1];
    candidatoSigma2[split] = sigma2[split-1];
  }
  else {
    candidatoSigma2[split-1] = sigma2[split-1] * epSigma2;
    candidatoSigma2[split] = sigma2[split-1] * (1-epSigma2);
  }

  candidatoBeta[(*r+1) * (split-1) + split - 1] = 0;
  candidatoQ[(*r+1) * (split-1) + split - 1] = 0;
  candidatoBeta[(*r+1) * split + split] = 0;
  candidatoQ[(*r+1) * split + split] = 0;
  candidatoBeta[(*r+1) * (split-1) + split] = gi0[0];
  candidatoQ[(*r+1) * (split-1) + split] = -gi0[0];
  candidatoBeta[(*r+1) * split + split-1] = gi0[1];
  candidatoQ[(*r+1) * split + split-1] = -gi0[1];


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
      /* if variances are different we include the prior */
      if (! *varEqual) {
	probSplit = probSplit + 
	  log(1.0 / (2 *sqrt(*maxVar * candidatoSigma2[i])));
      }
    }
    for (i=0; i<*r; ++i) {
      probSplit = probSplit - dnorm(mu[i], *muAlfa, *muBeta, 1);
      if (! *varEqual) {
	probSplit = probSplit - 
	  log(1.0 / (2 * sqrt(*maxVar * sigma2[i])));
      }
    }
    /*  take out zeros of the diagonal */
    for (i=0; i<((*r+1) * (*r+1)); ++i) {
      if ((i % (*r + 2)))
	probSplit = probSplit + dgamma(candidatoBeta[i], 1, 1, 1);
    }
    for (i=0; i<(*r * *r); ++i) {
      if ((i % (*r + 1)))
	probSplit = probSplit - dgamma(beta[i], 1, 1, 1);
    }
    for (k=0; k < *genome; ++k) {
      nn = index[k+1] - index[k];
      //P1(1633);
      double *yy; yy = Calloc(nn, double);
      double *xx; xx = Calloc(nn, double);
      //P2(1633);
      for (j=0; j < nn-1; ++j) {
	yy[j] = y[j + index[k]];
	xx[j] = x[j + index[k]];
      }
      xx[nn-1] = 0;
      yy[nn-1] = y[nn-1 + index[k]];

      normalNHHMMlikelihood(yy, &newState, xx, &nn, candidatoQ,  
			    statSplit, candidatoMu, candidatoSigma2, 
			    &loglikPartial);
#ifdef DEBUG
      Rprintf("\n       After likelihood\n"); fflush(stdout);
#endif

      loglikCandidate += loglikPartial;
      /*       delete [] xx; */
      /*       delete [] yy; */
      Free(xx);
      Free(yy);

    }
    probSplit = probSplit + heat * (loglikCandidate - *loglikLast);
    probSplit = probSplit + log(1-ps[*r]) - log(ps[*r-1]) + log(*r) - log(*r+1);
    probSplit = probSplit + log(*r + 1);
    probSplit = probSplit - log(2);
    probSplit = probSplit - dnorm(epMu, 0, *tauSplitMu, 1);
    probSplit = probSplit - dgamma(gi0[0], 1, 1, 1);
    probSplit = probSplit - dgamma(gi0[1], 1, 1, 1);

    /*  jacobian of the transformation */
    jacobian = log(2);
    if (! *varEqual) jacobian = jacobian + log(sigma2[split-1]);

    for(i=0; i < (split-1); i++) {
      jacobian = jacobian + log(beta[i* *r + (split-1)]);
    }
    for(i=split; i < *r; i++) {
      jacobian = jacobian + log(beta[i* *r + (split-1)]);
    }
    for (i=0; i< *r-1; ++i) {
      jacobian = jacobian + log((*tauSplitBeta/(1-ui[i])) + 
				(*tauSplitBeta/ui[i]));
    }
/*     probSplit = fmin(1.0, exp(probSplit + jacobian)); */
    probSplit = probSplit + jacobian;
    if(probSplit >= 0) probSplit = 1; else probSplit = exp(probSplit);
    CHECK_NUM(probSplit, Split);
/*     CHECK_NUM(jacobian, Split); */
/*     if (probSplit >1) probSplit = 1; */
    if (runif(0,1) <= probSplit) {
      *accepted =1;
      *loglikSplit = loglikCandidate;
      /*        Rprintf("acepto SPLIT\n");  */
      /*     printf("SI--------->%f %f\n", candidatoMu[split], candidatoMu[split+1]); */
    }
    else {
      /*       Rprintf("SPLIT/COMBINE\n"); */
      /*       Rprintf("NO->%f %f\n", candidatoMu[split], candidatoMu[split+1]); */
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
  Rprintf("\n Exiting split \n"); fflush(stdout);
#endif
}

/* tauSplitBeta is now the constant 1/(1-x) */
void Combine(double *y, double *x, int *varEqual, int *genome, int *index, 
	     double *mu, double *sigma2, 
	     double *beta, double *stat, double *statCombine, 
	     int *r, double *loglikLast, double *probK, double *ps, 
	     int *n, double *candidatoMu,
	     double *candidatoSigma2, 
	     double *candidatoBeta, double *loglikCombine,
	     double *muAlfa, double *muBeta,
	     double *tauSplitMu,
	     double *tauSplitBeta, int *accepted, 
	     double *maxVar, double heat) {

#ifdef DEBUG
  Rprintf("\n Entering combine \n"); fflush(stdout);
#endif

  
  double *candidatoQ; candidatoQ = Calloc((*r-1)*(*r-1), double);
  double probCombine;
  int combine;
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
    

/* common elements */

  for(i=0; i < (combine -1); ++i) {
    candidatoMu[i] = mu[i];
    candidatoSigma2[i] = sigma2[i];
    for (j=0; j< (combine-1); ++j) {
      candidatoBeta[i* (*r-1) + j] = beta[i* *r + j];
      candidatoQ[i* (*r-1) + j] = -candidatoBeta[i* (*r-1) + j];
    }
    for (j=(combine+1); j< *r; ++j) {
      candidatoBeta[i* (*r-1) + j-1] = beta[i* *r + j];
      candidatoQ[i* (*r-1) + j-1] = -candidatoBeta[i* (*r-1) + j-1];
      }
    }
  for(i=(combine+1); i < *r; ++i) {
    candidatoMu[i-1] = mu[i];
    candidatoSigma2[i-1] = sigma2[i];
    for (j=0; j< (combine-1); ++j) {
      candidatoBeta[(i-1)* (*r-1) + j] = beta[i* *r + j];
      candidatoQ[(i-1)* (*r-1) + j] = -candidatoBeta[(i-1)* (*r-1) + j];
    }
    for (j=(combine+1); j< *r; ++j) {
      candidatoBeta[(i-1)* (*r-1) + (j-1)] = beta[i* *r + j];
      candidatoQ[(i-1)* (*r-1) + (j-1)] = -candidatoBeta[(i-1)* (*r-1) + (j-1)];
    }
 }

  /* rows i,i0 */

  for(i=0; i < (combine-1); i++) {
    candidatoBeta[(combine-1) * (*r-1) + i] = 
      -log(exp(-beta[(combine-1) * *r + i] * (1.0 / *tauSplitBeta))  + 
	      exp(-beta[combine * *r + i] * (1.0 / *tauSplitBeta))) * *tauSplitBeta;
    candidatoQ[(combine-1) * (*r-1) + i] = - candidatoBeta[(combine-1) * (*r-1) + i];
    ui[i] = exp(-beta[(combine-1) * *r + i] * (1.0/ *tauSplitBeta)) / 
	(exp(-beta[(combine-1) * *r + i] * (1.0/ *tauSplitBeta))  + 
	 exp(-beta[combine * *r + i] * (1.0/ *tauSplitBeta)));


}

  for(i=(combine+1); i < *r; i++) {
    candidatoBeta[(combine-1) * (*r-1) + i-1] = 
      -log(exp(-beta[(combine-1) * *r + i] * (1.0/ *tauSplitBeta)) + 
	   exp(-beta[combine * *r + i] * (1.0/ *tauSplitBeta))) * *tauSplitBeta;
    candidatoQ[(combine-1) * (*r-1) + i-1] = - candidatoBeta[(combine-1) * (*r-1) + i-1];
    ui[i-2] = exp(-beta[(combine-1) * *r + i] * (1.0/ *tauSplitBeta)) / 
      (exp(-beta[(combine-1) * *r + i] * (1.0/ *tauSplitBeta))  + 
       exp(-beta[combine * *r + i] * (1.0/ *tauSplitBeta)));

  }

/* cols i0,j */

  for(i=0; i < (combine-1); i++) {
    candidatoBeta[i* (*r-1) + (combine-1)] = 
      beta[i* *r + (combine-1)] + beta[i* *r + combine];
    candidatoQ[i* (*r-1) + (combine-1)] = - candidatoBeta[i* (*r-1) + (combine-1)];
    epj[i] = beta[i* *r + (combine-1)] / 
      (beta[i* *r + (combine-1)] + beta[i* *r + combine]);

  }

  for(i=(combine+1); i < *r; i++) {
    candidatoBeta[(i-1)* (*r-1) + (combine-1)] = 
      beta[i* *r + (combine-1)] + beta[i* *r + combine];
    candidatoQ[(i-1)* (*r-1) + (combine-1)] = 
      - candidatoBeta[(i-1)* (*r-1) + (combine-1)];
    epj[i-2] = beta[i* *r + (combine-1)] / 
      (beta[i* *r + (combine-1)] + beta[i* *r + combine]);
    
  }

/* i0,j0 values */

  candidatoMu[combine-1] = (mu[combine-1] + mu[combine]) / 2.0;
  if (*varEqual) {
    candidatoSigma2[combine-1] = sigma2[combine-1];
  }
  else {
    candidatoSigma2[combine-1] = sigma2[combine-1] + sigma2[combine];
    if (candidatoSigma2[combine-1] > *maxVar) reachedMaxVar = 1;
  }
  candidatoBeta[(*r-1) * (combine-1) + combine - 1] = 0;
  candidatoQ[(*r-1) * (combine-1) + combine - 1] = 0;

  /* check beta > 0 */
  for(i=0; i< (*r-1) * (*r-1); ++i) {
    if (candidatoBeta[i] < 0) {
    /* As have zero prob on the prior, we stop the move */
      reachedMaxVar = 1;
    }
  }
    /*********************************************************/

  if (!reachedMaxVar) {
    /*  recover the other auxiliary variables */
    epMu = (mu[combine] - mu[combine-1]) / 2;
    epSigma2 = sigma2[combine-1] / (sigma2[combine-1] + sigma2[combine]);
    gi0[0] = beta[*r * (combine-1) + combine];
    gi0[1] = beta[*r * combine + combine -1];
    newState = *r - 1;
    /* we compute the same ratio as split and then compute its inverse */
    probCombine = log(probK[*r-1]) - log(probK[*r-2]);
    /*  sum of priors */
    for (i=0; i<*r; ++i) {
      probCombine = probCombine + dnorm(mu[i], *muAlfa, *muBeta, 1);
      if (! *varEqual) {
	probCombine = probCombine + 
	  log(1.0/(2*sqrt(*maxVar * sigma2[i])));
      }
    }
    for (i=0; i<(*r-1); ++i) {
      probCombine = probCombine - dnorm(candidatoMu[i], *muAlfa, *muBeta, 1);
      if (! *varEqual) {
	probCombine = probCombine - 
	  log(1.0/(2*sqrt(*maxVar * candidatoSigma2[i])));
      }
    }
    /*  take out zeros of the diagonal */
    for (i=0; i<(*r * *r); ++i) {
      if ((i % (*r + 1))) {
	probCombine = probCombine + dgamma(beta[i], 1, 1, 1);
      }
    }
    for (i=0; i<((*r-1) * (*r-1)); ++i) {
      if ((i % *r)) {
	probCombine = probCombine - dgamma(candidatoBeta[i], 1, 1, 1);
      }
    }
    for (k=0; k < *genome; ++k) {
      nn = index[k+1] - index[k];
      //P1(1893);
      double *yy; yy = Calloc(nn, double);
      double *xx; xx = Calloc(nn, double);
      //P2(1893);
      for (j=0; j < nn-1; ++j) {
	yy[j] = y[j + index[k]];
	xx[j] = x[j + index[k]];
      }
      xx[nn-1] = 0;
      yy[nn-1] = y[nn-1 + index[k]];
      normalNHHMMlikelihood(yy, &newState, xx, &nn, candidatoQ,  
			    statCombine, candidatoMu, candidatoSigma2, 
			    &loglikPartial);
      loglikCandidate += loglikPartial;
      Free(xx);
      Free(yy);
    }
    probCombine = probCombine + heat * (*loglikLast - loglikCandidate);
    probCombine = probCombine + log(1.0 -ps[*r-1]) - log(ps[*r-2]);
    probCombine = probCombine + log(*r-1) - log(*r);
    probCombine = probCombine + log(*r) - log(2);
    probCombine = probCombine - dnorm(epMu, 0, *tauSplitMu, 1);
    probCombine = probCombine - dgamma(gi0[0], 1, 1, 1);
    probCombine = probCombine - dgamma(gi0[1], 1, 1, 1);
    /*  jacobian of the transformation */
    /*  r-1 because of a 2 in the denominator and the split equivalent is r-1 */
    jacobian = log(2);
    if (! *varEqual) jacobian = jacobian + log(sigma2[combine-1]);
    for(i=0; i < (combine-1); i++) {
      jacobian = jacobian + log(candidatoBeta[i* (*r-1) + (combine-1)]);
    }
    for(i=(combine+1); i < *r; i++) {
      jacobian = jacobian + log(candidatoBeta[(i-1)* (*r-1) + (combine-1)]);
    }
    for (i=0; i< (*r-2); ++i) {
      jacobian = jacobian + log((*tauSplitBeta/(1-ui[i])) + 
				  (*tauSplitBeta/ui[i]));
      }

/*     probCombine = exp(probCombine + jacobian) ; */
/*     probCombine = fmin(1.0, 1.0 /probCombine); */
    probCombine = probCombine + jacobian ;
    if(probCombine <= 0) probCombine = 1; else probCombine = 1/exp(probCombine);

    CHECK_NUM(probCombine, Combine);
/*     if (probCombine >1) probCombine = 1; */
    if (runif(0,1) <= probCombine) {
      *accepted =1;
      *loglikCombine = loglikCandidate;
    }
  }

  Free(epj);
  Free(ui);
  Free(gi0);
  Free(candidatoQ);

#ifdef DEBUG
  Rprintf("\n Exiting combine \n"); fflush(stdout);
#endif

}

void SwapChains(double *y, double *x,  
		int *genome, int *index, 
		double *stat,
		int n, int NC, int *kMax, 
		double *muAlfa, double *muBeta,
		double *muCoupled, 
		double *sigma2Coupled, 
		double *betaCoupled, 
		int *rCoupled, 
		double *loglikLastCoupled, 
		double *heat, int *accepted, 
		int *triedChangeDim) {

  /* we only temper the likelihood, so the priors disappear */

  int i, j, k;
  double probSwap = 0;
  i = (int)rint(runif(1, NC));
  j = i;
  while (i == j) {
    j = (int)rint(runif(1, NC));
  }
  double *mui; mui = Calloc(rCoupled[i-1], double);
  double *sigma2i; sigma2i = Calloc(rCoupled[i-1], double);
  double *betai; betai = Calloc(rCoupled[i-1] * 
				rCoupled[i-1], double);
  double *qi; qi = Calloc(rCoupled[i-1] * 
			     rCoupled[i-1], double);
  double *stati; stati = Calloc(rCoupled[i-1], double);
  double loglikCandidatei=0;
  double *muj; muj = Calloc(rCoupled[j-1], double);
  double *sigma2j; sigma2j = Calloc(rCoupled[j-1], double);
  double *betaj; betaj = Calloc(rCoupled[j-1] * 
				rCoupled[j-1], double);
  double *qj; qj = Calloc(rCoupled[j-1] * 
			     rCoupled[j-1], double);
  double *statj; statj = Calloc(rCoupled[j-1], double);
  double loglikCandidatej=0;
  for (k=0; k < rCoupled[i-1]; k++) {
    mui[k] = muCoupled[((i-1) * *kMax * (*kMax+1)/2) + 
		       (rCoupled[i-1]* (rCoupled[i-1]-1)/2) + k];
    stati[k] = stat[(rCoupled[i-1] * (rCoupled[i-1]-1)/2) + k];
    sigma2i[k] = sigma2Coupled[((i-1) * *kMax * (*kMax+1)/2) + 
			       (rCoupled[i-1]* (rCoupled[i-1]-1)/2) + k];
  }
  for (k=0; k < rCoupled[j-1]; k++) {
    muj[k] = muCoupled[((j-1) * *kMax * (*kMax+1)/2) + 
		       (rCoupled[j-1]* (rCoupled[j-1]-1)/2) + k];
    statj[k] = stat[(rCoupled[j-1] * (rCoupled[j-1]-1)/2) + k];
    sigma2j[k] = sigma2Coupled[((j-1) * *kMax * (*kMax+1)/2) + 
			       (rCoupled[j-1]* (rCoupled[j-1]-1)/2) + k];
  }
  for (k=0; k < rCoupled[i-1] * rCoupled[i-1]; k++) {
    betai[k] = betaCoupled[((i-1) * *kMax * (*kMax+1) * (2* *kMax+1) / 6) + 
			   ((rCoupled[i-1]-1)* 
			    rCoupled[i-1] * (2*rCoupled[i-1]-1)/6) + k];
    qi[k] = -betai[k];
  }

  for (k=0; k < rCoupled[j-1] * rCoupled[j-1]; k++) {
    betaj[k] = betaCoupled[((j-1) * *kMax * (*kMax+1) * (2* *kMax+1) / 6) + 
			   ((rCoupled[j-1]-1)* 
			    rCoupled[j-1] * (2*rCoupled[j-1]-1)/6) + k];
    qj[k] = -betaj[k];
  }
  /* likelihood */
 
  loglikCandidatei = loglikLastCoupled[((i-1) * *kMax) + rCoupled[i-1]-1];
  loglikCandidatej = loglikLastCoupled[((j-1) * *kMax) + rCoupled[j-1]-1];
  probSwap += heat[i-1] * (loglikCandidatej - loglikCandidatei);
  probSwap += heat[j-1] * (loglikCandidatei - loglikCandidatej);
  CHECK_NUM(probSwap, SwapChains);
  if (probSwap >= 0) probSwap = 1;
  else probSwap = exp(probSwap);
  if (runif(0,1) <= probSwap) {
    /* we have swapped the cool chain */
    if ((i==1) | (j==1)) *accepted = 2;
    else *accepted = 1;
    /* swap parameters */
    for (k=0; k<rCoupled[j-1]; k++) {
      muCoupled[((i-1) * *kMax * (*kMax+1)/2) + 
		(rCoupled[j-1]* (rCoupled[j-1]-1)/2) + k] = muj[k];
      sigma2Coupled[((i-1) * *kMax * (*kMax+1)/2) + 
		(rCoupled[j-1]* (rCoupled[j-1]-1)/2) + k] = sigma2j[k];
    }
    for (k=0; k<rCoupled[i-1]; k++) {
      muCoupled[((j-1) * *kMax * (*kMax+1)/2) + 
		(rCoupled[i-1]* (rCoupled[i-1]-1)/2) + k] = mui[k];
      sigma2Coupled[((j-1) * *kMax * (*kMax+1)/2) + 
		(rCoupled[i-1]* (rCoupled[i-1]-1)/2) + k] = sigma2i[k];
    }
    for (k=0; k<rCoupled[j-1] * rCoupled[j-1]; k++) {
      betaCoupled[((i-1) * *kMax * (*kMax+1) * (2* *kMax+1) / 6) + 
		  ((rCoupled[j-1]-1)* 
		   rCoupled[j-1] * (2*rCoupled[j-1]-1)/6) + k] = 
	betaj[k];
    }
    for (k=0; k<rCoupled[i-1] * rCoupled[i-1]; k++) {
      betaCoupled[((j-1) * *kMax * (*kMax+1) * (2* *kMax+1) / 6) + 
		  ((rCoupled[i-1]-1)* 
		   rCoupled[i-1] * (2*rCoupled[i-1]-1)/6) + k] = 
	betai[k];
    }
    double Auxl1 = loglikLastCoupled[((i-1) * *kMax) + rCoupled[i-1]-1];
    loglikLastCoupled[((i-1) * *kMax) + rCoupled[j-1]-1] = 
      loglikLastCoupled[((j-1) * *kMax) + rCoupled[j-1]-1];
    loglikLastCoupled[((j-1) * *kMax) + rCoupled[i-1]-1] = Auxl1;
    int Auxr = rCoupled[i-1];
    rCoupled[i-1] = rCoupled[j-1];
    rCoupled[j-1] = Auxr;
  }
  /* we have tried a transdimensional move with the cool chain*/
  if (((i==1) | (j==1)) & (rCoupled[i-1] != rCoupled[j-1])) {
    *triedChangeDim = 1;
  }
  Free(mui);
  Free(sigma2i);
  Free(betai);
  Free(qi);
  Free(stati);
  Free(muj);
  Free(sigma2j);
  Free(betaj);
  Free(qj);
  Free(statj);
}

void MetropolisUpdate(double *y, double *x, int *varEqual, 
		      int *genome, int *index, double *mu, 
		      double *sigma2, double *beta, double *stat, int *r, int *n, 
		      double *muAlfa, double *muBeta, 
		      double *sigmaTauMu,
		      double *sigmaTauSigma2, double *sigmaTauBeta,
		      double *loglikLast, double *maxVar, 
		      double heat) {
#ifdef DEBUG
  Rprintf("\n Entering MetropolisUpdate\n"); fflush(stdout);
#endif


  int i, j, k;
  double acepProb;
  double loglikCandidate = 0;
  double loglikPartial = 0;
  double *q; q = Calloc(*r * *r, double);
  double *candidatoMu; candidatoMu = Calloc(*r, double);
  double *candidatoSigma2; candidatoSigma2 = Calloc(*r, double);
  double *candidatoBeta; candidatoBeta = Calloc(*r * *r, double);
  int reachedMaxVar=0;
  int nn;
  /* update mu */
  
  acepProb = 0;
  for(i=0; i<*r; ++i) {
    candidatoMu[i] = mu[i] + rnorm(0, *sigmaTauMu);
    acepProb = acepProb + (dnorm(candidatoMu[i], *muAlfa, *muBeta, 1) -
				   dnorm(mu[i], *muAlfa, *muBeta, 1));
  }
  for (i=0; i<(*r * *r); ++i) {
    q[i] = -beta[i];
  }

  for (k=0; k < *genome; ++k) {
    nn = index[k+1] - index[k];
    //P1(2119);
    double *yy; yy = Calloc(nn, double);
    double *xx; xx = Calloc(nn, double);
    //P2(2119);
    for (j=0; j < nn-1; ++j) {
      yy[j] = y[j + index[k]];
      xx[j] = x[j + index[k]];
    }
    xx[nn-1] = 0;
    yy[nn-1] = y[nn-1 + index[k]];
    normalNHHMMlikelihood(yy, r, xx, &nn, q,  stat, candidatoMu, sigma2, 
			  &loglikPartial);
    loglikCandidate += loglikPartial;
    Free(xx);
    Free(yy);
  }
      
  acepProb = acepProb + heat * (loglikCandidate - *loglikLast);
  if(acepProb >= 0) acepProb = 1; else acepProb = exp(acepProb);
  CHECK_NUM(acepProb, MetropolisUpdateMu);

/*   if (acepProb > 1) acepProb = 1; */
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
    if (*varEqual) {
      acepProb = acepProb + 
	(log(1.0/(2*sqrt(*maxVar * candidatoSigma2[0]))) - 
	 log(1.0/(2*sqrt(*maxVar * sigma2[0]))));
    }
    else {
      for(i=0; i<*r; ++i) {
	acepProb = acepProb + 
	  (log(1.0/(2*sqrt(*maxVar * candidatoSigma2[i]))) - 
	   log(1.0/(2*sqrt(*maxVar * sigma2[i]))));
      }
    }
    for (k=0; k < *genome; ++k) {
      nn = index[k+1] - index[k];
      //P1(2177);
      double *yy; yy = Calloc(nn, double);
      double *xx; xx = Calloc(nn, double);
      //P2(2177);
      for (j=0; j < nn-1; ++j) {
	yy[j] = y[j + index[k]];
	xx[j] = x[j + index[k]];
      }
      xx[nn-1] = 0;
      yy[nn-1] = y[nn-1 + index[k]];      
      normalNHHMMlikelihood(yy, r, xx, &nn, q, stat, mu, candidatoSigma2, 
			    &loglikPartial);
      loglikCandidate += loglikPartial;
      Free(xx);
      Free(yy);
    }
    acepProb = acepProb + heat * (loglikCandidate - *loglikLast);
    /* If variances are equal we should not do this r times */

    if (*varEqual) {
      acepProb = acepProb + log(candidatoSigma2[0]) - log(sigma2[0]);
    }
    else {
      for (i=0; i<*r; ++i) {
	acepProb = acepProb + log(candidatoSigma2[i]) - log(sigma2[i]);
      }
    }
    if(acepProb >= 0) acepProb = 1; else acepProb = exp(acepProb);
    CHECK_NUM(acepProb, MetropolisUpdateSigma);
/*     if (acepProb > 1) acepProb = 1; */
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
      if ((i % (*r+1))) {
	candidatoBeta[i] = exp(log(beta[i]) + rnorm(0, *sigmaTauBeta));
	q[i] = -candidatoBeta[i];
	acepProb = acepProb + (dgamma(candidatoBeta[i], 1, 1, 1) -
				       dgamma(beta[i], 1, 1, 1));
      }
      else {
	candidatoBeta[i] = beta[i];
	q[i] = -candidatoBeta[i];
      }
    }
    for (k=0; k < *genome; ++k) {
      nn = index[k+1] - index[k];
      //P1(2229);
      double *yy; yy = Calloc(nn, double);
      double *xx; xx = Calloc(nn, double);
      //P2(2229);
      for (j=0; j < nn-1; ++j) {
	yy[j] = y[j + index[k]];
	xx[j] = x[j + index[k]];
      }
      xx[nn-1] = 0;
      yy[nn-1] = y[nn-1 + index[k]];
      normalNHHMMlikelihood(yy, r, xx, &nn, q,  stat, mu, sigma2, 
			    &loglikPartial);
      loglikCandidate += loglikPartial;
      Free(xx);
      Free(yy);
    }
    acepProb = acepProb + heat * (loglikCandidate - *loglikLast);
/*     Rprintf("\n  acepProb 1 = %f\n", acepProb); */
/*     Rprintf("\n  acepProb 2 = %f\n", acepProb); */

    for (i=0; i< *r * *r; ++i) {
      if (i % (*r+1)) {
	acepProb = acepProb + log(candidatoBeta[i]) - log(beta[i]);
/* 	Rprintf("\n       candidatoBeta[i] = %f\n", candidatoBeta[i]); */
/* 	Rprintf("\n       beta[i] = %f\n", beta[i]); */
      }

    }
    if(acepProb >= 0) acepProb = 1; else acepProb = exp(acepProb);
    CHECK_NUM(acepProb, MetropolisUpdateBeta);
/*     if (acepProb > 1) acepProb = 1; */
    if (runif(0,1) < acepProb) {
      for (i=0; i< (*r * *r); ++i) beta[i] = candidatoBeta[i];
      *loglikLast = loglikCandidate;
    }
  }
  Free(q);
  Free(candidatoMu);
  Free(candidatoSigma2);
  Free(candidatoBeta);

#ifdef DEBUG
  Rprintf("\n Exiting MetropolisUpdate \n"); fflush(stdout);
#endif


}




void viterbi_nogenome(double *y, double *x,  
		      const int k, const int n, double *mu, double *sigma2,
		      double *beta, double *stat, int *states) {
  /* Taken from nnhl-o23.c with some modifications */
  int i, j, l;
  int nn_minus_1 = n - 1;
  int i_minus_1;
  double xx_i_1_1;
  double rowSumsQ = 0.0;
  double log_rowSumsQ = 0.0;
  
  int * restrict bAux; bAux = Calloc(k, int);
  double * restrict sd_j; sd_j = Calloc(k, double);
  double * restrict logstat_j; logstat_j = Calloc(k, double);
  double * restrict mAux; mAux = Calloc(k, double);
  int ** restrict b; b = Calloc(n, int*);
  double ** restrict m; m = Calloc(n, double*);
  double ** restrict log_Q; log_Q = Calloc(k, double*);
  double ** restrict beta_2D_tr; beta_2D_tr = Calloc(k, double*);

  for (int ii = 0; ii < k; ii++) {
    log_Q[ii] = Calloc(k, double);
    beta_2D_tr[ii] = Calloc(k, double);
  }
  for(int jj = 0; jj < n; jj++) {
    m[jj] = Calloc(k, double);
    b[jj] = Calloc(k, int);
  }
  
/*   we can traverse the array jumping or traverse beta jumping; */
/*   I guess both are equally bad. */
  for (int row = 0; row < k; row++) {
    for(int col = 0; col < k; col++) {
      beta_2D_tr[col][row] = beta[row * k + col];
      //FIXME: use an accumulator instead of multiplication?
      // beta[row + k_mult]; k_mult += k; initialize k_mult = 0 in row outer loop.
	}
  }
  for (int ii = 0; ii < k; ii++) {
    sd_j[ii] = sqrt(sigma2[ii]);
    logstat_j[ii] = log(stat[ii]);
  }

  
  /*  Forward recursion */
  for(j=0; j<k; j++) {
    m[0][j] = logstat_j[j] + dnorm(y[0], mu[j], sd_j[j], 1);
  }
  for(i=1; i < n; i++){
    xx_i_1_1 = x[i - 1] - 1;
    i_minus_1 = i - 1;
    for (j=0; j < k; j++) {
      /*	log_Q[j][j] = 0.0; */
      /*       Diagonals of logQ are known (1.0). Only do off-diagnoals. */
      rowSumsQ = 1.0;
      for(l = (j + 1); l < k; l++) {
	log_Q[j][l] = beta_2D_tr[j][l] * xx_i_1_1;
	rowSumsQ += exp(log_Q[j][l]);
      }
      for(l = 0; l < j; l++) {
	log_Q[j][l] = beta_2D_tr[j][l] * xx_i_1_1;
	rowSumsQ += exp(log_Q[j][l]);
      }
      log_rowSumsQ = log(rowSumsQ); 
      log_Q[j][j] = -(log_rowSumsQ);
      for(l = (j + 1); l < k; l++) {
	log_Q[j][l] -= log_rowSumsQ;
      }
      for(l = 0; l < j; l++) {
	log_Q[j][l] -= log_rowSumsQ;
      }
    }
    for(j=0; j < k; j++) {
      for (l=0; l < k; l++) {
	mAux[l] = m[i_minus_1][l] + log_Q[l][j];
	bAux[l] = l + 1;
      }
      revsort(mAux, bAux, k);
      b[i][j] = bAux[0];
      m[i][j] = mAux[0] +
	dnorm(y[i], mu[j], sd_j[j], 1);

    }
  }

  /*  Backward recursion */
  for (j=0; j < k; j++) {
    mAux[j] = m[nn_minus_1][j];
    bAux[j] = j + 1; 
  }
  revsort(mAux, bAux, k);
  states[nn_minus_1] = bAux[0];

  for (i = (n - 2); i >= 0; --i) {
    states[i] = b[(i + 1)][states[1 + i] - 1];
  }
 
  Free(bAux);
  Free(mAux);
  Free_2(m, n);
  Free_2_int(b, n);
  Free(logstat_j);
  Free(sd_j);
  Free_2(log_Q, k);
  Free_2(beta_2D_tr, k);
  
}


  

void viterbi_genome(double *y, double *x, int *genome, int *index, 
		    int *k, int *n, double *mu, double *sigma2,
		    double *beta, double *stat, int *states,
		    struct Sequence **arraySeqRefs,
		    int write_seq) {
  int i, j, g, nn, l, index_g, nn_minus_1, index_g_plus_1;
  double rowSumsQ = 0.0;
  double log_rowSumsQ = 0.0;
  double *sd_j; sd_j = Calloc(*k, double);
  double *logstat_j; logstat_j = Calloc(*k, double);
  double **log_Q; log_Q = Calloc(*k, double*);
  double **beta_2D_tr; beta_2D_tr = Calloc(*k, double*);

  for (int ii = 0; ii < *k; ii++) {
    log_Q[ii] = Calloc(*k, double);
    beta_2D_tr[ii] = Calloc(*k, double);
    sd_j[ii] = sqrt(sigma2[ii]);
    logstat_j[ii] = log(stat[ii]);
  }

/*   we can traverse the array jumping or traverse beta jumping; */
/*   I guess both are equally bad. */
  for (int row = 0; row < *k; row++) {
    for(int col = 0; col < *k; col++) {
      beta_2D_tr[col][row] = beta[row * *k + col];
    }
  }
  
  for (g=0; g < *genome; g++) {
    index_g = index[g];
    nn = index[g+1] - index_g;
    nn_minus_1 = nn - 1;
    int *bAux; bAux = Calloc(*k, int);
    double *yy; yy = Calloc(nn, double);
    double *xx; xx = Calloc(nn, double);
    double *mAux; mAux = Calloc(*k, double);
    double **m; m = Calloc(nn, double*);
    int **b; b = Calloc(nn, int*);
    for(int jj = 0; jj < nn; jj++) {
      m[jj] = Calloc(*k, double);
      b[jj] = Calloc(*k, int);
    }
    for (i=0; i < nn_minus_1; i++) {
      yy[i] = y[i + index_g];
      xx[i] = x[i + index_g];
    }
    xx[nn_minus_1] = 0;
    yy[nn_minus_1] = y[nn_minus_1 + index_g];
    
    /*  Forward recursion */
    for(j=0; j<*k; j++) {
      m[0][j] = logstat_j[j] + dnorm(yy[0], mu[j], sd_j[j], 1);
    }
    for(i=1; i < nn; i++){
      double xx_i_1_1 = xx[i - 1] - 1;
      int i_minus_1 = i - 1;
      for (j=0; j < *k; j++) {
	/*	log_Q[j][j] = 0.0; */
	rowSumsQ = 1.0;
	for(l = (j + 1); l < *k; l++) {
	  log_Q[j][l] = beta_2D_tr[j][l] * xx_i_1_1;
	  rowSumsQ += exp(log_Q[j][l]);
	}
	for(l = 0; l < j; l++) {
	  log_Q[j][l] = beta_2D_tr[j][l] * xx_i_1_1;
	  rowSumsQ += exp(log_Q[j][l]);
	}
	
	log_rowSumsQ = log(rowSumsQ);

	log_Q[j][j] = -(log_rowSumsQ);
	for(l = (j + 1); l < *k; l++) {
	  log_Q[j][l] -= log_rowSumsQ;
	}
	for(l = 0; l < j; l++) {
	  log_Q[j][l] -= log_rowSumsQ;
	}
      }
      for(j=0; j < *k; j++) {
	for (l=0; l < *k; l++) {
	  mAux[l] = m[i_minus_1][l] + log_Q[l][j];
	  bAux[l] = l + 1;
	}
	revsort(mAux, bAux, *k);
	b[i][j] = bAux[0];
	m[i][j] = mAux[0] +
	  dnorm(yy[i], mu[j], sd_j[j], 1);
      }
    }

    /*  Backward recursion */
    for (j=0; j < *k; j++) {
      mAux[j] = m[nn_minus_1][j];
      bAux[j] = j + 1; 
    }
    revsort(mAux, bAux, *k);
    states[index_g + nn_minus_1] = bAux[0];
    index_g_plus_1 = index_g + 1;
    for (i = nn-2; i >= 0; --i) {
      states[index_g + i] = b[(i + 1)][states[index_g_plus_1 + i] - 1];
    }
    Free(xx);
    Free(yy);
    Free(bAux);
    Free(mAux);
    Free_2(m, nn);
    Free_2_int(b, nn);
#ifdef DEBUGV
    Rprintf("\n **** Doing genome value %i\n", g);
    PR3(*k);
    PR3(index_g);
#endif
    if(write_seq) 
      viterbi_to_Sequence(arraySeqRefs + g, *k, nn, states + index_g); 
  }
  /*  Free(bAux_0); */
  Free(logstat_j);
  Free(sd_j);
  Free_2(log_Q, *k);
  Free_2(beta_2D_tr, *k);
    
}



void viterbi(double *y, double *x, int *genome, int *index, 
	     int *k, int *n, double *mu, double *sigma2,
	     double *beta, double *stat, int *states,
	     struct Sequence **seqRef, struct Sequence **arraySeqRefs,
	     int write_seq, unsigned int *viterbi_counts) {
 #ifdef DEBUGW 
  if(write_seq) {
    Rprintf("\n Entered viterbi\n");
    fflush(stdout);
  }
#endif 
  
  if((*genome) > 1) {
    viterbi_genome(y, x, genome, index, 
		   k, n, mu, sigma2,
		   beta, stat, states,
		   arraySeqRefs,
		   write_seq);
  } else {
    viterbi_nogenome(y, x, *k, *n, mu, sigma2,
		     beta, stat, states);
    if(write_seq) {
      viterbi_to_Sequence(seqRef, *k, *n, states); 
    }
  }
  if(write_seq) 
    viterbi_counts[(*k) - 1]++;
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
	    states, NULL, NULL, 0, NULL);
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

  int i, j, countrep, repeated;
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
    error("ERROR in wholeViterbi: Can't open file for writing\n");
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
	    states, NULL, NULL, 0, NULL);
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
	      states, NULL, NULL, 0, NULL);
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



// Was called from getEdges in R. But we should use the stored
// Viterbi sequences!!

/* void edges(double *y, double *x, int *genome, int *index, */
/* 	   int *k, int *n, int *N, double *mu, double *sigma2, */
/* 	   double *beta, double *stat, int *count_edge) { */
  
/*   int i, j; */
/*   int *states, *lastSeq; */
/*   double *AuxMu, *AuxSigma2, *AuxBeta; */
  
/*   states = (int *) R_alloc(*n, sizeof(int)); */
  
/*   lastSeq = (int *) R_alloc(*n, sizeof(int)); */
/*   AuxMu = (double *) R_alloc(*k, sizeof(double)); */
/*   AuxSigma2 = (double *) R_alloc(*k, sizeof(double)); */
/*   AuxBeta = (double *) R_alloc(*k * *k, sizeof(double)); */
  
/*   for (j=0; j < *n; j++) count_edge[j] = 0; /\*initialize *\/ */
  
/*   for(i=0; i < *N; i++) { */

/*     for(j=0; j < *k; j++) { */
/*       AuxMu[j] = mu[i* *k + j]; */
/*       AuxSigma2[j] = sigma2[i* *k + j]; */
/*     } */
/*     for(j=0; j < *k * *k; j++) { */
/*       AuxBeta[j] = beta[i * *k * *k + j]; */
/*     } */
/*     viterbi(y, x, genome, index, k, n, AuxMu, AuxSigma2, AuxBeta, stat, */
/* 	    states, NULL, NULL, 0, NULL); */
/*     for (j=1; j < *n; j++) { */
/*       if(states[j - 1] != states[j]) { */
/* 	count_edge[j] += 1; */
/*       } */
/*     } */
/*   } */
/* } */


void doBurnin(double *y, double *x, int *varEqual, int *genome, 
	      int *index, int *kMax, int *n, int *burnin,
	      int *probB,
	      int *probD, int *probS, int *probC, int *probE, double *probK,
	      double *pb, double *ps,
	      double *muAlfa, double *muBeta,
	      double *s1, double *s2, double *initMu, 
	      double *initSigma2, 
	      double *initBeta,
	      double *sigmaTauMu, double *sigmaTauSigma2,
	      double *sigmaTauBeta, double *tauSplitMu,
	      double *tauSplitBeta,
	      double *muCoupled, double *sigma2Coupled,
	      double *betaCoupled, int *rCoupled, 
	      double *loglikLastCoupled, double *stat, 
	      int *startK, int *RJ, double *maxVar, 
	      int *NC, double *heat)  {


  /*  Initializa RNG */
  GetRNGstate();
#ifdef DEBUG
  double dummy_random_number;
  dummy_random_number = runif(0, 1);
  PR(dummy_random_number);
#endif

/*   Rprintf("doBurnintausplit=%f\n", *tauSplitBeta); */

  int nc;
  double *OldStatCoupled; OldStatCoupled = Calloc(*kMax, double);
  double *OldMuCoupled; OldMuCoupled = Calloc(*kMax, double);
  double *OldSigma2Coupled; OldSigma2Coupled = Calloc(*kMax, double);
  double *OldBetaCoupled; OldBetaCoupled = Calloc(*kMax * *kMax, double);
  double *qCoupled; qCoupled = Calloc(*kMax * *kMax, double);
  
  int i,j,m;
  int t;
  int nn;
  int *accepted; accepted = Calloc(1, int);
  int *triedChangeDim; triedChangeDim = Calloc(1, int);
  /* index to permutations */
  int *indexPerm; indexPerm = Calloc(*kMax, int);
  /* new parameters (max limit) */

  double *NewStatCoupled; NewStatCoupled = Calloc(*kMax, double);
  double *NewMuCoupled; NewMuCoupled = Calloc(*kMax, double);
  double *NewSigma2Coupled; NewSigma2Coupled = Calloc(*kMax, double);
  double *NewBetaCoupled; NewBetaCoupled = Calloc(*kMax * *kMax, double);
    
  double *q; q = Calloc(*kMax * *kMax, double);

  /*          Initializations */
  *accepted = 0;
  *triedChangeDim = 0;
  for (i=0; i<*kMax; ++i) indexPerm[i] = i+1;

  /* FIXME  do this appropriately: fill with true junk zz*/
  int junkindex = 0;
  for(junkindex = 0; junkindex < ((*kMax) * (*kMax)); ++junkindex) {
    q[junkindex] = -99;
  }
/*   Rprintf("Overdispersed init\n"); */

  /*  Overdispersed init */
  for (nc=0; nc < *NC; nc++) {
    for (i=0; i<*kMax; ++i) {
      for (j=0; j<=i; ++j) {
	OldMuCoupled[j] = runif(-2, 2);
      }
      sigma2Coupled[(nc * *kMax * (*kMax+1)/2) + 
		    (i* (i+1)/2)] = runif(0, *maxVar);
      for (j=1; j<=i; ++j) {
	if (varEqual) {
	  sigma2Coupled[(nc * *kMax * (*kMax+1)/2) + 
			(i* (i+1)/2) + j] = 
	    sigma2Coupled[(nc * *kMax * (*kMax+1)/2) + 
			  (i* (i+1)/2)];
	}
	else {
	  sigma2Coupled[(nc * *kMax * (*kMax+1)/2) + 
			(i* (i+1)/2) + j] = runif(0, *maxVar);
	}
      }
      /*  relabelling */
      R_rsort(OldMuCoupled, i+1);
      for (j=0; j<=i; ++j) {
	muCoupled[(nc * *kMax * (*kMax+1)/2) + 
		  (i* (i+1)/2) + j] = OldMuCoupled[j];
      }
      for (j=0; j<((i+1)*(i+1)); ++j) {
	if (!(j % (i+2))) {
	  betaCoupled[(nc * *kMax * (*kMax+1) * (2* *kMax+1) / 6) + 
		      (i* (i+1) * (2*i+1)/6) + j] = 0;
	}
	else {
	  betaCoupled[(nc * *kMax * (*kMax+1) * (2* *kMax+1) / 6) + 
		      (i* (i+1) * (2*i+1)/6) + j] = runif(0,10);
	}
      }
    }
  }

  /* starting values supplied by user */
  /* we assume that thar sorted by mean value */
  if (*startK) {
    i = *startK - 1;
    for (j=0; j<=i; ++j) {
      muCoupled[i * (i+1)/2 + j] = initMu[j];
      sigma2Coupled[(i* (i+1)/2) + j] = initSigma2[j];
    }
    for (j=0; j<((i+1)*(i+1)); ++j) {
      betaCoupled[(i* (i+1) * (2*i+1)/6) + j] = initBeta[j];
    }
  }
  /* End starting values supplied by user */

  /*  loglik of start values */
  for (nc=0; nc < *NC; nc++) {
    for (i=1; i<=*kMax; i++) {
      for (j=0; j<i; j++) {
	OldStatCoupled[j] = stat[i * (i-1) /2 + j];
	OldMuCoupled[j] = muCoupled[(nc * *kMax * (*kMax+1)/2) + 
				    (i* (i-1)/2) + j];
	OldSigma2Coupled[j] = sigma2Coupled[(nc * *kMax * (*kMax+1)/2) + 
					    (i* (i-1)/2) + j];
      }
      for (j=0; j <i*i; ++j) {
	OldBetaCoupled[j] = betaCoupled[(nc * *kMax * (*kMax+1) * (2* *kMax+1) / 6) + 
					(i* (i-1) * (2*i-1)/6) + j];
	qCoupled[j] = -OldBetaCoupled[j];
      }
      for (m=0; m < *genome; ++m) {
	nn = index[m+1] - index[m];
	//P1(2851);
	double *yy; yy = Calloc(nn, double);
	double *xx; xx = Calloc(nn, double);
	//P2(2851);
	/*  Initialization */
	for (j=0; j < nn-1; ++j) {
	  yy[j] = -9999;
	  xx[j] = -9999;
	}
	double loglikPartial = 0;
	for (j=0; j < nn-1; ++j) {
	  yy[j] = y[j + index[m]];
	  xx[j] = x[j + index[m]];
	}
	xx[nn-1] = 0;
	yy[nn-1] = y[index[m + 1] - 1];
	  
	normalNHHMMlikelihood(yy, &i, xx, &nn, qCoupled, OldStatCoupled, 
			      OldMuCoupled, OldSigma2Coupled, 
			      &loglikPartial);
	loglikLastCoupled[(nc * *kMax) + (i-1)] += loglikPartial;
	Free(xx);
	Free(yy);
      }
    }
  }

  
  for (nc=0; nc < *NC; nc++) {
    rCoupled[nc] = (int)rint(runif(1, *kMax));
  }
  if (*startK!=0) {
    rCoupled[0] = *startK;
  }

  /*  Loop MCMC iterations */
  for(t=0; t<*burnin; ++t) {
    /* Allow R interrupts; check every 100 iterations */
    if (!(t % 100))
      R_CheckUserInterrupt(); 

    /* METROPOLIS UPDATE */
    
    for(nc=0; nc < *NC ; nc++) {
      /* old parameters */
      for (i=0; i<rCoupled[nc]; ++i) {
	OldStatCoupled[i] = stat[(rCoupled[nc] * (rCoupled[nc]-1) / 2) + i];
	OldMuCoupled[i] = muCoupled[(nc * *kMax * (*kMax+1)/2) + 
				    (rCoupled[nc]* (rCoupled[nc]-1)/2) + i];
	OldSigma2Coupled[i] = sigma2Coupled[(nc * *kMax * (*kMax+1)/2) + 
					    (rCoupled[nc]* (rCoupled[nc]-1)/2) + i];
      }
      for (i=0; i<rCoupled[nc] * rCoupled[nc]; ++i) {
	OldBetaCoupled[i] = betaCoupled[(nc * *kMax * (*kMax+1) * (2* *kMax+1) / 6) + 
					((rCoupled[nc]-1)* 
					 rCoupled[nc] * (2*rCoupled[nc]-1)/6) + i];
      }
      MetropolisUpdate(y, x, varEqual, genome, index, OldMuCoupled, 
		       OldSigma2Coupled, OldBetaCoupled, OldStatCoupled, &rCoupled[nc], 
		       n, muAlfa, muBeta, &sigmaTauMu[rCoupled[nc]-1], 
		       &sigmaTauSigma2[rCoupled[nc]-1], 
		       &sigmaTauBeta[rCoupled[nc]-1], 
		       &loglikLastCoupled[nc * *kMax + (rCoupled[nc]-1)], 
		       maxVar, heat[nc]);
      rsort_with_index(OldMuCoupled, indexPerm, rCoupled[nc]);
      for (i=0; i <rCoupled[nc]; ++i) {
	muCoupled[(nc * *kMax * (*kMax+1)/2) + 
		  (rCoupled[nc]* (rCoupled[nc]-1)/2) + i] = OldMuCoupled[i];
	sigma2Coupled[(nc * *kMax * (*kMax+1)/2) + 
		      (rCoupled[nc]* (rCoupled[nc]-1)/2) + i] = 
	  OldSigma2Coupled[indexPerm[i]-1];
      }
      for (i=0; i<rCoupled[nc]; ++i) {
	for (j=0; j<rCoupled[nc]; ++j) {
	  betaCoupled[(nc * *kMax * (*kMax+1) * (2* *kMax+1) / 6) + 
		      ((rCoupled[nc]-1)* rCoupled[nc] * 
		       (2*rCoupled[nc]-1)/6) + i* rCoupled[nc] + j] = 
	    OldBetaCoupled[(indexPerm[i]-1)*rCoupled[nc] + indexPerm[j]-1];
	  qCoupled[i*rCoupled[nc] +j] = 
	    -OldBetaCoupled[(indexPerm[i]-1)*rCoupled[nc] + indexPerm[j]-1];
	}
      }
      for(i=0; i<*kMax; ++i) indexPerm[i]= i+1;
    }
    /***********************************************************/
    /***********************************************************/
    /***********************************************************/

    if(*RJ) {
    /*  Birth or death */

    for (nc=0; nc < *NC; nc++) {
      if(runif(0,1) <= pb[rCoupled[nc]-1]) {
	for (i=0; i<rCoupled[nc]; ++i) {
	  OldStatCoupled[i] = stat[(rCoupled[nc] * (rCoupled[nc]-1) /2) + i];
	  OldMuCoupled[i] = muCoupled[(nc * *kMax * (*kMax+1)/2) +
				      (rCoupled[nc]* (rCoupled[nc]-1)/2) + i];
	  OldSigma2Coupled[i] = sigma2Coupled[(nc * *kMax * (*kMax+1)/2) +
					      (rCoupled[nc]* (rCoupled[nc]-1)/2) + i];
	}
	for (i=0; i<rCoupled[nc] * rCoupled[nc]; ++i) {
	  OldBetaCoupled[i] = betaCoupled[(nc * *kMax * (*kMax+1) * (2* *kMax+1) / 6) +
					  ((rCoupled[nc]-1)*
					   rCoupled[nc] * (2*rCoupled[nc]-1)/6) + i];
	}
	for (i=0; i<(rCoupled[nc]+1); ++i) {
	  NewStatCoupled[i] = stat[(rCoupled[nc] * (rCoupled[nc] + 1)/ 2) + i];
	}
	Birth(y, x, varEqual, genome, index, OldMuCoupled, OldSigma2Coupled,
	      OldBetaCoupled, OldStatCoupled, NewStatCoupled,
	      &rCoupled[nc],
	      &loglikLastCoupled[(nc * *kMax) + rCoupled[nc]-1], probK,
	      pb, muAlfa, muBeta, n, NewMuCoupled, NewSigma2Coupled,
	      NewBetaCoupled,
	      &loglikLastCoupled[(nc * *kMax) + rCoupled[nc]],
	      accepted, maxVar, s1, s2, heat[nc]);
	if (*accepted) {
	  if (nc==0) {
	    if (*accepted==1) probB[0] = probB[0] + 1;
	    else probB[1] = probB[1] + 1;
	  }
	  rCoupled[nc] = rCoupled[nc]+1;
	  *accepted = 0;
	  /* Relabelling */
	  rsort_with_index(NewMuCoupled, indexPerm, rCoupled[nc]);
	  for (i=0; i<rCoupled[nc]; ++i) {
	    muCoupled[(nc * *kMax * (*kMax+1)/2) +
		      (rCoupled[nc]* (rCoupled[nc]-1)/2) + i] =
	      NewMuCoupled[i];
	    sigma2Coupled[(nc * *kMax * (*kMax+1)/2) +
			  (rCoupled[nc]* (rCoupled[nc]-1)/2) + i] =
	      NewSigma2Coupled[indexPerm[i]-1];
	  }
	  for (i=0; i<rCoupled[nc]; ++i) {
	    for (j=0; j<rCoupled[nc]; ++j) {
	      betaCoupled[(nc * *kMax * (*kMax+1) * (2* *kMax+1) / 6) +
			  ((rCoupled[nc]-1)* rCoupled[nc] *
			   (2*rCoupled[nc]-1)/6) + i* rCoupled[nc] + j] =
		NewBetaCoupled[(indexPerm[i]-1)*rCoupled[nc] + indexPerm[j]-1];
	      qCoupled[i*rCoupled[nc] +j] =
		-NewBetaCoupled[(indexPerm[i]-1)*rCoupled[nc] + indexPerm[j]-1];
	    }
	  }
	  /* Recover original permutations */
	  for (i=0; i<*kMax; ++i) indexPerm[i] = i+1;
	}
      }
      else {
	for (i=0; i<rCoupled[nc]; ++i) {
	  OldStatCoupled[i] = stat[(rCoupled[nc] * (rCoupled[nc]-1)/2) + i];
	  OldMuCoupled[i] = muCoupled[(nc * *kMax * (*kMax+1)/2) +
				      (rCoupled[nc]* (rCoupled[nc]-1)/2) + i];
	  OldSigma2Coupled[i] = sigma2Coupled[(nc * *kMax * (*kMax+1)/2) +
					      (rCoupled[nc]* (rCoupled[nc]-1)/2) + i];
	}
	for (i=0; i<rCoupled[nc] * rCoupled[nc]; ++i) {
	  OldBetaCoupled[i] = betaCoupled[(nc * *kMax * (*kMax+1) * (2* *kMax+1) / 6) +
					  ((rCoupled[nc]-1)*
					   rCoupled[nc] * (2*rCoupled[nc]-1)/6) + i];
	}
	for (i=0; i < rCoupled[nc]-1; ++i) {
	  NewStatCoupled[i] = stat[((rCoupled[nc]-1) * (rCoupled[nc]-2) / 2)  + i];
	}
	Death(y, x, genome, index, OldMuCoupled, OldSigma2Coupled,
	      OldBetaCoupled,
	      OldStatCoupled, NewStatCoupled,
	      muAlfa, muBeta, &rCoupled[nc],
	      &loglikLastCoupled[(nc * *kMax) + rCoupled[nc]-1], probK,
	      pb, n, NewMuCoupled, NewSigma2Coupled,
	      NewBetaCoupled, &loglikLastCoupled[(nc * *kMax) + rCoupled[nc]-2],
	      accepted, s1, s2, heat[nc]);
	if (*accepted) {
	  if (nc==0) {
	    if (*accepted==1) probD[0] = probD[0] + 1;
	    else probD[1] = probD[1] + 1;
	  }
	  rCoupled[nc] = rCoupled[nc]-1;
	  *accepted = 0;
	  for (i=0; i<rCoupled[nc]; ++i) {
	    muCoupled[(nc * *kMax * (*kMax+1)/2) +
		      (rCoupled[nc]* (rCoupled[nc]-1)/2) + i] = NewMuCoupled[i];
	    sigma2Coupled[(nc * *kMax * (*kMax+1)/2) +
			  (rCoupled[nc]* (rCoupled[nc]-1)/2) + i] =
	      NewSigma2Coupled[i];
	  }
	  for (i=0; i<rCoupled[nc]; ++i) {
	    for (j=0; j<rCoupled[nc]; ++j) {
	      betaCoupled[(nc * *kMax * (*kMax+1) * (2* *kMax+1) / 6) +
			  ((rCoupled[nc]-1)* rCoupled[nc] *
			   (2*rCoupled[nc]-1)/6) + i* rCoupled[nc] + j] =
		NewBetaCoupled[i*rCoupled[nc] + j];
	      qCoupled[i*rCoupled[nc] +j] = -NewBetaCoupled[i*rCoupled[nc] + j];
	    }
	    
	  }
	}
      }
    }

    /********************************************************************/
    /********************************************************************/
    /********************************************************************/
    /********************************************************************/
    /*********************** SPLIT / COMBINE ****************************/
    for (nc=0; nc < *NC; nc++) {
      if(runif(0,1) <= ps[rCoupled[nc]-1]) {
	for (i=0; i<rCoupled[nc]; ++i) {
	  OldStatCoupled[i] = stat[(rCoupled[nc] * (rCoupled[nc]-1) /2) + i];
	  OldMuCoupled[i] = muCoupled[(nc * *kMax * (*kMax+1)/2) +
				      (rCoupled[nc]* (rCoupled[nc]-1)/2) + i];
	  OldSigma2Coupled[i] = sigma2Coupled[(nc * *kMax * (*kMax+1)/2) +
					      (rCoupled[nc]* (rCoupled[nc]-1)/2) + i];
	}
	for (i=0; i<rCoupled[nc] * rCoupled[nc]; ++i) {
	  OldBetaCoupled[i] = betaCoupled[(nc * *kMax * (*kMax+1) * (2* *kMax+1) / 6) +
					  ((rCoupled[nc]-1)*
					   rCoupled[nc] * (2*rCoupled[nc]-1)/6) + i];
	}
	for (i=0; i<(rCoupled[nc]+1); ++i) {
	  NewStatCoupled[i] = stat[(rCoupled[nc] * (rCoupled[nc] + 1)/ 2) + i];
	}
	Split(y, x, varEqual, genome, index, OldMuCoupled, OldSigma2Coupled,
	      OldBetaCoupled, OldStatCoupled, NewStatCoupled,
	      &rCoupled[nc], 
	      &loglikLastCoupled[(nc * *kMax) + rCoupled[nc]-1], 
	      probK, ps, n, NewMuCoupled,
	      NewSigma2Coupled, NewBetaCoupled, 
	      &loglikLastCoupled[(nc * *kMax) + rCoupled[nc]],
	      muAlfa, muBeta, tauSplitMu,
	      tauSplitBeta, accepted, maxVar, heat[nc]);

	if (*accepted) {
	  if (nc==0) *probS = *probS + 1;
	  rCoupled[nc] = rCoupled[nc]+1;
	  *accepted = 0;
	  /* Relabelling */
	  rsort_with_index(NewMuCoupled, indexPerm, rCoupled[nc]);
	  for (i=0; i<rCoupled[nc]; ++i) {
	    muCoupled[(nc * *kMax * (*kMax+1)/2) +
		      (rCoupled[nc]* (rCoupled[nc]-1)/2) + i] = 
	      NewMuCoupled[i];
	    sigma2Coupled[(nc * *kMax * (*kMax+1)/2) +
			  (rCoupled[nc]* (rCoupled[nc]-1)/2) + i] =
	      NewSigma2Coupled[indexPerm[i]-1];
	  }
	  for (i=0; i<rCoupled[nc]; ++i) {
	    for (j=0; j<rCoupled[nc]; ++j) {
	      betaCoupled[(nc * *kMax * (*kMax+1) * (2* *kMax+1) / 6) +
			  ((rCoupled[nc]-1)* rCoupled[nc] *
			   (2*rCoupled[nc]-1)/6) + i* rCoupled[nc] + j] =
		NewBetaCoupled[(indexPerm[i]-1)*rCoupled[nc] + indexPerm[j]-1];
	      qCoupled[i*rCoupled[nc] +j] =
		-NewBetaCoupled[(indexPerm[i]-1)*rCoupled[nc] + indexPerm[j]-1];
	    }
	  }
	  /* Recover original permutations */
	  for (i=0; i<*kMax; ++i) indexPerm[i] = i+1;
	}
      }
      else {
	for (i=0; i<rCoupled[nc]; ++i) {
	  OldStatCoupled[i] = stat[(rCoupled[nc] * (rCoupled[nc]-1)/2) + i];
	  OldMuCoupled[i] = muCoupled[(nc * *kMax * (*kMax+1)/2) +
				      (rCoupled[nc]* (rCoupled[nc]-1)/2) + i];
	  OldSigma2Coupled[i] = sigma2Coupled[(nc * *kMax * (*kMax+1)/2) +
					      (rCoupled[nc]* (rCoupled[nc]-1)/2) + i];
	}
	for (i=0; i<rCoupled[nc] * rCoupled[nc]; ++i) {
	  OldBetaCoupled[i] = betaCoupled[(nc * *kMax * (*kMax+1) * (2* *kMax+1) / 6) +
					  ((rCoupled[nc]-1)*
					   rCoupled[nc] * (2*rCoupled[nc]-1)/6) + i];
	}
	for (i=0; i < rCoupled[nc]-1; ++i) {
	  NewStatCoupled[i] = stat[((rCoupled[nc]-1) * (rCoupled[nc]-2) / 2)  + i];
	}
 	Combine(y, x, varEqual, genome, index, OldMuCoupled, OldSigma2Coupled, 
 		OldBetaCoupled, OldStatCoupled, NewStatCoupled, 
 		&rCoupled[nc], &loglikLastCoupled[(nc * *kMax) + rCoupled[nc]-1],
		probK, ps, n, NewMuCoupled, 
 		NewSigma2Coupled, NewBetaCoupled, 
		&loglikLastCoupled[(nc * *kMax) + rCoupled[nc]-2], muAlfa, 
 		muBeta, tauSplitMu, 
 		tauSplitBeta, accepted, maxVar, heat[nc]); 
	if (*accepted) {
	  if (nc==0)  *probC = *probC + 1;
	  rCoupled[nc] = rCoupled[nc]-1;
	  *accepted = 0;
	  for (i=0; i<rCoupled[nc]; ++i) {
	    muCoupled[(nc * *kMax * (*kMax+1)/2) +
		      (rCoupled[nc]* (rCoupled[nc]-1)/2) + i] = NewMuCoupled[i];
	    sigma2Coupled[(nc * *kMax * (*kMax+1)/2) +
			  (rCoupled[nc]* (rCoupled[nc]-1)/2) + i] =
	      NewSigma2Coupled[i];
	  }
	  for (i=0; i<rCoupled[nc]; ++i) {
	    for (j=0; j<rCoupled[nc]; ++j) {
	      betaCoupled[(nc * *kMax * (*kMax+1) * (2* *kMax+1) / 6) +
			  ((rCoupled[nc]-1)* rCoupled[nc] *
			   (2*rCoupled[nc]-1)/6) + i* rCoupled[nc] + j] =
		NewBetaCoupled[i*rCoupled[nc] + j];
	      qCoupled[i*rCoupled[nc] +j] = -NewBetaCoupled[i*rCoupled[nc] + j];
	    }
	    
	  }
	}
      }
    }
    }
    /* Exchange move */
    if (*NC > 1) {
      SwapChains(y, x, genome, index, stat, *n, *NC,
		 kMax, muAlfa, muBeta, muCoupled, sigma2Coupled, betaCoupled,
		 rCoupled, loglikLastCoupled, heat, accepted, 
		 triedChangeDim);
      if (*accepted==1) {
	probE[0] = probE[0] + 1;
      }
      if (*accepted==2) {
	/* We have swapped the state of the main chain */
	probE[1] = probE[1] + 1;
      }
      *accepted = 0;
      *triedChangeDim = 0;
    }
    /* End exchange move */
  }

  Free(OldMuCoupled);
  Free(OldSigma2Coupled);
  Free(OldBetaCoupled);
  Free(OldStatCoupled);
  Free(qCoupled);
  Free(accepted);
  Free(triedChangeDim);
  Free(q);
  Free(NewStatCoupled);
  Free(NewMuCoupled);
  Free(NewSigma2Coupled);
  Free(NewBetaCoupled);
  Free(indexPerm);
  
  PutRNGstate();
  }
  


void MetropolisSweep(double *y, double *x, int *varEqual, int *genome, 
		     int *index, int *kMax, int *n, int *burnin,
		     int *TOT, int *times, int *probB,
		     int *probD, int *probS, int *probC, int *probE,
		     double *probK,
		     double *pb, double *ps,
		     double *muAlfa, double *muBeta,
		     double *s1, double *s2, double *initMu, 
		     double *initSigma2, 
		     double *initBeta,
		     double *sigmaTauMu, double *sigmaTauSigma2,
		     double *sigmaTauBeta, double *tauSplitMu,
		     double *tauSplitBeta,
		     int *k, double *mu, double *sigma2,
		     double *beta, double *stat, 
		     int *startK, int *RJ, double *maxVar, 
		     double *probStates, double *loglik,
		     int *NC, double *deltaT,
		     int *write_seq,
		     char **filename,
		     int *num_sequences)  {

/*   Rprintf("firts tusplit=%f \n", *tauSplitBeta); */

  /*  Initialize RNG */
  GetRNGstate();
#ifdef DEBUG
  double dummy_random_number;
  dummy_random_number = runif(0, 1);
  PR(dummy_random_number);
#endif

  //  Rprintf("\n  +++++++ ENTERING MetropolisSweep in C +++++++ \n");
  /* Coupled parallel chains */
  double *muCoupled; muCoupled = Calloc(*NC * *kMax * (*kMax+1) / 2, double);
  double *sigma2Coupled; sigma2Coupled = Calloc(*NC * *kMax * (*kMax+1) / 2, double);
  double *betaCoupled; betaCoupled = Calloc(*NC * *kMax * (*kMax+1) * (2* *kMax+1) / 6, 
					    double);
  //  Rprintf("\n    ++++++++++ allocation 1 done +++++++\n");
  
  int nc;
  double *loglikLastCoupled; loglikLastCoupled = Calloc(*NC * *kMax, double);
  double *OldStatCoupled; OldStatCoupled = Calloc(*kMax, double);
  double *OldMuCoupled; OldMuCoupled = Calloc(*kMax, double);
  //  Rprintf("\n    ++++++++++ allocation 2 done +++++++\n");

  double *OldSigma2Coupled; OldSigma2Coupled = Calloc(*kMax, double);
  double *OldBetaCoupled; OldBetaCoupled = Calloc(*kMax * *kMax, double);
  double *qCoupled; qCoupled = Calloc(*kMax * *kMax, double);
  //  Rprintf("\n    ++++++++++ allocation 3 done +++++++\n");

  int *rCoupled; rCoupled = Calloc(*NC, int);
  double *heat; heat=Calloc(*NC, double);
  for (nc=0; nc < *NC; nc++) {
    heat[nc] = 1.0 / (1 + *deltaT * nc);
  }
  /* End variables for coupled parallel chains */

  //    Rprintf("\n    ++++++++++ allocation 4 done +++++++\n");

  int i,j,m;
  int t;
  int *states; states = Calloc(*n, int);
  int *accepted; accepted = Calloc(1, int);
  int *triedChangeDim; triedChangeDim = Calloc(1, int);
  int mainAccepted = 0;
  //  Rprintf("\n    ++++++++++ allocation 5 done +++++++\n");

  /* index to permutations */
  int *indexPerm; indexPerm = Calloc(*kMax, int);
  /* index to start every object */
  int *indexMu; indexMu = Calloc(*kMax, int); /* same for mu and sigma2 */
  int *indexBeta; indexBeta = Calloc(*kMax, int);
  int *indexStat; indexStat = Calloc(*kMax, int);
  /* new parameters (max limit) */

  //    Rprintf("\n    ++++++++++ allocation 6 done +++++++\n");

  double *NewStatCoupled; NewStatCoupled = Calloc(*kMax, double);
  double *NewMuCoupled; NewMuCoupled = Calloc(*kMax, double);
  double *NewSigma2Coupled; NewSigma2Coupled = Calloc(*kMax, double);
  double *NewBetaCoupled; NewBetaCoupled = Calloc(*kMax * *kMax, double);
  //  Rprintf("\n    ++++++++++ allocation 7 done +++++++\n");

  double *q; q = Calloc(*kMax * *kMax, double);

  unsigned int *viterbi_counts; 
  viterbi_counts = (unsigned int *) R_alloc(*kMax, sizeof(unsigned int));
  for(int dd = 0; dd < (*kMax); dd++) viterbi_counts[dd] = 0;
  //  Rprintf("\n    ++++++++++ allocation 8 done +++++++\n");

  unsigned int k_sum[(*kMax)];
  unsigned int Total_k; 
  for(int dd = 0; dd < (*kMax); dd++) k_sum[dd] = 0;
  Total_k = 0;

  
  /*   Genome deserves especial treatment. We could use a single entity
       because of pointer and array exchangeability. But to make things
       clear now, I keep them separate. */
  
  struct Sequence *CompactSeq;  CompactSeq = NULL;
  struct Sequence *arraySeqRefs[*genome];
  for (int jj = 0; jj < (*genome); jj++)
    arraySeqRefs[jj] = NULL;

  // FIXME: elimination of beta: keep only current value
  
  /*          Initializations */
  *accepted = 0;
  *triedChangeDim = 0;
  for (i=0; i<*kMax; ++i) indexPerm[i] = i+1;
  for (i=0; i<*kMax; ++i) states[i] = -99;

  /* FIXME  do this appropriately: fill with true junk zz*/
  int junkindex = 0;
  for(junkindex = 0; junkindex < ((*kMax) * (*kMax)); ++junkindex) {
    q[junkindex] = -99;
    NewBetaCoupled[junkindex] = -99;
  }
  for(junkindex = 0; junkindex < (*kMax) ; ++junkindex) {
    indexMu[junkindex] = -99;
    indexBeta[junkindex] = -99;
    indexStat[junkindex] = -99;
    NewStatCoupled[junkindex] = -99;
    NewMuCoupled[junkindex] = -99;
    NewSigma2Coupled[junkindex] = -99;
  }


  indexMu[0] = 0;
  indexBeta[0] = 0;
  indexStat[0] = 0;
  for(i=1; i<*kMax; ++i) {
    indexMu[i] = (*TOT)*i + indexMu[i-1];
    indexBeta[i] = (*TOT)*i*i + indexBeta[i-1];
    indexStat[i] = i + indexStat[i-1];
  }
  Rprintf("       Start burn-in\n");
  doBurnin(y=y, x=x, varEqual=varEqual, genome=genome, 
	   index=index, kMax=kMax, n=n, burnin=burnin,
	   probB=probB, 
	   probD=probD, probS=probS, probC=probC, probE=probE, probK=probK,
	   pb=pb, ps=ps, muAlfa=muAlfa, muBeta=muBeta,
	   s1=s1, s2=s2, initMu=initMu, 
	   initSigma2=initSigma2, 
	   initBeta=initBeta,
	   sigmaTauMu=sigmaTauMu, sigmaTauSigma2=sigmaTauSigma2,
	   sigmaTauBeta=sigmaTauBeta, tauSplitMu=tauSplitMu,
	   tauSplitBeta=tauSplitBeta,
	   muCoupled=muCoupled, sigma2Coupled=sigma2Coupled,
	   betaCoupled=betaCoupled, rCoupled=rCoupled, 
	   loglikLastCoupled=loglikLastCoupled, stat=stat, 
	   startK=startK, RJ=RJ, maxVar=maxVar, 
	   NC=NC, heat=heat);
  Rprintf("       End burn-in\n");

  /*  Take the last of the burn-in */
  for (i=0; i<*kMax; ++i) {
    for (j=0; j<=i; ++j) {
      mu[indexMu[i] + j] = muCoupled[(i*(i+1)/2) + j];
      sigma2[indexMu[i] + j] = sigma2Coupled[(i*(i+1)/2) + j];
    }
    for (j=0; j<((i+1)*(i+1)); ++j) {
      beta[indexBeta[i] +j] = betaCoupled[(i* (i+1) * (2*i+1)/6) + j];
    }
  }

  /*  loglik of start values */
  for (i=1; i<=*kMax;++i) {
    loglik[(i-1)* *TOT] = loglikLastCoupled[i-1];
    for (j=0; j<i; ++j) {
      OldStatCoupled[j] = stat[indexStat[i-1] + j];
      OldMuCoupled[j] = mu[indexMu[i-1] + j];
      OldSigma2Coupled[j] = sigma2[indexMu[i-1] + j];
    }
    for (j=0; j < i*i; ++j) {
      OldBetaCoupled[j] = beta[indexBeta[i-1] + j];
      q[j] = - OldBetaCoupled[j];
    }
    /*     viterbi*/
    if(i == 1) {
      add_constant_Sequence(genome, index, n, &CompactSeq,
			    arraySeqRefs, *write_seq, viterbi_counts);
    } else if (i >1) {
      viterbi(y, x, genome, index, &i, n, OldMuCoupled, 
	      OldSigma2Coupled, OldBetaCoupled, OldStatCoupled,
	      states, &CompactSeq, arraySeqRefs, *write_seq, viterbi_counts);
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
  
  k[0] = rCoupled[0];
  /*  Next free position */
  for (i=0; i <*kMax; ++i) times[i] = 1;

  /*  Loop MCMC iterations */
  for(t=1; t<*TOT; ++t) {

    /* Allow R interrupts; check every 100 iterations */
    if (!(t % 100))
      R_CheckUserInterrupt(); 

    /* METROPOLIS UPDATE */

    for(nc=0; nc < *NC ; nc++) {
      /* old parameters */
      for (i=0; i<rCoupled[nc]; ++i) {
	OldStatCoupled[i] = stat[indexStat[rCoupled[nc]-1] + i];
	OldMuCoupled[i] = muCoupled[(nc * *kMax * (*kMax+1)/2) + 
				    (rCoupled[nc]* (rCoupled[nc]-1)/2) + i];
	OldSigma2Coupled[i] = sigma2Coupled[(nc * *kMax * (*kMax+1)/2) + 
					    (rCoupled[nc]* (rCoupled[nc]-1)/2) + i];
      }
      for (i=0; i<rCoupled[nc] * rCoupled[nc]; ++i) {
	OldBetaCoupled[i] = betaCoupled[(nc * *kMax * (*kMax+1) * (2* *kMax+1) / 6) + 
					((rCoupled[nc]-1)* 
					 rCoupled[nc] * (2*rCoupled[nc]-1)/6) + i];
      }
      MetropolisUpdate(y, x, varEqual, genome, index, OldMuCoupled, 
		       OldSigma2Coupled, OldBetaCoupled, OldStatCoupled, &rCoupled[nc], 
		       n, muAlfa, muBeta, &sigmaTauMu[rCoupled[nc]-1], 
		       &sigmaTauSigma2[rCoupled[nc]-1], 
		       &sigmaTauBeta[rCoupled[nc]-1], 
		       &loglikLastCoupled[nc * *kMax + (rCoupled[nc]-1)], 
		       maxVar, heat[nc]);
      rsort_with_index(OldMuCoupled, indexPerm, rCoupled[nc]);
      for (i=0; i <rCoupled[nc]; ++i) {
	muCoupled[(nc * *kMax * (*kMax+1)/2) + 
		  (rCoupled[nc]* (rCoupled[nc]-1)/2) + i] = OldMuCoupled[i];
	sigma2Coupled[(nc * *kMax * (*kMax+1)/2) + 
		      (rCoupled[nc]* (rCoupled[nc]-1)/2) + i] = 
	  OldSigma2Coupled[indexPerm[i]-1];
      }
      for (i=0; i<rCoupled[nc]; ++i) {
	for (j=0; j<rCoupled[nc]; ++j) {
	  betaCoupled[(nc * *kMax * (*kMax+1) * (2* *kMax+1) / 6) + 
		      ((rCoupled[nc]-1)* rCoupled[nc] * 
		       (2*rCoupled[nc]-1)/6) + i* rCoupled[nc] + j] = 
	    OldBetaCoupled[(indexPerm[i]-1)*rCoupled[nc] + indexPerm[j]-1];
	  qCoupled[i*rCoupled[nc] +j] = 
	    -OldBetaCoupled[(indexPerm[i]-1)*rCoupled[nc] + indexPerm[j]-1];
	}
      }
      for(i=0; i<*kMax; ++i) indexPerm[i]= i+1;
    }
      
    /* Main chain */
    loglik[*TOT * (rCoupled[0]-1) + times[rCoupled[0]-1]] = 
      loglikLastCoupled[rCoupled[0]-1];
      
    /*  Auxiliaries for viterbi */
    for (i=0; i< rCoupled[0]; ++i) {
      OldStatCoupled[i] = stat[(rCoupled[0] * (rCoupled[0] -1 )/2) + i];
      OldMuCoupled[i] = 
	muCoupled[(rCoupled[0] * (rCoupled[0] -1) /2) + i];
      OldSigma2Coupled[i] = 
	sigma2Coupled[(rCoupled[0] * (rCoupled[0] -1) /2) + i];
    }
    for (i=0; i<rCoupled[0]; ++i) {
      for (j=0; j<rCoupled[0]; ++j) {
	OldBetaCoupled[i*rCoupled[0] +j] = 
	  betaCoupled[((rCoupled[0]-1)* rCoupled[0] * 
		       (2*rCoupled[0]-1)/6) + i* rCoupled[0] + j];
      }
    }
      
    /*  viterbi */
    if(rCoupled[0] == 1) {
      add_constant_Sequence(genome, index, n, &CompactSeq,
			    arraySeqRefs, *write_seq, viterbi_counts);
    } else if (rCoupled[0] >1) {
      viterbi(y, x, genome, index, &rCoupled[0], n, 
	      OldMuCoupled, OldSigma2Coupled, OldBetaCoupled, 
	      OldStatCoupled, states, &CompactSeq, arraySeqRefs, *write_seq, viterbi_counts);
      for (i=0; i < rCoupled[0]-1; ++i) {
	for (j=0; j < *n; ++j) {
	  if(states[j]==i+1) {
	    probStates[(*n * (i + (rCoupled[0]-1) * (rCoupled[0]-2) / 2)) + j] =
	      (1.0 / (double)(times[rCoupled[0]-1]+1)) +
	      probStates[(*n * (i + (rCoupled[0]-1) * (rCoupled[0]-2) / 2)) + j] *
	      (double)(times[rCoupled[0]-1]) /
	      (double)(times[rCoupled[0]-1]+1);
	  }
	  else {
	    probStates[(*n * (i + (rCoupled[0]-1) * (rCoupled[0]-2) / 2)) + j] =
	      probStates[(*n * (i + (rCoupled[0]-1) * (rCoupled[0]-2) / 2)) + j] *
	      (double)(times[rCoupled[0]-1]) /
	      (double)(times[rCoupled[0]-1]+1);
	  }
	}
      }
    }
  
    /*  save updates */
    for (i=0; i <rCoupled[0]; ++i) {
      mu[indexMu[rCoupled[0]-1] + (times[rCoupled[0]-1]*rCoupled[0]) + i] = 
	OldMuCoupled[i];
      sigma2[indexMu[rCoupled[0]-1] + (times[rCoupled[0]-1]*rCoupled[0]) + i] = OldSigma2Coupled[i];
    }
    for (i=0; i<rCoupled[0]; ++i) {
      for (j=0; j<rCoupled[0]; ++j) {
	beta[indexBeta[rCoupled[0]-1] + 
	     (times[rCoupled[0]-1]*rCoupled[0]*rCoupled[0]) + 
	     i*rCoupled[0] +j] = OldBetaCoupled[i*rCoupled[0] +j];
	q[i*rCoupled[0] +j] = -OldBetaCoupled[i*rCoupled[0] + j];
      }
    }
    times[rCoupled[0]-1]++;

    if(*RJ) {
    /* Coupled chains */
    for (nc=0; nc < *NC; nc++) {
      for (i=0; i<rCoupled[nc]; ++i) {
	OldStatCoupled[i] = stat[indexStat[rCoupled[nc]-1] + i];
	OldMuCoupled[i] = muCoupled[(nc * *kMax * (*kMax+1)/2) +
				    (rCoupled[nc]* (rCoupled[nc]-1)/2) + i];
	OldSigma2Coupled[i] = sigma2Coupled[(nc * *kMax * (*kMax+1)/2) +
					    (rCoupled[nc]* (rCoupled[nc]-1)/2) + i];
      }
      for (i=0; i<rCoupled[nc] * rCoupled[nc]; ++i) {
	OldBetaCoupled[i] = betaCoupled[(nc * *kMax * (*kMax+1) * (2* *kMax+1) / 6) +
					((rCoupled[nc]-1)*
					 rCoupled[nc] * (2*rCoupled[nc]-1)/6) + i];
      }
      if(runif(0,1) <= pb[rCoupled[nc]-1]) {
	for (i=0; i<(rCoupled[nc]+1); ++i) {
	  NewStatCoupled[i] = stat[indexStat[rCoupled[nc]] + i];
	}
	Birth(y, x, varEqual, genome, index, OldMuCoupled, OldSigma2Coupled,
	      OldBetaCoupled, OldStatCoupled, NewStatCoupled,
	      &rCoupled[nc],
	      &loglikLastCoupled[(nc * *kMax) + rCoupled[nc]-1], probK,
	      pb, muAlfa, muBeta, n, NewMuCoupled, NewSigma2Coupled, NewBetaCoupled,
	      &loglikLastCoupled[(nc * *kMax) + rCoupled[nc]],
	      accepted, maxVar, s1, s2, heat[nc]);
	if (*accepted) {
	  /* save parameters */
	  rCoupled[nc] = rCoupled[nc]+1;
	  rsort_with_index(NewMuCoupled, indexPerm, rCoupled[nc]);
	  for (i=0; i<rCoupled[nc]; ++i) {
	    muCoupled[(nc * *kMax * (*kMax+1)/2) +
		      (rCoupled[nc]* (rCoupled[nc]-1)/2) + i] =
	      NewMuCoupled[i];
	    sigma2Coupled[(nc * *kMax * (*kMax+1)/2) +
			  (rCoupled[nc]* (rCoupled[nc]-1)/2) + i] =
	      NewSigma2Coupled[indexPerm[i]-1];
	  }
	  for (i=0; i<rCoupled[nc]; ++i) {
	    for (j=0; j<rCoupled[nc]; ++j) {
	      betaCoupled[(nc * *kMax * (*kMax+1) * (2* *kMax+1) / 6) +
			  ((rCoupled[nc]-1)* rCoupled[nc] *
			   (2*rCoupled[nc]-1)/6) + i* rCoupled[nc] + j] =
		NewBetaCoupled[(indexPerm[i]-1)*rCoupled[nc] + indexPerm[j]-1];
	      qCoupled[i*rCoupled[nc] +j] =
		-NewBetaCoupled[(indexPerm[i]-1)*rCoupled[nc] + indexPerm[j]-1];
	    }
	  }
	  /* Recover original permutations */
	  for (i=0; i<*kMax; ++i) indexPerm[i] = i+1;
	  if (nc==0) {
	    mainAccepted = 1;
	    if (*accepted==1) probB[0] = probB[0] + 1;
	    else probB[1] = probB[1] + 1;
	  }
	  *accepted = 0;
	}
      }
      else {
	for (i=0; i < rCoupled[nc]-1; ++i) {
	  NewStatCoupled[i] = stat[indexStat[rCoupled[nc]-2] + i];
	}
	Death(y, x, genome, index, OldMuCoupled, OldSigma2Coupled,
	      OldBetaCoupled,
	      OldStatCoupled, NewStatCoupled,
	      muAlfa, muBeta, &rCoupled[nc],
	      &loglikLastCoupled[(nc * *kMax) + rCoupled[nc]-1], probK,
	      pb, n, NewMuCoupled, NewSigma2Coupled,
	      NewBetaCoupled, &loglikLastCoupled[(nc * *kMax) + rCoupled[nc]-2],
	      accepted, s1, s2, heat[nc]);
	if (*accepted) {
	  /* save parameters */
	  rCoupled[nc] = rCoupled[nc]-1;
	  for (i=0; i<rCoupled[nc]; ++i) {
	    muCoupled[(nc * *kMax * (*kMax+1)/2) +
		      (rCoupled[nc]* (rCoupled[nc]-1)/2) + i] =
	      NewMuCoupled[i];
	    sigma2Coupled[(nc * *kMax * (*kMax+1)/2) +
			  (rCoupled[nc]* (rCoupled[nc]-1)/2) + i] =
	      NewSigma2Coupled[i];
	  }
	  for (i=0; i<rCoupled[nc]; ++i) {
	    for (j=0; j<rCoupled[nc]; ++j) {
	      betaCoupled[(nc * *kMax * (*kMax+1) * (2* *kMax+1) / 6) +
			  ((rCoupled[nc]-1)* rCoupled[nc] *
			   (2*rCoupled[nc]-1)/6) + i* rCoupled[nc] + j] =
		NewBetaCoupled[i*rCoupled[nc] + j];
	      qCoupled[i*rCoupled[nc] +j] = -NewBetaCoupled[i*rCoupled[nc] + j];
	    }
	  }
	  if (nc==0) {
	    mainAccepted = 1;
	    if (*accepted==1) probD[0] = probD[0] + 1;
	    else probD[1] = probD[1] + 1;
	  }
	  *accepted = 0;
	}
      }
      /* save main parameters */
      if (mainAccepted) {
	mainAccepted = 0;
	loglik[*TOT * (rCoupled[0]-1) + times[rCoupled[0]-1]] =
	  loglikLastCoupled[rCoupled[0]-1];

	/*  Auxiliaries for viterbi */
	for (i=0; i< rCoupled[0]; ++i) {
	  OldStatCoupled[i] = stat[(rCoupled[0] * (rCoupled[0] -1 )/2) + i];
	  OldMuCoupled[i] = muCoupled[(rCoupled[0] * (rCoupled[0]-1)/2) +i];
	  OldSigma2Coupled[i] = sigma2Coupled[(rCoupled[0] * (rCoupled[0]-1)/2) +i];
	}
	for (i=0; i<rCoupled[0]; ++i) {
	  for (j=0; j<rCoupled[0]; ++j) {
	    OldBetaCoupled[i*rCoupled[0] +j] =
	      betaCoupled[((rCoupled[0]-1)* rCoupled[0] *
			   (2*rCoupled[0]-1)/6) + i* rCoupled[0] + j];
	  }
	}
	/*   viterbi  */
	if(rCoupled[0] == 1) {
	  add_constant_Sequence(genome, index, n, &CompactSeq,
				arraySeqRefs, *write_seq, viterbi_counts);
	} else if (rCoupled[0] >1) {
	  viterbi(y, x, genome, index, &rCoupled[0], n,
		  OldMuCoupled, OldSigma2Coupled, OldBetaCoupled,
		  OldStatCoupled, states, &CompactSeq, arraySeqRefs, *write_seq, viterbi_counts);
	  for (i=0; i < rCoupled[0]-1; ++i) {
	    for (j=0; j < *n; ++j) {
	      if(states[j]==i+1) {
		probStates[(*n * (i + (rCoupled[0]-1) * (rCoupled[0]-2) / 2)) + j] =
		  (1.0 / (double)(times[rCoupled[0]-1]+1)) +
		  probStates[(*n * (i + (rCoupled[0]-1) * (rCoupled[0]-2) / 2)) + j] *
		  (double)(times[rCoupled[0]-1]) /
		  (double)(times[rCoupled[0]-1]+1);
	      }
	      else {
		probStates[(*n * (i + (rCoupled[0]-1) * (rCoupled[0]-2) / 2)) + j] =
		  probStates[(*n * (i + (rCoupled[0]-1) * (rCoupled[0]-2) / 2)) + j] *
		  (double)(times[rCoupled[0]-1]) /
		  (double)(times[rCoupled[0]-1]+1);
	      }
	    }
	  }
	}
	/*  save updates */
	for (i=0; i <rCoupled[0]; ++i) {
	  mu[indexMu[rCoupled[0]-1] + (times[rCoupled[0]-1]*rCoupled[0]) + i] = OldMuCoupled[i];
	  sigma2[indexMu[rCoupled[0]-1] + (times[rCoupled[0]-1]*rCoupled[0]) + i] = OldSigma2Coupled[i];
	}
	for (i=0; i<rCoupled[0]; ++i) {
	  for (j=0; j<rCoupled[0]; ++j) {
	    beta[indexBeta[rCoupled[0]-1] +
		 (times[rCoupled[0]-1]*rCoupled[0]*rCoupled[0]) + i*rCoupled[0] +j] =
	      OldBetaCoupled[i*rCoupled[0] +j];
	    q[i*rCoupled[0] +j] = -OldBetaCoupled[i*rCoupled[0] + j];
	  }
	}
	times[rCoupled[0]-1]++;
      }
    }
    k[3*t-2] = rCoupled[0];

    /***********************************************************/
    /***********************************************************/
    /***********************************************************/
    /*  Split or combine */

    for (nc=0; nc < *NC; nc++) {
      for (i=0; i<rCoupled[nc]; ++i) {
	OldStatCoupled[i] = stat[indexStat[rCoupled[nc]-1] + i];
	OldMuCoupled[i] = muCoupled[(nc * *kMax * (*kMax+1)/2) +
				    (rCoupled[nc]* (rCoupled[nc]-1)/2) + i];
	OldSigma2Coupled[i] = sigma2Coupled[(nc * *kMax * (*kMax+1)/2) +
					    (rCoupled[nc]* (rCoupled[nc]-1)/2) + i];
      }
      for (i=0; i<rCoupled[nc] * rCoupled[nc]; ++i) {
	OldBetaCoupled[i] = betaCoupled[(nc * *kMax * (*kMax+1) * (2* *kMax+1) / 6) +
					((rCoupled[nc]-1)*
					 rCoupled[nc] * (2*rCoupled[nc]-1)/6) + i];
      }
      if(runif(0,1) <= ps[rCoupled[nc]-1]) {
	for (i=0; i<(rCoupled[nc]+1); ++i) {
	  NewStatCoupled[i] = stat[indexStat[rCoupled[nc]] + i];
	}
	Split(y, x, varEqual, genome, index, OldMuCoupled, OldSigma2Coupled,
	      OldBetaCoupled, OldStatCoupled, NewStatCoupled,
	      &rCoupled[nc], 
	      &loglikLastCoupled[(nc * *kMax) + rCoupled[nc]-1], 
	      probK, ps, n, NewMuCoupled,
	      NewSigma2Coupled, NewBetaCoupled, 
	      &loglikLastCoupled[(nc * *kMax) + rCoupled[nc]],
	      muAlfa, muBeta, tauSplitMu,
	      tauSplitBeta, accepted, maxVar, heat[nc]);
	if (*accepted) {
	  /* save parameters */
	  rCoupled[nc] = rCoupled[nc]+1;
	  rsort_with_index(NewMuCoupled, indexPerm, rCoupled[nc]);
	  for (i=0; i<rCoupled[nc]; ++i) {
	    muCoupled[(nc * *kMax * (*kMax+1)/2) +
		      (rCoupled[nc]* (rCoupled[nc]-1)/2) + i] = 
	      NewMuCoupled[i];
	    sigma2Coupled[(nc * *kMax * (*kMax+1)/2) +
			  (rCoupled[nc]* (rCoupled[nc]-1)/2) + i] =
	      NewSigma2Coupled[indexPerm[i]-1];
	  }
	  for (i=0; i<rCoupled[nc]; ++i) {
	    for (j=0; j<rCoupled[nc]; ++j) {
	      betaCoupled[(nc * *kMax * (*kMax+1) * (2* *kMax+1) / 6) +
			  ((rCoupled[nc]-1)* rCoupled[nc] *
			   (2*rCoupled[nc]-1)/6) + i* rCoupled[nc] + j] =
		NewBetaCoupled[(indexPerm[i]-1)*rCoupled[nc] + indexPerm[j]-1];
	      qCoupled[i*rCoupled[nc] +j] = 
		-NewBetaCoupled[(indexPerm[i]-1)*rCoupled[nc] + indexPerm[j]-1];
	    }
	  }
	  /* Recover original permutations */
	  for (i=0; i<*kMax; ++i) indexPerm[i] = i+1;
	  if (nc==0) {
	    mainAccepted = 1;
	    *probS = *probS + 1;
	  }
	  *accepted = 0;
	}
      }
      else {
	for (i=0; i < rCoupled[nc]-1; ++i) {
	  NewStatCoupled[i] = stat[indexStat[rCoupled[nc]-2] + i];
	}
	Combine(y, x, varEqual, genome, index, OldMuCoupled, OldSigma2Coupled, 
		OldBetaCoupled, OldStatCoupled, NewStatCoupled, 
		&rCoupled[nc], &loglikLastCoupled[(nc * *kMax) + rCoupled[nc]-1],
		probK, ps, n, NewMuCoupled, 
		NewSigma2Coupled, NewBetaCoupled, 
		&loglikLastCoupled[(nc * *kMax) + rCoupled[nc]-2], muAlfa, 
		muBeta, tauSplitMu, 
		tauSplitBeta, accepted, maxVar, heat[nc]); 
	if (*accepted) {

	  /* save parameters */
	  rCoupled[nc] = rCoupled[nc]-1;
	  for (i=0; i<rCoupled[nc]; ++i) {
	    muCoupled[(nc * *kMax * (*kMax+1)/2) +
		      (rCoupled[nc]* (rCoupled[nc]-1)/2) + i] = 
	      NewMuCoupled[i];
	    sigma2Coupled[(nc * *kMax * (*kMax+1)/2) +
			  (rCoupled[nc]* (rCoupled[nc]-1)/2) + i] =
	      NewSigma2Coupled[i];
	  }
	  for (i=0; i<rCoupled[nc]; ++i) {
	    for (j=0; j<rCoupled[nc]; ++j) {
	      betaCoupled[(nc * *kMax * (*kMax+1) * (2* *kMax+1) / 6) +
			  ((rCoupled[nc]-1)* rCoupled[nc] *
			   (2*rCoupled[nc]-1)/6) + i* rCoupled[nc] + j] =
		NewBetaCoupled[i*rCoupled[nc] + j];
	      qCoupled[i*rCoupled[nc] +j] = -NewBetaCoupled[i*rCoupled[nc] + j];
	    }
	  }
	  if (nc==0) {
	    mainAccepted = 1;
	    *probC = *probC + 1;
	  }
	  *accepted = 0;
	}
      }
      /* save main parameters */
      if (mainAccepted) {
	mainAccepted = 0;
	loglik[*TOT * (rCoupled[0]-1) + times[rCoupled[0]-1]] = 
	  loglikLastCoupled[rCoupled[0]-1];
	/*  Auxiliaries for viterbi */
	for (i=0; i< rCoupled[0]; ++i) {
	  OldStatCoupled[i] = stat[(rCoupled[0] * (rCoupled[0] -1 )/2) + i];
	  OldMuCoupled[i] = muCoupled[(rCoupled[0] * (rCoupled[0]-1)/2) +i];
	  OldSigma2Coupled[i] = sigma2Coupled[(rCoupled[0] * (rCoupled[0]-1)/2) +i];
	}
	for (i=0; i<rCoupled[0]; ++i) {
	  for (j=0; j<rCoupled[0]; ++j) {
	    OldBetaCoupled[i*rCoupled[0] +j] = 
	      betaCoupled[((rCoupled[0]-1)* rCoupled[0] * 
			   (2*rCoupled[0]-1)/6) + i* rCoupled[0] + j];
	  }
	}
	/*   viterbi  */
	if(rCoupled[0] == 1) {
	  add_constant_Sequence(genome, index, n, &CompactSeq,
				arraySeqRefs, *write_seq, viterbi_counts);
	} else if (rCoupled[0] >1) {
	  viterbi(y, x, genome, index, &rCoupled[0], n, 
		  OldMuCoupled, OldSigma2Coupled, OldBetaCoupled, 
		  OldStatCoupled, states, &CompactSeq, arraySeqRefs, *write_seq, viterbi_counts);
	  for (i=0; i < rCoupled[0]-1; ++i) {
	    for (j=0; j < *n; ++j) {
	      if(states[j]==i+1) {
		probStates[(*n * (i + (rCoupled[0]-1) * (rCoupled[0]-2) / 2)) + j] =
		  (1.0 / (double)(times[rCoupled[0]-1]+1)) +
		  probStates[(*n * (i + (rCoupled[0]-1) * (rCoupled[0]-2) / 2)) + j] *
		  (double)(times[rCoupled[0]-1]) /
		  (double)(times[rCoupled[0]-1]+1);
	      }
	      else {
		probStates[(*n * (i + (rCoupled[0]-1) * (rCoupled[0]-2) / 2)) + j] =
		  probStates[(*n * (i + (rCoupled[0]-1) * (rCoupled[0]-2) / 2)) + j] *
		  (double)(times[rCoupled[0]-1]) /
		  (double)(times[rCoupled[0]-1]+1);
	      }
	    }
	  }
	}
	/*  save updates */
	for (i=0; i <rCoupled[0]; ++i) {
	  mu[indexMu[rCoupled[0]-1] + (times[rCoupled[0]-1]*rCoupled[0]) + i] = OldMuCoupled[i];
	  sigma2[indexMu[rCoupled[0]-1] + (times[rCoupled[0]-1]*rCoupled[0]) + i] = OldSigma2Coupled[i];
	}
	for (i=0; i<rCoupled[0]; ++i) {
	  for (j=0; j<rCoupled[0]; ++j) {
	    beta[indexBeta[rCoupled[0]-1] + 
		 (times[rCoupled[0]-1]*rCoupled[0]*rCoupled[0]) + i*rCoupled[0] +j] = 
	      OldBetaCoupled[i*rCoupled[0] +j];
	    q[i*rCoupled[0] +j] = -OldBetaCoupled[i*rCoupled[0] + j];
	  }
	}
	times[rCoupled[0]-1]++;
      }
    }
    k[3*t-1] = rCoupled[0];
    }
    /* Exchange move */
    if (*NC > 1) {
      SwapChains(y, x, genome, index, stat, *n, *NC,
		 kMax, muAlfa, muBeta, muCoupled, sigma2Coupled, betaCoupled,
		 rCoupled, loglikLastCoupled, heat, accepted, 
		 triedChangeDim);

      if (*accepted==1) {
	probE[0] = probE[0] + 1;
	*accepted = 0;
      }
      if (*accepted==2) {
	/* We have swapped the state of the main chain */
	probE[1] = probE[1] + 1;
	*accepted = 0;
	for (i=0; i< rCoupled[0]; i++) {
	  NewStatCoupled[i] = stat[(rCoupled[0]* (rCoupled[0]-1)/2) + i];
	  NewMuCoupled[i] = muCoupled[(rCoupled[0]* (rCoupled[0]-1)/2) + i];
	  NewSigma2Coupled[i] = sigma2Coupled[(rCoupled[0]* (rCoupled[0]-1)/2) + i];
	}
	for (i=0; i< rCoupled[0]*rCoupled[0]; i++) {
	  NewBetaCoupled[i] = betaCoupled[((rCoupled[0]-1)* rCoupled[0]
					   * (2*rCoupled[0]-1)/6) + i];
	}
	/*   viterbi  */
	if(rCoupled[0] == 1) {
	  add_constant_Sequence(genome, index, n, &CompactSeq,
				arraySeqRefs, *write_seq, viterbi_counts);
	} else if (rCoupled[0] >1) {
	  viterbi(y, x, genome, index, &rCoupled[0], n, 
		  NewMuCoupled, NewSigma2Coupled, NewBetaCoupled, NewStatCoupled,
		  states, &CompactSeq, arraySeqRefs, *write_seq, viterbi_counts);
	    
	  for (i=0; i < rCoupled[0]-1; ++i) {
	    for (j=0; j < *n; ++j) {
	      if(states[j]==i+1) {
		probStates[(*n * (i + (rCoupled[0]-1) * (rCoupled[0]-2) / 2)) + j] =
		  (1.0 / (double)(times[rCoupled[0]-1]+1)) +
		  probStates[(*n * (i + (rCoupled[0]-1) * (rCoupled[0]-2) / 2)) + j] *
		  (double)(times[rCoupled[0]-1]) /
		  (double)(times[rCoupled[0]-1]+1);
	      }
	      else {
		probStates[(*n * (i + (rCoupled[0]-1) * (rCoupled[0]-2) / 2)) + j] =
		  probStates[(*n * (i + (rCoupled[0]-1) * (rCoupled[0]-2) / 2)) + j] *
		  (double)(times[rCoupled[0]-1]) /
		  (double)(times[rCoupled[0]-1]+1);
	      }
	    }
	  }
	}
	for (i=0; i<rCoupled[0]; ++i) {
	  mu[indexMu[rCoupled[0]-1] + (times[rCoupled[0]-1]*rCoupled[0]) + i] = NewMuCoupled[i];
	  sigma2[indexMu[rCoupled[0]-1] + (times[rCoupled[0]-1]*rCoupled[0]) + i] = NewSigma2Coupled[i];
	}
	for (i=0; i<rCoupled[0]*rCoupled[0]; ++i) {
	  beta[indexBeta[rCoupled[0]-1] + (times[rCoupled[0]-1]*rCoupled[0]*rCoupled[0]) + i] = NewBetaCoupled[i];
	}
	loglik[*TOT * (rCoupled[0]-1) + times[rCoupled[0]-1]] = 
	  loglikLastCoupled[rCoupled[0]-1];
	times[rCoupled[0]-1]++;
      }
      if (*triedChangeDim) {
	k[3*t] = rCoupled[0]; 
	*triedChangeDim = 0;
      }
    }
    /* End exchange move */
  }

  if (*write_seq){
    //Include  k counts in sequence file so we can compute probs. later.
    for(int ii = 0; ii < (2 * ((*TOT) - 1) + 1); ii++) {
      if(k[ii] > 0) {
	k_sum[k[ii] - 1]++;
	Total_k++;
      }
    }

    if((*genome) > 1) {
      char *filenames[(*genome)];
      filenames[0] = strtok(*filename, "\n");
      for(int gg = 1; gg < (*genome); gg++) {
	filenames[gg] = strtok(NULL, "\n");
      }
      for(int gg = 0; gg < (*genome); gg++) {
	num_sequences[gg] = 0;
	Sequence_to_gzfile(filenames[gg], arraySeqRefs[gg], num_sequences + gg,
			   viterbi_counts, k_sum, Total_k);
      }
    } else {
      *num_sequences = 0;
      Sequence_to_gzfile(*filename, CompactSeq, num_sequences,
			 viterbi_counts, k_sum, Total_k);
    }
  }


  Free(muCoupled);
  Free(sigma2Coupled);
  Free(betaCoupled);
  Free(loglikLastCoupled);
  Free(OldMuCoupled);
  Free(OldSigma2Coupled);
  Free(OldBetaCoupled);
  Free(OldStatCoupled);
  Free(qCoupled);
  Free(rCoupled);
  Free(heat);
  Free(states);
  Free(accepted);
  Free(triedChangeDim);
  Free(q);
  Free(indexMu);
  Free(indexBeta);
  Free(indexStat);
  Free(NewStatCoupled);
  Free(NewMuCoupled);
  Free(NewSigma2Coupled);
  Free(NewBetaCoupled);
  Free(indexPerm);
  
  PutRNGstate();
}

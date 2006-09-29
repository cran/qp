/*
  qp package - this R code implements part of the R qp package described
  in R. Castelo and A. Roverato. A robust procedure for Gaussian graphical
  model search from microarray data with p larger than n, Journal of Machine
  Learning Research, accepted for publication.
 
  Copyright (C) 2006 R. Castelo and A. Roverato
  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
 
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 
  You should have received a copy of the GNU General Public License
  along with this program; if not, you can obtain one via WWW at
  http://www.gnu.org/copyleft/gpl.html, or by writing to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/



#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Lapack.h>
#include <R_ext/RS.h>

/* datatype definitions */

typedef struct {
  unsigned int* e;
  unsigned int  n;
} uin_vector_t;

typedef struct tag_clique_t {
  unsigned int* vtc;
  unsigned int n;
  struct tag_clique_t* next;
} clique_t;
                                                                                                
typedef struct {
  clique_t* first;
  int       n;
} clique_set_t;

/* function prototypes */

uin_vector_t
new_uin_vector(int);
                                                                                                
void
destroy_uin_vector(uin_vector_t);

void
add_clique(clique_set_t* cset, uin_vector_t* clq);
                                                                                                
void
destroy_cliques(clique_set_t* cset);

#ifdef Win32
extern void R_ProcessEvents(void);
extern void R_FlushConsole(void);
#endif

static SEXP
qp_fast_search(SEXP S, SEXP N, SEXP q, SEXP T, SEXP significance);

static SEXP
qp_fast_edge_prob(SEXP S, SEXP NR, SEXP iR, SEXP jR, SEXP qR, SEXP TR, SEXP sigR);

static int
qp_edge_prob(double* S, int n_var, int N, int i, int j, int q, int T, double significance);

static void
sample_d_con_N_loop(int N_loop, int d_con, int v_i, int v_j, int n, int* y);

static SEXP
qp_fast_ci_test(SEXP S, SEXP NR, SEXP iR, SEXP jR, SEXP C);

static double
qp_ci_test(double* S, int n_var, int N, int i, int j, int* C, int q);

static SEXP
qp_fast_get_cliques(SEXP I);

void
cliques_extend(unsigned int n, unsigned int *imat, uin_vector_t* compsub, int* c,
               int* pos, uin_vector_t old, int ne, int ce, clique_set_t* clqlst);

static void
matprod(double *x, int nrx, int ncx, double *y, int nry, int ncy, double *z);
 
static void
matinv(double* inv, double* M, int n);


/* R-function register */

static R_CallMethodDef
callMethods[] = {
  {"qp_fast_search", (DL_FUNC) &qp_fast_search, 5},
  {"qp_fast_edge_prob", (DL_FUNC) &qp_fast_edge_prob, 7},
  {"qp_fast_ci_test", (DL_FUNC) &qp_fast_ci_test,5},
  {"qp_fast_get_cliques", (DL_FUNC) &qp_fast_get_cliques, 1},
  {NULL}
};

void
R_init_qp(DllInfo* info) {

  R_registerRoutines(info,NULL,callMethods,NULL,0);

  GetRNGstate(); /* initialize the R-builtin RNG */

}



/*
  FUNCTION: qp_fast_search
  PURPOSE: compute for each pair of vertices indexed by the rows (columns)
           of the matrix S, the number of times it would be removed out of
           S trials according to the variance-covariance matrix S and
           the p-value threshold thr_p
  RETURNS: matrix of counts of edge removal (null hypothesis acceptance) for
           each adjacency of the undirected graph
*/

static SEXP
qp_fast_search(SEXP S, SEXP N, SEXP q, SEXP T, SEXP significance) {
  int    n_var;
  SEXP   G;
  int    i,j,k,t,pct,ppct;

  PROTECT_INDEX Spi;

  PROTECT_WITH_INDEX(S,&Spi);

  /* number of variables equals number of rows */
  n_var = INTEGER(getAttrib(S,R_DimSymbol))[0];

  if (INTEGER(q)[0] < 0 || INTEGER(q)[0] > n_var-2)
    error("q=%d > n.var-2=%d",INTEGER(q)[0],n_var-2);

  REPROTECT(S = coerceVector(S,REALSXP),Spi);
  PROTECT(G   = allocMatrix(INTSXP,n_var,n_var));

  ppct = -1;
  t = n_var*(n_var-1)/2;
  k = 0;
  for (i = 0; i < n_var-1; i++) {
    for (j = 0; j <= i; j++)
      INTEGER(G)[i+j*n_var] = 0;
    for (j = i+1; j < n_var; j++) {
      INTEGER(G)[i+j*n_var] = qp_edge_prob(REAL(S),n_var,INTEGER(N)[0],i,j,INTEGER(q)[0],
                                           INTEGER(T)[0],REAL(significance)[0]);
      k++;
      pct = (int) ((k*100)/t);
      if (pct != ppct) {
        if (pct % 10 == 0)
          Rprintf("%d",pct);
        else
          Rprintf(".",pct);
#ifdef Win32
        R_FlushConsole();
        R_ProcessEvents();
#endif
        ppct = pct;
      }
    }
  }
  Rprintf("\n");

  for (j = 0; j < n_var; j++)
    INTEGER(G)[n_var-1+j*n_var] = 0;

  UNPROTECT(2);   /* S G */

  return G;
}



/*
  FUNCTION: qp_fast_ci_test
  PURPOSE: wrapper of the R-C interface for calling the function that performs
           a test for conditional independence between variables i and j
           given de conditioning set C
  RETURNS: a list with two members, the t-statistic value and the p-value
           on rejecting the null hypothesis of independence
*/

static SEXP
qp_fast_ci_test(SEXP S, SEXP NR, SEXP iR, SEXP jR, SEXP C) {
  int    N = INTEGER(NR)[0];
  int    n_var = INTEGER(getAttrib(S,R_DimSymbol))[0];
  int    q;
  int*   cond;
  int    i,j,k;
  double p_value;
  double t_value;
  SEXP   result;
  SEXP   result_names;
  SEXP   result_t_val;
  SEXP   result_p_val;

  PROTECT_INDEX Spi,Cpi;

  PROTECT_WITH_INDEX(S,&Spi);
  PROTECT_WITH_INDEX(C,&Cpi);

  REPROTECT(S = coerceVector(S,REALSXP),Spi);
  REPROTECT(C = coerceVector(C,INTSXP),Cpi);

  i = INTEGER(iR)[0] - 1;
  j = INTEGER(jR)[0] - 1;
  q = length(C);

  cond = Calloc(q, int);
  for (k=0;k<q;k++)
    cond[k] = INTEGER(C)[k]-1;

  t_value = qp_ci_test(REAL(S),n_var,N,i,j,cond,q);
  p_value = 2.0 * (1.0 - pt(fabs(t_value),N-q-2,1,0));

  PROTECT(result = allocVector(VECSXP,2));
  SET_VECTOR_ELT(result,0,result_t_val = allocVector(REALSXP,1));
  SET_VECTOR_ELT(result,1,result_p_val = allocVector(REALSXP,1));
  PROTECT(result_names = allocVector(STRSXP,2));
  SET_STRING_ELT(result_names,0,mkChar("t.value"));
  SET_STRING_ELT(result_names,1,mkChar("p.value"));
  setAttrib(result,R_NamesSymbol,result_names);
  REAL(VECTOR_ELT(result,0))[0] = t_value;
  REAL(VECTOR_ELT(result,1))[0] = p_value;

  UNPROTECT(4); /* S C result result_names */

  Free(cond);

  return result;
}



/*
  FUNCTION: qp_ci_test
  PURPOSE: perform a test for conditional independence between variables
           indexed by i and j given the conditioning set Q
  RETURNS: a list with two members, the t-statistic value and the p-value
           on rejecting the null hypothesis of independence
*/

static double
qp_ci_test(double* S, int n_var, int N, int i, int j, int* Q, int q) {
  int*    subvars;
  int     subn = q + 2;
  int     k,l;
  double* Mmar;
  double  S11;
  double* S12;
  double* S21;
  double* S22;
  double* S22inv;
  double* S22inv1col;
  double* tmpmat;
  double  tmpval;
  double  betahat;
  double  sigma;
  double  se;
  double  t_value;

  subvars     = Calloc(subn,int);
  Mmar        = Calloc(subn*subn,double);
  S12         = Calloc(subn,double);
  S21         = Calloc(subn,double);
  S22         = Calloc((subn-1)*(subn-1),double);
  S22inv      = Calloc((subn-1)*(subn-1),double);
  S22inv1col  = Calloc(subn-1,double);

  subvars[0] = i; /* order here is important, first variable i */
  subvars[1] = j; /* then variable j then the conditioning set */
  for (k=2;k<subn;k++)
    subvars[k] = Q[k-2];

  /* Mmar <- S[c(i, j, sp), c(i, j, sp)] 
     S11     <- Mmar[1,1]
     S12     <- Mmar[1,-1]
     S21     <- Mmar[-1,1]
     S22     <- Mmar[-1,-1] */
  for (k=0;k<subn;k++)
    for (l=0;l<subn;l++) {
      Mmar[k+l*subn] = S[subvars[k]+subvars[l]*n_var];
      if (k == 0 && l > 0)
        S12[l-1] = Mmar[k+l*subn];
      if (k > 0 && l == 0)
        S21[k-1] = Mmar[k+l*subn];
      if (k > 0 && l > 0)
        S22[k-1+(l-1)*(subn-1)] = Mmar[k+l*subn];
    }
  S11 = Mmar[0];

  /* S22inv  <- solve(S22) */
  matinv(S22inv,S22,subn-1);

  /* betahat <- S12 %*% S22inv[,1] */
  Memcpy(S22inv1col,S22inv,(size_t) (subn-1));
  matprod(S12,1,subn-1,S22inv1col,subn-1,1,&betahat);

  /* sigma   <- sqrt((S11 - S12 %*% S22inv %*% S21) * (N - 1) / (N - q - 2)) */
  tmpmat = Calloc(subn-1,double);
  matprod(S22inv,subn-1,subn-1,S21,subn-1,1,tmpmat);
  matprod(S12,1,subn-1,tmpmat,subn-1,1,&tmpval);
  Free(tmpmat);
  sigma = sqrt( (S11 - tmpval) * (N - 1) / (N - subn) );
  /* se      <- sigma * sqrt(S22inv[1,1] / (N - 1)) */
  se = sigma * sqrt(S22inv[0] / (N - 1));
  /* t.value <- betahat / se */
  t_value = betahat / se;

  Free(S22inv1col);
  Free(S22inv);
  Free(S22);
  Free(S21);
  Free(S12);
  Free(Mmar);
  Free(subvars);

  return t_value;
}



/*
  FUNCTION: qp_edge_prob
  PURPOSE: compute the probability of the edge as the number of tests
           that accept the null hypothesis of independence given the
           q-order conditionals
  RETURNS: number of edges that would be removed under the given
           threshold, i.e., number of times the null hypothesis has
           been accepted
*/

static int
qp_edge_prob(double* S, int n_var, int N, int i, int j, int q, int T,
             double significance) {
  double thr;
  int    subn = q + 2;
  int*   q_by_T_samples;
  int    k;
  int    nremedges = 0;

  thr = qt(1.0-(significance/2.0),N-subn,1,0);

  q_by_T_samples = Calloc(q * T, int);

  sample_d_con_N_loop(T,q,i,j,n_var,q_by_T_samples);

  for (k = 0; k < T; k++) {
    double t_value;

    t_value = qp_ci_test(S,n_var,N,i,j,(int*) (q_by_T_samples+k*q),q);

    if (fabs(t_value) < thr)
      nremedges++;
  }

  Free(q_by_T_samples);

  /* return the number of coefficients under the threshold,
     i.e., number of removed edges */

  return nremedges;
}



/*
  FUNCTION: qp_fast_edge_prob
  PURPOSE: wrapper of the R-C interface for calling the function that
           computes the probability of the edge as the number of tests
           that accept the null hypothesis of independence given the
           q-order conditionals
  RETURNS: number of edges that would be removed under the given
           threshold, i.e., number of times the null hypothesis has
           been accepted
*/

static SEXP
qp_fast_edge_prob(SEXP S, SEXP NR, SEXP iR, SEXP jR, SEXP qR, SEXP TR,
                  SEXP sigR) {
  int    i,j;
  int    N;
  int    q;
  int    T;
  int    n_var;
  double significance;
  SEXP   nremedges;

  PROTECT_INDEX Spi;

  PROTECT_WITH_INDEX(S,&Spi);

  REPROTECT(S = coerceVector(S,REALSXP),Spi);

  i = INTEGER(iR)[0] - 1;
  j = INTEGER(jR)[0] - 1;

  N = INTEGER(NR)[0];
  q = INTEGER(qR)[0];
  T = INTEGER(TR)[0];

  significance = REAL(sigR)[0];

  /* number of variables equals number of rows */
  n_var = INTEGER(getAttrib(S,R_DimSymbol))[0];

  if (i < 0 || i > n_var-1 || j < 0 || j > n_var-1)
    error("vertices of the selected edge (i,j) should lie in the range [1,nrow(M)]");

  PROTECT(nremedges = allocVector(INTSXP,1));

  INTEGER(nremedges)[0] = qp_edge_prob(REAL(S),n_var,N,i,j,q,T,significance);

  UNPROTECT(2); /* S nremedges */

  return nremedges;
}



/*
  FUNCTION: sample_d_con_N_loop
  PURPOSE: sample without replacement d_con elements from n, N_loop times. this
           is a re-make of the SampleNoReplace function of random.c specifically
           tailored to sample in one shot all we need
  RETURN: a vector with of the N_loop samples of d_con elements one after each
          other
*/

static void
sample_d_con_N_loop(int N_loop, int d_con, int v_i, int v_j, int n, int* y) {
  int  i;
  int  k;
  int  total_j = 0;
  int* x;
  int* z;

  k = 0;
  x = Calloc(n,int);
  z = Calloc(n,int);

  for (i = 0; i < n; i++) {             /* x is a working-only vector */
    x[i] = i;
    z[i] = i;                           /* maps each vertex into a proper place */
  }

  if (v_i < v_j) {                      /* we should take care that the mapping z   */
    z[v_i] = v_j != n-2 ? n-2 : n-1;    /* re-maps the v_i and v_j vertices to the  */
    z[v_j] = z[v_i] != n-1 ? n-1 : n-2; /* n-1 and n-2 properly when any of the two */
  } else {                              /* is smaller than n-2                      */
    z[v_j] = v_i != n-2 ? n-2 : n-1;
    z[v_i] = z[v_j] != n-1 ? n-1 : n-2;
  }

  for (i = 0; i < N_loop; i++) {
    int j;
    int m = n-2;                              /* we sample from n-2 elements */

    for (j = 0; j < d_con ; j++) {
      int r;

      r = (int) (((double) m) * unif_rand()); /* sample using R-builtin RNG */
      y[total_j + j] = x[r];
      x[r] = x[--m];                          /* sample without replacement */

    }

    for (j = total_j; j < total_j+d_con; j++) { /* replace again the sampled elements */
      x[y[j]] = y[j];                           /* for the next round of N_loop       */
      y[j] = z[y[j]];                           /* use the mapping z to avoid choosing v_i or v_j */
    }

    total_j += d_con;

  }

  Free(x);
  Free(z);
}



/*
  FUNCTION: new_uin_vector
  PURPOSE: create a new object of the type uin_vector_t
  RETURN: the created object
*/

uin_vector_t
new_uin_vector(int n) {
  uin_vector_t v;
                                                                                                
  v.e = Calloc(n,unsigned int);
  v.n = n;
                                                                                                
  return v;
}
                                                                                                


/*
  FUNCTION: destroy_uin_vector
  PURPOSE: destroys an object of the type uin_vector_t which
           in this case is just freeing the memory allocated
           for the elements of a vector
  RETURN: the created object
*/

void
destroy_uin_vector(uin_vector_t v) {
  if (v.e != NULL)
    Free(v.e);
}



/*
  FUNCTION: add_clique
  PURPOSE: add a clique of vertices, stored in a vector, to a
           set of cliques stored as an object of the type clique_set_t
  RETURN: none
*/

void
add_clique(clique_set_t* cset, uin_vector_t* clq) {
  clique_t* c;
  int       i;
                                                                                                
  c = Calloc(1,clique_t);
  c->next = NULL;
                                                                                                
  if (cset->n == 0) {
    cset->first = c;
  } else {
    clique_t* p;
                                                                                                
    p = cset->first;
    while (p->next != NULL)
      p = p->next;
    p->next = c;
  }
                                                                                                
  c->vtc = Calloc(clq->n,unsigned int);

  for (i=0;i<clq->n;i++)
    c->vtc[i] = clq->e[i];
  c->n = clq->n;
                                                                                                
  cset->n++;
}



/*
  FUNCTION: destroy_cliques
  PURPOSE: destroys an object of the type clique_set_t which
           consists of going through a dynamically linked list
           and freeing the memory allocated for each of its elements
  RETURN: none
*/

void
destroy_cliques(clique_set_t* cset) {
  clique_t* p;
                                                                                                
  if (cset->n == 0)
    return;
                                                                                                
  p = cset->first;
  while (p != NULL) {
    clique_t* tmp;
                                                                                                
    tmp = p->next;
    Free(p->vtc);
    Free(p);
    p = tmp;
  }

  cset->n = 0;
}



/*
  FUNCTION: qp_fast_get_cliques
  PURPOSE: finds the (maximal) cliques of an undirected graph, it implements
           the algorithm from

           Bron, C. and Kerbosch, J.
           Finding all cliques of an undirected graph
           Communications of the ACM, 16(9):575-577, 1973.
  RETURNS: a list of (maximal) cliques
*/

static SEXP
qp_fast_get_cliques(SEXP I) {
  int            n = INTEGER(getAttrib(I,R_DimSymbol))[0];
  unsigned int*  imat;
  clique_set_t   clqlst;
  uin_vector_t   all;
  uin_vector_t   compsub;
  int            c = 0;
  int            pos = -1;
  int            i;
  SEXP           clqlstR;

  PROTECT_INDEX Ipi;

  PROTECT_WITH_INDEX(I,&Ipi);

  REPROTECT(I = coerceVector(I,INTSXP),Ipi);

  /* copy the incidence matrix 'I' into a local version 'imat' of the
     imat to put the diagonal into ones */
  imat = Calloc(n*n,unsigned int);

  Memcpy(imat,INTEGER(I),(size_t) (n*n));
  for (i=0;i<n;i++)
    imat[i*n+i] = 1;

  UNPROTECT(1); /* I */

  all = new_uin_vector(n);
  compsub = new_uin_vector(n);
  for (i=0;i<n;i++) {
    all.e[i] = i+1;
    compsub.e[i] = -1;
  }
  compsub.n = 0;
  clqlst.n = 0;
  cliques_extend(n,imat,&compsub,&c,&pos,all,0,n,&clqlst);
                                                                                                
  destroy_uin_vector(all);
  destroy_uin_vector(compsub);
  Free(imat);

  PROTECT(clqlstR = allocVector(VECSXP,clqlst.n));

  if (clqlst.n > 0) {
    clique_t*      p;
    int            iclq;

    iclq = 0;
    p = clqlst.first;
    while (p != NULL) {
      SEXP clq;
      int  i;

      SET_VECTOR_ELT(clqlstR,iclq,clq = allocVector(INTSXP,p->n));
      for (i=0;i<p->n;i++)
        INTEGER(VECTOR_ELT(clqlstR,iclq))[i] = p->vtc[i];

      iclq++;
      p = p->next;
    }

  }

  UNPROTECT(1); /* clqlstR */

  destroy_cliques(&clqlst);

  return clqlstR;
}



/*
  FUNCTION: cliques_extend
  PURPOSE: recursive function called from qp_fast_get_cliques that
           finds the (maximal) cliques of an undirected graph, it
           implements the algorithm from

           Bron, C. and Kerbosch, J.
           Finding all cliques of an undirected graph
           Communications of the ACM, 16(9):575-577, 1973.
  RETURNS: none
*/

void
cliques_extend(unsigned int n, unsigned int *imat, uin_vector_t* compsub, int* c,
               int* pos, uin_vector_t old, int ne, int ce, clique_set_t* clqlst) {
  uin_vector_t new;
  unsigned int minnod;
  unsigned int fixp,nod;
  unsigned int i,j;
  unsigned int count;
  unsigned int newne,newce;
  unsigned int p,s,sel;

  new = new_uin_vector(ce);
  for (i=0;i<ce;i++)
    new.e[i] = i+1;

  minnod = ce;
  nod    = 0;
  i      = 1;
  fixp   = -999; /* this value should actually never be used !! */
  s      = -999; /* this value should actually never be used !! */

  while (i <= ce && minnod != 0) {
    p = old.e[i-1];
    count = 0;
    j = ne; j++;
    while (j <= ce && count < minnod) {
      if (!imat[(old.e[j-1]-1)*n+p-1]) {
        count++;
        *pos = j;
      }
      j++;
    }

    if (count < minnod) {
      fixp = p;
      minnod = count;
      if (i <= ne)
        s = *pos;
      else {
        s = i;
        nod = 1;
      }
    }
    i++;
  }

  nod = minnod + nod;
  while (nod >= 1) {
    p = old.e[s-1];
    old.e[s-1] = old.e[ne+1-1];
    sel = old.e[ne+1-1] = p;
    newne = 0;
    i = 1;
    while (i <= ne) {
      if (imat[(old.e[i-1]-1)*n+sel-1]) {
        newne++;
        new.e[newne-1] = old.e[i-1];
      }
      i++;
    }
    newce = newne;
    i = ne + 1;
    i++;
    while (i <= ce) {
      if (imat[(old.e[i-1]-1)*n+sel-1]) {
        newce++;
        new.e[newce-1] = old.e[i-1];
      }
      i++;
    }
    (*c)++;
    compsub->e[*c-1] = sel;
    compsub->n++;
    if (newce == 0)
      add_clique(clqlst,compsub);
    else if (newne < newce)
      cliques_extend(n,imat,compsub,c,pos,new,newne,newce,clqlst);
    (*c)--;
    compsub->n--;
    ne++;
    if (nod > 1) {
      s = ne;
      do s++;
      while (imat[(old.e[s-1]-1)*n+fixp-1]);
    }
    nod--;
  }

  destroy_uin_vector(new);
}



/*
  FUNCTION: matprod
  PURPOSE: multiply two matrices by using the LaPACK library that
           comes along with the R distribution, this code is taken from

           R-2.2.0/src/main/array.c
  RETURNS: none
*/

static void
matprod(double *x, int nrx, int ncx, double *y, int nry, int ncy, double *z) {
    char *transa = "N", *transb = "N";
    int i,  j, k;
    double one = 1.0, zero = 0.0, sum;
    Rboolean have_na = FALSE;

    if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0) {
        /* Don't trust the BLAS to handle NA/NaNs correctly: PR#4582
         * The test is only O(n) here
         */
        for (i = 0; i < nrx*ncx; i++)
            if (ISNAN(x[i])) {have_na = TRUE; break;}
        if (!have_na)
            for (i = 0; i < nry*ncy; i++)
                if (ISNAN(y[i])) {have_na = TRUE; break;}
        if (have_na) {
            for (i = 0; i < nrx; i++)
                for (k = 0; k < ncy; k++) {
                    sum = 0.0;
                    for (j = 0; j < ncx; j++)
                        sum += x[i + j * nrx] * y[j + k * nry];
                    z[i + k * nrx] = sum;
                }
        } else
            F77_CALL(dgemm)(transa, transb, &nrx, &ncy, &ncx, &one,
                            x, &nrx, y, &nry, &zero, z, &nrx);
    } else /* zero-extent operations should return zeroes */
        for(i = 0; i < nrx*ncy; i++) z[i] = 0;
}



/*
  FUNCTION: matinv
  PURPOSE: calculates de inverse of a matrix by using the LaPACK library
           that comes along with the R distribution, this code is taken from
           the function modLa_dgesv in file

           R-2.2.0/src/modules/lapack/Lapack.c
  RETURNS: none
*/

static void
matinv(double* inv, double* M, int n) {
  int     i,j;
  int     info;
  int*    ipiv;
  double* avals;
  double* work;
  double  anorm;
  double  rcond;
  double  tol = DBL_MIN;

  for (i=0;i<n;i++)
    for (j=0;j<n;j++)
      inv[i+j*n] = i == j ? 1.0 : 0.0; 

  ipiv = (int *) Calloc(n,double);
  avals = (double *) Calloc(n*n,double);
  Memcpy(avals,M,(size_t) (n*n));

  F77_CALL(dgesv)(&n,&n,avals,&n,ipiv,inv,&n,&info);
  if (info < 0)
    error("argument %d of Lapack routine %s had invalid value",-info, "dgesv");
  if (info > 0)
    error("Lapack routine dgesv: system is exactly singular");

  anorm = F77_CALL(dlange)("1", &n, &n, M, &n, (double*) NULL);

  work = (double *) Calloc(4*n,double);

  F77_CALL(dgecon)("1", &n, avals, &n, &anorm, &rcond, work, ipiv, &info);
  if (rcond < tol)
    error("system is computationally singular: reciprocal condition number = %g",rcond);

  Free(ipiv);
  Free(avals);
  Free(work);
}

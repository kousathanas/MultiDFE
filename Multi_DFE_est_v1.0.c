#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_linalg.h>
#include <getopt.h>
#include <string.h>

#include "libhdr"
#include "Library_DFE_v1.0.h"

/***
    version 1.0. Uploaded to bitbucket

    compile with:
    gcc -O3 -o Multi_DFE *.c -lm -lgsl -lgslcblas -w

***/


/*Function List*/
double  likelihood_function (const gsl_vector *v, void *params); 
double  likelihood_function_neutral (const gsl_vector *v, void *params); 
double calculate_lik(int N,float *p,double *FV0,double *FVS);
double calculate_lik_neutral(int N,float *p,double *FV0);
double simplex(float *par,double *in,double *out,double size1,double maxsize,int maxiter,int npar,char *function);
double my_minimization(double n2,float *par,int ranrep,gsl_rng *rgen1);
double my_minimization_neutral(double n2,float *par);
int my_N2_minimization(float *par,int gss_lower, int gss_upper);
double initialize_parameters_for_maximization(gsl_rng *rgen1,int number_routine,double *routine_parameters,int number_s,int number_prob,double *out,int ranrep,double n2);
void calculate_weighted_FVS(int n2,int nspikes, double *prob_vec, double **FVSX, double *FVS);
double calculate_FV(int model,int n1,int n2,int t2,double s0,double s,double beta_scaler,double mean_s,double f0,double *mean_FV);
double compute_mean_stats(double *s_evaluated_vec,int n_s_evaluated);
double compute_prop_muts(double lower, double upper);
double compute_fix_prob(double *s_evaluated_vec,int n_s_evaluated);
double jukescantor(double sites, int diffs);
double finalsize;

/*Global variable List*/
double finals,MAX_f0,finalLik=0;
double saveresults[100];
int npar=0;
int nspikes=0;
int selmode=0;
int sfsfold=0;
int MAX_t,finalN2=0;

int N1=100;//global parameter for N1
int en2=0;//global parameter for N2
double tau;
int nalleles=0;
int n2_step=0, conpop=0, output_egf_mod=0;

typedef struct{
  double *sval,*sprob;
  double alpha,beta,es,beta_scaler;
  double freq_zero,t_real;
  double maxlik;
  int N2;
  int actualreps;
  double prop[5];
  double *preFV0,*preFVS,*FV0,*FVS;
  double lethal_prop;
  double ne;
  double fix_prob;
  double mean_s,mean_s2,mean_h;
}save_mle_struct;

typedef struct{
  double s,s2,h;
}mean_stats_struct;

save_mle_struct save_mle;
mean_stats_struct mean_stats;
/******************************************************************************/
double  likelihood_function (const gsl_vector *v, void *params){

  int n2=en2;//set n2 to global parameter en2
  int i=0,j=0;
  int n2d=2*n2;
  int checker=0;
  float *p = (float *)params;

  double lik=0;
  double f0=0,t=0,mean_s=0,sum_prob=0;
  double gamma_beta=0;
  double beta_alpha=0,beta_beta=0;
  double beta_scaler=0;
  double logmean=0,logsigma=0;
  double mean_s_steps=0;
  f0 = gsl_vector_get(v, 0);

  if (!conpop){
    t=tau;
  }
  else{
    t=100;
  }
 
  double * spikes_vec = (double*) calloc (npar+1, sizeof(double));
  double * prob_vec = (double*) calloc (npar+1, sizeof(double));
  double * FV0 = (double*) calloc (maxnd+1, sizeof(double));
  double * FVS = (double*) calloc (maxnd+1, sizeof(double));
  double **FVSX= calloc(nspikes+1, sizeof(double *));

  /*calculate n_e*/  
  double n_e=0;
  n_e=calculate_ne(N1,n2,t);


  /*allocate memory for AFVs*/
  if (nspikes>0){
    for(i = 0; i <= nspikes; i++){
      FVSX[i] = calloc(maxnd+1,  sizeof(double));
    }
  }
        
  if (f0<0.0||f0>1.0){
    lik=undefined;
    goto End;
  }

  /*select selection mode*/
  switch (selmode) {
   
  case 0:case 1:/*spike model OR step model*/
    /*assign s values*/
    spikes_vec[0]=0.0;

    for(i = 1; i <= nspikes; i++)	{
      spikes_vec[i]=gsl_vector_get(v, i);

      if (spikes_vec[i]>0.0){
	lik=undefined;
	goto End;
      }
      if (spikes_vec[i-1]<spikes_vec[i]){
	lik=undefined;
	goto End;
      }
      //if (spikes_vec[i]/spikes_vec[i-1]>100){lik=undefined;goto End;} 
      //if (spikes_vec[i]*n_e<-2*n_e){lik=undefined;goto End;} 

    }

    //monitorinput();
    /*steps start from 0*/
    /*calculate FV0*/
    checker=calculate_FV(selmode,N1,n2,t,0.0,0.0,0.0,0.0,f0,FV0);
    if (checker==undefined){
      lik=undefined;
      goto End;
    }

    /*assign probability values*/
    sum_prob=0;
    for(i = 1; i < nspikes; i++){
      prob_vec[i]=gsl_vector_get(v, nspikes+i); 
      sum_prob+=prob_vec[i];
      if (prob_vec[i]<0.0){
	lik=undefined;
	goto End;
      }
    }
    prob_vec[nspikes]=1-sum_prob;
	
    if (nspikes==1){
      prob_vec[1]=1.0;
    }
    /*assign last element of probability vector*/
      
        
    if ((sum_prob+prob_vec[nspikes])>1.0){
      lik=undefined;
      goto End;
    }
    if (prob_vec[nspikes]<0.0){
      lik=undefined;
      goto End;
    }

    /*calculate FVSs*/
    for(i = 1; i <= nspikes; i++){
      checker=calculate_FV(selmode,N1,n2,t,spikes_vec[i-1],spikes_vec[i],0.0,spikes_vec[i],f0,FVSX[i]);  
      if (checker==undefined){
	lik=undefined;
	goto End;
      }

    }
   
    for (i=1;i<=nspikes;i++){
      for (j=0;j<n2d;j++){
	FVS[j]+=prob_vec[i] * FVSX[i][j];
      }
    }  
     
    break;
         
  case 2:/*gamma distribution model*/

    mean_s=gsl_vector_get(v,1);
    gamma_beta=gsl_vector_get(v,2);

    if (gamma_beta<0.0001||gamma_beta>100){lik=undefined;goto End;}
    if (mean_s>-1e-7){lik=undefined;goto End;}
  
    checker=calculate_FV(selmode,N1,n2,t,0.0,0.0,0.0,0.0,f0,FV0);
    if (checker==undefined){lik=undefined;goto End;}
  
    checker=calculate_FV(selmode,N1,n2,t,-gamma_beta/mean_s,mean_s,0.0,mean_s,f0,FVS);
    if (checker==undefined){lik=undefined;goto End;}
    break; 
         
  case 3:/*beta distribution model*/ 

    beta_alpha=gsl_vector_get(v,1);
    beta_beta=gsl_vector_get(v,2);
    beta_scaler=gsl_vector_get(v,3);
    mean_s=(beta_alpha/(beta_alpha+beta_beta));
    //mean_s=mean_s*beta_scaler;
      
    if (mean_s>0){mean_s=-mean_s;}
    //if (beta_alpha<0.0001||beta_beta<0.0001){lik=undefined;goto End;}

    //monitorinput();
    checker=calculate_FV(selmode,N1,n2,t,0.0,0.0,0.0,0.0,f0,FV0);
    if (checker==undefined){lik=undefined;goto End;}
    //monitorinput();
    checker=calculate_FV(selmode,N1,n2,t,beta_alpha,mean_s,beta_scaler,mean_s,f0,FVS);
    if (checker==undefined){lik=undefined;goto End;}
    // printf("beta_alpha %f beta_beta %f mean_s %f beta_scaler %f\n",beta_alpha,beta_beta,mean_s,beta_scaler);
    //monitorinput();     
   
    break;
      
  case 4:/*lognormal distribution model*/ 

    logmean=gsl_vector_get(v,1);
    logsigma=gsl_vector_get(v,2);
    mean_s=exp(logmean+pow(logsigma,2)/2);
      
    if (mean_s>0){mean_s=-mean_s;}
    if (logsigma<0){lik=undefined;goto End;}

    checker=calculate_FV(selmode,N1,n2,t,0.0,0.0,0.0,0.0,f0,FV0);
    if (checker==undefined){lik=undefined;goto End;}
  
    checker=calculate_FV(selmode,N1,n2,t,logmean,logsigma,0.0,mean_s,f0,FVS);
    if (checker==undefined){lik=undefined;goto End;}
    //printf("logmean %f logsigma %f mean_s %f\n",logmean,logsigma,mean_s);

   
    break;
      
  case 5://fixed effects mode (Wilson et al 2011;Plos Genetics)
    /*assign s values*/
    spikes_vec[1]=0.0;
    spikes_vec[2]=-1/n_e;
    spikes_vec[3]=-5/n_e;
    spikes_vec[4]=-10/n_e;
    spikes_vec[5]=-50/n_e;
    spikes_vec[6]=-1;

    /*calculate FV0*/
    checker=calculate_FV(selmode,N1,n2,t,0.0,0.0,0.0,0.0,f0,FV0);
    if (checker==undefined){lik=undefined;goto End;}

    /*assign probability values*/
    sum_prob=0;
    for(i = 1; i < nspikes; i++) {
      prob_vec[i]=gsl_vector_get(v, i); 
      sum_prob+=prob_vec[i];
      if (prob_vec[i]<0.0){lik=undefined;goto End;}	   
    }
    prob_vec[nspikes]=1-sum_prob;
        
    if ((sum_prob+prob_vec[nspikes])>1.0){lik=undefined;goto End;}
    if (prob_vec[nspikes]<0.0){lik=undefined;goto End;}

    /*calculate FVSs*/
    for(i = 1; i <= nspikes; i++) {
      checker=calculate_FV(0,N1,n2,t,spikes_vec[i-1],spikes_vec[i],0.0,spikes_vec[i],f0,FVSX[i]);  
      if (checker==undefined){lik=undefined;goto End;}
    }

   
    for (i=1;i<=nspikes;i++){
      for (j=0;j<n2d;j++){
	FVS[j]+=prob_vec[i] * FVSX[i][j];
      }
    }        
      
    break;
      
         
         
  }/*switch end*/////////////////////////////////////////////////////

  /* scale FV0 and FVS vectors */
  egf_scaling(n2,f0,FV0,FVS);
  egf_scaling(n2,f0,FV0,FV0);

  /*save expected FVS vector*/
  for(i = 0; i <= n2d; i++){  
    save_mle.preFV0[i]=FV0[i];
    save_mle.preFVS[i]=FVS[i];
  }

  /*fold FV0 and FVS vectors*/
  if (sfsfold)
    {
      fold_vector_double(FV0, n2d);
      fold_vector_double(FVS, n2d);
    } 

  //dumpvector(FVS, 0, n2*2, "FVS"); 
  
  /*calculate likelihood*/
  lik=calculate_lik(n2,p,FV0,FVS);

  /* FREE VECTORS */
 End:

  if (nspikes>0)
    {
      for(i = 0; i <= nspikes; i++){
	free(FVSX[i]);
      }
    }

  free(FV0);
  free(FVS);
  free(FVSX);
  free(spikes_vec);
  free(prob_vec);
  
  if (isnan(lik))
    {
      lik=undefined;
    }
  //printf("%f\t%f\t%d\n",lik,mean_s,n2);

  //monitorinput();

  return -lik;

}
/******************************************************************************/
double  likelihood_function_neutral (const gsl_vector *v, void *params){

  int n2=en2;//set n2 to global parameter en2
  int i=0,j=0;
  int n2d=2*n2;
  int checker=0;
  int *p = (int *)params;
  double lik=0;
  double f0=0,t=0,mean_s=0,sum_prob=0;
  double * FV0 = (double*) calloc (maxnd+1, sizeof(double));
  
  f0 = gsl_vector_get(v, 0);
  
  if (!conpop){
    t=gsl_vector_get(v, 1);
  }
  else{
    t=100;
  }

  checker=calculate_FV(selmode,N1,n2,t,0.0,0.0,0.0,0.0,f0,FV0);

  if (checker==undefined){lik=undefined;goto End;}
  /* scale FV0 and FVS vectors */
  egf_scaling(n2,f0,FV0,FV0);

  /*fold FV0 and FVS vectors*/
  if (sfsfold){
    fold_vector_double(FV0, n2d);
  } 

  /*calculate likelihood*/
  lik=calculate_lik_neutral(n2,p,FV0);

  /* FREE VECTORS */
 End:
  free(FV0);
  if (isnan(lik)){
    lik=undefined;
  }
  //printf("%d\n",lik);
  //monitorinput();

  return -lik;

}
/******************************************************************************/
double calculate_lik(int N,float *p,double *FV0,double *FVS){

  int i=0,j=0;
  double binomial_prob=0,sel_lik=0,neu_lik=0,prob=0,total_lik=0;
  int n2d=N*2;

  //dumpvector(FVS, 0, n2d, "FVSobs"); 
  //monitorinput();
  if (sfsfold){/*FOLDED*/
    for (i=0;i<=nalleles/2;i++){
      //printf("%d\n",nalleles/2);
	  
      sel_lik=0;
      neu_lik=0;

      for (j=0;j<=N;j++){
	prob=(double)j/n2d;
	if (!odd(nalleles)&&(i==nalleles/2)){
	  binomial_prob=gsl_ran_binomial_pdf (i,prob,nalleles);
	}
	else {
	  binomial_prob=gsl_ran_binomial_pdf (i,prob,nalleles)+ gsl_ran_binomial_pdf (nalleles-i,prob,nalleles);
	}

	sel_lik+=FVS[j]*binomial_prob;
	neu_lik+=FV0[j]*binomial_prob;
      }
      total_lik+=p[i]*log(sel_lik);
      total_lik+=p[i+nalleles]*log(neu_lik);
	  
      //printf("%f\t%f\n",p[i],p[i+nalleles]);
	  
    }
    // printf("OK");
    //monitorinput();	
  }else {/*NOT FOLDED*/

    for (i=0;i<nalleles;i++){
      sel_lik=0;
      neu_lik=0;
	
      for (j=0;j<=n2d;j++){
	prob=(double)j/n2d;

	binomial_prob=gsl_ran_binomial_pdf (i,prob,nalleles);
	sel_lik+=FVS[j]*binomial_prob;
	neu_lik+=FV0[j]*binomial_prob;
      }
      total_lik+=p[i]*log(sel_lik);
      total_lik+=p[i+nalleles]*log(neu_lik);
    }
  }//alleles

  return total_lik;
}
/******************************************************************************/
double calculate_lik_neutral(int N,float *p,double *FV0){

  int i=0,j=0;
  double binomial_prob=0,sel_lik=0,neu_lik=0,prob=0,total_lik=0;
  int n2d=N*2;

  if (sfsfold){ /*FOLDED*/
    for (i=0;i<=nalleles/2;i++){
      //printf("%d\n",nalleles/2);
	  
      sel_lik=0;
      neu_lik=0;

      for (j=0;j<=N;j++){
	prob=(double)j/n2d;
	if (!odd(nalleles)&&(i==nalleles/2)) {
	  binomial_prob=gsl_ran_binomial_pdf (i,prob,nalleles);
	}
	else {
	  binomial_prob=gsl_ran_binomial_pdf (i,prob,nalleles)+ gsl_ran_binomial_pdf (nalleles-i,prob,nalleles);
	}
	neu_lik+=FV0[j]*binomial_prob;
      }
      total_lik+=p[i+nalleles]*log(neu_lik);
	  
      //printf("%d\t%f\n",p[i+nalleles],p[i+nalleles]*log(neu_lik));
	  
    }
    // printf("OK");
    // monitorinput();	
  }else /*NOT FOLDED*/
    {
      for (i=0;i<nalleles;i++){
	sel_lik=0;
	neu_lik=0;
	
	for (j=0;j<=n2d;j++){
	  prob=(double)j/n2d;

	  binomial_prob=gsl_ran_binomial_pdf (i,prob,nalleles);
	  neu_lik+=FV0[j]*binomial_prob;
	}
	total_lik+=p[i+nalleles]*log(neu_lik);
      }

    } /*alleles*/

  return total_lik;
}

/******************************************************************************/
double simplex(float *par,double *in,double *out,double size1,double maxsize,int maxiter,int nopar,char *function){
  const gsl_multimin_fminimizer_type *T =
    gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;

  size_t iter = 0;
  int status=0;
  int i=0;
  double size=0;
  int evaluations=0;

  //printf("%f\n",par[0]);
  //printf("%f\n",par[1]);
  //printf("%f\n",par[2]);
  //monitorinput();
  /* Starting point */
  x = gsl_vector_alloc (nopar);
  for (i=0;i<nopar;i++){
    gsl_vector_set (x, i, in[i]);
  }

  /* Set initial step sizes to 1 */
  ss = gsl_vector_alloc (nopar);
  for (i=0;i<nopar;i++){
    gsl_vector_set (ss, i, in[i]*size1);
  }

  /* Initialize method and iterate */
  minex_func.n =nopar;

  if (function=="likelihood_function_neutral"){
    minex_func.f = likelihood_function_neutral;
  }

  if (function=="likelihood_function"){
    minex_func.f = likelihood_function;
  }
  minex_func.params =par;

  s = gsl_multimin_fminimizer_alloc (T, i);
  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

  do {
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);
    if (status){
      break;
    }
    size = gsl_multimin_fminimizer_size (s);
    status = gsl_multimin_test_size (size,maxsize);
    if (status == GSL_SUCCESS)
      {
	//printf ("converged to minimum at\n");

      }

    evaluations+=1;

  }while (status == GSL_CONTINUE && iter < maxiter);

  for (i=0;i<nopar;i++){
    out[i]=gsl_vector_get (s->x, i);
  }

  out[nopar]=s->fval;
  finalLik=s->fval;
  finalsize=size;

  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);

  return out[nopar];
}
/******************************************************************************/
double initialize_parameters_for_maximization(gsl_rng *rgen1,int number_routine,
					      double *routine_parameters,int number_s,int number_prob,double *out,int ranrep,double n2)
{

  int i=0,j=0,k=0,check=0;
  int total_par=0;
  double sum_prob=0;
  double * s_val = (double*) calloc ( number_s+5, sizeof(double));
  double power1=0;
  double ne=calculate_ne(N1,n2,tau);
  double norm=0;
  
  if (number_routine>0){
    for (i = 0; i < number_routine ; i++){
      total_par=i;
      out[i]=routine_parameters[i];
      //printf("%f" ,out[i]);
    }
  }
   
  if (number_s>1)
    {     
      check=0;
      while(check==0){
	power1=0;
	for (i = 1; i <= number_s ; i++){
	  norm=gsl_ran_gaussian(rgen1,0.1)+1;
	  s_val[i]=pow(ne, (double) (i/number_s)-norm);
	  //printf("%f,%f,\n", norm,s_val[i]);
	  //monitorinput();
	  s_val[i]=-s_val[i];
	  out[total_par+i]= s_val[i];
	}
	for (i = 1; i < number_s ; i++){
	  if (s_val[i+1]<s_val[i]){
	    check=1;
	  } else{
	    check=0;
	    break;
	  }
	}
      }//while(check==0;)

      total_par+=number_s;
    
    }

  if (number_s==1)
    {
      total_par++;
      do
	{
	  out[total_par]=gsl_rng_uniform (rgen1);
      
	}while (out[total_par]<0.0);
    
      out[total_par]=-out[total_par];

    }   
     
  if (number_prob>0)
    {   
      sum_prob=0;
      while (sum_prob==0.0||sum_prob>1.0)
	{
	  sum_prob=0;

	  for(i = 1; i <= number_prob; i++)
	    {
	      out[total_par+i]= 1.0/nspikes;
	      sum_prob+=out[total_par+i];
	      //printf("%f\n",out[total_par+i]);
	      //monitorinput();
	    }

	}
      total_par+=number_prob;
    }
  free(s_val);
          
  if ((total_par+1)!=npar){printf("something went wrong with the initialisation of the starting values\n");}

  return(1); 

}
/******************************************************************************/
double get_max_lik_results(int npar,double *results,int number_s,int number_prob){
  int i=0,j=0,params=0;
  double sum=0;
  double n_e=0;
  /*save f0*/
  save_mle.freq_zero=results[params];
  params++;
  /*save t2, if not a constant population*/
  if (!conpop)
    {
      save_mle.t_real=tau;
      n_e=calculate_ne(N1,en2,tau);
    }else{
    save_mle.t_real=100;
    n_e=100;
  }
  /*select s mode*/
  switch (selmode)
    {
    case 0:case 1:
      save_mle.sval=calloc (nspikes+1, sizeof(double));
      save_mle.sprob=calloc (nspikes+1, sizeof(double));

      save_mle.sval[0]=0;
      for (i=1; i<=number_s; i++){
	save_mle.sval[i]=results[params];
	params++;
      }
      sum=0;
      for (i=1; i<nspikes; i++){
	save_mle.sprob[i]=results[params];
	sum+=results[params];
	params++;
      }
      save_mle.sprob[nspikes]=1-sum;
      break;

    case 2:
      save_mle.es=results[params];params++;
      save_mle.beta=results[params];params++;
      save_mle.alpha=-(save_mle.beta/save_mle.es);
      break;

    case 3:
      save_mle.alpha=results[params];params++;
      save_mle.beta=results[params];params++; 
      save_mle.beta_scaler=results[params];params++; 
      break;

    case 4: 
      save_mle.alpha=results[params];params++;
      save_mle.beta=results[params];params++;  
      break;

    case 5:
      save_mle.sprob=calloc (nspikes+1, sizeof(double));
      save_mle.sval=calloc (nspikes+1, sizeof(double));
      save_mle.sval[1]=0;
      save_mle.sval[2]=-1/n_e;
      save_mle.sval[3]=-5/n_e;
      save_mle.sval[4]=-10/n_e;
      save_mle.sval[5]=-50/n_e;
      save_mle.sval[6]=-1;

      sum=0;
      for (i=1; i<nspikes; i++){
	save_mle.sprob[i]=results[params];
	sum+=results[params];
	params++;
      }
      save_mle.sprob[nspikes]=1-sum;
      break;

    }//end switch

  save_mle.maxlik=-results[npar];
  return(1); 
}
/******************************************************************************/
double my_minimization(double n2,float *par,int ranrep,gsl_rng  *rgen1)
{
  en2=n2; //set global parameter en2 to n2

  int i=0,j=0,z=0;
  int reps=0;
  double max_lik=undefined;
  
  //const gsl_rng_type * T1;
  //T1 = gsl_rng_taus;
  //gsl_rng  *rgen2 =gsl_rng_alloc(T1);
   
  double * in = (double*) calloc ( npar+1, sizeof(double));
  double * out = (double*) calloc (npar+1, sizeof(double));
  double max1=0,max2=0,diff_max=0;
  double size=0;

  /*set up number of "routine" parameters. Can add up more if needed (e.g. a distribution parameter)*/
  int number_routine=0;
  double * routine_parameters = (double*) calloc (1000, sizeof(double));
  save_mle.actualreps=0;
  do {
    number_routine=1;routine_parameters[0]=0.9;/*f0*/
    max1,max2,diff_max=0;
   
    /*select appropriate model*/

    switch (selmode){
    case 0:case 1:
      initialize_parameters_for_maximization(rgen1,number_routine,routine_parameters,nspikes,nspikes-1,in,ranrep,en2);    
      break;

    case 2:case 4:     
      if(selmode==2){routine_parameters[number_routine]=-gsl_rng_uniform(rgen1);number_routine+=1;}//Mean_s	
      if(selmode==4){routine_parameters[number_routine]=-1;number_routine+=1;}//logmean
      routine_parameters[number_routine]=gsl_rng_uniform(rgen1);number_routine+=1;//beta

      initialize_parameters_for_maximization(rgen1,number_routine,routine_parameters,0,0,in,ranrep,en2);
      break;   

    case 3:
      routine_parameters[number_routine]=gsl_rng_uniform(rgen1);number_routine+=1;//beta1
      routine_parameters[number_routine]=gsl_rng_uniform(rgen1);number_routine+=1;
      routine_parameters[number_routine]=1;number_routine+=1;
      initialize_parameters_for_maximization(rgen1,number_routine,routine_parameters,0,0,in,ranrep,en2);	
      break;

    case 5:
      initialize_parameters_for_maximization(rgen1,number_routine,routine_parameters,0,nspikes-1,in,ranrep,en2);
      break;

    }/*switch end*/
      

    /*Likelihood calculations*/
    max1= simplex(par,in,out,0.1,1e-6,250,npar,"likelihood_function");
    max2= simplex(par,out,out,0.05,1e-6,250,npar,"likelihood_function");
    diff_max=max1-max2;
    reps=0;

    //monitorinput();
    while (diff_max>1e-5){
      max1=max2;
      max2= simplex(par,out,out,0.05,1e-6,250,npar,"likelihood_function");
      diff_max=max1-max2;
      reps++;

      if (reps>4){
	break;
      }/*no more than ten calls to the simplex*/
    }

    if (max_lik<-out[npar]){

      for (i=0;i<npar;i++){
	saveresults[i]=out[i];
	//printf("%f\t",saveresults[i]);
      }
      printf("\n");
      get_max_lik_results(npar,out,nspikes,nspikes-1);
      save_mle.N2=n2;
      max_lik=-out[npar];
      for(i = 0; i <= n2*2; i++){ 
	save_mle.FV0[i]=save_mle.preFV0[i];
	save_mle.FVS[i]=save_mle.preFVS[i];
      }
    }
	
    save_mle.actualreps++;
    ranrep--;//repeat counter
  }while(ranrep>0);//end random starts
    
  free(routine_parameters);   
  free(in);
  free(out);
  return max_lik;

}
/******************************************************************************/
double my_minimization_neutral(double n2,float *par)
{

  en2=n2; //set global parameter en2 to n2

  int i=0,j=0,z=0;
  double max_lik=0;
  int npar0=2;
  double * in = (double*) calloc ( npar0+1, sizeof(double));
  double * out = (double*) calloc (npar0+1, sizeof(double));

  in[0]=0.9;
  in[1]=20;

  /*Likelihood calculations*/
  simplex(par,in,out,1,1e-6,250,npar0,"likelihood_function_neutral");
  max_lik=simplex(par,out,out,0.1,1e-6,250,npar0,"likelihood_function_neutral");

  tau=out[1];
 
  free(in);
  free(out);
  
  return -max_lik;
}


#define n_to_evaluate 7
/******************************************************************************/
int my_N2_minimization(float *par,int gss_lower, int gss_upper)
{

  int i=0, x[n_to_evaluate], ind=0, next_ind=0, prev_ind=0, highest_ind=0, n2_save=0, max_found=0,
    prev_n2_evaluated=0,n2=0;
  double lik=0, lik_prev=0, y=0, lik_vec[n_to_evaluate], highest_lik=0, n2_mle=0;
  y = (double)gss_lower;
  ind = nearest_n2_ind(y);
  x[0] = ind;

  for (i=1; i<n_to_evaluate - 1; i++){
    y = (int)(double)(gss_upper - gss_lower)*pow((double)i/(double)n_to_evaluate, 1.8);
    ind = nearest_n2_ind(y);
    x[i] = ind;
  }

  y = (double)gss_upper;
  ind = nearest_n2_ind(y);
  //   printf("y %lf ind %d\n", y, ind);
  x[n_to_evaluate - 1] = ind;
  prev_n2_evaluated = undefined_int;
  for (i=0; i<n_to_evaluate; i++){
    //      printf("x[%d] %d n2_evaluated_vec[x[i]] %d prev_n2_evaluated %d\n",
    //         i, x[i], n2_evaluated_vec[x[i]], prev_n2_evaluated);
    lik_vec[i] = undefined;
    if ((prev_n2_evaluated!=undefined_int)&&(n2_evaluated_vec[x[i]] <= prev_n2_evaluated))
      gabort("Invalid sequence of N2 in forward_search", i);
    prev_n2_evaluated = n2_evaluated_vec[x[i]];
  }
  //   
  lik_prev = undefined;
  for (i=0; i<n_to_evaluate; i++)
    {
      n2 = n2_evaluated_vec[x[i]];
      lik_vec[i] = my_minimization_neutral(n2,par);
      //      printf("lik_vec[%d] %lf\n", i, lik_vec[i]);
      //      
      if (lik_prev != undefined){
	if (lik_vec[i] < lik_prev){
	  break;
	}
      }
      n2_mle = n2;
      //assign_saved_mles(*n2);
      lik_prev = lik_vec[i];
    }
  highest_lik = lik_vec[0];
  highest_ind = 0;
  for (i=0; i<n_to_evaluate; i++){
    n2 = n2_evaluated_vec[x[i]];
    //printf("Index %d N2 %d logL %lf\n", i, n2, lik_vec[i]);
    if ((lik_vec[i] != undefined)&(lik_vec[i] >= highest_lik)){
      highest_lik = lik_vec[i];
      highest_ind = i;
    }
  }
  //printf("Initial search: highest logL %lf, index %d\n", highest_lik, highest_ind);
  // monitorinput();
  //   
  ind = x[highest_ind];
  if (highest_ind < n_to_evaluate-1){
    next_ind = x[highest_ind + 1];
  }
  else{
    next_ind = undefined_int;
  }
  //  printf("ind %d next_ind %d\n", ind, next_ind);
  //   
  n2_save = n2;
  max_found = 0;
  if (next_ind != undefined_int){
    for (i=ind+1; i<next_ind; i++){
      n2 = n2_evaluated_vec[i];
      lik =  my_minimization_neutral(n2,par);
      if (lik < highest_lik){
	n2 = n2_save;
	break;
      }
      else{
	highest_lik = lik;
	n2_save = n2;
	n2_mle = n2;
	//assign_saved_mles(*n2);
	max_found = 1;
      }
    }
  }
  else{
    //      printf("next_ind == undefined_int\n");
  }
  n2 = n2_save;
  //printf("Highest lik %lf n2 %d\n", highest_lik, *n2);
  //   
  if (!max_found) {
    ind = x[highest_ind];
    if (highest_ind > 0){
      prev_ind = x[highest_ind - 1];
    }
    else{
      prev_ind = undefined_int;
    }
    //      printf("ind %d prev_ind %d\n", ind, prev_ind);
    //      
    n2_save = n2;
    if (prev_ind != undefined_int){
      for (i=ind-1; i>prev_ind; i--){
	n2 = n2_evaluated_vec[i];
	lik =  my_minimization_neutral(n2,par);

	if (lik < highest_lik){
	  n2 = n2_save;
	  break;
	}
	else{
	  highest_lik = lik;
	  n2_save = n2;
	  n2_mle = n2;
	  //assign_saved_mles(*n2);
	}
      }
    }
    // printf("Highest lik %lf n2 %d\n", highest_lik, n2);
    //monitorinput();
    //      
  }
  n2 = n2_mle;
  return(n2);
}
/******************************************************************************/
void calculate_weighted_FVS(int n2,int nspikes, double *prob_vec, double **FVSX, double *FVS)
{
  int n2d=n2*2;
  int i=0,j=0;
  for (i=1;i<=nspikes;i++){
    for (j=0;j<n2d;j++){
      FVS[j]+=prob_vec[i] * FVSX[i][j];
    }
  }
}
/******************************************************************************/

double calculate_FV(int model,int n1,int n2,int t2,
		    double s0,double s, double beta_scaler,
		    double mean_s,double f0,double *mean_FV)
{
  
  static double egf_vec1_lower[maxnd+1], egf_vec2_lower[maxnd+1], 
    egf_vec1_upper[maxnd+1], egf_vec2_upper[maxnd+1],
    egf_vec1[maxnd+1], egf_vec2[maxnd+1],egf_vec[maxnd+1];

  double * density_vec = (double*) calloc (maxnd+1, sizeof(double));
 
  int i=0,j=0,file_size_bytes=0;
  int t2_upper=0;
  int t2_real=t2;
  int n1d=2*n1;
  int n2d=2*n2;

  char *buffer_p1_t2_lower, *buffer_p2_t2_lower,
    *buffer_p1_t2_upper, *buffer_p2_t2_upper, *buffer_const_pop;

  if (!conpop){
    get_upper_lower_int(t2_real, &t2_lower, &t2_upper, n_t2_evaluated, t2_evaluated_vec);
    if ((t2_lower==undefined_int)||(t2_upper==undefined_int)){
      return undefined;
    }

    //printf("\n%d",t2_lower);
    //printf("\n%f",s); 
      
    file_size_bytes=compute_file_size_bytes(n2);
    //Read buffers
          
    buffer_p1_t2_lower = (char*) malloc (file_size_bytes);
    buffer_p1_t2_upper = (char*) malloc (file_size_bytes);

    buffer_p2_t2_lower = (char*) malloc (file_size_bytes);
    buffer_p2_t2_upper = (char*) malloc (file_size_bytes);

    read_phase1_phase2_file_into_buffer(n1,1, n2, t2_lower,
					buffer_p1_t2_lower, file_size_bytes);
    read_phase1_phase2_file_into_buffer(n1,1, n2, t2_upper,
					buffer_p1_t2_upper, file_size_bytes);

    read_phase1_phase2_file_into_buffer(n1,2, n2, t2_lower,
					buffer_p2_t2_lower, file_size_bytes);
    read_phase1_phase2_file_into_buffer(n1,2, n2, t2_upper,
					buffer_p2_t2_upper, file_size_bytes);
   
    /* choose selection model*/
    if (mean_s!=0.0){
      /*selection different than 0*/
      
      switch (model){

      case 0:case 5:  
	set_up_density_vec_equal_effects(s, n_s_evaluated, density_vec);
	break;		      
      case 1: 
	set_up_density_vec_step_effects(s0,s,n_s_evaluated, density_vec);
	break;	
      case 2: 
	compute_gamma_densities(s0,-s*s0, s_evaluated_vec, n_s_evaluated, density_vec);
	break;	
      case 3: 
	compute_beta_densities(s0, ((s0*(1+s))/(-s)),beta_scaler,s_evaluated_vec, n_s_evaluated, density_vec);
	break;
      case 4: 
	compute_lognormal_densities(s0, s, s_evaluated_vec, n_s_evaluated, density_vec);
	break;
      } 

      for (i=0; i<n_s_evaluated; i++){
	get_average_egf_vecs_1_and_2(s_evaluated_vec[i],n1, n2, t2_real, t2_lower,
				     t2_upper, egf_vec, buffer_p1_t2_lower, buffer_p2_t2_lower, buffer_p1_t2_upper,
				     buffer_p2_t2_upper);
	for (j=0; j<=n2d; j++){
	  mean_FV[j] += egf_vec[j]*density_vec[i];
	}
      }  
    }
   
    if (s==0.0){
      get_average_egf_vecs_1_and_2(0.0,n1, n2, t2_real, t2_lower,
				   t2_upper, egf_vec, buffer_p1_t2_lower, buffer_p2_t2_lower, buffer_p1_t2_upper,
				   buffer_p2_t2_upper);
      for (j=0; j<=n2d; j++){
	mean_FV[j] = egf_vec[j];
      }
    }
        

    free(buffer_p1_t2_lower);
    free(buffer_p2_t2_lower);
    free(buffer_p1_t2_upper);
    free(buffer_p2_t2_upper);
  }   
  else{
    /*constant population*/
    if (2*n1 > maxnd){
      printf("ERROR: Value of 2*n1 %d exceeds max_nd %d\n", n1, maxnd);
      gabort("Program terminating", 0);
    }
    file_size_bytes = compute_file_size_bytes(n1);
    //      printf("const pop: file_size_bytes %d\n", file_size_bytes); monitorinput();
    buffer_const_pop = (char*) malloc (file_size_bytes);

    read_const_pop_file_into_buffer(n1, buffer_const_pop, file_size_bytes);
      
    if (mean_s!=0.0){
      // choose selection model
      switch (model){

      case 0:case 5:
	set_up_density_vec_equal_effects(s, n_s_evaluated, density_vec);
	break;		      
      case 1: 
	set_up_density_vec_step_effects(s0,s,n_s_evaluated, density_vec);
	break;	
      case 2: 
	compute_gamma_densities(s0, -s*s0, s_evaluated_vec, n_s_evaluated, density_vec);
	break;	
      case 3: 
	compute_beta_densities(s0, ((s0*(1+s))/(-s)),beta_scaler, s_evaluated_vec, n_s_evaluated, density_vec);
	break;
      case 4: 
	compute_lognormal_densities(s0,s,s_evaluated_vec, n_s_evaluated, density_vec );
	break;
      }    
        
      for (i=0; i<n_s_evaluated; i++){
	get_const_pop_egf_vec(s_evaluated_vec[i], n1, egf_vec, buffer_const_pop, 0);
	for (j=0; j<=n1d; j++) mean_FV[j] += egf_vec[j]*density_vec[i];
      }  
    }
   
    if (s==0.0){
      get_const_pop_egf_vec(0.0, n1, egf_vec, buffer_const_pop, 1);
      for (j=0; j<=n1d; j++) mean_FV[j] = egf_vec[j];
    }
    free(buffer_const_pop);
  }

  free(density_vec);

  return (1);
}
/******************************************************************************/
double compute_mean_stats(double *s_evaluated_vec,int n_s_evaluated)
{
  int i=0,j=0;
  double lower=0,upper=0,area=0,total_area=0,inc_x=0,inc_y=0,sum=0,sum2=0,hsum = 0;

  for (i=0; i<n_s_evaluated; i++){
    area=0;
    get_lower_upper(i, s_evaluated_vec, n_s_evaluated, &lower, &upper);
  
    switch(selmode){
    case 0:case 5:
      for (j=1;j<=nspikes;j++){
	if (upper>=-save_mle.sval[j]&&lower<=-save_mle.sval[j]){
	  area+=save_mle.sprob[j];
	}
      }//for end
      break;

    case 1:
      for (j=1;j<=nspikes;j++)
	{   
	  inc_x=gsl_cdf_flat_P(lower,-save_mle.sval[j-1],-save_mle.sval[j]);    
	  inc_y=gsl_cdf_flat_P(upper,-save_mle.sval[j-1],-save_mle.sval[j]);
 
	  area+=(inc_y-inc_x)*save_mle.sprob[j];
	}
      break;

    case 2:
      inc_x = gsl_cdf_gamma_P(lower, save_mle.beta,1/save_mle.alpha);
      inc_y = gsl_cdf_gamma_P(upper, save_mle.beta,1/save_mle.alpha);
      
      area = (inc_y - inc_x);
      break;

    case 3:
      inc_x = gsl_cdf_beta_P (lower, save_mle.alpha, save_mle.beta);
      inc_y = gsl_cdf_beta_P (upper, save_mle.alpha, save_mle.beta);
    
      area = (inc_y - inc_x);
      break;  
    
    case 4:
      inc_x = gsl_cdf_lognormal_P (lower, save_mle.alpha, save_mle.beta);
      inc_y = gsl_cdf_lognormal_P (upper, save_mle.alpha, save_mle.beta);
    
      area = (inc_y - inc_x);
      break; 
    
    }
    sum+=s_evaluated_vec[i]*area;
    sum2+=-pow(s_evaluated_vec[i],2)*area;
    if(s_evaluated_vec[i]!=0){hsum+=(1/s_evaluated_vec[i])*area;}
    
    total_area+=area;

    //printf("%f\t%f\t%f\t%f\t%f\n",lower,upper,area,total_area,s_evaluated_vec[i]);
  }    

  save_mle.mean_s=sum/total_area;
  save_mle.mean_s2=sum2/total_area;
  save_mle.mean_h=(1/hsum)/total_area;
  //printf("%f\t%f\n",total_area,save_mle.mean_s);
  return (1);
}

/******************************************************************************/
double compute_prop_muts(double lower, double upper)
{
  int i=0,n1=0,n2=0;
  double inc_x=0,inc_y=0,area=0,sum=0,t2=0,w1=0,w2=0,n_e=0;
  n1=N1;
  n2=n1;
  t2=100;
  if (!conpop){
    n2=save_mle.N2;
    t2=save_mle.t_real;
    n_e=calculate_ne(n1,save_mle.N2,save_mle.t_real);
  }
  else{
    n_e=n1;
  }

  //printf("%f\n",n_e);
  /*selection mode*/
  area=0;
  switch (selmode){
  case 0:case 5:
    for (i=1;i<=nspikes;i++){
      if (upper>-save_mle.sval[i]*n_e&&lower<=-save_mle.sval[i]*n_e){
	area+=save_mle.sprob[i];
      }
      //printf("%f %f %f %f %f %f %f %f\n",area,inc_x,inc_y,-lower,-upper, -save_mle.sval[i-1]*n_e,-save_mle.sval[i]*n_e,n_e);
      //monitorinput();
      //printf("%f \n",save_mle.sprob[i]);
      //monitorinput();
    }
    break;
  case 1:
    for (i=1;i<=nspikes;i++){
      inc_x=gsl_cdf_flat_P(lower, -save_mle.sval[i-1]*n_e, -save_mle.sval[i]*n_e);
      inc_y=gsl_cdf_flat_P(upper, -save_mle.sval[i-1]*n_e, -save_mle.sval[i]*n_e);
      area+=(inc_y-inc_x)*save_mle.sprob[i];
      //printf("%f %f %f %f %f %f %f %f\n",area,inc_x,inc_y,-lower,-upper, -save_mle.sval[i-1]*n_e,-save_mle.sval[i]*n_e,n_e);
      //monitorinput();

    }
    break;

  case 2:
    inc_x = gsl_cdf_gamma_P(lower/n_e, save_mle.beta,1/save_mle.alpha);
    inc_y = gsl_cdf_gamma_P(upper/n_e, save_mle.beta,1/save_mle.alpha);
    if (lower==0){area=inc_y;}else{area=(inc_y-inc_x);}
    area=(inc_y-inc_x);
    break;

  case 3:
    inc_x = gsl_cdf_beta_P (lower/n_e, save_mle.alpha, save_mle.beta);
    inc_y = gsl_cdf_beta_P (upper/n_e, save_mle.alpha, save_mle.beta);
    area=(inc_y-inc_x);
    break;

  case 4:
    inc_x = gsl_cdf_lognormal_P (lower/n_e, save_mle.alpha, save_mle.beta);
    inc_y = gsl_cdf_lognormal_P (upper/n_e, save_mle.alpha, save_mle.beta);
    area=(inc_y-inc_x);
    break;
  }

  return(area);
}
/******************************************************************************/
double compute_fix_prob(double *s_evaluated_vec,int n_s_evaluated)
{
  int i=0,j=0;
  int n1=0,n2=0,w1=0,w2=0;
  double t2=0,n_e=0;

  double * density_vec = (double*) calloc (maxnd+1, sizeof(double));
  double * step_density_vec = (double*) calloc (maxnd+1, sizeof(double));
 
  n1=N1;
  if (conpop){
    int n2=n1;
    t2=100;
    n_e=n1;
  }
  else {
    n_e=calculate_ne(n1,save_mle.N2,save_mle.t_real);
  }    

  double  fix_prob=0, lambda=0, s=0, 
    lower=0, upper=0, area=0, total_area=0, inc_x=0, inc_y=0,integral=0;

  total_area = 0;
  integral = 0;
  
  switch(selmode)
    {
    case 0:case 5:
      for (i=1; i<=nspikes; i++){
	set_up_density_vec_equal_effects(save_mle.sval[i],n_s_evaluated, step_density_vec);
	for (j=0; j<n_s_evaluated; j++){
	  density_vec[j]+=step_density_vec[j]*save_mle.sprob[i];
	} 

      }
      break;

    case 1:
      for (i=1; i<=nspikes; i++){
	set_up_density_vec_step_effects(save_mle.sval[i-1],save_mle.sval[i],n_s_evaluated, step_density_vec);
	for (j=0; j<n_s_evaluated; j++){
	  density_vec[j]+=step_density_vec[j]*save_mle.sprob[i];
	}
      }	
      break;	

    case 2: 
      compute_gamma_densities(save_mle.alpha,save_mle.beta, s_evaluated_vec, n_s_evaluated, density_vec);
      break;
	  	
    case 3: 
      compute_beta_densities(save_mle.alpha,save_mle.beta,save_mle.beta_scaler, s_evaluated_vec, n_s_evaluated, density_vec);
      break;
	  
    case 4: 
      compute_lognormal_densities(save_mle.alpha,save_mle.beta, s_evaluated_vec, n_s_evaluated, density_vec);
      break;
    } 
	
  for (i=0; i<n_s_evaluated; i++){

    s = s_evaluated_vec[i];
    if (s == 0){
      fix_prob = 0.5/(double)n_e;
    }
    else {
      fix_prob = kimura_fixation_prob(s, (double)n_e);
    }  
    integral += density_vec[i]*fix_prob;
    total_area+=density_vec[i];
    //printf("\ns %lf\n",s );
    //monitorinput();
  }//end for i
  //printf("\ntotal_area %lf\n",total_area );
  save_mle.fix_prob = 2*(double)n_e*integral;

  return(1);
}
/******************************************************************************/
double jukescantor(double sites, int diffs)
{
  double p=0, k=0;
  p = (double)diffs/(double)sites;
  k = -(3.0/4.0)*log(1.0 - 4.0*p/3.0);
  //printf("Jukescantor: p %lf k %lf\n", p, k);
  return k;
}

/******************************************************************************
 ******************************************************************************/
main(argc,argv)
int argc; char **argv;
{
  gsl_rng_env_setup(); 
  int i,j;
 
  char *sfs_filename;
    
  int fsim=0;
  int N2=0;
  int output=0;
  int gss_lower=2;
  int gss_upper=1000;
  nalleles=120;
  int ranrep;
  int acc=0;
  /*Read options*/
  
  while (1)
    {
      static struct option long_options[] =
	{
	  /* These options set a flag. */
	  // {"verbose", no_argument,       &verbose_flag, 1},
	  // {"brief",   no_argument,       &verbose_flag, 0},
	  /* These options don't set a flag.
	     We distinguish them by their indices. */
	  {"conpop",  required_argument,0, 'a'},
	  {"sfsfold",  required_argument, 0, 'b'},
	  {"selmode",  required_argument, 0, 'c'},
	  {"nspikes",    required_argument, 0, 'd'},
	  {"ranrep",    required_argument, 0, 'e'},
	  {"file",    required_argument, 0, 'f'},
	  {"fsim",    required_argument, 0, 'g'},
	  {"N2lower",    required_argument, 0, 'h'},
	  {"N2upper",    required_argument, 0, 'i'},
	  {"N2",    required_argument, 0, 'j'},
	  {"acc",    required_argument, 0, 'k'},
  
	};
      /* getopt_long stores the option index here. */
      int option_index = 0;
      int c=0;
      c = getopt_long_only (argc, argv, "",
			    long_options, &option_index);
     
      /* Detect the end of the options. */
      if (c == -1)
	break;
     
      switch (c)
	{
	case 0:
	  /* If this option set a flag, do nothing else now. */
	  if (long_options[option_index].flag != 0)
	    break;
	  printf ("option %s", long_options[option_index].name);
	  if (optarg)
	    printf (" with arg %s", optarg);
	  printf ("\n");
	  break;

	case 'a':
	  conpop=atoi(optarg);
	  break;     
	case 'b':
	  sfsfold=atoi(optarg);
	  break;     
	case 'c':
	  selmode=atoi(optarg);
	  break;     
	case 'd':
	  nspikes=atoi(optarg);
	  break;
	case 'e':
	  ranrep=atoi(optarg);
	  break;
	case 'f':
	  sfs_filename=optarg;
	  break;
	case 'g':
	  fsim=atoi(optarg);
	  break;
	case 'h':
	  gss_lower=atoi(optarg);	  
	  break;
	case 'i':
	  gss_upper=atoi(optarg);
	  break;
	case 'j':
	  N2=atoi(optarg);
	  break;
	case 'k':
	  acc=atoi(optarg);
	  break;
	case '?':
	  /* getopt_long already printed an error message. */
	  break;
     
	default:
	  abort ();
	}
    }

  /* Print any remaining command line arguments (not options). */
  if (optind < argc){
    printf ("non-option ARGV-elements: ");
    while (optind < argc)
      printf ("%s ", argv[optind++]);
    putchar ('\n');
  }

  /*set up random number generator
    the same random number generator is used for all operations to ensure that replicable results are produced.  
  */
  const gsl_rng_type * T1;
  T1 = gsl_rng_taus;
  gsl_rng  *rgen1 =gsl_rng_alloc(T1);

  /*parameterisation of selection*/
  npar=1;//f0 parameter
  
  switch (selmode)//add selection parameters
    {
    case 0:
      if (nspikes<=0)
	{
	  printf("\nFor models 0 (spike) and 1 (step) you have to have at least 1 spike! Exiting program...\n");
	  return(0);
	}
      npar+=nspikes;//add number of spikes
      npar+=nspikes-1;//add probability parameters-1
      break;
    case 1:  
      if (nspikes<=0){
	printf("\nFor models 0 (spike) and 1 (step) you have to have at least 1 spike! Exiting program...\n");
	return(0);
      }
      npar+=nspikes;//add number of spikes
      npar+=nspikes-1;//add probability parameters-1
      break;
    case 2:case 4:
      nspikes=0;
      npar+=2;// add alpha and beta parameters
      break;
    case 3:
      nspikes=0;
      npar+=3;//add shape 1,shape 2, beta_scaler parameters
      break;
    case 5:
      nspikes=6;
      npar+=nspikes-1;//add probability parameters-1
      break;
      
    default:
      printf("please select models 0-3\n");
      abort ();

    }

  /* Read the SFS*/ 
  float total_neu=0,total_sel=0;

  float *sfs_sel,*sfs_neu,*par;
  float *sel_sites,*sel_diff,*neu_sites,*neu_diff;

  int *nallelesp=&nalleles;
  get_nalleles(nallelesp,sfs_filename);
  
  sfs_sel = calloc (nalleles+1, sizeof(float));
  sfs_neu = calloc (nalleles+1, sizeof(float));

  /* read in the SFS
     FORMAT: 
     nalleles
     selected sfs[0..nalleles]
     neutral sfs[0..nalleles]
  */
  get_sfs(nalleles,sfs_sel,sfs_neu,sfs_filename);
  par = calloc (2*(nalleles+1), sizeof(float));

  /* get the directory path where the precomputed allele frequency vector tables are located*/
  get_data_path(data_path);

  /* Fold the SFS, if the user so chooses*/
  if (sfsfold){
    sfsfold_f(nalleles,sfs_sel,par,0);
    sfsfold_f(nalleles,sfs_neu,par,nalleles);
  }
  
  save_mle.preFV0 = (double*) calloc (maxnd+1, sizeof(double));
  save_mle.FV0 = (double*) calloc (maxnd+1, sizeof(double));

  save_mle.preFVS = (double*) calloc (maxnd+1, sizeof(double));
  save_mle.FVS = (double*) calloc (maxnd+1, sizeof(double));

  int sampleN,sampleS=0;
  for (i=0;i<nalleles;i++){
    sampleS=sampleS+(int)par[i];
    sampleN=sampleN+(int)par[i+nalleles];
    //printf("%f\t",par[i]);
    //printf("%f\n",par[i+nalleles]);
    // printf("%d\t",sampleS);
    //printf("%d\n",sampleN);
  }
 // monitorinput();

  /*load vectors to be evaluated*/

  set_up_file_name(N1, s_evaluated_vec_file_const, s_evaluated_vec_file);
  set_up_file_name(N1, s_range_file_const, s_range_file);
  get_s_evaluated_vec(s_evaluated_vec, &n_s_evaluated, &n_s_evaluated_file, 
		      s_evaluated_vec_file);
  get_s_ranges();

  if (!conpop){
    set_up_file_name(N1, n2_evaluated_vec_file_const, n2_evaluated_vec_file);
    set_up_file_name(N1, t2_evaluated_vec_file_const, t2_evaluated_vec_file);
    set_up_file_name(N1, phase_1_dir_const, phase_1_dir);
    set_up_file_name(N1, phase_2_dir_const, phase_2_dir);

    get_int_evaluated_vec(t2_evaluated_vec,&n_t2_evaluated, &t2_lower,
			  &t2_step,&t2_evaluated_vec_file);
    get_int_evaluated_vec(n2_evaluated_vec,&n_n2_evaluated, &n2_lower,
			  &n2_step,&n2_evaluated_vec_file);

    if (N2==0){
      /*launch N2 minimization*/

      N2=my_N2_minimization(par,gss_lower,gss_upper);
      my_minimization_neutral(N2,par);
      my_minimization(N2,par,ranrep,rgen1);
    }
    else{
      /*N2 is known*/
      save_mle.N2=N2;
      my_minimization_neutral(N2,par);
      my_minimization(N2,par,ranrep,rgen1);
    }
      
  }
 
  else{   /*constant population*/
    N2 = N1;      
    set_up_file_name(N1, "", const_pop_dir);
    my_minimization(N2,par,ranrep,rgen1);
      
  }
  save_mle.ne=calculate_ne(N1,save_mle.N2,save_mle.t_real);
  //monitorinput();
  compute_fix_prob(s_evaluated_vec, n_s_evaluated);

  ////////////////////////////////////////////////////////////////////////////////
  /*PRINTING results*/
  char filename_results[1000]="";
  strcat(filename_results,sfs_filename);
  FILE *fileresults= fopen(strcat(filename_results,".MAXL.out"), "a" );
 
  fprintf(fileresults,"seed:%lu\t",gsl_rng_default_seed);
  fprintf(fileresults,"acc:%lu\t",acc);
  fprintf(fileresults,"selmode:%d\t",selmode);
  fprintf(fileresults,"nspikes:%d\t",nspikes);
  fprintf(fileresults,"ranrep:%d\t",save_mle.actualreps);
  fprintf(fileresults,"L:%0.16f\t",save_mle.maxlik);
  fprintf(fileresults,"f0:%f\t",save_mle.freq_zero);


  if (!conpop)
    {
      fprintf(fileresults,"N2:%d\t",save_mle.N2); 
      fprintf(fileresults,"t:%f\t",save_mle.t_real);
      fprintf(fileresults,"Nw:%f\t",save_mle.ne);
    }
  switch (selmode)
    {
    case 0:case 1:case 5:
      for (i=1;i<=nspikes;i++){
	fprintf(fileresults,"s%d:%E\t",i,save_mle.sval[i]);
      }
      for (i=1;i<=nspikes;i++){
	fprintf(fileresults,"p%d:%E\t",i,save_mle.sprob[i]);
      }

      break;
    case 2:
      fprintf(fileresults,"E(s):%f\t",save_mle.es);
      fprintf(fileresults,"beta:%f\t",save_mle.beta);
      break;
    case 3:
      fprintf(fileresults,"alpha:%E\t",save_mle.alpha);
      fprintf(fileresults,"beta:%E\t",save_mle.beta);
      //fprintf(fileresults,"beta_scaler:%E\t",save_mle.beta_scaler);
      break;
 
    case 4:
      fprintf(fileresults,"logmean:%f\t",save_mle.alpha);
      fprintf(fileresults,"logsigma:%f\t",save_mle.beta);
      break;
    }//switch end

  /*calculate proportion of mutations in fixed Nes ranges*/
  double prop[4],Ne_s_lower=0,Ne_s_upper=0;
  Ne_s_lower=0;
  Ne_s_upper=0.1;

  for (i=0;i<=3;i++){
    prop[i]=compute_prop_muts(Ne_s_lower,Ne_s_upper);
    fprintf(fileresults,"Nes_%0.1f_%0.1f:%E \t",Ne_s_lower,Ne_s_upper,prop[i]);
    Ne_s_lower=Ne_s_upper;
    Ne_s_upper*=10;
  }
/*
  compute_mean_stats(s_evaluated_vec, n_s_evaluated);
  fprintf(fileresults,"mean:%E\tmean2:%E\tmeanH:%E\t",save_mle.mean_s,save_mle.mean_s2,save_mle.mean_h);
*/
  //printf("lethals:%f\t",save_mle.lethal_prop);
  //monitorinput();

  /*code for reading simulated means*/
  /*
  double MEAN_SIM[2], true_prop[4],true_fix_prob=0;
  if (fsim==1){
    //read simulated s and s^2 and harmonic mean
    char filename_mean[maxnd]="";
    strcat(filename_mean,sfs_filename);
    FILE *file1= fopen(strcat(filename_mean,".mean"), "r" );
    //printf("%s\n",filename_mean); 
    float read;
    for (i =0 ;i<=2; i++){
      fscanf (file1, "%f", &read);
      MEAN_SIM[i]=read;
    }
*/
/*
    fprintf(fileresults,"truemean:%E\ttruemean2:%E\ttruemeanH:%E\t",MEAN_SIM[0],MEAN_SIM[1],MEAN_SIM[2]);


    for(i = 0; i <= 3; i++){
      fscanf (file1, "%f", &read);
      true_prop[i]=read;
      fprintf(fileresults,"trueprop%d:%f\t",i,true_prop[i]);
    }
    
    fscanf (file1, "%f", &read);
    true_fix_prob=read;
    fprintf(fileresults,"true_fix_prob:%f\t",true_fix_prob);
    close(file1);
  }//if fsim=1 end
*/
  fprintf(fileresults,"fix_prob:%lf\t",save_mle.fix_prob);
  
  fprintf(fileresults,"\n");
  /*
  //save_mle.FVS[0]=save_mle.freq_zero;
  //dumpvector(save_mle.FV0, 0, N2*2, "FV0exp"); 
  //dumpvector(save_mle.FVS, 0, N2*2, "FVSexp"); 
  //printf("%d\n",sampleN);
  // printf("%d\n",sampleS);
  /
  */
  int * discreteN = (int*) calloc (nalleles+2, sizeof(int));
  int * discreteS = (int*) calloc (nalleles+2, sizeof(int));

  binomial_sampling(save_mle.N2,nalleles,sampleN,save_mle.FV0,discreteN);
  binomial_sampling(save_mle.N2,nalleles,sampleS,save_mle.FVS,discreteS);

  //////////////////////////// 
  char filename_SFS[maxnd]="";
  strcat(filename_SFS,sfs_filename);
  FILE *fileSFS= fopen(strcat(filename_SFS,".expSFS.out"), "a" );
  //fprintf(fileSFS,"%s\n0 0\n0 0\n%d\n",sfsname,nalleles);
  fprintf(fileSFS,"%d\n",nalleles);
  for (i=0;i<nalleles;i++)
    {
      fprintf(fileSFS,"%d ",discreteS[i]);
  
    }
  fprintf(fileSFS,"\n");
  for (i=0;i<nalleles;i++)
    {
      fprintf(fileSFS,"%d ",discreteN[i]);
  
    }
  
  fprintf(fileSFS,"\n\n");

  char filename_freq_SFS_sel[maxnd]="";
  strcat(filename_freq_SFS_sel,sfs_filename);
  FILE *freq_SFS_sel= fopen(strcat(filename_freq_SFS_sel,".obs-exp.neu-sel.SFS.out"), "w" );
  int *sfs_freqS,*sfs_freqN;

  sfs_freqS = calloc (nalleles, sizeof(int));
  sfs_freqN = calloc (nalleles, sizeof(int));
  
  sfsfold_f(nalleles,discreteS,sfs_freqS,0);
  sfsfold_f(nalleles,discreteN,sfs_freqN,0);
  
  fprintf(freq_SFS_sel,"no.derived_alleles\tneu_obs\tneu_exp\tsel_obs\tsel_exp\t\n");
  for (i=0;i<nalleles;i++)
    {
      fprintf(freq_SFS_sel,"%d\t%.0f\t%d\t%.0f\t%d\n",i,par[i+nalleles],sfs_freqN[i],par[i],sfs_freqS[i]);
  
    }
  fprintf(freq_SFS_sel,"\n");  


  //dumpvector(discreteN, 0, nalleles, "exp_neu"); 
  //dumpvector(discreteS, 0, nalleles, "exp_sel");

  free(discreteN);
  free(discreteS);
  
  free(save_mle.preFVS);
  free(save_mle.FVS);
  free(save_mle.sval);
  free(save_mle.sprob);
  free(sfs_sel);
  free(sfs_neu);
  free(par);
  gsl_rng_free (rgen1);
  return (0);
}


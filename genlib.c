/*GENLIB library routines */
/* Compile progs: gcc -O3 -o prog prog.c ~/genlib.c -lm -lgsl -lgslcblas */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define pi 3.14159265358979
#define infinity 999999
#define true 1
#define false 0
#define maxranges 5                    /* for numerical integration */
#define maxsimpsonpoints 1025          /* change libhdr if these changed */

#define MBIG 1000000000
#define MSEED 161803398
#define FAC (1.0/MBIG)


#define ib1 1
#define ib2 2
#define ib5 16
#define ib18 131072

#define undefined -99999999999999.9999999


typedef double (*fun1)();

FILE *seedfileptr, *tmoutfile, *fopen();


int inputcounter=0;

/* GSL random number definitions */
gsl_rng * gsl_rng_r;
unsigned long int seed;
int random_number_generator_initialised_flag = 0;

long sseed;  /* seed for the slow generator */
int tracelevel;


struct acc
{
   int n;
   double sum;
   double sumsq;
};

struct covacc
{
   int n;
   double sumxy, sumx, sumy, sumx2, sumy2;
};

gabort(char *s, int r)
{
   printf("ERROR %s %d\n", s, r);
   exit(1);
}

   
/* This is the random number generator ran3 from Numerical Recipes from
Knuth.  It's the fastest of the ``non-quick-and-dirty'' generators, and
is supposed to have decent properties */
double uniform_old()
{
	static int inext,inextp;
	static long ma[56];
	static int iff=0;
	long mj,mk;
        double res;
	int i,ii,k;
	
	if(sseed<0||iff==0)
	{
		iff=1;
		mj=MSEED-(sseed<0?-sseed:sseed);
		mj%=MBIG;
		ma[55]=mj;
		mk=1;
		for(i=1;i<=54;i++)
		{
			ii=(21*i)%55;
			ma[ii]=mk;
			mk=mj-mk;
			if(mk<0) mk+=MBIG;
			mj=ma[ii];
		}
		for(k=1;k<=4;k++)
			for(i=1;i<=55;i++)
			{
				ma[i] -=ma[1+(i+30)%55];
				if(ma[i]<0) ma[i]+=MBIG;
			}
		inext=0;
		inextp=31;
		sseed=1;
	}
	if(++inext==56) inext=1;
	if(++inextp==56) inextp=1;
	mj=ma[inext]-ma[inextp];
	if(mj<0) mj+=MBIG;
	ma[inext]=mj;
	sseed=mj;
        res = mj*FAC;
        if (res > 1.0) res = 1.0;
        else if (res < 0.0) res = 0.0;
	return res;
}



double uniform()
{
   double res;
   if (!random_number_generator_initialised_flag)
      gabort("uniform() called, but !random_number_generator_initialised_flag", 0);
   res = gsl_ran_flat(gsl_rng_r,0.0,1.0);
   return(res);
}


int irbit1(unsigned long *iseed)
/* routine to return random bits */
{
  unsigned long newbit;
  newbit = (*iseed &ib18) >> 17
  ^ (*iseed &ib5) >> 4
  ^ (*iseed &ib2) >> 1
  ^ (*iseed &ib1);
  *iseed = (*iseed<<1) | newbit;
  return (int)newbit;
}


initialise_uniform_generator()
{
/* create a generator chosen by the environment variable GSL_RNG_TYPE */
   const gsl_rng_type * T;
//   printf("initialise_uniform_generator\n");
   gsl_rng_env_setup();
   T = gsl_rng_default;
   gsl_rng_r = gsl_rng_alloc (T);
   random_number_generator_initialised_flag = 1;
}


getseed()
{
   FILE *seedfileptr;
   static int c = 1;
   if (c==1)         // first time routine called
   {
      initialise_uniform_generator();
      c = 0;
   }
   printf("Enter seed (-99 to read from file) ");
   manual: scanf("%lu", &seed);
   if (seed == -99)
   {
      seedfileptr = fopen("seedfile", "r");
      if (seedfileptr==0)
      {
         printf("No seedfile, enter seed please ");
         goto manual;
      }
      fscanf(seedfileptr, "%lu", &seed);
//      printf("Seed read %lu\n", seed);
      fclose(seedfileptr);
   }
   gsl_rng_set(gsl_rng_r, seed);
}

getseedquick()
{
   FILE *seedfileptr;
   static int c = 1;
   if (c==1)         // first time routine called
   {
      initialise_uniform_generator();
      c = 0;
   }
   seedfileptr = fopen("seedfile", "r");
   if (seedfileptr==0)
   {
      printf("No seedfile, enter seed please ");
      scanf("%lu", &seed);
   }
   else
   {
      fscanf(seedfileptr, "%lu", &seed);
      fclose(seedfileptr);
   }
//   printf("getseedquick: Seed read %lu\n", seed);
   gsl_rng_set(gsl_rng_r, seed);
}


terminate_uniform_generator()
{
   gsl_rng_free(gsl_rng_r);
}

writeseed()
{
   FILE *seedfileptr;
   double temp, x;
   x = uniform();
   temp = floor(x*100000000);
   seed = (unsigned long int)temp;
//   printf("Write seed %lu\n", seed);
   seedfileptr = fopen("seedfile", "w");
   fprintf(seedfileptr, "%lu\n", seed);
   fclose(seedfileptr);
   terminate_uniform_generator();
}




initacc(a)
struct acc *a;
{
   a -> n = 0;
   a -> sum = 0;
   a -> sumsq = 0;
}


accum(struct acc *a, double x)
{
   a->n = a->n + 1;
   a->sum = a->sum + x;
   a->sumsq = a->sumsq + x*x;
}

double accmean(struct acc *a)
{
   return(a->sum/(double)(a->n));
}

double variance(struct acc *a)
{
   double num, denom;
   if (a->n == 0) return((double)-infinity);
   num = a->sumsq - (a->sum)*(a->sum)/((double)(a->n));
   denom = (double)(a->n) - (double)1;
   return(num/denom);
}


/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! normal
! ------
! Returns 2 random long reals from the standard normal distribution.
!
! Parameters: mu = mean
!           sdev = standard deviation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/



double normal(double mu, double sdev)
{
  double u1, u2, r;
  u1= uniform();
  if (u1==0.0) u1 = 0.00001;            /*prevent fatal error*/
  if (u1==1.0) u1 = .999999;
  u2= uniform();
  r = sqrt (-2.0*log(u1)) * cos(2.0*pi*u2);
  return(r*sdev + mu);
}


/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! quick normal
! ------------
! Returns two uniform long reals from the normal distribution.
!
! Parameters: x1, x2 -> return values
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
quicknormal (double *x1, double *x2)
{
  double u1, u2, r;
  u1= uniform();
  if (u1==0.0) u1 = 0.00001;            /*prevent fatal error*/
  if (u1==1.0) u1 = .999999;
  u2 = uniform();
  *x1 = sqrt (-2.0*log(u1)) * cos(2.0*pi*u2);
  *x2 = sqrt (-2.0*log(u1)) * sin(2.0*pi*u2);
}



/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! wishart
! -------
!
! Returns two uniform real numbers from the wishart (bivariate) distribution.
! These are obtained by squaring correlated uniform numbers from the bivariate
! normal distribution.
!
! Parameters : z1 = return value for symmetrical distribution (i.e. metric)
!              z2 =  ,,     ,,    ,, negative sided distribution (i.e. fitness)
!              epsilon1 = sqrt(E(z1*z1))
!              epsilon2 = sqrt(E(z2*z2))
!              rho = correlation of absolute values.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
wishart(double *z1, double *z2, double epsilon1, double epsilon2, double rho)
{
   double sq;
   rho = sqrt(rho);
   sq = 1/sqrt(3.0);
   quicknormal(z1, z2);
   *z2 = *z1*rho + *z2*sqrt(1.0 - rho*rho);
   *z1 = sq*epsilon1*(*z1)*(*z1);
   *z2 = -sq*epsilon2*(*z2)*(*z2);
   if (uniform() > 0.5) *z1 = -(*z1);
}


/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! bivnormal
! -------
!
! Returns two real numbers from the bivariate normal distribution.
!
! Parameters : z1 = return value for symmetrical distribution (i.e. metric)
!              z2 =  ,,     ,,    ,, negative sided distribution (i.e. fitness)
!              epsilon1 = sqrt(E(z1*z1))
!              epsilon2 = sqrt(E(z2*z2))
!              rho = correlation of absolute values.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
bivnormal(double *z1, double *z2, double epsilon1, double epsilon2, double rho)
{
   rho = sqrt(rho);
   quicknormal(z1, z2);
   *z2 = *z1*rho + *z2*sqrt(1.0 - rho*rho);
   *z1 = epsilon1*(*z1);
   *z2 = -epsilon2*(*z2);
   if (*z2 > 0.0) *z2 = -*z2;
}



cubenormal(z1, z2, epsilon1, epsilon2, rho)
double *z1, *z2, epsilon1, epsilon2, rho;
{
   double sq;
   sq = 1.0/sqrt(16.0);
   rho = sqrt(rho);
   quicknormal(z1, z2);
   *z2 = *z1*rho + *z2*sqrt(1.0 - rho*rho);
   *z1 = sq*epsilon1*(*z1)*(*z1)*(*z1);
   *z2 = -sq*epsilon2*(*z2)*(*z2)*(*z2);
   if (*z2 > 0.0) *z2 = -*z2;
}



/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! doublegamma
! -----------
!
! Returns uniform long real from the double gamma distribution.
!
! Parameter: epsilon = sqrt(E(a*a))
!            proportion positive = proportion of distribution +ve
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
double doublegamma (epsilon, proppositive)
double epsilon, proppositive;
{
  double r, sigma;
  sigma = sqrt(epsilon*1.0/sqrt(3.0));
  r= normal(0.0, sigma);
  r= r*r;
  if (uniform() > proppositive) r= -r;       /*uniformly allocate sign*/
  return(r);
}


double gamdev(int ia, double epsilon)
/* returns random number from gamma distribution with integer parameter
taken from Numerical Recipes in C */
{
   int j;
   double x;
   if ((ia<1)||(ia>=6)) gabort("Parameter to gamdev() invalid ", (double)ia);
   x = 1.0;
   for (j=1; j<=ia; j++)
   {
      x *= uniform();
   }
   x = -log(x);
   x = (epsilon*x)/sqrt((double)(ia*(ia+1)));
   return x;
}



initcovacc(struct covacc *a)
{
   a->n = 0;
   a->sumx = 0;
   a->sumy = 0;
   a->sumxy = 0;
   a->sumx2 = 0;
   a->sumy2 = 0;
}


covaccum(struct covacc *a, double x, double y)

{
   a->n = a->n + 1;
   a->sumx = a->sumx + x;
   a->sumy = a->sumy + y;
   a->sumxy = a->sumxy + x*y;
   a->sumx2 = a->sumx2 + x*x;
   a->sumy2 = a->sumy2 + y*y;
}



double covariance(struct covacc *a)
{
   if (a->n < 2)
   {
      return(-infinity);
   }
   return((a->sumxy - (a->sumx*a->sumy/(double)a->n))/((double)a->n - 1.0));
}



double correlation(struct covacc *a)
{
   double num, denom, cov, v1, v2;
   if (a->n < 2)
   {
      return(-infinity);
   }
   cov = covariance(a);
   num = a->sumx2 - (a->sumx)*(a->sumx)/((double)(a->n));
   denom = (double)(a->n) - 1.0;
   v1 = num/denom;
   num = a->sumy2 - (a->sumy)*(a->sumy)/((double)(a->n));
   denom = (double)(a->n) - 1.0;
   v2 = num/denom;
   return(cov/sqrt(v1*v2));
}


double se(a)
struct acc *a;
{
   if (a->n < 1) return(-infinity);
   if (variance(a)/a->n < 0) return(-infinity);
   return(sqrt(variance(a)/a->n));
}


printmse(s, r)
char *s;
struct acc *r;
{

   printf("%s %f +/- %f\n", s, accmean(r), se(r));
}


/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! generate poisson table
! ----------------------
!
! generates a table of probabilities for each whole number given a poisson
! distribution.
!
! Parametes: mean = mean of poisson.
!            last number = max. elements in table
!            table -> table of probabilities.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

generatepoissontable (double mean, int *lastnumber, double *table, int max)
{
   int j;
   double prev;
   prev = exp(-mean);
   if (prev == 0) gabort("Insufficient resolution in generate poisson table", 0);
   table[0] = prev;
   for (j = 1; j <= max; j++)
   {
      *lastnumber = j;
      prev = prev*mean/(double)j;
      table[j] = table[j-1] + prev;
      if (1 - table[j] < 0.000001) break;
   }
}




/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! poisson
! -------
!
! Looks up a table of probabilities of the poisson distribution already
! set up.  Returns the index of this probability.  This discrete vaue
! will come from the posson distribution.
!
! Parameters: last number = last entry in the table
!             table      -> table of probabilities
!
! Returns   : integer from poisson distribution.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

int poisson (int lastnumber, double *table)
{
   int j;
   float u;
   u = uniform();
   for (j = 0; j<=lastnumber; j++)
   {
      if ((double)u < table[j]) return(j);
   }
   return(lastnumber);
}


/*-------------------------------------------------------------------*/
/* binarysearchcumulative                                            */
/* ----------------------                                            */

/* Assumes array contains a cumulative density function of a discrete*/
/* distribution.  Returns an                                         */
/* integer randomly drawn from the density function.                 */

/* Parameters array -> cdf.                                          */
/*            size = no. of elements in array [1..size]              */
/*-------------------------------------------------------------------*/

int binarysearchcumulative(array, size)
double *array;
int size;
{
   int cur, pre, post;
   double r;
   r = uniform();
/*   printf("Uniform %f\n", r);*/
   pre = 1;
   post = size;
   cur = (size + 1)/2;
   do
   {
      if (array[cur] < r) pre = cur; else post = cur;
      cur = (pre + post)/2;
      if ((cur==post) || (cur==pre))
      {
         if (r < array[pre])
         {
/*            printf("Returning pre %d\n", pre);*/
            return(pre);
         }
         else
         {
/*            printf("Returning post %d\n", post);*/
            return(post);
         }
      }
   } while (size > 0);
}



/*----------------------------------------------------------------------*/
/*									*/
/* discrete 								*/
/* --------								*/
/*									*/
/* Returns random integer in range 1..n with equal probability.          */
/* Parameter: n = max. integer to return.				*/
/*									*/
/*----------------------------------------------------------------------*/
int discrete(int n)
{
   int res;
   res = (int)(uniform()*(double)n) + 1;
   if (res<1) res = 1;
   else if (res>n) res = n;
   return(res);
}



   
/*----------------------------------------------------------------------*/
/*                            						*/
/* samplewithoutreplacement						*/
/* ------------------------						*/
/*									*/
/* Returns a random integer from those present in array, which is valid */
/* from 1 to limit.  Overwrites sampled integer with array[limit],      */
/* then decrements limit.						*/
/* Parameters: limit -> no. valid integers in array                     */
/*             array -> array[1..limit] of integers to sample           */
/* 									*/
/*----------------------------------------------------------------------*/
int samplewithoutreplacement(limit, array)
int *limit, *array;
{
   int index, res;
   index = discrete(*limit);
   res = array[index];
   array[index] = array[*limit];
   *limit = *limit - 1;
   return(res);
}
   
/* trap - asks for input to stop program */
trap()
{
   int i;
   printf("Enter any int to continue\n");
   scanf("%d", &i);
}


/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! getint;
! -------
!
! Displays a prompt and reads in an integer.  Checks the range of the integer
! and stops the program if it is out of range.
!
! Parameters: s -> string prompt
!             i -> integer variable to input.
!           min =  max value
!           max = min value.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
getint(char *s, int *i, int min, int max)
{
   printf("%s", s);
   scanf("%d", i);
   if ((min != -infinity) && (*i < min))
   {
      printf("Value too small. Reading %s. Terminating.", s);
      printf("\n");
      printf("Min allowable: ");
      printf("%d", min);
      printf("\n");
      gabort("", 0);
   } 
   if ((max != infinity) && (*i > max)) 
   {
      printf("Value too large. Reading %s. Terminating.", s);
      printf("\n");
      printf("Max allowable: ");
      printf("%d", max);
      printf("\n");
      gabort("", 0);
   }
}

/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! getint and skip;
! ----------------
!
! Displays a prompt and reads in an integer, skipping
! to the end of the line.  Checks the range of the integer
! and stops the program if it is out of range.
!
! Parameters: s -> string prompt
!             i -> integer variable to input.
!           min =  max value
!           max = min value.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
getintandskip(char *s, int *i, int min, int max)
{
   char ch;
   printf(s);
   printf(" ");
   scanf("%d", i);
   do
   {
     ch = getchar();
   }
   while (ch != '\n');
   if ((min != -infinity) && (*i < min))
   {
      printf("Value too small. Reading %s. Terminating.", s);
      printf("\n");
      printf("Min allowable: ");
      printf("%d", min);
      printf("\n");
      gabort("", 0);
   } 
   if ((max != infinity) && (*i > max)) 
   {
      printf("Value too large. Reading %s. Terminating.", s);
      printf("\n");
      printf("Max allowable: ");
      printf("%d", max);
      printf("\n");
      gabort("", 0);
   }
}


/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
printreal(s,r,  x, y)
char *s;
double r;
int x, y;
{
   printf("%s %f\n", s, r);

}



/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! printint;
! ---------
!
! Prints out the integer with leading spaces x after the string s.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
printint (s, i, x)
char *s;
int i, x;
{
   printf("%s %d", s, i);
   printf("\n");
}


/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! print tail
! ----------
!
! Prints a line of dashes.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/

printtail()
{
   int j;
   printf("\n");
   for (j=1; j<=80; j++)
   {
      printf("-");
   }
   printf("\n\n");
}




/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! cointoss
! --------
!
! Uses random number to return false with probability 0.5 otherwise
! returns a non zero integer.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
int cointoss()
{
   if (uniform() > 0.5) return(0);
   return(1);
}



/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! printline
! ----------
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
printline (s)
char *s;
{
   printf(s);
   printf("\n");
}

spaces(n)
int n;
{
   int i;
   for (i=1; i<=n; i++) printf(" ");
}

/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! getreal and skip;
! ----------------
!
! Displays a prompt and reads in an real, skipping
! to the end of the line.  Checks the range of the integer
! and stops the program if it is out of range.
!
! Parameters: s -> string prompt
!             i -> integer variable to input.
!           min =  max value
!           max = min value.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
getrealandskip(char *s, double *i, double min, double max)
{
   char ch;
   printf(s);
   printf(" ");
   scanf("%lf", i);
   do
   {
     ch = getchar();
   }
   while (ch != '\n');
   if ((min != -infinity) && (*i < min))
   {
      printf("Value too small. Reading %s. Terminating.", s);
      printf("\n");
      printf("Min allowable: ");
      printf("%f", min);
      printf("\n");
      gabort("", 0);
   } 
   if ((max != infinity) && (*i > max)) 
   {
      printf("Value too large. Reading %s. Terminating.", s);
      printf("\n");
      printf("Max allowable: ");
      printf("%f", max);
      printf("\n");
      gabort("", 0);
   }
}

double calculaterealmean(double *a, int n)
{
   double sum;
   int i;
   sum = 0.0;
   for (i=1; i<=n; i++)
   {
      sum = sum + a[i];
   }
   return(sum/(double)n);
}
   
tracestart()
{
printf("Enter trace level ");
scanf("%d", &tracelevel);
}

trace(s, i)
{
   if (tracelevel==0) return(0);
   printf("%s %d\n", s, i);
}


outputrep(int j, int howoften)
{
  if ((double)j/(double)howoften  - floor((double)j/(double)howoften) == 0.0)
     printf("Completed iteration %d\n", j);
}


outputrepdot(int j, int reps)
{
  int howoften;
  howoften = reps/10;
  if ((double)j/(double)howoften  - floor((double)j/(double)howoften) == 0.0)
  {
     printf("%d.", j/howoften);
     if (j==reps) printf("\n");
     fflush(stdout);
  }
}



monitorinput()
{
   inputcounter--;
   if (inputcounter <= 0)
   {
      printf("Enter int to continue ");
      scanf("%d", &inputcounter);
   }
}


double normalheight(double x, double mu, double var)
{
   double res;
   res  = (1.0/sqrt(2.0*pi*var))*exp(-(x-mu)*(x-mu)/(2.0*var));
   return(res);
}

/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! simpson
! -------
!
! Computes the area under the function using simpsons rule for the set of
! points given.
!
! Parameters: x -> array of pairs of points.
!              points = no. of points.
!                   a = lower value
!                   v = upper value
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/

double simpson (double *x, int points, double a, double b)
{
   int j, coeff;
   double  res, h;
   if (2*(points/2) == points)
   {
      printf("FATAL ERROR: Even no. of points for SIMPSON.\n");
      monitorinput();
   }
   h = (b-a) / ((double)points-1.0);
   res = 0.0;
   for (j = 1; j<=points; j++)
   {
      if ((j == points) || (j==1 ))
      {
         coeff = 1;
      }
      else if (j == 2)
      {
         coeff = 4;
      }
      else if (coeff == 2) coeff = 4; else coeff = 2;
/*      %if trace level > 0 %then printstring("Coeff, x(j) ") %and write(coeff, 3) %and print(x(j), 2, 4) %and newline*/
      res = res + (double)coeff*x[j];
   }
   return((h/3.0) * res);
}




int odd(int i)
{
   int x;
   x = i/2;
   if ((double)(i)/2.0 - (double)x == 0.0) return(false); else return(true);
}




/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! solve quadratic
! ---------------
!
! Returns the roots of quadratic of parameters given.
!
! Parameters: a, b, c = parameters of quadratic.
!             root1, root2 = roots of quadratic.
!
! Returns : true => real roots.
!          false => no real root.
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
int solvequadratic(double a, double b, double c, double *r1, double *r2)
{
   double x;
   x = b*b - 4.0*a*c;
   if (x < 0.0) return(false);
   *r1 = (-b + sqrt(x))/(2.0*a);
   *r2 = (-b - sqrt(x))/(2.0*a);
   return(true);
}

   



int skiptoendofline(FILE *fptr)
{
   int status;
   char ch;
   do
   {
      status = fscanf(fptr, "%c", &ch);
      if (status==EOF) break;
/*      printf("skiptoendofline ch %c\n", ch);*/
   }
   while (ch!='\n');
   return(status);
}



/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!
! aitken
! ------
!
! Uses Aitken's formula to predict the asymptote from the three last points
! in the array given.
!
! Parameters: x-> array of points.
!               t = max. point in array.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
double aitken (double *a, int t)
{
   double x;
   x = a[t-2] - 2.0*a[t-1] + a[t];
   if (x == 0) return(a[t]);
   return(-((a[t-1] - a[t])*(a[t-1] - a[t]))/x + a[t]);
}

int factorial(int i)
{
   int res;
   res = 1;   
   for (;;)
   {
      if (i==0) break;
      res = res*i;
      i--;
   }
   return(res);
}




double doublefactorial(double x)
{
   double res;
   x = floor(x);
   res = 1.0;   
   for (;;)
   {
      if (x==0.0) break;
      res = res*x;
      x = x - 1.0;
   }
   return(res);
}



double logdoublefactorial(double x)
{
   double res;
   x = floor(x);
   res = 0.0;
   for (;;)
   {
      if (x==0.0) break;
      res = res + log(x);
      x = x - 1.0;
   }
   return(res);
}




/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! gammaalpha
! ----------
!
! Returns the value of the parameter alpha of the gamma function.
!
! 
! Parameters: beta = parameter of gamma func.

!            epsilon = sqrt(E(a*a))
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
double gammaalpha( beta, epsilon)
double beta, epsilon;
{
   double res;
   res = sqrt(beta*(beta + 1.0))/epsilon;
   return(res);
}




double asymptote_gamma(double z)
{
   double res;
   res = (z - 0.5)*log(z) -z +0.5*log(2.0*pi) + 1.0/(12.0*z) 
        - 1.0/(360.0*z*z*z) + 1.0/(1260.0*z*z*z*z*z) - 1.0/(1680.0*z*z*z*z*z*z*z);
   return(exp(res));
}
   
    
   
double series_gamma(double z)
/* return gamma(z) using series expansion 6.1.34 of Abramowitz and Stegun */
{
   #define maxc 26
   double c[maxc+1], res;
   int i;
   if (z > 2.0) return(asymptote_gamma(z));
   c[1] = 1.0;
   c[2] = 0.5772156649015329;
   c[3] =-0.6558780715202538;
   c[4] =-0.0420026350340952;
   c[5] = 0.1665386113822915;
   c[6] =-0.0421977345555443;
   c[7] =-0.0096219715278770;
   c[8] = 0.0072189432466630;
   c[9] =-0.0011651675918591;
   c[10]=-0.0002152416741149;
   c[11]= 0.0001280502823882;
   c[12]=-0.0000201348547807;
   c[13]=-0.0000012504934821;
   c[14]= 0.0000011330272320;
   c[15]=-0.0000002056338417;
   c[16]= 0.0000000061160950;
   c[17]= 0.0000000050020075;
   c[18]=-0.0000000011812746;
   c[19]= 0.0000000001043427;
   c[20]= 0.0000000000077823;
   c[21]=-0.0000000000036968;
   c[22]= 0.0000000000005100;
   c[23]=-0.0000000000000206;
   c[24]=-0.0000000000000054;
   c[25]= 0.0000000000000014;
   c[26]= 0.0000000000000001;
   res = 0.0;
   for (i=1; i<=26; i++)
   {
      res += c[i]*pow(z, (double)i);
/*      printf("gamma(z) %15.12f\n", 1.0/res); */
   }
   return(1.0/res);
}





/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! gammaht
! -------
!
! Returns the height at point a of a gamma distribution.
!
! Parameters: epsilon, beta = gamma parameters.
!             a = value to evaluate.
!
! Returns gamma fn.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
double gammaht(double epsilon, double beta, double a)
{
   double gammabeta, res, a1, a2, a3, alpha;
   if (a==0.0)
   {
      printf("Illegal value (0.0) for gammaht()\n");
      monitorinput();
   }
   if (beta < 0.0) gabort("gammaht: invalid beta\n", beta);
   alpha = gammaalpha(beta, epsilon);
   gammabeta = series_gamma(beta);
   a1 = pow(alpha, beta);
   a2 = exp(-alpha*a);
   a3 = pow(a, beta - 1.0);
   res = a1*a2*a3/gammabeta;
   return(res);
}


/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! half points
! -----------
!
! This routine is called to check that the simpson numerical integration
! is converging satisfactorily.  It halves the number of points in
! the array specified.
!

! Parameters: a -> array of points.
!              points -> no. of points.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
halfpoints(double *a, int *points)
{
   int ptr1, ptr2, newpoints;
   if (odd(*points)==false)
   {
      printf("Half points: ERROR: no. of points for simpson must be odd.35\n");
      monitorinput();
   }
   ptr1 = 1;
   ptr2 = 1;
   do
   {
      a[ptr1] = a[ptr2];
      newpoints = ptr1;
      ptr1++;
      ptr2 = ptr2 + 2;
   }
   while (ptr2 <= *points);
   *points = newpoints;
}



setuppoints(int pts, double *points, double lower, double upper, fun1 fn)
{
      int i;
      double interval, x, a;
      interval = (upper - lower)/((double)pts - 1.0);
      x = lower;
      for (i = 1; i<=pts; i++)
      {
         a = fn(x);
         points[i] = a;
         x = x + interval;
      }

}


/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!
! numerical integ
! ---------------
!
! Numerically integrates the function fn and returns the result in res.
! Integrates over a number of ranges with different numbers of points per
! range.  Points must be a number like 65, 129, 257, 513... in each
! range.  Halves number of points to do integration 3 times as a check on
! convergence
!
! Parameters:
!
!
! resvec: contains vector of result for different numbers of points
! start, finish: values limiting the function
! fn(x): is the function to integrate
! ranges: no. of ranges of points to evaluate.
! points per range: no of points in each range.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/

double numericalinteg(double *resvec, double start, double finish, fun1 fn, int ranges,
   int *pointsperrange)
{
   double lower, upper, r;
   double res[maxranges+2][4];
/*        The first dimension is contains the integral for each range, the
       last element containing the total.  The second dimension contains
       results for each halving of the points.
*/
   double points[maxsimpsonpoints+1];
   int i, j, k, fac, p;
   static int tmoutfileflag = 0;           
   if (tmoutfileflag==0)
   {
/*      tmoutfile = fopen("numericalinteg.out", "w");*/
      tmoutfileflag = 1;
   }
/*   printf("Start %f, finish %f, ranges %d\n", start, finish, ranges);*/
   upper = start;
   for (i = 1;  i<= ranges; i++)
   {
      lower = upper;
      upper = lower + (finish - start)/(double)ranges;
      p = pointsperrange[i];
      if (p > maxsimpsonpoints) gabort("too many points for simpson",p);
      setuppoints(p, points, lower, upper, fn);
      for (j = 1; j<=3; j++)
      {
         r = simpson(points, p, lower, upper);
         res[i][j] = r;
         if (j != 3) halfpoints(points, &p);
      }
   }
   for (i = 1; i<=3; i++)
   {
      res[ranges+1][i] = 0.0;
      for (j = 1; j<=ranges; j++)
         res[ranges+1][i] = res[ranges+1][i] + res[j][i];
      resvec[i] = res[ranges+1][i];
   }
/*   fprintf(tmoutfile, "Results of numerical integration\n");
   for (i=1; i<=ranges; i++) fprintf(tmoutfile, " Pts   Range %d", i);
   fprintf(tmoutfile, "     All");
   fprintf(tmoutfile, "\n");
*/
   fac = 1;
/*   for (i = 1; i<=3; i++)
   {
      for (j = 1; j<=ranges+1; j++)
      {
         if (j!=ranges+1) fprintf(tmoutfile, "%4d", (pointsperrange[j]-1)/fac + 1);
         fprintf(tmoutfile, "%10.4f", res[j][i]);
      }
      fac = fac*2;
      fprintf(tmoutfile, "\n");
   }
*/
   return(res[ranges+1][1]);
}



FILE *openforread(char *str)
{
   FILE *f;
   f = fopen(str, "r");
   if (f==0)
   {
      printf("ERROR: File %s not found.\n", str);
      exit(0);
   }
   else //printf("Opened file %s for read.\n", str);
   return(f);
}

FILE *openforwrite(char *str, char *mode)
{
   FILE *f;
   f = fopen(str, mode);
   if (f==0)
   {
      printf("ERROR: Unable to open file %s for write.\n", str);
      exit(0);
   }
   else //printf("Opened file %s for write, ", str);
   if (mode[0]=='a') printf("append mode.\n");
   else if (mode[0]=='w') printf("overwrite mode.\n");
   else gabort("Invalid mode given in openforwrite\n", 0);
   return(f);
}


FILE *openforreadsilent(char *str)
{
   FILE *f;
   f = fopen(str, "r");
   if (f==0)
   {
      printf("ERROR: File %s not found.\n", str);
      exit(0);
   }
   return(f);
}

FILE *openforwritesilent(char *str, char *mode)
{
   FILE *f;
   f = fopen(str, mode);
   if (f==0)
   {
      printf("ERROR: Unable to open file %s for write.\n", str);
      exit(0);
   }
   if (mode[0]=='a')
   {
   }
   else if (mode[0]=='w')
   {
   }
   else gabort("Invalid mode given in openforwrite\n", 0);
   return(f);
}


 
printmsefile(f, s, r)
FILE *f;
char *s;
struct acc *r;
{

   fprintf(f, "%s %f +/- %f\n", s, accmean(r), se(r));
}




quadratic_regression(double n, double sumx, double sumx2, double sumy, double sumx3,
  double sumxy, double sumx4, double sumx2y, double *b1, double *b2, double *b3)
{
/* Solve system of 3 simultaneous equations:

a*b1 + b*b2 + c*b3 = d         (1)
e*b1 + f*b2 + g*b3 = h         (2)
i*b1 + j*b2 + k*b3 = l         (3)

The actual parameters are appropriate for solving a quadratic regression -

y = b3*x^2 + b2*x + b3

See assignments below.

*/
   double a, b, c, d, e, f, g, h, i, j, k, l;
   double z1, z2, z3, z4, z5, z6;
   a = n;
   b = sumx;
   c = sumx2;
   d = sumy;
   e = b;
   f = c;
   g = sumx3;
   h = sumxy;
   i = c;
   j = g;
   k = sumx4;
   l = sumx2y;
/*   printf("\n%lf %lf %lf %lf\n", a, b, c, d);
   printf("%lf %lf %lf %lf\n", e, f, g, h);
   printf("%lf %lf %lf %lf\n\n", i, j, k, l);
*/

/* (1) - (2) -> b1*z1 + b2*z2 = z3  (4)     */
   z1 = a/c - e/g;
   z2 = b/c - f/g;
   z3 = d/c - h/g;

/* (1) - (3) -> b1*z4 + b2*z5 = z6  (5)     */
   z4 = a/c - i/k;
   z5 = b/c - j/k;
   z6 = d/c - l/k;


/* (4) - (5) */

   *b1 = (z3/z2 - z6/z5)/(z1/z2 - z4/z5);

   *b3 = ((d - a*(*b1))/b - (h - e**b1)/f) / (c/b - g/f);
  
   *b2 = (d - c*(*b3) - a*(*b1))/b;
}



int findtext(char *s, FILE *inptr)
{
   int len, i, curchar, dum;
   char c;
   curchar = 0;
   len = 0;
   for (i=0; i<=100; i++)
   {
      if (s[i] == '\0') break;
      len++;
   }
   for (;;)
   {
      c = getc(inptr);
      if (c == EOF) return(false);
      if (c==s[curchar])
      {
         curchar++;
         if (curchar==len) return(true);
      }
      else curchar = 0;
   }
}

/* FIND1DMAX */
/* Finds the maximum of 3 points using a quadratic approximation */

find1dmax(double *xvec, double *yvec, double *xmax, double *ymax)
{
   int i;
   double a, b, c, x1, x2, x3, y1, y2, y3, temp;
   double num, denom;
/*   for (i=1; i<=3; i++)
   {
      printf("xvec[%1d] %lf yvec[%1d] %lf\n", i, xvec[i], i, yvec[i]);
   }
*/
   x1 = xvec[1];
   x2 = xvec[2];
   x3 = xvec[3];
   y1 = yvec[1];
   y2 = yvec[2];
   y3 = yvec[3];
   if (y1==y3) y1+= 0.000000000000001;
   if (y3==y2) y2+= 0.0000000000000001;
/*   printf("x1 %f x2 %f x3 %f y1 %f y2 %f y3 %f\n", x1, x2, x3, y1, y2, y3);*/
   num = (y1 - y3)*(x1 - x2)/(x1 - x3) - y1 + y2;
   denom = (x1*x1 - x3*x3)*(x1 - x2)/(x1 - x3) - (x1*x1 - x2*x2);
   a = num/denom;
   num = y1 - y2 - a*(x1*x1 - x2*x2);
   b = num/(x1 - x2);
   c = y1 - a*x1*x1 - b*x1;
   temp = -b/(2.0*a);
   *xmax = temp;
   *ymax = a*temp*temp + b*temp + c;
}



/* FIND2DMAX */
/* Finds the maximum of a 3x3 grid using a sequential quadratic approximation */

double find2dmax(double *xvec, double *yvec, double zmatrix[4][4], double *xmax, double *ymax)
{
   int i, j;
   double vec[4], vec1[4], temp, max, max1, max2;
/*   for (i=1; i<=3; i++)
   {
      for (j=1; j<=3; j++)
      {
         if (zmatrix[i][j]==undefined) return(undefined);
      }
   }
*/
/* find the maximum for columns */
   for (i=1; i<=3; i++)
   {
      for (j=1; j<=3; j++)
      {
         vec[j] = zmatrix[i][j];
      }
      find1dmax(yvec, vec, &temp, &max);
/*      printf("Find1smax: zmax %lf ymax %lf\n", max, temp);*/
      vec1[i] = max;
   }
   find1dmax(xvec, vec1, xmax, &max1);
/*   printf("Columns: max %lf xmax %lf\n", max1, *xmax);*/
/* find the maximum for rows */
   for (i=1; i<=3; i++)
   {
      for (j=1; j<=3; j++)
      {
         vec[j] = zmatrix[j][i];
      }
      find1dmax(xvec, vec, &temp, &max);
/*      printf("Find1dmax: zmax %lf ymax %lf\n", max, temp);*/
      vec1[i] = max;
   }
   find1dmax(yvec, vec1, ymax, &max2);
/*   printf("Rows: max %lf ymax %lf\n", max2, *ymax);*/
   max = (max1 + max2)/2.0;
   return(max);
}


int bnldev(double pp, int n)
{
   int j, bnl;
   double p;
   p = (pp < 0.5 ? pp : 1.0 - pp);
   bnl = 0;
   for (j=1; j<=n; j++)
      if (uniform() < p) ++bnl;
   if (p!=pp) bnl = n - bnl;
   return bnl;
}


void metropolis(double y[], int ndim, double(*func)(double []),
   int limit, double delta[], double prevres, double *temperature)
{
   #define step_change 0.2
   int i, j, param, accept;
   double step, res, diff, u;
   param = 0;
   *temperature = 1.0;
   for (i=1; i<=limit; i++)
   {
/*      printf("%d ", i);*/
      param++;
      if (param > ndim) param = 1;
      step = normal(0.0, delta[param]);
/*      printf("param %d delta[param] %lf, step %lf\n", param, delta[param], step);*/
      y[param] += step;
      res = (*func)(y);
      if (res > prevres)
      {
         diff = prevres - res;
         u = uniform();
/*         printf("diff %lf, exp(diff) %lf, u %lf\n", diff, exp(diff), u);*/
         if (u < pow(exp(diff), 1.0/(*temperature)))
         {
            accept = true;
/*            printf("Accepted\n");*/
         }
         else
         {
            accept = false;
/*            printf("Rejected\n");*/
         }
      }
      else
      {
         accept = true;
/*         printf("Accepted unconditionally\n");*/
      }
      if (accept == true)      /* Note assume to be minimizing not maximizing */
      {
         prevres = res;      /* accept the change */
         delta[param] *= (1.0 + step_change); /* Increase the step size */
      }
      else
      {
         y[param] -=step;     /* reject the change */
         delta[param] *= (1.0 - step_change); /* Decrease the step size */
      }
/*      printf("Dump of delta vector: ");
      for (j=1; j<=ndim; j++)
      {
         printf("%lf ", delta[j]);
      }
      printf("\n"); monitorinput();
*/
      *temperature -= 2.0/(double)limit;
      if (*temperature < 0.001) *temperature = 0.001;
   }
}


/* A slower way of generating poissons that the tabular method in poisson */
int genpoisson(double xm)
{
   static double sq, alxm, g, oldm = (-1.0);
   double em, t, y;
//   if (xm>300) gabort("genpoisson: xm too large", 0);
   if (xm>=200.0)   /* Use normal generator for very high xm */
   {
      return normal(0.0, sqrt(xm)) + xm;
   }
   if (xm!=oldm)
   {
      oldm = xm;
      g = exp(-xm);
   }
   em = -1;
   t = 1.0;
   do
   {
      ++em;
      t *= uniform();
   }
   while (t > g);
   return em;
}



void qsortreals (double *a, int from, int to)
   {
      int l, v;
      double d;
      if (from >= to) return;
      l = from;
      v = to;
      d = a[v];
      for (;;)
      {
         while ((l < v) && (a[l] <= d))
         {
            l = l + 1 ;
         }
         if (l==v) break;
         a[v] = a[l];
         while ((v > l) && (a[v] >= d))
         {
            v = v - 1;
         }
         if (v==l) break;
         a[l] = a[v];
      }

/*      ! now l=v*/
      a[v] = d;
      l = l-1;
      v = v+1;
      qsortreals(a, from, l);
      qsortreals(a, v, to);
}



#define maxofrs 100

FILE *openfilecheck2dir(char *str)
/*
 Attempts to open file str, if fails first time attempts to open
../str or ../../str
*/
{
   FILE *f;
   static char s[maxofrs];
   int i;
   f = fopen(str, "r");
   if (f==0)
   {
      s[0] = '.';
      s[1] = '.';
      s[2] = '/';
      i=0;
      for (;;)
      {
         if (i+3 > maxofrs) gabort("String too long", i);
         s[i+3] = str[i];
         if (str[i] == '\0') break;
         i++;
      }
      f = fopen(s, "r");
      if (f==0)
      {
         s[0] = '.';
         s[1] = '.';
         s[2] = '/';
         s[3] = '.';
         s[4] = '.';
         s[5] = '/';
         i=0;
         for (;;)
         {
            if (i+6 > maxofrs) gabort("String too long", i);
            s[i+6] = str[i];
            if (str[i] == '\0') break;
            i++;
         }
         f = fopen(s, "r");
         if (f==0)
         {
            printf("ERROR: File %s not found.\n", str);
            exit(0);
         }
         else  printf("Opened file %s for read.\n", s);
      }
      else  printf("Opened file %s for read.\n", s);
   }
   else printf("Opened file %s for read.\n", str);
   return(f);
}



float rtsafe(void (*funcd)(float, float *, float *), float x1, float x2, float xacc)
{
   #define MAXIT 100
   int j;
   float df, dx, dxold, f, fh, fl;
   float temp, xh, xl, rts;
   (*funcd)(x1, &fl, &df);
   (*funcd)(x2, &fh, &df);
   if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
      gabort("Rtsafe: roots must be bracketed", 0);
   if (fl == 0.0) return(x1);
   if (fh == 0.0) return(x2);
   if (fl < 0.0)
   {
      xl = x1;
      xh = x2;
   }
   else
   {
      xh = x1;
      xl = x2;
   }
   rts = 0.5*(x1+x2);
   dxold = fabs(x2-x1);
   dx = dxold;
   (*funcd)(rts, &f, &df);
   for (j=1; j<=MAXIT; j++)
   {
      if ((((rts-xh)*df-f)*((rts-xl)*df-f) >=0.0)
        || (fabs(2.0*f) > fabs(dxold*df)))
      {
         dxold = dx;
         dx = 0.5*(xh-xl);
         rts = xl +dx;
         if (xl==rts) return rts;
      }
      else
      {
         dxold = dx;
         dx = f/df;
         temp = rts;
         rts -= dx;
         if (temp==rts) return (rts);
      }
      if (fabs(dx) < xacc) return rts;
      (*funcd)(rts, &f, &df);
      if (f < 0.0)
         xl = rts;
      else
         xh = rts;
   }
   printf("warning: rtsafe: Maximum interations exceeded\n"); monitorinput();
   return (rts);
}


/* Algorith to generate bivariate normal deviates provided by Ian
   White
*/

bivariatenormal(double *z1, double *z2, double var1, double var2, double cov)
{
   double s;
   quicknormal(z1, z2);     /* uncorrelated normal deviates mean 0 SD 1 */
   *z1 *= sqrt(var1);       /* Scale Z1 according to the standard deviation */
   s = var2 - cov*cov/var1;
   if (s<0) gabort("Invalid parameters for bivnormal", s);
   *z2 = (cov/var1)*(*z1) + sqrt(s)*(*z2);
}


get_silent_mode(int argc, char *argv[], int *silent_mode)
{
   int n_command_line_param, i;
   *silent_mode = 0;
   n_command_line_param = argc;
//   printf("n_command_line_param %d\n", n_command_line_param);
   for (i=1; i<n_command_line_param; i++)
   {
//      printf("argv[%d] %s\n", i, argv[i]);
   }
   if (n_command_line_param==2)
   {
      if (strcmp(argv[1], "-s")==0)
      {
//         printf("Silent mode switched on\n");
         *silent_mode = 1;
      }
   }
}



make_numbered_file_name(char *outfilename, char *filestr, int ind)
{
   int len, i;
   int i1;
   char ch;
   len = strlen(filestr);
   for (i=0; i<len; i++) outfilename[i] = filestr[i];
   if (ind > 99) 
      gabort("makefilename: file index too large", ind);
   if (ind > 9)
   {
      i1 = ind/10;
      ch = i1 + '0';
      outfilename[len] = ch;
//      printf("Char 1 ind %d i1 %d ch %c\n", ind, i1, ch); monitorinput();
      i1 = ind - (ind/10)*10;
      ch = i1 + '0';
//      printf("Char 2 ind %d i1 %d ch %c\n", ind, i1, ch); monitorinput();
      outfilename[len+1] = ch;
      outfilename[len+2] = '\0';
   }
   else
   {
      ch = ind + '0';
//      printf("Char 1 ind %d i1 %d ch %c\n", ind, i1, ch); monitorinput();
      outfilename[len] = ch;
      outfilename[len+1] = '\0';
   }
//   printf("makefilename: results %s\n", outfilename);
}


// --------------------------------------------------------------------
// Routine for grid searching a function f that takes integer values as
// it's paramters
// --------------------------------------------------------------------

double fill_out_grid(double **grid, int x0, int x1, int y0, int y1,
   int step_size_x, int step_size_y,
   double (*f)(int, int), int *max_x, int *max_y, int *nfunc)
{
   double max_fun, f_result;
   int i, j, evaluated;
   *max_x = -1;
   *max_y = -1;
   max_fun = undefined;
   i = x0;
   for (;;)
   {
      j = y0;
      for (;;)
      {
         evaluated = 0;
         if (grid[i][j]==undefined)
         {
            f_result = (*f)(i, j);
            (*nfunc)++;
            grid[i][j] = f_result;
            evaluated = 1;
         }
         else f_result = grid[i][j];
//         printf("i %d j %d f_result %lf ", i, j, f_result);
//         if (evaluated == 1) printf("evaluated\n");
//         else printf("not evaluated\n");
         if ((f_result > max_fun)||(max_fun==undefined))
         {
            max_fun = f_result;
            *max_x = i;
            *max_y = j;
         }
         if (j==y1) break;
         j += step_size_y;
         if (j > y1) j = y1;
      }
      if (i==x1) break;
      i += step_size_x;
      if (i > x1) i = x1;
   }
   return max_fun;
}


int_to_string(int v, char *str, int max_len)
{
   int divisor = 10000000, i1, ind = 0, digit_encountered = 0;
   char ch;
   if (v > divisor) gabort("int_to_string: integer parameter to large", v);
   for (;;)
   {
      i1 = v/divisor;
//      printf("i1 %d divisor %d\n", i1, divisor);
      if (i1 >= 0)
      {
         ch = i1 + '0';
         if (ch!='0') digit_encountered = 1;
//         printf("ch %c\n", ch);
         if (digit_encountered)
         {
            str[ind] = ch;
            ind++;
            if (ind == max_len) gabort("int_to_string: string too long", ind);
         }
         v = v - i1*divisor;
      }
      divisor /= 10;
      if (divisor == 1) break;
//      monitorinput();
   }
   ch = v + '0';
//   printf("ch %c\n", ch);
   str[ind] = ch;
   ind++;
   if (ind == max_len) gabort("int_to_string: string too long", ind);
   str[ind] = '\0';
}


// bilinear_interpolation:
// First 4 parameters are positions on x and y axes
// Second 4 are function values
// x and y are positions to be evaluated.

double bilinear_interpolation(double x1, double x2, double y1, double y2,
   double x1y1, double x1y2, double x2y1, double x2y2,
   double x, double y)
{
   double res;
   res = (x2 - x)*(y2 - y)*x1y1 +
         (x - x1)*(y2 - y)*x2y1 +
         (x2 - x)*(y - y1)*x1y2 +
         (x - x1)*(y - y1)*x2y2;

   res /= ((x2 - x1)*(y2 - y1));
   return res;
}

// --------------------
// kimura_fixation_prob
// --------------------
// Returns fixation probability of an additive mutation, assuming that s
// is small in a diploid population
//
// Parameters
// ---------
// s - selection coefficient, difference between the homozygotes.
// N - Population effective size.
//

double kimura_fixation_prob(double s, int N)
{
   double num, denom;
   if ((s < 0.0)&&((double)N*s > -0.0000001))  // Trap numerical problems
   {
      num = 1.0;
      denom = 2.0*(double)N;
   }
   else if ((s > 0.0)&&((double)N*s < 0.0000001))
   {
      num = 1.0;
      denom = 2.0*(double)N;
   }
   else
   {
      num = 1 - exp(-s);
      denom = 1 - exp(-2*(double)N*s);
   }
//   printf("s %14.14lf N %d num %lf denom %lf num/denom %lf\n",
///   s, N, num, denom, num/denom);
//   monitorinput();
   return num/denom;
}

// --------------------
// kimura_fixation_prob
// --------------------
// Returns fixation probability of an additive mutation, assuming that s
// is small in a diploid population
//
// Parameters
// ---------
// s - selection coefficient, difference between the homozygotes.
// N - Population effective size.
//

double kimura_fixation_prob_X_loci(double s, int N)
{
   double num, denom;
   if ((s < 0.0)&&((double)N*s > -0.0000001))  // Trap numerical problems
   {
      num = 1.0;
      denom = 2.0*(double)N;
   }
   else if ((s > 0.0)&&((double)N*s < 0.0000001))
   {
      num = 1.0;
      denom = 2.0*(double)N;
   }
   else
   {
      num = 1 - exp(-s);
      denom = 1 - exp(-2*(double)N*s);
   }
//   printf("s %14.14lf N %d num %lf denom %lf num/denom %lf\n",
///   s, N, num, denom, num/denom);
//   monitorinput();
   return num/denom;
}


//----------------------------------------
// normal_cumulative_probability_function
// ---------------------------------------
//
// Returns the area under the normal pdf from -infinity to x.
// From algorithm described in Abramowitz and Stegun.

double normal_cumulative_probability_function(const double x)
{
  const double b1 =  0.319381530;
  const double b2 = -0.356563782;
  const double b3 =  1.781477937;
  const double b4 = -1.821255978;
  const double b5 =  1.330274429;
  const double p  =  0.2316419;
  const double c  =  0.39894228;

  if(x >= 0.0) {
      double t = 1.0 / ( 1.0 + p * x );
      return (1.0 - c * exp( -x * x / 2.0 ) * t *
      ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
  }
  else {
      double t = 1.0 / ( 1.0 - p * x );
      return ( c * exp( -x * x / 2.0 ) * t *
      ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
    }
}

//----------------------------------------
// dump_matrix
//----------------------------------------
//
// Prints out the matrix X with 1..rows, 1..cols rows and columns, respectively.
//
// s -> string to print
// X -> matrix, usually allocated using dmatrix(), dimension rows, cols
// rows = rows, 1 relative
// cols = cols, 1 relative

dump_matrix(char *s, double **X, int rows, int cols)
{
   int i, j;
   printf ("%s\n", s);
   for (i = 1; i <=rows; i++)
   {
      for (j = 1; j <= cols; j++)
      {
         printf("%lf ", X[i][j]);
      }
      printf("\n");
   }
}

//----------------------------------------
// compute_matrix_transpose
//----------------------------------------
//
// Generates the transpose of X.
//
// X -> matrix, usually allocated using dmatrix(), dimension rows, cols
// rows = rows of X, 1 relative
// cols = cols of X, 1 relative
// XT -> result matrix, allocated using dmatrix(), dimension cols, rows

compute_matrix_transpose(double **X, double **XT, int rows, int cols)
{
   int i, j;
   for (i=1; i<=rows; i++)
   {
      for (j=1; j<=cols; j++)
      {
         XT[j][i] = X[i][j];
      }
   }
}

//----------------------------------------
// multiply_matrices
//----------------------------------------
//
// Carry out matrix operation XY, putting result in RES.
//
// RES -> result matrix dimension rows_X, cols_X, usually allocated using
//          dmatrix()
// X -> matrix dimension rows_X, cols_X, usually allocated using dmatrix()
// Y -> matrix dimension rows_Y, cols_Y, usually allocated using dmatrix()
// rows_X = rows of X, 1 relative
// cols_X = cols of X, 1 relative
// rows_Y = rows of Y, 1 relative
// cols_Y = cols of Y, 1 relative

multiply_matrices(double **RES, double **X, double **Y, int rows_X, int cols_X,
   int rows_Y, int cols_Y)
// Matrix X has dimension rows_X, cols_X and Y has dimension rows_Y, cols_Y. The
// result has dimension rows_X, cols_Y.
{
   int i, j, k;
   double element;
   if (cols_X != rows_Y)
   {
      gabort("multiply_matrices: cols_X != rows_Y", 0);
   }
//   printf("multiply_matrices: rows_X %d cols_X %d\n", rows_X, cols_X);
//   dump_matrix("multiply_matrices: Matrix X", X, rows_X, cols_X);
//   monitorinput();
//   printf("multiply_matrices: rows_Y %d cols_Y %d\n", rows_Y, cols_Y);
//   dump_matrix("multiply_matrices: Matrix Y", Y, rows_Y, cols_Y);
//   monitorinput();
   for (i=1; i<=rows_X; i++)
   {
      for (j=1; j<=cols_Y; j++)
      {
         element = 0.0;
         for (k=1; k<=cols_X; k++)
         {
//            printf("X[%d][%d] %lf X[%d][%d] %lf\n",
//                 i, k, X[i][k], k, j, Y[k][j]);
            element += X[i][k]*Y[k][j];
//            printf("i %d j %d k %d element %lf\n", i, j, k, element);
//            monitorinput();
         }
//         printf("RES[%d][%d] %lf\n", i, j, element); monitorinput();
         RES[i][j] = element;
      }
   }
}

// ----------------------------------------------
// normal_approx_to_binomial
// ----------------------------------------------
//
// Uses the normal approximation to calculate the probability of i successes
// from n trials if the probabilituy of one sucess is p

double normal_approx_to_binomial(int i, int n, double p)
{
   double m, s, z_upper, z_lower, area_lower, area_upper;
   m = p*(double)n;
   s = sqrt((1.0-p)*p*(double)n);
//   printf("m %lf s %lf\n", m, s);
   z_lower = ((double)i - 0.5 - m)/s;
   area_lower = normal_cumulative_probability_function(z_lower);
   z_upper = ((double)i + 0.5 - m)/s;
   area_upper = normal_cumulative_probability_function(z_upper);
//   printf("z_lower %lf area_lower %lf z_upper %lf area_upper %lf\n",
//      z_lower, area_lower, z_upper, area_upper);
   return(area_upper - area_lower);
}

// ----------------------------------------------

// loadparam
// paramvec = vector with 2 elements: element 0 = 0 => parameter not variable
//                                    element 0 = 1 => parameter variable
//                                    element 1 location to get parameter value
// curpos = index - 1 from which to load parameter
// p: matrix into which parameter loaded
//

loadparam(double *paramvec, int *curpos, double **p)
{
   if (paramvec[0]!=0.0)
   {
      *curpos = *curpos + 1;
      p[1][*curpos] = paramvec[1];
   }
}

// ----------------------------------------------
// unloadparam
// paramvec = vector with 2 elements: element 0 = 0 => parameter not variable
//                                    element 0 = 1 => parameter variable
// curpos = index -1 in which to unload parameter
// x = vector containing parameter value

unloadparam(double *paramvec, int *curpos, double x[])
{
   if (paramvec[0]!=0.0)
   {
      *curpos = *curpos + 1;
      paramvec[1] = x[*curpos];
   }
}

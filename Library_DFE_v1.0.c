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


/******************************************************************************/
int get_average_egf_vecs_1_and_2(double s,int n1, int n2, double t2_real,
				 int t2_lower, int t2_upper, double *egf_vec, 
				 char *buffer_p1_n2_t2_lower,
				 char *buffer_p2_n2_t2_lower, 
				 char *buffer_p1_n2_t2_upper,
				 char *buffer_p2_n2_t2_upper)
{
  static double egf_vec1_lower[maxnd+1], egf_vec2_lower[maxnd+1], 
    egf_vec1_upper[maxnd+1], egf_vec2_upper[maxnd+1],
    egf_vec1[maxnd+1], egf_vec2[maxnd+1];
  double total_density;
  int i, res;
  int n1d=2*n1;
  int n2d=2*n2;

  if (s<-1) {
    /* 
       The contribution of strongly selected mutations
       is assumed to be inversely proportional to s, as
       expected from mutation selection balance, with a
       contribution from phase 2 mutations only.
    */

      for (i=0; i<=2*n2; i++)
	{
	  egf_vec1[i] = 0;
	  egf_vec2[i] = 0;
	}
      egf_vec2[1] = -2/s;
    }
   else {
      res =  get_binary_egf_vec(buffer_p1_n2_t2_lower, n2, t2_lower, 
				s, egf_vec1_lower);
      if (res==0) return 0;
      res =  get_binary_egf_vec(buffer_p2_n2_t2_lower, n2, t2_lower, 
				s, egf_vec2_lower);
      if (res==0) return 0;
      res =  get_binary_egf_vec(buffer_p1_n2_t2_upper, n2, t2_upper, 
				s, egf_vec1_upper);
      if (res==0) return 0;
      res =  get_binary_egf_vec(buffer_p2_n2_t2_upper, n2, t2_upper, 
				s, egf_vec2_upper);
      if (res==0) return 0;

      compute_weighted_average_egf_vec(t2_real, t2_lower, t2_upper,
				       egf_vec1_lower, egf_vec1_upper, egf_vec1, 2*n2);
      //   monitorinput();
      compute_weighted_average_egf_vec(t2_real, t2_lower, t2_upper,
				       egf_vec2_lower, egf_vec2_upper, egf_vec2, 2*n2);
    }

  total_density = compute_weighted_average(egf_vec, egf_vec1, egf_vec2, n1d, 2*n2);
  if (s==0)
    {
      total_density_s0 = total_density;
      //      printf("s=0: total_density %lf\n", total_density);
    }

  return (1);
}

/******************************************************************************/
void get_upper_lower_int(double parm, int *lower, int *upper, int n_int_evaluated, int *int_evaluated_vec)
{
  int i;
  *lower = undefined_int;
  *upper = undefined_int;
  //   printf("get_upper_lower_int parm %lf n_int_evaluated %d\n",
  //      parm, n_int_evaluated);
  for (i=0; i<n_int_evaluated; i++)
    {
      //      printf("int_evaluated_vec[%d] %d\n", i, int_evaluated_vec[i]);
      if ((double)int_evaluated_vec[i] >= parm)
	{
	  *upper = int_evaluated_vec[i];
	  if (i != 0) *lower = int_evaluated_vec[i-1];
	  return(1);
	}
    }
}
/******************************************************************************/
void get_upper_lower_double(double parm, double*lower, double *upper,
		       int n_int_evaluated, double *double_evaluated_vec)
{
  int i;
  *lower = undefined;
  *upper = undefined;
  //   printf("get_upper_lower_int parm %lf n_int_evaluated %d\n",
  //      parm, n_int_evaluated);
  for (i=0; i<n_int_evaluated; i++)
    {
      //      printf("int_evaluated_vec[%d] %d\n", i, int_evaluated_vec[i]);
      if (double_evaluated_vec[i] >= parm)
	{
	  *upper =double_evaluated_vec[i];
	  if (i != 0) *lower = double_evaluated_vec[i-1];
	  return(1);
	}
    }
}
/******************************************************************************/
void get_s_ranges()
{
  int i, stat, ind, i1;
  double d1, d2, d3;
  FILE *inf, *fopen();
  inf = openforreadsilent(s_range_file);
  stat = fscanf(inf, "%d", &nranges);
  if (stat!=1) gabort("get_s_ranges: input error 1", 0);
  if (nranges >max_s_ranges)
    gabort("get_s_ranges: nranges >max_s_ranges", nranges);
  for (i=1; i<=nranges; i++)
    {
      stat = fscanf(inf, "%d %lf %lf %lf %d", &ind, &d1, &d2, &d3, &i1);
      if (stat!=5) gabort("get_s_ranges: input error 2", 0);
      if (ind!=i) gabort("get_s_ranges: ind!=i", i);
      s_ranges[i].s_lower = d1;
      s_ranges[i].s_upper = d2;
      s_ranges[i].s_step = d3;
      s_ranges[i].s_lower_index = i1;
    }

}
/******************************************************************************/
void get_int_evaluated_vec(int *int_evaluated_vec, int *n_int_evaluated,
		      int *int_lower, int *int_step, 
		      char *int_evaluated_vec_file)
{
  int i, stat;
  FILE *inf, *fopen();
  int int_input, first_int, second_int, int_upper;
  inf = openforreadsilent(int_evaluated_vec_file);
  i = 0;
  for (;;)
    {
      stat = fscanf(inf, "%d", &int_input);
      if (stat==EOF) break;
      if (stat!=1) gabort("get_int_evaluated_vec: read error 1", stat);
      if (*n_int_evaluated >= max_n_int_evaluated)
	gabort("Value of n_int_evaluated too large: get_int_evaluated_vec",
	       get_int_evaluated_vec);
      int_evaluated_vec[*n_int_evaluated] = int_input;
      //      printf("int_evaluated_vec[%d] %d\n",
      //         *n_int_evaluated, int_evaluated_vec[*n_int_evaluated]);
      (*n_int_evaluated)++;
      if (i==0)
	{
	  first_int = int_input;
	  *int_lower = first_int;
	}
      if (i==1)
	{
	  second_int = int_input;
	}
      int_upper = int_input;
      i++;
    }
  fclose (inf);
  *int_step = second_int - first_int;
  //printf("int_lower %d int_upper %d int_step %d\n",
  //  *int_lower, int_upper, *int_step);

}
/******************************************************************************/
void scale_vector(double *vec, int n2, double total_density_s0)
{
  int i;
  double tot = 0;
  //   printf("scale_vector: n2 %d total_density_s0 %lf\n", n2, total_density_s0);
  //   monitorinput();
  for (i=1; i<n2; i++)
    {
      vec[i] /= total_density_s0;
      tot += vec[i];
    }
  vec[0] = 1.0 - tot;        // Put all "missing" density in the lost mutant class
}
/******************************************************************************/
int get_binary_egf_vec(char *buffer, int n, int t2, double s, double *egf_vec)
{
  int i, j, file_size_bytes, n_read, t2_read, offset, tab_num, tab_size,
    size_of_float;
  float s_read, freq_read;
  char *ptr, ch;
  // Each expected gene freq table has two integers (n, t2), a float (s), then
  //   2*n -1 floats (frequencies)
  tab_size = (2*sizeof(int)+sizeof(float) + sizeof(float)*(2*n -1));
  file_size_bytes = n_s_evaluated_file*tab_size;
  //   printf("tab_size %d file_size_bytes %d\n", tab_size, file_size_bytes);
  tab_num =  compute_s_tab_num(s);
  //   printf("tab_num assuming ranges %d\n", tab_num);
  //   monitorinput();
  if (tab_num>=n_s_evaluated_file) return 0;
  //   printf("s %lf s_lower %lf s_step %lf (s - s_lower)/s_step %lf tab_num %d\n",
  //        s, s_lower, s_step, (s - s_lower)/s_step, tab_num);
  //   monitorinput();
  offset = tab_num*tab_size;
  assign_int_from_read_buff(&n_read, buffer, offset);
  //   printf("n_read %d\n", n_read);
  if (n_read!=n)
    {
      // printf("n %d n_read %d\n", n, n_read);
      gabort("get_binary_egf_vec: read error 1", n_read);
    }
  offset += sizeof(int);
  assign_int_from_read_buff(&t2_read, buffer, offset);
  //printf("t2 %d t2_read %d\n", t2, t2_read);
  if (t2_read!=t2) gabort("get_binary_egf_vec: read error 2", t2_read);
  offset += sizeof(int);
  s_read = assign_float_from_read_buff(buffer, offset);
  //   printf("s_read %f\n", s_read);
  size_of_float = sizeof(float);
  for (i=1; i<2*n; i++)
    {
      offset += size_of_float;
      //    assign_float_from_read_buff - in line code for speed
      ptr = &freq_read;
      for (j=0; j<size_of_float; j++)
	{
	  //      printf("assign_float_from_read_buff: ch %d\n", ch);
	  *(ptr + j) = buffer[offset + j];
	}
      //      printf("freq_read %f\n", freq_read);
      egf_vec[i] = (double)freq_read;
    }
  //   dumpvector(egf_vec, 1, 2*n -1, "egf_vec:");
  //   monitorinput();
  return 1;
}
/******************************************************************************/
void fold_vector_double(double *vec, int counter)
{
  int i, j;
  j = 1;
  for (i=counter - 1; i > counter/2; i--)
    {
      vec[j] += vec[i];
      vec[i] = 0;
      j++;

    }
}
/******************************************************************************/
void fold_vector_int(int *vec, int counter)
{
  int i, j;
  j = 1;
  for (i=counter - 1; i > counter/2; i--)
    {
      vec[j] += vec[i];
      vec[i] = 0;
      j++;

    }
}
/******************************************************************************/
int nearest_n2_ind(int n2_real)
{
  int i, ind = undefined_int;
  double min = undefined, diff, res;
  for (i=0; i<n_n2_evaluated; i++)
    {
      diff = n2_evaluated_vec[i] - n2_real;
      if (diff < 0) diff = -diff;

      if (min==undefined)
	{
	  min = diff;
	  ind = i;
	}
      else
	{
	  if (min > diff)
	    {
	      min = diff;
	      ind = i;
	    }
	}
    }
  if (ind == undefined_int) gabort("nearest_n2: ind == undefined_int", ind);
  res = n2_evaluated_vec[ind];;
  //   printf("nearest_n2: n2_real %lf result %lf\n", n2_real, res);
  //   monitorinput();
  return ind;
}
/******************************************************************************/
void read_phase1_phase2_file_into_buffer(int n1, int phase, int n2, int t2,
				    char *buffer,
				    int file_size_bytes)
{
  static char n1_str[max_str_len], n2_str[max_str_len], t2_str[max_str_len];
  char file_str[max_str_len];
  int n2d;
  FILE *inf, *fopen();

  //   printf("read_phase1_phase2_file_into_buffer: phase %d n2 %d t2 %d\n",
  //      phase, n2, t2);
  n2d = 2*n2;
  int_to_string(n1, n1_str, max_str_len);

  int_to_string(n2, n2_str, max_str_len);

  int_to_string(t2, t2_str, max_str_len);

  if (strlen(n1_str) + strlen(n2_str) + strlen(t2_str) >= max_str_len)
    gabort("read_phase1_phase2_file_into_buffer: string too long", max_str_len);
  if (phase==1) 
    {
      strcpy(file_str, phase_1_dir);
      strcat(file_str, "/n1:");
      strcat(file_str, n1_str);
      strcat(file_str, ":");
    }
  else
    {
      strcpy(file_str, phase_2_dir);
      strcat(file_str, "/");
    }
  strcat(file_str, "n2:");
  strcat(file_str, n2_str);
  strcat(file_str, ":t2:");
  strcat(file_str, t2_str);

  inf = openforreadsilent(file_str);
  fread(buffer, 1, file_size_bytes, inf);
  fclose(inf);
}
/******************************************************************************/
void read_const_pop_file_into_buffer(int n1, char *buffer, int file_size_bytes)
{
  char file_str[max_str_len];
  FILE *inf, *fopen();

  //   printf("read_const_pop_file_into_buffer: n1 %d file_size_bytes %d\n",
  //      n1, file_size_bytes);
  strcpy(file_str, const_pop_dir);
  strcat(file_str, const_pop_file);
  //   printf("read_const_pop_file_into_buffer: file_str %s\n", file_str);
  //   monitorinput();
  inf = openforreadsilent(file_str);
  fread(buffer, 1, file_size_bytes, inf);
  fclose(inf);
}
/******************************************************************************/
void compute_weighted_average_egf_vec(double t2_real, int t2_lower, int t2_upper,
				 double *egf_vec_lower, double *egf_vec_upper,
				 double *egf_vec, int n2d)
{
  int i;
  for (i=0; i<=n2d; i++)
    {
      egf_vec[i] = egf_vec_lower[i] + (egf_vec_upper[i] - egf_vec_lower[i])*(t2_real - (double)t2_lower)/(double)(t2_upper - t2_lower);
    }
}
/******************************************************************************/
double compute_weighted_average(double *weighted_vec, double *phase1_mut_vec,
				double *phase2_mut_vec, int n1d, int n2d)
{
  int i;
  double total = 0;
  for (i=1; i<n2d; i++)
    {
      weighted_vec[i]= (n1d*phase1_mut_vec[i] + n2d*phase2_mut_vec[i])/
	(n1d + n2d);
      total += weighted_vec[i];
    }
  //   dumpvector(weighted_vec, 0, n2d, "weighted average vector:"); monitorinput();
  return(total);
}
/******************************************************************************/
int get_const_pop_egf_vec(double s, int n1, double *egf_vec, char *buffer_const_pop,
			  int assign_tds0_mode)
{
  double total_density;
  int i, res;
  if (s<-1)          // The contribution of strongly selected mutations
    // is assumed to be inversely proportional to s, as
    // expected from mutation selection balance, with a
    // contribution from phase 2 mutations only.
    {
      for (i=0; i<=2*n1; i++)
	{
	  egf_vec[i] = 0;
	}
      egf_vec[1] = -2/s;
      //      printf("get_const_pop_egf_vec s %lf\n", s);
      //      dumpvector(egf_vec, 0, 2*n1, "get_const_pop_egf_vec:");
      //      monitorinput();
    }
  else
    {
      res =  get_binary_egf_vec(buffer_const_pop, n1, 0, s, egf_vec);
      if (res==0) return 0;
      //      dumpvector(egf_vec, 0, 2*n1, "get_const_pop_egf_vec:");
      //      monitorinput();
    }
  total_density = compute_weighted_average(egf_vec, egf_vec, egf_vec, 2*n1, 2*n1);
  if (assign_tds0_mode)
    {
      total_density_s0 = total_density;
    }
  //   printf("assign_tds0_mode %d: total_density_s0 %lf\n", assign_tds0_mode, total_density_s0);
  return (1);
}
/******************************************************************************/
void assign_int_from_read_buff(int *res, char *buffer, int start)
{
  int temp, i;
  char *ptr, ch;
  ptr = &temp;
  for (i=0; i<sizeof(int); i++)
    {
      ch = buffer[start+i];

      *(ptr + i) = ch;
    }

  *res = temp;
}
/******************************************************************************/
float assign_float_from_read_buff(char *buffer, int start)
{
  float temp;
  int i;
  char *ptr, ch;
  ptr = &temp;
  for (i=0; i<sizeof(float); i++)
    {
      ch = buffer[start+i];
      //      printf("assign_float_from_read_buff: ch %d\n", ch);
      *(ptr + i) = ch;
    }
  //   printf("assign_float_from_read_buff: temp %f\n", temp);
  return temp;
}
/******************************************************************************/
int compute_s_tab_num(double s)
{
  double s_lower, s_upper, s_step;
  int i, tab_num;
  for (i=1; i<=nranges; i++)
    {
      s_lower = s_ranges[i].s_lower;
      s_upper = s_ranges[i].s_upper;
      if ((s >= s_lower)&&(s <= s_upper))
	{
	  s_step = s_ranges[i].s_step;
	  tab_num = (int)(floor((s - s_lower)/s_step + 0.0000000001));
	  return(tab_num + s_ranges[i].s_lower_index - 1); 
	}
    }
  gabort("compute_s_tab_num: failure to find tab_num", 0);
  return(1);
}
/******************************************************************************/
int compute_file_size_bytes(int n)
{
  int tab_size, file_size_bytes;
  tab_size = (2*sizeof(int)+sizeof(float) + sizeof(float)*(2*n - 1));
  file_size_bytes = n_s_evaluated_file*tab_size;
  return (file_size_bytes);
}
/******************************************************************************/
void set_up_file_name(int n1, char *froot, char *file_name)
{
  int len;
  static char n1_str[max_str_len];
  int_to_string(n1, n1_str, max_str_len);
  //   printf("set_up_file_name n1 %d froot %s n1_str %s\n", n1, froot, n1_str);
  //   monitorinput();
  len = strlen(data_path);
  len += strlen(n1_str) + 1;
  len += strlen(froot);
  if (len >= max_str_len) gabort("set_up_file_name: string too long (1)", len);
  strcpy(file_name, data_path);
  strcat(file_name, "n1_");
  strcat(file_name, n1_str);
  strcat(file_name, "/");
  strcat(file_name, froot);
  //   printf("set_up_file_name n1 file_name %s\n", file_name);
  //   monitorinput();
}
/******************************************************************************/
void get_data_path(char *data_path)
{
  FILE *inf, *fopen();
  int stat;
  inf = openforread("directory_config.dat");
  stat = fscanf(inf, "%s", data_path);
  if (stat!=1) gabort("get_data_path: read error", 1);
  //   printf("get_data_path: data_path %s\n", data_path); monitorinput();
}
/******************************************************************************/
void get_s_evaluated_vec(double *s_evaluated_vec, int *n_s_evaluated,
		    int *n_s_evaluated_file, char *s_evaluated_vec_file)
{
  int stat;
  FILE *inf, *fopen();
  double s;
  *n_s_evaluated = 0;
  inf = openforreadsilent(s_evaluated_vec_file);
  s = -100.0;                                 // Set up large |s| values specially
  while (s < -1)
    {
      s_evaluated_vec[*n_s_evaluated] = s;
      (*n_s_evaluated)++;
      if (s > -20) s = s + 1.0; else s = s + 5.0;
    }
  *n_s_evaluated_file = 0;
  for (;;)
    {
      stat = fscanf(inf, "%lf", &s);
      if (stat==EOF) break;
      if (stat!=1) gabort("get_s_limits: read error 1", stat);
      if (*n_s_evaluated >= max_n_s_evaluated)
	gabort("Value of n_s_evaluated too large: get_s_evaluated_vec",
	       get_s_evaluated_vec);
      s_evaluated_vec[*n_s_evaluated] = s;
      (*n_s_evaluated)++;
      (*n_s_evaluated_file)++;
    }
  fclose (inf);
  //   printf("n_s_evaluated %d n_s_evaluated_file %d\n", *n_s_evaluated,
  //      *n_s_evaluated_file);
  //   monitorinput();
}
/******************************************************************************/
void get_lower_upper(int ind, double *s_evaluated_vec, int n_s_evaluated,
		double *lower, double *upper)
{
  double x, temp_r;
  if ((ind<0)||(ind>=n_s_evaluated))
    gabort("get_lower_upper: ind out of range", ind);
  x = s_evaluated_vec[ind];
  if (ind<n_s_evaluated-1)
    *upper = (x + s_evaluated_vec[ind+1])/2.0;
  else
    *upper = x;
  if (ind!=0)
    *lower = (x + s_evaluated_vec[ind-1])/2.0;
  else
    *lower = x;
  if (x <=0.0)
    {
      temp_r = *lower;
      *lower = *upper;
      *upper = temp_r;
      if (*lower!=0.0) *lower = -*lower;
      if (*upper!=0.0) *upper = -*upper;
    }
}
/******************************************************************************/
void dumpvector(double *v, int min, int max, char *s)
{
  int i, control = 0;
  //   printf("\n");
  printf("%s", s); printf("\n");
  for (i = min; i<max; i++)
    {
      printf("%3d %10.8lf ", i, v[i]);
      control++;
      if (control==5)
	{
	  control = 0;
	  printf("\n");
	}
    }
  printf("\n");
}
/******************************************************************************/
void dumpmatrix(double **m, int s1,
	   int rows, int s2, int cols, char *s)
{
  int i, j;
  printf("\n%s\n", s);
  for (i = s1; i<=rows; i++)
    {
      for (j = s2; j<=cols; j++)
	{
	  printf(" %2.3f", m[i][j]);
	}
      printf("\n");
    }
  printf("\n");
}
/******************************************************************************/
double selfn(double q, double s, double h)
{
  double num, denom;
  num = s*q*(1-q)*(q + h*(1-2*q));
  denom = 1 + 2.0*h*s*q*(1-q) + s*q*q;
  if (denom== 0) denom = .00000000000001;
  return( q + (num/denom));
}
/******************************************************************************/
void setuprow(double **a, double p, int row, int n)
{

  double z;         /*probability of zero failures*/
  double x, y, temp1, temp2;
  int j;

  temp1 = 1.0 - p;
  if (temp1 == 0.0)  temp1 = .000000000001;            /*prevent zero divide*/
  temp1 = p/temp1;       /*First prevent negative log*/
  if (temp1 <=0)   temp1 = .000000000001;
  x = log(temp1);    /*used for later multiplication of series*/
  temp2 = 1.0-p;
  if (temp2<=0)   temp2 = .000000000001;
  z = (double)n*log(temp2);
  if (z > 150)   z = 150;
  a[row][0] = exp(z);
  y = a[row][0];        /*add up for later*/
  for (j = 1;j<=n; j++)
    {

      z = z + log(((double)n + 1.0 - (double)j)/(double)j) + x;
      /*             cumulative binomial coeff.*/
      a[row][j] = exp(z);
      y = y + a[row][j];          /*add up*/
    }
  if (y == 0.0)   y = .0000000000001;
  for (j=0; j<=n; j++) a[row][j] = a[row][j]/y;
}
/******************************************************************************/
void setupmatrix(double **a, int n1, int n2, double s, double h)
{
  int i;
  double p;


  for (i = 1; i<=n1 - 1; i++)             /*do variable rows*/
    {
      p = (double)i/(double)n1;
      p = selfn(p, s,h);
      setuprow(a, p, i, n2);
    }
  for (i = 0; i<=n2; i++)
    {              /*do constant edges*/
      a[0][i] = 0;
      a[n1][i] = 0;
    }
  a[0][0] = 1.0;                          /*do absorbing corners*/
  a[n1][n2] = 1.0;

}
/******************************************************************************/
void matrixinvert(int N, double **a,double *inverter)
{

  int i, j, lotkin_signum;
  gsl_matrix  *lotkin_a,*lotkin_inv;
  gsl_vector *x, *lotkin_b, *lotkin_x;

  gsl_permutation *lotkin_perm;

  /* allocate a, x, b */
  lotkin_a = gsl_matrix_alloc(N-1, N-1);
  lotkin_inv = gsl_matrix_alloc(N-1, N-1);
  x = gsl_vector_alloc(N-1);
  lotkin_b = gsl_vector_alloc(N-1);

  lotkin_x = gsl_vector_alloc(N-1);

  gsl_matrix *identity=gsl_matrix_alloc(N-1, N-1);


  for (i = 0; i<=N-2; i++)
    {
      for (j = 0; j<=N-2; j++)
	{
          gsl_matrix_set( identity, i, i,1);
	}
      ;
    }

  for (i = 0; i<=N-2; i++)
    {
      for (j = 0; j<=N-2; j++)
	{
	  gsl_matrix_set(lotkin_a,i,j,gsl_matrix_get(identity,i,j)-a[i+1][j+1]);
	}

    }


  /* LU decomposition and forward&backward substition */

  //gsl_blas_dgemv(CblasNoTrans, 1.0, lotkin_a, x, 0.0, lotkin_b);

  lotkin_perm = gsl_permutation_alloc(N-1);
  gsl_linalg_LU_decomp(lotkin_a, lotkin_perm, &lotkin_signum);
  gsl_linalg_LU_solve(lotkin_a, lotkin_perm, lotkin_b, lotkin_x);
  gsl_linalg_LU_invert (lotkin_a,lotkin_perm, lotkin_inv);

  /*apparently the inversion adds +1 to the first element of the vector
    so I substract it*/

  gsl_matrix_set(lotkin_inv,0,0,gsl_matrix_get(lotkin_inv,0,0)-1);

  //printf("\n%s\n", "Equilibrium");

  for (i = 0; i<1; i++)
    {
      for (j = 0; j<=N-2; j++)
	{

	  // printf("%d: %2.3f  ",j,gsl_matrix_get(lotkin_inv, i, j));
	}
      //printf("\n");
    }


  for (i=0;i<=N-2;i++){

    inverter[i+1]=gsl_matrix_get(lotkin_inv,0,i);

  }

  gsl_matrix_free(lotkin_a);
  gsl_matrix_free(lotkin_inv);
  gsl_vector_free(x);
  gsl_vector_free(lotkin_b);
  gsl_vector_free(lotkin_x);
  gsl_matrix_free(identity);
  gsl_permutation_free(lotkin_perm);
}
/******************************************************************************/
void tmiterate(double **a,double *mut1,
	  int t, int n1, int n2,int decay,double *sumf)
{

  int k, i, j = 0;
  double z;

  double * steadystatefreq = (double*) calloc (maxnd+1, sizeof(double));
  double * mut2 = (double*) calloc (maxnd+1, sizeof(double));


  for (k =1 ;k<=t; k++)
    {                 /*t iterations*/

      /* Increment steadystatefreq. distribution*/
      for (i = 0; i<=n1; i++)
	{
	  steadystatefreq[i] = steadystatefreq[i] + mut1[i];
	}

      /* Perform matrix multiplication*/
      for (i = 0; i<=n2; i++)
	{
	  z = 0;
	  for (j = 0; j<=n1; j++)
	    {
	      z = z + a[j][i]*mut1[j];
	    }
	  mut2[i] = z;
	}

      /* Copy result of multiplication*/
      for (i=0; i<=n2; i++) {mut1[i] = mut2[i];

	if (decay==0){
	  sumf[i]+=mut1[i];}
	if(decay==1){sumf[i]=mut1[i];}
      }


    }

  free(steadystatefreq);
  free(mut2);
}
/******************************************************************************/
void tmiterate2(double **a,double *mut1,
	   int t, int n1, int n2,int decay,double **a1)
{

  double *  steadystatefreq = (double*) calloc (maxnd+1, sizeof(double));
  double * mut2 = (double*) calloc (maxnd+1, sizeof(double));


  int k, i, j = 0;
  double z;



  for (k = 1; k<=t; k++){                 /*t iterations*/

    /* Increment steadystatefreq. distribution*/
    for (i = 0; i<=n1; i++)
      {
	steadystatefreq[i] = steadystatefreq[i] + mut1[i];
      }

    /* Perform matrix multiplication*/
    for (i = 0; i<=n2; i++)
      {
	z = 0;
	for (j = 0; j<=n1; j++)
	  {
            z = z + a[j][i]*mut1[j];
	  }
	mut2[i] = z;
      }

    /* Copy result of multiplication*/
    for (i=0; i<=n2; i++) {mut1[i] = mut2[i];

      if (decay==0){

	a1[k][i]+=mut1[i];
	if (k>1){a1[k][i]+=a1[k-1][i];}
      }

      if(decay==1){
          
	a1[k][i]=mut1[i];
	//if (k>0){ a1[k][i]+=a1[k-1][i];}
      }
    }

    //dumpvector(a1[k],0,n2,"matrix");
  }

  free(steadystatefreq);
  free(mut2);

}
/******************************************************************************/
void eqf_using_matrix_inversion(int n1,double s,double **a,double *egf_out)
{
  setupmatrix(a, n1, n1, s, H);
  matrixinvert(n1,a,egf_out);
}
/******************************************************************************/
void vector_average(int n1,int n2,double *fv1,double *fv2, double *fv_out)
{
  int i;
  for (i=0;i<n2;i++)
    {

      fv_out[i]=((n1*fv1[i])+(n2*fv2[i]))/(n1+n2);

    }
}
/******************************************************************************/
void vector_s_average(int n,double P1,double *fv1,double *fv2, double *fv_out)
{
  int i;
  for (i=0;i<n;i++)
    {

      fv_out[i]=(P1*fv1[i])+((1-P1)*fv2[i]);

    }
}
/******************************************************************************/
void output_sfs_to_file_thanasis_format(int n,int sample1, int *sfs1,int *sfs2,char *sfs_filename)
{


	  FILE *file= fopen(sfs_filename, "w" );
	  int i;

	  //selected sites first
	  fprintf(file,"%d\n",sample1);

	  for (i=1;i<n;i++){

	    fprintf(file,"%d ",sfs2[i]);

	  }

	  fprintf(file,"\n");
	  //neutral sites second
	  for (i=1;i<n;i++){

	    fprintf(file,"%d ",sfs1[i]);


	  }

	  fclose(file);
}
/******************************************************************************/
void output_sfs_to_file_peter_format(int n,int sample1, int sample2, int *sfs1,int *sfs2,char *sfs_filename)
{

  FILE *file= fopen(sfs_filename, "w" );
  int i;
  fprintf(file,"1\n");
  fprintf(file,"%d\n",n);

  //selected sites first
  fprintf(file,"%d\n",sample2);

  for (i=1;i<n;i++){

    fprintf(file,"%d ",sfs2[i]);

  }

  fprintf(file,"\n");
  //neutral sites second

  fprintf(file,"%d\n",sample1);
  for (i=1;i<n;i++){

    fprintf(file,"%d ",sfs1[i]);


  }

  fclose(file);

}
/******************************************************************************/
void output_sfs_to_file_peter_format2(int n,int sample1, int sample2, int *sfs1,int *sfs2,char *sfs_filename)
{

  FILE *file= fopen(sfs_filename, "w" );
  int i;
  fprintf(file,"name\n");

  fprintf(file,"0 0\n0 0\n");
  fprintf(file,"%d\n",n);
  for (i=0;i<n;i++){

    fprintf(file,"%d ",sfs1[i]);

  }

  fprintf(file,"0\n");
  //fprintf(file,"%d\n",sample1);
  for (i=0;i<n;i++){

    fprintf(file,"%d ",sfs2[i]);


  }
fprintf(file,"0");
  fclose(file);

}
/******************************************************************************/
void get_sfs(int *nalleles,float *sfs_sel,float *sfs_neu,char *sfs_filename)
{

  int i=0;
  int gene_number=0,sum=0;
  FILE *file= fopen(sfs_filename, "r" );

  fscanf (file, "%d", nalleles);

  //printf("%s %d %d %d %d %d\n",name,sel_sites,sel_diff,neu_sites,neu_diff,*nalleles);

  fscanf (file, "%f", &sfs_sel[0]);
  //printf("%d\n",sfs_sel[0]);
  sum=0;
  for (i =1 ;i<=*nalleles; i++)
    {
      fscanf (file, "%f", &sfs_sel[i]);
      //printf("%f\n",sfs_sel[i]);
      sum+=sfs_sel[i];
    }

//printf("next\n");
  fscanf (file, "%f", &sfs_neu[0]);
  //printf("%d\n",sfs_neu[0]);
  sum=0;
  for (i =1 ;i<=*nalleles; i++)
    {
      fscanf (file, "%f", &sfs_neu[i]);
      //printf("%d\n",sfs_neu[i]);
      sum+=sfs_neu[i];
    }

  //monitorinput();
  fclose(file);

}
/******************************************************************************/
void get_sfs_peter1(int *total_neu,int *total_sel,int *nalleles,float *sfs_sel,float *sfs_neu,char *sfs_filename)
{
  int i=0;
  int gene_number=0,sum=0;
  
  FILE *file= fopen(sfs_filename, "r" );
  fscanf (file, "%d", &gene_number);
  fscanf (file, "%d", nalleles);

  //printf("%d %d\n",read_sfs,*nalleles);

  fscanf (file, "%f", &sfs_sel[0]);
  //printf("%d\n",sfs_sel[0]);
  sum=0;
  for (i =1 ;i<=*nalleles; i++)
    {
      fscanf (file, "%f", &sfs_sel[i]);
      //printf("%f\n",sfs_sel[i]);
      sum+=sfs_sel[i];
    }

//printf("next\n");
  fscanf (file, "%f", &sfs_neu[0]);
  //printf("%d\n",sfs_neu[0]);
  sum=0;
  for (i =1 ;i<=*nalleles; i++)
    {
      fscanf (file, "%f", &sfs_neu[i]);
      //printf("%d\n",sfs_neu[i]);
      sum+=sfs_neu[i];
    }

  //monitorinput();
  fclose(file);
}
/******************************************************************************/
void get_sfs_peter2(char *sfsname,float *total_neu,float *total_sel,int *nalleles,float *sfs_sel,float *sfs_neu,char *sfs_filename,float *sel_sites,float *sel_diff,float *neu_sites,float *neu_diff)
{

  int i=0;
  int gene_number=0,sum=0;
  FILE *file= fopen(sfs_filename, "r" );
  
  fscanf (file, "%s", sfsname);
  fscanf (file, "%f", sel_sites);
  fscanf (file, "%f", sel_diff);
  fscanf (file, "%f", neu_sites);
  fscanf (file, "%f", neu_diff);
  fscanf (file, "%d", nalleles);

  //printf("%s %d %d %d %d %d\n",name,sel_sites,sel_diff,neu_sites,neu_diff,*nalleles);

  fscanf (file, "%f", &sfs_sel[0]);
  //printf("%d\n",sfs_sel[0]);
  sum=0;
  for (i =1 ;i<=*nalleles; i++)
    {
      fscanf (file, "%f", &sfs_sel[i]);
      //printf("%f\n",sfs_sel[i]);
      sum+=sfs_sel[i];
    }

//printf("next\n");
  fscanf (file, "%f", &sfs_neu[0]);
  //printf("%d\n",sfs_neu[0]);
  sum=0;
  for (i =1 ;i<=*nalleles; i++)
    {
      fscanf (file, "%f", &sfs_neu[i]);
      //printf("%d\n",sfs_neu[i]);
      sum+=sfs_neu[i];
    }

  //monitorinput();
  fclose(file);
}

/******************************************************************************/
double sfsfold_f(int nalleles,float *sfs,float *sfs_folded,int index)
{

int i=0;

      for (i=0;i<=nalleles/2;i++)
	{
	  sfs_folded[index+i]=sfs[i]+sfs[nalleles-i]; 
	  
	  if(i==nalleles/2&&(!odd(nalleles)))
	  {
	  sfs_folded[index+i]=sfs[i];
	  }
	  
          //printf("%d:%d\n",sfs[i],sfs_folded[index+i]);
	}
      return(1);
}

/******************************************************************************/
int binarysearchcumulative_from_zero(double *array, int size)
{
  int cur=0, pre=0, post=0;
  double r=0;
  r = uniform();
  /*   printf("Uniform %f\n", r);*/
  pre = 0;
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
/******************************************************************************/
double randomly_select_allele_frequency(double *egf_vec, double *temp1, int nd)
{
  double u=0, freq=0;
  int ind=0;
  ind = binarysearchcumulative_from_zero(egf_vec, nd-1);
  //   printf("ind %d\n", ind);
  (temp1[ind])++;
  freq = (double)ind/(double)nd;
  //   printf("ind %d freq %lf\n", ind, freq); monitorinput();
  return(freq);
}
/******************************************************************************/
void set_up_cumulative_allele_freq_vec(int nd, double *egf_vec, double *cum_vec,
				  char *str)
{
  int i;
  for (i=0; i<nd; i++)
    {
      cum_vec[i] = egf_vec[i];
    }
  for (i=1; i<=nd; i++)
    {
      cum_vec[i] += cum_vec[i-1];
    }
}
/******************************************************************************/
void generate_simulated_data(double *egf_vec, int nsites, int sample_size, int nd,
			int *sfs)
{
  int i, n_seg;
  double allele_freq;
  struct acc a;
  initacc(&a);
  for (i=1; i<=nsites; i++)
    {
      allele_freq = 
	randomly_select_allele_frequency(egf_vec, allele_freq_sampled, nd);
      n_seg = bnldev(allele_freq, sample_size);
      accum(&a, (double)n_seg);
      (sfs[n_seg])++;
    }
}

/******************************************************************************/
void egf_scaling_s(int n1,double f0,double *v_s_in, double *v_s_out)
{

  double sumx,sumy;
  int i;
  sumx=0;
  sumy=0;

  for (i=1;i<n1;i++){

    sumx+=v_s_in[i];

  }
//printf("%f\n",sumx);
  for (i=1;i<n1;i++){
    v_s_out[i]/=sumx;
    sumy+=v_s_out[i];
  }

  v_s_out[0]=(1-sumy);
}
/******************************************************************************/
void egf_scaling_f0(int n1,double f0,double *fv)
{

  int i;
  for (i=0;i<n1;i++){
    fv[i]*=(1-f0);
  }
  fv[0]+=f0;
}
/******************************************************************************/
void egf_scaling(int N,double f0,double *FV0,double *FVS)
{
  int n2d=N*2;
  egf_scaling_s(n2d,f0,FV0, FVS);
  egf_scaling_f0(n2d,f0,FVS);
}
/******************************************************************************/
void binomial_sampling(int n1,int alleles,int sample, double *invy,int *discrete)
{

  int s1=0,success=0;
  int i=0;
  int n1d=2*n1;
  double prob=0;


  const gsl_rng_type * T;
  T = gsl_rng_taus;
  gsl_rng  *rgen =gsl_rng_alloc(T);

  /*
    based on equilibrium frequency vector (use it as a probability vector),
    I randomly generate numbers from 0 to n1d-1 for sample sites.
  */


  for (i=0;i<sample;i++)
    {

      gsl_ran_discrete_t *r= gsl_ran_discrete_preproc (n1d,invy);
      s1=gsl_ran_discrete (rgen, r);
      prob=(double)(s1)/n1d;

      success=gsl_ran_binomial(rgen,prob,alleles-1);
      discrete[success]++;

      gsl_ran_discrete_free(r);
    }
  gsl_rng_free (rgen);
  

}
/******************************************************************************/
int set_up_density_vec_equal_effects(double s, int n_s_evaluated, double *density_vec)
{
   double s_lower=0, s_upper=0, p_upper=0;
   int i=0, s_ind=0;
   s_lower = undefined;
   s_upper = undefined;
//   printf("set_up_gamma_density_vec_equal_effects s %lf\n", s);
//   monitorinput();
   for (i=0; i<n_s_evaluated; i++)
   {
      density_vec[i] = 0.0;
   }
   for (i=1; i<n_s_evaluated; i++)
   {
//      printf("%d %lf  ", i, s_evaluated_vec[i]);
//      printf("\n");
      if (s_evaluated_vec[i] >= s)
      {
         s_upper = s_evaluated_vec[i];
         s_ind = i;
         if (i != 0 ) s_lower = s_evaluated_vec[i-1];
         break;
      }
   }
   if ((s_lower==undefined)||(s_upper==undefined)) return 0;
   p_upper = (s - s_lower)/(s_upper - s_lower);
//   printf("set_up_density_vec_equal_effects: s_lower %lf s_upper %lf\n",
//      s_lower, s_upper);
//   printf("p_upper %lf\n", p_upper);
//   monitorinput();
   density_vec[s_ind] = p_upper;
   if (s_ind==1)
      density_vec[s_ind] += (1.0 - p_upper);
   else
      density_vec[s_ind-1] = (1.0 - p_upper);
//   for (i=1; i<n_s_evaluated; i++)
//   {
//      printf("s %lf density_vec[%d] %lf\n", s_evaluated_vec[i], i,
//         density_vec[i]);
//   }
//   monitorinput();
return (1);
}
/******************************************************************************/
int set_up_density_vec_step_effects(double alpha, double beta, 
				    int n_s_evaluated, double *density_vec)
{
 if (alpha<0){alpha=-alpha;}
  if (beta<0){beta=-beta;}
  int i=0,s_ind=0;
  double lower=0, upper=0, area=0, total_area = 0.0, x=0, diff=0, step=0, inc_x=0, inc_y=0, gf=0;
  double mean=0;
  
  for (i=0; i<n_s_evaluated; i++)
    {
 
      x = s_evaluated_vec[i];

      get_lower_upper(i, s_evaluated_vec, n_s_evaluated, &lower, &upper);
      
      /*calculate cumulative uniform distribution for lower and upper*/
      inc_x=gsl_cdf_flat_P(lower,alpha,beta);    
      inc_y=gsl_cdf_flat_P(upper,alpha,beta);
 
      area=inc_y-inc_x;
      density_vec[i] = area;
      total_area+=area;
      //printf("x %f inc_x %f inc_y %f area %f\n",x,inc_x,inc_y,area);
    }
//printf("alpha:%f beta:%f total_area:%f\n",alpha,beta, total_area);
//monitorinput();
// area =gsl_cdf_flat_P(-s_evaluated_vec[0], alpha, beta);
// density_vec[0]+=1.0-area;
 
 //if (density_vec[0]>0){printf("%f\n",density_vec[0]);}
  return (1);
}
/******************************************************************************/
double compute_beta_densities(double alpha,double beta,double beta_scaler, double *s_evaluated_vec,int n_s_evaluated, double *density_vec)
{
  int i=0;
  double lower=0, upper=0, area=0, total_area = 0.0, x=0, diff=0, step=0, inc_x=0, inc_y=0, gf=0;

  if (alpha <= 0) return undefined;
  if (beta <= 0) return undefined;
  
  for (i=0; i<n_s_evaluated; i++)
    {
      x = s_evaluated_vec[i];
      get_lower_upper(i, s_evaluated_vec, n_s_evaluated, &lower, &upper);

      inc_x = gsl_cdf_beta_P (lower, alpha, beta);
      inc_y = gsl_cdf_beta_P (upper, alpha, beta);
      
      area = (inc_y - inc_x);
      density_vec[i] = area;

      //printf("x %f inc_x %f inc_y %f area %f\n",x,inc_x,inc_y,area);
    
    }
   area = gsl_cdf_beta_P(-s_evaluated_vec[0], alpha, beta);
   density_vec[0]+=1.0-area;


  return (1);
}
/******************************************************************************/
double compute_gamma_densities(double alpha,double beta, double *s_evaluated_vec, int n_s_evaluated, double *density_vec)
{
  int i=0;
  double lower=0, upper=0, area=0, total_area = 0.0, x=0, diff=0, step=0, inc_x=0, inc_y=0, gf=0;
  
  if (alpha <= 0) return undefined;

  for (i=0; i<n_s_evaluated; i++)
    {
      x = s_evaluated_vec[i];
      get_lower_upper(i, s_evaluated_vec, n_s_evaluated, &lower, &upper);

      inc_x = gsl_cdf_gamma_P(lower, beta,1/alpha);
      inc_y = gsl_cdf_gamma_P(upper, beta,1/alpha);
      
      area = (inc_y - inc_x);
      density_vec[i] = area;
      total_area+=area;
      //printf("x%f 1/alpha %f beta %f area %f\n",x,1/alpha,beta,area);
      
    }

  area = gsl_cdf_gamma_P(-s_evaluated_vec[0], beta,1/alpha);
  density_vec[0]+=1.0-area;
  //printf("x%f inc_x %f inc_y %f area %f\n",x,inc_x,inc_y,area);


  return (1);

}

/******************************************************************************/
double compute_lognormal_densities(double alpha, double beta,double *s_evaluated_vec,int n_s_evaluated, double *density_vec)
{
  int i=0;
  double lower=0, upper=0, area=0, total_area = 0.0, x=0, diff=0, step=0, inc_x=0, inc_y=0, gf=0;
  

  for (i=0; i<n_s_evaluated; i++)
    {
      x = s_evaluated_vec[i];
      get_lower_upper(i, s_evaluated_vec, n_s_evaluated, &lower, &upper);

      inc_x = gsl_cdf_lognormal_P(lower, alpha,beta);
      inc_y = gsl_cdf_lognormal_P(upper, alpha,beta);
      
      area = (inc_y - inc_x);
      density_vec[i] = area;
      total_area+=area;
      //printf("x%f 1/alpha %f beta %f area %f\n",x,1/alpha,beta,area);
      
    }

  area = gsl_cdf_lognormal_P(-s_evaluated_vec[0], alpha,beta);
  density_vec[0]+=1.0-area;
  //printf("x%f inc_x %f inc_y %f area %f\n",x,inc_x,inc_y,area);


  return (1);

}



/******************************************************************************/
double calculate_ne(int n1,int n2,int t)
{
double w1=0,w2=0,n_e=0;
double N1=(double)n1;
double N2=(double)n2;
double t2=(double)t;
/*compute weighted Ne*/
w1=pow((1-1/(2*N2)),t2);
w2=(1-exp(-t2/(2*N2)));

n_e=(N1*w1+N2*w2)/(w1+w2);

return n_e;
}

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

/******************************************************************************/
#define H 0.5
#define maxnd 10000

#define neutral 0
#define selected 1
#define max_n_int_evaluated 1000
#define max_n_s_evaluated 2000

#define undefined -99999999999.99
#define undefined_int 999999
#define max_str_len 100
#define n_decimal_places 6
#define max_s_ranges 100  

#define s_evaluated_vec_file_const "s_evaluated.dat"
#define t2_evaluated_vec_file_const "t2_evaluated.dat"
#define n2_evaluated_vec_file_const "n2_evaluated.dat"
#define s_range_file_const "s_ranges_evaluated.dat"
#define phase_1_dir_const "n1_n2_t2_s_vec_dir"
#define phase_2_dir_const "n2_t2_s_vec_dir"
#define const_pop_file "control_generate_allele_distrib_n1.out.bin"
/******************************************************************************/
/*Function List*/
double compute_weighted_average(double *weighted_vec, double *phase1_mut_vec,
				double *phase2_mut_vec, int n1d, int n2d);
float assign_float_from_read_buff(char *buffer, int start);
int get_average_egf_vecs_1_and_2(double s,int n1, int n2, double t2_real,
				 int t2_lower, int t2_upper, double *egf_vec, char *buffer_p1_n2_t2_lower,
				 char *buffer_p2_n2_t2_lower, char *buffer_p1_n2_t2_upper,
				 char *buffer_p2_n2_t2_upper);
				 get_upper_lower_int(double parm, int *lower, int *upper, int n_int_evaluated, int *int_evaluated_vec);
get_upper_lower_double(double parm, double*lower, double *upper,int n_int_evaluated, double *double_evaluated_vec);
get_s_ranges();
get_int_evaluated_vec(int *int_evaluated_vec, int *n_int_evaluated,
		      int *int_lower, int *int_step, char *int_evaluated_vec_file);
scale_vector(double *vec, int n2, double total_density_s0);
int get_binary_egf_vec(char *buffer, int n, int t2, double s, double *egf_vec);

fold_vector_double(double *vec, int counter);
fold_vector_int(int *vec, int counter);

int nearest_n2_ind(int n2_real);
read_phase1_phase2_file_into_buffer(int n1, int phase, int n2, int t2, char *buffer,
				    int file_size_bytes);
read_const_pop_file_into_buffer(int n1, char *buffer, int file_size_bytes);
compute_weighted_average_egf_vec(double t2_real, int t2_lower, int t2_upper,
				 double *egf_vec_lower, double *egf_vec_upper,
				 double *egf_vec, int n2d);
double compute_weighted_average(double *weighted_vec, double *phase1_mut_vec,
				double *phase2_mut_vec, int n1d, int n2d);
int get_const_pop_egf_vec(double s, int n1, double *egf_vec, char *buffer_const_pop,
			  int assign_tds0_mode);
assign_int_from_read_buff(int *res, char *buffer, int start);
float assign_float_from_read_buff(char *buffer, int start);
int compute_s_tab_num(double s);
int compute_file_size_bytes(int n);
set_up_file_name(int n1, char *froot, char *file_name);
get_data_path(char *data_path);
get_s_evaluated_vec(double *s_evaluated_vec, int *n_s_evaluated,
   int *n_s_evaluated_file, char *s_evaluated_vec_file);
get_lower_upper(int ind, double *s_evaluated_vec, int n_s_evaluated,
   double *lower, double *upper);
dumpvector(double *v, int min, int max, char *s);
dumpmatrix(double **m, int s1,
        int rows, int s2, int cols, char *s);
double selfn(double q, double s, double h);
setuprow(double **a, double p, int row, int n);
setupmatrix(double **a, int n1, int n2, double s, double h);
matrixinvert(int N, double **a,double *inverter);
tmiterate(double **a,double *mut1,
        int t, int n1, int n2,int decay,double *sumf);
tmiterate2(double **a,double *mut1,
        int t, int n1, int n2,int decay,double **a1);
eqf_using_matrix_inversion(int n1,double s,double **a,double *egf_out);
vector_average(int n1,int n2,double *fv1,double *fv2, double *fv_out);
vector_s_average(int n,double P1,double *fv1,double *fv2, double *fv_out);

output_sfs_to_file_thanasis_format(int n, int *sfs1,int *sfs2,char *sfs_filename);
output_sfs_to_file_peter_format(int n,int sample1, int sample2, int *sfs1,int *sfs2,char *sfs_filename);
get_sfs(int alleles,int *par,char *sfs_filename);

get_sfs_peter1(int *total_neu,int *total_sel,int *nalleles,float *sfs_sel,float *sfs_neu,char *sfs_filename);
get_sfs_peter2(char *sfsname,float *total_neu,float *total_sel,int *nalleles,float *sfs_sel,float *sfs_neu,char *sfs_filename,float *sel_sites,float *sel_diff,float *neu_sites,float *neu_diff);

double sfsfold_f(int nalleles,float *sfs,float *sfs_folded,int index);


int binarysearchcumulative_from_zero(double *array, int size);
double randomly_select_allele_frequency(double *egf_vec, double *temp1, int nd);
set_up_cumulative_allele_freq_vec(int nd, double *egf_vec, double *cum_vec,
				  char *str);
generate_simulated_data(double *egf_vec, int nsites, int sample_size, int nd,
			int *sfs);
			
egf_scaling_s(int n1,double f0,double *v_s_in, double *v_s_out);
egf_scaling_f0(int n1,double f0,double *fv);
egf_scaling(int N,double f0,double *FV0,double *FVS);
binomial_sampling(int n1,int alleles,int sample, double *invy,int *discrete);


int set_up_density_vec_equal_effects(double s, int n_s_evaluated, double *gamma_density_vec);
int set_up_density_vec_step_effects(double alpha, double beta, int n_s_evaluated, double *gamma_density_vec);
double compute_beta_densities(double alpha,double beta,double beta_scaler, double *s_evaluated_vec,int n_s_evaluated, double *density_vec);
double compute_gamma_densities(double alpha,double beta, double *s_evaluated_vec, int n_s_evaluated, double *density_vec);				
double compute_lognormal_densities(double alpha, double beta,double *s_evaluated_vec,int n_s_evaluated, double *density_vec);
			       				
double calculate_ne(int n1,int n2,int t);
			       			
/******************************************************************************/
/*Structure List*/
struct rangestruct
{
  double s_lower, s_upper, s_step;
  int s_lower_index;
};

struct rangestruct s_ranges[max_s_ranges+1];
/******************************************************************************/
/*Global variable List*/
double allele_freq_sampled[maxnd+1];
double total_density_s0;
char s_evaluated_vec_file[max_str_len], t2_evaluated_vec_file[max_str_len],
  n2_evaluated_vec_file[max_str_len], s_range_file[max_str_len],
  phase_1_dir[max_str_len], phase_2_dir[max_str_len],
  const_pop_dir[max_str_len];
char data_path[max_str_len], file_label_str[max_str_len];
int nranges, trace_level, n_sfs;
int n_s_evaluated, n_s_evaluated_file, gss_n2;
int t2_evaluated_vec[max_n_int_evaluated], n_t2_evaluated, t2_lower;
int t2_step, verbose_mode, neutrals_only_mode;
int n2_evaluated_vec[max_n_int_evaluated], n_n2_evaluated, n2_lower;
double s_evaluated_vec[max_n_s_evaluated];
/******************************************************************************/

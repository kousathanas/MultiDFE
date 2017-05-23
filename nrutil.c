#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#define NR_END 1
#define FREE_ARG char*

void nrerror(char error_txt[])
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_txt);
	fprintf(stderr,"...now exiting to system\n");
	exit(1);
}

float *vector(long nl, long nh)
{
	float *v;
	
	v=(float *)malloc((size_t)((nh-nl+1+NR_END)*sizeof(float)));
	if(!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}

int *ivector(long nl, long nh)
{
	int *v;
	
	v=(int *)malloc((size_t)((nh-nl+1+NR_END)*sizeof(int)));
	if(!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}

unsigned char *cvector(long nl, long nh)
{
	unsigned char *v;
	
	v=(unsigned char *)malloc((size_t)((nh-nl+1+NR_END)*sizeof(unsigned char)));
	if(!v) nrerror("allocation failure in cvector()");
	return v-nl+NR_END;
}

unsigned long *lvector(long nl, long nh)
{
	unsigned long *v;
	
	v=(unsigned long *)malloc((size_t)((nh-nl+1+NR_END)*sizeof(unsigned long)));
	if(!v) nrerror("allocation failure in lvector()");
	return v-nl+NR_END;
}

double *dvector(long nl, long nh)
{
	double *v;
	
	v=(double *)malloc((size_t)((nh-nl+1+NR_END)*sizeof(double)));
	if(!v) nrerror("allocation failure in dvector()");
	return v-nl+NR_END;
}

float **matrix(long nrl,long nrh,long ncl,long nch)
{
	long i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;
	
	m=(float **)malloc((size_t)((nrow+NR_END)*sizeof(float*)));
	if(!m) nrerror("allocation failure 1 in matrix()");
	m+=NR_END;
	m-=nrl;
	m[nrl]=(float *)malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
	if(!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl]+=NR_END;
	m[nrl]-=ncl;
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	return(m);
}

double **dmatrix(long nrl,long nrh,long ncl,long nch)
{
	long i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;
	
	m=(double **)malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if(!m) nrerror("allocation failure 1 in dmatrix()");
	m+=NR_END;
	m-=nrl;
	m[nrl]=(double *)malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if(!m[nrl]) nrerror("allocation failure 2 in dmatrix()");
	m[nrl]+=NR_END;
	m[nrl]-=ncl;
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	return(m);
}

int **imatrix(long nrl,long nrh,long ncl,long nch)
{
	long i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;
	
	m=(int **)malloc((size_t)((nrow+NR_END)*sizeof(int*)));
	if(!m) nrerror("allocation failure 1 in imatrix()");
	m+=NR_END;
	m-=nrl;
	m[nrl]=(int *)malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
	if(!m[nrl]) nrerror("allocation failure 2 in imatrix()");
	m[nrl]+=NR_END;
	m[nrl]-=ncl;
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	return(m);
}

float **submatrix(float **a,long oldrl,long oldrh,long oldcl,long oldch,long newrl,long newcl)
{
	long i,j,nrow=oldrh-oldrl+1,ncol=oldch-oldcl+1;
	float **m;
	
	m=(float **)malloc((size_t)((nrow+NR_END)*sizeof(float*)));
	if(!m) nrerror("allocation failure in submatrix()");
	m+=NR_END;
	m-=newrl;
	for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;
	return(m);
}

float **convert_matrix(float *a,long nrl,long nrh,long ncl,long nch)
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;
	
	m=(float **)malloc((size_t)((nrow+NR_END)*sizeof(float*)));
	if(!m) nrerror("allocation failure in convert_matrix()");
	m+=NR_END;
	m-=nrl;
	m[nrl]=a-ncl;
	for(i=1,j=nrl+1;i<=nrow;i++,j++) m[j]=m[i-1]+ncol;
	return(m);
}

float ***f3tensor(long nrl,long nrh,long ncl,long nch,long ndl,long ndh)
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	float ***t;
	
	t=(float ***)malloc((size_t)((nrow+NR_END)*sizeof(float**)));
	if(!t) nrerror("allocation failure 1 in f3tensor()");
	t+=NR_END;
	t-=nrl;
	t[nrl]=(float **)malloc((size_t)((nrow*ncol+NR_END)*sizeof(float*)));
	if(!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
	t[nrl]+=NR_END;
	t[nrl]-=ncl;
	t[nrl][ncl]=(float *)malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(float)));
	if(!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
	t[nrl][ncl]+=NR_END;
	t[nrl][ncl]-=ndl;
	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++)
	{
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol+ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}
	return(t);
}

void free_vector(float *v,long nl,long nh)
{
	free((FREE_ARG)(v+nl-NR_END));
}

void free_ivector(int *v,long nl,long nh)
{
	free((FREE_ARG)(v+nl-NR_END));
}

void free_cvector(unsigned char *v,long nl,long nh)
{
	free((FREE_ARG)(v+nl-NR_END));
}

void free_lvector(unsigned long *v,long nl,long nh)
{
	free((FREE_ARG)(v+nl-NR_END));
}

void free_dvector(double *v,long nl,long nh)
{
	free((FREE_ARG)(v+nl-NR_END));
}

void free_matrix(float **m,long nrl,long nrh,long ncl,long nch)
{
	free((FREE_ARG)(m[nrl]+ncl-NR_END));
	free((FREE_ARG)(m+nrl-NR_END));
}

void free_dmatrix(double **m,long nrl,long nrh,long ncl,long nch)
{
	free((FREE_ARG)(m[nrl]+ncl-NR_END));
	free((FREE_ARG)(m+nrl-NR_END));
}

void free_imatrix(int **m,long nrl,long nrh,long ncl,long nch)
{
	free((FREE_ARG)(m[nrl]+ncl-NR_END));
	free((FREE_ARG)(m+nrl-NR_END));
}

void free_submatrix(float **b,long nrl,long nrh,long ncl,long nch)
{
	free((FREE_ARG)(b+nrl-NR_END));
}

void free_convert_matrix(float **b,long nrl,long nrh,long ncl,long nch)
{
	free((FREE_ARG)(b+nrl-NR_END));
}

void free_f3tensor(float ***t,long nrl,long nrh,long ncl,long nch,long ndl,long ndh)
{
	free((FREE_ARG)(t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG)(t[nrl]+ncl-NR_END));
	free((FREE_ARG)(t+nrl-NR_END));
}


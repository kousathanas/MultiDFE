// generate_allele_distrib_n1_n2.c: generates the allele frequency vector
// from a starting frequency vector that has bee read in after t2 generations
// incorporating a change of population size from n1 to n2.
#include "libhdr"
#include "bctmat01.hd"

#define maxt2 1000

output_startingfreq(int n1, int n2, double selcoeff, double *startingfreq, int t2,
   char *file_name)
{
   FILE *outf, *fopen();
   int i, count;
   outf = openforwritesilent(file_name, "a");
   fprintf(outf, "%d %d %lf\n2", n2, t2, selcoeff);
//   printf("output_startingfreq: %d %d %lf\n", n2, t2, selcoeff); monitorinput();
   count = 0;
   for (i=0; i<=2*n2; i++)
   {
      fprintf(outf, "%10.8lf", startingfreq[i]);
      count++;
      if (count==5)
      {
         count = 0;
         if (i!= 2*n2) fprintf(outf, "\n");
      }      
      else fprintf(outf, " ");
   }
   fprintf(outf, "\n");
   fclose(outf);
}

output_steadystatefreq(int n2, double selcoeff,
   double steadystatefreq[4][maxnd+1], int t2,
   char *file_name)
{
   FILE *outf, *fopen();
   int i, count;
   outf = openforwritesilent(file_name, "a");
   fprintf(outf, "%d %d %lf\n2", n2, t2, selcoeff);
   count = 0;
   for (i=0; i<=2*n2; i++)
   {
      fprintf(outf, "%10.8lf", steadystatefreq[3][i]);
      count++;
      if (count==5)
      {
         count = 0;
         if (i!= 2*n2) fprintf(outf, "\n");
      }      
      else fprintf(outf, " ");
   }
   fprintf(outf, "\n");
   fclose(outf);
}

read_t2_vector(char *t2_file_name, int *t2_vector, int *nt2)
{
   FILE *infile, *fopen();
   int stat, i1;
   infile = openforreadsilent(t2_file_name);
   *nt2 = 0;
   for (;;)
   {
      stat = fscanf(infile, "%d", &i1);
      if (stat==EOF) break;
      if (stat!=1) gabort("read_t2_vector: read error 1", 1);
      t2_vector[*nt2] = i1;
      (*nt2)++;
      if (*nt2 >= maxt2) gabort("read_t2_vector: nt2 >= maxt2", *nt2);
   }
//   printf("read_t2_vector: elements read %d\n", *nt2);
}


int t_in_t2_vector(int t, int *t2_vector, int nt2)
{
   int i;
   for (i=0; i<nt2; i++)
   {
      if (t == t2_vector[i]) return 1;
   }
   return 0;
}


tmiterate1(double a[maxnd+1][maxnd+1], double s, double h, double *mut1,
   double steadystatefreq[4][maxnd+1], char *t2_file_name, int t_max,
   int n1, int n2, int output_file, char *file_name)
{
   int k, i, j, nt2;
   double z, gf, res;
   static double u[maxnd+1], mut2[maxnd+1], temp1[4];
   static int t2_vector[maxt2];
   int n1d, n2d;
   n1d = 2*n1;
   n2d = 2*n2;
// Ensure that outputs are made 1 generation apart to account for the 1
// generation change in population size for the case of phase 1 vec
   if (output_file>0)
   {
      read_t2_vector(t2_file_name, t2_vector, &nt2);
      t_max = t2_vector[nt2-1];
   }
//   dumpvector(mut1, 0, n1d, "tmiterate1: Vector mut1 - starting state");
//   monitorinput();
   for (i=0; i<=n1d; i++)
   {
      u[i]=((double)i/n1d);
   }
   for (i=0; i<=n1d; i++) steadystatefreq[3][i] = 0;
   for (k = 1; k<=t_max; k++)                  /*t iterations*/
   {
      for (i = 0; i<=n1d; i++)
      {
         gf = mut1[i]*u[i];
      }
/* Increment steadystatefreq. distribution*/
      for (i = 0; i<=n1d; i++)
      {
         steadystatefreq[3][i] = steadystatefreq[3][i] + mut1[i];
      }
/* Assign last but 2 and last but 1 generations as well */
      if ((k == t_max - 1) || (k == t_max - 2))
      {
         for (i = 1; i<=n1d - 1; i++)
         {
            steadystatefreq[t_max - k][i] = steadystatefreq[3][i];
         }
      }

/* Perform matrix multiplication*/
      for (i = 0; i<=n2d; i++)
      {
         z = 0;
         for (j = 0; j<=n1d; j++)
         {
            z = z + a[j][i]*mut1[j];
         }
         mut2[i] = z;
      }

/* Copy result of multiplication*/
      for (i=0; i<=n2d; i++) mut1[i] = mut2[i];
//      printf("tmiterate: t %d ", k);
//      dumpvector(mut1, 0, n2d, "Vector mut1");
//      monitorinput();
      if (output_file>0)
      {
//         printf("output starting frequency vec at t %d\n", k);
         if ((output_file==1)&&(t_in_t2_vector(k+1, t2_vector, nt2)))
         {
            output_startingfreq(n1, n2, s, mut1, k+1, file_name);
         }
         else if ((output_file==2)&&(t_in_t2_vector(k, t2_vector, nt2)))
         {
            output_steadystatefreq(n2, s, steadystatefreq, k, file_name);
         }
      }
   }

   for (i=0; i<=n2d; i++)
   {
      for (j=1; j<=3; j++)
      {
         temp1[j] = steadystatefreq[j][i];
//         printf("temp1[%d] %lf\n", j, temp1[j]);
      }
      res = aitken(temp1, 3);
//      printf("res %lf\n", res);
//      monitorinput();
      steadystatefreq[0][i] = res;
   }
//   dumpvector(mut1, 0, n2d,
//      "genefreq. probability vector (mut1) - after t generations"); monitorinput();
}


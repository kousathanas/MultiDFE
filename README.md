#Steps for installing and running MultiDFE

1) Unpack, cd to the unpacked directory and then compile with this command:

gcc -O3 -o MultiDFE *.c  -lm -lgsl -lgslcblas -w

Alternatively, you can also type  
make

2) Now you should have an executable Multi_DFE.

To run the program you need a folder with the lookup tables for the likelihood calculation. The program will not run without specifying the folder where these tables are located. These tables are quite large in size (4.5 Gbytes). 

These lookup tables can be downloaded here :  http://datashare.is.ed.ac.uk/bitstream/handle/10283/2730/data.tar.gz?sequence=1&isAllowed=y

If you have trouble downloading them make sure to contact me (akousath *at* pasteur.fr) to help you.

After you download these tables, you need to specify in the file "directory_config.dat" the location of the tables. The directory specification should branch so the contents of the specified folder are these 3 folders: n1_20, n1_100, n1_1000. In the present implemention only lookup tables in folder n1_100 will be used.

3) If you have specified the lookup tables correctly then you can run the program with the following command:

./MultiDFE -conpop 0 -sfsfold 1 -selmode 2 -nspikes 0 -ranrep 1 -file example.sfs

the example.sfs file has an example of the input SFS (see at the bottom of README for explanation)  

###Explanation of arguments:  
-conpop    0/1    model population size change (2-epoch change/constant)  
-sfsfold   0/1    fold sfs (no/yes)  
-selmode (0/1/2/3/4/5)    selection model (spike/step/gamma/beta/lognormal/6-fixed-spikes)  
-nspikes (0...Inf)    Specify how many spikes you want only if selmode=(0/1)  
-ranrep (0..Inf)    Specify how many replicate runs you want from random starting values (recommended 5-10 for spike/step models).  

#########################################################################################
#########################################################################################
#INPUT files

Please check the example.sfs file which looks like this:  

20  
994838 2351 788 478 271 233 165 139 120 85 69 70 72 61 46 58 31 35 40 50 0  
972726 10843 4558 2650 1700 1231 887 789 705 534 498 442 385 354 323 296 271 294 274 240 0 
  
You have three lines:  
Line 1: Sample size  
Line 2&3: the SFS, space or tab (white space) separated  
Line2: SFS for selected class  
Line3: SFS for neutral class  

the SFSs should always be of length sample size +1  
if you have a sample size of 20, then the SFS contains 21 values (this is because the 0th column, ie. unmutated class is included)  

I recommend that you use always the folded spectrum so have option -sfsfold set to 1 (this will always fold the SFS).  

#########################################################################################  
#########################################################################################  
#OUTPUT files  

After the program finishes running you will get an output file that has the same name as the input but with the addition of .MAXL.out suffix  

so in this case it will be example.sfs.MAXL.out  

The output is appended to this file as lines that look like this:  

seed:0  acc:0   selmode:2       nspikes:0       ranrep:1        L:-217141.5448478720500134      f0:0.929167     N2:279  t:122.716564    Nw:172.576286   E(s):-1255.653717       beta:0.136181   Nes_0.0_0.1:1.114742E-01        Nes_0.1_1.0:4.105598E-02        Nes_1.0_10.0:5.617678E-02     Nes_10.0_100.0:7.686491E-02     mean:-1.192397E+01      mean2:-6.341868E+02     meanH:-7.711070E-03     fix_prob:0.141050  

###explanation of OUTPUT:  
  
seed: the seed used for the random generator. You can change it by adding GSL_RNG_SEED= before the command like this:  
GSL_RNG_SEED=1 ./MultiDFE -conpop 0 -sfsfold 1 -selmode 2 -nspikes 0 -ranrep 1 -file example.out  

L: the log-likelihood of the model given the data  

parameter estimates:  
f0: mutation parameter  
N2,t: demographic parameter  

Nw: the weighted effective pop.size over the 2-epochs  

depending on the selection model you chose, the corresponding parameter estimates will be printed. For example if you chose gamma (-selmode 2) then you will get :  
E(s):-1255.653717       beta:0.136181  

if you chose a 2-spike model (-selmode 0 -nspikes 2),  you will get:  
s1:-6.995761E-04        s2:-7.590491E-01        p1:1.738205E-01 p2:8.261795E-01  

for the sel.coeffs and probabilities for each spike.  

then you also have several statistics printed such as:  
prob. density in 4 Nes ranges:  
Nes_0.0_0.1:1.114742E-01        Nes_0.1_1.0:4.105598E-02        Nes_1.0_10.0:5.617678E-02     Nes_10.0_100.0:7.686491E-02 

please note that Nes >100 is not printed. It is equal to 1-sum(of other range proportions)   

mean, squared mean E(s^2) (mean2) and harmonic mean effect of a new mutation (obtained through integration of the DFE):  
mean:-1.192397E+01      mean2:-6.341868E+02     meanH:-7.711070E-03 

the average fixation probability of a new mutation (this can be used together with divergence to calculate alpha and omega_a as described in the DFE and faster-X genetics papers)
fix_prob:0.141050



#########################################################################################  
#GENERAL ADVICE 

###For getting mean (Nes):

To get an estimate for E(Nes) you should multiply the entries Nw and E(s). It's the standard for comparison with other methods such as DFE_alpha. For the example: Nw:172.576286 E(s):-1255.653717 : Nw E(s)=−216696

Having large values for E(Nes) is not a problem and is usually expected because there is not much power to estimate it for relatively small samples.

###To get alpha and omega_a:
you should use the formulas 10 and 11 from Kousathanas and Keightley (2013). The average fixation probability is given by MultiDFE and to obtain dN and dS you just need to divide the number of nonsynonymous and synonymous substitutions over the number for sites, respectively. A simple Jukes-Cantor correction can be applied with an R-function:

```
####################################   
***div.jukes***
calculates Jukes-Cantor divergence.
Input: x<-total sites,y<-site diffs
####################################   
div.jukes<-function(x,y)
{
d<-vector(length=length(x));
for (i in 1:length(x))
{
if (y[i]<=0){d[i]=NA;next;}
p=y[i]/x[i]
if ((1-(4/3)*p)<0){d[i]=NA;next;}
d[i]=(-3/4)*log(1-(4/3)*p)
}

return(d)
}
```

### Bad performance of the fixed 6-spikes and the beta models:
usually the fixed 6-spikes and beta models have bad performance compared to the other models because they cannot predict mutations with very strong effects. The beta is constrained to s=[0,1], and there is a similar constraint for the 6-fixed spikes model. Of course s cannot be greater than 1, at least theoretically. But in the inference method s can be greater than 1 to be able to use sufficiently large Nes values to be fitted to a realistic dataset.

This can be understood by the fact that MultiDFE uses evolutionary scaling. It models the evolutionary process in smaller population sizes (say N=100) and scales the parameter Nes to fit data that could have been generated for a much larger population size. To do this it has to accept s larger than 1.

### Increasing no.of steps/spikes and log-likelihood values
In respect to the log-likelihood values, as you increase the number of parameters, they are not expected to get worse (only better or at least as good as a nested model with less parameters). Getting a lower log-likelihood when increasing the number of parameters is a sign of non-convergence due to the algorithm getting stuck in a local logL maximum. In such a case you could exclude the models that did not converge, or put a star next to the likelihood of these models and mention under the table "not converged".
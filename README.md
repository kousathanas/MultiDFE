Steps for running MultiDFEest

1) Unpack, cd to the unpacked directory and then compile with this command:

gcc -O3 -o MultiDFE *.c  -lm -lgsl -lgslcblas -w

2) Now you should have an executable Multi_DFE.

To run the program you need a folder with the lookup tables for the likelihood calculation. The program will not run without specifying the folder where these tables are located. These tables are quite large in size (5.3 Gbytes for N1=100). These lookup tables can be downloaded at prof. Peter Keigtley's website: http://www.homepages.ed.ac.uk/pkeightl//dfe_alpha/dfe-alpha-download.html along with the original DFE-alpha program.

After you download these tables, you need to specify in the file "directory_config.dat" the location of the tables. The directory specification should branch so the contents of the specified folder are these 3 folders: n1_20, n1_100, n1_1000

3) If you have specified the lookup tables correctly then you can run the program with the following command:

./Multi_DFE -N1 100 -conpop 0 -sfsfold 1 -selmode 2 -nspikes 0 -ranrep 1 -file example.sfs

the example.sfs file has an example of the input SFS (see at the bottom of README for explanation)

Explanation of arguments:
-conpop    0/1    model population size change (constant/2-epoch change)
-sfsfold    0/1    fold sfs (no/yes)
-selmode (0/1/2/3/4/5)    selection model (spike/step/gamma/beta/lognormal/6-fixed-spikes)
-nspikes (0...Inf)    Specify how many spikes you want only if selmode=(0/1)
-ranrep (0..Inf)    Specify how many replicate runs you want from random starting values (recommended 5-10 for spike/step models).

#########################################################################################
#########################################################################################
#OUTPUT FILE format
After the program finishes running you will get an output file that has the same name as the input but with the addition of .MAXL.out suffix

so in this case it will be example.sfs.MAXL.out

The output is appended to this file as lines that look like this:

seed:0  acc:0   selmode:2       nspikes:0       ranrep:1        L:-217141.5448478720500134      f0:0.929167     N2:279  t:122.716564    Nw:172.576286   E(s):-1255.653717       beta:0.136181   Nes_0.0_0.1:1.114742E-01        Nes_0.1_1.0:4.105598E-02        Nes_1.0_10.0:5.617678E-02     Nes_10.0_100.0:7.686491E-02     mean:-1.192397E+01      mean2:-6.341868E+02     meanH:-7.711070E-03     fix_prob:0.141050

explanation of OUTPUT:
seed: the seed used for the random generator. You can change it by adding GSL_RNG_SEED= before the command.
GSL_RNG_SEED=1 ./Multi_DFE -N1 100 -conpop 0 -sfsfold 1 -selmode 2 -nspikes 0 -ranrep 1 -file example.ou

L: the loglikelihood of the model given the data

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

mean, squared mean E(s^2) (mean2) and harmonic mean effect of a new mutation (obtained through integration of the DFE):
mean:-1.192397E+01      mean2:-6.341868E+02     meanH:-7.711070E-03 

the average fixation probability of a new mutation (this can be used together with divergence to calculate alpha and omega_a as described in the DFE and faster-X genetics papers)
fix_prob:0.141050


#########################################################################################
#########################################################################################
#INPUT FILE format

Please check the example.sfs file which looks like this:

name
0 0
0 0
20
994838 2351 788 478 271 233 165 139 120 85 69 70 72 61 46 58 31 35 40 50 0
972726 10843 4558 2650 1700 1231 887 789 705 534 498 442 385 354 323 296 271 294 274 240 0

You have six lines:
Line 1: Just put a string for name, it doesnt matter
Line 2&3: this was supposed to contain the divergence data. Since my program only infers the DFE, these are not even used so just leave it at 0 for now.
Line 4: Sample size
Line 5&6: the SFS, space or tab (white space) separated
Line5: SFS for selected class
Line6: SFS for neutral class

the SFSs should always be of length sample size +1
if your you have a sample size of 20, then the SFS contains 21 values (this is because the 0th column, ie. unmutated class is included)
Basically you only need to vary lines 4,5,6 between datasets. Lines 1,2,3 dont matter.

I recommend that you use always the folded spectrum so have option -sfsfold set to 1 (this will always fold the SFS).
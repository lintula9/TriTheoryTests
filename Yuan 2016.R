#The code is to evaluate the cutoff values RMSEA_e in equivalence testing
#corresponding to the conventional cutoff values of RMSEA=.01, .05, .08, and .10, respectively;
#as described in the article
#"Confirm Structural Equation Models by Equivalence Testing with Adjusted Fit Indices"
#by Yuan, Chan, Marcoulides and Bentler;

#needed inputs are degrees of freedom (df) and sample size (N);

df=23; N=145; n=N-1;

RMSEA_e01=exp(
  1.34863-.51999*log(df)+.01925*log(df)*log(df)-.59811*log(n)+.00902*sqrt(n)+.01796*log(df)*log(n)
);
#corresponding to R-square=.9997;

RMSEA_e05=exp(
  2.06034-.62974*log(df)+.02512*log(df)*log(df)-.98388*log(n)
  +.05442*log(n)*log(n)-.00005188*n+.05260*log(df)*log(n)
);
#corresponding to R-square=.9996;

RMSEA_e08=exp(
  2.84129-.54809*log(df)+.02296*log(df)*log(df)-.76005*log(n)
  +.10229*log(n)*log(n)-1.11167*(n^.2)+.04845*log(df)*log(n)
);
#corresponding to R-square=.9977;

RMSEA_e10=exp(
  2.36352-.49440*log(df)+.02131*log(df)*log(df)-.64445*log(n)
  +.09043*log(n)*log(n)-1.01634*(n^.2)+.04422*log(df)*log(n)
);
#corresponding to R-square=.9955;

cutoff=cbind(RMSEA_e01, RMSEA_e05, RMSEA_e08, RMSEA_e10);
cutoff_3=round(cutoff,3);
print(cutoff);

cat('--excellent--', cutoff_3[1], '--close--', cutoff_3[2], '--fair--', cutoff_3[3], '--mediocre--', cutoff_3[4], '--poor--',"\n")
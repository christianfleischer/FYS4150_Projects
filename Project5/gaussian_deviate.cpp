// function for gaussian random numbers
double gaussian_deviate(long *);
// ran2 for uniform deviates, initialize with negative seed.

// random numbers with gaussian distribution
double gaussian_deviate(long * idum)
{
static int iset = 0;
static double gset;
double fac, rsq, v1, v2;
if ( idum < 0) iset =0;
if (iset == 0) {
do {
v1 = 2.*ran2(idum) -1.0;
v2 = 2.*ran2(idum) -1.0;
rsq = v1*v1+v2*v2;

} while (rsq >= 1.0 || rsq == 0.);
fac = sqrt(-2.*log(rsq)/rsq);
gset = v1*fac;
iset = 1;
return v2*fac;
} else {
iset =0;
return gset;
}
} // end function for gaussian deviates

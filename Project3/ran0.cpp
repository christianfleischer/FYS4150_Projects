
/*
** The function
**           ran0()
** is an "Minimal" random number generator of Park and Miller
** (see Numerical recipe page 279). Set or reset the input value
** idum to any integer value (except the unlikely value MASK)
** to initialize the sequence; idum must not be altered between
** calls for sucessive deviates in a sequence.
** The function returns a uniform deviate between 0.0 and 1.0.
*/

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

double ran0(long *idum)
{
long     k;
double   ans;

*idum ^= MASK;
k = (*idum)/IQ;
*idum = IA*(*idum - k*IQ) - IR*k;
if(*idum < 0) *idum += IM;
ans=AM*(*idum);
*idum ^= MASK;
return ans;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef MASK

// End: function ran0()

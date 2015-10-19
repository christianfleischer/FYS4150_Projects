#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <stdio.h>
#include <ctime>
#include <gauleg.cpp>
#include <gauss_laguerre.cpp>
#include <ran0.cpp>
//#include <omp.h>
//#define NUM_THREADS 4

using namespace std;

double int_function_a(double x1, double y1, double z1, double x2, double y2, double z2);
double int_function_b(double r1, double r2, double theta1, double theta2, double phi1, double phi2);
double brute_force_MC(double *);
double exponential_MC(double *);

int main()
{
    // Integration limits a and b:
         double a = -2;     //-2
         double b = 2;      //2
         int Nleg = 51;     //51 for Legendre to 3 digits (takes a long time).
         int Nlag = 16;     //16 for Laguerre to 3 digits.
         double const  pi = 3.14159265359;
    // The closed form solution which is known for the integral we want to solve:
         double ClosedForm = 5*pi*pi/(16*16);
    // Arrays for integration points and weights using Legendre polynomials:
         double *x = new double [Nleg];
         double *w = new double [Nleg];
    // Set up the mesh points and weights:
         gauleg(a,b,x,w, Nleg);

    // Gauss-Laguerre is old-fashioned translation of F77 --> C++
    // arrays start at 1 and end at n.
    // Arrays for integration points and weights using Laguerre polynomials:
         double *xGLaguerre = new double [Nlag+1];
         double *wGLaguerre = new double [Nlag+1];
    // Set up the mesh points and weights and the power of x^alf:
         double alf = 2.0;
         gauss_laguerre(xGLaguerre,wGLaguerre, Nlag, alf);
    // rrays for integration points and weights using Legendre polynomials for the angular variables:
         double *xt = new double [Nlag];
         double *wt = new double [Nlag];
         double *xp = new double [Nlag];
         double *wp = new double [Nlag];
    // Set up the mesh points and weights:
         gauleg(0, pi, xt, wt, Nlag);
         gauleg(0, 2*pi, xp, wp, Nlag);
         
         //I attempted to parallelize the Gauss-Legendre method, but didn't manage to get it to work. 
         //omp_set_num_threads(NUM_THREADS);

         //cout << "  The number of threads available    = " << omp_get_max_threads ( )  <<  endl;
         //cout << "  The number of processors available = " << omp_get_num_procs ( ) << endl;

    // Start timing Gauss-Legendre:
         clock_t start1, finish1;
         start1 = clock();
    // Evaluate the integral with the Gauss-Legendre method:
    // Note that we initialize the sum.
         double int_gauss = 0.;
         double percent = 0.; //See how close the Gauss-Legendre method is to finshing (useful for those 15 min calculations).
    // Six-double loops:
         //# pragma omp parallel for reduction(+:int_gauss) private (i, x, w)
         for (int i=0;i<Nleg;i++){
                 cout << percent << "% " << endl;
                 percent += 100./Nleg;
                 for (int j = 0;j<Nleg;j++){
                 for (int k = 0;k<Nleg;k++){
                 for (int l = 0;l<Nleg;l++){
                 for (int m = 0;m<Nleg;m++){
                 for (int n = 0;n<Nleg;n++){
            int_gauss+=w[i]*w[j]*w[k]*w[l]*w[m]*w[n]
           *int_function_a(x[i],x[j],x[k],x[l],x[m],x[n]);
                    }}}}}
            }
         cout << percent << "%" << endl;
    // Finish timing Gauss-Legendre:
         finish1 = clock();
         double ComputationTimeLeg = ((finish1-start1)/(double) CLOCKS_PER_SEC);
         cout << endl
              << "N: " << Nleg << endl
              << "Gauss-Legendre: " << int_gauss << endl
              << "Closed Form: " << ClosedForm << endl
              << "Computation time: " << ComputationTimeLeg << endl;

    // Start timing Gauss-Laguerre
         clock_t start2, finish2;
         start2 = clock();
    // Evaluate the integral with the Gauss-Laguerre method:
    // Note that we initialize the sum.
         double int_gausslag = 0.;
         for ( int i = 1;  i <= Nlag; i++){
                 for (int j = 1;j<=Nlag;j++){
                 for (int k = 0;k<Nlag;k++){
                 for (int l = 0;l<Nlag;l++){
                 for (int m = 0;m<Nlag;m++){
                 for (int n = 0;n<Nlag;n++){
            int_gausslag += wGLaguerre[i]*wGLaguerre[j]*wt[k]*wt[l]*wp[m]*wp[n]
           *sin(xt[k])*sin(xt[l])*int_function_b(xGLaguerre[i],xGLaguerre[j],xt[k],xt[l],xp[m],xp[n]);
                    }}}}}
            }

         int_gausslag /= pow(4,5);  //Need to divide by this factor due to the variable change r = u/4.
    // Finish timing Gauss-Laguerre:
         finish2 = clock();
         double ComputationTimeLag = ((finish2-start2)/(double) CLOCKS_PER_SEC);
         cout << endl
              << "N: " << Nlag << endl
              << "Gauss-Laguerre: " << int_gausslag << endl
              << "Closed Form: " << ClosedForm << endl
              << "Computation time: " << ComputationTimeLag << endl;

    // Test that we get (to computer precision) exactly pi^2/64 without the 1/sqrt factor gvien in int_function_b:
         double int_test = 0.;
         for ( int i = 1;  i <= Nlag; i++){
                 for (int j = 1;j<=Nlag;j++){
                 for (int k = 0;k<Nlag;k++){
                 for (int l = 0;l<Nlag;l++){
                 for (int m = 0;m<Nlag;m++){
                 for (int n = 0;n<Nlag;n++){
            int_test += wGLaguerre[i]*wGLaguerre[j]*wt[k]*wt[l]*wp[m]*wp[n]
           *sin(xt[k])*sin(xt[l]);
                 }}}}}
            }

         int_test /= pow(4,6);
         double tolerance = pow(10, -10);
         if (abs(int_test - pi*pi/64) < tolerance){
             cout << "TEST SUCCESSFUL" << endl;
         }

         int nBrute = 2*pow(10, 7);   //Number of simulations for brute force MC.
         double xMCbrute[6], fxBrute;
         double int_mcBrute = 0.; double varianceBrute = 0.;
         double sum_sigmaBrute = 0. ; long idum = -1 ;
         double length = 2.; //We fix the max size of the box to L=2.
         double jacobidetBrute = pow((2*length),6);

    // Start timing brute force Monte Carlo:
         clock_t start3, finish3;
         start3 = clock();
         for (int i = 1; i <= nBrute; i++){
             for (int j = 0; j < 6; j++){
                 xMCbrute[j] = -length + 2*length*ran0(&idum);
             }
             fxBrute = brute_force_MC(xMCbrute);
             int_mcBrute += fxBrute;
             sum_sigmaBrute += fxBrute*fxBrute;
         }
         int_mcBrute = int_mcBrute/((double) nBrute );
         sum_sigmaBrute = sum_sigmaBrute/((double) nBrute );
         varianceBrute=sum_sigmaBrute - int_mcBrute*int_mcBrute;
    // Finish timing brute force Monte Carlo:
         finish3 = clock();
         double ComputationTimeMCbrute = ((finish3-start3)/(double) CLOCKS_PER_SEC);

         cout << endl
              << "N: " << nBrute << endl
              << "Montecarlo (BF): " << jacobidetBrute*int_mcBrute << endl
              << "Sigma: " << jacobidetBrute*sqrt(varianceBrute/((double) nBrute)) << endl
              << "Computation time: " << ComputationTimeMCbrute << endl;

         int nImportance = 7*pow(10,4);   //Number of simulations for importance sampling MC.
         double xMCimportance[6], fxImportance;
         double int_mcImportance = 0.; double varianceImportance = 0.;
         double sum_sigmaImportance = 0; long idum2 = -1;
         double jacobidetImportance = (pow(pi,4)/4);

    // Start timing importance sampling Monte Carlo:
         clock_t start4, finish4;
         start4 = clock();
         for (int i = 1; i <= nImportance; i++){
             xMCimportance[0] = (-1./4)*log(1-ran0(&idum2));
             xMCimportance[1] = (-1./4)*log(1-ran0(&idum2));
             xMCimportance[2] = pi*ran0(&idum2);
             xMCimportance[3] = pi*ran0(&idum2);
             xMCimportance[4] = 2*pi*ran0(&idum2);
             xMCimportance[5] = 2*pi*ran0(&idum2);
             fxImportance = exponential_MC(xMCimportance);
             int_mcImportance += fxImportance;
             sum_sigmaImportance += fxImportance*fxImportance;
         }
         int_mcImportance = int_mcImportance/((double) nImportance);
         sum_sigmaImportance = sum_sigmaImportance/((double) nImportance);
         varianceImportance = sum_sigmaImportance - int_mcImportance*int_mcImportance;
    // Finish timing importance sampling Monte Carlo:
         finish4 = clock();
         double ComputationTimeMCimportance = ((finish4-start4)/(double) CLOCKS_PER_SEC);

         cout << endl
              << "N: " << nImportance << endl
              << "Montecarlo (IS): " << jacobidetImportance*int_mcImportance << endl
              << "Sigma: " << jacobidetImportance*sqrt(varianceImportance/((double) nImportance)) << endl
              << "Computation time: " << ComputationTimeMCimportance << endl;
}
// This function defines the function to integrate for Gauss-Legendre:
double int_function_a(double x1, double y1, double z1, double x2, double y2, double z2)
{
   double alpha = 2.;
// Evaluate the different terms of the exponential:
   double exp1=-2*alpha*sqrt(x1*x1+y1*y1+z1*z1);
   double exp2=-2*alpha*sqrt(x2*x2+y2*y2+z2*z2);
// Evaluate the denominator:
   double deno=sqrt(pow((x1-x2),2)+pow((y1-y2),2)+pow((z1-z2),2));
// Cheating!!
  if(deno <pow(10.,-6.)) { return 0;}
  else return exp(exp1+exp2)/deno;
} // End of function to evaluate.

// This function defines the denominator in the function to integrate for Gauss-Laguerre:
double int_function_b(double r1, double r2, double theta1, double theta2, double phi1, double phi2)
{
    double cosbeta = cos(theta1)*cos(theta2) + sin(theta1)*sin(theta2)*cos(phi1 - phi2);
// Evaluate the denominator:    
    double denom = sqrt(r1*r1 + r2*r2 - 2*r1*r2*cosbeta);
    if(denom <pow(10.,-6)||isnan(denom)) { return 0;}
    else return 1./denom;
} // End of function the denominator.

// This function defines the function to integrate for brute force Monte Carlo:
double brute_force_MC(double *x)
{    
    double alpha = 2.;
    double x1 = x[0]; double y1 = x[1]; double z1 = x[2];
    double x2 = x[3]; double y2 = x[4]; double z2 = x[5];
// Evaluate the different terms of the exponential:    
    double exp1 = -2*alpha*sqrt(x1*x1+y1*y1+z1*z1);
    double exp2=-2*alpha*sqrt(x2*x2+y2*y2+z2*z2);
// Evaluate the denominator:    
    double deno=sqrt(pow((x1-x2),2)+pow((y1-y2),2)+pow((z1-z2),2));
    return exp(exp1+exp2)/deno;
} // End of function for the integrand.

// This function defines the function to integrate for importance sampling Monte Carlo:
double exponential_MC(double *x)
{
    double r1 = x[0]; double r2 = x[1]; double theta1 = x[2];
    double theta2 = x[3]; double phi1 = x[4]; double phi2 = x[5];
    double cosbeta = cos(theta1)*cos(theta2) + sin(theta1)*sin(theta2)*cos(phi1 - phi2);
// Evaluate the denominator:    
    double denom = sqrt(r1*r1 + r2*r2 - 2*r1*r2*cosbeta);
    if(denom <pow(10.,-6)||isnan(denom)) { return 0;}
    else return r1*r1*r2*r2*sin(theta1)*sin(theta2)/denom;
} // End of function for the integrand.

#include <iostream>
#include <cmath>
#include <armadillo>

using namespace std;
using namespace arma;


void InitializeOneElectron (int N, mat &A, double rhoMin, double rhoMax, mat &SaveEigenvector)
{
    double h = (rhoMax-rhoMin)/N;

    //Set up the vector rho and the matrix A:
    vec rho = rhoMin + linspace(0, N, N+1)*h;
    vec V = rho%rho;

    SaveEigenvector.col(0) = rho.subvec(1, N-1);    //Saves the rho vector for output.

    double Constant = 1/(h*h);
    A.diag(0)  =  2*Constant + V.subvec(1,N-1);     //Set d_i elements in A
    A.diag(1)  = -Constant*ones(N-2);               //Set e_i elements in A
    A.diag(-1) = A.diag(1);                         //Set e_i elements in A

    return;
}

void InitializeTwoElectrons(int N, mat &A, double rhoMin, double rhoMax, mat &SaveEigenvector, double omega_r)
{
    double h = (rhoMax-rhoMin)/N;

    //Set up the vector rho and the matrix A:
    vec rho = rhoMin + linspace(0, N, N+1)*h;

    //No repulsion:
    //vec V = omega_r*omega_r*(rho%rho);

    //With Repulsion:
    vec V = omega_r*omega_r*(rho%rho) + 1.0/rho;

    SaveEigenvector.col(0) = rho.subvec(1, N-1);    //Saves the rho vector for output.

    double Constant = 1/(h*h);
    A.diag(0)  =  2*Constant + V.subvec(1,N-1);     //Set d_i elements in A
    A.diag(1)  = -Constant*ones(N-2);               //Set e_i elements in A
    A.diag(-1) = A.diag(1);                         //Set e_i elements in A

    return;
}

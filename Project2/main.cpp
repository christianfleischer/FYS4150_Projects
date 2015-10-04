#include <iostream>
#include <armadillo>
#include <time.h>
#include <initialize_electrons.cpp>
#include <jacobi_rotations.cpp>

using namespace std;
using namespace arma;

int main()
{
    int N = 200;
    double rhoMin = 0;
    double rhoMax = 5;
    double omega_r = 0.01;

    mat A = zeros(N-1,N-1);
    mat SaveEigenvector = zeros(N-1, 2);

    //Set initial condtions and set up matrix:
    InitializeOneElectron(N, A, rhoMin, rhoMax, SaveEigenvector);
    //InitializeTwoElectrons(N, A, rhoMin, rhoMax, SaveEigenvector, omega_r);

    clock_t start1, finish1;
    start1 = clock();

    //Finding eigenvalues and eigenvectors using armadillo:
    vec ArmadilloEigenvalues;
    mat Eigenvectors;
    eig_sym(ArmadilloEigenvalues, Eigenvectors, A);

    finish1 = clock();
    double ComputationTimeArma = ((finish1-start1)/(double) CLOCKS_PER_SEC);

    //Finding eigenvalues using Jacobi rotations:
    int NumberOfRotations = 0;      //Counter for similarity transformations.
    double Tolerance = 1E-8;        //Tolerance for when the matrix should be considered as diagonal.
    mat absA = abs(A);              //Creating a matrix for finding the max value element in A.
    absA.diag(0) = zeros(N-1);      //Remove the diagonal elements since we're interested in the max off-diagonal element.

    uword RowIndexMax;      //Row index for the maximum element.
    uword ColIndexMax;      //Column index for the maximum element.
    double MaxValue = absA.max(RowIndexMax, ColIndexMax);   //Value of the max element.

    clock_t start2, finish2;
    start2 = clock();

    //Diagonalize A with Jacobi rotations:
    JacobiRotations(MaxValue, A, absA, NumberOfRotations, Tolerance, RowIndexMax, ColIndexMax, N);
    //Eigenvalues listed on the diagonal of A sorted from smallest to largest:
    vec JacobiEigenvalues = sort(A.diag(0));

    finish2 = clock();
    double ComputationTimeJacobi = ((finish2-start2)/(double) CLOCKS_PER_SEC);

    SaveEigenvector.col(1) = Eigenvectors.col(0);
    //SaveEigenvector.save("omega0_01withrepulsion.dat", raw_ascii);

    cout << "Eigenvalues, Jacobi:" << endl;
    cout << "1:     " << JacobiEigenvalues(0) << endl;
    cout << "2:     " << JacobiEigenvalues(1) << endl;
    cout << "3:     " << JacobiEigenvalues(2) << endl;
    cout << "Eigenvalues, Armadillo:" << endl;
    cout << "1:     " << ArmadilloEigenvalues(0) << endl;
    cout << "2:     " << ArmadilloEigenvalues(1) << endl;
    cout << "3:     " << ArmadilloEigenvalues(2) << endl;
    cout << endl;
    cout << "Number of similarity transformations: " << NumberOfRotations << endl;
    cout << "Computation time (sec):" << endl;
    cout << "Jacobi: " << ComputationTimeJacobi << "    " << "Aramadillo: " << ComputationTimeArma << endl;

    return 0;
}

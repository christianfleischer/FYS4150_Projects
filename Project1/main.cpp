#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <time.h>
#include <armadillo>

using namespace arma;
using namespace std;

int main()
{
    //Creating the grid points vector x:
    int n = 10;
    double h = 1./(n+1);
    vec x = linspace(h, 1-h, n);

    //Creating a vector for the right hand side of the equation Av = ~b:
    double h2 = h*h;
    vec f = 100.*exp(-10.*x);     //Source term vector.
    vec b = zeros(n);
    b.fill(h2);
    b = b%f;
    vec original_b = b;     //Copy of the original RHS vector for use with LU decomposition later.

    //Three vectors corresponding to the non-zero diagonals in the tridiagonal matrix:
    vec Aa = -1.*ones(n);
    vec Ab = 2.*ones(n);
    vec Ac = Aa;

    //Start timer for the tridiagonal method:
    clock_t startTri, finishTri;
    startTri = clock();
    //Boundary conditions:
    Aa(0) = 0;
    Ac(n-1) = 0;

    //Forward substitution:
    Ac(0) = Ac(0)/Ab(0);
    b(0) = b(0)/Ab(0);
    int i;

    for (i = 1; i<=n-1; i++){
        Ac(i) = Ac(i) / (Ab(i) - Aa(i)*Ac(i-1));
        b(i) = (b(i) - Aa(i)*b(i-1)) / (Ab(i) - Aa(i)*Ac(i-1));
    }

    //Backward substitution:
    vec v = zeros(n);   //Vector for the solution.
    v(n-1) = b(n-1);

    for (i=n-1; i>=1; i--) {
        v(i-1) = b(i-1) - Ac(i-1)*v(i);
    }

    //Stop timer for the tridiagonal method and calculate computation time:
    finishTri = clock();
    double ComputationTimeTri = ((finishTri - startTri)/(double) CLOCKS_PER_SEC);

    //ofstream ofile;
    //char *outfilename;
    //outfilename = "tridiag_solution_n10.dat";
    //ofile.open(outfilename);
    //ofile << setiosflags(ios::showpoint | ios::uppercase);
    // //ofile << setprecision(8) << "x:" << setw(21) << "v(x):" << endl;
    //ofile << scientific << setprecision(8) << 0. << setw(18) << 0. << endl;
    //for (i=1; i <= n-1; i++) {
    //    ofile << scientific << setprecision(8) << x(i) << setw(18) << v(i) << endl;
    //}
    //ofile << scientific << setprecision(8) << 1. << setw(18) << 0. << endl;
    //ofile.close();

    //Save the solution vector v(x) and the grid points vector x as a data file,
    //which is used for plotting with the Python program (plot_data.py):
    mat data = zeros(n+2, 2);
    for (i = 1; i<=n; i++) {
        data(i,0) = x(i-1);
        data(i,1) = v(i-1);
    }
    data(n+1,0) = 1;
    data.save("tridiag_solution_n10.dat", raw_ascii);

    //Calculating the relative error of numerical solution (tridiagonal method):
    vec u = 1. - (1-exp(-10))*x - exp(-10*x);    //Exact solution.
    vec rel_error = log10(abs((v-u)/u));
    rel_error.save("rel_error_n10.dat", raw_ascii);
    cout << "For n=" << n << " we have a max value of the relative error: rel_error(max)=" << max(rel_error) << endl;

    //Solving the same problem with LU decomposition to compare the calculation time with that of the tridiagonal method:
    //Creating the tridiagonal matrix:
    mat A = 2.*eye(n, n);        //Creates a nxn matrix elements on the main diagonal equal to 2, all other elemets are 0.
    A.diag(-1) = -1.*ones(n-1);    //Makes all elements on the diagonal below the main diagonal equal to -1.
    A.diag(1) = A.diag(-1);     //Same as above for the diagonal above the main diagonal.
    mat L, U, P;                //Declaring vectors to be used in armadillo's LU decomposition.

    //Start timer for the LU decomposition method:
    clock_t startLU, finishLU;
    startLU = clock();

    //Calculating the solution vector with armadillo LU decomposition:
    lu(L, U, P, A);
    vec y = solve(L, original_b);
    vec vLU = solve(U, y);          //Solution vector obtained with LU.

    //Stop timer for the LU decomposition method and calculate computation time:
    finishLU = clock();
    double ComputationTimeLU = ((finishLU - startLU)/(double) CLOCKS_PER_SEC);

    //Save data from LU method to file:
    mat dataLU = zeros(n+2, 2);
    for (i = 1; i<=n; i++) {
        dataLU(i,0) = x(i-1);
        dataLU(i,1) = vLU(i-1);
    }
    dataLU(n+1,0) = 1;
    dataLU.save("LU_solution_n10.dat", raw_ascii);

    //Compare computation times in terminal:
    cout << "Computation time for n=" << n << endl;
    cout << scientific << "Tridiagonal Solver: " << ComputationTimeTri << "s" << endl;
    cout << scientific << "LU with armadillo: " << ComputationTimeLU << "s" << endl;

    return 0;
}


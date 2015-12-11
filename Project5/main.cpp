#include <iostream>
#include <armadillo>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include "ran2.cpp"
#include "gaussian_deviate.cpp"

using namespace std;
using namespace arma;

void ForwardEuler(int, double, double, double, vec&, vec&, vec&);
void BackwardEuler(int, double, double, double, vec&, vec&, vec&);
void CrankNicolson(int, double, double, double, vec&, vec&, vec&);
vec TridiagSolver(int, vec a, vec b, vec c, vec v);
void MonteCarlo(double D, double d);
void MonteCarloSampling(int, double, int, double, double, vec&);
void SaveToFile(char*, vec, vec, double, int);

int main()
{
    double D = 1;   //Diffusion coefficient
    double d = 1;
    double T = 1;

    int Nx = 100;
    int Nt = 50000;

    double dt = 1./Nt;
    double dx = 1./Nx;
    //double dt = 0.5*dx*dx;  //Maximum stable dt for given Nx.
    //double dt = 1.*dx*dx;  //-> alpha = 1 to test that forward Euler is unstable for this alpha.

    vec x = linspace(0, d, Nx+1);

    double alpha = D*dt/(dx*dx);

    vec vFE = zeros(Nx+1);
    vec vNewFE = zeros(Nx+1);

    //Set boundary conditions
    vFE(0) = 1;
    vNewFE(0) = 1;

    vec vBE = vFE;
    vec vNewBE = vNewFE;
    vec vCN = vFE;
    vec vNewCN = vNewFE;

    ForwardEuler(Nx, T, dt, alpha, x, vFE, vNewFE);
    BackwardEuler(Nx, T, dt, alpha, x, vBE, vNewBE);
    CrankNicolson(Nx, T, dt, alpha, x, vCN, vNewCN);
    MonteCarlo(D, d);

    //cout << alpha << endl;
    //cout << dt << endl;

    return 0;
}

void ForwardEuler(int Nx, double T, double dt, double alpha, vec& x, vec& v, vec& vNew)
{
    //Set initial condition (v(0) = 1, v(Nx) = 0)
    //for (int i=1; i<Nx; i++) v(i) = x(i) - 1;

    if (alpha>0.5) cout << "WARNING: Unstable alpha" << endl;

    double t = dt;
    do{
        for (int i=1; i<Nx; i++){
            vNew(i) = v(i) + alpha*(v(i+1) - 2*v(i) + v(i-1));
        }
        swap(v, vNew);  //Pointer swap.
        t = t + dt;

        if (T/10.-dt/2. < t && t < T/10.+dt/2.){
            char Filename[100];
            sprintf(Filename, "FE_dx%.3f_t%.3f.txt", 1./Nx, t);
            SaveToFile(Filename, v, x, t, Nx);
        }
    }while(t<=T);

    char FilenameSteady[100];
    sprintf(FilenameSteady, "FE_dx%.3f_t%.3f.txt", 1./Nx, T);
    SaveToFile(FilenameSteady, v, x, T, Nx);

    for (int i=0; i<=Nx; i++){
        cout << v(i) << "     " << 1-x(i) << endl;
    }

}

void BackwardEuler(int Nx, double T, double dt, double alpha, vec& x, vec& v, vec& vNew)
{
    vec Aa = -alpha*ones(Nx+1);     //Off diagonal elements of the tridiagonal matrix.
    vec Ab = (1+2*alpha)*ones(Nx+1);    //Diagonal elements of the tridiagonal matrix.
    vec Ac = Aa;

    //Boundary conditions:
    Aa(0) = 0;
    Aa(Nx) = 0;
    Ac(0) = 0;
    Ac(Nx) = 0;

    double t = dt;
    do{
        vNew = TridiagSolver(Nx, Aa, Ab, Ac, v);
        swap(v, vNew);  //Pointer swap.
        t = t + dt;

        if (T/10.-dt/2. < t && t < T/10.+dt/2.){
            char Filename[100];
            sprintf(Filename, "BE_dx%.3f_t%.3f.txt", 1./Nx, t);
            SaveToFile(Filename, v, x, t, Nx);
        }
    }while(t <= T);

    char FilenameSteady[100];
    sprintf(FilenameSteady, "BE_dx%.3f_t%.3f.txt", 1./Nx, T);
    SaveToFile(FilenameSteady, v, x, T, Nx);
    cout << v << endl;
}

void CrankNicolson(int Nx, double T, double dt, double alpha, vec& x, vec& v, vec& vNew)
{
    vec Aa = -alpha*ones(Nx+1);     //Off diagonal elements of the tridiagonal matrix.
    vec Ab = (2+2*alpha)*ones(Nx+1);    //Diagonal elements of the tridiagonal matrix.
    vec Ac = Aa;
    vec vTilde = v;

    //Boundary conditions:
    Aa(0) = 0;
    Aa(Nx) = 0;
    Ac(0) = 0;
    Ac(Nx) = 0;

    //The solver for Crank-Nicolson has two steps.
    double t = dt;
    do{
        //First step: Calculate vTilde = (2*I-alpha*B)*v in the same way as for forward Euler:
        for (int i=1; i<Nx; i++){
            vTilde(i) = alpha*v(i-1) + (2 - 2*alpha)*v(i) + alpha*v(i+1);
        }

        //Second step: Solve the matrix equation (2*I + alpha*B)*vNew = vTilde in the same way as for backward Euler:
        vNew = TridiagSolver(Nx, Aa, Ab, Ac, vTilde);

        swap(v, vNew);  //Pointer swap.
        t = t + dt;

        if (T/10.-dt/2. < t && t < T/10.+dt/2.){
            char Filename[100];
            sprintf(Filename, "CN_dx%.3f_t%.3f.txt", 1./Nx, t);
            SaveToFile(Filename, v, x, t, Nx);
        }
    }while(t <= T);

    char FilenameSteady[100];
    sprintf(FilenameSteady, "CN_dx%.3f_t%.3f.txt", 1./Nx, T);
    SaveToFile(FilenameSteady, v, x, T, Nx);
    cout << v << endl;


}

vec TridiagSolver(int Nx, vec Aa, vec Ab, vec Ac, vec v){
    //Forward substitution
    vec temp(Nx+1);
    double Abtemp = Ab(1);

    vec vNew = v;
    vNew(1) = v(1)/Abtemp;

    for(int i=1; i < Nx ; i++) {
        temp(i) = Ac(i-1)/Abtemp;
        Abtemp = Ab(i)-Aa(i)*temp(i);
        vNew(i) = (v(i) - Aa(i)*vNew(i-1))/Abtemp;
    }

    //Backward substitution
    for(int i=Nx-1 ; i >= 1 ; i--) {
        vNew(i) -= temp(i+1)*vNew(i+1);
    }
    return vNew;
}

void MonteCarlo(double D, double d){
    double NumWalkers = 100000;  //Total amount of walkers.
    double prob = 0.5;
    double T = 0.01;
    double dt = 0.00001;

    bool Gauss = true;//false;
    double l0;
    double xi = 1;

    if (Gauss == true){
        long idumG = -1;
        xi = abs(gaussian_deviate(&idumG));
    }
    l0 = sqrt(2*D*dt)*xi;
    //cout << xi << endl;
    int Nx = round(d/l0)+1;

    vec x = linspace(0, d, Nx+1);
    vec CountPosition = zeros(Nx+1);

    MonteCarloSampling(NumWalkers, T, Nx, dt, prob, CountPosition);

    char Filename[100];
    if (Gauss == true){
        sprintf(Filename, "MC_Gauss_t%.3f.txt", T);
    }
    else{
        sprintf(Filename, "MC_t%.3f.txt", T);
    }

    ofstream outfile(Filename);
    for (int i=0; i<=Nx; i++){
        outfile << CountPosition(i)/NumWalkers << "  " << x(i) << endl;
    }
    outfile.close();
}

void MonteCarloSampling(int NumWalkers, double T, int Nx, double dt, double prob, vec& CountPosition){
    long idum = -1;
    CountPosition(0) = NumWalkers;

    for (int NT=1; NT<=NumWalkers; NT++){
        int Position = 0;

        for (double t=dt; t<=T; t+=dt){

            if (ran2(&idum) <= prob){
                Position -= 1;
                if (Position <= 0){
                    break;
                }
            }
            else{
                Position += 1;
                if (Position>=Nx){
                    break;
                }
            }
            CountPosition(Position) += 1;
        }
    }
}

void SaveToFile(char* Filename, vec v, vec x, double t, int Nx){
    ofstream outfile(Filename);
    for (int i=0; i<=Nx; i++){
        outfile << v(i) << "  " << x(i) << "  " << t << endl;
    }
    outfile.close();
}

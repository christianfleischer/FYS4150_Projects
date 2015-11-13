#include <iostream>
#include <cmath>
#include <cstring>
#include <random>
#include <fstream>
#include <mpi.h>
#include <chrono>

using namespace std;

inline int periodic(int i, int limit, int add) {
    return (i+limit+add) % (limit);
}
void initialize(int NumSpins, int **SpinMatrix, double& E, double& M);
void initializeRandom(int NumSpins, int **SpinMatrix, double& E, double& M);
void closedForm (int NumSpins, double T, double& Z, double& expectE,
                 double& expectM, double& CV, double& chi);
void MetropolisAlgo (int NumSpins, int **SpinMatrix,
                     double& E, double& M, double *w,
                     int& AcceptedConfigCounter);
void CalculateExpectationValuesMC (double T, int NumSpins, int Cycles,
                                   int** SpinMatrix, bool &RandomInitialSpins,
                                   double *w, int my_rank, ofstream &EnergiesFile,
                                   ofstream &ExpValsMCFile);
void CalculateExpectationValuesTemp (double T, int NumSpins, int Cycles,
                                     int** SpinMatrix, bool &RandomInitialSpins,
                                     double *w, int myloop_begin, int myloop_end,
                                     int my_rank, ofstream &ExpValsTempFile);
void output(int NumSpins, int CycleMC, double T, double *average, int AcceptedConfigCounter, ofstream &ofile);
void compare(int NumSpins, int CycleMC, double T, double *average);

int main(int argc, char *argv[])
{
    int NumSpins = 20;      // Number of spins in one dimension (L).

    double T = 1.0;        // Temperature for calculations with constant temperature.
    double StartTemp = 2.0; // Start temperature for calculations with varying temperature.
    double EndTemp = 2.4;   // End temperature for calculations with varying temperature.
    double dT = 0.05;       // Temperature step for calculations with varying temperature.

    int NumCyclesMC = 100000;   // Number of Monte Carlo cycles.
    bool RandomInitialSpins = false;    // Switch for initial cofiguration (all spins up if false, random config if true).

    // MPI initializations
    int numprocs, my_rank;
    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

    /*
    Determine number of intervall which are used by all processes
    myloop_begin gives the starting point on process my_rank
    myloop_end gives the end point for summation on process my_rank
    */
    int no_intervalls = NumCyclesMC/numprocs;
    int myloop_begin = my_rank*no_intervalls + 1;
    int myloop_end = (my_rank+1)*no_intervalls;
    if ( (my_rank == numprocs-1) &&( myloop_end < NumCyclesMC) ) myloop_end = NumCyclesMC;

    // Broadcast common variables to all nodes
    MPI_Bcast (&NumSpins, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&StartTemp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&EndTemp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&dT, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Allocate memory for spin matrix:
    int **SpinMatrix;
    SpinMatrix = new int*[NumSpins];

    for (int i=0; i<NumSpins; i++){
        SpinMatrix[i] = new int[NumSpins];
    }


    // Setup array for possible energy changes:
    double w[17];
    for( int dE =-8; dE <= 8; dE++) w[dE+8] = 0;
    for( int dE =-8; dE <= 8; dE+=4) w[dE+8] = exp(-dE/T);

    // Calculating expectation values etc. as functions of MC cycles:
    // Create datafiles:
    char EnergiesFileName[100];
    char ExpValsMCFileName[100];
    sprintf(EnergiesFileName, "Energies_T%.1f_L%d_R%d.txt", T, NumSpins, RandomInitialSpins);
    sprintf(ExpValsMCFileName, "ExpValsMC_T%.1f_L%d_R%d.txt", T, NumSpins, RandomInitialSpins);
    ofstream EnergiesFile(EnergiesFileName);
    ofstream ExpValsMCFile(ExpValsMCFileName);

    // Do calculation:
    CalculateExpectationValuesMC(T, NumSpins, NumCyclesMC, SpinMatrix,
                                 RandomInitialSpins, w, my_rank,
                                 EnergiesFile, ExpValsMCFile);

    EnergiesFile.close();
    ExpValsMCFile.close();

    // Calculating expectation values etc. as functions of temperature:
    // Create datafiles:
    char ExpValsTempFileName[100];
    sprintf(ExpValsTempFileName, "ExpValsTemp_L%d_R%d.txt", NumSpins, RandomInitialSpins);
    ofstream ExpValsTempFile(ExpValsTempFileName);

    // Do calculation with a loop over the temperature interval:
    for (T=StartTemp; T<=EndTemp; T+=dT){
        cout << "Current temperature: " << T << endl;

        // Need to re-initialize w for each temperature (due to the dependancy on T):
        for(int dE = -8; dE <= 8; dE++) w[dE+8] = 0;
        for(int dE = -8; dE <= 8; dE+=4) w[dE+8] = exp(-dE/T);

        CalculateExpectationValuesTemp(T, NumSpins, NumCyclesMC, SpinMatrix,
                                       RandomInitialSpins, w, myloop_begin,
                                       myloop_end, my_rank, ExpValsTempFile);
    }

    ExpValsTempFile.close();

    MPI_Finalize();

    return 0;
}

void initialize (int NumSpins, int **SpinMatrix, double& E, double& M)
{
    // Initialize spin matrix:
    for (int i=0; i<NumSpins; i++){
        for (int j=0; j<NumSpins; j++){
            SpinMatrix[i][j] = 1;
        }
    }

    // Initialize energy and magnetization:
    for (int i=0; i<NumSpins; i++){
        for (int j=0; j<NumSpins; j++){
            M += SpinMatrix[i][j];
            E -= SpinMatrix[i][j]*
                    (SpinMatrix[periodic(i, NumSpins, -1)][j] +
                     SpinMatrix[i][periodic(j, NumSpins, -1)]);
        }
    }
}

void initializeRandom (int NumSpins, int **SpinMatrix, double& E, double& M)
{
    // Random number generator:
    default_random_engine RNG;
    uniform_real_distribution<double> distribution(0.0,1.0);
    double RandNum;

    // Initialize spin matrix:
    for (int i=0; i<NumSpins; i++){
        for (int j=0; j<NumSpins; j++){
            RandNum = distribution(RNG);
            if (RandNum>0.5) SpinMatrix[i][j] = 1;
            else SpinMatrix[i][j] = -1;
        }
    }

    // Initialize energy and magnetization:
    for (int i=0; i<NumSpins; i++){
        for (int j=0; j<NumSpins; j++){
            M += SpinMatrix[i][j];
            E -= SpinMatrix[i][j]*
                    (SpinMatrix[periodic(i, NumSpins, -1)][j] +
                     SpinMatrix[i][periodic(j, NumSpins, -1)]);
        }
    }
}

void closedForm (int NumSpins, double T, double& Z, double& expectE,
                 double& expectM, double& CV, double& chi)
{
    // Closed form expressions:
    int N = NumSpins*NumSpins;
    Z = 12 + 4*cosh(8/T);
    expectE = -(8*sinh(8/T)/(3 + cosh(8/T)))/N;
    expectM = 8/Z*(exp(8/T)+2)/N;
    CV = 64/T*((1 + 3*cosh(8/T))/(Z/4*Z/4))/N;
    chi = (32/T*((exp(8/T) + 1)/Z) - 64/T*((exp(8/T) + 2)/Z)*((exp(8/T) + 2)/Z))/N;
}

void MetropolisAlgo (int NumSpins, int **SpinMatrix,
                     double& E, double& M, double *w,
                     int& AcceptedConfigCounter)
{
    // Random number generator:
    unsigned Seed = std::chrono::system_clock::now().time_since_epoch().count();
    default_random_engine RNG(Seed);
    uniform_real_distribution<double> distribution(0.0,1.0);

    //AcceptedConfigCounter = 0;

    // Loop over all spins:
    for (int i=0; i<NumSpins; i++){
        for (int j=0; j<NumSpins; j++){
            // Find random position
            int ix = (int) (distribution(RNG)*(double)NumSpins);
            int iy = (int) (distribution(RNG)*(double)NumSpins);
            int deltaE = 2*SpinMatrix[iy][ix]*
                    (SpinMatrix[iy][periodic(ix, NumSpins, -1)]+
                     SpinMatrix[periodic(iy, NumSpins, -1)][ix]+
                     SpinMatrix[iy][periodic(ix, NumSpins, 1)] +
                     SpinMatrix[periodic(iy, NumSpins, 1)][ix]);
            // Metropolis test
            if (distribution(RNG) <= w[deltaE+8]){
                SpinMatrix[iy][ix] *= -1;   // Flip one spin and accept the new spin config.
                // Update the energy and the magnetization
                M += (double) 2*SpinMatrix[iy][ix];
                E += (double) deltaE;
                AcceptedConfigCounter += 1;
            }
        }
    }
}

void CalculateExpectationValuesMC (double T, int NumSpins, int Cycles,
                                   int** SpinMatrix, bool &RandomInitialSpins,
                                   double *w, int my_rank, ofstream &EnergiesFile,
                                   ofstream &ExpValsMCFile)
{
    double Average[5];              // Array for expectation values.
    int AcceptedConfigCounter = 0;  // Counter for accepted configurations.
    double E = 0;                   // Energy
    double M = 0;                   // Magnetization

    int TotalNumSpins = NumSpins*NumSpins;

    for( int i=0; i<5; i++) Average[i] = 0.;

    // Initialize:
    if (RandomInitialSpins){initializeRandom(NumSpins, SpinMatrix, E, M);}
    else {initialize(NumSpins, SpinMatrix, E, M);}

    // Loop over MC cycles:
    for (int Cycle=1; Cycle<=Cycles; Cycle++){
        cout << "Cycle: " << Cycle << endl;
        MetropolisAlgo(NumSpins, SpinMatrix, E, M, w, AcceptedConfigCounter);
        EnergiesFile << E/TotalNumSpins << endl;

        // Update expectation values
        Average[0] += E; Average[1] += E*E;
        Average[2] += M; Average[3] += M*M;
        Average[4] += fabs(M);

        // Save calculated results to file:
        if (my_rank==0){output(NumSpins, Cycle, T, Average, AcceptedConfigCounter, ExpValsMCFile);}
    }
    // Compare calculated results and closed form results:
    if (NumSpins==2){
        compare(NumSpins, Cycles, T, Average);
    }
}

void CalculateExpectationValuesTemp (double T, int NumSpins, int Cycles,
                                     int** SpinMatrix, bool &RandomInitialSpins,
                                     double *w, int myloop_begin, int myloop_end,
                                     int my_rank, ofstream &ExpValsTempFile)
{
    double Average[5];              // Array for expectation values calculated by a single node.
    double TotalAverage[5];         // Array for total expectation values from all nodes.
    int AcceptedConfigCounter = 0;  // Counter for accepted configurations.
    double E = 0;                   // Energy
    double M = 0;                   // Magnetization

    for( int i=0; i<5; i++) Average[i] = 0.;
    for( int i=0; i<5; i++) TotalAverage[i] = 0.;

    // Initialize:
    if (RandomInitialSpins){initializeRandom(NumSpins, SpinMatrix, E, M);}
    else {initialize(NumSpins, SpinMatrix, E, M);}

    // Loop over Monte Carlo cycles:
    for (int Cycle=myloop_begin; Cycle<=myloop_end; Cycle++){
        MetropolisAlgo(NumSpins, SpinMatrix, E, M, w, AcceptedConfigCounter);

        // Update expectation values
        Average[0] += E; Average[1] += E*E;
        Average[2] += M; Average[3] += M*M;
        Average[4] += fabs(M);
    }

    // Combine the calculated data from all nodes to find the total expectation values
    for (int i=0; i<5; i++){
        MPI_Reduce(&Average[i], &TotalAverage[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }

    // Master node calls output function to save results to file:
    if (my_rank==0){output(NumSpins, Cycles, T, TotalAverage, AcceptedConfigCounter, ExpValsTempFile);}
}

void output(int NumSpins, int CycleMC, double T, double *average, int AcceptedConfigCounter, ofstream &ofile)
{
    double norm = 1/((double) (CycleMC)); // divide by total number of cycles.
    double expectE = average[0]*norm;     // Expectation value of energy.
    double expectE2 = average[1]*norm;    // Expectation value of E*E.
    double expectM = average[2]*norm;     // Expectation value of magnetization.
    double expectM2 = average[3]*norm;    // Expectation value of M*M.
    double expectAbsM = average[4]*norm;  // Expectation value of absolute value of magnetization.

    double varianceE = expectE2- expectE*expectE;           // Variance of energy.
    double varianceAbsM = expectM2-expectAbsM*expectAbsM;   // Variance of magnetization.

    double CV = varianceE/T/T;            // Specific heat.
    double chi = varianceAbsM/T;          // Susceptibility.

    // all expectation values are per spin, divide the values we're interested in by 1/NumSpins/NumSpins
    expectE /= NumSpins*NumSpins;
    expectAbsM /= NumSpins*NumSpins;
    varianceE /= NumSpins*NumSpins;
    CV /= NumSpins*NumSpins;
    chi /= NumSpins*NumSpins;

    // Write results to file:
    ofile << T << "  " << CycleMC << "  " << AcceptedConfigCounter*norm << "  "
          << expectE << "  " << expectAbsM << "  " << varianceE << "  "
          << CV << "  " << chi << endl;
}

void compare(int NumSpins, int CycleMC, double T, double *average)
{
    double norm = 1/((double) (CycleMC)); // divided by total number of cycles
    double expectE = average[0]*norm;     // Expectation value of energy.
    double expectE2 = average[1]*norm;    // Expectation value of E*E.
    double expectM = average[2]*norm;     // Expectation value of magnetization.
    double expectM2 = average[3]*norm;    // Expectation value of M*M.
    double expectAbsM = average[4]*norm;  // Expectation value of absolute value of magnetization.

    double varianceE = expectE2- expectE*expectE;           // Variance of energy.
    double varianceAbsM = expectM2-expectAbsM*expectAbsM;   // Variance of magnetization.

    double CV = varianceE/T/T;            // Specific heat.
    double chi = varianceAbsM/T;          // Susceptibility.

    // all expectation values are per spin, divide the values we're interested in by 1/NumSpins/NumSpins
    expectE /= NumSpins*NumSpins;
    expectAbsM /= NumSpins*NumSpins;
    varianceE /= NumSpins*NumSpins;
    CV /= NumSpins*NumSpins;
    chi /= NumSpins*NumSpins;

    // Calculate closed form vlaues:
    double Z_CF, expectE_CF, expectM_CF, CV_CF, chi_CF;
    closedForm(NumSpins, T, Z_CF, expectE_CF, expectM_CF, CV_CF, chi_CF);

    cout << "Values from closed form expressions and numerical calculations" << endl << "(all values are per spin):"
         << endl << "Value:   " << "Closed Form:   " << "Numerical:"
         << endl << "<E>      " << expectE_CF << "       " << expectE
         << endl << "<M>      " << expectM_CF << "       " << expectAbsM
         << endl << "CV       " << CV_CF << "      " << CV
         << endl << "chi      " << chi_CF << "     " << chi << endl;
}

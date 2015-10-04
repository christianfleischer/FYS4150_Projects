#include <iostream>
#include <cmath>
#include <armadillo>

using namespace std;
using namespace arma;

void JacobiRotations (double MaxValue,mat &A,mat &absA,int &NumberOfRotations,double Tolerance,uword RowIndexMax,uword ColIndexMax,int N)
{
    uword k = RowIndexMax;
    uword l = ColIndexMax;
    double s, c, t, tau;
    double A_kk, A_ll, A_ik, A_il;

    while (MaxValue>Tolerance){
        tau = (A(l,l) - A(k,k))/(2*A(k,l));

        //Computing t = tan(theta), choosing the smallest root:
        if (tau > 0){
            t = 1.0/(tau + sqrt(1.0 + tau*tau));
        }
        else {
            t = -1.0/(-tau + sqrt(1.0 + tau*tau));
        }

        c = 1.0/sqrt(1 + t*t);
        s = c*t;

        //Using the similarity transformation on A and replacing A with the transformed matrix:
        A_kk = A(k,k);
        A_ll = A(l,l);

        A(k,k) = c*c*A_kk - 2*c*s*A(k,l) + s*s*A_ll;
        A(l,l) = s*s*A_kk + 2*c*s*A(k,l) + c*c*A_ll;
        A(k,l) = 0;
        A(l,k) = 0;

        for (int i=0; i<N-1; i++){
            if (i != k && i!= l){
                A_ik = A(i,k);
                A_il = A(i,l);
                A(i,k) = c*A_ik - s*A_il;
                A(k,i) = A(i,k);
                A(i,l) = c*A_il + s*A_ik;
                A(l,i) = A(i,l);
            }
        }

        //Find the new largest element in A:
        absA = abs(A);
        absA.diag(0) = zeros(N-1);
        MaxValue = absA.max(k,l);

        //Count the number of similarity transformations:
        NumberOfRotations++;

    }
    return;
}

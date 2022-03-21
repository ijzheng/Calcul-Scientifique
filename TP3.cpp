#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>
#include <ctime>
#include <math.h>


using namespace std;


double *zero(int N)
{
    double *b = new double [N];
    for (int i=0;i<N;i++){
        b[i]=0;
    }
    return b;
}

double *one(int N)
{
    double *b = new double [N];
    for (int i=0;i<N;i++){
        b[i]=1;
    }
    return b;
}


double *linspace(double a, double b, int N) {
	double dx = (b - a) / (N - 1);
	double *v = zero(N);
	for (int i = 0; i < N; i++) {
		v[i] = a+i*dx;
	}
	return v;
}

double **Zero(int N)
{
    double **A = new double *[N];
    for (int i=0;i<N;i++){
        A[i] = new double [N];
        for(int j=0;j<N;j++){
            A[i][j]=0;
            
        }
    }
    return A;
}

double ps(double*u,double*v,int N)
{
    double r=0;
    for(int i=0;i<N;i++){
        r+=u[i]*v[i];
    }
    return r;
}


double **pm(double **A,double **B,int N){
    double **C = Zero(N);
    for (int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            for (int k=0;k<N;k++){
                C[i][j] += A[i][k]*B[k][j];
            }
            
        }
    }
    return C;
}

double *pmv(double **A,double *v,int N){
    double *u = zero(N);
    for (int i=0;i<N;i++){
        u[i] = ps(A[i],v,N);
    }
    return u;
}

double **am(double **A, double **B,int N){
    double **C = Zero(N);
    for (int i=0;i<N;i++){
        for (int j=0;j<N;j++){
            C[i][j] = A[i][j]+B[i][j];
        }
    }
    return C;
}

double *av(double *u,double *v,int N){
    double *r = zero(N);
    for (int i=0;i<N;i++){
        r[i] = u[i]+v[i]; 
    }
    return r;
}

double **prm(double **A, double b, int N){
    for (int i=0;i<N;i++){
        for (int j=0;j<N;j++){
            A[i][j] *= b;
        }
    }
    return A;
}

double *prv(double *v,double a,int N){
    for (int i=0;i<N;i++){
        v[i] *= a;
    }
    return v;
}

double norme(double *u, int N){
    double r = 0;
    for(int i=0;i<N;i++){
        r += u[i]*u[i];
    }
    r = sqrt(r);
    return r;
}


void print_m(double**A,int N)
{
    for( int i =0; i<N ; i++){
        cout<<"ligne"<<i<<": " ;
        for ( int j =0; j<N ; j++){
            cout<<A[i][j]<<" " ;
        }
    cout<<endl ;
    }
}


void print_v(double*A,int N)
{
    cout<<"[";
    for( int i =0; i<N-1 ; i++){
        cout<<A[i]<<"," ;
    }
    cout<<A[N-1]<<"]";
}



double *jacobi(double **A,double *b, double *x0,int N){
    int maxiter = 10000;
    double tol = 0.00000001;
    double **S = Zero(N);
    double **T = Zero(N);
    double *x = x0;
    for (int i=0;i<N;i++){
        for (int j=0;j<N;j++){
            if (i==j){
                S[i][j] = A[i][j];
            }
            else{
                T[i][j] = A[i][j]*(-1);
            }
        }
    }
    double **invS = Zero(N);
    for(int i=0;i<N;i++){
        invS[i][i] = 1/S[i][i];
    }
    double *xx = zero(N);
    for (int i=0;i<maxiter;i++){
        xx = av(pmv(pm(invS,T,N),x,N),pmv(invS,b,N),N);
        if ((norme(av(xx,prv(x,-1,N),N),N))-tol < 0){
            break;
        }
        x = prv(xx,1,N);
        if (i == maxiter-1){
            return zero(N);
        }
    }
    double *sol = xx;
    return sol;
}

double f(double x){
    double r = 4*x*x*x+ 3*x*x + 2*x+1;
    return r;
}

double *moment(double a,double b,int N){
    double *x = linspace(a,b,N);
    double *M = one(N-2);
    double h = (b-a)/(N-1);
    double *y = zero(N-2);
    double *yy = zero(N);
    yy[0] = f(a);
    yy[N-1] = f(b);
    for (int i=1;i<N-1;i++){
        y[i-1] = f(x[i]);
    }
    double  **A = Zero(N-2);
    for (int i=0;i<N-2;i++){
        for (int j=0;j<N-2;j++){
            if (i==j){
                A[i][j] = 4;
            }
            if(i-j==1||i-j==-1){
                A[i][j] = 1;
            }
        }
    }
    double **B = Zero(N-2);
    for (int i=0;i<N-2;i++){
        for (int j=0;j<N-2;j++){
            if (i==j){
                B[i][j] = -12/h/h;
            }
            if(i-j==1||i-j==-1){
                B[i][j] = 6/h/h;
            }
        }
    }
    double *u = pmv(B,y,N-2);
    double *R = jacobi(A,u,M,N-2);
    double *sol = zero(N);
    for (int i=1;i<N-1;i++){
        sol[i] = R[i-1];
    }
    return sol;
}








int main(){
    int N = 10;
    int a = -5;
    int b = 5;
    double g = 0.001;
    double *sol = moment(a, b, N);
    print_v(sol,N);
}
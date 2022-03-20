#include <iostream>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>
#include <ctime>

using namespace std;


double *zero(int N)
{
    double *b = new double [N];
    for (int i=0;i<N;i++){
        b[i]=0;
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



double **jacobi(double **A,int N){
    double **L = prm(A,2,N);
    return L;
}

//double *moment(f, double a,double b,int N){

//}








int main(){
    int N = 10;
    double **A = Zero(N);
    double *u = zero(N);
    double *v = zero(N);
    double *r = av(u,v,N);
    print_v(A[2],N);
}
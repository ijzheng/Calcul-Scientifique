#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <time.h>
#include <ctime>

using namespace std;

double *zero(int N) {
    double *b = new double [N];
    for (int i=0;i<N;i++){
        b[i]=0;
    }
    return b;
}
void print_v(double*A,int N)
{
    cout<<"[";
    for( int i =0; i<N-1 ; i++){
        cout<<A[i]<<"," ;
    }
    cout<<A[N-1]<<"]";
}

double *euler_explicite(int T, int N){
    double *y = zero(N+1);
    y[0]=1;
    double h = T/N;
    for (int i=0;i<N;i++){
        double ti=i*h;
        y[i+1]=(1-2*h*ti)*y[i];
    }
    return y;
}
double *euler_implicite(int T, int N){
    double *y =zero(N+1);
    double h = T/N;
    y[0]=1;
    for (int i=0;i<N;i++){
        double ti1=(i+1)*h;
        y[i+1]=y[i]* (1/(1+2*h*ti1));
    }
    return y;
}
double *crank_nicholson(int T,int N){
    double h=T/N;
    double *y =zero(N+1);
    y[0]=1;
    for (int i=0;i<N;i++){
        double ti=i*h;
        double ti1=(i+1)*h;
        y[i+1]=y[i]*(1+ti)/(1+h*ti1);
    }
    return y;

}
double *Heun(int T,int N){
    double h=T/N;
    double *y =zero(N+1);
    y[0]=1;
    for (int i=0;i<N;i++){
        double ti=i*h;
        double ti1=(i+1)*h;
        y[i+1]=y[i]*(1-(ti+ti1)*h)+2*h*h*ti*ti1;
    }
    return y;

}


int main(){
    int T=300;
    int N=200;
    //print_v(euler_explicite(T,N),N+1);
    //print_v(euler_implicite(T,N),N+1);
    //print_v(Heun(T,N),N+1);
    print_v(crank_nicholson(T,N),N+1);

}

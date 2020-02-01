#include <iostream>
#include <cmath>
using namespace std;
/*
 * abs(a-b)
 * abs((a-b)/a)
 * abs(m-manterior)
 * abs((m-manterior)/m)
 * */
double f(double x){return 5*pow(x,2)+6*x-2;}
double f1(double x){return pow(x,2)-x-1.2;}
double f1l(double x){return 2*x-1;}
double g(double x){return sqrt(x+1.2);}
double gl(double x){return 1.0/(2*sqrt(x-1.2));}
void Bisseccao(double a, double b, double prec){
    double m=0, mseguinte; int i=0;
    do{
        m=(a+b)/2.0;
        if (f(a)*f(m)<0)
            b = m;
        else
            a = m;
        mseguinte=(a+b)/2.0;
        cout << "Iteracao n " << ++i << " - a: " << a << "   b: " << b << endl;
    }while(abs((mseguinte-m))>prec);
}

void Newton(double xo){
    int i=1; double nextx;
    do{
        nextx = xo - f1(xo)/f1l(xo);
        cout << "x: " << xo << "   f(x): " << f1(xo) << "   fl(x): " << f1l(xo) << "  it n: " << i << endl;
        xo = nextx;
    }while(++i != 6);
}

void Picard(double xo){
    int i=1; double nextx;
    do{
        nextx = g(xo);
        cout << "x: " << xo << "   g(x): " << g(xo) << "  it n: " << i << endl;
        xo = nextx;
    }while(++i != 6);
}
void Corda(double a, double b, double prec){
    double m=0, mseguinte; int i=0;
    do{
        m=(a*f(b) - b*f(a))/(f(b)-f(a));
        if (f(a)*f(m)<0)
            b=m;
        else
            a=m;
        mseguinte = (a*f(b) - b*f(a))/(f(b)-f(a));
        cout << "Iteracao n " << ++i << " - a: " << a << "   b: " << b << " mseg-m: " <<abs(m-mseguinte) << " prec: " << prec <<endl;
    }while(abs(mseguinte-m)>prec);
}


int main() {
    cout << "Hello, World!" << endl;
    /*Bisseccao(-2,-1,pow(10,-4));*/
    Corda(-2,-1,pow(10,-4));
    /*Newton(4);
    Picard(4);*/
    return 0;
}
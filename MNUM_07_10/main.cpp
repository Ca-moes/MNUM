#include <iostream>
#include <math.h>
//#include <cmath>
using namespace std;

double f3(double x){return exp((x-5.0)/2.0);}
double f2(double x) {return 2*log(x)+ 5;}
double f(double x){return exp(0.7*x)-pow(x,2)-0.5;}
double fl(double x){ return (2.0/3.0)*M_PI*(9.0-x)*x - M_PI*pow(x,2.0)/3.0;}
double g1(double x, double y) {return sqrt((x*y+5.0*x-1.0)/2.0);}
double g2(double x, double y) {return sqrt(x+3*log(x));}
double g1dx(double x, double y){return (0.3535533905932738 *(y + 5.0))/(sqrt(x*y + 5.0*x - 1.0));}
double g1dy(double x, double y){return (0.3535533905932738*x)/(sqrt(x*y + 5.0*x - 1.0));}
double g2dx(double x, double y){return (3/x +1)/(2 * sqrt(3*log(x)+x));}
double g2dy(double x, double y){return 0;}

void Bisseccao(double a, double b, double erro){
    // intervalo [a,b]
    double m=0, mant;
    int i =0;
    cout << "Inicial - " << "  a: " << a << "   b: " << b << "   m: " << (a+b)/2 << "  f(a): " << f(a) << "   f(b): " << f(b) << "   f(m): " << f((a+b)/2) << endl;
    do{
        mant = m;
        m = (a+b)/2;
        if (f(a)*f(m)<0)
            b = m;
        else
            a = m;
        cout << "Iteracao n" << ++i << "  a: " << a << "   b: " << b << "   m: " << (a+b)/2 << "  f(a): " << f(a) << "   f(b): " << f(b) << "   f(m): " << f((a+b)/2) << endl;
    }while(abs(m - mant) >= erro);
    cout << endl;
}

void Corda(double a, double b, double erro){
    double w = 0, want;
    int i=0;
    do{
        want = w;
        w = (a*f(b) - b*f(a))/(f(b)-f(a));
        if (f(a)*f(w)<0)
            b = w;
        else
            a = w;
        cout << "Iteracao n" << ++i << "  a: " << a << "   b: " << b << endl;
    }while(abs(w - want) >= erro);
    cout << endl;
}

void Tangente(double a, double erro){
    double nextx = a ;
    int i=0;
    //cout << "Iteracao n" << ++i << "  x: " << a << endl;
    do{
        a = nextx;
        nextx = a - f(a)/fl(a);
        cout << "Iteracao n" << ++i << "  x: " << a << endl;
    }while(abs(nextx - a) >= erro);
    cout << endl;
}

void Picard(double a, double erro){
    double nextx = a, currentx;
    //cout << "current x: " << currentx << "    nextx: " << nextx << endl;
    do{
        currentx = nextx;
        nextx = f3(currentx);
        cout << "current x: " << currentx << "    nextx: " << nextx << endl;
    }while(abs(nextx - currentx) >= erro);
}
void Picard_sistema(double x, double y, double erro){
    double currentx, currenty, nextx = x, nexty = y;
    //out << "current x: " << currentx << "    nextx: " << nextx << "   current y: " << currenty << "    nexty: " << nexty << endl;
    do{
        currentx = nextx;
        currenty = nexty;
        nextx = g1(currentx, currenty);
        nexty = g2(currentx, currenty);
        cout << "current x: " << currentx << "    nextx: " << nextx << "   current y: " << currenty << "    nexty: " << nexty << endl;
    }while(abs(nextx - currentx) >= erro && abs(nexty - currenty) >= erro);
}
void Newton_sistema(double x, double y, double erro){
    double nextx = x - (g1(x,y) * g2dy(x,y) - g2(x,y) * g1dy(x,y))/(g1dx(x,y) * g2dy(x,y) - g2dx(x,y) * g1dy(x,y));
    double nexty = y - (g2(x,y) * g1dx(x,y) - g1(x,y) * g2dx(x,y))/(g1dx(x,y) * g2dy(x,y) - g2dx(x,y) * g1dy(x,y));
    int i=0;
    cout << "Iteracao n" << ++i << "  x: " << x << "   y: " << y << endl;
    while(abs(nextx - x) >= erro && abs(nexty - y) >= erro){
        x = nextx;
        y = nexty;
        nextx = x - (g1(x,y) * g2dy(x,y) - g2(x,y) * g1dy(x,y))/(g1dx(x,y) * g2dy(x,y) - g2dx(x,y) * g1dy(x,y));
        nexty = y - (g2(x,y) * g1dx(x,y) - g1(x,y) * g2dx(x,y))/(g1dx(x,y) * g2dy(x,y) - g2dx(x,y) * g1dy(x,y));
        cout << "Iteracao n" << ++i << "  x: " << x << "   y: " << y << endl;
    }
    cout << endl;
}
int main() {
    Bisseccao(-1.0,0,pow(10,-5));
    /*Corda(1.9,2.5,pow(10,-5));
    Tangente(1.9, pow(10,-5));
    Picard(5, pow(10,-5));
    Picard_sistema(4,4, pow(10,-3));
    Newton_sistema(4, 4,pow(10,-3));*/
    double a =5*cos(10);
    return 0;
}
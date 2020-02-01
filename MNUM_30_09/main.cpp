#include <iostream>
#include <math.h>

using namespace std;
double f(double x){
    return x - 2*log(x) -5;
}

void Corda(double a, double b, double error){
    double xn = (f(a)*b - f(b)-a) / (f(a)-f(b)), oldxn;
    cout << "error -> " << error << endl;
    do {
        cout << "a: " << a << "   b: " << b << "   f(a): " << f(a) << "  f(b): " << f(b);
        if (f(a) * f(xn) < 0)
            b = xn;
        else
            a = xn;
        oldxn = xn;
        xn = (f(a)*b - f(b)-a) / (f(a)-f(b));
        cout << "    Interval: xn:  " << xn << " xnanterior:" << oldxn << "    diff= " << abs(xn-oldxn) << endl;
    }while(abs(xn-oldxn)>error);

    if(b > a)
        cout << "\n\nINTERVALO -> " << a << " - " << b <<endl;
    else
        cout << "\n\nINTERVALO -> " << b << " - " << a <<endl;
}

void Bi(double a, double b,double error){
    double xn = (a+b)/2, oldxn;
    cout << "error -> " << error << endl;
    do {
        cout << "a: " << a << "   b: " << b << "   f(a): " << f(a) << "  f(b): " << f(b);
        if (f(a) * f(xn) < 0)
            b = xn;
        else
            a = xn;
        oldxn = xn;
        xn = (a+b)/2;
        cout << "    Interval: xn:  " << xn << " xnanterior:" << oldxn << "    diff= " << abs(xn-oldxn) << endl;
    }while(abs(xn-oldxn)>error);

    if(b > a)
        cout << "\n\nINTERVALO -> " << a << " - " << b <<endl;
    else
        cout << "\n\nINTERVALO -> " << b << " - " << a <<endl;
}

void Picard(double a){

}
int main() {
    cout << f(0.05) << "" << f(2) << endl;
    Bi(0.05,2,pow(10,-3));
    Corda(0.05,2,pow(10,-3));
    return 0;
}
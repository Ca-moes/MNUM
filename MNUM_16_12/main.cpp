#include <iostream>
#include "cmath"

using namespace std;

const double B = (sqrt(5)-1.0)/2.0;
const double A = pow((sqrt(5)-1.0)/2.,2);

double f(double x){return pow(2*x+1,2) - 5*cos(10*x);}
double f2(double x, double y){return y*y - 2*x*y - 6*y + 2*x*x + 12;}
double dfdx(double x, double y){return -2*y + 4*x;}
double dfdy(double x,double y){return 2*y - 2*x -6;}

void Regra_Aurea_min(){
    double x1 = -1, x2 = 0, x3, x4;
    do {
        x3=x1 + A*(x2-x1); x4 = x1 + B*(x2-x1);
        f(x3) < f(x4) ? x2 = x4 : x1 = x3;  // para achar máximo é trocar o < para >
    }while (abs(x2-x1)>=pow(10,-3));
    f(x1) < f(x2) ? x1=x1 : x1=x2;
    cout << "Valor x onde esta minimo -> " << x1 << endl;
}
void Gradiente(){
    // fazer deltaf á parte, calcular derivadas e por em cima.
    double y=0, x=0, h=1, xn = 0, yn=0;
    bool diminui = false;
    do{
        x = xn;
        y = yn;
        xn = x - h*(dfdx(x,y));
        yn = y - h*(dfdy(x,y));
/*        if (f2(xn,yn) < f2(x,y)){
            h = h*2;
        } else{
            h = h/2.0;
            xn = x;
            yn = y;
        }*/
        f2(xn,yn) < f2(x,y) ? h *= 2 : h /= 2;
    }while(abs(xn-x)>=pow(10,-3) || abs(yn-y)>=pow(10,-3) || xn==0 || yn==0);
    cout << "x -> " << x << "\ty -> " << y;
}

void Quadratica(){
    // calcular hessiana
    double y=0, x=0, h=1.0/4.0, xn = 0, yn=0;
    do{
        x = xn;
        y = yn;
        xn = x - h*(dfdx(x,y));
        yn = y - h*(dfdy(x,y));
        f2(xn,yn) < f2(x,y) ? h *= 2 : h /= 2;
    }while(abs(xn-x)>=pow(10,-3) || abs(yn-y)>=pow(10,-3));
    cout << "x -> " << x << "\ty -> " << y;
}
int main() {
    //Regra_Aurea_min();
    Gradiente();
    return 0;
}

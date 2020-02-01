#include <iostream>
#include <cmath>
using namespace std;

double Gx1(double x2, double x3){return (7-x2-x3)/3.0;}
double Gx2(double x1, double x3){return (4-x1-(2*x3))/4.0;}
double Gx3(double x2){return (5-(2*x2))/5.0;}
double f(double x){return sin(x)/pow(x,2);}

void Gauss_Jacobi(double vix1, double vix2, double vix3){
    //vix1 - valor inicial x1
    cout << "METODO GAUSS-JACOBI\n";
    double x1 = Gx1(vix2,vix3), x2 = Gx2(vix1,vix3), x3 = Gx3(vix2), nx1, nx2, nx3;
    int i = 0;
    cout << "Iteracao n " << ++i << "  - x1: " << x1 << "   x2: " << x2 << "    x3: " << x3 << endl;
    do{
        nx1 = Gx1(x2,x3); nx2 = Gx2(x1,x3); nx3 = Gx3(x2);
        cout << "Iteracao n " << ++i << "  - x1: " << nx1 << "   x2: " << nx2 << "    x3: " << nx3 << endl;
        if(abs(nx1-x1)<pow(10,-3) && abs(nx2-x2)<pow(10,-3) &&abs(nx3-x3)<pow(10,-3)) break;
        x1 = nx1; x2 = nx2; x3 = nx3;
    }while(true);
}

void Gauss_Seidel(double vix1, double vix2, double vix3){
    //vix1 - valor inicial x1
    cout << "METODO GAUSS-SEIDEL\n";
    double x1 = Gx1(vix2,vix3), x2 = Gx2(vix1,vix3), x3 = Gx3(vix2), nx1, nx2, nx3;
    int i = 0;
    cout << "Iteracao n " << ++i << "  - x1: " << x1 << "   x2: " << x2 << "    x3: " << x3 << endl;
    do{
        nx1 = Gx1(x2,x3); nx2 = Gx2(nx1,x3); nx3 = Gx3(nx2);
        cout << "Iteracao n " << ++i << "  - x1: " << nx1 << "   x2: " << nx2 << "    x3: " << nx3 << endl;
        if(abs(nx1-x1)<pow(10,-3) && abs(nx2-x2)<pow(10,-3) &&abs(nx3-x3)<pow(10,-3)) break;
        x1 = nx1; x2 = nx2; x3 = nx3;
    }while(true);
    cout << "\n";
}

void Trapezios(int n, double x0, double xn, double h){
    // usar função f em cima
    double hinicial = h, h_2 = h/2.0, soma = 0;
    for (int i = 1; i <= n-1; ++i) {
        soma = soma + f(h*i);
    }
    soma = soma * 2;
    soma = soma + f(x0) + f(xn);
    cout << "Resultado: "<< h_2 * soma << endl;
}
int main() {
    /*Gauss_Jacobi(0,0,0);
    Gauss_Seidel(0,0,0);*/
    Trapezios(4,M_PI_2, M_PI,(M_PI - M_PI_2)/4);
    return 0;
}
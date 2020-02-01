#include <iostream>
#include <cmath>
using namespace std;

double y(double valor){return sin(valor);}
double y3d(double x, double y){return sin(x+y);}

void Simpson(int n){
    // h/3 * (y0 + 4*(\sum_(i=1,3,5)^2x-1 yi) + 2*(\sum_(i=2)^2x-2 yi) + yn)
    double h = M_PI/(2.0*n);
    double tempsoma=0;
    for (int i = 1; i <= 2*n-1; i = i+2) {
        cout << "valor a somar : "<< y(i*h) << "    valor de tempsoma : " << tempsoma << endl;
        tempsoma = tempsoma + 4*y(i * h);
    }
    for (int j = 2; j <= 2*n-2; j = j+2) {
        cout << "2 valor a somar : "<< y(j*h) << "    valor de tempsoma : " << tempsoma << endl;
        tempsoma = tempsoma + 2*y(j * h);
    }
    tempsoma += y(0);
    tempsoma += y(M_PI);
    tempsoma = tempsoma * (h/3);
    cout << tempsoma << endl;
}

void Simpson3D(){
    int n = 2;
    double hx = (M_PI/2.0)/(n);
    double hy = (M_PI/2.0)/(n);
    double htot = (hx-hy)/9;
    double tempsoma=0;
    for (int i = 0; i <= 2; i+=2) {
        for (int j = 0; j <= 2 ; j+=2) {
            cout << "valor a somar : "<< y3d(i*hx,j*hy) << "    valor de tempsoma : " << tempsoma << endl;
            tempsoma = tempsoma + y3d(i*hx,j*hy);
        }
    }
    for (int i = 0; i <= 2; i++) {
        for (int j = 0; j <= 2 ; j++) {
            if((i==0 && j==1) || (i==2 && j == 1)){
                cout << "valor a somar : "<< 4*y3d(i*hx,j*hy) << "    valor de tempsoma : " << tempsoma << endl;
                tempsoma = tempsoma + 4*y3d(i*hx,j*hy);
            }
            if(i==1 && j!=1){
                cout << "valor a somar : "<< 4*y3d(i*hx,j*hy) << "    valor de tempsoma : " << tempsoma << endl;
                tempsoma = tempsoma + 4*y3d(i*hx,j*hy);
            }
        }
    }
    cout << "valor a somar : "<< 16*y3d(1*hx,1*hy) << "    valor de tempsoma : " << tempsoma << endl;
    tempsoma = tempsoma + 16*y3d(1*hx,1*hy);

    tempsoma = tempsoma * htot;
    cout << tempsoma << endl;
}
int main() {
    //Simpson(4);
    Simpson3D();
    std::cout << "Hello, World!" << std::endl;
    return 0;
}
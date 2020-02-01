#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;
//Pergunta 1
double dz(double t, double y, double z){
    return 0.5 + t*t + t*z;
}
double dy(double t, double y, double z){
    return z;
}
void Euler(){
    double h=0.25, t=0, y=0, z=1, znext, ynext, tnext; int i=0;
    cout << i++ << "   " << t << "   " << y << endl;
    znext = z + h*dz(t,y,z);
    ynext = y + h*dy(t,y,z);
    tnext = t + h;
    z = znext; y = ynext; t = tnext;
    cout << i++ << "   " << t << "   " << y << endl;
    znext = z + h*dz(t,y,z);
    ynext = y + h*dy(t,y,z);
    tnext = t + h;
    z = znext; y = ynext; t = tnext;
    cout << i++ << "   " << t << "   " << y << endl;
}
void RK4(){
    double h=0.25, t=0, y=0, z=1, znext=z, ynext=y, tnext=t; int i=0;
    double dely1, dely2,dely3,dely4,delz1,delz2,delz3,delz4;
    cout << i++ << "   " << t << "   " << y << endl;
    delz1 = h*dz(t,y,z);
    dely1 = h*dy(t,y,z);
    delz2 = h*dz(t + h/2, y+dely1/2.0, z+delz1/2.0);
    dely2 = h*dy(t + h/2, y+dely1/2.0, z+delz1/2.0);
    delz3 = h*dz(t + h/2, y+dely2/2.0, z+delz2/2.0);
    dely3 = h*dy(t + h/2, y+dely2/2.0, z+delz2/2.0);
    delz4 = h*dz(t+h, y+dely3, z+delz3);
    dely4 = h*dy(t+h, y+dely3, z+delz3);
    znext = z + delz1/6.0 + delz2/3.0 + delz3/3.0 + delz4/6.0;
    ynext = y + dely1/6.0 + dely2/3.0 + dely3/3.0 + dely4/6.0;
    tnext = t + h;
    z = znext; y = ynext; t = tnext;
    cout << i++ << "   " << t << "   " << y << endl;
    delz1 = h*dz(t,y,z);
    dely1 = h*dy(t,y,z);
    delz2 = h*dz(t + h/2, y+dely1/2.0, z+delz1/2.0);
    dely2 = h*dy(t + h/2, y+dely1/2.0, z+delz1/2.0);
    delz3 = h*dz(t + h/2, y+dely2/2.0, z+delz2/2.0);
    dely3 = h*dy(t + h/2, y+dely2/2.0, z+delz2/2.0);
    delz4 = h*dz(t+h, y+dely3, z+delz3);
    dely4 = h*dy(t+h, y+dely3, z+delz3);
    znext = z + delz1/6.0 + delz2/3.0 + delz3/3.0 + delz4/6.0;
    ynext = y + dely1/6.0 + dely2/3.0 + dely3/3.0 + dely4/6.0;
    tnext = t + h;
    z = znext; y = ynext; t = tnext;
    cout << i++ << "   " << t << "   " << y << endl;
}
//Pergunta 3
double Z(double x, double y){return 3*x*x - x*y + 11*y + y*y - 8*x;}
double dZdx(double x, double y){return 6*x -y -8;}
double dZdy(double x, double y){return -x +11 +2*y;}
void grad(){
    double h=0.5, x=2, y=2, z = Z(x,y), gradx=dZdx(x,y), grady = dZdy(x,y);
    int i=0;
    cout << i++ << "\t" << x << "\t   " << z << "\t" << gradx <<"\n\t" << y << "\t\t\t" << grady << "\n\n";
    x = x - h*gradx;
    y = y - h*grady;
    z = Z(x,y);
    cout << i++ << "\t" << x << "\t   " << z << "\n\t" << y << endl;
}
//Pergunta 4
double y(double x){return exp(x*1.5);}
double trapezios(double h){
    double x0 = 1, x2n = 1.5, y0 = y(x0), y2n = y(x2n), sumpar=0, sumimpar=0;
    bool impar = true;
    for (double i = x0+h; i < x2n; i = i+h) {
        impar ? sumimpar += y(i) : sumpar += y(i);
        impar = !impar;
    }
    return (h/3.0) * (y0 + y2n + 4*sumimpar + 2*sumpar);
}
void trapezios3times(){
    double h = 0.125;
    cout << h << "\t" << trapezios(h) << endl;
    cout << h/2 << "\t" << trapezios(h/2) << endl;
    cout << h/4 << "\t" << trapezios(h/4) << endl;
    cout << "QC=  " << (trapezios(h/2) - trapezios(h))/(trapezios(h/4) - trapezios(h/2)) << "  ~ 16" << endl;
    cout << "err= " << trapezios(h/4) - trapezios(h/2) << endl;
}
//Pergunta 5
double newtf(double x){return (x-3.7)+pow(cos(x+1.2),3);}
double newtfl(double x){return 1 - 3*sin(x+1.2)*pow(cos(x+1.2),2);}
void Newton(){
    double x = 3.8;
    cout << x - newtf(x)/newtfl(x);
}
int main() {
    cout << setprecision(5) << fixed;
    Newton();
    return 0;
}

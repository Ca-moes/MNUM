#include <iostream>
#include <cmath>
using namespace std;

double flinha(double x, double y){return pow(x,2)+pow(y,2);}
double f(double x, double y){return (x*pow(y,2) + (pow(x,3)/3.0)); }

void Euler(double lowerlim,double upperlim, double x, double y, double h){
    cout << "EULER\n x     y" << endl;
    int numiter = upperlim/h + h;
    do{
        cout << x << "   " << y << endl;
        y = y + h*flinha(x,y);
        x = x + h;
        numiter--;
    }while(numiter);
}

void RK2(double lowerlim,double upperlim, double x, double y, double h){

    cout << "RK2\n x     y" << endl;
    int numiter = upperlim/h + h;
    do{
        cout << x << "   " << y << endl;
        y = y + h*flinha(x+(h/2.0),y + (h/2.0)*flinha(x,y));
        x = x + h;
        numiter--;
    }while(numiter);
}

void RK4(double lowerlim,double upperlim, double x, double y, double h){

    cout << "RK4\n x     y" << endl;
    int numiter = upperlim/h + h;
    do{
        cout << x << "   " << y << endl;
        y = y + ((1.0/6.0)*(h*flinha(x,y))) + ((1.0/3.0)*(h*flinha(x+(h/2),y+(h*flinha(x,y)/(2.0))))) + ((1.0/3.0) * h * flinha(x+h/2.0,y+(h*flinha(x+(h/2),y+(h*flinha(x,y)/(2.0))))/(2.0))) + (1.0/6.0)*(h*flinha(x+h,y+(flinha(x+h/2.0,y+(h*flinha(x+(h/2),y+(h*flinha(x,y)/(2.0))))/(2.0)))));
        x = x + h;
        if (x == 1) break;
        numiter--;
    }while(numiter);
}

int main() {
    Euler (0,1.5,0,0,0.1);
    cout << "\n\n";
    RK2(0,1.5,0,0,0.1);
    cout << "\n\n";
    RK4(0,1.5,0,0,0.1);
    cout << "Hello, World!" << endl;
    return 0;
}
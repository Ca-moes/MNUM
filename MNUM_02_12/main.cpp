#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

double frety(double x, double y, double z){return x*y*z + pow(x,2)/2.0; }
double fretz(double x, double y, double z){return x*y + pow(x,2)*z/2.0; }
double fretyli(double x, double y, double z){return z*y + x;}
double fretzli(double x, double y, double z){return z*x + y;}

void Euler(double lowerlim,double upperlim, double x, double y, double z, double h){
    vector<double> valy, valz;
    double xini = x, yini = y, zini = z;
    for (int i = 0; i < 2; ++i) {
        cout << "EULER\n h = " << h << "\n x     y       z" << endl;
        double xn, yn, zn;
        int numiter = upperlim/h + h;
        do{
            xn = x + h;
            yn = y + h * fretyli(x,y,z);
            zn = z + h * fretzli(x,y,z);
            cout << xn << "   " << yn << "   " << zn << endl;
            x = xn;
            y = yn;
            z = zn;
            numiter--;
        }while(numiter);

    }




    cout << "EULER\n h = " << h << "\n x     y       z" << endl;
    double xn, yn, zn;
    double xini = x, yini = y, zini = z;
    int numiter = upperlim/h + h;
    do{
        xn = x + h;
        yn = y + h * fretyli(x,y,z);
        zn = z + h * fretzli(x,y,z);
        cout << xn << "   " << yn << "   " << zn << endl;
        x = xn;
        y = yn;
        z = zn;
        numiter--;
    }while(numiter);


    h = h/2;
    x = xini; y = yini; z = zini;
    numiter = upperlim/h + h;
    cout << "EULER\n h = " << h << "\n x      y     z" << endl;
    do{
        xn = x + h;
        yn = y + h * fretyli(x,y,z);
        zn = z + h * fretzli(x,y,z);
        cout << xn << "   " << yn << "   " << zn << endl;
        x = xn;
        y = yn;
        z = zn;
        numiter--;
    }while(numiter);

    h = h/2;
    x = xini; y = yini; z = zini;
    numiter = upperlim/h + h;
    cout << "EULER\n h = " << h << "\n x       y       z" << endl;
    do{
        xn = x + h;
        yn = y + h * fretyli(x,y,z);
        zn = z + h * fretzli(x,y,z);
        cout << xn << "   " << yn << "   " << zn << endl;
        x = xn;
        y = yn;
        z = zn;
        numiter--;
    }while(numiter);
}

/*void RK2(double lowerlim,double upperlim, double x, double y, double h){

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
}*/

int main() {
    Euler(0,0.5,0,1,1,0.05);
    cout << "\n\n";
    /*RK2(0,1.5,0,0,0.1);
    cout << "\n\n";
    RK4(0,1.5,0,0,0.1);
    cout << "Hello, World!" << endl;*/
    return 0;
}
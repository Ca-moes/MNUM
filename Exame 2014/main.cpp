#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
using namespace std;
//Pergunta 2
double f(double x){return -x + 40*cos(sqrt(x))+3;}
double fl(double x){return -(20*sin(sqrt(x)))/(sqrt(x)) - 1;}
double nextx(double x){return x - f(x)/fl(x);}
void resposta2(){
    double x=1.7, gx = f(x);
    cout << x << "\t" << gx << endl;
    x = nextx(x);
    gx = f(x);
    cout << x << "\t" << gx << endl;
    x = nextx(x);
    gx = f(x);
    cout << x << "\t" << gx << endl;
}
//Pergunta 5
double y(double x)
{   //cout << "5+cos(x) = " << 5.0*cos(x) << "\tsin(x) = " << sin(x) << "\t5.0*cos(x) - sin(x) = " << 5.0*cos(x) - sin(x) << endl;
    return (5.0*cos(x)) - sin(x);
}
double B = (sqrt(5)-1)/2.0, A=pow(B,2);
void resposta5(){
    double x1=2, x2=4, x3=A*(x2-x1)+x1, x4 = B*(x2-x1)+x1;
    cout << x1 << "  " << x2 << "  " << x3 << "  " << x4 << "  " << y(x1) << "  " << y(x2) << "  " << y(x3) << "  " << y(x4) << endl;
    if (y(x3) < y(x4)){
        x2=x4;
        x4=x3;
        x3=A*(x2-x1)+x1;
    }
    else{
        x1=x3;
        x3=x4;
        x4 = B*(x2-x1)+x3;
    }
    cout << x1 << "  " << x2 << "  " << x3 << "  " << x4 << "  " << y(x1) << "  " << y(x2) << "  " << y(x3) << "  " << y(x4) << endl;
    if (y(x3) < y(x4)){
        x2=x4;
        x4=x3;
        x3=A*(x2-x1)+x1;
    }
    else{
        x1=x3;
        x3=x4;
        x4 = B*(x2-x1)+x3;
    }
    cout << x1 << "  " << x2 << "  " << x3 << "  " << x4 << "  " << y(x1) << "  " << y(x2) << "  " << y(x3) << "  " << y(x4) << endl;

}
//Pergunta 7
double g(double x){return pow(4*x*x*x - x + 3,1.0/4.0);}
int main() {
    cout << setprecision(5) << fixed;
    cout << "3.5" << endl << g(3.5) << endl << g(g(3.5))<< endl;
    return 0;
}

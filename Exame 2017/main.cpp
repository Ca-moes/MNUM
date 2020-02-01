#include <iostream>
#include <cmath>
using namespace std;

//Pergunta 1
const double B = (sqrt(5)-1.0)/2.0, A = pow(B,2);
double f1(double x){return pow((x-4),2) + pow(x,4);}
void Aurea(){
    double x1 = 0, x2 = 2, x3 = A*(x2-x1)+x1, x4=B*(x2-x1)+x1;
    cout << x1 << "\t" << x2 << endl;
    do {
        if (f1(x3) < f1(x4)) {
            x2 = x4;
            x4 = x3;
            x3 = B * (x4 - x1) + x1;
        } else {
            x1 = x3;
            x3 = x4;
            x4 = B * (x2 - x3) + x3;
        }
        cout << x1 << "\t" << x2 << endl;
    }while(abs(x2-x1)>0.001);
}
// Pergunta 2
double f2(double x){return sqrt(1+pow(2.5,2)*exp(5*x));}
double trapezio(double h){
    double x0 = 0, xn = 1, y0 = f2(x0), yn = f2(xn),sum = 0;
    for (double i = h; i < xn; i = i + h)
        sum = sum + f2(i);
    double result = h/2.0 * (y0+yn+2*sum);
    return result;
}
void Trapezios3times(){
    double h = 0.125, l, ll, lll;
    l = trapezio(h);
    cout << h << "\t" << l << endl;
    ll = trapezio(h/2);
    cout << h/2 << "\t" << ll << endl;
    lll = trapezio(h/4);
    cout << h/4 << "\t" << lll << endl;
    cout << "QC= " << (ll-l)/(lll-ll) << endl;
    cout << "erro= " << (lll-ll)/3.0;
}
double simpson(double h){
    double x0 = 0, x2n = 1, y0 = f2(x0), y2n = f2(x2n),sumpar = 0, sumimpar=0;
    bool impar = true;
    for (double i = h; i < x2n; i = i + h){
        if(impar){
            sumimpar = sumimpar + f2(i);
            impar = false;
        }
        else{
            sumpar = sumpar + f2(i);
            impar = true;
        }
    }

    double result = h/3.0 * (y0+y2n+4*sumimpar+2*sumpar);
    return result;
}
void simpson3times(){
    double h = 0.125, l, ll, lll;
    l = simpson(h);
    cout << h << "\t" << l << endl;
    ll = simpson(h/2);
    cout << h/2 << "\t" << ll << endl;
    lll = simpson(h/4);
    cout << h/4 << "\t" << lll << endl;
    cout << "QC= " << (ll-l)/(lll-ll) << endl;
    cout << "erro= " << (lll-ll)/15.0;
}
//Pergunta 3
// 1. [-6,-4] & [1,3]
// 2. derivadas têm de ser <1 nos intervalos, 1) é <1 no int. negativo, 2) é <1 nos positivo
double xn(double x){return log(5+x);}
double f3(double x){return exp(x)-x-5;}
double f3l(double x){return exp(x)-1;}
void PicardoPeano(){
    double xnext=1, x; int i = 0;
    do{
        x = xnext;
        xnext = xn(x);
        cout << "iteracao " << ++i << "\txnext= " << xnext  << endl;
    }while(abs(xnext-x) > 0.0001);
}
void Newton(){
    double xnext=1, x; int i = 0;
    do{
        x = xnext;
        xnext = x - f3(x)/f3l(x);
        cout << "iteracao " << ++i << "\txnext= " << xnext  << endl;
    }while(abs(xnext-x) > 0.0001);
}
//Pergunta 4
double dC(double C, double T){return -exp(-0.5/(T+273.0)) * C;}
double dT(double C, double T){return 30*exp(-0.5/(T+273.0)) * C - 0.5*(T-20);}
double Euler(double dx, int niter){
    double Tnext, Cnext, xnext, C=2.5, T=25, x=0; int i=0;
    do{
        Tnext = T + dx*dT(C,T);
        Cnext = C + dx*dC(C,T);
        xnext = x+dx;
        cout << "iter " << ++i << "\tx=" << xnext << "\tC=" << Cnext << "\tT=" << Tnext << endl;
        T=Tnext;
        C=Cnext;
        x=xnext;
    }while(i < niter);
    return T;
}
void RK4(){
    double dx=0.25, Tnext, Cnext, xnext, C=2.5, T=25, x=0, dC1, dC2, dC3, dC4, dT1, dT2, dT3, dT4; int i=0;
    do{
        dT1 = dx*dT(C,T);
        dC1 = dx*dC(C,T);

        dT2 = dx*dT(C+(dC1/2.0),T+(dT1/2.0));
        dC2 = dx*dC(C+dC1/2.0,T+dT1/2.0);

        dT3 = dx*dT(C+(dC2/2.0), T+(dT2/2.0));
        dC3 = dx*dC(C+dC2/2.0, T+dT2/2.0);

        dT4 = dx*dT(C+dC3, T+dT3);
        dC4 = dx*dC(C+dC3, T+dT3);

        Tnext = T + dT1/6 + dT2/3 + dT3/3 + dT4/6;
        Cnext = C + dC1/6 + dC2/3 + dC3/3 + dC4/6;
        xnext = x+dx;
        cout << "iter " << ++i << "\tx=" << xnext << "\tC=" << Cnext << "\tT=" << Tnext << endl;
        T=Tnext;
        C=Cnext;
        x=xnext;
    }while(i <2);
}
void Euler3times(){
    double th = Euler(0.25, 2), thl = Euler(0.25/2.0, 4), thll = Euler(0.25/4.0, 8);
    cout << th << "\n" << thl << "\n" << thll << "\n" << (thl-th)/(thll-thl) << "\n" << thll-thl;

}
//Pergunta 5
double wdx(double x, double y){return -1.1*y + 14*x -8;}
double wdy(double x, double y){return -1.1*x + 12;}
void result(){
    double x = 3 - 0.1*wdx(3,1), y=1-0.1*wdy(3,1);
    cout << -1.1*x*y + 12*y + 7*x*x - 8*x;
}
int main() {
    cout.precision(7);
    //Trapezios3times();
    //simpson3times();
    //PicardoPeano();
    //Euler3times();
    return 0;
}

#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

void print(vector< vector<double> > A) {
    int n = A.size();
    for (int i=0; i<n; i++) {
        for (int j=0; j<n+1; j++) {
            cout << A[i][j] << "\t";
            if (j == n-1) {
                cout << "| ";
            }
        }
        cout << "\n";
    }
    cout << endl;
}
vector<double> gauss(vector< vector<double> > A) {
    int n = A.size();

    for (int i=0; i<n; i++) {
        // Search for maximum in this column
        double maxEl = abs(A[i][i]);
        int maxRow = i;
        for (int k=i+1; k<n; k++) {
            if (abs(A[k][i]) > maxEl) {
                maxEl = abs(A[k][i]);
                maxRow = k;
            }
        }

        // Swap maximum row with current row (column by column)
        for (int k=i; k<n+1;k++) {
            double tmp = A[maxRow][k];
            A[maxRow][k] = A[i][k];
            A[i][k] = tmp;
        }

        // Make all rows below this one 0 in current column
        for (int k=i+1; k<n; k++) {
            double c = -A[k][i]/A[i][i];
            for (int j=i; j<n+1; j++) {
                if (i==j) {
                    A[k][j] = 0;
                } else {
                    A[k][j] += c * A[i][j];
                }
            }
        }
    }

    // Solve equation Ax=b for an upper triangular matrix A
    vector<double> x(n);
    for (int i=n-1; i>=0; i--) {
        x[i] = A[i][n]/A[i][i];
        for (int k=i-1;k>=0; k--) {
            A[k][n] -= A[k][i] * x[i];
        }
    }
    return x;
}

// gauss jacobi e seidel
double Gx1(double x2, double x3){return (20+x2-2*x3)/4.0;}
double Gx2(double x1, double x3){return (25-x1-2*x3)/8.0;}
double Gx3(double x1, double x2){return (-10-3*x1+x2)/5.0;}
//trapezios e simpson
double f(double x){return sin(x)/(x*x);}
//trapezios e simpson 3d
double f3(double x, double y){ return exp(y-x);}
// euler & Runga-Kutta's
double dy(double x, double y){return pow(x,2)+pow(y,2);}
double y(double x, double y){return (x*pow(y,2) + (pow(x,3)/3.0)); }

double dz(double t, double z){return 2 + t*t + t*z;}
//euler e runga 3 vars
double dy3vars(double x, double y, double z){return z*y+x;}
double dz3vars(double x, double y, double z){return z*x+y;}

void Gauss_Jacobi(double vix1, double vix2, double vix3){
    //vix1 - valor inicial x1
    cout << "METODO GAUSS-JACOBI\n";
    double x1 = Gx1(vix2,vix3), x2 = Gx2(vix1,vix3), x3 = Gx3(vix1, vix2), nx1, nx2, nx3;
    int i = 0;
    cout << "Iteracao n " << ++i << "  - x1: " << x1 << "   x2: " << x2 << "    x3: " << x3 << endl;
    do{
        nx1 = Gx1(x2,x3); nx2 = Gx2(x1,x3); nx3 = Gx3(x1, x2);
        cout << "Iteracao n " << ++i << "  - x1: " << nx1 << "   x2: " << nx2 << "    x3: " << nx3 << endl;
        if(abs(nx1-x1)<pow(10,-3) && abs(nx2-x2)<pow(10,-3) &&abs(nx3-x3)<pow(10,-3)) break;
        x1 = nx1; x2 = nx2; x3 = nx3;
    }while(true);
}
void Gauss_Seidel(double vix1, double vix2, double vix3){
    //vix1 - valor inicial x1
    cout << "METODO GAUSS-SEIDEL\n";
    double x1 = Gx1(vix2,vix3), x2 = Gx2(x1,vix3), x3 = Gx3(x1, x2), nx1, nx2, nx3;
    int i = 0;
    cout << "Iteracao n " << ++i << "  - x1: " << x1 << "   x2: " << x2 << "    x3: " << x3 << endl;
    do{
        nx1 = Gx1(x2,x3); nx2 = Gx2(nx1,x3); nx3 = Gx3(nx1, nx2);
        cout << "Iteracao n " << ++i << "  - x1: " << nx1 << "   x2: " << nx2 << "    x3: " << nx3 << endl;
        if(abs(nx1-x1)<pow(10,-3) && abs(nx2-x2)<pow(10,-3) &&abs(nx3-x3)<pow(10,-3)) break;
        x1 = nx1; x2 = nx2; x3 = nx3;
    }while(true);
    cout << "\n";
}
double Trapezio(int n, double a, double b){
    double h = (b-a)/n, area = f(a) + f(b);
    for (int i = 1; i < n; i++){ area += 2* f(a + i*h);}
    area *= (h/2);
    cout << "For " << n << " intervals, between " << a << " and " << b << ": area: " << area << endl;
    return area;
}
void Trapezio_QC_ERROR(int n, double a, double b){
    double a1 = Trapezio(n, a, b);
    double a2 = Trapezio(2*n, a, b);
    double a3 = Trapezio(4*n, a, b);
    cout << "QC: " << (a2 - a1) / (a3 - a2) << " ~4" << endl;
    cout << "ERROR: " << (a3 - a2)/ 3.0 << endl;
}
double Simpson(int n, double x0, double xn){
    double h = (xn-x0)/n, sumS= f(x0) + f(xn), solS;
    for (int i = 1; i < n; i++) {
        if (i % 2 == 0)
            sumS += 2 * f(x0 + i * h);
        else
            sumS += 4 * f(x0 + i * h);
    }
    solS = h / 3 * sumS;
    cout << "Simpson - For " << n << " intervals, between " << x0 << " and " << xn << ": area: " << solS << endl;
    return solS;
}
void Simpson_QC_ERROR(int n, double x0, double xn){
    double a1 = Simpson(n, x0, xn);
    double a2 = Simpson(2*n, x0, xn);
    double a3 = Simpson(4*n, x0, xn);
    cout << "QC: " << (a2 - a1) / (a3 - a2) << " ~16" << endl;
    cout << "ERROR: " << (a3 - a2)/ 15.0 << endl;
}
double Trapezio_3(double x0, double xn, double y0, double yn, double nx, double ny){
    double hx = (xn-x0)/nx;
    double hy = (yn-y0)/ny;
    double E0 = f3(x0, y0) + f3(x0, yn) + f3(xn, y0) + f3(xn, yn);
    double E1 = f3(x0 + hx, y0) + f3(x0, y0 + hy) + f3(x0 + hx, yn) + f3(xn, y0 + hy);
    double E2 = f3(x0 + hx, y0 + hy);
    return ((hx * hy) / 4.0) * (E0 + 2*E1 + 4*E2);
}
double Simpson_3(double x0, double xn, double y0, double yn, double nx, double ny){
    double hx = (xn-x0)/nx;
    double hy = (yn-y0)/ny;
    double E0 = f3(x0, y0) + f3(x0, yn) + f3(xn, y0) + f3(xn, yn);
    double E1 = f3(x0 + hx, y0) + f3(x0, y0 + hy) + f3(x0 + hx, yn) + f3(xn, y0 + hy);
    double E2 = f3(x0 + hx, y0 + hy);
    return ((hx * hy) / 9.0) * (E0 + 4*E1 + 16*E2);
}

double Euler(double xini,double xfin, double x0, double y0, double h, double n = 0){
    cout << "EULER        com h = " << h << endl;
    cout << 0 << ": x = " << x0 << "\t\ty = " << y0 << endl;
    double xn = x0,yn=y0;
    if (n==0)
        n = abs((xfin-xini)/h);
    for (int i = 0; i < n; i++) {
        x0 = xn;
        y0 = yn;
        xn = x0 + h;
        yn = y0 + h * dy(x0, y0);
        cout << i +1 << ": x = " << xn << "\t\ty = " << yn << endl;
    }
    return yn;
}
void Euler_QC_ERROR(double xini,double xfin, double x0, double y0, double h, double n = 0){
    double a1 = Euler(xini,xfin,x0,y0,h, n);
    double a2 = Euler(xini,xfin,x0,y0,h/2, n);
    double a3 = Euler(xini,xfin,x0,y0,h/4, n);
    cout << "QC: " << (a2 - a1) / (a3 - a2) << " ~2" << endl;
    cout << "ERROR: " << (a3 - a2) << endl;
}
vector<double> Euler3vars(double xini,double xfin, double x0, double y0, double z0, double h, double n = 0){
    cout << "EULER 3vars       com h = " << h << endl;
    double xn = x0,yn=y0, zn = z0;
    if (n==0)
        n = abs((xfin-xini)/h);
    for (int i = 0; i < n; i++) {
        x0 = xn;
        y0 = yn;
        z0 = zn;
        xn = x0 + h;
        yn = y0 + h * dy3vars(x0, y0, z0);
        zn = z0 + h * dz3vars(x0, y0, z0);
        cout << i +1 << ": x = " << xn << "\t\ty = " << yn << "\t\tz = " << zn << endl;
    }
    vector<double> result = {yn, zn};
    return result;
}
void Euler3vars_QC_ERROR(double xini,double xfin, double x0, double y0, double z0, double h, double n = 0){
    vector<double> a1 = Euler3vars(xini,xfin,x0,y0,z0,h, n);
    vector<double> a2 = Euler3vars(xini,xfin,x0,y0,z0,h/2, n);
    vector<double> a3 = Euler3vars(xini,xfin,x0,y0,z0,h/4, n);
    cout << "a1[0] = " << a1[0] << "\ta2[0] = " << a2[0]<< "\ta3[0] = " << a3[0] << endl;
    cout << "a1[1] = " << a1[1] << "\ta2[1] = " << a2[1]<< "\ta3[1] = " << a3[1] << endl;
    cout << "QCy: " << (a2[0] - a1[0]) / (a3[0] - a2[0]) << " ~2" << endl;
    cout << "ERRORy: " << (a3[0] - a2[0]) << endl;
    cout << "QCz: " << (a2[1] - a1[1]) / (a3[1] - a2[1]) << " ~2" << endl;
    cout << "ERRORz: " << (a3[1] - a2[1]) << endl;
}

double Runge_Kutta2(double xini,double xfin, double x0, double y0, double h, double n = 0){
    cout << "RK2        com h = " << h << endl;
    double xn = x0,yn=y0;
    if (n==0)
        n = abs((xfin-xini)/h);
    for (int i = 0; i < n; i++) {
        x0 = xn;
        y0 = yn;
        xn = x0 + h;
        yn = y0 + h * dy(x0+(h/2.0), y0+(h/2.0)*dy(x0,y0));
        cout << i +1 << ": x = " << xn << "\t\ty = " << yn << endl;
    }
    return yn;
}
void Runge_Kutta2_QC_ERRO(double xini,double xfin, double x0, double y0, double h, double n = 0){
    double a1 = Runge_Kutta2(xini,xfin,x0,y0,h, n);
    double a2 = Runge_Kutta2(xini,xfin,x0,y0,h/2, n);
    double a3 = Runge_Kutta2(xini,xfin,x0,y0,h/4, n);
    cout << "QC: " << (a2 - a1) / (a3 - a2) << " ~4" << endl;
    cout << "ERROR: " << (a3 - a2)/ 3.0 << endl;
}
double Runge_Kutta4(double xini,double xfin, double x0, double y0, double h, double n = 0){
    cout << "RK4        com h = " << h << endl;
    double xn = x0,yn=y0, dy1, dy2, dy3, dy4;
    if (n==0)
        n = abs((xfin-xini)/h);
    for (int i = 0; i < n; i++) {
        x0 = xn;
        y0 = yn;
        xn = x0 + h;
        dy1 = h * dy(x0, y0);
        dy2 = h * dy(x0 + h/2.0, y0 + dy1/2.0);
        dy3 = h * dy(x0 + h/2.0, y0 + dy2/2.0);
        dy4 = h * dy(x0 + h, y0 + dy3);
        yn = y0 + dy1/6.0 + dy2/3.0 + dy3/3.0 + dy4/6.0;
        cout << i +1 << ": x = " << xn << "\t\ty = " << yn << endl;
    }
    return yn;
}
void Runge_Kutta4_QC_ERRO(double xini,double xfin, double x0, double y0, double h, double n = 0){
    double a1 = Runge_Kutta4(xini,xfin,x0,y0,h, n);
    double a2 = Runge_Kutta4(xini,xfin,x0,y0,h/2, n);
    double a3 = Runge_Kutta4(xini,xfin,x0,y0,h/4, n);
    cout << "QC: " << (a2 - a1) / (a3 - a2) << " ~16" << endl;
    cout << "ERROR: " << (a3 - a2)/ 15.0 << endl;
}
vector<double> Runge_Kutta2_3vars(double xini,double xfin, double x0, double y0, double z0, double h, double n = 0){
    cout << "RK2 3vars       com h = " << h << endl;
    double xn = x0,yn=y0,zn = z0;
    if (n==0)
        n = abs((xfin-xini)/h);
    for (int i = 0; i < n; i++) {
        x0 = xn;
        y0 = yn;
        z0 = zn;
        xn = x0 + h;
        yn = y0 + h * dy3vars(x0+(h/2.0), y0+(h/2.0)*dy3vars(x0,y0,z0), z0+(h/2.0)*dz3vars(x0,y0,z0));
        zn = z0 + h * dz3vars(x0+(h/2.0), y0+(h/2.0)*dy3vars(x0,y0,z0), z0+(h/2.0)*dz3vars(x0,y0,z0));
        cout << i +1 << ": x = " << xn << "\t\ty = " << yn << endl;
    }
    vector<double> res = {yn,zn};
    return res;
}
void Runge_Kutta2_3vars_QC_ERRO(double xini,double xfin, double x0, double y0, double z0, double h, double n = 0){
    vector<double> a1 = Runge_Kutta2_3vars(xini,xfin,x0,y0,z0,h, n);
    vector<double> a2 = Runge_Kutta2_3vars(xini,xfin,x0,y0,z0,h/2, n);
    vector<double> a3 = Runge_Kutta2_3vars(xini,xfin,x0,y0,z0,h/4, n);
    cout << "a1[0] = " << a1[0] << "\ta2[0] = " << a2[0]<< "\ta3[0] = " << a3[0] << endl;
    cout << "a1[1] = " << a1[1] << "\ta2[1] = " << a2[1]<< "\ta3[1] = " << a3[1] << endl;
    cout << "QCy: " << (a2[0] - a1[0]) / (a3[0] - a2[0]) << " ~4" << endl;
    cout << "ERRORy: " << (a3[0] - a2[0]) / 3.0 << endl;
    cout << "QCz: " << (a2[1] - a1[1]) / (a3[1] - a2[1]) << " ~4" << endl;
    cout << "ERRORz: " << (a3[1] - a2[1]) / 3.0 << endl;
}
vector<double> Runge_Kutta4_3vars(double xini,double xfin, double x0, double y0, double z0, double h, double n = 0){
    cout << "RK4        com h = " << h << endl;
    double xn = x0,yn=y0,zn=z0, dy1, dy2, dy3, dy4, dz1, dz2, dz3, dz4;
    vector<double> res;
    if (n==0)
        n = abs((xfin-xini)/h);
    for (int i = 0; i < n; i++) {
        x0 = xn;
        y0 = yn;
        z0 = zn;
        xn = x0 + h;
        dy1 = h * dy3vars(x0, y0, z0);
        dz1 = h * dz3vars(x0, y0, z0);
        dy2 = h * dy3vars(x0 + h/2.0, y0 + dy1/2.0, z0 + dz1/2.0);
        dz2 = h * dz3vars(x0 + h/2.0, y0 + dy1/2.0, z0 + dz1/2.0);
        dy3 = h * dy3vars(x0 + h/2.0, y0 + dy2/2.0, z0 + dz2/2.0);
        dz3 = h * dz3vars(x0 + h/2.0, y0 + dy2/2.0, z0 + dz2/2.0);
        dy4 = h * dy3vars(x0 + h, y0 + dy3, z0 + dz3);
        dz4 = h * dz3vars(x0 + h, y0 + dy3, z0 + dz3);
        yn = y0 + dy1/6.0 + dy2/3.0 + dy3/3.0 + dy4/6.0;
        zn = z0 + dz1/6.0 + dz2/3.0 + dz3/3.0 + dz4/6.0;
        cout << i +1 << ": x = " << xn << "\t\ty = " << yn << "\t\tz = " << zn << endl;
    }
    res.push_back(yn);
    res.push_back(zn);
    return res;
}
void Runge_Kutta4_3vars_QC_ERRO(double xini,double xfin, double x0, double y0, double z0, double h, double n = 0){
    vector<double> a1 = Runge_Kutta4_3vars(xini,xfin,x0,y0,z0,h, n);
    vector<double> a2 = Runge_Kutta4_3vars(xini,xfin,x0,y0,z0,h/2, n);
    vector<double> a3 = Runge_Kutta4_3vars(xini,xfin,x0,y0,z0,h/4, n);
    cout << a1.at(0) << " " << a2.at(0) << " " << a3.at(0) << endl;
    cout << a1.at(1) << " " << a2.at(1) << " " << a3.at(1) << endl;
    cout << "QCy: " << (a2.at(0) - a1.at(0)) / (a3.at(0) - a2.at(0)) << " ~16" << endl;
    cout << "ERROR: " << (a3.at(0) - a2.at(0))/ 15.0 << endl;
    cout << "QCz: " << (a2.at(1) - a1.at(1)) / (a3.at(1) - a2.at(1)) << " ~16" << endl;
    cout << "ERROR: " << (a3.at(1) - a2.at(1))/ 15.0 << endl;
}

void EstabExtern2017Anri(){
    vector<vector<double>> matrix = { {18,-1,1},{3,-5,4},{6,8,29} };
    double coefficient_error = 0.1;
    vector<double> sol_list = { 0.552949, -0.15347, -0.10655 };
    vector<double> err_list = { 0, 0, 0 };
    double div = 0, mul = 0;
    for (int i = 0; i < matrix.size(); i++)
    {
        err_list[i] = coefficient_error * (1 - sol_list[0] - sol_list[1] - sol_list[2]);
    }
    for (int i = 0; i < matrix.size(); i++)
    {
        div = matrix[i][i];
        for (int j = i; j < matrix.size(); j++)
        {
            matrix[i][j] /= div;
        }
        err_list[i] /= div;
        for (int j = i + 1; j < matrix.size(); j++)
        {
            mul = matrix[j][i];
            for (int k = i; k < matrix.size(); k++)
            {
                matrix[j][k] -= matrix[i][k] * mul;
            }
            err_list[j] -= err_list[i] * mul;
        }
    }
    for (int i = 0; i < matrix.size(); i++)
    {
        for (int j = 0; j < matrix.size(); j++)
        {
            cout << matrix[i][j] << "\t";
        }
        cout << "|" << err_list[i] << endl;
    }
    for (int i = matrix.size() + 1; i >= 0; i--)
    {
        for (int j = i + 1; j < matrix.size(); j++)
        {
            err_list[i] -= matrix[i][j] * err_list[j];
        }
    }
    cout << "Delta x: " << err_list[0] << endl;
    cout << "Delta y: " << err_list[1] << endl;
    cout << "Delta z: " << err_list[2] << endl;
}
double EulerTest2017(double tini,double tfin, double t0, double y0, double z0, double h, double n = 0){
    cout << "EULERTest        com h = " << h << endl;
    cout << 0 << ": t = " << t0 << "\t\ty = " << y0 <<"\t\tz = " << z0 << endl;
    double tn = t0,yn=y0, zn=z0;
    if (n==0)
        n = abs((tfin-tini)/h);
    for (int i = 0; i < n; i++) {
        t0 = tn;
        z0 = zn;
        y0 = yn;
        tn = t0 + h;
        zn = z0 + h*(2 + t0*t0 + t0*z0);
        yn = y0 + h * z0;
        cout << i +1 << ": t = " << tn << "\t\ty = " << yn << "\t\tz = " << zn << endl;
    }
    return yn;
}
double Gx1t(double x2, double x3, double x4){return -(-25 + 0.5*x2 + 3*x3 + 0.25*x4)/6.0;}
double Gx2t(double x1, double x3, double x4){return -(-10 + 1.2*x1 + 0.25*x3 + 0.2*x4)/3.0;}
double Gx3t(double x1, double x2, double x4){return (7 + 1*x1 - 0.25*x2 - 2*x4)/4.0;}
double Gx4t(double x1, double x2, double x3){return -(12 + 2*x1 + 4*x2 + 1*x3)/8.0;}
void Gauss_SeidelTest2017(double vix1, double vix2, double vix3, double vix4){
    //vix1 - valor inicial x1
    cout << "METODO GAUSS-SEIDELTest\n";
    double x1 = Gx1t(vix2,vix3,vix4), x2 = Gx2t(x1,vix3,vix4), x3 = Gx3t(x1, x2, vix4), x4 = Gx4t(x1,x2,x3), nx1, nx2, nx3, nx4;
    int i = 0;
    cout << "Iteracao n " << ++i << "  - x1: " << x1 << "\tx2: " << x2 << "\tx3: " << x3 << "\tx4: " << x4 << endl;
    do{
        nx1 = Gx1(x2,x3); nx2 = Gx2(nx1,x3); nx3 = Gx3(nx1, nx2);
        cout << "Iteracao n " << ++i << "  - x1: " << nx1 << "   x2: " << nx2 << "    x3: " << nx3 << endl;
        if(abs(nx1-x1)<pow(10,-3) && abs(nx2-x2)<pow(10,-3) &&abs(nx3-x3)<pow(10,-3)) break;
        x1 = nx1; x2 = nx2; x3 = nx3;
    }while(true);
    cout << "\n";
}
void Gauss_Seidel2017Anri(){
    vector<double> result = { 2.12687,2.39858,3.99517,-3.73040 };
    vector<vector<double>> matrix = { {6,0.5,3,0.25}, {1.2,3,0.25,0.2}, {-1,0.25,4,2}, {2,4,1,8} };
    vector<double> sol_list = { 25,10,7,-12 };
    for (int i = 0; i < matrix.size(); i++)
    {
        result[i] = sol_list[i];
        for (int j = 0; j < matrix.at(i).size(); j++)
        {
            if (j != i)
                result[i] -= result[j]*matrix[i][j];
        }
        result[i] /= matrix[i][i];
    }
    cout << "Results: ";
    for (double i : result)
        cout << i << "; ";
}
void Simpson2017Anri(){
    vector<double> results = { 0,0,0 };
    double result, h = 1, err, dx = 0.25;
    int counter;
    vector<double> list = { 1.04,0.37,0.38,1.49,1.08,0.13,0.64,0.84,0.12 };
    for (int j = 0; j < 3; j++)
    {
        result = list[0];
        counter = 1;
        for (int i = h/dx; i < list.size() - h/dx; i += h/dx)
        {
            if (counter % 2 == 0)
                result += list[i] * 2;
            else
                result += list[i] * 4;
            counter++;
        }
        result += list[list.size() - 1];
        result *= h / 3;
        cout << "Step: " << h << "; Result: " << result << endl;
        results[j] = result;
        h /= 2;
    }
    err = (results[2] - results[1])/15;
    cout << "Error: " << err;
}
void CubTrapezios2017Anri(){
    double sum = 0;
    vector<vector<double>> xy = { {1.1,1.4,7.7},{2.1,3.1,2.2},{7.3,1.5,1.2} };
    {
        for (int i = 0; i < xy.size(); i++)
        {
            for (int j = 0; j < xy.size(); j++)
            {
                if ((i == 0 || i == xy.size() - 1) && (j == 0 || j == xy.size() - 1))
                    sum += xy[i][j];
                else if (i == 0 || j == 0 || i == xy.size() - 1 || j == xy.size() - 1)
                    sum += 2 * xy[i][j];
                else
                    sum += 4 * xy[i][j];
            }
        }
    }
    sum /= 4;
    cout << "Sum: " << sum;
}
void Runge_Kutta2017Anri(){
    double h = 0.25;
    double t = 1;
    double y = 1;
    double z = 0;
    double z1, z2, z3, z4, dz1, dz2, dz3, dz4;
    for (int i = 0; i < 3; i++)
    {
        cout << "Iteracao " << i << endl;
        cout << "t: " << t << endl;
        cout << "y: " << y << endl;
        z1 = (h * z);
        dz1 = h * dz(t, z);
        z2 = (h * (z + dz1 / 2));
        dz2 = (h * (dz(t + h/2, z + dz1 / 2)));
        z3 = (h * (z + dz2 / 2));
        dz3 = (h * (dz(t + h/2, z + dz2 / 2)));
        z4 = (h * (z + dz3));
        dz4 = (h * (dz(t + h, z + dz3)));
        y += z1 / 6 + z2 / 3 + z3 / 3 + z4 / 6;
        z += dz1 / 6 + dz2 / 3 + dz3 / 3 + dz4 / 6;
        t += h;
    }

}

double EulerTest2016(double tini,double tfin, double t0, double C0, double T0, double h, double n = 0){
    cout << "EULERTest2016        com h = " << h << endl;
    cout << 0 << ": t = " << t0 << "\t\tC = " << C0 <<"\t\tT = " << T0 << endl;
    double tn = t0, Cn=C0, Tn=T0;
    if (n==0)
        n = abs((tfin-tini)/h);
    for (int i = 0; i < n; i++) {
        t0 = tn;
        T0 = Tn;
        C0 = Cn;
        tn = t0 + h;
        Cn = C0 + h*(-exp(-0.5/(T0+273)) * C0);
        Tn = T0 + h * (20*(exp(-0.5/(T0+273)) * C0) - 0.5 * (T0-20));
        cout << i +1 << ": t = " << tn << "\t\tC = " << Cn << "\tT = " << Tn << endl;
    }
    return Cn;
}
double dC(double T, double C){return -exp(-0.5 / (T + 273))*C;}
double dT(double T, double C){return 20*exp(-0.5 / (T + 273))*C - 0.5*(T - 20);}
void Euler2016Anri(){
    vector<double> sol_list = { 0,0,0 };
    cout << "a) EULER:" << endl;
    double t = 0, C = 1, T = 15, deltaC, deltaT, h = 0.25;
    for (int i = 0; i < 3; i++)
    {
        cout << "Iteracao " << i << endl;
        cout << "t: " << t << endl;
        cout << "C: " << C << endl;
        cout << "T: " << T << endl;
        deltaC = h * dC(T, C);
        deltaT = h * dT(T, C);
        if (i != 2)
        {
            C += deltaC;
            T += deltaT;
            t += h;
        }
    }
    sol_list[0] = C;
    cout << endl;
    cout << "c) QC and Error on Euler:" << endl;
    double qc, err;
    h /= 2;
    for (int i = 0; i < 2; i++)
    {
        t = 0;
        C = 1;
        T = 15;
        for (int j = 0; j < (0.5/h) + 1; j++)
        {
            deltaC = h * dC(T, C);
            deltaT = h * dT(T, C);
            if (j != (0.5/h))
            {
                C += deltaC;
                T += deltaT;
                t += h;
            }
        }
        sol_list[i + 1] = C;
        h /= 2;
    }
    for (int i = 0; i < 3; i++)
    {
        cout << "Iteracao " << i << " : " << sol_list[i] << endl;
    }
    qc = (sol_list[1] - sol_list[0]) / (sol_list[2] - sol_list[1]);
    err = abs(sol_list[2] - sol_list[1]);
    cout << "QC: " << qc << endl;
    cout << "Epsilon: " << err << endl;
}
void RK42016Anri(){
    vector<double> sol_list = { 0,0,0 };
    cout << "a) EULER:" << endl;
    double t = 0, C = 1, T = 15, deltaC, deltaT, h = 0.25;
    cout << "b) RUNGE KUTTA:" << endl;
    t = 0;
    C = 1;
    T = 15;
    double dC1, dC2, dC3, dC4, dT1, dT2, dT3, dT4;
    for (int i = 0; i < 3; i++)
    {
        cout << "Iteracao " << i << endl;
        cout << "t: " << t << endl;
        cout << "C: " << C << endl;
        cout << "T: " << T << endl;
        dC1 = h * dC(T, C);
        dT1 = h * dT(T, C);
        dC2 = h * dC(T + dT1 / 2, C + dC1 / 2);
        dT2 = h * dT(T + dT1 / 2, C + dC1 / 2);
        dC3 = h * dC(T + dT2 / 2, C + dC2 / 2);
        dT3 = h * dT(T + dT2 / 2, C + dC2 / 2);
        dC4 = h * dC(T + dT3, C + dC3);
        dT4 = h * dT(T + dT3, C + dC3);
        C += dC1/6 + dC2/3 + dC3/3 + dC4/6;
        T += dT1/6 + dT2/3 + dT3/3 + dT4/6;
        t += h;
    }
}
void Simpson2016Anri(){
    vector<double> f = { 0.18,0.91,0.83,1.23,0.88,1.37,0.80,1.34,0.43 }, sums = { 0,0,0 };
    double h = 0.8, x = 0, sum = 0, err;
    int counter;
    for (int i = 0; i < 3; i++)
    {
        x = 0;
        counter = 1;
        sum = f[x / 0.2];
        x = counter*h;
        while (x < 1.5)
        {
            if (counter % 2 == 0)
                sum += 2 * f[x / 0.2];
            else
                sum += 4 * f[x / 0.2];
            counter++;
            x = counter*h;
        }
        sum += f[x / 0.2];
        cout << (h/3)*sum << endl;
        sums[i] = (h/3)*sum;
        h /= 2;
    }
    err = (sums[2] - sums[1]) / 15;
    cout << "Error: " << err;
}
void Euler22016Anri(){
    double T = 10, t = 5, h = 0.4;
    for (int i = 0; i < 2; i++)
        T += h*(-0.25*(T - 42));
    cout << T;
}
void Gauss_Sidel2016Anri(){
    double x0 = 0, y0 = 0, z0 = 0, t0 = 0;
    double x = (2.5 - 0.5*y0 - 3 * z0 - 0.25*t0)/6;
    double y = (3.8 - 1.2*x - 0.25 * z0 - 0.2*t0)/3;
    double z = (10 + x - 0.25 * y - 2*t0)/4;
    double t = (7 - 2*x - 4 * y - z)/8;
    cout << "x: " << x << endl << "y: " << y << endl << "z: " << z << endl << "t: " << t << endl;
}

int main() {
    /*Gauss_Jacobi(0,0,0);
    Gauss_Seidel(0,0,0);
    Simpson_QC_ERROR(8,M_PI_2,M_PI);
    cout << Trapezio_3(0,0.5,0,0.5,2,2) << endl;
    cout << Simpson_3(0,0.5,0,0.5,2,2) << endl;
    Euler_QC_ERROR(0,1.4,0,0,0.1);
    Runge_Kutta2_QC_ERRO(0,1.4,0,0,0.1);
    Runge_Kutta4_QC_ERRO(0,1.4,0,0,0.1/8);
    Euler3vars_QC_ERROR(0,0.5,0,1,1,0.05);
    Runge_Kutta2_3vars_QC_ERRO(0,0.5,0,1,1,0.05);
    Runge_Kutta4_3vars_QC_ERRO(0,0.5,0,1,1,0.05);*/

    //EulerTest2017(1,1.5,1,1,0,0.25);
    //Gauss_SeidelTest2017(2.12687, 2.39858, 3.99517, -3.73040);
    //Runge_Kutta4Test2017();
    //Gauss_Seidel2017Anri();
    //EulerTest2016(0,0.5,0,1,15, 0.25);
    //Simpson2017Anri();
    //CubTrapezios2017Anri();
    Euler2016Anri();


    //   Para Gauss - Usar com eliminação de Guass e Khaletsy
    // n = numero de variaveis
    /*int n = 3;
    vector<double> line(n+1,0);
    vector< vector<double> > A(n,line);
    // A:
    A[0][0] = 4;
    A[0][1] = -1;
    A[0][2] = 2;
    A[1][0] = 1;
    A[1][1] = 8;
    A[1][2] = 2;
    A[2][0] = 3;
    A[2][1] = -1;
    A[2][2] = 5;
    // B:
    A[0][3] = 20;
    A[1][3] = 25;
    A[2][3] = -10;
    // Print input
    print(A);
    // Calculate solution
    vector<double> x(n);
    x = gauss(A);
    // Print result
    cout << "Result:\t";
    for (int i=0; i<n; i++) {
        cout << x[i] << " ";
    }
    cout << endl;*/
}
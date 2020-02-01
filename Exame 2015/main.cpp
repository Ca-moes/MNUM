#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
using namespace std;
typedef vector<vector<double>> Matrix;

//Pergunta 1
double dTdt(double T){return -0.25*(T-37);}
void Euler(){
    double T=3, t=5, Ta=37, h=0.4, tnext, Tnext;
    Tnext = T + dTdt(T)*h;
    tnext = t + h;
    cout << "T: " << Tnext << "\tt: " << tnext << endl;
    T = Tnext;
    t = tnext;
    Tnext = T + dTdt(T)*h;
    tnext = t + h;
    cout << "T: " << Tnext << "\tt: " << tnext << endl;
}
//Pergunta 3
Matrix solveMatrix(Matrix matrix){
    double aux;
    for (int i = 0; i < matrix.size(); i++) {
        aux = matrix[i][i];

        for (int j = 0; j < matrix[i].size(); j++) {
            matrix[i][j] /= aux;
        }

        for (int j = i+1; j < matrix.size(); j++) {
            aux = matrix[j][i];
            for (int k = i; k < matrix[j].size(); k++) {
                matrix[j][k] -= matrix[i][k] * aux;
            }
        }
    }
    return matrix;
}
void printMatrix(Matrix matrix){
    for (int i = 0; i < matrix.size(); i++) {
        for (int j = 0; j < matrix[i].size(); ++j) {
            cout << matrix[i][j] << "\t";
        }
        cout << endl;
    }
}
vector<double> resultsMatrix(Matrix matrix){
    vector<double> results = { matrix[0][matrix[0].size()-1], matrix[1][matrix[0].size()-1], matrix[2][matrix[0].size()-1]};
    for (int i = matrix.size() - 1; i >= 0 ; i--) {
        for (int j = i+1; j < matrix.size(); ++j) {
            results[i] -= matrix[i][j] * results[j];
        }
    }
    return results;
}
void answer(){
    Matrix matrix = { {1.0,1.0 / 2.0,1.0 / 3.0,-1.0},{1.0/2.0,1.0 / 3.0,1.0 / 4.0,1.0},{1.0 / 3.0,1.0 / 4.0,1.0 / 5.0,1.0} };
    Matrix solvedmatrix = solveMatrix(matrix);

    printMatrix(solvedmatrix);
    vector<double> results = resultsMatrix(solvedmatrix);

    for (int i = 0; i < results.size(); i++) {
        cout << "Resultado x" << i+1 << " = " << results[i]<<endl;
    }

    double residuo = 0.05;
    vector<double> resi = {0,0,0};

    for (int i = 0; i < resi.size(); i++) {
        resi[i] = residuo;
        for (int j = 0; j < results.size(); j++) {
            resi[i] -= residuo * results[j];
        }
    }

    for (int i = 0; i < 3; i++) {
        matrix[i][matrix[i].size()-1]=resi[i];
    }

    printMatrix(solveMatrix(matrix));

    results = resultsMatrix(solveMatrix(matrix));

    for (int i = 0; i < results.size(); i++) {
        cout << "Resultado x" << i+1 << " = " << results[i] << endl;
    }
}
//Pergunta 4
double xnext(double x){return 2*log(2*x);}
void thiss(){cout << xnext(1.1) << endl << xnext(1.1)-1.1;}
//Pergunta 5
double y(double x){return sqrt(1+pow(2.5*exp(x*2.5),2));}
double trap(double h){
    double a=0, b=1, y0=y(a), yn = y(b), sum=0;
    for (double i = h; i < b; i=i+h) {
        sum = sum + y(i);
    }
    return h/2.0 * (y0+yn+2*sum);
}
void trap3(){
    double h=0.125;
    cout << h << endl;
    cout << h/2 << endl;
    cout << h/4 << endl;
    cout << trap(h) << endl;
    cout << trap(h/2) << endl;
    cout << trap(h/4) << endl;
    cout << (trap(h/2) - trap(h))/(trap(h/4) - trap(h/2))<< endl;
    cout << (trap(h/4)-trap(h/2))/3.0<< endl;
}
double simp(double h){
    double a=0, b=1, y0=y(a), yn = y(b), sumpar=0, sumimpar=0; bool impar=true;
    for (double i = h; i < b; i=i+h) {
        if (impar){
            sumimpar = sumimpar + y(i);
            impar = false;
        } else{
            sumpar = sumpar + y(i);
            impar = true;
        }

    }
    return h/3.0 * (y0+yn+4*sumimpar+2*sumpar);
}
void simp3(){
    double h=0.125;
    cout << h << endl;
    cout << h/2 << endl;
    cout << h/4 << endl;
    cout << simp(h) << endl;
    cout << simp(h/2) << endl;
    cout << simp(h/4) << endl;
    cout << (simp(h/2) - simp(h))/(simp(h/4) - simp(h/2)) << " ~ 16" << endl;
    cout << (simp(h/4)-simp(h/2))/15.0<< endl;
}
int main() {
    cout << setprecision(5) << fixed;
    //
    answer();
}

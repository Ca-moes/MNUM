#include <iostream>
#include <iomanip>
#include <math.h>
using namespace std;

int main() {
    double anterior = exp(1) - 1;
    for (int i = 0; i < 26; ++i) {
        cout <<setw(2) << i << " - " << setprecision(10) << anterior << endl;
        //cout << "anterior : " << anterior << "       i : " << i << "         - 1.0";
        anterior = anterior * (double(i) +1) - 1.0;
    }
    return 0;
}
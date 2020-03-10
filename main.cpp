//
//  main.cpp
//  SLAU. Gauss
//
//  Created by Паршиков Мирон on 29.03.2018.
//  Copyright © 2018 Паршиков Мирон. All rights reserved.
//

#include <iostream>
using namespace std;
class matrix{
    int m;
    int n;
    double **A;
    double *B;
protected:
    double sum_for_line (const double* line){
        double sum;
        return sum;
    }
public:
    matrix (){
        cout << "m = ";
        cin >> m;
        cout << "n = ";
        cin >> n;
        A = new double*[m];
        B = new double[m];
        for (int i = 0; i < m; i++) {
            A[i] = new double[n];
        }
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n+1; j++) {
                if (j != n) {
                    printf("a[%d][%d] = ", i+1, j+1);
                    cin >> A[i][j];
                }
                if (j == n) {
                    printf("b[%d] = ", i+1);
                    cin >> B[i];
                }
            }
        }
    }
    matrix (int m1, int n1){
        m = m1;
        n = n1;
        A = new double*[m];
        for (int i = 0; i < m; i++) {
            A[i] = new double[n];
        }
    }
    ~matrix (){
        for (int i = 0; i < m; i++) {
            delete [] A[i];
        }
    }
    int det () const {
        //определитель матрицы будем раскладывать по 1-ой строке => в индексации Си - по нулевой
        int d = 0;
        if (n == 1) {
            d = A[0][0];
            return d;
        }
        else {
            for (int j = 0; j < n; j++) {
                //создаем дополнительный минор на каждом шаге, пока матрица не будет иметь размер 1х1
                matrix minor(n-1, n-1);
                for (int i1 = 0+1; i1 < n; i1++) {
                    for (int j1 = 0; j1 < j; j1++) {
                        minor.A[i1-1][j1] = A[i1][j1];
                    }
                    for (int j1 = j+1; j1 < n; j1++) {
                        minor.A[i1-1][j1-1] = A[i1][j1];
                    }
                }
                if ((0 + j) % 2 == 0){
                    d = d + A[0][j]*minor.det(); //
                }
                else{
                    d = d - A[0][j]*minor.det();
                }
            }
            return d;
        }
    }
    
    void gauss_for_slau(){
        if (det() != 0) {
            for (int i = 0; i < n-1; i++) {
                if (A[i][i] == 0) {
                    for (int j = 1; j < n; j++) {
                        if (A[j][i] != 0) {
                            int help;
                            for (int k = 0; k < n; k++) {
                                help = A[j][k];
                                A[j][k] = A[i][k];
                                A[i][k] = help;
                            }
                            break;
                        }
                    }
                }
                for (int k = i+1; k < n; k++) {
                    int help = A[k][i];
                    for (int j = 0; j < n; j++) {
                        A[k][j] = A[k][j] - (A[i][j]/A[i][i])*help;
                    }
                    B[k] = B[k] - (B[i]/A[i][i])*help;
                }
            }
            for (int i = 0; i < n; i++) {
                double sum = 0;
                for (int j = n-i; j < n; j++) {
                    sum = sum + A[n-i-1][j]*B[j];
                }
                B[n-i-1]=(B[n-i-1]-sum)/A[n-i-1][n-i-1];
            }
        }
    }
    void output(){
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                cout << A[i][j] << " ";
            }
            cout << B[i];
            cout << endl;
        }
    }
};
int main(int argc, const char * argv[]) {
    matrix A;
    printf("Det = %d\n", A.det());
    A.gauss_for_slau();
    A.output();
    return 0;
}

#include <iostream>
#include <cmath>
#include <fstream>

double f(double x){
    return sin(x) + log(x);
}

void addToFile(double *x, double *y){
    std::ofstream file("input.txt");
    for (int i = 0; i < 10; ++i) {
        file << x[i] << " " << y[i] << "\n";
    }
    file.close();
}

void addFuncToFile(){
    double x[11]{1,2,3,4,5,6,7,8,9,10};
    double y[11];
    for (int i = 1; i <= 10; ++i) {
        y[i] = f(x[i]);
    }
    addToFile(x,y);
    return;
}

double *getMatrix(double *x, int n){
    double *xy = new double [n * n];
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            xy[i * n + j] = pow(x[i], j);
        }
    }
    return xy;
}

int readFromFile(double *x, double *y){
    std::ifstream file("input.txt");
    if (!file){
      return 1;
    }
    for (int i = 0; i < 10; ++i) {
        file >> x[i];
        file >> y[i];
    }
    file.close();
    return 0;
}
double* solveGaussian(double* A, double* y, int n) {
    double* x = new double[n];

    for (int k = 0; k < n; ++k) {
        int maxRow = k;
        for (int i = k + 1; i < n; ++i) {
            if (std::fabs(A[i * n + k]) > std::fabs(A[maxRow * n + k])) {
                maxRow = i;
            }
        }

        for (int j = 0; j < n; ++j) {
            std::swap(A[k * n + j], A[maxRow * n + j]);
        }
        std::swap(y[k], y[maxRow]);

        double pivot = A[k * n + k];
        if (std::fabs(pivot) < 1e-9) {
            throw std::runtime_error("ERROR");
        }
        for (int j = k; j < n; ++j) {
            A[k * n + j] /= pivot;
        }
        y[k] /= pivot;

        for (int i = k + 1; i < n; ++i) {
            double factor = A[i * n + k];
            for (int j = k; j < n; ++j) {
                A[i * n + j] -= factor * A[k * n + j];
            }
            y[i] -= factor * y[k];
        }
    }

    for (int i = n - 1; i >= 0; --i) {
        x[i] = y[i];
        for (int j = i + 1; j < n; ++j) {
            x[i] -= A[i * n + j] * x[j];
        }
    }

    return x;
}


double LX(double *a, double x){
    double res = 0;
    for (int i = 0; i < 10; ++i) {
        res += a[i] * pow(x,i);
    }
    return res;
}



int main(){
//    addFuncToFile();
    double *x = new double[10];
    double *y = new double [10];
    if (readFromFile(x,y)){
        return 1;
    }
    double *m =  getMatrix(x, 10);
    for (int i = 0; i < 10; ++i) {
        for (int j = 0; j < 10; ++j) {
            std::cout << m[i * 10 + j] << " ";
        }
        std::cout << y[i] << "\n";
    }
    double *a = solveGaussian(m,y,10);
    for (int i = 0 ; i < 10; i++) std::cout<<a[i] << " ";
    std::ofstream
    for (double temp = 0.; temp <= 3; temp += 0.1){

    }
//    for (int i = 0; i < 9;i++) {
//        std::cout << i * 10 + i << "\t";
//        for (int j = i + 1; j < 10; ++j) {
//            for (int k = 0; k < 10; ++k) {
//                std::cout << k;
//            }
//            std::cout << j << " ";
//        }
//        std::cout<< "\n";
////        i+=10+i;
////        double k = m[i*10] / m[i*10+10];
////        for (int j = 0; j < 9; ++j) {
////            std::cout << i << " " << j << "\n";
////            m[(i+1)*10+j] -= k * m[i*10+j];
////        }
////        y[i+1] -= k * y[i];
//    }
////    for (int i = 0; i < 10; ++i) {
////        for (int j = 0; j < 10; ++j) {
////            std::cout << m[i * 10 + j] << " ";
////        }
////        std::cout << y[i] << "\n";
////    }


    


}


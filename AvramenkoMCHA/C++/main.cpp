#include <iostream>
#include "gauss.h"

int main()
{
    int k = 1;
    std::vector<double> b = {4.2, 4.2, 4.2, 4.2, 4.2};

    std::vector<std::vector<double>> C = {
        {0.2, 0.0, 0.2, 0.0, 0.0},
        {0.0, 0.2, 0.0, 0.2, 0.0},
        {0.2, 0.0, 0.2, 0.0, 0.2},
        {0.0, 0.2, 0.0, 0.2, 0.0},
        {0.0, 0.0, 0.2, 0.0, 0.2}};
    std::vector<std::vector<double>> D = {
        {2.33, 0.81, 0.67, 0.92, -0.53},
        {-0.53, 2.33, 0.81, 0.67, 0.92},
        {0.92, -0.53, 2.33, 0.81, 0.67},
        {0.67, 0.92, -0.53, 2.33, 0.81},
        {0.81, 0.67, 0.92, -0.53, 2.33}};

    std::vector<std::vector<double>> A = sumMatrix(multiplyMatrixByNumber(k, C), D);
    std::vector<double> resGauss = solveGauss(A, b);
    std::vector<double> resChooseColumn = solveChooseByColumn(A, b);
    std::vector<double> resChooseMatrix = solveChooseByMatrix(A, b);

    std::cout << "x1: ";
    for (int i = 0; i < 5; i++)
    {
        std::cout << resGauss[i] << "    " ;
    }
    std::cout << "\n\nx2: ";
    for (int i = 0; i < 5; i++)
    {
        std::cout << resChooseColumn[i] << "    "  ;
    }
    std::cout << "\n\nx3: ";
    for (int i = 0; i < 5; i++)
    {
        std::cout << resChooseMatrix[i] << "    "  ;
    }
    
    std::vector<double> b_g = checkSolution(A, resGauss);
    std::vector<double> b_cc = checkSolution(A, resChooseColumn);
    std::vector<double> b_cm = checkSolution(A, resChooseMatrix);

    std::vector<double> resid1 = checkSolution(A, resGauss);
    std::vector<double> resid2 = checkSolution(A, resChooseColumn);
    std::vector<double> resid3 = checkSolution(A, resChooseMatrix);

    std::cout << "\n\nresid1: ";
    for (int i = 0; i < 5; i++)
    {
        std::cout << fabs(b_g[i] - b[i]) << "    " ;
        resid1.push_back(fabs(b_g[i] - b[i]));
    }
    std::cout << "\n\nresid2: ";
    for (int i = 0; i < 5; i++)
    {
        std::cout << fabs(b_cc[i] - b[i])<< "    " ;
        resid2.push_back(fabs(b_cc[i] - b[i]));
    }
    std::cout << "\n\nresid3: ";
    for (int i = 0; i < 5; i++)
    {
        std::cout << fabs(b_cm[i] - b[i]) << "    " ;
        resid3.push_back(fabs(b_cm[i] - b[i]));
    }

    std::cout << "\n\nnorm1: ";
    std::cout << norm(resid1);
    std::cout << "\n\nnorm2: ";
    std::cout << norm(resid2);
    std::cout << "\n\nnorm3: ";
    std::cout << norm(resid3);
}
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

    for (int i = 0; i < 5; i++)
    {
        std::cout << resGauss[i] << '\t' << resChooseColumn[i] << '\t' << resChooseMatrix[i];
        std::cout << '\n';
    }
}
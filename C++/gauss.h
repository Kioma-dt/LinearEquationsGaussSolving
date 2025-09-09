#pragma once
#include <vector>
#include <cmath>
#include <stdexcept>

std::vector<std::vector<double>> sumMatrix(std::vector<std::vector<double>> A, std::vector<std::vector<double>> B);
std::vector<std::vector<double>> multiplyMatrixByNumber(double num, std::vector<std::vector<double>> A);

std::vector<double> solveGauss(std::vector<std::vector<double>> A, std::vector<double> b);

std::vector<double> solveChooseByColumn(std::vector<std::vector<double>> A, std::vector<double> b);

std::vector<double> solveChooseByMatrix(std::vector<std::vector<double>> A, std::vector<double> b);

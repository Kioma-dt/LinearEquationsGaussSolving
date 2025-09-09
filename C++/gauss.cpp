#include "gauss.h"

std::vector<std::vector<double>> sumMatrix(std::vector<std::vector<double>> A, std::vector<std::vector<double>> B)
{
    std::vector<std::vector<double>> C;

    for (int i = 0; i < A.size(); i++)
    {
        std::vector<double> row;
        for (int j = 0; j < A[i].size(); j++)
        {
            row.push_back(A[i][j] + B[i][j]);
        }
        C.push_back(row);
    }

    return C;
}

std::vector<std::vector<double>> multiplyMatrixByNumber(double num, std::vector<std::vector<double>> A)
{
    std::vector<std::vector<double>> C;

    for (int i = 0; i < A.size(); i++)
    {
        std::vector<double> row;
        for (int j = 0; j < A[i].size(); j++)
        {
            row.push_back(num * A[i][j]);
        }
        C.push_back(row);
    }

    return C;
}

void nullColumn(std::vector<std::vector<double>> &A, std::vector<double> &b, int i)
{
    for (int j = i + 1; j < A.size(); j++)
    {
        double q = A[j][i] / A[i][i];

        for (int k = i; k < A.size(); k++)
        {
            A[j][k] -= q * A[i][k];
        }
        b[j] -= q * b[i];
    }
}

void swapRows(std::vector<std::vector<double>> &A, std::vector<double> &b, int i, int j)
{
    A[i].swap(A[j]);
    std::swap(b[i], b[j]);
}

void swapColumns(std::vector<std::vector<double>> &A, int i, int j, std::vector<double> &variablesSubstitution)
{
    for (int k = 0; k < A.size(); k++)
    {
        std::swap(A[k][i], A[k][j]);
    }

    std::swap(variablesSubstitution[i], variablesSubstitution[j]);
}

std::pair<int, int> findNotZeroCoefficient(std::vector<std::vector<double>> A, int i)
{
    for (int j = i; j < A.size(); j++)
    {
        for (int k = i; k < A.size(); k++)
        {
            if (A[j][k] != 0)
            {
                return std::make_pair(j, k);
            }
        }
    }

    return std::make_pair(-1, -1);
}

void placeNotZeroCoefficient(std::vector<std::vector<double>> &A, std::vector<double> &b, int i, std::vector<double> &variablesSubstitution)
{
    auto [firstIndex, secondIndex] = findNotZeroCoefficient(A, i);

    if (firstIndex < 0)
    {
        for (int j = i; j < b.size(); j++)
        {
            if (b[j] != 0)
            {
                throw std::runtime_error("No solutions");
            }
        }

        throw std::runtime_error("Infinite number of solutions");
    }

    if (firstIndex != i)
    {
        swapRows(A, b, i, firstIndex);
    }

    if (secondIndex != i)
    {
        swapColumns(A, i, secondIndex, variablesSubstitution);
    }
}

std::vector<double> backSubstitution(std::vector<std::vector<double>> A, std::vector<double> b, std::vector<double> variablesSubstitution)
{
    std::vector<double> res(A.size());
    for (int i = A.size() - 1; i >= 0; i--)
    {
        double temp = 0;
        for (int j = i + 1; j < A.size(); j++)
        {
            temp += A[i][j] * res[variablesSubstitution[j]];
        }
        res[variablesSubstitution[i]] = ((b[i] - temp) / A[i][i]);
    }
    return res;
}

std::vector<double> solveGauss(std::vector<std::vector<double>> A, std::vector<double> b)
{
    std::vector<double> variablesSubstitution{0, 1, 2, 3, 4};
    for (int i = 0; i < A.size() - 1; i++)
    {
        if (A[i][i] == 0)
        {
            placeNotZeroCoefficient(A, b, i, variablesSubstitution);
        }

        nullColumn(A, b, i);
    }

    auto res = backSubstitution(A, b, variablesSubstitution);
    return res;
};

std::vector<double> solveChooseByColumn(std::vector<std::vector<double>> A, std::vector<double> b)
{
    std::vector<double> variablesSubstitution{0, 1, 2, 3, 4};
    for (int i = 0; i < A.size() - 1; i++)
    {
        double maxElement = fabs(A[i][i]);
        int maxInColumnIndex = i;

        for (int j = i + 1; j < A.size(); j++)
        {
            if (fabs(A[j][i]) > maxElement)
            {
                maxInColumnIndex = j;
                maxElement = fabs(A[j][i]);
            }
        }

        if (maxInColumnIndex != i)
        {
            swapRows(A, b, i, maxInColumnIndex);
        }

        if (A[i][i] == 0)
        {
            placeNotZeroCoefficient(A, b, i, variablesSubstitution);
        }

        nullColumn(A, b, i);
    }

    auto res = backSubstitution(A, b, variablesSubstitution);
    return res;
}

std::vector<double> solveChooseByMatrix(std::vector<std::vector<double>> A, std::vector<double> b)
{
    std::vector<double> variablesSubstitution{0, 1, 2, 3, 4};
    for (int i = 0; i < A.size() - 1; i++)
    {
        double maxElement = fabs(A[i][i]);
        int maxElementFirstIndex = i;
        int maxElementSecondIndex = i;

        for (int j = i; j < A.size(); j++)
        {
            for (int k = i; k < A.size(); k++)
            {
                if (fabs(A[j][k]) > maxElement)
                {
                    maxElement = fabs(A[j][k]);
                    maxElementFirstIndex = j;
                    maxElementSecondIndex = k;
                }
            }
        }

        if (maxElementFirstIndex != i)
        {
            swapRows(A, b, i, maxElementFirstIndex);
        }

        if (maxElementSecondIndex != i)
        {
            swapColumns(A, i, maxElementSecondIndex, variablesSubstitution);
        }

        if (A[i][i] == 0)
        {
            for (int j = i; j < b.size(); j++)
            {
                if (b[j] != 0)
                {
                    throw std::runtime_error("No solutions");
                }
            }
            throw std::runtime_error("Infinite number of solutions");
        }

        nullColumn(A, b, i);
    }

    auto res = backSubstitution(A, b, variablesSubstitution);
    return res;
}
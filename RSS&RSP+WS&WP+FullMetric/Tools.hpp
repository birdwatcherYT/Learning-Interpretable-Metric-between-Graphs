#ifndef TOOLS_HPP
#define TOOLS_HPP

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <Eigen/Dense>
#include <numeric>
//#include <SymEigsSolver.h>

using namespace std;
using namespace Eigen;

pair<MatrixXd, VectorXi> dataLoad(const string &filename);
VectorXd projectPositive(const VectorXd &v);
pair<double, VectorXd> smallestEigh(const MatrixXd& A, double eps = 1e-11);
pair<double, VectorXd> smallestEigh(const MatrixXd& A, const VectorXd& _x, double eps = 1e-11);
MatrixXd projectSemidefinite(const MatrixXd &M);
MatrixXd projectOneNegative(const MatrixXd &M);
vector<size_t> randomIndex(size_t n);

template <typename T>
vector<size_t> argsort(const vector<T> &v) {
    // initialize original index locations
    vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0);
    // sort indexes based on comparing values in v
    sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
    return idx;
}

#endif /* TOOLS_HPP */

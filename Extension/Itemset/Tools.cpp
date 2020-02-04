#include "Tools.hpp"

pair<MatrixXd, VectorXi> dataLoad(const string &filename) {
    ifstream ifs(filename);
    if (ifs.fail()) {
        cerr << "File open failed." << endl;
        exit(1);
    }
    vector<int> y_temp;
    vector< vector<double> > X_temp;
    string line;
    while (getline(ifs, line)) {
        stringstream ss(line);
        // y
        int label;
        ss >> label;
        y_temp.push_back(label);
        // X
        double value;
        vector<double> row;
        while (ss >> value)
            row.push_back(value);
        X_temp.push_back(row);
    }
    Index n = X_temp.size(), dim = X_temp[0].size();
    MatrixXd X(n, dim);
    VectorXi y(n);
    for (Index i = 0; i < n; ++i) {
        y.coeffRef(i) = y_temp[i];
        for (Index j = 0; j < dim; ++j)
            X.coeffRef(i, j) = X_temp[i][j];
    }
    return pair<MatrixXd, VectorXi>(X, y);
}

VectorXd projectPositive(const VectorXd &v) {
    return v.array().max(0);
}

MatrixXd projectSemidefinite(const MatrixXd &M) {
    SelfAdjointEigenSolver<MatrixXd> eig(M);
    const MatrixXd& V = eig.eigenvectors();
    VectorXd lams = eig.eigenvalues();
    bool negative = false;
    Index n = lams.size();
    for (Index i = 0; i < n; ++i)
        if (lams.coeff(i) < 0)
            lams.coeffRef(i) = 0, negative = true;
    return negative ? (V * lams.asDiagonal() * V.transpose()) : M;
}

vector<size_t> randomIndex(size_t n){
    vector<size_t> idx(n);
    iota(idx.begin(), idx.end(), 0);
    random_shuffle(idx.begin(), idx.end());
    return idx;
}

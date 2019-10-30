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

pair<double, VectorXd> smallestEigh(const MatrixXd& A, double eps) {
    return smallestEigh(A, VectorXd::Random(A.rows()), eps);
}

pair<double, VectorXd> smallestEigh(const MatrixXd& A, const VectorXd& _x, double eps) {
    double px, pp, pAx, pAp, xAx, a, b, c, lambda, lambda_prev, alpha, gSquaredNorm, gSquaredNorm_prev;
    VectorXd x = _x;
    x.normalize();
    VectorXd Ax = A * x;
    lambda = xAx = x.dot(Ax);
    VectorXd g = (Ax - lambda * x), g_prev;
    VectorXd p = -g;
    gSquaredNorm = g.squaredNorm();
    while (true){
        lambda_prev = lambda, gSquaredNorm_prev = gSquaredNorm, g_prev = g;
        pp = p.squaredNorm(), px = p.dot(x), pAx = p.dot(Ax), pAp = p.dot(A * p);
        a = px * pAp - pp * pAx, b = pAp - pp * xAx, c = pAx - px * xAx;
        alpha = (b > 0) ? ( -2 * c / (b + sqrt(b * b - 4 * a * c)) ) : ( (-b + sqrt(b * b - 4 * a * c)) / (2 * a) );
        x += alpha * p;
        x.normalize();
        Ax = A * x;
        lambda = xAx = x.dot(Ax);
        g = (Ax - lambda * x);
        gSquaredNorm = g.squaredNorm();
        if (gSquaredNorm <= eps || fabs(lambda - lambda_prev) <= eps * fabs(lambda))
            return pair<double, VectorXd>(lambda, x);
        if (fabs(g.dot(g_prev)) >= 0.2 * gSquaredNorm)
            p = -g;
        else
            p = -g + (gSquaredNorm / gSquaredNorm_prev) * p;
    }
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

MatrixXd projectOneNegative(const MatrixXd &M) {
    const pair<double, VectorXd>& lam_v = smallestEigh(M);
    const double& lam = lam_v.first;
    const VectorXd& v = lam_v.second;
    /*
    DenseGenMatProd<double> _M(M);
    SymEigsSolver<double, SMALLEST_ALGE, DenseGenMatProd<double> > eig(&_M, 1, M.rows());
    eig.init();
    eig.compute();
    lam = eig.eigenvalues()(0);
    v = eig.eigenvectors().col(0);
     */
    return (lam < 0) ? (M - lam * (v * v.transpose())) : M;
}

vector<size_t> randomIndex(size_t n){
    vector<size_t> idx(n);
    iota(idx.begin(), idx.end(), 0);
    random_shuffle(idx.begin(), idx.end());
    return idx;
}

#include "SequenceMetric.hpp"
using namespace std;

const double L = 1;
const double U = 0.95;
const double eta=1;
const int splitNum = 100;
const int maxloop = INT_MAX;
const double eps = 1e-6;
const int freq = 1;
const double R=pow(0.01, 1.0/(splitNum-1));

string filename="./data/";
const int minsup=1;
int maxpat=0;
const double train=0.6;
const int Ksample=10;
string kernelPath="../sequenceKernel/";

int main(int argc, char** argv) {
    if (argc==4){
        srand(atoi(argv[1]));
        filename+=string(argv[2]);
        maxpat=atoi(argv[3]);
    }else{
        cerr<<"Usage: ./run srand filename maxpat"<<endl;
        exit(1);
    }
    SequenceMetric m(filename, minsup, maxpat, train, Ksample);

    // string name=string(argv[2]);
    // double validMaxF=0;
    // string argMax="";
    // for (int i=3;i<=5;++i){
        // cout<<"kernel"<<i<<endl;
        // string path=kernelPath+name+to_string(i);
        // double f_measure=m.kernel(path);
        // if (validMaxF<f_measure)
            // validMaxF=f_measure, argMax=path;
    // }
    // m.regularizationPath(L, U, eta, splitNum, maxloop, eps, freq, R, argMax);
    // return 0;

    m.regularizationPath(L, U, eta, splitNum, maxloop, eps, freq, R);
    return 0;
}

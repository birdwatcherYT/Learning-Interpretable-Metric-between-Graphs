#ifndef ITEMSETMETRIC_HPP
#define ITEMSETMETRIC_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <map>
#include <vector>
#include <cmath>
#include <ctime>
#include <climits>
#include <cfloat>
#include "tree.hh"
#include "Tools.hpp"

typedef std::vector<int> Pattern;
typedef std::vector< std::pair<int, int> > ProjectDB;
typedef std::vector< std::vector<bool> > vecVec;
typedef std::vector< std::vector<int> > Transact;
typedef std::vector< std::pair< Pattern, int> > PatternSupportList;
typedef std::vector< Pattern> PatternList;

std::ostream& operator<<(std::ostream& os, const PatternSupportList& patternSupportList);

template<class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T> vec) {
    os << "[ ";
    for (const T &e : vec)
        os << e << ", ";
    os << "]";
    return os;
}


// this が 引数を含むかどうか
// 含む:1, 同じ:0, 含まない:-1
int contain(const Pattern &pattern1, const Pattern &pattern2);

class ItemsetMetric {
private:
    Transact TRAIN;
    Transact TEST;
    int minsup;
    int maxpat;
    int Ksample;
    std::vector<int> y_train;
    std::vector<int> y_test;
    
    int trainSize;
    int testSize;

    std::set<int> classLabel;
    std::vector< std::vector<int> > sameClass;
    std::vector< std::vector<int> > diffClass;
    std::vector<size_t> trainIndex;
    std::vector<size_t> testIndex;

    class Save {
    public:
        std::vector<bool> x;
        Pattern pattern;
        ProjectDB pdb;
        bool nextCheck=false;
        double screeningScore=DBL_MAX;
        double pruningScore=DBL_MAX;
        Save(){}
        Save(const std::vector<bool> &x, const Pattern &pattern, const ProjectDB &pdb){
            this->x=x, this->pattern=pattern, this->pdb=pdb;
        }
    };
    tree<Save> Tree;
    int visit;

    tree<Save>::iterator createRoot();
    void createChildren(const tree<Save>::iterator &node);
    void __getMaxValue(const tree<Save>::iterator &node, double& maxval, const std::vector< std::vector<double> >& zS, const std::vector< std::vector<double> >& zD);
    void __safePatternPruning(double lambda_prev, double lambda, const std::vector< std::vector<double> >& zS, const std::vector< std::vector<double> >& zD, double z2, double epsilon, const tree<Save>::iterator &node, vecVec& Xt, PatternList& patternList);
    double __workingSetPruning(double lambda, const std::vector< std::vector<double> >& zS, const std::vector< std::vector<double> >& zD, const tree<Save>::iterator &node, vecVec& Xt, PatternList& patternList);
    void __getTestX(const ProjectDB& pdb, vecVec& Xt, Pattern& pattern, const PatternList& patternList);
    void setKsample(const std::string &kernel);
    void nearest(const std::string &mat);
public:
    ItemsetMetric(const std::string& filename, int minsup, int maxpat, double train, int Ksample);
    double getMaxValue(const std::vector< std::vector<double> >& zS, const std::vector< std::vector<double> >& zD) ;
    void safePatternPruning(double lambda_prev, double lambda, const std::vector< std::vector<double> >& zS, const std::vector< std::vector<double> >& zD, double z2, double epsilon, vecVec& Xt, PatternList& patternList);
    double workingSetPruning(double lambda, const std::vector< std::vector<double> >& zS, const std::vector< std::vector<double> >& zD, vecVec& Xt, PatternList& patternList);
    void regularizationPath(double L=1, double U=0.95, double eta=1, int splitNum = 100, int maxloop = 100000, double eps = 1e-6, int freq = 5, double R=0.96, const std::string &kernel="");
    vecVec getTestX(const PatternList& patternList); 
    double error(const std::vector< std::vector<double> > &distance, const int K, int from, int to) const;
    double f_measure(const std::vector< std::vector<double> > &distance, const int K, int from, int to, bool macro=true) const;
    void run(const vecVec &active_Xt, const vecVec &test_Xt, const vector<double>& m, double lambda, double L, double U, double eta, double eps, unsigned loopMax);
    double kernel(const std::string &mat);
};

#endif /* ITEMSETMETRIC_HPP */


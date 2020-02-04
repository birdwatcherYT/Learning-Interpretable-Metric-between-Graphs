#include "ItemsetMetric.hpp"

using namespace std;

std::ostream& operator<<(std::ostream& os, const PatternSupportList& patternSupportList) {
    os << "[";
    for (const pair<Pattern, int>& p : patternSupportList) {
        os << "([ ";
        for (const int& element : p.first)
            os << element << " ";
        os << "], " << p.second << "), ";
    }
    os << "]" << endl;
    return os;
}

// this が 引数を含むかどうか
// 含む:1, 同じ:0, 含まない:-1
int contain(const Pattern &pattern1, const Pattern &pattern2) {
    size_t size = pattern2.size();
    if (pattern1.size() < size )
        return -1;
    for (size_t i=0;i<size;++i)
        if (pattern2[i] != pattern1[i])
            return -1;
    if (pattern1.size() == size )
        return 0;
    return 1;
}


ItemsetMetric::ItemsetMetric(const std::string& filename, int minsup, int maxpat, double train, int Ksample) {
    this->minsup = minsup;
    this->maxpat = maxpat;
    this->Ksample=Ksample;

    ifstream ifs(filename);
    if (ifs.fail()) {
        cerr << "ファイルオープンに失敗" << endl;
        exit(1);
    }
    //行ごと
    string line;
    double value;
    Transact temp;
    std::vector<int> y_temp;
    while (getline(ifs, line)) {
        stringstream ss(line);
        vector<int> items;
        ss >> value;
        y_temp.push_back(value);
        int item;
        while (ss >> item)
            items.push_back(item);
        temp.push_back(items);
    }

    this->trainSize=y_temp.size()* train;
    this->testSize=y_temp.size() - this->trainSize ;
    vector<size_t> random = randomIndex(y_temp.size());

    TRAIN.resize(trainSize); y_train.resize(trainSize); trainIndex.resize(trainSize);
    for(int i=0;i<trainSize;++i){
        TRAIN[i]=temp[random[i]];
        y_train[i]=y_temp[random[i]];
        trainIndex[i]=random[i];
        classLabel.insert(y_train[i]);
    }
    TEST.resize(testSize); y_test.resize(testSize); testIndex.resize(testSize);
    for(int i=0;i<testSize;++i){
        TEST[i]=temp[random[trainSize+i]];
        y_test[i]=y_temp[random[trainSize+i]];
        testIndex[i]=random[trainSize+i];
        classLabel.insert(y_test[i]);
    }
}

void ItemsetMetric::setKsample(const string &kernel){
    map<int, vector<int> > y2same;
    map<int, vector<int> > y2diff;
    for (int y : classLabel){
        for (int i=0;i<trainSize;++i){
            if (y_train[i]==y)
                y2same[y].push_back(i);
            else
                y2diff[y].push_back(i);
        }
        cout << "#y2same = " << y2same[y].size() << ", #y2diff = " << y2diff[y].size() <<endl;
    }
    sameClass.resize(trainSize);
    diffClass.resize(trainSize);

    if (kernel != ""){
        nearest(kernel);
    }else{
        for (int i=0;i<trainSize;++i){
            vector<size_t> randomSame = randomIndex(y2same[y_train[i]].size());
            vector<size_t> randomDiff = randomIndex(y2diff[y_train[i]].size());
            bool collision = false;
            sameClass[i].reserve(Ksample); diffClass[i].reserve(Ksample);
            for (int j=0;j<Ksample;++j){
                diffClass[i].push_back(y2diff[y_train[i]][randomDiff[j]]);
                if (i==y2same[y_train[i]][randomSame[j]]){
                    collision=true;
                    continue;
                }
                sameClass[i].push_back(y2same[y_train[i]][randomSame[j]]);
            }
            if (collision)
                sameClass[i].push_back(y2same[y_train[i]][randomSame[Ksample]]);
        }
    }
}

void ItemsetMetric::nearest(const std::string &mat){
    ifstream ifs(mat);
    if (ifs.fail()){
        cerr<<"file open error"<<endl;
        exit(1);
    }
    string str;
    vector< vector<double> > GraphKernel;
    while (getline(ifs, str)){
        stringstream ss(str);
        vector<double> k;
        while(!ss.eof()){
            double value;
            ss >> value;
            k.push_back(value);
        }
        GraphKernel.push_back(k);
    }

    vector< vector<double> > distance(trainSize);
    #pragma omp parallel for
    for (int I=0; I < trainSize; ++I){
        int i=trainIndex[I];
        distance[I].resize(trainSize);
        for (int J=0; J < trainSize; ++J){
            int j=trainIndex[J];
            distance[I][J] = GraphKernel[i][i]-2*GraphKernel[i][j]+GraphKernel[j][j];
        }
    }

    for (size_t i=0;i<trainSize;++i){
        vector<size_t> sortIndex = argsort(distance[i]);
        sameClass[i].reserve(Ksample); diffClass[i].reserve(Ksample);
        for(size_t j:sortIndex){
            if (y_train[i]!=y_train[j]){//diff class
                if (diffClass[i].size() < Ksample )
                    diffClass[i].push_back(j);
            }else{//same class
                if (i!=j && sameClass[i].size() < Ksample )
                    sameClass[i].push_back(j);
            }
            if (sameClass[i].size() >= Ksample && diffClass[i].size() >= Ksample )
                break;
        }
    }
}


tree<ItemsetMetric::Save>::iterator ItemsetMetric::createRoot(){
    if (Tree.empty()){
        ProjectDB pdb;
        for (int i = 0; i < trainSize; ++i)
            pdb.push_back(pair<int, int>(i, -1));
        tree<Save>::iterator root=Tree.insert(Tree.begin(), Save(vector<bool>(), Pattern(), pdb));
        createChildren(root);
    }
    return Tree.begin();
}

void ItemsetMetric::createChildren(const tree<Save>::iterator &node){
    if (!node->nextCheck){
        node->nextCheck=true;
        if (node->pattern.size() >= maxpat)
            return;
        map<int, ProjectDB> counter;
        for (const pair<int, int>& p : node->pdb) {
            int id = p.first;
            for (int j = p.second + 1; j < TRAIN[id].size(); ++j)
                counter[TRAIN[id][j]].push_back(pair<int, int>(id, j));
        }
        Pattern pattern(node->pattern);
        for (const auto& it : counter) {
            vector<bool> x(trainSize, 0);
            int support = 0, oldID = -1;
            for (const pair<int, int>& p : it.second) {
                int id = p.first;
                if (oldID != id) {
                    ++support;
                    x[id] = 1;
                    oldID = id;
                }
            }
            if (support < minsup)
                continue;
            pattern.push_back(it.first);
            Tree.append_child(node, Save(x, pattern, it.second));
            pattern.pop_back();
        }
    }
}

double ItemsetMetric::getMaxValue(const std::vector< std::vector<double> >& zS, const std::vector< std::vector<double> >& zD)  {
    visit=0;
    double maxval = 0;
    tree<Save>::iterator root=createRoot();
    for (tree<Save>::sibling_iterator it=Tree.begin(root);it!=Tree.end(root);++it)
        __getMaxValue(it, maxval, zS, zD);
    cout<<"getMaxValue:visit = "<<visit<<endl;
    return maxval;
}

void ItemsetMetric::__getMaxValue(const tree<Save>::iterator &node, double& maxval, const std::vector< std::vector<double> >& zS, const std::vector< std::vector<double> >& zD){
    ++visit;

    double dot = 0;
    double dotSPPC = 0;
    for (int i = 0; i < TRAIN.size(); ++i) {
        double tmp0=0, tmp1=0;
        for (int a=0,j=sameClass[i][a];a<Ksample;++a,j=sameClass[i][a]){
            bool cij = (node->x[i] != node->x[j]);
            if (cij)
                dot -= zS[i][a];
            if (node->x[i] && !node->x[j])
                tmp1 -= zS[i][a];
        }
        for (int a=0,l=diffClass[i][a];a<Ksample;++a,l=diffClass[i][a]){
            bool cil = (node->x[i] != node->x[l]);
            if (cil)
                dot += zD[i][a];
            if (node->x[l])
                tmp0 += zD[i][a];
            if (node->x[i])
                tmp1 += zD[i][a];
        }
        dotSPPC += fmax(tmp0, tmp1);
    }

    if (dotSPPC <= maxval)
        return;
    maxval = fmax(maxval, dot);

    createChildren(node);

    for (tree<Save>::sibling_iterator it=Tree.begin(node);it!=Tree.end(node);++it)
        __getMaxValue(it, maxval, zS, zD);
    return;
}

void ItemsetMetric::safePatternPruning(double lambda_prev, double lambda, const std::vector< std::vector<double> >& zS, const std::vector< std::vector<double> >& zD, double z2, double epsilon_d, vecVec& Xt, PatternList& patternList) {
    clock_t start=clock();
    Xt.clear();
    patternList.clear();
    visit=0;

    tree<Save>::iterator root=createRoot();

    for (tree<Save>::sibling_iterator it=Tree.begin(root);it!=Tree.end(root);++it)
        __safePatternPruning(lambda_prev, lambda, zS, zD, z2, epsilon_d, it, Xt, patternList);
    clock_t end=clock();
    cout<<"safePatternPruning:visit = "<<visit<<endl;
    cout<<"safePatternPruning:time = "<<(double)(end-start)/CLOCKS_PER_SEC<<endl;
}

void ItemsetMetric::__safePatternPruning(double lambda_prev, double lambda, const std::vector< std::vector<double> >& zS, const std::vector< std::vector<double> >& zD, double z2, double epsilon_d, const tree<Save>::iterator &node, vecVec& Xt, PatternList& patternList) {

    if (node->pruningScore<=lambda){
        return;
    }else if (node->screeningScore>lambda){
        ++visit;
        double dot = 0, h_norm = 0;
        double dotSPPC = 0, h_normSPPC = 0;
        for (int i = 0; i < trainSize; ++i) {
            double tmp0=0, tmp1=0;
            for (int a=0,j=sameClass[i][a];a<Ksample;++a,j=sameClass[i][a]){
                bool cij = (node->x[i] != node->x[j]);
                if (cij)
                    dot -= zS[i][a], ++h_norm;
                if (node->x[i] && !node->x[j])
                    tmp1 -= zS[i][a];
                if (node->x[i] || node->x[j])
                    ++h_normSPPC;
            }
            for (int a=0,l=diffClass[i][a];a<Ksample;++a,l=diffClass[i][a]){
                bool cil = (node->x[i] != node->x[l]);
                if (cil)
                    dot += zD[i][a], ++h_norm;
                if (node->x[l])
                    tmp0 += zD[i][a];
                if (node->x[i])
                    tmp1 += zD[i][a];
                if (node->x[i] || node->x[l])
                    ++h_normSPPC;
            }
            dotSPPC += fmax(tmp0, tmp1);
        }
        node->pruningScore = lambda_prev*(2*epsilon_d*sqrt(h_normSPPC)+sqrt(h_normSPPC*z2)+dotSPPC)/(2*lambda_prev+sqrt(h_normSPPC*z2)-dotSPPC);
        if (node->pruningScore <= lambda)
            return;
        if (dotSPPC <= lambda)// WS
            return;
        node->screeningScore = lambda_prev*(2*epsilon_d*sqrt(h_norm)+sqrt(h_norm*z2)+dot)/(2*lambda_prev+sqrt(h_norm*z2)-dot);
        if (dot > lambda && node->screeningScore > lambda) {
            Xt.push_back(node->x);
            patternList.push_back(node->pattern);
        }
    }
    
    createChildren(node);
    
    for (tree<Save>::sibling_iterator it=Tree.begin(node);it!=Tree.end(node);++it)
        __safePatternPruning(lambda_prev, lambda, zS, zD, z2, epsilon_d, it, Xt, patternList);
    return;
}

double ItemsetMetric::workingSetPruning(double lambda, const std::vector< std::vector<double> >& zS, const std::vector< std::vector<double> >& zD, vecVec& Xt, PatternList& patternList) {
    clock_t start=clock();
    Xt.clear();
    patternList.clear();
    visit=0;

    double m_alpha_lambda_eta2=0;
    tree<Save>::iterator root=createRoot();

    for (tree<Save>::sibling_iterator it=Tree.begin(root);it!=Tree.end(root);++it)
        m_alpha_lambda_eta2+=__workingSetPruning(lambda, zS, zD, it, Xt, patternList);
    clock_t end=clock();
    cout<<"workingSetPruning:visit = "<<visit<<endl;
    cout<<"workingSetPruning:time = "<<(double)(end-start)/CLOCKS_PER_SEC<<endl;
    return m_alpha_lambda_eta2;
}

double ItemsetMetric::__workingSetPruning(double lambda, const std::vector< std::vector<double> >& zS, const std::vector< std::vector<double> >& zD, const tree<Save>::iterator &node, vecVec& Xt, PatternList& patternList) {
    
    double m_alpha_lambda_eta2=0;
    if (node->pruningScore<=lambda){
        return 0;
    }else if (node->screeningScore>lambda){
        ++visit;
        double dot = 0;
        double dotSPPC = 0;
        for (int i = 0; i < trainSize; ++i) {
            double tmp0=0, tmp1=0;
            for (int a=0,j=sameClass[i][a];a<Ksample;++a,j=sameClass[i][a]){
                bool cij = (node->x[i] != node->x[j]);
                if (cij)
                    dot -= zS[i][a];
                if (node->x[i] && !node->x[j])
                    tmp1 -= zS[i][a];
            }
            for (int a=0,l=diffClass[i][a];a<Ksample;++a,l=diffClass[i][a]){
                bool cil = (node->x[i] != node->x[l]);
                if (cil)
                    dot += zD[i][a];
                if (node->x[l])
                    tmp0 += zD[i][a];
                if (node->x[i])
                    tmp1 += zD[i][a];
            }
            dotSPPC += fmax(tmp0, tmp1);
        }
        if (dotSPPC <= lambda)
            return 0;
        if (dot > lambda) {
            m_alpha_lambda_eta2+=(dot-lambda)*(dot-lambda);
            Xt.push_back(node->x);
            patternList.push_back(node->pattern);
        }
    }
    
    createChildren(node);
    
    for (tree<Save>::sibling_iterator it=Tree.begin(node);it!=Tree.end(node);++it)
        m_alpha_lambda_eta2+=__workingSetPruning(lambda, zS, zD, it, Xt, patternList);
    return m_alpha_lambda_eta2;
}

void ItemsetMetric::regularizationPath(double L, double U, double eta, int splitNum, int maxloop, double eps, int freq, double R, const std::string &kernel) {
    setKsample(kernel);
    
    cout<<"L = "<<L<<", U = "<<U<<", eta = "<<eta<<endl;
    cout<<"trainSize = "<<trainSize<<", testSize = "<<testSize<<", Ksample = "<<Ksample<<endl;
    cout<<"maxpat = "<<maxpat<<", minsup = "<<minsup<<", R = "<<R<<endl;

    vector<double> m;
    vector<double> m_prev(m);
    vector<double> grad;
    vector<double> grad_prev(grad);


    double zD1=0, zS1=0, z2=0;
    double loss = 0;
    vector< vector<double> > zD(trainSize);
    vector< vector<double> > zS(trainSize);
    for (int i = 0; i < trainSize; ++i){
        zS[i].resize(Ksample,0);
        zD[i].resize(Ksample,2*L);
    }
    zD1=2*L*Ksample*trainSize;
    z2=2*L*2*L*Ksample*trainSize;
    loss=L*L*Ksample*trainSize;
    
    double lambda_prev = getMaxValue(zS, zD); //LAMBDA_MAX
    cout << "lambda_max = " << lambda_prev << endl;
    double m1 = 0, m2 = 0;
    map<Pattern, double> mSaveDict;
    vecVec Xt;
    PatternList pat;
    size_t dim = 0;
    double epsilon_d=0;
    
    double validMaxF=0;
    double correspondTestF;

    vecVec bestActive_Xt, bestTestXt;
    vector<double> bestActive_m;
    double bestLambda(0);

    vector<double> hMaxWS;
    vector< vector<double> > zSMaxWS=zS, zDMaxWS=zD;
    double z2MaxWS=z2, zS1MaxWS=zS1, zD1MaxWS=zD1;
    vector< vector<double> > zSMax=zS, zDMax=zD;
    double z2Max=z2;

    for (int idx = 1; idx < splitNum; ++idx) {
        double lambda = R * lambda_prev;

        double alpha = 1e-6;
        double primal_prev = DBL_MAX;
        vector<size_t> index;
        double dualMaxWS = -DBL_MAX;
        double dualMax = -DBL_MAX;

        cout << "********************************" << endl;
        cout << "[" << idx << "] lambda = " << lambda << endl;
        clock_t start = clock();

        cout<<"rrpb_d = "<<sqrt(z2Max)*((lambda_prev-lambda)/(2*lambda_prev))+epsilon_d<<endl;
        safePatternPruning(lambda_prev, lambda, zSMax, zDMax, z2Max, epsilon_d, Xt, pat);
        m.clear();
        for (const Pattern &it : pat)
            m.push_back(mSaveDict[it]);
        dim = m.size();
        index.resize(dim);
        for (size_t k = 0; k < dim; ++k)
            index[k] = k;

        cout<<"dim = "<<dim<<endl;
        for (int loop = 0; loop < maxloop; ++loop) {
            cout << "---------- loop  = " << loop << "----------" << endl;
            vector<double> h(dim, 0);
            vector<double> h_norm(dim, 0);
            //-----------------------------
            zD1 = zS1 = z2 = m1 = m2 = loss = 0;
            for (int i = 0; i < trainSize; ++i) {
                for (int a=0,j=sameClass[i][a];a<Ksample;++a,j=sameClass[i][a]){
                    vector<bool> cij(dim);
                    double value = 0;
                    for (size_t k : index) {
                        cij[k] = (Xt[k][i] != Xt[k][j]);
                        if (cij[k])
                            value += m[k], ++h_norm[k];
                    }
                    double hinge = fmax(0,-U+value);
                    for (size_t k : index)
                        h[k] -= cij[k]*2*hinge;
                    loss += hinge*hinge, zS1 += 2*hinge, z2 += (2*hinge)*(2*hinge), zS[i][a] = 2*hinge;
                }
                for (int a=0, l=diffClass[i][a];a<Ksample;++a, l=diffClass[i][a]){
                    vector<bool> cil(dim);
                    double value = 0;
                    for (size_t k : index) {
                        cil[k] = (Xt[k][i] != Xt[k][l]);
                        if (cil[k])
                            value += m[k], ++h_norm[k];
                    }
                    double hinge = fmax(0,L-value);
                    for (size_t k : index)
                        h[k] += cil[k]*2*hinge;
                    loss += hinge*hinge, zD1 += 2*hinge, z2 += (2*hinge)*(2*hinge), zD[i][a] = 2*hinge;
                }
            }
            for (size_t k : index)
                m1 += m[k], m2 += m[k]*m[k];
            //-----------------------------
            
            double primal = loss + lambda * (m1+0.5*eta*m2);
            if (primal > primal_prev) {
                alpha *= 0.1, --loop;
                cout << "alpha = " << alpha << endl;
                for (size_t k : index) {
                    m[k] = m_prev[k] - alpha * grad[k];
                    if (m[k] < 0) m[k] = 0;
                }
                continue;
            }
            double m_alpha2 = 0;
            for (size_t k : index)
                if (h[k]>lambda)
                    m_alpha2 += (h[k]-lambda)*(h[k]-lambda);
            double dual = -0.25 * z2 + L*zD1-U*zS1-0.5*m_alpha2/(lambda*eta);

            if (dualMaxWS<dual)
                dualMaxWS = dual, hMaxWS=h, zSMaxWS=zS, zDMaxWS=zD, z2MaxWS=z2, zS1MaxWS=zS1, zD1MaxWS=zD1;
            
            double gap = primal - dualMaxWS;
            cout << "primal = " << primal << ", dual = " << dualMaxWS << ", gap = " << gap << endl;
            if (gap <= eps * primal) {// WS内での収束
                cout<<"***converge***"<<endl;

                vecVec Xt_temp;
                PatternList pat_temp;
                
                gap = primal - dualMax;
                if (!(gap <= eps * primal)){// これまでのDualの値で収束しなかったとき
                    m_alpha2=workingSetPruning(lambda, zSMaxWS, zDMaxWS, Xt_temp, pat_temp);
                    dual = -0.25 * z2MaxWS + L*zD1MaxWS-U*zS1MaxWS-0.5*m_alpha2/(lambda*eta);
                    if (dualMax<dual)
                        dualMax = dual, zSMax=zSMaxWS, zDMax=zDMaxWS;
                    gap = primal - dualMax;
                }
                cout << "primal = " << primal << ", dual = " << dualMax << ", gap = " << gap << endl;
                if (gap <= eps * primal){
                    // convergence
                } else {
                    mSaveDict.clear();
                    for (int k=0;k<m.size();++k)
                        mSaveDict[pat[k]] = m[k];
                    Xt=Xt_temp, pat=pat_temp;
                    m.clear();
                    for (const Pattern &it : pat)
                        m.push_back(mSaveDict[it]);
                    dim = m.size();
                    index.resize(dim);
                    for (size_t k = 0; k < dim; ++k)
                        index[k] = k;
                    cout << "dim = " << index.size() << endl;
                    
                    primal_prev = DBL_MAX;
                    dualMaxWS = -DBL_MAX;

                    continue;
                }
                if (gap<0)
                    gap=0;
                epsilon_d = 2*sqrt(gap);
                vector<double> active_m;
                PatternList active_pat;
                vecVec active_Xt;
                for (size_t k : index){
                    if (m[k]!=0){
                        active_m.push_back(m[k]);
                        active_pat.push_back(pat[k]);
                        active_Xt.push_back(Xt[k]);
                    }
                }
                int activeSize  = active_m.size();
                cout<< "active = " << activeSize <<endl;

                mSaveDict.clear();
                for (int k=0;k<active_m.size();++k )
                    mSaveDict[active_pat[k]] = active_m[k];
                cout << "m = " << active_m <<endl;
                cout << "pat = "<< active_pat<<endl;
                cout << "time = " << (double) (clock() - start) / CLOCKS_PER_SEC << endl;
                
                clock_t s=clock();
                vecVec testXt = getTestX(active_pat);
                // cout << testXt <<endl;
                cout << "testTime = " << (double) (clock() - s) / CLOCKS_PER_SEC << endl;

                vector< vector<double> >distance(testSize);
                #pragma omp parallel for
                for (int i=0;i<testSize;++i){
                    distance[i].resize(trainSize,0);
                    for (int j=0;j<trainSize;++j)
                        for (int k=0;k<activeSize;++k)
                            if (testXt[k][i] != active_Xt[k][j])
                                distance[i][j] += active_m[k];
                }

                //////---------------------
                for (int k=1;k<50;k+=2){
                    double eV=error(distance, k, 0, testSize/2);
                    double eT=error(distance, k, testSize/2, testSize);
                    double fmacroV=f_measure(distance, k, 0, testSize/2);
                    double fmacroT=f_measure(distance, k, testSize/2, testSize);
                    double fmicroV=f_measure(distance, k, 0, testSize/2, false);
                    double fmicroT=f_measure(distance, k, testSize/2, testSize, false);
                    cout<<k<<"-nn: ValidError = "<< eV <<", TestError = "<< eT<< endl;
                    cout<<k<<"-nn: Valid-Fmacro = "<< fmacroV <<", Test-Fmacro = "<< fmacroT<< endl;
                    cout<<k<<"-nn: Valid-Fmicro = "<< fmicroV <<", Test-Fmicro = "<< fmicroT << endl;
                    if (validMaxF<fmicroV)
                        validMaxF=fmicroV, correspondTestF=fmicroT, bestActive_Xt=active_Xt, bestTestXt=testXt, bestActive_m=active_m, bestLambda=lambda;
                }
                // run(active_Xt, testXt, active_m, lambda, L, U, eta, eps, maxloop);
                break;
            }
            if (loop % freq == 0) {
                double r = 2*sqrt(gap);
                cout<<"r_dgb = "<<r<<endl;
                if (loop != 0) {
                    vector<size_t> newIndex;
                    newIndex.reserve(index.size());
                    for (size_t k : index) {
                        double ub = hMaxWS[k] + r * sqrt(h_norm[k]);
                        if (ub > lambda)
                            newIndex.push_back(k);
                    }
                    index = newIndex;
                }
                cout << "dim = " << index.size() << endl;
            }

            grad.resize(h.size());
            for (size_t k : index)
                grad[k] = -h[k] + lambda*(1+eta*m[k]);
            if (primal_prev != DBL_MAX) {
                double dot = 0, g_sqr = 0, m_sqr = 0;
                for (size_t k : index) {
                    double m_diff = m[k] - m_prev[k], g_diff = grad[k] - grad_prev[k];
                    dot += m_diff*g_diff, g_sqr += g_diff*g_diff, m_sqr += m_diff*m_diff;
                }
                alpha = 0.5 * fabs(dot / g_sqr + m_sqr / dot);
                if (std::isnan(alpha)) {
                    cout << "alpha is nan." << endl;
                    // exit(1);
                    alpha=1e-6;
                }
            }
            primal_prev = primal, m_prev = m, grad_prev = grad;
            //update m
            for (size_t k : index) {
                m[k] -= alpha * grad[k];
                if (m[k] < 0) m[k] = 0;
            }
        }
        lambda_prev = lambda;
    }
    cout<<"validMaxF = "<<validMaxF<<", correspondTestF = "<<correspondTestF<<endl;
    run(bestActive_Xt, bestTestXt, bestActive_m, bestLambda, L, U, eta, eps, maxloop);
}

vecVec ItemsetMetric::getTestX(const PatternList& patternList) {
    ProjectDB pdb;
    for (int i = 0; i < testSize; ++i)
        pdb.push_back(pair<int, int>(i, -1));

    vecVec Xt;

    Pattern pattern;
    __getTestX(pdb, Xt, pattern, patternList);
    for (int j=0, count=patternList.size() - Xt.size(); j < count; ++j )
        Xt.push_back(vector<bool>(TEST.size(), 0));
    return Xt;
}

void ItemsetMetric::__getTestX(const ProjectDB& pdb, vecVec& Xt, Pattern& pattern, const PatternList& patternList) {
    if (pattern.size() > maxpat)
        return;
    bool flag=true;
    for (int i=Xt.size();i<patternList.size();++i){
        if (contain(patternList[i], pattern) >= 0){
            flag = false;
            break;
        }
    }
    if (flag)
        return;
    map<int, ProjectDB> counter;
    for (const pair<int, int>& p : pdb) {
        int id = p.first;
        for (int j = p.second + 1; j < TEST[id].size(); ++j)
            counter[TEST[id][j]].push_back(pair<int, int>(id, j));
    }
    for (const auto& it : counter) {
        pattern.push_back(it.first);
        for (int i=Xt.size();i<patternList.size();++i ){
            if (contain(patternList[i], pattern) == 0){
                for (int j=0, count=i - Xt.size(); j < count;++j )
                    Xt.push_back(vector<bool>(TEST.size(), 0));
                vector<bool> x(testSize, 0);
                for (const pair<int, int>& p : it.second)
                    x[p.first] = 1;
                Xt.push_back(x);
                break;
            }
        }
        __getTestX(it.second, Xt, pattern, patternList);
        pattern.pop_back();
    }
}

// distance : test x train
double ItemsetMetric::error(const std::vector< std::vector<double> > &distance, const int K, int from, int to) const{
    double err=0;
    #pragma omp parallel 
    {
        #pragma omp for reduction(+:err)
        for (size_t i=from;i<to;++i){
            vector<size_t> sortIndex = argsort(distance[i]);
            map<int,int> classMap;
            for (size_t k = 0; k < K; ++k)
                classMap[y_train[sortIndex[k]]]++;
            int max=-1, argmax=0;
            for (const pair<int,int> &p : classMap){
                if (p.second > max)
                    max = p.second, argmax = p.first;
            }
            err += (y_test[i] != argmax);
        }
    }
    return err / (to-from);
}

// distance : test x train
double ItemsetMetric::f_measure(const std::vector< std::vector<double> > &distance, const int K, int from, int to, bool macro) const{
    int classNum = classLabel.size();
    vector<int> TP(classNum, 0), TP_FP(classNum, 0), TP_FN(classNum, 0);
    
    #pragma omp declare reduction(+ : std::vector<int> : \
                          std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<int>())) \
                initializer(omp_priv = omp_orig)
    #pragma omp parallel for reduction(+:TP,TP_FP,TP_FN)
    for (size_t i=from;i<to;++i){
        vector<size_t> sortIndex = argsort(distance[i]);
        map<int,int> classMap;
        for (size_t k = 0; k < K; ++k)
            classMap[y_train[sortIndex[k]]]++;
        int max=-1, argmax=0;
        for (const pair<int,int> &p : classMap){
            if (p.second > max)
                max = p.second, argmax = p.first;
        }
        int j=0;
        for (int c : classLabel){
            // F-measure
            TP_FN[j] += (y_test[i]==c);
            TP_FP[j] += (argmax==c);
            TP[j] += (argmax==c && y_test[i]==c);
            ++j;
        }
    }
    
    if (macro){
        double f_macro=0;
        for (int j=0;j<classNum;++j){
            double precision = (double) TP[j]/TP_FP[j];
            double recall = (double) TP[j]/TP_FN[j];
            f_macro += 2*precision*recall/(recall+precision);
        }
        f_macro /= classNum;
        return f_macro;
    }
    int TP_sum=accumulate(TP.begin(), TP.end(), 0), TP_FP_sum=accumulate(TP_FP.begin(), TP_FP.end(), 0), TP_FN_sum=accumulate(TP_FN.begin(), TP_FN.end(), 0);
    double precision = (double) TP_sum/TP_FP_sum;
    double recall = (double) TP_sum/TP_FN_sum;
    double f_micro = 2*precision*recall/(recall+precision);
    return f_micro;
}

void ItemsetMetric::run(const vecVec &active_Xt, const vecVec &test_Xt, const vector<double>& m, double lambda, double L, double U, double eta, double eps, unsigned loopMax) {

    int dim = active_Xt.size(), trainSize=active_Xt[0].size(), testSize=test_Xt[0].size();
    cout<<"Full Matrix Dim = "<<dim<<endl;
    MatrixXd X_train = MatrixXd::Zero(trainSize, dim);
    for (int i=0;i<X_train.rows();++i)
        for (int k=0;k<dim;++k)
            X_train(i,k)=active_Xt[k][i];
    MatrixXd X_test = MatrixXd::Zero(testSize, dim);
    for (int i=0;i<X_test.rows();++i)
        for (int k=0;k<dim;++k)
            X_test(i,k)=test_Xt[k][i];

    double primal_prev = DBL_MAX;
    MatrixXd M = MatrixXd::Zero(dim, dim);
    for (int k=0;k<dim;++k)
        M(k,k)=m[k];
    MatrixXd M_prev = M;

    MatrixXd Grad = MatrixXd::Zero(dim, dim);
    MatrixXd Grad_prev = Grad;

    double zD1, zS1, z2, m1, m2, loss;
    vector< vector<double> > zS(trainSize);
    vector< vector<double> > zD(trainSize);
    for (int i = 0; i < trainSize; ++i){
        zS[i].resize(Ksample,0);
        zD[i].resize(Ksample,0);
    }
    double alpha = 1e-6;
    double dualMax=-DBL_MAX;
    
    clock_t start = clock();
    for (unsigned loop = 0; loop < loopMax; ++loop) {
        cout << "-----" << "loop = " << loop << "-----" << endl;

        MatrixXd H = MatrixXd::Zero(dim,dim);
        zD1 = zS1 = z2 = m1 = m2 = loss = 0;
        
        #pragma omp declare reduction(+ : Eigen::MatrixXd : omp_out=omp_out+omp_in) initializer(omp_priv = omp_orig)
        #pragma omp parallel for reduction(+:H, loss, zS1, zD1, z2)
        for (int i = 0; i < X_train.rows(); ++i) {
            for (int a=0,j=sameClass[i][a];a<Ksample;++a,j=sameClass[i][a]){
                VectorXd cij = X_train.row(i)-X_train.row(j);
                double value = cij.dot(M*cij);
                double hinge = fmax(0,-U+value);
                H -= 2*hinge*(cij*cij.transpose());
                loss += hinge*hinge, zS1 += 2*hinge, z2 += (2*hinge)*(2*hinge), zS[i][a] = 2*hinge;
            }
            for (int a=0,l=diffClass[i][a];a<Ksample;++a,l=diffClass[i][a]){
                VectorXd cil = X_train.row(i)-X_train.row(l);
                double value = cil.dot(M*cil);
                double hinge = fmax(0,L-value);
                H += 2*hinge*(cil*cil.transpose());
                loss += hinge*hinge, zD1 += 2*hinge, z2 += (2*hinge)*(2*hinge), zD[i][a] = 2*hinge;
            }
        }
        
        m2 = M.squaredNorm();
        m1 = M.trace();

        double primal = loss + lambda * (m1+0.5*eta*m2);
        if (primal > primal_prev) {
            cout << "primal = " << primal  << endl;
            alpha *= 0.1;
            cout << "alpha = " << alpha << endl;
            --loop;
            M = projectSemidefinite(M_prev - alpha * Grad);
            continue;
        }
        double m_alpha2 = projectSemidefinite(H-lambda*MatrixXd::Identity(dim,dim)).squaredNorm();
        double dual = -0.25 * z2 + L*zD1-U*zS1-0.5*m_alpha2/(lambda*eta);
        if (dualMax<dual)
            dualMax=dual;
        double gap = primal - dualMax;
        cout << "primal = " << primal  << ", dual = " << dualMax <<  ", gap = " << gap << endl;
        if (gap <= eps * fabs(primal))
            break;

        Grad = -H + lambda * (MatrixXd::Identity(dim,dim)+eta*M);

        if (loop != 0){
            const MatrixXd &M_diff = M-M_prev;
            const MatrixXd &G_diff = Grad-Grad_prev;
            const double dot = M_diff.cwiseProduct(G_diff).sum();
            const double G_sqr = G_diff.cwiseProduct(G_diff).sum();
            const double M_sqr = M_diff.cwiseProduct(M_diff).sum();
            alpha = 0.5 * fabs(dot/G_sqr + M_sqr/dot);
            if (std::isnan(alpha)){
                cout<<"alpha is nan."<<endl;
                exit(1);
            }
        }
        primal_prev = primal, M_prev = M, Grad_prev = Grad;
        M = projectSemidefinite(M - alpha * Grad);
    }
    clock_t end = clock();
    cout << "Fulltime = " << (double) (end - start) / CLOCKS_PER_SEC << endl;

    vector< vector<double> >distance(testSize);
    #pragma omp parallel for
    for (int i=0;i<testSize;++i){
        distance[i].resize(trainSize,0);
        for (int j=0;j<trainSize;++j){
            VectorXd cij = X_test.row(i)-X_train.row(j);
            distance[i][j] = cij.dot(M*cij);
        }
    }
    for (int k=1;k<50;k+=2){
        double eV=error(distance, k, 0, testSize/2);
        double eT=error(distance, k, testSize/2, testSize);
        double fmacroV=f_measure(distance, k, 0, testSize/2);
        double fmacroT=f_measure(distance, k, testSize/2, testSize);
        double fmicroV=f_measure(distance, k, 0, testSize/2, false);
        double fmicroT=f_measure(distance, k, testSize/2, testSize, false);
        cout<<k<<"-nn: FullValidError = "<< eV <<", FullTestError = "<< eT<< endl;
        cout<<k<<"-nn: FullValid-Fmacro = "<< fmacroV <<", FullTest-Fmacro = "<< fmacroT<< endl;
        cout<<k<<"-nn: FullValid-Fmicro = "<< fmicroV <<", FullTest-Fmicro = "<< fmicroT << endl;
    }
    cout<<"M2 = "<<m2<<endl;
}

double ItemsetMetric::kernel(const std::string &mat){
    ifstream ifs(mat);
    if (ifs.fail()){
        cerr<<"file open error"<<endl;
        exit(1);
    }
    string str;
    vector< vector<double> > GraphKernel;
    while (getline(ifs, str)){
        stringstream ss(str);
        vector<double> k;
        while(!ss.eof()){
            double value;
            ss >> value;
            k.push_back(value);
        }
        GraphKernel.push_back(k);
    }

    vector< vector<double> > distance(testSize);
    #pragma omp parallel for
    for (int I=0; I < testSize; ++I){
        int i=testIndex[I];
        distance[I].resize(trainSize);
        for (int J=0; J < trainSize; ++J){
            int j=trainIndex[J];
            distance[I][J] = GraphKernel[i][i]-2*GraphKernel[i][j]+GraphKernel[j][j];
        }
    }

    //////---------------------
    double validMaxF=0, correspondTestF;
    for (int k=1;k<50;k+=2){
        double eV=error(distance, k, 0, testSize/2);
        double eT=error(distance, k, testSize/2, testSize);
        double fmacroV=f_measure(distance, k, 0, testSize/2);
        double fmacroT=f_measure(distance, k, testSize/2, testSize);
        double fmicroV=f_measure(distance, k, 0, testSize/2, false);
        double fmicroT=f_measure(distance, k, testSize/2, testSize, false);
        cout<<k<<"-nn: KernelValidError = "<< eV <<", KernelTestError = "<< eT<< endl;
        cout<<k<<"-nn: KernelValid-Fmacro = "<< fmacroV <<", KernelTest-Fmacro = "<< fmacroT<< endl;
        cout<<k<<"-nn: KernelValid-Fmicro = "<< fmicroV <<", KernelTest-Fmicro = "<< fmicroT << endl;
        if (validMaxF<fmicroV)
            validMaxF=fmicroV, correspondTestF=fmicroT;
    }
    cout<<"KernelValidMaxF = "<<validMaxF<<", KernelCorrespondTestF = "<<correspondTestF<<endl;
    return validMaxF;
}


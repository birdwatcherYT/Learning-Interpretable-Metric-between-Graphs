#include "gspan.h"

using namespace std;

namespace GSPAN {

    gSpan::gSpan(const string& filename, unsigned int minsup, unsigned int maxpat, double train, int Ksample, bool ignoreEdgeLabel, int VertexLabelDegreeQuantization, bool directed) {
        this->minsup = minsup;
        this->maxpat = maxpat;
        this->directed = directed;
        this->Ksample=Ksample;
        ifstream fp(filename);
        read(fp, train, ignoreEdgeLabel, VertexLabelDegreeQuantization);
        /* Do single node handling, as the normal gspan DFS code based processing
         * cannot find subgraphs of size |subg|==1.  Hence, we find frequent node
         * labels explicitly.
         */
        for (unsigned int id = 0; id < TRAIN.size(); ++id) {
            for (unsigned int nid = 0 ; nid < TRAIN[id].size() ; ++nid) {
                if (singleVertexTRAIN[id][TRAIN[id][nid].label] == 0) {
                    // number of graphs it appears in
                    singleVertexLabelTRAIN[TRAIN[id][nid].label] += 1;// support
                }
                singleVertexTRAIN[id][TRAIN[id][nid].label] += 1;//count
            }
        }
        for (unsigned int id = 0; id < TEST.size(); ++id) {
            for (unsigned int nid = 0 ; nid < TEST[id].size() ; ++nid) {
                if (singleVertexTEST[id][TEST[id][nid].label] == 0) {
                    // number of graphs it appears in
                    singleVertexLabelTEST[TEST[id][nid].label] += 1;// support
                }
                singleVertexTEST[id][TEST[id][nid].label] += 1;//count
            }
        }
        for (map<unsigned int, unsigned int>::iterator it = singleVertexLabelTRAIN.begin () ; it != singleVertexLabelTRAIN.end () ; ++it) {
            if (it->second < minsup)
                continue;

            vector<bool> x(TRAIN.size(),0);
            for (map<unsigned int, map<unsigned int, unsigned int> >::iterator cur = singleVertexTRAIN.begin(); cur != singleVertexTRAIN.end(); ++cur)
                x[cur->first] = (cur->second[it->first]>0);

            DFS_CODE.createSingleVertex(it->first);
            singleTree.push_back(Save(x, DFS_CODE, NULL));
            DFS_CODE.clear();
        }
        cout<<"filename = "<<filename<<endl;
    }

    std::istream &gSpan::read(std::istream &is, double train, bool ignoreEdgeLabel, int VertexLabelDegreeQuantization) {
        Graph g(directed);
        std::vector < Graph > temp;
        while (true) {
            g.read(is, ignoreEdgeLabel);
            if (g.empty()) break;
            temp.push_back(g);
        }

        if (VertexLabelDegreeQuantization>0){
            int maxDegree=0;
            for (Graph &g : temp)
                for (Vertex &v : g)
                    if (maxDegree<v.edge.size())
                        maxDegree=v.edge.size();
            cout<<"maxDegree = "<<maxDegree<<endl;
            for (Graph &g : temp)
                for (Vertex &v : g)
                    v.label=v.edge.size()/VertexLabelDegreeQuantization*VertexLabelDegreeQuantization;
        }

        this->trainSize=temp.size()* train;
        this->testSize=temp.size() - this->trainSize ;
        vector<size_t> random = randomIndex(temp.size());

        TRAIN.resize(trainSize); trainIndex.resize(trainSize);
        for(int i=0;i<trainSize;++i){
            TRAIN[i]=temp[random[i]];
            trainIndex[i]=random[i];
            classLabel.insert(TRAIN[i].y);
        }
        TEST.resize(testSize); testIndex.resize(testSize);
        for(int i=0;i<testSize;++i){
            TEST[i]=temp[random[trainSize+i]];
            testIndex[i]=random[trainSize+i];
            classLabel.insert(TEST[i].y);
        }
        return is;
    }
    
    void gSpan::setKsample(const string &kernel){
        map<int, vector<int> > y2same;
        map<int, vector<int> > y2diff;
        for (int y : classLabel){
            for (int i=0;i<trainSize;++i){
                if (TRAIN[i].y==y)
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
                vector<size_t> randomSame = randomIndex(y2same[TRAIN[i].y].size());
                vector<size_t> randomDiff = randomIndex(y2diff[TRAIN[i].y].size());
                bool collision = false;
                sameClass[i].reserve(Ksample); diffClass[i].reserve(Ksample);
                for (int j=0;j<Ksample;++j){
                    diffClass[i].push_back(y2diff[TRAIN[i].y][randomDiff[j]]);
                    if (i==y2same[TRAIN[i].y][randomSame[j]]){
                        collision=true;
                        continue;
                    }
                    sameClass[i].push_back(y2same[TRAIN[i].y][randomSame[j]]);
                }
                if (collision)
                    sameClass[i].push_back(y2same[TRAIN[i].y][randomSame[Ksample]]);
            }
        }
    }

    void gSpan::nearest(const std::string &mat){
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
                if (TRAIN[i].y!=TRAIN[j].y){//diff class
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

    std::map<unsigned int, unsigned int> gSpan::support_counts(Projected &projected) {
        std::map<unsigned int, unsigned int> counts;

        for (Projected::iterator cur = projected.begin(); cur != projected.end(); ++cur)
            counts[cur->id] += 1;

        return counts;
    }

    unsigned int gSpan::support(Projected &projected) {
        unsigned int oid = 0xffffffff, size = 0;

        for (Projected::iterator cur = projected.begin(); cur != projected.end(); ++cur) {
            if (oid != cur->id)
                ++size;
            oid = cur->id;
        }
        return size;
    }

    void gSpan::mining(vecVec& Xt, PatternSupportList& patternSupportList) {
        Xt.clear();
        patternSupportList.clear();
        visit=0;

        /////////////////////////////////////////////////////////////
        for (map<unsigned int, unsigned int>::iterator it = singleVertexLabelTRAIN.begin () ; it != singleVertexLabelTRAIN.end () ; ++it) {
            if (it->second < minsup)
                continue;
            ++visit;

            DFS_CODE.createSingleVertex(it->first);
            patternSupportList.push_back(pair<DFSCode, int>(DFS_CODE, it->second));
            DFS_CODE.clear();

            vector<bool> x(TRAIN.size(),0);
            for (map<unsigned int, map<unsigned int, unsigned int> >::iterator cur = singleVertexTRAIN.begin(); cur != singleVertexTRAIN.end(); ++cur)
                x[cur->first] = (cur->second[it->first]>0);
            Xt.push_back(x);
        }
        /////////////////////////////////////////////////////////////

        EdgeList edges;
        Projected_map3 root;

        for (unsigned int id = 0; id < TRAIN.size(); ++id) {
            Graph &g = TRAIN[id];
            for (unsigned int from = 0; from < g.size(); ++from) {
                if (get_forward_root(g, g[from], edges)) {
                    for (EdgeList::iterator it = edges.begin(); it != edges.end(); ++it)
                        root[g[from].label][(*it)->elabel][g[(*it)->to].label].push(id, *it, 0);
                }
            }
        }

        for (Projected_iterator3 fromlabel = root.begin(); fromlabel != root.end(); ++fromlabel) {
            for (Projected_iterator2 elabel = fromlabel->second.begin(); elabel != fromlabel->second.end(); ++elabel) {
                for (Projected_iterator1 tolabel = elabel->second.begin(); tolabel != elabel->second.end(); ++tolabel) {
                    /* Build the initial two-node graph.  It will be grown
                     * recursively within project.
                     */
                    DFS_CODE.push(0, 1, fromlabel->first, elabel->first, tolabel->first);
                    __mining(tolabel->second, Xt, patternSupportList);
                    DFS_CODE.pop();
                }
            }
        }
        cout<<"mining:visit = "<<visit<<endl;
    }

    /* Recursive subgraph mining function (similar to subprocedure 1
     * Subgraph_Mining in [Yan2002]).
     */
    void gSpan::__mining(Projected &projected, vecVec& Xt, PatternSupportList& patternSupportList) {
        /* Check if the pattern is frequent enough.
         */
        unsigned int sup = support(projected);
        if (sup < minsup || DFS_CODE.nodeCount() > maxpat)
            return;
        /* The minimal DFS code check is more expensive than the support check,
         * hence it is done now, after checking the support.
         */
        if (!is_min())
            return;

        ++visit;
        patternSupportList.push_back(pair<DFSCode, int>(DFS_CODE, sup));

        vector<bool> x(TRAIN.size(),0);
        for (Projected::iterator cur = projected.begin(); cur != projected.end(); ++cur)
            x[cur->id] = 1;
        Xt.push_back(x);

        /* We just outputted a frequent subgraph.  As it is frequent enough, so
         * might be its (n+1)-extension-graphs, hence we enumerate them all.
         */
        const RMPath &rmpath = DFS_CODE.buildRMPath();
        int minlabel = DFS_CODE[0].fromlabel;
        int maxtoc = DFS_CODE[rmpath[0]].to;

        Projected_map3 new_fwd_root;
        Projected_map2 new_bck_root;
        EdgeList edges;

        /* Enumerate all possible one edge extensions of the current substructure.
         */
        for (unsigned int n = 0; n < projected.size(); ++n) {

            unsigned int id = projected[n].id;
            PDFS *cur = &projected[n];
            History history(TRAIN[id], cur);

            // XXX: do we have to change something here for directed edges?

            // backward
            for (int i = (int) rmpath.size() - 1; i >= 1; --i) {
                Edge *e = get_backward(TRAIN[id], history[rmpath[i]], history[rmpath[0]], history);
                if (e)
                    new_bck_root[DFS_CODE[rmpath[i]].from][e->elabel].push(id, e, cur);
            }

            // pure forward
            // FIXME: here we pass a too large e->to (== history[rmpath[0]]->to
            // into get_forward_pure, such that the assertion fails.
            //
            // The problem is:
            // history[rmpath[0]]->to > TRAIN[id].size()
            if (get_forward_pure(TRAIN[id], history[rmpath[0]], minlabel, history, edges))
                for (EdgeList::iterator it = edges.begin(); it != edges.end(); ++it)
                    new_fwd_root[maxtoc][(*it)->elabel][TRAIN[id][(*it)->to].label].push(id, *it, cur);

            // backtracked forward
            for (int i = 0; i < (int) rmpath.size(); ++i)
                if (get_forward_rmpath(TRAIN[id], history[rmpath[i]], minlabel, history, edges))
                    for (EdgeList::iterator it = edges.begin(); it != edges.end(); ++it)
                        new_fwd_root[DFS_CODE[rmpath[i]].from][(*it)->elabel][TRAIN[id][(*it)->to].label].push(id, *it, cur);
        }

        /* Test all extended substructures.
         */
        // backward
        for (Projected_iterator2 to = new_bck_root.begin(); to != new_bck_root.end(); ++to) {
            for (Projected_iterator1 elabel = to->second.begin(); elabel != to->second.end(); ++elabel) {
                DFS_CODE.push(maxtoc, to->first, -1, elabel->first, -1);
                __mining(elabel->second, Xt, patternSupportList);
                DFS_CODE.pop();
            }
        }

        // forward
        for (Projected_riterator3 from = new_fwd_root.rbegin(); from != new_fwd_root.rend(); ++from) {
            for (Projected_iterator2 elabel = from->second.begin(); elabel != from->second.end(); ++elabel) {
                for (Projected_iterator1 tolabel = elabel->second.begin(); tolabel != elabel->second.end(); ++tolabel) {
                    DFS_CODE.push(from->first, maxtoc + 1, -1, elabel->first, tolabel->first);
                    __mining(tolabel->second, Xt, patternSupportList);
                    DFS_CODE.pop();
                }
            }
        }
        return;
    }
    
    tree<gSpan::Save>::iterator gSpan::createRoot(){
        if (Tree.empty()){
            EdgeList edges;
            tree<Save>::iterator root=Tree.insert(Tree.begin(), Save());
            
            Projected_map3 &Root=root->new_fwd_root;
            for (unsigned int id = 0; id < TRAIN.size(); ++id) {
                Graph &g = TRAIN[id];
                for (unsigned int from = 0; from < g.size(); ++from) {
                    if (get_forward_root(g, g[from], edges)) {
                        for (EdgeList::iterator it = edges.begin(); it != edges.end(); ++it)
                            Root[g[from].label][(*it)->elabel][g[(*it)->to].label].push(id, *it, 0);
                    }
                }
            }
            root->nextCheck=true;

            for (Projected_iterator3 fromlabel = Root.begin(); fromlabel != Root.end(); ++fromlabel) {
                for (Projected_iterator2 elabel = fromlabel->second.begin(); elabel != fromlabel->second.end(); ++elabel) {
                    for (Projected_iterator1 tolabel = elabel->second.begin(); tolabel != elabel->second.end(); ++tolabel) {
                        DFS_CODE.push(0, 1, fromlabel->first, elabel->first, tolabel->first);

                        Projected *projected=&(tolabel->second);
                        vector<bool> x(TRAIN.size(), 0);
                        unsigned int sup = 0;
                        for (Projected::iterator cur = projected->begin(); cur != projected->end(); ++cur)
                            if (!x[cur->id])
                                ++sup, x[cur->id] = 1;
                        if (sup >= minsup && DFS_CODE.nodeCount() <= maxpat && is_min())
                            Tree.append_child(root, Save(x,DFS_CODE,projected));

                        DFS_CODE.pop();
                    }
                }
            }
        }
        return Tree.begin();
    }
    
    void gSpan::createChildren(const tree<gSpan::Save>::iterator &node){
        if (!node->nextCheck){
            const RMPath &rmpath = DFS_CODE.buildRMPath();
            int minlabel = DFS_CODE[0].fromlabel;
            int maxtoc = DFS_CODE[rmpath[0]].to;
            
            EdgeList edges;
            for (unsigned int n = 0; n < node->projected->size(); ++n) {
                unsigned int id = (*node->projected)[n].id;
                PDFS *cur = &(*node->projected)[n];
                History history(TRAIN[id], cur);

                for (int i = (int) rmpath.size() - 1; i >= 1; --i) {
                    Edge *e = get_backward(TRAIN[id], history[rmpath[i]], history[rmpath[0]], history);
                    if (e)
                        node->new_bck_root[DFS_CODE[rmpath[i]].from][e->elabel].push(id, e, cur);
                }

                if (get_forward_pure(TRAIN[id], history[rmpath[0]], minlabel, history, edges))
                    for (EdgeList::iterator it = edges.begin(); it != edges.end(); ++it)
                        node->new_fwd_root[maxtoc][(*it)->elabel][TRAIN[id][(*it)->to].label].push(id, *it, cur);

                for (int i = 0; i < (int) rmpath.size(); ++i)
                    if (get_forward_rmpath(TRAIN[id], history[rmpath[i]], minlabel, history, edges))
                        for (EdgeList::iterator it = edges.begin(); it != edges.end(); ++it)
                            node->new_fwd_root[DFS_CODE[rmpath[i]].from][(*it)->elabel][TRAIN[id][(*it)->to].label].push(id, *it, cur);
            }
            
            node->nextCheck=true;
            for (Projected_iterator2 to = node->new_bck_root.begin(); to != node->new_bck_root.end(); ++to) {
                for (Projected_iterator1 elabel = to->second.begin(); elabel != to->second.end(); ++elabel) {
                    DFS_CODE.push(maxtoc, to->first, -1, elabel->first, -1);

                    Projected *projected=&(elabel->second);
                    vector<bool> x(TRAIN.size(), 0);
                    unsigned int sup = 0;
                    for (Projected::iterator cur = projected->begin(); cur != projected->end(); ++cur)
                        if (!x[cur->id])
                            ++sup, x[cur->id] = 1;
                    if (sup >= minsup && DFS_CODE.nodeCount() <= maxpat && is_min())
                        Tree.append_child(node, Save(x,DFS_CODE,projected));

                    DFS_CODE.pop();
                }
            }

            for (Projected_riterator3 from = node->new_fwd_root.rbegin(); from != node->new_fwd_root.rend(); ++from) {
                for (Projected_iterator2 elabel = from->second.begin(); elabel != from->second.end(); ++elabel) {
                    for (Projected_iterator1 tolabel = elabel->second.begin(); tolabel != elabel->second.end(); ++tolabel) {
                        DFS_CODE.push(from->first, maxtoc + 1, -1, elabel->first, tolabel->first);

                        Projected *projected=&(tolabel->second);
                        vector<bool> x(TRAIN.size(), 0);
                        unsigned int sup = 0;
                        for (Projected::iterator cur = projected->begin(); cur != projected->end(); ++cur)
                            if (!x[cur->id])
                                ++sup, x[cur->id] = 1;
                        if (sup >= minsup && DFS_CODE.nodeCount() <= maxpat && is_min())
                            Tree.append_child(node, Save(x,DFS_CODE,projected));

                        DFS_CODE.pop();
                    }
                }
            }
        }
    }

    double gSpan::getMaxValue(const std::vector< std::vector<double> >& zS, const std::vector< std::vector<double> >& zD) {
        visit=0;
        double maxval = 0;
        /////////////////////////////////////////////////////////////
        for (const Save &node : singleTree) {
            ++visit;

            double dot = 0;
            for (int i = 0; i < TRAIN.size(); ++i) {
                for (int a=0,j=sameClass[i][a];a<Ksample;++a,j=sameClass[i][a]){
                    bool cij = (node.x[i] != node.x[j]);
                    if (cij)
                        dot -= zS[i][a];
                }
                for (int a=0,l=diffClass[i][a];a<Ksample;++a,l=diffClass[i][a]){
                    bool cil = (node.x[i] != node.x[l]);
                    if (cil)
                        dot += zD[i][a];
                }
            }
            maxval = fmax(maxval, dot);
        }
        /////////////////////////////////////////////////////////////

        tree<Save>::iterator root=createRoot();

        for (tree<Save>::sibling_iterator it=Tree.begin(root);it!=Tree.end(root);++it)
            __getMaxValue(it, maxval, zS, zD);

        cout<<"getMaxValue:visit = "<<visit<<endl;
        return maxval;
    }

    void gSpan::__getMaxValue(const tree<Save>::iterator &node, double& maxval, const std::vector< std::vector<double> >& zS, const std::vector< std::vector<double> >& zD) {
        ++visit;

        DFS_CODE=node->dfscode;

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

    void gSpan::safePatternPruning(double lambda_prev, double lambda, const std::vector< std::vector<double> >& zS, const std::vector< std::vector<double> >& zD, double z2, double epsilon_d, vecVec& Xt, std::vector<DFSCode>& patternList) {
        clock_t start=clock();
        Xt.clear();
        patternList.clear();
        visit=0;

        /////////////////////////////////////////////////////////////
        for (Save &node : singleTree) {
            if (node.screeningScore<=lambda)
                continue;
            ++visit;
            
            double dot = 0, h_norm = 0;
            for (int i = 0; i < TRAIN.size(); ++i) {
                for (int a=0,j=sameClass[i][a];a<Ksample;++a,j=sameClass[i][a]){
                    bool cij = (node.x[i] != node.x[j]);
                    if (cij)
                        dot -= zS[i][a], ++h_norm;
                }
                for (int a=0,l=diffClass[i][a];a<Ksample;++a,l=diffClass[i][a]){
                    bool cil = (node.x[i] != node.x[l]);
                    if (cil)
                        dot += zD[i][a], ++h_norm;
                }
            }
            node.screeningScore=lambda_prev*(2*epsilon_d*sqrt(h_norm)+sqrt(h_norm*z2)+dot)/(2*lambda_prev+sqrt(h_norm*z2)-dot);
            if (node.screeningScore > lambda) {
                Xt.push_back(node.x);
                patternList.push_back(node.dfscode);
            }
        }
        /////////////////////////////////////////////////////////////

        tree<Save>::iterator root=createRoot();

        for (tree<Save>::sibling_iterator it=Tree.begin(root);it!=Tree.end(root);++it)
            __safePatternPruning(lambda_prev, lambda, zS, zD, z2, epsilon_d, it, Xt, patternList);
        clock_t end=clock();
        cout<<"safePatternPruning:visit = "<<visit<<endl;
        cout<<"safePatternPruning:time = "<<(double)(end-start)/CLOCKS_PER_SEC<<endl;
    }

    void gSpan::__safePatternPruning(double lambda_prev, double lambda, const std::vector< std::vector<double> >& zS, const std::vector< std::vector<double> >& zD, double z2, double epsilon_d, const tree<Save>::iterator &node, vecVec& Xt, std::vector<DFSCode>& patternList) {
        
        DFS_CODE=node->dfscode;

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
            node->screeningScore = lambda_prev*(2*epsilon_d*sqrt(h_norm)+sqrt(h_norm*z2)+dot)/(2*lambda_prev+sqrt(h_norm*z2)-dot);
            if (node->screeningScore > lambda) {
                Xt.push_back(node->x);
                patternList.push_back(node->dfscode);
            }
        }
        
        createChildren(node);
        
        for (tree<Save>::sibling_iterator it=Tree.begin(node);it!=Tree.end(node);++it)
            __safePatternPruning(lambda_prev, lambda, zS, zD, z2, epsilon_d, it, Xt, patternList);
        return;
    }

    void gSpan::regularizationPath(double L, double U, double eta, int splitNum, int maxloop, double eps, int freq, double R, const string &kernel) {
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
        map<string, double> mSaveDict;
        vecVec Xt;
        vector<DFSCode> pat;
        size_t dim = 0;
        double epsilon_d=0;
        
        double validMaxF=0;
        double correspondTestF;

        vecVec bestActive_Xt, bestTestXt;
        vector<double> bestActive_m;
        double bestLambda(0);

        vector<double> hMax;
        vector< vector<double> > zSMax=zS, zDMax=zD;
        double z2Max=z2;
        
        for (int idx = 1; idx < splitNum; ++idx) {
            double lambda = R * lambda_prev;

            double alpha = 1e-6;
            double primal_prev = DBL_MAX;
            vector<size_t> index;
            double dualMax = -DBL_MAX;

            cout << "********************************" << endl;
            cout << "[" << idx << "] lambda = " << lambda << endl;
            clock_t start = clock();

            cout<<"rrpb_d = "<<sqrt(z2Max)*((lambda_prev-lambda)/(2*lambda_prev))+epsilon_d<<endl;
            safePatternPruning(lambda_prev, lambda, zSMax, zDMax, z2Max, epsilon_d, Xt, pat);
            m.clear();
            for (const DFSCode &it : pat)
                m.push_back(mSaveDict[it.toString()]);
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

                if (dualMax<dual)
                    dualMax = dual, hMax = h, zSMax=zS, zDMax=zD, z2Max=z2;
                
                double gap = primal - dualMax;
                cout << "primal = " << primal << ", dual = " << dualMax << ", gap = " << gap << endl;
                if (gap <= eps * primal) {
                    if (gap<0)
                        gap=0;
                    epsilon_d = 2*sqrt(gap);
                    vector<double> active_m;
                    vector<DFSCode> active_pat;
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
                        mSaveDict[active_pat[k].toString()] = active_m[k];
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

                    break;
                }
                if (loop % freq == 0) {
                    double r = 2*sqrt(gap);
                    cout<<"r_dgb = "<<r<<endl;
                    if (loop != 0) {
                        vector<size_t> newIndex;
                        newIndex.reserve(index.size());
                        for (size_t k : index) {
                            double ub = hMax[k] + r * sqrt(h_norm[k]);
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
                if (loop != 0) {
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
    }

    std::ostream& operator<<(std::ostream& os, const PatternSupportList & patternSupportList) {
        os << "[";
        for (const std::pair<DFSCode, int>& p : patternSupportList)
            os << p << ", ";
        os << "]" << std::endl;
        return os;
    }

    std::ostream& operator<<(std::ostream& os, const std::pair<DFSCode, int>& p) {
        os << "([ " << p.first << " ], " << p.second << ")";
        return os;
    }


    vecVec gSpan::getTestX(const std::vector< DFSCode >& dfscodes) {
        EdgeList edges;
        Projected_map3 root;
        vecVec Xt;

        /////////////////////////////////////////////////////////////
        for (map<unsigned int, unsigned int>::iterator it = singleVertexLabelTEST.begin () ; it != singleVertexLabelTEST.end () ; ++it) {

            DFS_CODE.createSingleVertex(it->first);
            for (int i=Xt.size();i<dfscodes.size();++i ){
                if (dfscodes[i].contain(DFS_CODE) == 0){
                    for (int j=0, count=i - Xt.size(); j < count;++j )
                        Xt.push_back(vector<bool>(TEST.size(), 0));
                    vector<bool> x(TEST.size(),0);
                    for (map<unsigned int, map<unsigned int, unsigned int> >::iterator cur = singleVertexTEST.begin(); cur != singleVertexTEST.end(); ++cur)
                        x[cur->first] = (cur->second[it->first]>0);
                    Xt.push_back(x);
                    break;
                }
            }
            DFS_CODE.clear();
        }
        /////////////////////////////////////////////////////////////

        for (unsigned int id = 0; id < TEST.size(); ++id) {
            Graph &g = TEST[id];
            for (unsigned int from = 0; from < g.size(); ++from) {
                if (get_forward_root(g, g[from], edges)) {
                    for (EdgeList::iterator it = edges.begin(); it != edges.end(); ++it)
                        root[g[from].label][(*it)->elabel][g[(*it)->to].label].push(id, *it, 0);
                }
            }
        }

        for (Projected_iterator3 fromlabel = root.begin(); fromlabel != root.end(); ++fromlabel) {
            for (Projected_iterator2 elabel = fromlabel->second.begin(); elabel != fromlabel->second.end(); ++elabel) {
                for (Projected_iterator1 tolabel = elabel->second.begin(); tolabel != elabel->second.end(); ++tolabel) {
                    DFS_CODE.push(0, 1, fromlabel->first, elabel->first, tolabel->first);
                    __getTestX(tolabel->second, Xt, dfscodes);
                    DFS_CODE.pop();
                }
            }
        }
        for (int j=0, count=dfscodes.size() - Xt.size(); j < count;++j )
            Xt.push_back(vector<bool>(TEST.size(), 0));
        return Xt;
    }

    void gSpan::__getTestX(Projected &projected, vecVec& Xt, const std::vector< DFSCode >& dfscodes) {
        if (DFS_CODE.nodeCount() > maxpat)
            return;
        bool flag=true;
        for (int i=Xt.size();i<dfscodes.size();++i){
            if (dfscodes[i].contain(DFS_CODE) >= 0){
                flag = false;
                break;
            }
        }
        if (flag)
            return;
        if (!is_min())
            return;

        for (int i=Xt.size();i<dfscodes.size();++i ){
            if (dfscodes[i].contain(DFS_CODE) == 0){
                for (int j=0, count=i - Xt.size(); j < count;++j )
                    Xt.push_back(vector<bool>(TEST.size(), 0));
                vector<bool> x(TEST.size(), 0);
                for (Projected::iterator cur = projected.begin(); cur != projected.end(); ++cur) 
                    x[cur->id] = 1;
                Xt.push_back(x);
                break;
            }
        }

        const RMPath &rmpath = DFS_CODE.buildRMPath();
        int minlabel = DFS_CODE[0].fromlabel;
        int maxtoc = DFS_CODE[rmpath[0]].to;

        Projected_map3 new_fwd_root;
        Projected_map2 new_bck_root;
        EdgeList edges;

        for (unsigned int n = 0; n < projected.size(); ++n) {

            unsigned int id = projected[n].id;
            PDFS *cur = &projected[n];
            History history(TEST[id], cur);

            for (int i = (int) rmpath.size() - 1; i >= 1; --i) {
                Edge *e = get_backward(TEST[id], history[rmpath[i]], history[rmpath[0]], history);
                if (e)
                    new_bck_root[DFS_CODE[rmpath[i]].from][e->elabel].push(id, e, cur);
            }

            if (get_forward_pure(TEST[id], history[rmpath[0]], minlabel, history, edges))
                for (EdgeList::iterator it = edges.begin(); it != edges.end(); ++it)
                    new_fwd_root[maxtoc][(*it)->elabel][TEST[id][(*it)->to].label].push(id, *it, cur);

            for (int i = 0; i < (int) rmpath.size(); ++i)
                if (get_forward_rmpath(TEST[id], history[rmpath[i]], minlabel, history, edges))
                    for (EdgeList::iterator it = edges.begin(); it != edges.end(); ++it)
                        new_fwd_root[DFS_CODE[rmpath[i]].from][(*it)->elabel][TEST[id][(*it)->to].label].push(id, *it, cur);
        }

        for (Projected_iterator2 to = new_bck_root.begin(); to != new_bck_root.end(); ++to) {
            for (Projected_iterator1 elabel = to->second.begin(); elabel != to->second.end(); ++elabel) {
                DFS_CODE.push(maxtoc, to->first, -1, elabel->first, -1);
                __getTestX(elabel->second, Xt, dfscodes);
                DFS_CODE.pop();
            }
        }

        for (Projected_riterator3 from = new_fwd_root.rbegin(); from != new_fwd_root.rend(); ++from) {
            for (Projected_iterator2 elabel = from->second.begin(); elabel != from->second.end(); ++elabel) {
                for (Projected_iterator1 tolabel = elabel->second.begin(); tolabel != elabel->second.end(); ++tolabel) {
                    DFS_CODE.push(from->first, maxtoc + 1, -1, elabel->first, tolabel->first);
                    __getTestX(tolabel->second, Xt, dfscodes);
                    DFS_CODE.pop();
                }
            }
        }
        return;
    }

    // distance : test x train
    double gSpan::error(const std::vector< std::vector<double> > &distance, const int K, int from, int to) const{
        double err=0;
        #pragma omp parallel 
        {
            #pragma omp for reduction(+:err)
            for (size_t i=from;i<to;++i){
                vector<size_t> sortIndex = argsort(distance[i]);
                map<int,int> classMap;
                for (size_t k = 0; k < K; ++k)
                    classMap[TRAIN[sortIndex[k]].y]++;
                int max=-1, argmax=0;
                for (const pair<int,int> &p : classMap){
                    if (p.second > max)
                        max = p.second, argmax = p.first;
                }
                err += (TEST[i].y != argmax);
            }
        }
        return err / (to-from);
    }

    // distance : test x train
    double gSpan::f_measure(const std::vector< std::vector<double> > &distance, const int K, int from, int to, bool macro) const{
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
                classMap[TRAIN[sortIndex[k]].y]++;
            int max=-1, argmax=0;
            for (const pair<int,int> &p : classMap){
                if (p.second > max)
                    max = p.second, argmax = p.first;
            }
            int j=0;
            for (int c : classLabel){
                // F-measure
                TP_FN[j] += (TEST[i].y==c);
                TP_FP[j] += (argmax==c);
                TP[j] += (argmax==c && TEST[i].y==c);
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

    vector< vector<double> > gSpan::euclid() {
        visit=0;
        vector<Graph> temp(TEST);
        temp.insert(temp.end(),TRAIN.begin(),TRAIN.end() );
        
        vector< vector<double> > Euclid(testSize);
        for (vector<double>& v : Euclid)
            v.resize(trainSize, 0);

        /////////////////////////////////////////////////////////////
        for (map<unsigned int, unsigned int>::iterator it = singleVertexLabelTRAIN.begin () ; it != singleVertexLabelTRAIN.end () ; ++it) {
            // trainのサポートによる枝刈り
            if (it->second < minsup)
                continue;
            ++visit;

            vector<bool> xTRAIN(TRAIN.size(),0);
            for (map<unsigned int, map<unsigned int, unsigned int> >::iterator cur = singleVertexTRAIN.begin(); cur != singleVertexTRAIN.end(); ++cur)
                xTRAIN[cur->first] = (cur->second[it->first]>0);
            vector<bool> xTEST(TEST.size(),0);
            for (map<unsigned int, map<unsigned int, unsigned int> >::iterator cur = singleVertexTEST.begin(); cur != singleVertexTEST.end(); ++cur)
                xTEST[cur->first] = (cur->second[it->first]>0);

            for (size_t i = 0; i < testSize; ++i)
                for (size_t j = 0; j < trainSize; ++j)
                    Euclid[i][j] += (xTEST[i]!=xTRAIN[j]);
        }
        /////////////////////////////////////////////////////////////

        EdgeList edges;
        Projected_map3 root;

        for (unsigned int id = 0; id < temp.size(); ++id) {
            Graph &g = temp[id];
            for (unsigned int from = 0; from < g.size(); ++from) {
                if (get_forward_root(g, g[from], edges)) {
                    for (EdgeList::iterator it = edges.begin(); it != edges.end(); ++it)
                        root[g[from].label][(*it)->elabel][g[(*it)->to].label].push(id, *it, 0);
                }
            }
        }

        for (Projected_iterator3 fromlabel = root.begin(); fromlabel != root.end(); ++fromlabel) {
            for (Projected_iterator2 elabel = fromlabel->second.begin(); elabel != fromlabel->second.end(); ++elabel) {
                for (Projected_iterator1 tolabel = elabel->second.begin(); tolabel != elabel->second.end(); ++tolabel) {
                    /* Build the initial two-node graph.  It will be grown
                     * recursively within project.
                     */
                    DFS_CODE.push(0, 1, fromlabel->first, elabel->first, tolabel->first);
                    __euclid(tolabel->second, temp, Euclid);
                    DFS_CODE.pop();
                }
            }
        }
        cout<<"euclid:visit = "<<visit<<endl;
        return Euclid;
    }

    void gSpan::__euclid(Projected &projected, vector<Graph> &temp, vector< vector< double > >& Euclid) {
        if (DFS_CODE.nodeCount() > maxpat)
            return;
        if (!is_min())
            return;

        vector<bool> x(temp.size());
        for (Projected::iterator cur = projected.begin(); cur != projected.end(); ++cur)
            x[cur->id] = 1;

        int sup=0;
        // trainのサポートによる枝刈り
        for(int i=0;i<trainSize;++i )
            sup += x[testSize+i];
        if (sup < minsup)
            return;

        ++visit;
        for (size_t i = 0; i < testSize; ++i)
            for (size_t j = 0; j < trainSize; ++j)
                Euclid[i][j] += (x[i]!=x[testSize+j]);

        const RMPath &rmpath = DFS_CODE.buildRMPath();
        int minlabel = DFS_CODE[0].fromlabel;
        int maxtoc = DFS_CODE[rmpath[0]].to;

        Projected_map3 new_fwd_root;
        Projected_map2 new_bck_root;
        EdgeList edges;

        /* Enumerate all possible one edge extensions of the current substructure.
         */
        for (unsigned int n = 0; n < projected.size(); ++n) {

            unsigned int id = projected[n].id;
            PDFS *cur = &projected[n];
            History history(temp[id], cur);

            // XXX: do we have to change something here for directed edges?

            // backward
            for (int i = (int) rmpath.size() - 1; i >= 1; --i) {
                Edge *e = get_backward(temp[id], history[rmpath[i]], history[rmpath[0]], history);
                if (e)
                    new_bck_root[DFS_CODE[rmpath[i]].from][e->elabel].push(id, e, cur);
            }

            // pure forward
            // FIXME: here we pass a too large e->to (== history[rmpath[0]]->to
            // into get_forward_pure, such that the assertion fails.
            //
            // The problem is:
            // history[rmpath[0]]->to > temp[id].size()
            if (get_forward_pure(temp[id], history[rmpath[0]], minlabel, history, edges))
                for (EdgeList::iterator it = edges.begin(); it != edges.end(); ++it)
                    new_fwd_root[maxtoc][(*it)->elabel][temp[id][(*it)->to].label].push(id, *it, cur);

            // backtracked forward
            for (int i = 0; i < (int) rmpath.size(); ++i)
                if (get_forward_rmpath(temp[id], history[rmpath[i]], minlabel, history, edges))
                    for (EdgeList::iterator it = edges.begin(); it != edges.end(); ++it)
                        new_fwd_root[DFS_CODE[rmpath[i]].from][(*it)->elabel][temp[id][(*it)->to].label].push(id, *it, cur);
        }

        /* Test all extended substructures.
         */
        // backward
        for (Projected_iterator2 to = new_bck_root.begin(); to != new_bck_root.end(); ++to) {
            for (Projected_iterator1 elabel = to->second.begin(); elabel != to->second.end(); ++elabel) {
                DFS_CODE.push(maxtoc, to->first, -1, elabel->first, -1);
                __euclid(elabel->second, temp, Euclid);
                DFS_CODE.pop();
            }
        }

        // forward
        for (Projected_riterator3 from = new_fwd_root.rbegin(); from != new_fwd_root.rend(); ++from) {
            for (Projected_iterator2 elabel = from->second.begin(); elabel != from->second.end(); ++elabel) {
                for (Projected_iterator1 tolabel = elabel->second.begin(); tolabel != elabel->second.end(); ++tolabel) {
                    DFS_CODE.push(from->first, maxtoc + 1, -1, elabel->first, tolabel->first);
                    __euclid(tolabel->second, temp, Euclid);
                    DFS_CODE.pop();
                }
            }
        }
        return;
    }

    double gSpan::kernel(const std::string &mat){
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
}

#ifndef GSPAN_H
#define GSPAN_H

#include <iostream>
#include <map>
#include <vector>
#include <set>
#include <algorithm>
#include <iterator>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <cstring>
#include <string>
#include <sstream>
#include <cassert>
#include <cfloat>
#include <numeric>
#include <climits>
#include "tree.hh"

namespace GSPAN {

    template <class T> inline void _swap(T &x, T &y) {
        T z = x;
        x = y;
        y = z;
    }

    struct Edge {
        int from;
        int to;
        int elabel;
        unsigned int id;

        Edge() : from(0), to(0), elabel(0), id(0) {
        };
    };

    class Vertex {
    public:
        typedef std::vector<Edge>::iterator edge_iterator;

        int label;
        std::vector<Edge> edge;

        void push(int from, int to, int elabel) {
            edge.resize(edge.size() + 1);
            edge[edge.size() - 1].from = from;
            edge[edge.size() - 1].to = to;
            edge[edge.size() - 1].elabel = elabel;
            return;
        }
    };

    class Graph : public std::vector<Vertex> {
    private:
        unsigned int edge_size_;
    public:
        typedef std::vector<Vertex>::iterator vertex_iterator;

        Graph(bool _directed) {
            directed = _directed;
        };
        bool directed;

        double y;

        unsigned int edge_size() {
            return edge_size_;
        }

        unsigned int vertex_size() {
            return (unsigned int) size();
        } // wrapper
        void buildEdge();
        std::istream &read(std::istream &, bool ignoreEdgeLabel); // read
        std::ostream &write(std::ostream &); // write

        void check(void);

        Graph() : edge_size_(0), directed(false) {
        };
    };

    class DFS {
    public:
        int from;
        int to;
        int fromlabel;
        int elabel;
        int tolabel;

        friend bool operator==(const DFS &d1, const DFS &d2) {
            return (d1.from == d2.from && d1.to == d2.to && d1.fromlabel == d2.fromlabel
                    && d1.elabel == d2.elabel && d1.tolabel == d2.tolabel);
        }

        friend bool operator!=(const DFS &d1, const DFS &d2) {
            return (!(d1 == d2));
        }

        DFS() : from(0), to(0), fromlabel(0), elabel(0), tolabel(0) {
        };
    };

    typedef std::vector<int> RMPath;

    struct DFSCode : public std::vector <DFS> {
    private:
        RMPath rmpath;
    public:
        /////////////////////////////////////////////////////////////
        const static int NOTHING = -1;
        /////////////////////////////////////////////////////////////
        const RMPath& buildRMPath();

        /* Convert current DFS code into a graph.
         */
        bool toGraph(Graph &);

        /* Return number of nodes in the graph.
         */
        unsigned int nodeCount(void);

        /////////////////////////////////////////////////////////////
        void createSingleVertex(unsigned int label) {
            resize(1);
            DFS &d = (*this)[0];

            d.from = 0;
            d.to = NOTHING;
            d.fromlabel = label;
            d.elabel = NOTHING;
            d.tolabel = NOTHING;
        }
        /////////////////////////////////////////////////////////////

        void push(int from, int to, int fromlabel, int elabel, int tolabel) {
            resize(size() + 1);
            DFS &d = (*this)[size() - 1];

            d.from = from;
            d.to = to;
            d.fromlabel = fromlabel;
            d.elabel = elabel;
            d.tolabel = tolabel;
        }

        void pop() {
            resize(size() - 1);
        }
        int contain(const DFSCode &dfscode)const;

        std::ostream &write(std::ostream &)const; // write
        std::ostream &writeCode(std::ostream &os);
        friend std::ostream& operator<<(std::ostream &os, const DFSCode &dfs);
        std::string toString()const;
    };

    struct PDFS {
        unsigned int id; // ID of the original input graph
        Edge *edge;
        PDFS *prev;

        PDFS() : id(0), edge(0), prev(0) {
        };
    };

    class History : public std::vector<Edge*> {
    private:
        std::vector<int> edge;
        std::vector<int> vertex;

    public:

        bool hasEdge(unsigned int id) {
            return (bool)edge[id];
        }

        bool hasVertex(unsigned int id) {
            return (bool)vertex[id];
        }
        void build(Graph &, PDFS *);

        History() {
        };

        History(Graph& g, PDFS *p) {
            build(g, p);
        }

    };

    class Projected : public std::vector<PDFS> {
    public:

        void push(int id, Edge *edge, PDFS *prev) {
            resize(size() + 1);
            PDFS &d = (*this)[size() - 1];
            d.id = id;
            d.edge = edge;
            d.prev = prev;
        }
    };

    typedef std::vector <Edge*> EdgeList;

    bool get_forward_pure(Graph&, Edge *, int, History&, EdgeList &);
    bool get_forward_rmpath(Graph&, Edge *, int, History&, EdgeList &);
    bool get_forward_root(Graph&, Vertex&, EdgeList &);
    Edge *get_backward(Graph&, Edge *, Edge *, History&);

    typedef std::vector< std::vector<bool> > vecVec;
    typedef std::vector< std::pair< DFSCode, int> > PatternSupportList;

    class gSpan {
    private:

        typedef std::map<int, std::map <int, std::map <int, Projected> > > Projected_map3;
        typedef std::map<int, std::map <int, Projected> > Projected_map2;
        typedef std::map<int, Projected> Projected_map1;
        typedef std::map<int, std::map <int, std::map <int, Projected> > >::iterator Projected_iterator3;
        typedef std::map<int, std::map <int, Projected> >::iterator Projected_iterator2;
        typedef std::map<int, Projected>::iterator Projected_iterator1;
        typedef std::map<int, std::map <int, std::map <int, Projected> > >::reverse_iterator Projected_riterator3;
        
        class Save {
        public:
            std::vector<bool> x;
            DFSCode dfscode;
            Projected *projected=NULL;
            Projected_map3 new_fwd_root;
            Projected_map2 new_bck_root;
            bool nextCheck=false;
            double screeningScore=DBL_MAX;
            double pruningScore=DBL_MAX;
            Save(){}
            Save(const std::vector<bool> &x, const DFSCode &dfscode, Projected *projected){
                this->x=x, this->dfscode=dfscode, this->projected=projected;
            }
        };
        tree<Save> Tree;
        std::vector<Save> singleTree;

        std::vector < Graph > TEST;
        std::vector < Graph > TRAIN;
        DFSCode DFS_CODE;
        DFSCode DFS_CODE_IS_MIN;
        Graph GRAPH_IS_MIN;
        // std::vector < std::vector < int > > Euclid;

        ////////////////////////////////////////////////////////////////////////////////
    public:
        int trainSize ;
        int testSize ;
        int Ksample ;
    private:
        std::set<int> classLabel;
        std::vector< std::vector<int> > sameClass;
        std::vector< std::vector<int> > diffClass;
        std::vector<size_t> trainIndex;
        std::vector<size_t> testIndex;
        ////////////////////////////////////////////////////////////////////////////////

        unsigned int minsup;
        unsigned int maxpat;
        bool directed;
        int visit;
        /////////////////////////////////////////////////////////////
        /* Singular vertex handling stuff
         * [graph][vertexlabel] = count.
         */
        std::map<unsigned int, std::map<unsigned int, unsigned int> > singleVertexTRAIN;
        std::map<unsigned int, unsigned int> singleVertexLabelTRAIN;
        std::map<unsigned int, std::map<unsigned int, unsigned int> > singleVertexTEST;
        std::map<unsigned int, unsigned int> singleVertexLabelTEST;
        /////////////////////////////////////////////////////////////

        bool is_min();
        bool project_is_min(Projected &);

        std::map<unsigned int, unsigned int> support_counts(Projected &projected);
        unsigned int support(Projected&);
        void __getMaxValue(const tree<Save>::iterator &node, double& maxval, const std::vector< std::vector<double> >& zS, const std::vector< std::vector<double> >& zD);
        double __workingSetPruning(double lambda, const std::vector< std::vector<double> >& zS, const std::vector< std::vector<double> >& zD, const tree<Save>::iterator &node, vecVec& Xt, std::vector<DFSCode>& patternList);
        std::istream &read(std::istream &, double train, bool ignoreEdgeLabel, int VertexLabelDegreeQuantization);
        tree<Save>::iterator createRoot();
        void createChildren(const tree<Save>::iterator &node);

        void __getTestX(Projected &projected, vecVec& Xt, const std::vector< DFSCode >& dfscodes) ;

        void setKsample(const std::string &kernel);
        void nearest(const std::string &mat);
    public:

        gSpan(const std::string& filename, unsigned int minsup, unsigned int maxpat, double train, int Ksample, bool ignoreEdgeLabel, int VertexLabelDegreeQuantization, bool directed = false);
        double getMaxValue(const std::vector< std::vector<double> >& zS, const std::vector< std::vector<double> >& zD);
        double workingSetPruning(double lambda, const std::vector< std::vector<double> >& zS, const std::vector< std::vector<double> >& zD, vecVec& Xt, std::vector<DFSCode>& patternList);
        void regularizationPath(double L=1, double U=0.95, double eta=1, int splitNum = 100, int maxloop = 100000, double eps = 1e-6, int freq = 5, double R=0.96, const std::string &kernel="");
        vecVec getTestX(const std::vector< DFSCode >& dfscodes); 
        double error(const std::vector< std::vector<double> > &distance, const int K, int from, int to) const;
        double f_measure(const std::vector< std::vector<double> > &distance, const int K, int from, int to, bool macro=true) const;
        double kernel(const std::string &mat);
    };

    template<class T>
    std::ostream& operator<<(std::ostream& os, const std::vector< std::vector<T> >& mat) {
        os << "[" << std::endl;
        for (const std::vector<T>& v : mat)
            os << v <<std::endl;
        os << "]" << std::endl;
        return os;
    }

    template<class T>
    std::ostream& operator<<(std::ostream& os, const std::vector<T> &vec) {
        os << "[ ";
        for (const T &e : vec)
            os << e << ", ";
        os << "]";
        return os;
    }

    template<class T1, class T2>
    double operator*(const std::vector<T1>& v1, const std::vector<T2> &v2) {
        double dot = 0;
        for (int i = 0, size = v1.size(); i < size; ++i)
            dot += v1[i] * v2[i];
        return dot;
    }

    std::ostream& operator<<(std::ostream& os, const PatternSupportList& patternSupportList);
    std::ostream& operator<<(std::ostream& os, const std::pair<DFSCode, int>& p);
    std::vector<size_t> randomIndex(size_t n);
    template <typename T>
    std::vector<size_t> argsort(const std::vector<T> &v) {
        // initialize original index locations
        std::vector<size_t> idx(v.size());
        std::iota(idx.begin(), idx.end(), 0);
        // sort indexes based on comparing values in v
        std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
        return idx;
    }
    
};

#endif // GSPAN_H

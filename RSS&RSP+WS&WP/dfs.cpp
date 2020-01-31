#include "gspan.h"

using namespace std;

namespace GSPAN {


    bool DFSCode::toGraph(Graph &g) {
        g.clear();

        for (DFSCode::iterator it = begin(); it != end(); ++it) {
            g.resize(std::max(it->from, it->to) + 1);

            if (it->fromlabel != -1)
                g[it->from].label = it->fromlabel;
            if (it->tolabel != -1)
                g[it->to].label = it->tolabel;

            g[it->from].push(it->from, it->to, it->elabel);
            if (!g.directed)
                g[it->to].push(it->to, it->from, it->elabel);
        }
        g.buildEdge();
        return true;
    }

    unsigned int DFSCode::nodeCount(void) {
        unsigned int nodecount = 0;

        for (DFSCode::iterator it = begin(); it != end(); ++it)
            nodecount = std::max(nodecount, (unsigned int) (std::max(it->from, it->to) + 1));

        return nodecount;
    }

    std::ostream &DFSCode::write(std::ostream &os) const {
        if (empty())return os;
        /////////////////////////////////////////////////////////////
        if (size()==1 && (*this)[0].to == NOTHING && (*this)[0].elabel == NOTHING && (*this)[0].tolabel == NOTHING){
            os << "(" << (*this)[0].fromlabel << ")";
            return os;
        }
        /////////////////////////////////////////////////////////////

        os << "(" << (*this)[0].fromlabel << ") " << (*this)[0].elabel << " (0f" << (*this)[0].tolabel << ")";

        for (unsigned int i = 1; i < size(); ++i) {
            if ((*this)[i].from < (*this)[i].to)
                os << " " << (*this)[i].elabel << " (" << (*this)[i].from << "f" << (*this)[i].tolabel << ")";
            else
                os << " " << (*this)[i].elabel << " (b" << (*this)[i].to << ")";
        }
        return os;
    }

    std::ostream &DFSCode::writeCode(std::ostream &os) {
        if (empty()) return os;
        /////////////////////////////////////////////////////////////
        if (size()==1 && (*this)[0].to == NOTHING && (*this)[0].elabel == NOTHING && (*this)[0].tolabel == NOTHING){
            os << "v " << (*this)[0].from << " " << (*this)[0].fromlabel << endl;
            return os;
        }
        /////////////////////////////////////////////////////////////

        os << "v " << (*this)[0].from << " " << (*this)[0].fromlabel << endl
                << "e " << (*this)[0].from << " " << (*this)[0].to << " " << (*this)[0].elabel << endl
                << "v " << (*this)[0].to << " " << (*this)[0].tolabel << endl;
        for (unsigned int i = 1; i < size(); ++i) {
            if ((*this)[i].from < (*this)[i].to) {
                os << "e " << (*this)[i].from << " " << (*this)[i].to << " " << (*this)[i].elabel << endl
                        << "v " << (*this)[i].to << " " << (*this)[i].tolabel << endl;
            } else {
                os << "e " << (*this)[i].from << " " << (*this)[i].to << " " << (*this)[i].elabel << endl;
            }
        }
        os << endl;
        return os;
    }

    std::ostream& operator<<(std::ostream& os, const DFSCode& dfs) {
        dfs.write(os);
        return os;
    }

    std::string DFSCode::toString() const {
        std::stringstream ss;
        this->write(ss);
        return std::string(ss.str());
    }

    // this が 引数を含むかどうか (最小DFSコード)
    // 含む:1, 同じ:0, 含まない:-1
    int DFSCode::contain(const DFSCode &dfscode) const{
        size_t size = dfscode.size();
        if (this->size() < size )
            return -1;
        for (size_t i=0;i<size;++i)
            if (dfscode[i] != (*this)[i])
                return -1;
        if (this->size() == size )
            return 0;
        return 1;
    }
}

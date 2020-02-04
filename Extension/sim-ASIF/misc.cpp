#include "gspan.h"

namespace GSPAN {

    const RMPath &DFSCode::buildRMPath() {
        rmpath.clear();

        int old_from = -1;

        for (int i = size() - 1; i >= 0; --i) {
            if ((*this)[i].from < (*this)[i].to && // forward
                    (rmpath.empty() || old_from == (*this)[i].to)) {
                rmpath.push_back(i);
                old_from = (*this)[i].from;
            }
        }

        return rmpath;
    }

    void History::build(Graph &graph, PDFS *e) {
        // first build history
        clear();
        edge.clear();
        edge.resize(graph.edge_size());
        vertex.clear();
        vertex.resize(graph.size());

        if (e) {
            push_back(e->edge);
            edge[e->edge->id] = vertex[e->edge->from] = vertex[e->edge->to] = 1;

            for (PDFS *p = e->prev; p; p = p->prev) {
                push_back(p->edge); // this line eats 8% of overall instructions(!)
                edge[p->edge->id] = vertex[p->edge->from] = vertex[p->edge->to] = 1;
            }
            std::reverse(begin(), end());
        }
    }

    bool get_forward_rmpath(Graph &graph, Edge *e, int minlabel, History& history, EdgeList &result) {
        result.clear();
        assert(e->to >= 0 && e->to < graph.size());
        assert(e->from >= 0 && e->from < graph.size());
        int tolabel = graph[e->to].label;

        for (Vertex::edge_iterator it = graph[e->from].edge.begin(); it != graph[e->from].edge.end(); ++it) {
            int tolabel2 = graph[it->to].label;
            if (e->to == it->to || minlabel > tolabel2 || history.hasVertex(it->to))
                continue;
            if (e->elabel < it->elabel || (e->elabel == it->elabel && tolabel <= tolabel2))
                result.push_back(&(*it));
        }

        return (!result.empty());
    }

    bool get_forward_pure(Graph &graph, Edge *e, int minlabel, History& history, EdgeList &result) {
        result.clear();

        assert(e->to >= 0 && e->to < graph.size());
        /* Walk all edges leaving from vertex e->to.
         */
        for (Vertex::edge_iterator it = graph[e->to].edge.begin();
                it != graph[e->to].edge.end(); ++it) {
            /* -e-> [e->to] -it-> [it->to] */
            assert(it->to >= 0 && it->to < graph.size());
            if (minlabel > graph[it->to].label || history.hasVertex(it->to))
                continue;
            result.push_back(&(*it));
        }

        return (!result.empty());
    }

    bool get_forward_root(Graph &g, Vertex &v, EdgeList &result) {
        result.clear();
        for (Vertex::edge_iterator it = v.edge.begin(); it != v.edge.end(); ++it) {
            assert(it->to >= 0 && it->to < g.size());
            if (v.label <= g[it->to].label)
                result.push_back(&(*it));
        }
        return (!result.empty());
    }

    Edge *get_backward(Graph &graph, Edge* e1, Edge* e2, History& history) {
        if (e1 == e2)
            return 0;

        assert(e1->from >= 0 && e1->from < graph.size());
        assert(e1->to >= 0 && e1->to < graph.size());
        assert(e2->to >= 0 && e2->to < graph.size());

        for (Vertex::edge_iterator it = graph[e2->to].edge.begin(); it != graph[e2->to].edge.end(); ++it) {
            if (history.hasEdge(it->id))
                continue;

            if ((it->to == e1->from) &&
                    ((e1->elabel < it->elabel) ||
                    (e1->elabel == it->elabel) &&
                    (graph[e1->to].label <= graph[e2->to].label)
                    )) {
                return &(*it);
            }
        }

        return 0;
    }

    std::vector<size_t> randomIndex(size_t n){
        std::vector<size_t> idx(n);
        std::iota(idx.begin(), idx.end(), 0);
        std::random_shuffle(idx.begin(), idx.end());
        return idx;
    }
}


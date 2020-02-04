#include "matchingWeight.hpp"

using namespace std;

// 最小費用流を求める関数
MatchingWeight::COST MatchingWeight::MinCostFlow(Graph &G, int s, int t, FLOW inif) {
    COST dist[G.V];
    int prevv[G.V];
    int preve[G.V];

    COST res = 0;
    FLOW f = inif;
    while (f > 0) {
        fill(dist, dist + G.V, INF);
        dist[s] = 0;
        while (true) {
            bool update = false;
            for (int v = 0; v < G.V; ++v) {
                if (dist[v] == INF) continue;
                for (int i = 0; i < G[v].size(); ++i) {
                    Edge &e = G[v][i];
                    if (e.cap > 0 && dist[e.to] > dist[v] + e.cost+1e-12) {
                        dist[e.to] = dist[v] + e.cost;
                        prevv[e.to] = v;
                        preve[e.to] = i;
                        update = true;
                    }
                }
            }
            if (!update) break;
        }

        if (dist[t] == INF) return INF;

        FLOW d = f;
        for (int v = t; v != s; v = prevv[v]) {
            d = min(d, G[prevv[v]][preve[v]].cap);
        }
        f -= d;
        res += dist[t] * d;
        for (int v = t; v != s; v = prevv[v]) {
            Edge &e = G[prevv[v]][preve[v]];
            Edge &re = G.redge(e);
            e.cap -= d;
            re.cap += d;
        }
    }
    return res;
}


MatchingWeight::COST MatchingWeight::distance(int NUM_WORKER, int NUM_JOB, const vector< tuple<size_t, size_t, double> > &edges){
	
    // グラフの定義 (ノード数を引数に)
    Graph G(NUM_WORKER + NUM_JOB + 2);             // +2 は S, T の分

    // スーパーノード S, T の index
    int S_node = NUM_WORKER + NUM_JOB;
    int T_node = NUM_WORKER + NUM_JOB + 1;
	
	for (const auto &t: edges)
		G.addedge(std::get<0>(t), std::get<1>(t) + NUM_WORKER, 1, std::get<2>(t));
	
    for (int i = 0; i < NUM_WORKER; ++i) {
        // S から従業員 i へと、容量 1, コスト 0 の枝を張る
        G.addedge(S_node, i, 1, 0);
    }

    for (int j = 0; j < NUM_JOB; ++j) {
        // 仕事 j から T へと、容量 1, コスト 0 の枝を張る
        G.addedge(j + NUM_WORKER, T_node, 1, 0);
    }

    // 最小費用流を求める
    // COST res = MinCostFlow(G, S_node, T_node, NUM_JOB);
    COST res = MinCostFlow(G, S_node, T_node, NUM_WORKER);
	return res;
}

#include "matching.hpp"

using namespace std;

void Matching::dibfs(Graph &G, int s, int *level) {
	for (int i = 0; i < G.V; ++i) level[i] = -1;
	level[s] = 0;
	queue<int> que;
	que.push(s);
	while (!que.empty()) {
		int v = que.front();
		que.pop();
		for (int i = 0; i < G[v].size(); ++i) {
			Edge &e = G[v][i];
			if (level[e.to] < 0 && e.cap > 0) {
				level[e.to] = level[v] + 1;
				que.push(e.to);
			}
		}
	}
}

Matching::FLOW Matching::didfs(Graph &G, int v, int t, FLOW f, int *level, int* iter) {
	if (v == t) return f;
	for (int &i = iter[v]; i < G[v].size(); ++i) {
		Edge &e = G[v][i], &re = G.redge(e);
		if (level[v] < level[e.to] && e.cap > 0) {
			FLOW d = didfs(G, e.to, t, min(f, e.cap), level, iter);
			if (d > 0) {
				e.cap -= d;
				re.cap += d;
				return d;
			}
		}
	}
	return 0;
}

// 最大流を求めるメイン関数
Matching::FLOW Matching::Dinic(Graph &G, int s, int t) {
	// 最大流を求めるサブルーチンたち
	int level[G.V];
	int iter[G.V];
	FLOW res = 0;
	while (true) {
		dibfs(G, s, level);
		if (level[t] < 0) return res;
		memset(iter, 0, sizeof(iter));
		FLOW flow;
		while ((flow = didfs(G, s, t, INF, level, iter)) > 0) {
			res += flow;
		}
	}
}

Matching::FLOW Matching::matching(int srcNum, int sinkNum, const std::vector< std::pair<size_t, size_t> > &edges){
	Graph G(srcNum + sinkNum + 2);             // +2 は S, T の分
    // スーパーノード S, T の index
    int S_node = srcNum + sinkNum;
    int T_node = srcNum + sinkNum + 1;

	for (const pair<size_t,size_t> p: edges)
		G.addedge(p.first, p.second + srcNum, 1);
    for (int i = 0; i < srcNum; ++i) {
        // S から男 i へと、容量 1 の枝を張る
        G.addedge(S_node, i, 1);
    }

    for (int j = 0; j < sinkNum; ++j) {
        // 女 j から T へと、容量 1 の枝を張る
        G.addedge(j + srcNum, T_node, 1);
    }

    // 最大流を求める
	return Dinic(G, S_node, T_node);
}

#ifndef MATCHING_H
#define MATCHING_H

#include <iostream>
#include <vector>
#include <map>
#include <queue>
#include <cstring>
#include <climits>

namespace Matching{
	using namespace std;
	typedef int FLOW;                // フローを表す型、今回は int 型
	const FLOW INF = INT_MAX;      // 十分大きい値

	// グラフの辺の構造体
	struct Edge {
		int rev, from, to;
		FLOW cap, icap;
		Edge(int r, int f, int t, FLOW c) : rev(r), from(f), to(t), cap(c), icap(c) {}
		friend ostream& operator << (ostream& s, const Edge& E) {
			if (E.cap > 0) return s << E.from << "->" << E.to << '(' << E.cap << ')';
			else return s;
		}
	};

	// グラフ構造体
	struct Graph {
		int V;
		vector<Edge> *list;

		Graph(int n = 0) : V(n) {list =new vector<Edge>[V]; for (int i = 0; i < V; ++i) list[i].clear(); }
		void init(int n = 0) { V = n; for (int i = 0; i < V; ++i) list[i].clear(); }
		void resize(int n = 0) { V = n; }
		void reset() { for (int i = 0; i < V; ++i) for (int j = 0; j < list[i].size(); ++j) list[i][j].cap = list[i][j].icap; }
		inline vector<Edge>& operator [] (int i) { return list[i]; }
		~Graph(){delete []list;}
		Edge &redge(Edge e) {
			if (e.from != e.to) return list[e.to][e.rev];
			else return list[e.to][e.rev + 1];
		}

		void addedge(int from, int to, FLOW cap) {
			list[from].push_back(Edge((int)list[to].size(), from, to, cap));
			list[to].push_back(Edge((int)list[from].size() - 1, to, from, 0));
		}
	};
	void dibfs(Graph &G, int s, int *level);
	FLOW didfs(Graph &G, int v, int t, FLOW f, int *level, int* iter);
	// 最大流を求めるメイン関数
	FLOW Dinic(Graph &G, int s, int t);
	FLOW matching(int srcNum, int sinkNum, const std::vector< std::pair<size_t, size_t> > &edges);
};

#endif // MATCHING_H

#ifndef MATCHINGWEIGHT_H
#define MATCHINGWEIGHT_H

#include <iostream>
#include <vector>
#include <tuple>
#include <cmath>

namespace MatchingWeight{
	using namespace std;

	typedef int FLOW;                // フローを表す型、今回は int 型
	typedef double COST;                // コストを表す型、今回は int 型
	const COST INF = INFINITY;      // 十分大きい値

	// グラフの辺の構造体
	struct Edge {
		int rev, from, to;
		FLOW cap, icap;
		COST cost;
		Edge(int r, int f, int t, FLOW ca, COST co) : rev(r), from(f), to(t), cap(ca), icap(ca), cost(co) {}
	};

	// グラフ構造体
	struct Graph {
		int V;
		vector<Edge> *list;

		Graph(int n = 0) : V(n) {list=new vector<Edge>[V]; for (int i = 0; i < V; ++i) list[i].clear(); }
		void init(int n = 0) { V = n; for (int i = 0; i < V; ++i) list[i].clear(); }
		void resize(int n = 0) { V = n; }
		void reset() { for (int i = 0; i < V; ++i) for (int j = 0; j < list[i].size(); ++j) list[i][j].cap = list[i][j].icap; }
		inline vector<Edge>& operator [] (int i) { return list[i]; }
		~Graph(){delete []list;}
		Edge &redge(Edge &e) {
			if (e.from != e.to) return list[e.to][e.rev];
			else return list[e.to][e.rev + 1];
		}

		void addedge(int from, int to, FLOW cap, COST cost) {
			list[from].push_back(Edge((int)list[to].size(), from, to, cap, cost));
			list[to].push_back(Edge((int)list[from].size() - 1, to, from, 0, -cost));
		}
	};
	// 最小費用流を求める関数
	COST MinCostFlow(Graph &G, int s, int t, FLOW inif);
	COST distance(int NUM_WORKER, int NUM_JOB, const vector< tuple<size_t, size_t, double> > &edges);
};

#endif // MATCHINGWEIGHT_H

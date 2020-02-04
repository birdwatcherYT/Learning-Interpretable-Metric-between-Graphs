#ifndef WL_HPP
#define WL_HPP

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
#include "gspan.h"
#include "matchingWeight.hpp"

namespace GSPAN{
	class Graph;
}

template<class T> 
int index(const std::vector<T> &v, const T& e, const int from=0){
	for (int i=from;i<v.size();++i){
		if (v[i]==e)
			return i;
	}
	return -1;
}


class WL{
	typedef std::pair<int, std::vector<int> > PairLabel;
	const GSPAN::Graph *graph;
	std::vector< int > WLlabeling; // node ID to PairLabel
	std::vector< PairLabel > labelId2pair; // label id to PairLabel
	std::vector< std::vector< int > > history;
	
public:
	WL(const GSPAN::Graph &g);
	void updateWLlabeling(int loop=1);
	double contains(const WL& pattern, 
			const std::vector< std::vector<double> > &vertexLabelDistance, const std::map<int, int> &vertexLabel2idx, double threshold)const;
};

#endif // WL_HPP

#include "WL.hpp"

WL::WL(const GSPAN::Graph &g){
	graph=&g;
	WLlabeling.resize(g.size() );
	for (int i=0;i<g.size();++i){
		const GSPAN::Vertex &v=g[i];
		PairLabel pairLabel( v.label, std::vector<int>() );
		int idx=index(labelId2pair, pairLabel);
		if (idx<0){
			idx=labelId2pair.size();
			labelId2pair.push_back(pairLabel);
		}
		WLlabeling[i]=idx;
	}
	history.push_back(WLlabeling);
}

// update WLlabeling
void WL::updateWLlabeling(int loop){
	for (int _=0;_<loop;++_){
		std::vector< int > WLlabelingTemp(graph->size()); // node ID to PairLabel (list) while updating
		int from = labelId2pair.size();
		// for each vertex
		for (int i=0;i<graph->size();++i){
			std::vector<int> v;
			v.reserve((*graph)[i].edge.size());
			for ( int j=0; j<(*graph)[i].edge.size(); ++j ){
				if ( (*graph)[i].edge[j].from==i ) {
					v.push_back(WLlabeling[(*graph)[i].edge[j].to]);
				} else if ( (*graph)[i].edge[j].to==i ) {
					v.push_back(WLlabeling[(*graph)[i].edge[j].from]);
				}else{
					std::cerr<< "graph error" <<std::endl;
					std::exit(1);
				}
			}
			PairLabel pairLabel( WLlabeling[i], v );
			int idx=index(labelId2pair, pairLabel, from);
			if (idx<0){
				idx=labelId2pair.size();
				labelId2pair.push_back(pairLabel);
			}
			WLlabelingTemp[i]=idx;
		}
		WLlabeling=WLlabelingTemp;
		history.push_back(WLlabeling);
	}
}

// return (pattern in this)
bool WL::contains(const WL& pattern) const {
	// pattern -> this
	std::vector< std::set<int> > containMap(pattern.labelId2pair.size(), std::set<int>());
	// for each step
	for (size_t h=0;h< history.size(); ++h ){
		
		std::vector< std::pair<size_t,size_t> > edges;
		for (size_t i = 0; i < pattern.history[h].size(); ++i){ //***
			int patLab = pattern.history[h][i]; 
			const PairLabel &patPair = pattern.labelId2pair[patLab];
			for (size_t j=0; j < history[h].size();++j ){ //***
				int thisLab = history[h][j]; 
				const PairLabel &thisPair = labelId2pair[thisLab];
				if (h==0){
					if (patPair.first!=thisPair.first)
						continue;
				}else{
					// patLab in thisLab ----------------------------
					// first
					if (containMap[patPair.first].count(thisPair.first)==0)
						continue;
					// second
					std::vector< std::pair<size_t,size_t> > _edges;
					for (size_t _i=0;_i<patPair.second.size();++_i){
						for (size_t _j=0;_j<thisPair.second.size();++_j){
							if (containMap[patPair.second[_i]].count(thisPair.second[_j])!=0)
								_edges.push_back(std::pair<size_t,size_t>(_i,_j));
						}
					}
					if (Matching::matching(patPair.second.size(), thisPair.second.size(), _edges) < patPair.second.size())
						continue;
					//-----------------------------------------------
				}
				containMap[patLab].insert(thisLab);
				edges.push_back(std::pair<size_t,size_t>(i,j));
			}
			if (containMap[patLab].empty())
				return false;
		}
		if (Matching::matching(pattern.history[h].size(), history[h].size(), edges) < pattern.history[h].size())
			return false;
	}
	return true;
}

#include "dbscan.h"
#include <pcl/kdtree/kdtree_flann.h>

namespace EXX{

	template <typename T>
	void dbscan::cluster( flann::Matrix<T> &features, double eps, int minPts, std::vector<std::vector<int> > &c ){

		if (c.size() > 0){ c.clear(); }

		// Build flann nearest neighbour search.
		std::vector< std::vector<int> > indices;
        std::vector<std::vector<T> > dists;
        std::vector< std::vector<int> > indices_i;
        std::vector<std::vector<T> > dists_i;
		flann::Index<flann::L2<T> > index(features, flann::KDTreeIndexParams(4));
		flann::Matrix<float> current(new float[1*features.rows], 1, features.rows);
        index.buildIndex();
		float radius = 1.0;
        // index.radiusSearch(features, indices, dists, radius, flann::SearchParams(128));

		std::vector<bool> visited(features.rows, false);
		std::vector<bool> cluster_member(features.rows, false);
		std::vector<int> cluster;

		int t = 0;
		int numPts = 0;
		for (size_t i = 0; i < features.rows; ++i ){
			if ( !visited.at(i) ){ 
				visited.at(i) = true;

				for (size_t j = 0; j < features.cols; ++j){
					current[0][j] = features[i][j];
				}

				indices.clear();
				dists.clear();
				numPts = index.radiusSearch(current, indices, dists, radius, flann::SearchParams(128));

				if( numPts >= minPts ){
					
					cluster.clear();
					cluster.push_back(i);


					for (size_t k = 0; k < indices[0].size(); ++k ){
					// for ( auto k : indices[0] ){
						t = indices[0].at(k);
						if ( !visited.at(t) ){
							visited.at(t) = true;

							for (size_t j = 0; j < features.cols; ++j){
								current[0][j] = features[t][j];
							}

							if ( index.radiusSearch(current, indices, dists, radius, flann::SearchParams(128)) > minPts ){
								indices[0].push_back(t);
							}
						}
						if ( !cluster_member.at(t) ){
							cluster_member.at(t) = true;
							cluster.push_back(t);
						}
					}
				}
				if ( cluster.size() > 0 ){
					c.push_back(cluster);
				}

			}
		}
	}
}
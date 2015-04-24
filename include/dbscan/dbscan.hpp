template <typename T>
void EXX::dbscan::cluster( flann::Matrix<T> &features, std::set<int> walls, std::set<int> floors, T eps, int minPts, std::vector<std::vector<int> > &c ){

	if (c.size() > 0){ c.clear(); }

	int nn = 10;

	// Build flann nearest neighbour search.
	// std::vector< std::vector<int> > indices(1);
    // std::vector<std::vector<T> > dists(1);
    std::vector< std::vector<int> > indices_i(1);
    std::vector<std::vector<T> > dists_i(1);
	flann::Index<flann::L2<T> > index(features, flann::KDTreeIndexParams(4));
	flann::Matrix<T> current(new T[1*features.cols], 1, features.cols);
	flann::Matrix<T> dists(new T[current.rows*nn], current.rows, nn);
    flann::Matrix<int> indices(new int[current.rows*nn], current.rows, nn);
    index.buildIndex();
    // index.radiusSearch(features, indices, dists, radius, flann::SearchParams(128));

	std::vector<bool> visited(features.rows, false);
	std::vector<bool> cluster_member(features.rows, false);
	std::vector<int> neighborPts;
	std::vector<int> cluster;

	int t = 0;
	int numPts = 0;
	for (size_t i = 0; i < features.rows; ++i ){

		if (walls.find(i) != walls.end() || floors.find(i) != floors.end()){
			visited.at(i) = true;
			cluster_member.at(i) = true;
			continue;
		}

		if ( !visited.at(i) ){ 
			visited.at(i) = true;
			cluster.clear();
			cluster.push_back(i);

			for (size_t j = 0; j < features.cols; ++j){
				current[0][j] = features[i][j];
			}

			numPts = index.radiusSearch(current, indices, dists, eps, flann::SearchParams(128));
			if( numPts >= minPts ){
				
				neighborPts.clear();

				for (size_t k = 0; k < indices.cols; ++k){
					if ( indices[0][k] == -1 ){ break; }
					neighborPts.push_back(indices[0][k]);
				}

				for (size_t k = 0; k < neighborPts.size(); ++k ){
					std::cout << "inni" << std::endl;
					t = neighborPts.at(k);
					if ( !visited.at(t) ){
						visited.at(t) = true;
						for (size_t j = 0; j < features.cols; ++j){
							current[0][j] = features[t][j];
						}

						if ( index.radiusSearch(current, indices, dists, eps, flann::SearchParams(128)) > minPts ){
							for (size_t k = 0; k < indices.cols; ++k){
								if ( indices[0][k] == -1 ){ break; }
								neighborPts.push_back(indices[0][k]);
							}
						}
					}
					if ( !cluster_member.at(t) ){
						cluster_member.at(t) = true;
						cluster.push_back(t);
					}
				}
			}
			c.push_back(cluster);

		}
	}
}
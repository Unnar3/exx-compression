#ifndef DBSCAN_H
#define DBSCAN_H

#include <pcl/kdtree/kdtree_flann.h>
#include <vector>

namespace EXX{
class dbscan{
public:
	template <typename T>
	void cluster( flann::Matrix<T> &features, std::set<int> walls, std::set<int> floors, T eps, int minPts, std::vector<std::vector<int> > &c);
};
}
#include "dbscan.hpp"
#endif
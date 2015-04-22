#ifndef DBSCAN_H
#define DBSCAN_H

#include <pcl/kdtree/kdtree_flann.h>
#include <vector>

namespace EXX{
class dbscan{
public:
	template <typename T>
	void cluster( flann::Matrix<T> &features, double eps, int minPts, std::vector<std::vector<int> > &c);
};
}
#endif
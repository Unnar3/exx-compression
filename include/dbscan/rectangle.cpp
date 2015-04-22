#include "rectangle.h"
#include <pcl/kdtree/kdtree_flann.h>
void Rectangle::set_values (int x, int y) {
  width = x;
  height = y;
}
void Rectangle::cluster( flann::Matrix<double> &features, double eps, int minPts, std::vector<std::vector<int> > &co ){
	std::cout << "hahah" << std::endl;

	if (co.size() > 0){ co.clear(); }

	// Build flann nearest neighbour search.
	std::vector< std::vector<int> > indices;
    std::vector<std::vector<double> > dists;
    std::vector< std::vector<int> > indices_i;
    std::vector<std::vector<double> > dists_i;
}
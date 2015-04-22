#include <pcl/kdtree/kdtree_flann.h>
class Rectangle {
    int width, height;
  public:
    void set_values (int,int);
    int area() {return width*height;}

    template <typename T>
	void cluster( flann::Matrix<T> &features, double eps, int minPts, std::vector<std::vector<int> > &c);
};
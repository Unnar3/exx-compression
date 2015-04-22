#include <pcl/kdtree/kdtree_flann.h>
class Rectangle {
    int width, height;
  public:
    void set_values (int,int);
    int area() {return width*height;}

    void cluster( flann::Matrix<double> &features, double eps, int minPts, std::vector<std::vector<int> > &co);
};
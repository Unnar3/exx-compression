#ifndef PLANE_FEATURES_H
#define PLANE_FEATURES_H

#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/ModelCoefficients.h>
#include <pcl/surface/gp3.h>

#include <moment_of_inertia/moment_of_inertia_estimation.h>
#include <Eigen/Dense>
#include <vector>

typedef pcl::PointXYZRGB PointT;
typedef pcl::PointCloud<PointT> PointCloudT;
typedef pcl::PointXYZRGBA PointTA;
typedef pcl::PointCloud<PointTA> PointCloudTA;
typedef pcl::ModelCoefficients ModelCoeffT;
typedef std::vector<PointCloudT::Ptr> vPointCloudT;
typedef std::vector<PointCloudTA::Ptr> vPointCloudTA;


namespace EXX{

struct planeDescriptor{
    std::vector<float> momentOfInertia;
    double boundingBoxArea;
    double hullArea;
    double hullBoundingBoxRatio;
    double widthLengthRatio;
    double distFromCentroid;
    Eigen::Vector4d normal;
};

class planeFeatures{

public:
    planeFeatures(){};
    ~planeFeatures(){};

    // Calculate features for each plane in planes.
    static void calculateFeatures(vPointCloudT planes, vPointCloudT hulls, std::vector<Eigen::Vector4d> normal, std::vector<planeDescriptor> *vPlaneDescriptor);
    static void matchSimilarFeatures(std::vector<planeDescriptor> descriptor, std::vector<std::set<int> > *sets);

private:
	// Takes in two vectors, returns angle between them in range 0 to 1  
    // where 1 means no difference, pi/2 is considered maximum angle and pi as 0.
    static double angleBetweenVectors(Eigen::Vector4d a, Eigen::Vector4d b);

    // Returns RGB value from red to green depending on value, low value results in blue,
    // high value results in red.
    static void getValueBetweenTwoFixedColors(double value, int *red, int *green, int *blue);

    // Takes in min and max points for a bounding box and returns the biggest area
    // and the ratio of length and with for the area.
    static void getBiggestCubeArea(PointT minPoint, PointT maxPoint, double *area, double *WLRatio);

};

}
#endif
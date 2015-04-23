#ifndef PLANE_FEATURES_H
#define PLANE_FEATURES_H

#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/ModelCoefficients.h>
#include <pcl/surface/gp3.h>
#include <pcl/visualization/cloud_viewer.h>

#include <moment_of_inertia/moment_of_inertia_estimation.h>
#include "dbscan/dbscan.h"
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
    Eigen::Vector3f massCenter;
    double boundingBoxArea;
    double hullArea;
    double hullBoundingBoxRatio;
    double widthLengthRatio;
    double distFromCentroid;
    Eigen::Vector4d normal;
};

struct featureSet{
    flann::Matrix<double> features;
    flann::Matrix<int> indices;
    std::set<std::set<int> > objects;
    std::set<int> walls;
    std::set<int> floors;
    std::vector<Eigen::Vector3f> position;
    int objectSize;
};

class planeFeatures{
    boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer;
    bool viewerIsSet;
    double bigPlaneSize;
    double centroid_to_wall_distance_;

public:
    planeFeatures() : 
        viewerIsSet(false), 
        bigPlaneSize(0.8),
        centroid_to_wall_distance_(0.35f) {};
    ~planeFeatures(){};

    // Calculate features for each plane in planes.
    void loadFeatures(const vPointCloudT &planes, const std::vector<Eigen::Vector4d> &normal, const std::vector<int> &normalInd, featureSet &fSet);
    // Nearest neighbour search using flann library
    void matchFeatures(featureSet &fSet, std::vector<std::vector<int> > &c);
    void matchFeatures(featureSet &fSet);
    // Create clusters from nearest neighbour search and load into sets
    void groupFeatures(featureSet &fSet);
    // Detects wall segments that were wrongly classified.
    void improveWalls(const vPointCloudT &planes, featureSet &fSet);
    void improveWalls(const vPointCloudT &planes, featureSet &fSet, std::vector<std::vector<int> > &c);
    // Find repeating objects from set of planes
    void groupObjects(featureSet &fSet);

    void calculateFeatures(vPointCloudT planes, vPointCloudT hulls, std::vector<Eigen::Vector4d> normal, std::vector<int> normalInd, std::vector<planeDescriptor> *vPlaneDescriptor);
    void matchSimilarFeatures(std::vector<planeDescriptor> descriptor, std::vector<std::set<int> > *sets);
    void findRepeatingObjects(vPointCloudT planes, std::vector<planeDescriptor> desc, std::vector<std::set<int> > sets, std::vector< std::set<int>> *objects);    
    void setViewer(boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer);

    //Setter methods
    void setCentroidToWallDistance( float dist ){ centroid_to_wall_distance_ = dist; }
    
private:
	// Takes in two vectors, returns angle between them in range 0 to 1  
    // where 1 means no difference, pi/2 is considered maximum angle and pi as 0.
    static double angleBetweenVectors(Eigen::Vector4d a, Eigen::Vector4d b);

    // Takes in plane, returns angle between vector and floor  
    static double angleToFloor(Eigen::Vector4d a);

    // Returns RGB value from red to green depending on value, low value results in blue,
    // high value results in red.
    static void getValueBetweenTwoFixedColors(double value, int *red, int *green, int *blue);

    // Takes in min and max points for a bounding box and returns the biggest area
    // and the ratio of length and with for the area.
    static void getBiggestCubeArea(PointT minPoint, PointT maxPoint, double *area, double *WLRatio);

    // Returns the euclidean distance between two points.
    static double euclideanDistance(Eigen::Vector3f a, Eigen::Vector3f b);

    // Prints the planeDescriptor in a good way.
    static void printDescriptor(planeDescriptor descr);

    static float calculatePolygonAreaMine (const PointCloudT &polygon);
    static float calculatePolygonAreaMineInv (const PointCloudT &polygon);

};

}
#endif
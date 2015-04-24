#ifndef COMPRESSION_H
#define COMPRESSION_H

#include <moment_of_inertia/moment_of_inertia_estimation.h>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/ModelCoefficients.h>
#include <pcl/surface/gp3.h>
#include <Eigen/Dense>

typedef pcl::PointXYZRGB PointT;
typedef pcl::PointCloud<PointT> PointCloudT;
typedef pcl::PointXYZRGBA PointTA;
typedef pcl::PointCloud<PointTA> PointCloudTA;
typedef pcl::ModelCoefficients ModelCoeffT;
typedef std::vector<PointCloudT::Ptr> vPointCloudT;
typedef std::vector<PointCloudTA::Ptr> vPointCloudTA;

namespace EXX{

struct planesAndCoeffs{
	std::vector<PointCloudT::Ptr> cloud;
	std::vector<ModelCoeffT::Ptr> coeff;
};

struct cloudMesh{
	PointCloudT::Ptr cloud;
	pcl::PolygonMesh mesh;
};

struct densityDescriptor{
	float voxel_res;
	float seed_res;
	float rw_max_dist;
	float gp3_search_rad;
};

class compression{

	// Filtering
	float v_leaf_size_;
	float sv_voxel_res_, sv_seed_res_, sv_color_imp_, sv_spatial_imp_;
	// RANSAC
	int max_ite_, min_inliers_, max_number_planes_; float dist_thresh_;
	// EUCLIDEAN CLUSTERING
	float ec_cluster_tolerance_; int ec_min_cluster_size_;
	// CONCAVE HULLS
	float hull_alpha_, rw_hull_eps_, rw_hull_max_dist_;
	bool simplify_hulls_;
	// GREEDY PROJECTION TRIANGULATION
	float gp3_search_rad_, gp3_Mu_, gp3_max_surface_angle_, gp3_min_angle_, gp3_max_angle_;
	int gp3_max_nearest_neighbours_, gp3_Ksearch_;

	PointCloudT::Ptr cloud_;
	std::vector<PointCloudT::Ptr > planes_;
	std::vector<PointCloudT::Ptr > c_planes_;
	std::vector<PointCloudT::Ptr > sv_planes_;
	std::vector<PointCloudT::Ptr > hulls_;
	std::vector<PointCloudT::Ptr > rw_hulls_;
	std::vector<pcl::ModelCoefficients::Ptr> coeffs_;
	std::vector<cloudMesh > cloud_mesh_;

public:
	compression() : 
		v_leaf_size_(0.02f),
		sv_voxel_res_(0.1f), 
		sv_seed_res_(0.3f), 
		sv_color_imp_(0.5f), 
		sv_spatial_imp_(0.1f),
		max_ite_(100), 
		min_inliers_(200),
		max_number_planes_(100), 
		dist_thresh_(0.04f), 
		ec_cluster_tolerance_(0.05f),
		ec_min_cluster_size_(100),
		hull_alpha_(0.1f), 
		rw_hull_eps_(0.04f),
		rw_hull_max_dist_(0.3f), 
		simplify_hulls_(true),
		gp3_search_rad_(0.3f), 
		gp3_Mu_(3.0f), 
		gp3_max_surface_angle_(M_PI/4.0f),
		gp3_min_angle_(M_PI/20),  
		gp3_max_angle_(2.0f*M_PI/2.5f),
		gp3_max_nearest_neighbours_(100), 
		gp3_Ksearch_(20)
	{ }
	~compression(){ };

	void triangulate();
	void triangulatePlanes();
	
	// Filtering Methods 
	void voxelGridFilter(PointCloudT::Ptr cloud, PointCloudT::Ptr out_cloud);
	void superVoxelClustering(vPointCloudT *cloud, vPointCloudT *out_vec, std::vector<densityDescriptor> &dDesc);

	// RANSAC METHODS
	void extractPlanesRANSAC(PointCloudT::Ptr cloud, planesAndCoeffs *pac);
	void projectToPlane(planesAndCoeffs *pac);
	static void projectToPlaneS(PointCloudT::Ptr cloud, Eigen::Vector4d coeff);
	static void projectToPlaneS(PointCloudT::Ptr cloud, ModelCoeffT::Ptr coeff);

	// EUCLIDIAN CLUSTERING OF PLANES
	void euclideanClusterPlanes(vPointCloudT* cloud, vPointCloudT* out_vec, std::vector<int> *normalIndex);

	// CONCAVE HULLS
	void planeToConcaveHull(vPointCloudT *planes, vPointCloudT *hulls);
	void planeToConvexHull(vPointCloudT &planes, vPointCloudT &hulls, std::vector<double> &area);
	void reumannWitkamLineSimplification(vPointCloudT* hulls, vPointCloudT* s_hulls, std::vector<densityDescriptor> &dDesc);

	void getPlaneDensity( vPointCloudT &planes, vPointCloudT &hulls, std::vector<densityDescriptor> &dDesc);

	// TRIANGULATION
	void greedyProjectionTriangulation(PointCloudT::Ptr nonPlanar, vPointCloudT *planes, vPointCloudT *hulls, std::vector<cloudMesh> *cm);
	void greedyProjectionTriangulationPlanes(PointCloudT::Ptr nonPlanar, vPointCloudT *planes, vPointCloudT *hulls, std::vector<cloudMesh> *cm,std::vector<densityDescriptor> &dDesc);
	void improveTriangulation(std::vector<cloudMesh> &cm, vPointCloudT &planes, vPointCloudT &hulls);

	// SET METHODS
	void setVoxelLeafSize(float leaf){ v_leaf_size_ = leaf; }
	void setSVVoxelResolution(float res){ sv_voxel_res_ = res; }
	void setSVSeedResolution(float res){ sv_seed_res_ = res; }
	void setSVColorImportance(float imp){ sv_color_imp_ = imp; }
	void setSVSpatialImportance(float imp){ sv_spatial_imp_ = imp; }
	void setRANSACMaxIteration(int ite){ max_ite_ = ite; }
	void setRANSACMinInliers(int inl){ min_inliers_ = inl; }
	void setRANSACDistanceThreshold(double thresh){ dist_thresh_ = thresh; }
	void setECClusterTolerance(double dist){ ec_cluster_tolerance_ = dist; }
	void setECMinClusterSize(int size){ ec_min_cluster_size_ = size; }
	void setHULLAlpha(double alpha){ hull_alpha_ = alpha; }
	void setRWHUllEps(double eps){ rw_hull_eps_ = eps; }
	void setRWHullMaxDist(double dist){ rw_hull_max_dist_ = dist; }
	void setSimplifyHulls(bool simp){ simplify_hulls_ = simp; }
	void setGP3SearchRad(double rad){ gp3_search_rad_ = rad; }
	void setGP3Mu(double mu){ gp3_Mu_ = mu; }
	void setGP3MaxSurfaceAngle(double angle){ gp3_max_surface_angle_ = angle; }
	void setGP3MinAngle(double angle){ gp3_min_angle_ = angle; }
	void setGP3MaxAngle(double angle){ gp3_max_angle_ = angle; }
	void setGP3MaxNearestNeighbours(int K){ gp3_max_nearest_neighbours_ = K; }
	void setGP3Ksearch(int K){ gp3_Ksearch_ = K; }


	// RETURN METHODS
	PointCloudT::Ptr returnCloud() { return cloud_; }
	std::vector<PointCloudT::Ptr > returnPlanes() { return planes_; }
	std::vector<PointCloudT::Ptr > returnECPlanes() { return c_planes_; }
	std::vector<PointCloudT::Ptr > returnSuperVoxelPlanes() { return sv_planes_; }
	std::vector<PointCloudT::Ptr > returnHulls() { return hulls_; }
	std::vector<PointCloudT::Ptr > returnRWHulls() { return rw_hulls_; }
	std::vector<cloudMesh > returnCloudMesh() { return cloud_mesh_; }

	// SAVE METHODS
	void saveCloud(std::string path="./", std::string name = "cmprs_cloud"){ savePCD(cloud_, path+name); }
	void savePlanes(std::string name = "cmprs_planes"){ savePCD(planes_, name); }
	// void saveSVPlanes(std::string name = "cmprs_sv_cloud"){ savePCD(sv_planes_, name); }
	void saveHulls(std::string name = "cmprs_hulls"){ savePCD(hulls_, name); }
	void saveRWHulls(std::string name = "cmprs_rw_hulls"){ savePCD(rw_hulls_, name); }
private:

	PointCloudT::Ptr superVoxelClustering_s(PointCloudT::Ptr cloud, densityDescriptor &dDesc);
	double pointToLineDistance(PointT current, PointT next, PointT nextCheck);
	double distBetweenPoints(PointT a, PointT b);
	PointCloudT::Ptr planeToConcaveHull_s(PointCloudT::Ptr cloud);
	void planeToConvexHull_s(const PointCloudT::Ptr cloud, PointCloudT::Ptr out, double &area);
	PointCloudT::Ptr reumannWitkamLineSimplification_s(PointCloudT::Ptr cloud, densityDescriptor &dDesc);
	void savePCD(PointCloudT::Ptr cloud, std::string name);
	void savePCD(std::vector<PointCloudT::Ptr> cloud, std::string name);
	void saveVTK(pcl::PolygonMesh mesh, std::string name);
	cloudMesh greedyProjectionTriangulation_s(PointCloudT::Ptr cloud, float gp3_rad);
	PointCloudT::Ptr PointRGBAtoRGB( PointCloudTA::Ptr cloudRGBA );

};

}
#endif
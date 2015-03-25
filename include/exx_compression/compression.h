
#ifndef COMPRESSION_H
#define COMPRESSION_H

#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/ModelCoefficients.h>
#include <pcl/surface/gp3.h>

typedef pcl::PointXYZRGB PointT;
typedef pcl::PointCloud<PointT> PointCloudT;

namespace EXX{

struct cloudMesh{
	PointCloudT::Ptr cloud;
	pcl::PolygonMesh mesh;
};

class compression{

	// Filtering
	float v_leaf_size_;
	float sv_voxel_res_, sv_seed_res_, sv_color_imp_, sv_spatial_imp_;
	// RANSAC
	int max_ite_, min_inliers_; double dist_thresh_;
	// CONCAVE HULLS
	double hull_alpha_, rw_hull_eps_;
	// GREEDY PROJECTION TRIANGULATION
	double gp3_search_rad_, gp3_Mu_, gp3_max_surface_angle_, gp3_min_angle_, gp3_max_angle_;
	int gp3_max_nearest_neighbours_, gp3_Ksearch_;

	PointCloudT::Ptr cloud_;
	std::vector<PointCloudT::Ptr > planes_;
	std::vector<PointCloudT::Ptr > sv_planes_;
	std::vector<PointCloudT::Ptr > hulls_;
	std::vector<PointCloudT::Ptr > rw_hulls_;
	std::vector<pcl::ModelCoefficients::Ptr> coeffs_;
	std::vector<pcl::PolygonMesh> meshes_;
	std::vector<cloudMesh > cloud_mesh_;
	pcl::PolygonMesh mesh_;

public:
	compression() : 
		v_leaf_size_(0.02f),
		sv_voxel_res_(0.1f), sv_seed_res_(0.3f), sv_color_imp_(0.5f), sv_spatial_imp_(0.1f),
		max_ite_(100), min_inliers_(200), dist_thresh_(0.04),
		hull_alpha_(0.1), rw_hull_eps_(0.02),
		gp3_search_rad_(0.3), gp3_Mu_(3.0), gp3_max_surface_angle_(M_PI/4),
		gp3_min_angle_(M_PI/20),  gp3_max_angle_(2*M_PI/2.5),
		gp3_max_nearest_neighbours_(100), gp3_Ksearch_(20)
	{ }
	~compression(){ };


	void setInputCloud(PointCloudT::Ptr inputCloud);

	void triangulate();
	void triangulatePlanes();
	
	// Filtering Methods 
	void voxelGridFilter();
	void voxelGridFilter(float leaf_size);
	void superVoxelClustering();
	void superVoxelClustering(float voxel_res, float seed_res, float color_imp, float spatial_imp);

	// RANSAC METHODS
	void extractPlanesRANSAC();
	void extractPlanesRANSAC(int max_iterations, int min_inliers, double distance_threshold);

	// CONCAVE HULLS (DEPENDS ON THAT PLANES HAVE BEEN IDENTIFIEDL)
	void planeToConcaveHull();
	void planeToConcaveHull(double alpha);
	void reumannWitkamLineSimplification();
	void reumannWitkamLineSimplification(double eps);

	// TRIANGULATION
	void greedyProjectionTriangulation();
	void greedyProjectionTriangulationPlanes();

	// SET METHODS
	void setVoxelLeafSize(float leaf){ v_leaf_size_ = leaf; }
	void setSVVoxelResolution(float res){ sv_voxel_res_ = res; }
	void setSVSeedResolution(float res){ sv_seed_res_ = res; }
	void setSVColorImportance(float imp){ sv_color_imp_ = imp; }
	void setSVSpatialImportance(float imp){ sv_spatial_imp_ = imp; }
	void setRANSACMaxIteration(int ite){ max_ite_ = ite; }
	void setRANSACMinInliers(int inl){ min_inliers_ = inl; }
	void setRANSACDistanceThreshold(double thresh){ dist_thresh_ = thresh; }
	void setHULLAlpha(double alpha){ hull_alpha_ = alpha; }
	void setRWHUllEps(double eps){ rw_hull_eps_ = eps; }
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
	std::vector<PointCloudT::Ptr > returnSuperVoxelPlanes() { return sv_planes_; }
	std::vector<PointCloudT::Ptr > returnHulls() { return hulls_; }
	std::vector<PointCloudT::Ptr > returnRWHulls() { return rw_hulls_; }
	std::vector<pcl::PolygonMesh > returnMeshes() { return meshes_; }
	pcl::PolygonMesh returnMesh() { return mesh_; }
	std::vector<cloudMesh > returnCloudMesh() { return cloud_mesh_; }

	// SAVE METHODS
	void saveCloud(std::string path="./", std::string name = "cmprs_cloud"){ savePCD(cloud_, path+name); }
	void savePlanes(std::string name = "cmprs_planes"){ savePCD(planes_, name); }
	void saveSVPlanes(std::string name = "cmprs_sv_cloud"){ savePCD(sv_planes_, name); }
	void saveHulls(std::string name = "cmprs_hulls"){ savePCD(hulls_, name); }
	void saveRWHulls(std::string name = "cmprs_rw_hulls"){ savePCD(rw_hulls_, name); }
	void saveMesh(std::string path="./", std::string name="cmprs_mesh"){ saveVTK(mesh_, path+name); }
private:

	PointCloudT::Ptr superVoxelClustering_s(PointCloudT::Ptr cloud, float voxel_res, float seed_res, float color_imp, float spatial_imp);
	void projectToPlane();
	double pointToLineDistance(PointT current, PointT next, PointT nextCheck);
	PointCloudT::Ptr planeToConcaveHull_s(PointCloudT::Ptr cloud, double alpha);
	PointCloudT::Ptr reumannWitkamLineSimplification_s(PointCloudT::Ptr cloud, double eps);
	void savePCD(PointCloudT::Ptr cloud, std::string name);
	void savePCD(std::vector<PointCloudT::Ptr> cloud, std::string name);
	void saveVTK(pcl::PolygonMesh mesh, std::string name);
	cloudMesh greedyProjectionTriangulation_s(PointCloudT::Ptr cloud);

};

}
#endif
#include <plane_features/plane_features.h>
#include <moment_of_inertia/moment_of_inertia_estimation.h>
#include <exx_compression/compression.h>
#include <pcl/common/common.h>
#include <pcl/kdtree/kdtree_flann.h>


namespace EXX{


    void planeFeatures::loadFeatures(const vPointCloudT &planes, const std::vector<Eigen::Vector4d> &normal, const std::vector<int> &normalInd, featureSet &fSet){
        
        int nof = 4;
        // features = flann::Matrix<double>(new double[planes.size()*nof], planes.size(), nof );
        fSet.features = flann::Matrix<double>(new double[planes.size()*nof], planes.size(), nof );
        pcl::MomentOfInertiaEstimation<PointT> feature_extractor;
        PointT min_point_OBB;
        PointT max_point_OBB;
        PointT position_OBB;
        Eigen::Matrix3f rotational_matrix_OBB;  
        Eigen::Vector3f mass_center;

        double area;
        double wlRatio;
        double averageNumberOfPointsPerArea = 0;
        double floorAngle;
        // calculate features for each plane.
        for (size_t i = 0; i < planes.size(); ++i){
            
            feature_extractor.setInputCloud (planes[i]);
            feature_extractor.compute ();
            feature_extractor.getMassCenter( mass_center );
            feature_extractor.getOBB (min_point_OBB, max_point_OBB, position_OBB, rotational_matrix_OBB);
            Eigen::Vector3f position (position_OBB.x, position_OBB.y, position_OBB.z);
            Eigen::Quaternionf quat (rotational_matrix_OBB);
            planeFeatures::getBiggestCubeArea(min_point_OBB, max_point_OBB, &area, &wlRatio);

            if (viewerIsSet == true){                
                if (area > bigPlaneSize){
                    viewer->addCube (min_point_OBB.x, max_point_OBB.x, min_point_OBB.y, max_point_OBB.y, min_point_OBB.z, max_point_OBB.z, 1.0, 1.0, 0.0, "a"+std::to_string(i));
                    viewer->addCube (position, quat, max_point_OBB.x - min_point_OBB.x, max_point_OBB.y - min_point_OBB.y, max_point_OBB.z - min_point_OBB.z, "o"+std::to_string(i));
                }
            }
            averageNumberOfPointsPerArea += planes[i]->points.size() / area;


            // Find angle to floor and convert to range [0,1] where 0 means 90 degrees and 
            // 1 means zero or 180 degrees.
            floorAngle = planeFeatures::angleToFloor( normal.at( normalInd.at(i) ));
            floorAngle = std::abs(floorAngle - M_PI/2);
            
            // Check if it is a wall and set to inf;
            if ( area > bigPlaneSize &&  floorAngle < M_PI/5 ){
                fSet.features[i][0] = std::log(0);
                fSet.features[i][1] = std::log(0);
                fSet.features[i][2] = std::log(0);
                fSet.features[i][3] = std::log(0);
                fSet.walls.insert(i);
            }
            // Check if floor and set to inf
            else if ( mass_center(2) < 0.2 && floorAngle > M_PI/2 - M_PI/5 ){
                fSet.features[i][0] = std::log(0);
                fSet.features[i][1] = std::log(0);
                fSet.features[i][2] = std::log(0);
                fSet.features[i][3] = std::log(0);
                fSet.floors.insert(i);
            }
            else{
                fSet.features[i][0] = 2 * std::log( area );
                fSet.features[i][1] = std::log( wlRatio );
                fSet.features[i][2] = 2 * floorAngle / (M_PI/2);
                fSet.features[i][3] = 4 * mass_center(2);
            }
            fSet.position.push_back(position);
        }
    }

    void planeFeatures::matchFeatures(featureSet &fSet, std::vector<std::vector<int> > &c){
        dbscan scan;
        scan.cluster(fSet.features, fSet.walls, fSet.floors, 0.5, 2, c);
    }

    void planeFeatures::matchFeatures(featureSet &fSet){

        // Maximum neighbours found using radius search;
        int nn = 30;
        fSet.indices = flann::Matrix<int>(new int[fSet.features.rows*nn], fSet.features.rows, nn);

        flann::Matrix<double> dists(new double[fSet.features.rows*nn], fSet.features.rows, nn);

        // construct an randomized kd-tree index using 4 kd-trees
        flann::Index<flann::L2<double> > index(fSet.features, flann::KDTreeIndexParams(4));
        index.buildIndex();
        double radius = 1.0;
        std::vector<int> num_inliers;
        index.radiusSearch(fSet.features, fSet.indices, dists, radius, flann::SearchParams(128));
    }

    void planeFeatures::groupFeatures(featureSet &fSet){
        
        // load the matrix into sets
        std::set<int> set;
        for (size_t j = 0; j < fSet.indices.rows; ++j){
            for (size_t i = 0; i < fSet.indices.cols; ++i){
                if ( fSet.indices[j][i] == -1 ){ 
                    break; 
                }
                set.insert( fSet.indices[j][i] );
            }
            fSet.objects.insert(set);
            set.clear();
        }

        // Make sure sets don't intersect
        std::set<int > tmp;
        std::set<std::set<int> >::iterator it = fSet.objects.begin();
        std::set<std::set<int> >::iterator ite;
        for ( ; it != fSet.objects.end(); ++it){
            if (it->size() == 0){
                //sets.remove( it );
                continue;
            }
            for ( ite = it; ite != fSet.objects.end(); ++ite){
                if ( ite->size() == 0 || ite == it){
                    continue;
                }
                tmp.clear();
                if (it->size() >= ite->size()){
                    std::set_difference(ite->begin(), ite->end(), it->begin(), it->end(), 
                        std::inserter(tmp, tmp.begin()));
                    (*ite).swap(tmp);
                } else {
                    std::set_difference(it->begin(), it->end(), ite->begin(), ite->end(), 
                        std::inserter(tmp, tmp.begin()));
                    (*it).swap(tmp);
                    if(it->size() == 0){
                        break;
                    }
                }
            }
        }
    }

    void planeFeatures::improveWalls(const vPointCloudT &planes, featureSet &fSet){
        PointCloudT::Ptr cloud (new PointCloudT ());
        for (auto i : fSet.walls){
            *cloud += *planes.at(i);
        }
        pcl::KdTreeFLANN<PointT> kdtree;
        kdtree.setInputCloud (cloud);
        PointT searchPoint;

        std::set<std::set<int> > new_object;
        std::set<int> single_object;
        std::vector<int> pointIdx;
        std::vector<float> pointRadius;
        std::vector<int> newWalls;

        for ( auto set : fSet.objects ){
            single_object.clear();
            for ( auto i : set ){
                searchPoint.x = fSet.position.at(i)(0);
                searchPoint.y = fSet.position.at(i)(1);
                searchPoint.z = fSet.position.at(i)(2);
                if (kdtree.radiusSearch (searchPoint, centroid_to_wall_distance_, pointIdx, pointRadius) > 0){
                    fSet.walls.insert(i);
                    std::cout << "found probable wall." << std::endl;
                } else {
                    single_object.insert(i);
                }
            }
            new_object.insert(single_object);
            fSet.objectSize += single_object.size();
        }
        fSet.objects = new_object;
    }

    void planeFeatures::improveWalls(const vPointCloudT &planes, featureSet &fSet, std::vector<std::vector<int> > &c){
        PointCloudT::Ptr cloud (new PointCloudT ());
        for (auto i : fSet.walls){
            *cloud += *planes.at(i);
        }
        pcl::KdTreeFLANN<PointT> kdtree;
        kdtree.setInputCloud (cloud);
        PointT searchPoint;

        std::vector<std::vector<int> > new_object;
        std::vector<int> single_object;
        std::vector<int> pointIdx;
        std::vector<float> pointRadius;
        std::vector<int> newWalls;

        for ( auto vec : c ){
            single_object.clear();
            for ( auto i : vec ){
                searchPoint.x = fSet.position.at(i)(0);
                searchPoint.y = fSet.position.at(i)(1);
                searchPoint.z = fSet.position.at(i)(2);
                if (kdtree.radiusSearch (searchPoint, centroid_to_wall_distance_, pointIdx, pointRadius) > 0){
                    fSet.walls.insert(i);
                    std::cout << "found probable wall." << std::endl;
                } else {
                    single_object.push_back(i);
                }
            }
            new_object.push_back(single_object);
            fSet.objectSize += single_object.size();
        }
        c = new_object;
    }

    void planeFeatures::groupObjects(featureSet &fSet){


        // load center point of all planes into flann::Matrix
        int nn = 10;
        // flann::Matrix<int> indices (new int[fSet.objectSize*nn], fSet.objectSize, nn);
        // flann::Matrix<float> dists(new float[fSet.objectSize*nn], fSet.objectSize, nn);
        std::vector< std::vector<int> > indices;
        std::vector<std::vector<float> > dists;
        flann::Matrix<float> features(new float[fSet.objectSize*3], fSet.objectSize, 3);
        // construct an randomized kd-tree index using 4 kd-trees  
        std::map< int, std::pair<int, int> > maps;

        int k = 0;
        int j = 0;
        for ( auto set : fSet.objects){
            for ( auto i : set){    
                features[k][0] = fSet.position.at(i)(0);
                features[k][1] = fSet.position.at(i)(1);
                features[k][2] = fSet.position.at(i)(2);
                maps[k] = std::make_pair(j, i);
                k++;
            }
            j++;
        }
        flann::Index<flann::L2<float> > index(features, flann::KDTreeIndexParams(4));
        index.buildIndex();
        float radius = 1.0;
        index.radiusSearch(features, indices, dists, radius, flann::SearchParams(128));

        for (auto vec : indices){

        }
    }















    void planeFeatures::calculateFeatures(vPointCloudT planes, vPointCloudT hulls, std::vector<Eigen::Vector4d> normal, std::vector<int> normalInd, std::vector<planeDescriptor> *vPlaneDescriptor){
        // Check to see if they are the same size
        if ( planes.size() != hulls.size() ){
            // TODO: break in some way.
        }

        flann::Matrix<double> dataset(new double[planes.size()*3], planes.size(), 3 );

        pcl::MomentOfInertiaEstimation<PointT> feature_extractor;
        PointT min_point_OBB;
        PointT max_point_OBB;
        PointT position_OBB;
        Eigen::Matrix3f rotational_matrix_OBB;  
        // std::vector<float> moment_of_inertia;
        Eigen::Vector3f mass_center;

        planeDescriptor pDescriptor;
        double area;
        double wlRatio;
        bool hasAdded = false;
        double averageNumberOfPointsPerArea = 0;
        // calculate features for each plane.
        for (size_t i = 0; i < planes.size(); ++i){
            
            feature_extractor.setInputCloud (planes[i]);
            feature_extractor.compute ();
            //feature_extractor.getMomentOfInertia (moment_of_inertia);
            feature_extractor.getMassCenter( mass_center );
            feature_extractor.getOBB (min_point_OBB, max_point_OBB, position_OBB, rotational_matrix_OBB);
            Eigen::Vector3f position (position_OBB.x, position_OBB.y, position_OBB.z);
            Eigen::Quaternionf quat (rotational_matrix_OBB);
            
            pDescriptor.massCenter = mass_center;
            planeFeatures::getBiggestCubeArea(min_point_OBB, max_point_OBB, &area, &wlRatio);
            if (viewerIsSet == true){                
                if (area > bigPlaneSize){
                    viewer->addCube (min_point_OBB.x, max_point_OBB.x, min_point_OBB.y, max_point_OBB.y, min_point_OBB.z, max_point_OBB.z, 1.0, 1.0, 0.0, "a"+std::to_string(i));
                    viewer->addCube (position, quat, max_point_OBB.x - min_point_OBB.x, max_point_OBB.y - min_point_OBB.y, max_point_OBB.z - min_point_OBB.z, "o"+std::to_string(i));
                } else {
                    viewer->addCube (min_point_OBB.x, max_point_OBB.x, min_point_OBB.y, max_point_OBB.y, min_point_OBB.z, max_point_OBB.z, 1.0, 0.0, 1.0, "a"+std::to_string(i));
                }
                hasAdded = true;
            }
            pDescriptor.boundingBoxArea = area;
            pDescriptor.hullArea = double( pcl::calculatePolygonArea(*hulls[i]) );
            pDescriptor.hullBoundingBoxRatio = pDescriptor.hullArea / pDescriptor.boundingBoxArea;
            pDescriptor.widthLengthRatio = wlRatio;
            pDescriptor.normal = normal.at(normalInd.at(i));
            vPlaneDescriptor->push_back(pDescriptor);
            // printDescriptor(pDescriptor);
            averageNumberOfPointsPerArea += planes[i]->points.size() / area;
        }
        std::cout << "average number of points: " << averageNumberOfPointsPerArea / planes.size() << std::endl;
    }

    void planeFeatures::matchSimilarFeatures(std::vector<planeDescriptor> desc, std::vector<std::set<int> > *sets){
        
        double eDist;
        // eigen vDescriptor;
        // Eigen::Vector5d vDescriptor;
        Eigen::VectorXd vDescriptor(5);
        Eigen::MatrixXd mDescriptor( desc.size(), desc.size() );
        int red, green, blue;
        double angle, norm;

        std::set<int> set;
        std::set<int> walls;
        std::set<int> floors;
        std::set<int> filled;

        double floorAngle;

        for ( size_t i = 0; i < desc.size(); ++i){
            
            floorAngle = planeFeatures::angleToFloor( desc[i].normal );
            // find walls.
            if ( desc[i].boundingBoxArea > bigPlaneSize &&  std::abs( floorAngle - M_PI/2 ) < M_PI/5  ){
                walls.insert( i );
                filled.insert( i );
                std::cout << "Found big wall: " << desc[i].boundingBoxArea << std::endl;
                continue;
            }
            else if ( desc[i].massCenter(2) < 0.2 && std::abs( floorAngle - M_PI/2 ) > M_PI/2 - M_PI/5 ){
                floors.insert( i );
                filled.insert( i );
                std::cout << "Found floor: " << std::endl;
                continue;
            }

            for (size_t j = 0; j < desc.size(); ++j){
                eDist = 0;
                // Check distance 
                vDescriptor(0) = ( std::log( desc[i].boundingBoxArea ) - std::log( desc[j].boundingBoxArea ));
                vDescriptor(1) = 0; //( std::log( desc[i].hullArea ) - std::log( desc[j].hullArea ));
                vDescriptor(2) = 0; //( std::log( desc[i].hullBoundingBoxRatio ) - std::log( desc[j].hullBoundingBoxRatio ));
                vDescriptor(3) = 2 * ( std::log( desc[i].widthLengthRatio ) - std::log( desc[j].widthLengthRatio ));
                vDescriptor(4) = 2 * std::log( planeFeatures::angleBetweenVectors( desc[i].normal, desc[j].normal ));
                norm = vDescriptor.norm();
                mDescriptor(i,j) = norm;
                // planeFeatures::getValueBetweenTwoFixedColors(norm, &red, &green, &blue);
            }
        }

        // std::cout << "similarity matrix is:" << std::endl;
        // std::cout << mDescriptor << std::endl;

        // Loop through all columns and group together similar planes.
        // TODO: Improve.
        for (size_t i = 0; i < mDescriptor.cols(); ++i){
            if ( filled.count(i) ){
                continue;
            }
            set.clear();
            // std::cout << "set for col: " << i << std::endl;
            for (size_t j = i+1; j < mDescriptor.rows(); ++j){
                if (mDescriptor(i,j) < 1.0 && !filled.count(j) ){
                    set.insert( j );
                    // std::cout << j << "  ";
                    filled.insert(j);
                }
            }
            if ( !set.empty() ){
                set.insert( i );
                sets->push_back( set );
            }
            // std::cout << "" << std::endl;
        }
        sets->push_back( floors );
        sets->push_back( walls );
    }   


    void planeFeatures::findRepeatingObjects(vPointCloudT planes, std::vector<planeDescriptor> desc, std::vector<std::set<int> > sets, std::vector< std::set<int>> *objects){

        double size = sets.size()-2;
        Eigen::MatrixXd matching = Eigen::MatrixXd::Zero(size, size);
        double hmm;
        for (int iset = 0; iset < size; ++iset){
            for ( auto i : sets.at(iset) ){
                if ( desc.at(i).boundingBoxArea > bigPlaneSize ){ continue; }
                
                for (int jset = 0; jset < size; ++jset){
                    for ( auto j : sets.at(jset) ){
                        if ( j == i ) { continue; }
                        if ( desc.at(j).boundingBoxArea > bigPlaneSize ){ continue; }
                        hmm = planeFeatures::euclideanDistance( desc.at(i).massCenter, desc.at(j).massCenter);
                        std::cout << hmm << std::endl;
                        if ( hmm > 1.0 ){ continue; }
                        
                        if ( iset > jset){
                            matching(iset,jset) += 1; 
                        } else {
                            matching(jset,iset) += 1; 
                        }
                    }
                }
            }
        }

        std::cout << "Object matching" << std::endl;
        std::cout << matching << std::endl;
        std::set<int> s;
        Eigen::MatrixXd::Index maxRow, maxCol;
        // for ( size_t i = 0; i < matching.rows(); ++i ){
        //     for ( size_t j = 0; j < matching.cols(); ++j  ){
        //         if (matching(i,j) > 14){
        //             s.clear();
        //             s.insert(i);
        //             s.insert(j);
        //             objects->push_back( s );
        //         }
        //     }
        // }
        matching.maxCoeff(&maxRow, &maxCol);
        s.clear();
        s.insert(maxRow);
        s.insert(maxCol);
        objects->push_back( s );

    }

    void planeFeatures::setViewer(boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer){
        this->viewer = viewer;
        viewerIsSet = true;
    }

	// Takes in two vectors, returns angle between them in range 0 to 1  
    // where 1 means no difference, pi/2 is cnsidered maximum angle and pi as 0.
    double planeFeatures::angleBetweenVectors(Eigen::Vector4d a, Eigen::Vector4d b){
        if (a == b) { return 1; }
        
        Eigen::Vector3d an(a(0),a(1),a(2));
        Eigen::Vector3d bn(b(0),b(1),b(2));
        Eigen::Vector3d cn(0,0,1);

        // Norm the vectors before calculation.
        an = an / an.norm();
        bn = bn / bn.norm();
        cn = cn / cn.norm();

        double ang1 = std::abs( std::acos( an.dot(cn)) - M_PI/2 ) / (M_PI/2); 
        double ang2 = std::abs( std::acos( bn.dot(cn)) - M_PI/2 ) / (M_PI/2); 

        double ret = std::abs( ang1 - ang2);

        // std::cout << "angle1: " << ang1*90 << "  ";
        // std::cout << "angle2: " << ang2*90 << "  ";
        // std::cout << "return: " << 2* std::log(1 - ret) << "  ";
        // std::cout << "an: " << an.transpose() << "  ";
        // std::cout << "bn: " << bn.transpose() << std::endl;

        return 1 - ret;
    }

    double planeFeatures::angleToFloor( Eigen::Vector4d a ){
        
        Eigen::Vector3d an(a(0),a(1),a(2));
        Eigen::Vector3d fn(0,0,1);
        an = an / an.norm();
        fn = fn / fn.norm();

        return std::abs( std::acos( an.dot(fn) ) );
    }

    // Returns RGB value from red to green depending on value, low value results in blue,
    // high value results in red.
    void planeFeatures::getValueBetweenTwoFixedColors(double value, int *red, int *green, int *blue)
    {
        if (value > 10.0 && value == INFINITY){ 
            value = 1.0; 
        } else {
            value = value / 10.0;
        }

        int aR = 0;   int aG = 0; int aB=255;  // RGB for our 1st color (blue in this case).
        int bR = 255; int bG = 0; int bB=0;    // RGB for our 2nd color (red in this case).

        *red   = double(bR - aR) * value + aR;      // Evaluated as -255*value + 255.
        *green = double(bG - aG) * value + aG;      // Evaluates as 0.
        *blue  = double(bB - aB) * value + aB;      // Evaluates as 255*value + 0.
    }

    void planeFeatures::getBiggestCubeArea(PointT minPoint, PointT maxPoint, double *area, double *WLRatio){
        double x = std::abs( maxPoint.x - minPoint.x );
        double y = std::abs( maxPoint.y - minPoint.y );
        double z = std::abs( maxPoint.z - minPoint.z );
        // if (x > z && y > z){
            *area = x*y;
            *WLRatio = x/y;
        // } else if (x > y && z > y) {
        //     *area = x*z;
        //     *WLRatio = x/z;
        // }
        // else {
        //     *area = y*z;
        //     *WLRatio = y/z;
        // }
        if( *WLRatio > 1 ){
            *WLRatio = 1 / *WLRatio;
        }
    }

    double planeFeatures::euclideanDistance(Eigen::Vector3f a, Eigen::Vector3f b){
        a(0) -= b(0);
        a(1) -= b(1);
        a(2) -= b(2);

        return a.norm();
    }


    void planeFeatures::printDescriptor(planeDescriptor descr){
        std::cout << "Descriptor" << std::endl;
        std::cout << " " << std::endl;
        std::cout << "Bounding Box Area:\t" << descr.boundingBoxArea << std::endl;
        std::cout << "Hull Area:        \t" << descr.hullArea << std::endl;
        std::cout << "Hull Bounding box Ratio:\t" << descr.hullBoundingBoxRatio << std::endl;
        std::cout << "Width Length Ratio:\t" << descr.widthLengthRatio << std::endl;
        std::cout << "Normal:           \t" << descr.normal.transpose() << std::endl;
        std::cout << " " << std::endl;

    }

    float planeFeatures::calculatePolygonAreaMine (const PointCloudT &polygon) 
    {
        float area = 0.0f;
        int num_points = polygon.size ();
        int j = 0;
        Eigen::Vector3f va,vb,res;

        res(0) = res(1) = res(2) = 0.0f;
        for (int i = 0; i < num_points; ++i) 
        {
            j = (i + 1) % num_points;
            va = polygon[i].getVector3fMap ();
            vb = polygon[j].getVector3fMap ();
            res += va.cross (vb);
        }
        area = res.norm ();
        return (area*0.5);
    }

    float planeFeatures::calculatePolygonAreaMineInv (const PointCloudT &polygon) 
    {
        float area = 0.0f;
        int num_points = polygon.size ();
        int j = 0;
        Eigen::Vector3f va,vb,res;

        res(0) = res(1) = res(2) = 0.0f;
        for (int i = num_points-1; i >= 0; --i) 
        {
            j = (i + 1) % num_points;
            va = polygon[i].getVector3fMap ();
            vb = polygon[j].getVector3fMap ();
            res += va.cross (vb);
        }
        area = res.norm ();
        return (area*0.5);
    }

}
#include <plane_features/plane_features.h>
#include <pcl/common/common.h>


namespace EXX{


    void planeFeatures::calculateFeatures(vPointCloudT planes, vPointCloudT hulls, std::vector<Eigen::Vector4d> normal, std::vector<int> normalInd, std::vector<planeDescriptor> *vPlaneDescriptor){
        // Check to see if they are the same size
        if ( planes.size() != hulls.size() ){
            // TODO: break in some way.
        }

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
                if (area > 0.8){
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
            pDescriptor.normal = normal.at(     normalInd.at(i) );
            vPlaneDescriptor->push_back(pDescriptor);
            // printDescriptor(pDescriptor);
        }
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
            if ( desc[i].boundingBoxArea > 0.8 &&  std::abs( floorAngle - M_PI/2 ) < M_PI/5  ){
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
                vDescriptor(3) = ( std::log( desc[i].widthLengthRatio ) - std::log( desc[j].widthLengthRatio ));
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

}
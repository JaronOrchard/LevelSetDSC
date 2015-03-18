//
//  Point cloud fitting velocity function for the Deformable Simplicial Complex
//  Created by Jared Saul (2015)

#pragma once

#include "velocity_function.h"

/**
 A velocity function which moves the interface vertices towards a point cloud.
 */
class PointCloudFunc: public DSC::VelocityFunc<> {
    
private:

#ifdef _WIN32
    const std::string obj_path = "data\\";
    const std::string log_path = "LOG\\";
#else
    const std::string obj_path = "./data/";
    const std::string log_path = "./LOG/";
#endif

    std::vector<vec3> point_cloud;
    

    // Variables:
    std::string point_cloud_file_name = "new\\teapot.obj";
    real scale_target = 1.1; // Outermost edge of imported .obj should reach this
    real alpha = 0.2;


    void import_point_cloud() {
        point_cloud.clear();
        std::cout << "Loading points from " << point_cloud_file_name << "...";
        std::string filename = obj_path + point_cloud_file_name;
        std::ifstream fin(filename.data());
        if (fin.fail()) {
            std::cout << "FAILED." << std::endl;
            std::cout << " - Error opening " << point_cloud_file_name << std::endl;
            std::cout << " - Using default sample points." << std::endl;
            point_cloud.push_back(vec3(0, 1.5, 0));
            point_cloud.push_back(vec3(0, -1.5, 0));
            point_cloud.push_back(vec3(1.5, 0, 0));
            point_cloud.push_back(vec3(0, 0, 1.5));
        } else {
            std::string temp;
            fin >> temp;
            while (!fin.eof()) {
                if (temp == "v") { // Vertex
                    real x, y, z; // The (x,y,z) coordinates of a vertex.
                    fin >> x >> y >> z;
                    point_cloud.push_back(vec3(x, y, z));
                }
                fin >> temp;
            }
            fin.close();
            std::cout << "done." << std::endl;
            std::cout << " - Read in " << point_cloud.size() << " points." << std::endl;
        }   
    }

    void normalize_point_cloud() {
        std::cout << "Normalizing point cloud..." << std::endl;
        if (point_cloud.empty()) {
            std::cout << "FAILED." << std::endl;
            std::cout << " - Error: Point cloud is empty." << std::endl;
            return;
        }
        real min_x = INFINITY;
        real min_y = INFINITY;
        real min_z = INFINITY;
        real max_x = -INFINITY;
        real max_y = -INFINITY;
        real max_z = -INFINITY;
        for (size_t i = 0; i < point_cloud.size(); i++) {
            min_x = std::min(min_x, point_cloud[i][0]);
            min_y = std::min(min_y, point_cloud[i][1]);
            min_z = std::min(min_z, point_cloud[i][2]);
            max_x = std::max(max_x, point_cloud[i][0]);
            max_y = std::max(max_y, point_cloud[i][1]);
            max_z = std::max(max_z, point_cloud[i][2]);
        }
        
        std::cout << "Centering point cloud on the origin..." << std::endl;
        real x_shift = -((min_x + max_x) / 2);
        real y_shift = -((min_y + max_y) / 2);
        real z_shift = -((min_z + max_z) / 2);
        if (x_shift != 0) {
            std::cout << " - Shifting x values by " << x_shift << "." << std::endl;
            for (size_t i = 0; i < point_cloud.size(); i++) {
                point_cloud[i][0] += x_shift;
            }
        }
        if (y_shift != 0) {
            std::cout << " - Shifting y values by " << y_shift << "." << std::endl;
            for (size_t i = 0; i < point_cloud.size(); i++) {
                point_cloud[i][1] += y_shift;
            }
        }
        if (z_shift != 0) {
            std::cout << " - Shifting z values by " << z_shift << "." << std::endl;
            for (size_t i = 0; i < point_cloud.size(); i++) {
                point_cloud[i][2] += z_shift;
            }
        }
        std::cout << "Scaling point_cloud to a max of " << scale_target << "..." << std::endl;
        real x_range = (max_x - min_x) / 2;
        real y_range = (max_y - min_y) / 2;
        real z_range = (max_z - min_z) / 2;
        real max_range = std::max(x_range, std::max(y_range, z_range));
        real scale_factor = scale_target / max_range;
        std::cout << " - Scaling points by " << scale_factor << "..." << std::endl;
        for (size_t i = 0; i < point_cloud.size(); i++) {
            for (size_t j = 0; j < 3; j++) {
                point_cloud[i][j] *= scale_factor;
            }
        }
    }

    real get_angular_defect(DSC::DeformableSimplicialComplex<>& dsc, is_mesh::NodeKey nodeKey) {
        real pi = 3.141592654;
        int numFaces = 0;
        real totalAngles = 0;
        for (auto f : dsc.get_faces(nodeKey)) {
            if (dsc.get(f).is_interface()) {
                real angle;
                numFaces++;
                auto faceNodes = dsc.get_nodes(f);
                vec3 nodePosOrig = dsc.get_pos(nodeKey);
                vec3 nodePos1 = dsc.get(faceNodes[0]).get_pos();
                vec3 nodePos2 = dsc.get(faceNodes[1]).get_pos();
                vec3 nodePos3 = dsc.get(faceNodes[2]).get_pos();
                vec3 len12 = nodePos2 - nodePos1;
                vec3 len13 = nodePos3 - nodePos1;
                vec3 len23 = nodePos3 - nodePos2;
                if (nodePosOrig == nodePos1) { // Node 1 is Node nodeKey
                    if (length(len12) >= length(len13) && length(len12) >= length(len23)) {
                        angle = acos(length(len13) / length(len12)); // 1-2 is hypotenuse
                    } else if (length(len13) >= length(len23)) {
                        angle = acos(length(len12) / length(len13)); // 1-3 is hypotenuse
                    } else {
                        angle = pi - acos(length(len12) / length(len23)) - acos(length(len13) / length(len23)); // 2-3 is hypotenuse
                    }
                }
                else if (nodePosOrig == nodePos2) { // Node 2 is Node nodeKey
                    if (length(len12) >= length(len13) && length(len12) >= length(len23)) {
                        angle = acos(length(len23) / length(len12)); // 1-2 is hypotenuse
                    } else if (length(len13) >= length(len23)) {
                        angle = pi - acos(length(len12) / length(len13)) - acos(length(len23) / length(len13)); // 1-3 is hypotenuse
                    } else {
                        angle = acos(length(len12) / length(len23)); // 2-3 is hypotenuse
                    }
                }
                else { // Node 3 is Node nodeKey
                    if (length(len12) >= length(len13) && length(len12) >= length(len23)) {
                        angle = pi - acos(length(len13) / length(len12)) - acos(length(len23) / length(len12)); // 1-2 is hypotenuse
                    } else if (length(len13) >= length(len23)) {
                        angle = acos(length(len23) / length(len13)); // 1-3 is hypotenuse
                    } else {
                        angle = acos(length(len13) / length(len23)); // 2-3 is hypotenuse
                    }
                }
                totalAngles += angle;
            }
        }
        return (pi * 2 - totalAngles);
    }

public:
    /**
     Creates a velocity function which moves the interface vertices towards a point cloud.
     */
    PointCloudFunc(real velocity, real accuracy, int max_time_steps = 500) :
        VelocityFunc<>(velocity, accuracy, max_time_steps)
    {
        import_point_cloud();
        normalize_point_cloud();
        std::cout << std::endl;
    }

    /**
     Returns the name of the velocity function.
     */
    virtual std::string get_name() const
    {
        return std::string("POINT CLOUD FITTING");
    }
    
    /**
     Computes the motion of each interface vertex and stores the new position in new_pos in the simplicial complex class.
     */
    virtual void deform(DSC::DeformableSimplicialComplex<>& dsc)
    {
        // Formula:
        //   Speed(x):     alpha * (Normal [dot] (p - x))
        //                    where p is the closest point in the point cloud to x
        //   Movement(x):  speed(x) * Normal

        auto init_time = std::chrono::system_clock::now();
        vec3 new_pos;
        vec3 p_minus_x;
        real closest_dist;
        vec3 closest_point;
        vec3 point_normal;
        real dot_product;
        for(auto nit = dsc.nodes_begin(); nit != dsc.nodes_end(); nit++)
        {
            if(dsc.is_movable(nit.key()))
            {
                // Find the closest point in the point cloud:
                closest_point = point_cloud[0];
                closest_dist = (closest_point - nit->get_pos()).length();
                for (size_t i = 1; i < point_cloud.size(); i++) {
                    if ((point_cloud[i] - nit->get_pos()).length() < closest_dist) {
                        closest_point = point_cloud[i];
                        closest_dist = (closest_point - nit->get_pos()).length();
                    }
                }
                // Calculate the movement vector and set destination:
                point_normal = dsc.get_normal(nit.key());
                p_minus_x = closest_point - nit->get_pos();
                dot_product = point_normal[0] * p_minus_x[0] + point_normal[1] * p_minus_x[1] + point_normal[2] * p_minus_x[2];
                //std::cout << "POINT: " << nit->get_pos() << std::endl;
                //std::cout << "SPEED: " << alpha*dot_product << std::endl;
                //std::cout << "MOVEMENT: " << (alpha*dot_product)*point_normal << std::endl;
                real speed = alpha * dot_product;

                real angular_defect_constant = 0.0025;
                real angular_defect = get_angular_defect(dsc, nit.key());
                if (abs(angular_defect) > 0.001) {
                    speed -= (angular_defect_constant * angular_defect);
                }

                new_pos = (speed * point_normal) + nit->get_pos();
                dsc.set_destination(nit.key(), new_pos);
            }
        }
        update_compute_time(init_time);
        VelocityFunc::deform(dsc);
    }
};

//
//  Point cloud fitting velocity function for the Deformable Simplicial Complex
//  Created by Jared Saul (2015)

#pragma once

#include "velocity_function.h"

/**
 A velocity function which moves the interface vertices towards a point cloud.
 */
class PointCloudFunc: public DSC::VelocityFunc<> {
    
    
public:
    /**
     Creates a velocity function which moves the interface vertices towards a point cloud.
     */
    PointCloudFunc(real velocity, real accuracy, int max_time_steps = 500) :
        VelocityFunc<>(velocity, accuracy, max_time_steps)
    {
        
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

        real alpha = 0.1;

        std::vector<vec3> point_cloud;
        vec3 high_point(0, 1.5, 0);
        vec3 low_point(0, -1.5, 0);
        vec3 x_point(1.5, 0, 0);
        vec3 z_point(0, 0, 1.5);
        point_cloud.push_back(high_point);
        point_cloud.push_back(low_point);
        point_cloud.push_back(x_point);
        point_cloud.push_back(z_point);

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
                new_pos = ((alpha * dot_product) * point_normal) + nit->get_pos();
                dsc.set_destination(nit.key(), new_pos);
            }
        }
        update_compute_time(init_time);
        VelocityFunc::deform(dsc);
    }
};

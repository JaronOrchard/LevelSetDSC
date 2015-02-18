//
//  Level set method for the Deformable Simplicial Complex
//  Created by Jared Saul (2015)

#pragma once

#include "velocity_function.h"

/**
 A velocity function which moves the interface vertices towards a point cloud.
 */
class LevelSetFunc: public DSC::VelocityFunc<> {
    
    
public:
    /**
     Creates a velocity function which moves the interface vertices towards a point cloud.
     */
    LevelSetFunc(real velocity, real accuracy, int max_time_steps = 500) :
        VelocityFunc<>(velocity, accuracy, max_time_steps)
    {
        
    }

    /**
     Returns the name of the velocity function.
     */
    virtual std::string get_name() const
    {
        return std::string("LEVEL SET MOTION");
    }
    
    /**
     Computes the motion of each interface vertex and stores the new position in new_pos in the simplicial complex class.
     */
    virtual void deform(DSC::DeformableSimplicialComplex<>& dsc)
    {
        auto init_time = std::chrono::system_clock::now();
        vec3 new_pos;
        for(auto nit = dsc.nodes_begin(); nit != dsc.nodes_end(); nit++)
        {
            if(dsc.is_movable(nit.key()))
            {
                new_pos = nit->get_pos() + 0.1*VELOCITY * dsc.get_normal(nit.key());
                dsc.set_destination(nit.key(), new_pos);
            }
        }
        update_compute_time(init_time);
        VelocityFunc::deform(dsc);
    }
};

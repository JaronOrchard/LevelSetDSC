## PointCloudDSC
### An implementation of the Deformable Simplicial Complex (DSC) method modified to reconstruct surfaces from point clouds
#### Modified by Jared Saul under the direction of Professor John C. Hart at the University of Illinois at Urbana-Champaign

---
### How to use
1. Set the `model_file_name` in `user_interface.h` to "cube" and the `point_cloud_file_name` in `point_cloud_function.h` to point at whatever .obj file you wish to reconstruct
2. Tweak additional settings in `point_cloud_function.h` as desired
3. Start the program
4. Hit `4` to load the point cloud and set it as the new velocity function
5. Hit `m` to move 1 time step or `M` to move multiple time steps (typically 50 or 100, defined in `user_interface.cpp`)
6. `i` exports the mesh to the data directory, while `I` exports speed and curvature data for each vertex position.  Other added commands are listed at the bottom of the terminal window

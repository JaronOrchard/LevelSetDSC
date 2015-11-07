//
//  Deformabel Simplicial Complex (DSC) method
//  Copyright (C) 2013  Technical University of Denmark
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  See licence.txt for a copy of the GNU General Public License.

#include "user_interface.h"
#include "rotate_function.h"
#include "average_function.h"
#include "normal_function.h"
#include "point_cloud_function.h"

#include <iostream>
#include <iomanip>
#include <ctime>
#include <chrono>

using namespace DSC;
using std::cout;
using std::endl;
using std::vector;

void display_(){
    UI::get_instance()->display();
}

void keyboard_(unsigned char key, int x, int y){
    UI::get_instance()->keyboard(key, x, y);
}

void reshape_(int width, int height){
    UI::get_instance()->reshape(width, height);
}

void visible_(int v){
    UI::get_instance()->visible(v);
}

void animate_(){
    UI::get_instance()->animate();
}

UI* UI::instance = NULL;

UI::UI(int &argc, char** argv)
{
    instance = this;

    glutInit(&argc, argv);
#ifdef _WIN32
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH | GLUT_MULTISAMPLE);
#else
    glutInitDisplayMode(GLUT_3_2_CORE_PROFILE | GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH | GLUT_MULTISAMPLE);
#endif
    glutCreateWindow("");
    
    glutDisplayFunc(display_);
    glutKeyboardFunc(keyboard_);
	glutIgnoreKeyRepeat(true);
    glutVisibilityFunc(visible_);
    glutReshapeFunc(reshape_);
	glutIdleFunc(animate_);

#ifndef __APPLE__
	glewExperimental = GL_TRUE;  // See http://www.opengl.org/wiki/OpenGL_Loading_Library
	GLint GlewInitResult = glewInit();
	if (GlewInitResult != GLEW_OK) {
		printf("ERROR: %s\n", glewGetErrorString(GlewInitResult));
	}
    check_gl_error(); // Catches a GL_INVALID_ENUM error. See http://www.opengl.org/wiki/OpenGL_Loading_Library
#endif
    
    // Read input
    std::string motion = "";
    real discretization = 2.5;
    real velocity = 5.;
    real accuracy = 0.25;
    
    if(argc == 2)
    {
        model_file_name = std::string(argv[1]);
    }
    else if(argc > 2)
    {
        for(int i = 0; i < argc; ++i)
        {
            std::string str(argv[i]);
            if (str == "nu") {
                velocity = std::atof(argv[i+1]);
            }
            else if (str == "delta") {
                discretization = std::atof(argv[i+1]);
            }
            else if (str == "alpha") {
                accuracy = std::atof(argv[i+1]);
            }
            else if (str == "model") {
                model_file_name = argv[i+1];
            }
            else if (str == "motion") {
                motion = argv[i+1];
            }
        }
    }
    painter = std::unique_ptr<Painter>(new Painter(light_pos));
    load_model(model_file_name, discretization);
    
    if(motion.empty())
    {
        vel_fun = std::unique_ptr<VelocityFunc<>>(new VelocityFunc<>(velocity, accuracy, 500));
        start("");
    }
    else {
        keyboard(*motion.data(), 0, 0);
    }
    
	glutReshapeWindow(WIN_SIZE_X, WIN_SIZE_Y);
    check_gl_error();
}

void UI::load_model(const std::string& file_name, real discretization)
{
    std::cout << "\nLoading " << obj_path + file_name + ".dsc" << std::endl;
    dsc = nullptr;
    std::vector<vec3> points;
    std::vector<int>  tets;
    std::vector<int>  tet_labels;
    is_mesh::import_tet_mesh(obj_path + file_name + ".dsc", points, tets, tet_labels);
    
    dsc = std::unique_ptr<DeformableSimplicialComplex<>>(new DeformableSimplicialComplex<>(points, tets, tet_labels));
    
    vec3 p_min(INFINITY), p_max(-INFINITY);
    for (auto nit = dsc->nodes_begin(); nit != dsc->nodes_end(); nit++) {
        for (int i = 0; i < 3; i++) {
            p_min[i] = Util::min(nit->get_pos()[i], p_min[i]);
            p_max[i] = Util::max(nit->get_pos()[i], p_max[i]);
        }
    }
    
    vec3 size = p_max - p_min;
    real var = Util::max(Util::max(size[0], size[1]), size[2]);
    real dist = 1.2*var;
    eye_pos = {dist, var, dist};
    camera_pos = {var, var, -dist};
    per_step_screenshot_pos = {var*0.1, var*0.2, -dist*0.6};
    light_pos = {0., 0., dist};
    
    painter->update(*dsc);
    std::cout << "Loading done" << std::endl << std::endl;

    /*
    Degenerate edge quality : 0.1
    Minimum edge quality : 0.5
    Degenerate face quality : 0.0005
    Minimum face quality : 0.015
    Degenerate tet quality : 0.02
    Minimum tet quality : 0.3
    Minimum edge length : 0
    Maximum edge length : 2
    Minimum face area : 0.2
    Maximum face area : 5
    Minimum tet volume : 0.2
    Maximum tet volume : Infinity
    */
    std::cout << "Setting new parameters..." << std::endl << std::endl;
    //parameters old_pars = { 0.1, 0.5, 0.0005, 0.015, 0.02, 0.3, 0., 2., 0.2, 5., 0.2, INFINITY };
    parameters new_pars = { 0.1, 0.5, 0.0005, 0.015, 0.03, 0.2, 0., 2., 0.2, 5., 0.2, INFINITY };
    //parameters april_pars = { 0.1, 0.5, 0.0005, 0.015, 0.03, 0.2, 0., 2., 0.1, 3., 0.1, INFINITY };
    dsc->set_parameters(new_pars);
    
}

void UI::update_title()
{
    std::ostringstream oss;
    oss << "3D DSC\t" << vel_fun->get_name() << ", Time step " << vel_fun->get_time_step();
    oss << " (Nu = " << vel_fun->get_velocity() << ", Delta = " << dsc->get_avg_edge_length() << ", Alpha = " << vel_fun->get_accuracy() << ")";
    std::string str(oss.str());
    glutSetWindowTitle(str.c_str());
}

void UI::display()
{
    if (glutGet(GLUT_WINDOW_WIDTH) != WIN_SIZE_X || glutGet(GLUT_WINDOW_HEIGHT) != WIN_SIZE_Y) {
        return;
    }
    GLfloat timeValue = glutGet(GLUT_ELAPSED_TIME)*0.0002;
    //vec3 ep( eye_pos[0] * sinf(timeValue), eye_pos[1] * cosf(timeValue) , eye_pos[2] * cosf(timeValue));
    vec3 ep(eye_pos[0] * sinf(timeValue), eye_pos[1] * cosf(timeValue) * 0.1, eye_pos[2] * cosf(timeValue));
    painter->set_view_position(ep);
    painter->draw();
    glutSwapBuffers();
    update_title();
    check_gl_error();
}

void UI::reshape(int width, int height)
{
    WIN_SIZE_X = width;
    WIN_SIZE_Y = height;
    painter->reshape(width, height);
}

void UI::animate()
{
    if(CONTINUOUS)
    {
        std::cout << "\n***************TIME STEP " << vel_fun->get_time_step() + 1 <<  " START*************\n" << std::endl;
        vel_fun->take_time_step(*dsc);
        painter->update(*dsc);
        if(RECORD && basic_log)
        {
            painter->set_view_position(camera_pos);
            painter->save_painting(basic_log->get_path(), vel_fun->get_time_step());
            basic_log->write_timestep(*vel_fun, *dsc);
        }
        if (vel_fun->is_motion_finished(*dsc))
        {
            stop();
            if (QUIT_ON_COMPLETION) {
                exit(0);
            }
        }
        std::cout << "\n***************TIME STEP " << vel_fun->get_time_step() <<  " STOP*************\n" << std::endl;
    }
    glutPostRedisplay();
}

void UI::keyboard(unsigned char key, int x, int y) {
    switch(key) {
        case '\033':
            stop();
            exit(0);
            break;
        case '0':
            stop();
            QUIT_ON_COMPLETION = false;
            RECORD = false;
            vel_fun = std::unique_ptr<VelocityFunc<>>(new VelocityFunc<>(vel_fun->get_velocity(), vel_fun->get_accuracy(), 500));
            start("");
            break;
        case '1':
            stop();
            QUIT_ON_COMPLETION = true;
            RECORD = true;
            vel_fun = std::unique_ptr<VelocityFunc<>>(new RotateFunc(vel_fun->get_velocity(), vel_fun->get_accuracy()));
            start("rotate");
            break;
        case '2':
            stop();
            QUIT_ON_COMPLETION = true;
            RECORD = true;
            vel_fun = std::unique_ptr<VelocityFunc<>>(new AverageFunc(vel_fun->get_velocity(), vel_fun->get_accuracy()));
            start("smooth");
            break;
        case '3':
            stop();
            QUIT_ON_COMPLETION = true;
            RECORD = true;
            vel_fun = std::unique_ptr<VelocityFunc<>>(new NormalFunc(vel_fun->get_velocity(), vel_fun->get_accuracy()));
            start("expand");
            break;
        case '4':
            stop();
            QUIT_ON_COMPLETION = true;
            RECORD = true;
            vel_fun = std::unique_ptr<VelocityFunc<>>(new PointCloudFunc(vel_fun->get_velocity(), vel_fun->get_accuracy()));
            start("point_cloud");
            break;
        case ' ':
            if(!CONTINUOUS)
            {
                std::cout << "MOTION STARTED" << std::endl;
                if(RECORD && basic_log)
                {
                    painter->set_view_position(camera_pos);
                    painter->save_painting(basic_log->get_path(), vel_fun->get_time_step());
                }
            }
            else {
                std::cout << "MOTION PAUSED" << std::endl;
            }
            CONTINUOUS = !CONTINUOUS;
            break;
        case 'c': {
            int node_count = 0;
            int face_count = 0;
            for (auto nit = dsc->nodes_begin(); nit != dsc->nodes_end(); nit++)
            {
                if (dsc->is_movable(nit.key()))
                {
                    node_count++;
                }
            }
            for (auto fit = dsc->faces_begin(); fit != dsc->faces_end(); fit++) {
                if (dsc->get(fit.key()).is_interface()) {
                    face_count++;
                }
            }
            std::cout << "INTERFACE NODE COUNT: " << node_count << std::endl;
            std::cout << "INTERFACE FACE COUNT: " << face_count << std::endl;
            std::cout << "AVERAGE EDGE LENGTH: " << dsc->get_avg_edge_length() << std::endl;
            std::cout << std::endl;
            vel_fun->analyze_result(*dsc);
            std::cout << std::endl;
            break;
        }
        case 'C': {
            vel_fun->print_quality_measure(*dsc);
            break;
        }
        case 'n': {
            // Calculate the curvature for all nodes using Implicit Fairing's equation 14.
            int i = 0;
            for (auto nit = dsc->nodes_begin(); nit != dsc->nodes_end(); nit++)
            {
                if (dsc->is_movable(nit.key()))
                {
                    // To perform Implicit Fairing's equation 14, do this for each vertex:
                    // 1) For each *neighboring* vertex/edge...
                    // 2) Find angles alpha and beta
                    // 3) Multiply sum of cotangents by vector x_j - x_i
                    // 4) Sum all triangles' areas (or sum as you go and divide by 2)
                    // 5) Divide by 4A and take length of normal as curvature
                    vec3 sum(0, 0, 0);
                    real totalArea = 0;
                    for (auto e : dsc->get_edges(nit.key())) {
                        if (dsc->get(e).is_interface()) { // For each interface edge around the vertex...
                            vector<vec3> sixFaceNodes;
                            for (auto f : dsc->get(e).face_keys()) {
                                if (dsc->get(f).is_interface()) { // For both interface faces around that edge...
                                    auto faceNodes = dsc->get_nodes(f);
                                    vec3 nodePos1 = dsc->get(faceNodes[0]).get_pos();
                                    vec3 nodePos2 = dsc->get(faceNodes[1]).get_pos();
                                    vec3 nodePos3 = dsc->get(faceNodes[2]).get_pos();
                                    sixFaceNodes.push_back(nodePos1);
                                    sixFaceNodes.push_back(nodePos2);
                                    sixFaceNodes.push_back(nodePos3);
                                }
                            }
                            // The sixFaceNodes vector should now have 6 node positions in it from the two interface faces around edge e.
                            // Of [0,1,2] and [3,4,5], two are shared between the sets, and one is unique in both set.
                            // We can find the alpha and beta nodes by finding the odd node out in both sets of three.
                            int alphaIdx, betaIdx, otherIdx1, otherIdx2;
                            if (sixFaceNodes[0] != sixFaceNodes[3] && sixFaceNodes[0] != sixFaceNodes[4] && sixFaceNodes[0] != sixFaceNodes[5]) { alphaIdx = 0; otherIdx1 = 1; otherIdx2 = 2; }
                            else if (sixFaceNodes[1] != sixFaceNodes[3] && sixFaceNodes[1] != sixFaceNodes[4] && sixFaceNodes[1] != sixFaceNodes[5]) { alphaIdx = 1; otherIdx1 = 0; otherIdx2 = 2; }
                            else if (sixFaceNodes[2] != sixFaceNodes[3] && sixFaceNodes[2] != sixFaceNodes[4] && sixFaceNodes[2] != sixFaceNodes[5]) { alphaIdx = 2; otherIdx1 = 0; otherIdx2 = 1; }
                            if (sixFaceNodes[3] != sixFaceNodes[0] && sixFaceNodes[3] != sixFaceNodes[1] && sixFaceNodes[3] != sixFaceNodes[2]) { betaIdx = 3; }
                            else if (sixFaceNodes[4] != sixFaceNodes[0] && sixFaceNodes[4] != sixFaceNodes[1] && sixFaceNodes[4] != sixFaceNodes[2]) { betaIdx = 4; }
                            else if (sixFaceNodes[5] != sixFaceNodes[0] && sixFaceNodes[5] != sixFaceNodes[1] && sixFaceNodes[5] != sixFaceNodes[2]) { betaIdx = 5; }
                            
                            // x_i is the original node, and x_j is one of the
                            // two remaining shared nodes of the current edge...
                            vec3 x_i = dsc->get_pos(nit.key());
                            vec3 x_j;
                            if (x_i == sixFaceNodes[otherIdx1]) {
                                x_j = sixFaceNodes[otherIdx2];
                            } else {
                                x_j = sixFaceNodes[otherIdx1];
                            }
                            vec3 edgeVec = x_j - x_i;

                            // Compute angles; cotangent is tan(M_PI_2 - angle) (http://stackoverflow.com/questions/3738384/stable-cotangent):
                            real alphaAngle = Util::angle<real, vec3>(sixFaceNodes[alphaIdx], sixFaceNodes[otherIdx1], sixFaceNodes[otherIdx2]);
                            real betaAngle = Util::angle<real, vec3>(sixFaceNodes[betaIdx], sixFaceNodes[otherIdx1], sixFaceNodes[otherIdx2]);
                            edgeVec *= (tan(M_PI_2 - alphaAngle) + tan(M_PI_2 - betaAngle));
                            sum += edgeVec;

                            // Add area of both faces to total area so far:
                            totalArea += Util::area<real, vec3>(sixFaceNodes[0], sixFaceNodes[1], sixFaceNodes[2]);
                            totalArea += Util::area<real, vec3>(sixFaceNodes[3], sixFaceNodes[4], sixFaceNodes[5]);
                        }
                    }
                    totalArea /= 2; // (We added each face's area twice on purpose, so divide by 2 here)
                    sum /= (4 * totalArea);
                    cout << "NODE: " << ++i << " | k = " << sum.length() << endl;
                }
            }
            break;
        }
        case 'p': { // Split all interface faces whose area is greater than the average
            vel_fun->split_larger_than_average_interface_faces(*dsc);
            break;
        }
        case 'o':
            vel_fun->print_face_speed_stats(*dsc, false);
            break;
        case 'O':
            vel_fun->print_face_speed_stats(*dsc, true);
            break;
        case 'm':
            std::cout << "MOVE" << std::endl;
            vel_fun->take_time_step(*dsc);
            painter->update(*dsc);
            break;
        case 'M': {
            int steps_to_take = 200;
            bool save_obj_and_screenshot_after_every_step = true;
            std::cout << "MOVE (" << steps_to_take << " steps)" << std::endl;
            real total_time = 0.;
            for (int i = 0; i < steps_to_take; i++) {
                std::cout << "Time step " << (i + 1) << "/" << steps_to_take << ", " << total_time << " sec elapsed so far" << std::endl;
                vel_fun->take_time_step(*dsc);
                total_time += vel_fun->get_deform_time();
                if (i % 4 == 3) {
                    vel_fun->print_face_speed_stats(*dsc, true);
                }
                //if (i % 100 == 74) {
                //    vel_fun->split_larger_than_average_interface_faces(*dsc);
                //}
                if (save_obj_and_screenshot_after_every_step) {
                    painter->update(*dsc);
                    // Export screenshot:
                    std::cout << "TAKING SCREEN SHOT " << (i+1) << std::endl;
                    painter->set_view_position(per_step_screenshot_pos);
                    painter->save_painting("LOG_PAPER", i+1);
                    // Export .obj mesh:
                    std::cout << "EXPORTING SURFACE MESH " << (i+1) << std::endl;
                    std::string filename("LOG_PAPER/mesh" + std::string(Util::concat4digits("_", i+1)) + ".obj");
                    std::vector<vec3> points;
                    std::vector<int> faces;
                    dsc->extract_surface_mesh(points, faces);
                    is_mesh::export_surface_mesh(filename, points, faces);
                }
            }
            painter->update(*dsc);
            std::cout << "Move complete.  Total elapsed time: " << total_time << " sec" << std::endl;
            break;
        }
        case 'r':
            std::cout << "RELOAD MODEL" << std::endl;
            load_model(model_file_name, dsc->get_avg_edge_length());
            break;
        case 't':
            std::cout << "TEST VELOCITY FUNCTION" << std::endl;
            vel_fun->test(*dsc);
            painter->update(*dsc);
            break;
        case '\t':
            painter->switch_display_type();
            painter->update(*dsc);
            break;
        case 's':
            std::cout << "TAKING SCREEN SHOT" << std::endl;
            painter->set_view_position(camera_pos);
            painter->save_painting("LOG");
            break;
        case 'e':
        {
            std::cout << "EXPORTING MESH" << std::endl;
            std::string filename("data/mesh.dsc");
            std::vector<vec3> points;
            std::vector<int> tets;
            std::vector<int> tet_labels;
            dsc->extract_tet_mesh(points, tets, tet_labels);
            is_mesh::export_tet_mesh(filename, points, tets, tet_labels);
        }
            break;
        case 'i':
        {
            std::cout << "EXPORTING SURFACE MESH" << std::endl;
            std::string filename("data/mesh.obj");
            std::vector<vec3> points;
            std::vector<int> faces;
            dsc->extract_surface_mesh(points, faces);
            is_mesh::export_surface_mesh(filename, points, faces);
        }
            break;
        case 'I': {
            std::cout << "EXPORTING SPEED AND CURVATURE DATA" << std::endl;
            vel_fun->write_data_files(*dsc);
        }
            break;
        case '+':
        {
            real velocity = std::min(vel_fun->get_velocity() + 1., 100.);
            vel_fun->set_velocity(velocity);
            update_title();
        }
            break;
        case '-':
        {
            real velocity = std::max(vel_fun->get_velocity() - 1., 0.);
            vel_fun->set_velocity(velocity);
            update_title();
        }
            break;
        case '.':
        {
            real discretization = std::min(dsc->get_avg_edge_length() + 0.5, 100.);
            dsc->set_avg_edge_length(discretization);
            update_title();
        }
            break;
        case ',':
        {
            real discretization = std::max(dsc->get_avg_edge_length() - 0.5, 1.);
            dsc->set_avg_edge_length(discretization);
            update_title();
        }
            break;
        case '<':
        {
            real accuracy = std::min(vel_fun->get_accuracy() + 1., 100.);
            vel_fun->set_accuracy(accuracy);
            update_title();
        }
            break;
        case '>':
        {
            real accuracy = std::max(vel_fun->get_accuracy() - 1., 1.);
            vel_fun->set_accuracy(accuracy);
            update_title();
        }
            break;
    }
}

void UI::visible(int v)
{
    if(v==GLUT_VISIBLE)
        glutIdleFunc(animate_);
    else
        glutIdleFunc(0);
}

void UI::stop()
{
    if(RECORD && basic_log)
    {
        basic_log->write_message("MOTION STOPPED");
        basic_log->write_log(*dsc);
        basic_log->write_log(*vel_fun);
        basic_log->write_timings(*vel_fun);
        
        std::vector<vec3> points;
        std::vector<int> faces;
        std::vector<int> tets;
        std::vector<int> tet_labels;
        dsc->extract_tet_mesh(points, tets, tet_labels);
        is_mesh::export_tet_mesh(basic_log->get_path() + std::string("/mesh.dsc"), points, tets, tet_labels);
        points.clear();
        dsc->extract_surface_mesh(points, faces);
        is_mesh::export_surface_mesh(basic_log->get_path() + std::string("/mesh.obj"), points, faces);
        basic_log = nullptr;
    }
    
    CONTINUOUS = false;
    painter->update(*dsc);
    update_title();
    glutPostRedisplay();
}

void UI::start(const std::string& log_folder_name)
{
    if(RECORD)
    {
        basic_log = std::unique_ptr<Log>(new Log(log_path + log_folder_name));
        painter->set_view_position(camera_pos);
        painter->save_painting(log_path, vel_fun->get_time_step());
        basic_log->write_message(vel_fun->get_name().c_str());
        basic_log->write_log(*vel_fun);
        basic_log->write_log(*dsc);
    }
    
    painter->update(*dsc);
    update_title();
    glutPostRedisplay();
}

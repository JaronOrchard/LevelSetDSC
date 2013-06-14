//
//  user_interface.cpp
//  3D_DSC
//
//  Created by Asger Nyman Christiansen on 2/21/13.
//  Copyright (c) 2013 Asger Nyman Christiansen. All rights reserved.
//

#include "user_interface.h"

void _check_gl_error(const char *file, int line)
{
    GLenum err (glGetError());
    
    while(err!=GL_NO_ERROR) {
        std::string error;
        
        switch(err) {
            case GL_INVALID_OPERATION:      error="INVALID_OPERATION";      break;
            case GL_INVALID_ENUM:           error="INVALID_ENUM";           break;
            case GL_INVALID_VALUE:          error="INVALID_VALUE";          break;
            case GL_OUT_OF_MEMORY:          error="OUT_OF_MEMORY";          break;
            case GL_INVALID_FRAMEBUFFER_OPERATION:  error="INVALID_FRAMEBUFFER_OPERATION";  break;
        }
        
        std::cerr << "GL_" << error.c_str() <<" - "<<file<<":"<<line<<std::endl;
        err=glGetError();
    }
}

#define check_gl_error() _check_gl_error(__FILE__,__LINE__)


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

void motion_(int x, int y){
    UI::get_instance()->motion(x,y);
}

void mouse_(int button, int state, int x, int y)
{
    UI::get_instance()->mouse(button, state, x, y);
}

void animate_(){
    UI::get_instance()->animate();
}

UI* UI::instance = NULL;

UI::UI(int &argc, char** argv)
{
    instance = this;
	WIN_SIZE_X = 1000;
    WIN_SIZE_Y = 1000;

    glutInit(&argc, argv);
    glutInitWindowSize(WIN_SIZE_X,WIN_SIZE_Y);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH | GLUT_MULTISAMPLE);
    glutCreateWindow("");
    
    glutDisplayFunc(display_);
    glutKeyboardFunc(keyboard_);
	glutIgnoreKeyRepeat(true);
    glutVisibilityFunc(visible_);
    glutReshapeFunc(reshape_);
	glutMouseFunc(mouse_);
	glutMotionFunc(motion_);
	glutIdleFunc(animate_);
    
    glEnable(GL_MULTISAMPLE);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glEnable(GL_POINT_SMOOTH);
    glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
    
	glEnable(GL_DEPTH_TEST);
    glDepthMask(GL_TRUE);
    
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    
	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
    
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glShadeModel(GL_FLAT);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
    
	glEnable (GL_BLEND);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
        
    dsc = (DSC<GELTypes>*) NULL;
    
    if(argc > 1)
    {
        QUIT_ON_COMPLETION = true;
        CONTINUOUS = true;
        RECORD = true;
        
        Util::ArgExtracter ext(argc, argv);
        ext.extract("nu", VELOCITY);
        ext.extract("delta", DISCRETIZATION);
        ext.extract("alpha", ACCURACY);
    }
    else {
        VELOCITY = 1.;
        DISCRETIZATION = 25.;
        ACCURACY = 5.;
        
        CONTINUOUS = false;
        RECORD = true;
        QUIT_ON_COMPLETION = false;
    }
    update_title();
    check_gl_error();
}

void UI::update_title()
{
    std::ostringstream oss;
    oss << "3D DSC\t";
    if(dsc)
    {
        oss << dsc->get_title();
    }
    oss << " (Nu = " << VELOCITY << ", Delta = " << DISCRETIZATION << ", Alpha = " << ACCURACY << ")";
    std::string str(oss.str());
    glutSetWindowTitle(str.c_str());
}

void UI::display()
{
    if (glutGet(GLUT_WINDOW_WIDTH) != WIN_SIZE_X || glutGet(GLUT_WINDOW_HEIGHT) != WIN_SIZE_Y) {
        return;
    }
    draw();
    update_title();
    
    if(dsc && CONTINUOUS)
    {
        bool finished = dsc->take_time_step();
        if (finished)
        {
            stop();
            if (QUIT_ON_COMPLETION) {
                exit(0);
            }
        }
    }
    
    draw();
    update_title();
    check_gl_error();
}

void UI::reshape(int width, int height)
{
    WIN_SIZE_X = width;
    WIN_SIZE_Y = height;
    
    if (dsc) {
        view_ctrl->reshape(WIN_SIZE_X,WIN_SIZE_Y);
    }
}

void UI::animate()
{
    glutPostRedisplay();
}

void UI::mouse(int button, int state, int x, int y)
{
    if (dsc) {
        CGLA::Vec2i pos(x,y);
        if (state == GLUT_DOWN)
        {
            if (button == GLUT_LEFT_BUTTON)
                view_ctrl->grab_ball(GLGraphics::ROTATE_ACTION,pos);
            else if (button == GLUT_MIDDLE_BUTTON)
                view_ctrl->grab_ball(GLGraphics::ZOOM_ACTION,pos);
            else if (button == GLUT_RIGHT_BUTTON)
                view_ctrl->grab_ball(GLGraphics::PAN_ACTION,pos);
        }
        else if (state == GLUT_UP)
            view_ctrl->release_ball();
    }
}

void UI::motion(int x, int y)
{
    if (dsc) {
        CGLA::Vec2i pos(x,y);
        view_ctrl->roll_ball(pos);
    }
}

void UI::keyboard(unsigned char key, int x, int y) {
    switch(key) {
        case '\033':
            stop();
            exit(0);
            break;
        case '0':
            stop();
            break;
        case '1':
            motion1();
            break;
        case '2':
            motion2();
            break;
        case '3':
            motion3();
            break;
        case '4':
            motion4();
            break;
        case '5':
            motion5();
            break;
        case '6':
            motion6();
            break;
        case ' ':
            if(!CONTINUOUS)
            {
                std::cout << "MOTION STARTED" << std::endl;
            }
            else {
                std::cout << "MOTION PAUSED" << std::endl;
            }
            CONTINUOUS = !CONTINUOUS;
            break;
        case 'm':
            if(dsc)
            {
                std::cout << "MOVE" << std::endl;
                dsc->take_time_step();
            }
            break;
        case 't':
            if(dsc)
            {
                std::cout << "TEST" << std::endl;
                dsc->test();
            }
            break;
        case '\t':
            if(dsc)
            {
                switch_display_type();
            }
            break;
        case 's':
            if(dsc)
            {
                std::cout << "TAKING SCREEN SHOT" << std::endl;
                Painter<GELTypes>::save_painting(WIN_SIZE_X, WIN_SIZE_Y, "LOG");
            }
            break;
        case '+':
            if(!dsc)
            {
                VELOCITY = std::min(VELOCITY + 1., 100.);
                update_title();
            }
            break;
        case '-':
            if(!dsc)
            {
                VELOCITY = std::max(VELOCITY - 1., 0.);
                update_title();
            }
            break;
        case '.':
            if(!dsc)
            {
                DISCRETIZATION = std::min(DISCRETIZATION + 0.5, 100.);
                update_title();
            }
            break;
        case ',':
            if(!dsc)
            {
                DISCRETIZATION = std::max(DISCRETIZATION - 0.5, 1.);
                update_title();
            }
            break;
        case '<':
            if(!dsc)
            {
                ACCURACY = std::min(ACCURACY + 1., 100.);
                update_title();
            }
            break;
        case '>':
            if(!dsc)
            {
                ACCURACY = std::max(ACCURACY - 1., 1.);
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

void UI::draw()
{
    Painter<GELTypes>::begin();
    if (dsc)
    {
        view_ctrl->set_gl_modelview();
        
        DeformableSimplicialComplex<GELTypes> *complex = dsc->get_complex();
        Painter<GELTypes>::draw_complex(complex);
        if(RECORD && CONTINUOUS)
        {
            Painter<GELTypes>::save_painting(WIN_SIZE_X, WIN_SIZE_Y, dsc->get_log_path(), dsc->get_time_step());
        }
    }
    Painter<GELTypes>::end();
}

void UI::stop()
{
    if(dsc)
    {
        draw();
        delete dsc;
        delete view_ctrl;
        dsc = (DSC<GELTypes>*) NULL;
    }
}


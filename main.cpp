#include <iostream>
#include "points_cloud.h"

#define GLUT_DISABLE_ATEXIT_HACK
#include "C:\Program Files\CodeBlocks\MinGW\x86_64-w64-mingw32\include\GL\glu.h"
#include "C:\Program Files\CodeBlocks\MinGW\x86_64-w64-mingw32\include\GL\glut.h"
#include <windows.h>

double *alpha, *beta;

bool is_one(std::vector<double> vec){
    return vec[0]==1;
}

double line_x(double t, double v){

    return t*std::cos(t);
}

double line_y(double t, double v){

    return t*std::sin(t);
}

double line_z(double t, double v){

    return std::sin(t);
}

PointCloud* example_1(){
    gluOrtho2D(0, 200.0, 0, 200.0);
    PointCloud *pc = new PointCloud();
    std::vector<double> vec{};
    for(double i = 50; i < 150; i += 2){
        for(double j = 50; j < 150; j += 2){
            if( std::abs(i - 100)<15 && std::abs(j - 100)<15 ){
                vec.push_back(i);
                vec.push_back(j);
                pc->add_point(vec);
                vec = {};
            }

            if( std::abs(i-50)<35 && std::abs(j-50)<35 ){
                vec.push_back(i);
                vec.push_back(j);
                pc->add_point(vec);
                vec = {};
            }
        }
    }
    return pc;
}


PointCloud* example_2(){
    gluOrtho2D(0, 500.0, 0, 500.0);
    PointCloud *pc = new PointCloud();
    std::vector<double> vec{};
    for(double i = 50; i < 150; i += 2){
        for(double j = 50; j < 150; j += 2){
            if( std::abs(i - 100)<60  && std::abs(i - 100)>30){
                vec.push_back(i);
                vec.push_back(j);
                pc->add_point(vec);
                vec = {};
            }
            if(std::sqrt((i-100)*(i-100)+(j-100)*(j-100)) == 40){
                vec.push_back(i);
                vec.push_back(j);
                pc->add_point(vec);
                vec = {};
            }
        }
    }
    return pc;
}

PointCloud* henon(){
    gluOrtho2D(-1.5, 1.5, -1, 1.0);
    PointCloud *pc = new PointCloud();
    std::vector<double> vec{};
    double x_p{}, y_p{}, x_n{},y_n{};

    for(double i = 0; i < 1; i += 0.01){
        for(double j = 0; j < 1; j += 0.01){
            x_p = i;
            y_p = j;
            for(int k = 0; k < 15; ++k){
                x_n = 1 - 1.4*x_p*x_p+y_p;
                y_n = 0.3*x_p;
                x_p = x_n;
                y_p = y_n;
            }
            vec.push_back(x_n);
            vec.push_back(y_n);
            pc->add_point(vec);
            vec = {};
        }
    }
    return pc;
}

PointCloud* tinkerbell(){
    gluOrtho2D(-1.8, 1.8, -2, 1.2);
    PointCloud *pc = new PointCloud();
    std::vector<double> vec{};
    double x_p{}, y_p{}, x_n{},y_n{};

    for(double i = -0.5; i < 0.5; i += 0.008){
        for(double j = -0.5; j < 0.5; j += 0.008){
            x_p = i;
            y_p = j;
            for(int k = 0; k < 15; ++k){
                x_n = x_p*x_p - y_p*y_p+0.9*x_p-0.6013*y_p;
                y_n = 2*x_p*y_p+2*x_p+0.5*y_p;
                x_p = x_n;
                y_p = y_n;
            }
            vec.push_back(x_n);
            vec.push_back(y_n);
            pc->add_point(vec);
            vec = {};
        }
    }
    return pc;
}

void display() {
    gluOrtho2D(-20, 20, -20, 50);
	glClear(GL_COLOR_BUFFER_BIT);

    PointCloud *pc = txt_to_points_cloud("lorenz_x_z.txt");

    std::vector<HOMOLOGY_TYPE> homo = pc->homo_dim(*alpha,*beta);

    glColor3f(0.0, 0.0, 1.0);

	glBegin(GL_POINTS);
        for(int i = 0; i < pc->cloud_size(); ++i){

                switch(homo[i]){
                    case Z_Z:
                        glColor3f(1,0,1);//pink
                        break ;
                    case Z:
                        glColor3f(0.0, 0.0, 1.0);//blue
                        break ;
                    case Z2:
                        glColor3f(0.0, 1.0, 0.0);//green
                        break ;
                    case Z3:
                        glColor3f(1.0, 0.0, 0.0);//red
                        break ;
                    case ANOTHER:
                        glColor3f(0.5, 0.5, 0.0);//gray
                        break ;
                    case ISOLATED:
                        glColor3f(0.9, 0.7, 0);//yellow
                        break ;
                    }
                glVertex2d((pc->get_points())[i][0], (pc->get_points())[i][1]);
            }
	glEnd();
	glFlush();
}

void myinit() {
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glColor3f(1.0, 0.0, 0.0);
	glPointSize(1.2);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
}

int main(int argc, char** argv)
{
    double a, b;

    std::cin>>a>>b;

    alpha = &a;
    beta = &b;

    glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(0, 0);
	glutCreateWindow("Points");
	glutDisplayFunc(display);

	myinit();
	glutMainLoop();

    return 0;
}

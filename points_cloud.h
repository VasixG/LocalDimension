#ifndef POINTS_CLOUD_H_INCLUDED
#define POINTS_CLOUD_H_INCLUDED

#include "vector"
#include <algorithm>
#include <cmath>
#include <ostream>
#include <iostream>
#include <fstream>
#include <string>

enum ERR_CODE{DIM_ERR, EMPTY};
enum HOMOLOGY_TYPE{ISOLATED, Z_Z, Z, Z2, Z3, ANOTHER};

typedef bool (*filt)(std::vector<double> point);
typedef double (*xfunc)(double t, double v);
typedef double (*yfunc)(double t, double v);
typedef double (*zfunc)(double t, double v);

double inner_prod(const std::vector<double> v1, const std::vector<double> v2){

    //if(v1.size() != v2.size())throw DIM_ERR;
    double dotprod{0};
    for(auto it1 = v1.begin(), it2 = v2.begin(); it1 != v1.end(); ++it1, ++it2){
        dotprod+=(*it1)*(*it2);
    }
    return dotprod;
}

std::vector<double> operator - (std::vector<double> v1, std::vector<double> v2){

    if(v1.size() != v2.size())throw DIM_ERR;
    std::vector<double> v{};
    for(auto it1 = v1.begin(), it2 = v2.begin(); it1 != v1.end(); ++it1, ++it2){
        v.push_back((*it1)-(*it2));
    }
    return v;
}

std::vector<double> operator - (std::vector<double> v1){

    std::vector<double> v{};
    for(auto it1 = v1.begin(); it1 != v1.end(); ++it1){
        v.push_back(-(*it1));
    }
    return v;
}

std::vector<double> operator + (std::vector<double> v1, std::vector<double> v2){

    if(v1.size() != v2.size())throw DIM_ERR;
    std::vector<double> v{0};
    for(auto it1 = v1.begin(), it2 = v2.begin(); it1 != v1.end(); ++it1, ++it2){
        v.push_back((*it1)+(*it2));
    }
    return v;
}

double norm(std::vector<double> vec){
    return std::sqrt(inner_prod(vec,vec));
}

double distance(std::vector<double> p1, std::vector<double> p2){
    return norm(p1-p2);
}

class Dict{
public:
    Dict(std::vector<double>* _vect, double _value): vect{_vect}, value{_value}{}
    ~Dict(){}

    double get_value(){
        return value;
    }

    std::vector<double>* get_vec(){
        return vect;
    }

    bool operator < (const Dict& dict1){
        if(value < dict1.value){
            return true;
        }
        return false;
    }

private:
    std::vector<double>* vect;
    double value;

};

class Homopoints{
public:
    Homopoints(std::vector<double>* _vect, HOMOLOGY_TYPE _value): vect{_vect}, value{_value}{}
    ~Homopoints(){}

    double get_value(){
        return value;
    }

    std::vector<double>* get_vec(){
        return vect;
    }


private:
    std::vector<double>* vect;
    HOMOLOGY_TYPE value;

};


class PointCloud{
public:
    PointCloud(std::vector<std::vector<double>> _points = {{}}):points{_points}{
        int e_dim = points[0].size();
        for(auto it = points.begin(); it != points.end(); ++it){
            if(it->size() != e_dim)std::cout<<"DIM_ERR \n";
        }
    };

    virtual ~PointCloud(){}


    int eucl_dim() const{
        if(is_empty())return 0;
        return points[0].size();
    }


    int cloud_size() const{
        return points.size();
    }


    bool is_empty() const{
        return (points.size()==1 && points[0].size() == 0);
    }


    std::vector<std::vector<double>> get_points() const{
        return points;
    }


    PointCloud& add_point(std::vector<double> p){
        if((p.size() != eucl_dim()) && (!is_empty()))throw DIM_ERR;
        if(is_empty())points.pop_back();
        points.push_back(p);
        return *this;
    }


    PointCloud& filt_points(filt f){
        auto newEnd = std::remove_if(points.begin(), points.end(), f);

        points.erase(newEnd, points.end());
        return *this;
    }

    PointCloud* eps_area_of_point(double eps, std::vector<double> point){
        PointCloud *area = new PointCloud();

        for(auto it = points.begin(); it != points.end(); ++it)if(norm(point-(*it))<eps) area->add_point((*it));

        return area;
    }

    friend std::ostream& operator << (std::ostream &stream, PointCloud& v){
        stream<<"[";
        int i = 0, j = 0;

        for(; i<v.cloud_size() - 1; ++i){
            stream<<"{ ";
            for(j = 0; j < v.eucl_dim() - 1; ++j){
                stream<<v.points[i][j]<<", ";
            }
            stream<<v.points[i][j];
            stream<<"}, ";
        }
        stream<<"{ ";
        for(int j = 0; j < v.eucl_dim() - 1; ++j){
            stream<<v.points[i][j]<<", ";
        }
        stream<<v.points[i][j];
        stream<<"} ";

        stream<<"]";
        return stream;
    }

    double stoh_dim_of_point(double eps, std::vector<double> point){

        double prod_s{1};

        for(auto it = points.begin(); it != points.end(); ++it){
            if((*it) != point){ prod_s *= norm((*it)-point);}
        }

        std::cout<<"Euclid dim: "<< eucl_dim()<<"\n";
        std::cout<<"r^k: "<< std::pow(eps, eucl_dim())<<"\n";
        std::cout<<"prod{s}: "<< prod_s<<"\n";
        std::cout<<"log(prod_s/r^k): "<< std::log(prod_s/(std::pow(eps, eucl_dim())))<<"\n";


        return -((double)(eucl_dim()))/(std::log(prod_s/(std::pow(eps, eucl_dim()))));
    }

    HOMOLOGY_TYPE homo_dim_of_the_point(double alpha, double beta, std::vector<double> point){
        int holes{0};
        std::vector<double>id = {1, 0};
        std::vector<Dict>dict_vec_angle{};
        for(auto i = points.begin(); i != points.end(); ++i){

            std::vector<double> cp{(*i) - point};

            if(std::abs(cp[0]) <= (alpha+beta)/2 && std::abs(cp[1]) <= (alpha+beta)/2){

                double norm_cp{norm(cp)};
                double phi{0};
                if(norm_cp <= (alpha+beta)/2 && norm_cp >= (beta-alpha)/2){
                    if(cp[1] < 0){
                        phi = std::acos(inner_prod(-cp, id)/norm_cp) + M_PI;
                    }else{
                        phi = std::acos(inner_prod(cp, id)/norm_cp);
                    }

                    dict_vec_angle.push_back(Dict(&(*i), phi));

                }
            }
        }

        if(dict_vec_angle.size() == 0){
            return ISOLATED;
        }

        std::sort(dict_vec_angle.begin(), dict_vec_angle.end());

        for(auto j = dict_vec_angle.begin();j != dict_vec_angle.end()-1; ++j){
            if(distance(*(j->get_vec()),*((j+1)->get_vec()) ) > alpha )++holes;
        }
        if(distance(*(dict_vec_angle.begin()->get_vec()),*((dict_vec_angle.end()-1)->get_vec())) > alpha )++holes;

        switch(holes){
        case 0:return Z_Z;
        case 1: return Z;
        case 2: return Z2;
        case 3: return Z3;
        }

        return ANOTHER;
    }


    std::vector<HOMOLOGY_TYPE> homo_dim(double alpha, double beta){
        std::vector<HOMOLOGY_TYPE> homo_points{};

        for(auto i = points.begin(); i != points.end(); ++i){
            homo_points.push_back(homo_dim_of_the_point(alpha, beta, (*i)));
        }

        return homo_points;
    }


private:
    std::vector<std::vector<double>> points;


};


std::vector<double> curve_point2D(double t, xfunc x_t, yfunc y_t){
    return std::vector<double>{x_t(t, 0), y_t(t, 0)};
}


PointCloud* parametrizationCurve2D(xfunc x_t, yfunc y_t, double start_t, double end_t, double step_t){
    PointCloud *pc = new PointCloud();

    for(double t = start_t; t < end_t; t+=step_t){
        pc->add_point(curve_point2D(t, x_t, y_t));
    }

    return pc;
}


std::vector<double> curve_point3D(double t, xfunc x_t, yfunc y_t, zfunc z_t){
    return std::vector<double>{x_t(t, 0), y_t(t, 0), z_t(t, 0)};
}


PointCloud* parametrizationCurve3D(xfunc x_t, yfunc y_t, zfunc z_t, double start_t, double end_t, double step_t){
    PointCloud *pc = new PointCloud();

    for(double t = start_t; t < end_t; t+=step_t){
        pc->add_point(curve_point3D(t, x_t, y_t, z_t));
    }

    return pc;
}

std::vector<double> surface_point3D(double t, double v, xfunc x_t, yfunc y_t, zfunc z_t){
    return std::vector<double>{x_t(t, v), y_t(t, v), z_t(t, v)};
}


PointCloud* parametrizationSurface3D(xfunc x_t, yfunc y_t, zfunc z_t, double start_t, double end_t, double step_t, double start_v, double end_v, double step_v){
    PointCloud *pc = new PointCloud();

    for(double t = start_t; t < end_t; t+=step_t){
        for(double v = start_v; v < end_v; v += step_v){
            pc->add_point(surface_point3D(t, v, x_t, y_t, z_t));
        }
    }

    return pc;
}

PointCloud* txt_to_points_cloud(std::string path){
    PointCloud *pc = new PointCloud();
    std::ifstream txt;
    txt.open(path);

    std::string mystring;
    std::vector<double> vect{};
    int k = 0;
    if (txt.is_open() ) {
        while ( txt.good() ){
         txt >> mystring;
         if(k == 0){
            vect.push_back(std::stod(mystring));
            ++k;
         }else{
            vect.push_back(std::stod(mystring));
            pc->add_point(vect);
            vect = {};
            k = 0;
         }
        }
    }


    return pc;
}


#endif // POINTS_CLOUD_H_INCLUDED

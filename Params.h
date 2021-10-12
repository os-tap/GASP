#pragma once

#include <cassert>

#define M_PI 3.14159265358979323846
#include <cmath>

#include <string>
#include <iostream>
#include <fstream>

#include <fmt/core.h>
#include <nlohmann/json.hpp>

#include <SPLINTER/datatable.h>
#include <SPLINTER/bspline.h>
#include <SPLINTER/bsplinebuilder.h>

#include "Point.hpp"


namespace ps {
    class Params;
}

using json = nlohmann::json;

class ps::Params {

public:

    Params();
    void Read();
    void Load(json j);
    json GetFromFile();
    void Print();

    const int screen_width = 400;
    const int screen_height = 700;
    const int screen_bottom_gap = 0;

    std::string csv_folder;

    double screen_proportion, screen_x_proportion, screen_y_proportion;

    double area_beg, area_end, area_size, area_center, area_height;
    double stream_beg, stream_end, stream_width, stream_radius;
    double L, scale, DSR;

    double burn_fix_a, burn_fix_b;

    double burn_radius, burn_radius_2;
    double burn_radius_cross, burn_radius_2_cross;

    double base_speed, burn_speed, const_speed;
    double const_speed_x, const_speed_z;

    double emitter_begin;
    int base_particles;
    double particles_dist;
    double base_particles_NR2;

    int iterations;
    double iterate_burn_radius;
    double iterate_speed;
    double iterate_const;
    double iterate_const_x, iterate_const_z;
    int iterate_particles;

    int burn_time, sage_time, wave_time, burn_fade;


    int frontline_window_steps;
    double frontline_window_size;

    int frontline_spline_steps;
    double frontline_spline_alpha;
    bool move_normal, move_speed;
    bool curve_by_radius;

    int frontline_stencil_h, frontline_radius_h;


    double refract_coef, refract_offset;

    int svm_count, print_count, display_count;

    bool frontline_cross_chunk, calc_cross;
    double frontline_cross_multipler, frontline_cross_area, frontline_cross_radius, frontline_cross_radius_2;
    double frontline_cross_border;
    double frontline_cross_spline_alpha;


    bool scale_burn;


    double curve_start, curve_end;
    double curve_burn_coef;


    int grid_curve_calc;
    double grid_curve_area;
    bool refill;
    double grid_count_gap;

    SPLINTER::BSpline curve_spline{ 1 };

    void set_curve_spline(SPLINTER::BSpline spline, double _start, double _end) {
        curve_spline = spline;
        curve_start = _start;
        curve_end = _end;
    }
    void clear_curve_spline() {
        set_curve_spline(SPLINTER::BSpline{ 1 }, area_end, area_beg);
        /*curve_spline = SPLINTER::BSpline {1};
        curve_start = area_end;
        curve_end = area_beg;*/
    }

    double get_curvature(const double x) const
    {
        if (x < curve_start || x > curve_end) return 0;

        //std::cout << curve_start;

        SPLINTER::DenseVector xd(1);
        xd(0) = x;
        double c = curve_spline.eval(xd);

        return c;

    }



    /*double get_curvature(const double x) const {
        if (x < curve_start || x > curve_end) return 0;
        //double c = frontline_curve[lround((x - curve_start) / curve_width * (frontline_curve.size() - 1))];
        return sqrt(frontline_curve[lround((x - curve_start) / curve_width * (frontline_curve.size() - 1))]);
    }*/



    //std::string swarm_params() const;
//    std::string frontline_params() const;


    double make_radius_cross_fix(double radius) const {
        return radius / (1 - burn_fix_a / pow(base_particles_NR2* radius* radius, burn_fix_b));
    }
    double get_burn_radius(double x_cord) const {

        //return burn_radius_cross;
        return burn_radius_cross * (1 + get_curvature(x_cord) * /*system_speed(x_cord)*/  curve_burn_coef * scale_burn);


        /*double br = burn_radius_cross;
        for (const auto x_kink : frontline_kinks) {
            if (scale_burn && fabs(x_kink - x_cord) < scale_burn_size) {
                br *= 1 + pow(1 - fabs(x_kink - x_cord) / scale_burn_size, scale_burn_pow) * particle_speed(x_cord) * particle_speed(x_cord) * scale_burn_multipler;
                break;
            }
        }
        return br;*/


        /*double br = burn_radius_cross;
        if (scale_burn && fabs(x_cord) < scale_burn_size) {
            br *= 1 + pow(1 - fabs(x_cord) / scale_burn_size, scale_burn_pow) * scale_burn_multipler;
        }
        return br;*/
        //return burn_radius_cross;
        //double ir = 0.2;
        //return burn_radius_cross * (1 + (fabs(x_cord) < 0.2) * pow(1 - fabs(x_cord) / 0.2 * 1, 2) * 1.5);
        //return burn_radius_cross * scale_burn * (1 + (fabs(x_cord) < scale_burn_size) * pow(1 - fabs(x_cord) / scale_burn_size, scale_burn_pow) * scale_burn_multipler);
        //return burn_radius_cross * scale_burn * (1 + (fabs(x_cord) < scale_burn_size) * sqrt(1 - fabs(x_cord) / scale_burn_size) * scale_burn_multipler);
        //return scale_burn ? burn_radius_cross * (1 + (fabs(x_cord) < 0.2)) : burn_radius_cross;
    }


    double stream_function(double x) const {
        return (this->*stream_function_p)(x);
    }

    void SetStream(int i) {
        assert(i && i <= 4);
        stream_function_p = streams[i];
        Read();
    }

private:
    typedef double(Params::* stream)(double) const;
    stream streams[5]{
            0,
            &Params::linear_stream,
            &Params::log_stream,
            &Params::x2_stream,
            &Params::const_stream,
    };
    stream stream_function_p = &Params::x2_stream;

    double linear_stream(const double x) const {
        return 1 - fabs(from_center(x)) / stream_radius;
    }
    double log_stream(const double x) const {
        return log(stream_radius + 1 - fabs(from_center(x))) / log(stream_radius + 1);
    }
    double x2_stream(const double x) const {
        return 1 - from_center(x) * from_center(x) / stream_radius / stream_radius;
    }
    double const_stream(const double x) const  {
        return 0.5;
    }
    double rising_stream(const double x) const {
        return burn_speed * burn_radius + (x - stream_beg) / stream_width;
    }


public:
    double from_center(const double x) const{
        return x - area_center;
    }
    double screen_to_area_x(const int x) const {
        return x * screen_x_proportion + area_beg;
    }
    double screen_to_area_y(const int y) const {
        return area_height - y * screen_y_proportion;
    }
    double center_percentage(const double x) const {
        return 1 - fabs(1 - (x - area_beg) / area_size * 2);
    }
    double system_speed(const double x) const {
        return base_speed * stream_function(x) + const_speed;
    }
    double burn_speed_fu(const double x) const {
        return burn_speed * stream_function(x);
    }
    double profile_speed(const double x) const {
        return iterate_speed * stream_function(x);
    }
    double particle_speed(const double x) const {
        return profile_speed(x) + iterate_const;
    }


    double particle_speed_x(const double x, const double z) const {
        //return iterate_const_x * z;
        return iterate_const_x * (1 - 2 * z) * (x - x * x);
    }
    double particle_speed_z(const double x, const double z) const {
        //return iterate_const_z * -x;
        return -iterate_const_z * (1 - 2 * x) * (z - z * z) + particle_speed(x);
    }

    double refract_func(const double x) const {
        return x * refract_coef + refract_offset;
    }




};


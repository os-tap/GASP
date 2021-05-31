#pragma once

#include <vector>
#include <string>

#include <fmt/core.h>
#include <iostream>
#include <fstream>


//#include <SPLINTER/datatable.h>
//#include <SPLINTER/bspline.h>
//#include <SPLINTER/bsplinebuilder.h>

#include "Spline.h"


#include "Params.h"
#include "Particle.h"


namespace ps {

    class Frontline
    {
        const Params* P;
    public:
        Frontline(const Params& P_) : P(&P_) { Init(); }
        void Init();


        void Calc(const std::vector <Particle*>& particle_list);
        void WindowMiddle(const std::vector <Particle*>& particle_list);
        void SplineSmooth(double alpha);
        void SplineSmooth2(double alpha);
        void FivePointStencil(int h_div);
        void CalcNormal();
        void CalcError();
        void CalcRadius(int h_div);
        void CalcCurve();
        void SetCrosses(const std::vector <double>& crosses);

        void BuildCurvatureSpline();
        void BuildCrossesSpline(const std::vector <double>& crosses);
        double get_curvature(const double x) const;
        double getZ(const double x);
        double getZ2(const double x);


        Spline spl;

        //SPLINTER::BSpline curve_spline{1};
        //SPLINTER::BSpline front_spline{1};
        //SPLINTER::DataTable spline_samples;

        double spline_start, spline_end;
        double curve_start, curve_end;

        void Calc2(const std::vector <Particle*>& particle_list);
        void Print(unsigned num);


        std::vector<Point> spline_points;

        struct window_info {
            double sum = 0;
            unsigned count = 0;
        };
        std::vector<Point> window_points;
        std::vector<window_info> window_infos;

        struct flow_point {
            double x, Vx, Vx2;
        };
        std::vector <flow_point> flow_points;


        struct analys_point {
            double div=0, diff2=0, div2=0, cross=0, raw_cross, r=0, c=0, k=0, cx=0, cz=0;
        };
        std::vector <analys_point> analys_points;
        std::vector <double> curvature{0};




        double area_start, area_end, area_size;
        int window_steps=0;
        double window_size, window_radius;
        double window_step_size, window_step_start, window_step_end;


        int spline_steps=0;
        double spline_step_size;



        //double *x{ 0 }, *Vx{ 0 }, * Vx2{ 0 };

        /*
        int steps = 0;
        double window = 0;
        unsigned h;
        double radius;
        double steps_start;
        double steps_end;
        double steps_area;
        double step_size;

        double w_percent;

        double avg = 0, _avg = 0;
        double deviation = 0, error = 0;
        */

    };

}
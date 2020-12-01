#pragma once

#include "Params.h"
#include "Particle.h"
#include <list>
#include <vector>
#include <string>
#include <fmt/core.h>
#include <iostream>
#include <fstream>


namespace ps {

    class Frontline
    {
        const Params* P;
    public:
        Frontline(const Params& P_) : P(&P_) { Init(); }
        void Init(int steps, double window, double
        tart, double area_end);
        void Init();


        void Calc(const std::vector <Particle*>& particle_list);
        void CalcError();
        void CalcRadius();

        void Calc2(const std::vector <Particle*>& particle_list);
        void FivePointStencil();
        void Print(unsigned num);

        struct flow_point {
            double x, Vx, Vx2;
        };
        flow_point* points{ 0 };

        struct front_line_point {
            double x=0, z=0, div=0, sum=0, diff2=0, div2=0, cross=0, r=0;
            unsigned count = 0;
        };

        Point* frontline_coords{0};

        std::vector <front_line_point> front_line_points;

        //double *x{ 0 }, *Vx{ 0 }, * Vx2{ 0 };

        int steps = 0;
        double window = 0;
        unsigned h;
        double area_start, area_end;
        double radius;
        double steps_start;
        double steps_end;
        double steps_area;
        double step_size;

        double w_percent;

        double avg = 0, _avg = 0;
        double deviation = 0, error = 0;

    public:
        ~Frontline() {
            delete[]points;
        }
    };

}
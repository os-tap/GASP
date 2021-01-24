#include "Frontline.h"

#include <SPLINTER/datatable.h>
#include <SPLINTER/bspline.h>
#include <SPLINTER/bsplinebuilder.h>

namespace ps {


    void Frontline::Init() {

        area_start = P->area_beg;
        area_end = P->area_end;
        area_size = area_end - area_start;

        /*   WINDOW MIDDLE PARAMS   */

        if (window_steps != P->frontline_window_steps)
        {
            window_steps = P->frontline_window_steps;
            window_steps_2 = window_steps * window_steps;
            window_points.resize(window_steps_2);
        }

        window_radius = (area_size) / (window_steps + 1.) * P->frontline_window_size;
        window_size = window_radius * 2;

        window_step_size = (area_size - window_size) / (window_steps - 1.);
        for (int i = 0; i < window_steps; ++i) {
            for (int j = 0; j < window_steps; ++j) {
                auto& wp = window_points[j * window_steps + i];
                wp.x = area_start + window_radius + window_step_size * i;
            }
        }
        for (int j = 0; j < window_steps; ++j) {
            for (int i = 0; i < window_steps; ++i) {
                auto& wp = window_points[j * window_steps + i];
                wp.y = area_start + window_radius + window_step_size * j;
            }
        }
        for (int i = 0; i < window_steps_2; i++)
        {
            auto& wp = window_points[i];
            wp.r = sqrt(wp.x * wp.x + wp.y * wp.y);
        }
        window_area_radius = area_size / 2 - window_radius;

        /*   SPLINE PARAMS   */

        if (spline_steps != P->frontline_spline_steps)
        {
            spline_steps = P->frontline_window_steps;
            spline_points.resize(window_steps_2);
            flow_points.resize(window_steps_2);
            analys_points.resize(window_steps_2);
        }

        spline_step_size = (area_size - window_size) / (spline_steps - 1.);

        /*for (int i = 0; i < length; i++)
        {
            double r = 0;
            int r_steps = r * r / 2;
            double deg_step = M_PI * 2 / r_steps;
            double deg = 0;
            for (int i = 0; i < r_steps; i++)
            {
                double x = r * sin(deg);
                double y = r * cos(deg);
                deg += deg_step;
                //spline_points.push_back();
            }
        }*/



        /*for (int i = 0; i < spline_steps; ++i) {
            spline_points[i].x = area_start + window_radius + spline_step_size * i;
            
            flow_points[i].Vx = P->system_speed(spline_points[i].x) * P->burn_speed;
            flow_points[i].Vx2 = flow_points[i].Vx * flow_points[i].Vx;

            flow_points[i].x = spline_points[i].x;

            analys_points[i] = {};
        }*/

        for (int i = 0; i < window_steps_2; ++i) {
            spline_points[i].x = window_points[i].x;
            spline_points[i].y = window_points[i].y;
        }


    }


    void Frontline::Calc(const std::vector <Particle*>& particle_list) {

        WindowMiddle(particle_list);
        //SplineSmooth(P->frontline_spline_alpha);
        //CalcError();
        //FivePointStencil(P->frontline_stencil_h);
        //CalcNormal();
        //CalcRadius(P->frontline_radius_h);
    }

    void Frontline::WindowMiddle(const std::vector <Particle*>& particle_list) {
        for (auto& sp : spline_points) {
            sp.z = 0;
        }

        for (auto& wp : window_points) {
            wp.z = wp.sum = wp.count = 0;
        }

        for (const auto& particle : particle_list) {

            int beg_x = (int)ceil((particle->x - area_start - window_size) / window_step_size);
            beg_x *= beg_x > 0;
            int end_x = (int)floor((particle->x - area_start) / window_step_size);
            if (end_x > window_steps - 1) end_x = window_steps - 1;
            //end_x = end_x * (end_x < window_steps) + (end_x > window_steps - 1) * (window_steps - 1);

            int beg_y = (int)ceil((particle->y - area_start - window_size) / window_step_size);
            beg_y *= beg_y > 0;
            int end_y = (int)floor((particle->y - area_start) / window_step_size);
            if (end_y > window_steps - 1) end_y = window_steps - 1;
            //end_y = end_y * (end_y < window_steps) + (end_y > window_steps - 1) * (window_steps - 1);

            for (int i = beg_x; i <= end_x; ++i) {
                for (int j = beg_y; j <= end_y; ++j) {
                    auto &wp = window_points[j * window_steps + i];
                    wp.sum += particle->z;
                    wp.count++;
                }
            }

        }
        avg = 0;
        for (int i = 0; i < window_steps_2; ++i) {
            if (window_points[i].count) {
                window_points[i].z = spline_points[i].z = window_points[i].sum / window_points[i].count;
                //window_points[i].z -= P->system_speed(window_points[i].x) / 2;
                avg += window_points[i].z;
            }
        }
        avg /= window_steps_2;


    }

    void Frontline::SplineSmooth(double alpha) {
        SPLINTER::DataTable samples;
        SPLINTER::DenseVector xy(2);

        for (const auto& wp : window_points) {
            //if (wp.z && wp.r < window_area_radius) 
            {
                xy(0) = wp.x;
                xy(1) = wp.y;
                samples.addSample(xy, wp.z);
            }
        }
        
        if (samples.cbegin() == samples.cend()) {
            for (auto& sp : spline_points) {
                sp.z = 0;
            }
            return;
        };

        SPLINTER::BSpline pspline = SPLINTER::BSpline::Builder(samples)
            .degree(3)
            .smoothing(SPLINTER::BSpline::Smoothing::PSPLINE)
            .alpha(P->frontline_spline_alpha)
            .build();


        //SPLINTER::DenseVector xd(1);

        /*for (size_t i = 0; i < spline_steps; i++)
        {
            spline_points[i].x = flow_points[i].x;
            xy(0) = spline_points[i].x;
            spline_points[i].z = pspline.eval(xy);
        }*/
        for (auto& sp : spline_points) {
            xy(0) = sp.x;
            xy(1) = sp.y;
            sp.z = pspline.eval(xy);
        }

    }


    void Frontline::CalcNormal()
    {
        double gap = P->burn_radius / 2;
        for (size_t i = 0; i < spline_steps; i++)
        {
            if (spline_points[i].z && analys_points[i].div)
            {
                double d = ((int)(analys_points[i].div > 0) * 2 - 1);
                double k = 1. / analys_points[i].div;
                double dx = d * gap / sqrt(k * k + 1);
                double dy = dx * k;
                spline_points[i].x += dx;
                spline_points[i].z -= dy;
            }
            else {
                spline_points[i].z = 0;
            }
        }
    }

    void Frontline::CalcError()
    {
        int active_window_points = 0;
        deviation = 0;
        double d;
        for (int i = 0; i < window_steps_2; ++i) {
            if (window_points[i].count) {
                d = window_points[i].z - avg;
                deviation += d * d;
                ++active_window_points;
            }
        }
        deviation /= (double)(active_window_points - 1);
        error = sqrt(deviation) / P->stream_width;
    }




    //void Frontline::Calc2(const std::vector <Particle*>& particle_list) {

    //    for (int i = 0; i < steps; ++i) {
    //        //Vx[i] = P->system_speed(frontline_window_points[i].x) * P->burn_speed;
    //        //Vx2[i] = Vx[i] * Vx[i];
    //        frontline_window_points[i].sum = frontline_window_points[i].count = 0;
    //        frontline_window_points[i].z = 0;
    //    }

    //    for (const auto& particle : particle_list) {

    //        int w = (int)floor((particle->_x() - area_start) * w_percent);
    //        if (particle->_z() < frontline_window_points[w].z)
    //        {
    //            frontline_window_points[w].count = 1;
    //            frontline_window_points[w].z = particle->_z();
    //            //frontline_window_points[w].x = particle->_x();
    //        }

    //    }

    //    /*for (int i = 0; i < steps; ++i) {
    //        frontline_window_points[i].z = frontline_window_points[i].sum / frontline_window_points[i].count;
    //    }*/

    //    FivePointStencil(P->frontline_stencil_h);
    //}

    void Frontline::Print(unsigned num) {
        std::ofstream csv(P->csv_folder + "line.csv." + std::to_string(num));

        //std::string output = P->frontline_params();
        //std::string output = "x,z,wx,wz,div,div2,Vx2,diff2,cross,r";
        std::string output = "x,y,z,sz";

        for (int i = 0; i < window_steps_2; ++i) {
            //if (window_points[i].r < window_area_radius) 
            {
                output += fmt::format("\n{},{},{},{}",
                    window_points[i].x,
                    window_points[i].y,
                    window_points[i].z,
                    spline_points[i].z
                );
            }

                /*output += fmt::format("\n{},{},{},{},{},{},{},{},{},{}",
                    spline_points[i].x,
                    spline_points[i].z,
                    window_points[i].x,
                    window_points[i].z,
                    analys_points[i].div,
                    analys_points[i].div2,
                    flow_points[i].Vx2,
                    analys_points[i].diff2,
                    analys_points[i].cross,
                    analys_points[i].r
                );*/
                //output+= '\n' + std::to_string(frontline_window_points[i].x) + ',' + std::to_string(frontline_window_points[i].y) + ',' + std::to_string(frontline_window_points[i].div);
            
        }
        csv << output;
        csv.close();
    }

    void Frontline::FivePointStencil(int h_div) {
        assert(h_div >= 1);


        for (int i = h_div*2; i < spline_steps - h_div * 2; ++i) {


            analys_points[i].div =
                8 * (spline_points[i + h_div].z - spline_points[i - h_div].z)
                - spline_points[i + h_div * 2].z
                + spline_points[i - h_div * 2].z;

            analys_points[i].div /= h_div * spline_step_size * 12;

            analys_points[i].div2 = analys_points[i].div * analys_points[i].div + 1;

            analys_points[i].diff2 =
                - spline_points[i + h_div * 2].z
                + 16 * spline_points[i + h_div].z
                - 30 * spline_points[i].z
                + 16 * spline_points[i - h_div].z
                - spline_points[i - h_div * 2].z;

            analys_points[i].diff2 /= 12 * h_div * h_div * spline_step_size * spline_step_size;



        }

        //exit(0);
    }


    void Frontline::CalcRadius(int h_div) {
        
        for (int i = h_div; i < spline_steps - h_div; ++i) {
            const auto& A = spline_points[i];
            const auto& B = spline_points[i - h_div];
            const auto& C = spline_points[i + h_div];

            Point M[2], Cntr{0}; double H1, H2, G;

            M[0].x = B.x - A.x;
            M[0].z = B.z - A.z;

            M[1].x = C.x - A.x;
            M[1].z = C.z - A.z;

            H1 = M[0].x * (A.x + B.x) + M[0].z * (A.z + B.z);
            H2 = M[1].x * (A.x + C.x) + M[1].z * (A.z + C.z);

            G = (M[0].x * (C.z - B.z) - M[0].z * (C.x - B.x)) * 2;

            Cntr.x = (M[1].z * H1 - M[0].z * H2) / G;
            Cntr.z = (M[0].x * H2 - M[1].x * H1) / G;

            double dx = (A.x - Cntr.x);
            double dy = (A.z - Cntr.z);

            analys_points[i].r = 1. / sqrt(dx*dx + dy*dy);
        }


    }


}


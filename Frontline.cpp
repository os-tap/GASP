#include "Frontline.h"


namespace ps {


    void Frontline::Init() {

        area_start = P->area_beg;
        area_end = P->area_end;
        area_size = area_end - area_start;

        /*   WINDOW MIDDLE PARAMS   */

        if (window_steps != P->frontline_window_steps)
        {
            window_steps = P->frontline_window_steps;
            window_points.resize(window_steps);
            flow_points.resize(window_steps);
        }

        window_radius = (area_size) / (window_steps + 1.) * P->frontline_window_size;
        window_size = window_radius * 2;

        first_point = area_start + window_radius;
        last_point = area_end - window_radius;

        window_step_size = (area_size - window_size) / (window_steps - 1.);
        for (int i = 0; i < window_steps; ++i) {
            window_points[i].x = flow_points[i].x = area_start + window_radius + window_step_size * i;

            flow_points[i].Vx = P->system_speed(flow_points[i].x) * P->burn_speed;
            flow_points[i].Vx2 = flow_points[i].Vx * flow_points[i].Vx;
        }


        /*   SPLINE PARAMS   */

        if (spline_steps != P->frontline_spline_steps)
        {
            spline_steps = P->frontline_spline_steps;
            spline_points.resize(spline_steps);
            analys_points.resize(spline_steps);
            curvature.resize(spline_steps);
        }

        spline_step_size = (area_size - window_size) / (spline_steps - 1.);

        for (int i = 0; i < spline_steps; ++i) {
            spline_points[i].x = area_start + window_radius + spline_step_size * i;

            //analys_points[i] = {};
        }


    }


    void Frontline::Calc(const std::vector <Particle*>& particle_list) {

        /*for (int i = 0; i < spline_steps; ++i) {
            analys_points[i] = {};
        }*/
        for (auto& ap : analys_points) ap = {};
        for (double& c : curvature) c = 0;

        WindowMiddle(particle_list);
        SplineSmooth(P->frontline_spline_alpha);
        //CalcError();
        FivePointStencil(P->frontline_stencil_h);

        if(P->move_normal) CalcNormal();

        CalcRadius(P->frontline_radius_h);
        CalcCurve();

        for (size_t i = 0; i < spline_steps; i++)
        {
            curvature[i] = P->curve_by_radius ? analys_points[i].r : analys_points[i].c;
        }
    }

    void Frontline::WindowMiddle(const std::vector <Particle*>& particle_list) {

        for (auto& wp : window_points) {
            wp.z = wp.sum = wp.count = 0;
        }

        for (const auto& particle : particle_list) {

            int beg_i = (int) ceil((particle->x - area_start - window_size ) / window_step_size);
            int end_i =	(int)floor((particle->x - area_start) / window_step_size);
            beg_i *= beg_i > 0;
            end_i = end_i * (end_i < window_steps) + (end_i > window_steps - 1) * (window_steps - 1);
            //if (end_i > window_steps - 1) end_i = window_steps - 1;

            for (int i = beg_i; i <= end_i; ++i) {
                window_points[i].sum += particle->z;
                window_points[i].count++;
            }

        }

        for (int i = 0; i < window_steps; ++i) {
            if (window_points[i].count) {
                window_points[i].z += window_points[i].sum / window_points[i].count;
            }
        }

        if (P->move_speed)
        {
            for (int i = 0; i < window_steps; ++i) {
                if (window_points[i].count) {
                    window_points[i].z -= P->system_speed(window_points[i].x) / 2;
                }
            }
        }


    }

    void Frontline::SplineSmooth(double alpha) {
        SPLINTER::DataTable samples;

        double spline_start = area_start;
        double spline_end = area_end;

        for (const auto& wp : window_points) {
            if (wp.z) samples.addSample(wp.x, wp.z);
        }

        for (size_t i = window_points.size()-2; i ; i--) {
            if (window_points[i].z) {
                spline_end = window_points[i].x;
                break;
            }
        }
        for (size_t i = 1; i < window_points.size(); i++) {
            if (window_points[i].z) {
                spline_start = window_points[i].x;
                break;
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


        SPLINTER::DenseVector xd(1);
        for (size_t i = 0; i < spline_steps; i++)
        {
            spline_points[i].x = flow_points[i].x;
            if (spline_points[i].x >= spline_start && spline_points[i].x <= spline_end)
            {
                xd(0) = spline_points[i].x;
                spline_points[i].z = pspline.eval(xd);
            }
            
        }
        /*for (auto& sp : spline_points) {
            xd(0) = sp.x;
            sp.z = pspline.eval(xd);
        }*/

    }


    void Frontline::CalcNormal()
    {
        //double gap = P->burn_radius / 2;
        for (size_t i = 0; i < spline_steps; i++)
        {
            if (spline_points[i].z && analys_points[i].div)
            {
                double gap = P->get_burn_radius(spline_points[i].x) / 2;
                double d = ((int)(analys_points[i].div > 0) * 2 - 1);
                double k = 1. / analys_points[i].div;
                analys_points[i].k = k;
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
        /*int active_window_points = 0;
        deviation = 0;
        double d;
        for (int i = 0; i < steps; ++i) {
            if (frontline_window_points[i].count) {
                d = frontline_window_points[i].z - avg;
                deviation += d * d;
                ++active_window_points;
            }
        }
        deviation /= (double)(active_window_points - 1);
        error = sqrt(deviation) / P->stream_width;*/
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
        std::string output = "x,z,wx,wz,div,div2,Vx2,diff2,cross,raw_cross,r,c,k";

        for (int i = 0; i < spline_steps; ++i) {

            if (spline_points[i].z)
            {

                output += fmt::format("\n{},{},{},{},{},{},{},{},{},{},{},{},{}",
                    spline_points[i].x,
                    spline_points[i].z,
                    window_points[i].x,
                    window_points[i].z,
                    analys_points[i].div,
                    analys_points[i].div2,
                    flow_points[i].Vx2,
                    analys_points[i].diff2,
                    analys_points[i].cross,
                    analys_points[i].raw_cross,
                    analys_points[i].r,
                    analys_points[i].c,
                    analys_points[i].k
                );
                //output+= '\n' + std::to_string(frontline_window_points[i].x) + ',' + std::to_string(frontline_window_points[i].y) + ',' + std::to_string(frontline_window_points[i].div);
            }
        }
        csv << output;
        csv.close();
    }
    
    void Frontline::FivePointStencil(int h_div) {
        assert(h_div >= 1);

        int start = 0;
        int end = spline_steps;
        while (start < end && spline_points[start].z <= 0) ++start;
        while (end > start && spline_points[end-1].z <= 0) --end;


        for (int i = start + h_div*2; i < end - h_div * 2; ++i) {


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

        int start = 0;
        int end = spline_steps;
        while (start < end && spline_points[start].z <= 0) ++start;
        while (end > start && spline_points[end - 1].z <= 0) --end;
        
        for (int i = start + h_div + P->frontline_stencil_h; i < end - h_div - P->frontline_stencil_h; ++i) {
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
    void Frontline::CalcCurve() {
        int start = 0;
        int end = spline_steps;
        while (start < end && spline_points[start].z <= 0) ++start;
        start += P->frontline_stencil_h * 2;
        while (end > start && spline_points[end - 1].z <= 0) --end;
        end -= P->frontline_stencil_h*2;

        for (int i = start; i < end; ++i) {
            double ddy = (1 + analys_points[i].div * analys_points[i].div) / analys_points[i].diff2;
            double x1 = spline_points[i].x - analys_points[i].div * ddy;
            double z1 = spline_points[i].z + ddy;
            analys_points[i].c = 1./ sqrt
                ((x1 - spline_points[i].x) * (x1 - spline_points[i].x) + (z1 - spline_points[i].z) * (z1 - spline_points[i].z));
        }

        
        bool scale_area = false;
        double scale_start;
        kinks.clear();
        for (int i = 0; i < spline_points.size(); i++)
        {
            if (!scale_area && analys_points[i].c >= P->scale_burn_condition) {
                scale_area = true;
                scale_start = spline_points[i].x;
            }
            else if (scale_area && analys_points[i].c < P->scale_burn_condition) {
                scale_area = false;
                kinks.push_back(scale_start + (spline_points[i - 1].x - scale_start)/2);
                //std::cout << "scale\n";
            }
        }

    }

    void Frontline::SetCrosses(const std::vector<double>& crosses)
    {
        if (crosses.size() != analys_points.size())
        {
            std::cout << "\n\ncrosses.size() != analys_points.size()";
            exit(22);
        }
        for (size_t i = 0; i < spline_steps; i++)
        {
            analys_points[i].cross = crosses[i];
        }
    }

    void Frontline::BuildCurvatureSpline(const std::vector<double>& crosses)
    {
        if (crosses.size() != analys_points.size())
        {
            std::cerr << "\n\ncrosses.size() != analys_points.size()";
            exit(22);
        }

        for (size_t i = 0; i < spline_steps; i++)
        {
            analys_points[i].raw_cross = crosses[i];
        }

        SPLINTER::DataTable samples;

        for (size_t i = 0; i < spline_steps; i++) {
            curve_start = spline_points[i].x;
            if (crosses[i]) break;
        }

        for (size_t i = spline_steps - 1; i>=0; i--) {
            curve_end = spline_points[i].x;
            if (crosses[i]) break;
        }

        //std::cout << '\n' << curve_start;



        for (int i = 0; i < spline_steps; ++i) {
            if (crosses[i]) samples.addSample(spline_points[i].x, P->frontline_cross_border - crosses[i]);
        }

        if (samples.cbegin() == samples.cend()) return;

    
        curve_spline = SPLINTER::BSpline::Builder(samples)
            .degree(3)
            .smoothing(SPLINTER::BSpline::Smoothing::PSPLINE)
            .alpha(P->frontline_cross_spline_alpha)
            .build();



        for (size_t i = 0; i < spline_steps; i++)
        {
            analys_points[i].cross = get_curvature(spline_points[i].x);
        }
    }

    double Frontline::get_curvature(const double x) const
    {
        if (x < curve_start || x > curve_end) return 0;

        SPLINTER::DenseVector xd(1);
        xd(0) = x;
        double c = curve_spline.eval(xd);

        return c * (c > 0);

    }


}


#pragma once

#include <vector>

#include <SPLINTER/datatable.h>
#include <SPLINTER/bspline.h>
#include <SPLINTER/bsplinebuilder.h>

#include "Point.h"

namespace ps {
    class PSpline;
}
class ps::PSpline
{

    SPLINTER::BSpline spline{1};
    SPLINTER::DataTable spline_samples;

    double spline_start, spline_end;

    std::unordered_map<double, double> cache;

public:
    PSpline() {}

    void Fit(const std::vector<Point>& points, double _alpha) {

        spline_samples = SPLINTER::DataTable{};


        for (const auto& wp : points) {
            if (wp.z) spline_samples.addSample(wp.x, wp.z);
        }

        if (spline_samples.cbegin() == spline_samples.cend()) {
            return;
        }


        spline_start = points.front().x;
        spline_end = points.back().x;
        for (size_t i = 0; i < points.size(); i++) {
            spline_start = points[i].x;
            if (points[i].z) {
                break;
            }
        }
        for (size_t i = points.size()-1; i>0 ; i--) {
            spline_end = points[i].x;
            if (points[i].z) {
                break;
            }
        }
    

        spline = SPLINTER::BSpline::Builder(spline_samples)
            .degree(3)
            .smoothing(SPLINTER::BSpline::Smoothing::PSPLINE)
            .alpha(_alpha)
            .build();


    }
    double eval(const double x) const
    {
        if (x < spline_start || x > spline_end) return 0;

        SPLINTER::DenseVector xd(1);
        xd(0) = x;
        return spline.eval(xd);

    }

    ~PSpline() { }
};


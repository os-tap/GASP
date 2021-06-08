#pragma once
#pragma once

#include <vector>
#include <unordered_map>

#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>

#include "Point.h"

namespace ps {
    class BSpline;
}
class ps::BSpline
{

    size_t n = 0;
    size_t k = 4;
    size_t ncoeffs = 12;
    size_t nbreak = ncoeffs - 2;

    /* use uniform breakpoints on [0, 15] */
    /* allocate a cubic bspline workspace (k = 4) */
    gsl_bspline_workspace* bw;
    gsl_vector* B;
    gsl_matrix* X;
    gsl_vector* y, * w;

    gsl_multifit_linear_workspace* mw;
    gsl_vector* c;
    gsl_matrix* cov;

    std::unordered_map<double, double> cache;

public:
    BSpline() {
        Init(4,12);
    }

    void Init(size_t k, size_t coef) {
        ncoeffs = coef;
        nbreak = ncoeffs - 2;

        gsl_bspline_free(bw);
        gsl_vector_free(B);
        gsl_vector_free(c);
        gsl_matrix_free(cov);

        bw = gsl_bspline_alloc(k, nbreak);
        B = gsl_vector_alloc(ncoeffs);

        c = gsl_vector_alloc(ncoeffs);
        cov = gsl_matrix_alloc(ncoeffs, ncoeffs);


        gsl_matrix_free(X);
        gsl_multifit_linear_free(mw);
        X = gsl_matrix_alloc(n, ncoeffs);
        mw = gsl_multifit_linear_alloc(n, ncoeffs);
    }

    void Fit(const std::vector<Point>& points, double start, double end) {

        gsl_bspline_knots_uniform(start, end, bw);

        if (n != points.size()) {

            gsl_matrix_free(X);
            gsl_vector_free(y);
            gsl_vector_free(w);
            gsl_multifit_linear_free(mw);

            n = points.size();
            X = gsl_matrix_alloc(n, ncoeffs);
            mw = gsl_multifit_linear_alloc(n, ncoeffs);

            y = gsl_vector_alloc(n);
            w = gsl_vector_alloc(n);
        }

        for (int i = 0; i < n; ++i)
        {
            gsl_vector_set(y, i, points[i].z);
            double sigma = 0.1 * points[i].z;
            gsl_vector_set(w, i, 1.0 / (sigma * sigma));

            /* compute B_j(xi) for all j */
            gsl_bspline_eval(points[i].x, B, bw);

            /* fill in row i of X */
            for (int j = 0; j < ncoeffs; ++j)
            {
                double Bj = gsl_vector_get(B, j);
                gsl_matrix_set(X, i, j, Bj);
            }
        }

        double chisq;
        gsl_multifit_wlinear(X, w, y, c, cov, &chisq, mw);


    }
    double eval(const double xi) {
        double yi, yerr;
        gsl_bspline_eval(xi, B, bw);
        gsl_multifit_linear_est(B, c, cov, &yi, &yerr);
        return yi;
    }

    ~BSpline() {
        gsl_bspline_free(bw);
        gsl_vector_free(B);
        gsl_vector_free(y);
        gsl_matrix_free(X);
        gsl_vector_free(c);
        gsl_vector_free(w);
        gsl_matrix_free(cov);
        gsl_multifit_linear_free(mw);
    }
};


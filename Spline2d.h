//#pragma once
//#include <gsl/gsl_math.h>
//#include <gsl/gsl_interp2d.h>
//#include <gsl/gsl_spline2d.h>
//
//#include <iostream>;
//
//namespace ps {
//
//class Spline2d
//{
//
//    gsl_spline2d* spline;
//    size_t nx, ny;
//    double* xa, * ya, * za;
//    gsl_interp_accel* xacc;
//    gsl_interp_accel* yacc;
//
//public:
//    void init(double *_xa, double *_ya, size_t _nx, size_t _ny) {
//
//        nx = _nx;
//        ny = _ny;
//        xa = _xa;
//        ya = _ya;
//        za = new double[nx * ny];
//
//        spline = gsl_spline2d_alloc(gsl_interp2d_bilinear, nx, ny);
//        xacc = gsl_interp_accel_alloc();
//        yacc = gsl_interp_accel_alloc();
//
//    }
//
//    int set(size_t xi, size_t yi, double z) {
//        return gsl_spline2d_set(spline, za, xi, yi, z);
//    }
//    void fit() {
//        gsl_spline2d_init(spline, xa, ya, za, nx, ny);
//        //std::cout << "\nxmin: " << spline->interp_object.xmin;
//        //std::cout << "\nxmax: " << spline->interp_object.xmax;
//    }
//    double eval(double x, double y) {
//        if (x < spline->interp_object.xmin || x > spline->interp_object.xmax || y > spline->interp_object.ymax || y < spline->interp_object.ymin) return 0;
//        return gsl_spline2d_eval(spline, x, y, xacc, yacc);
//    }
//
//    ~Spline2d()
//    {
//        delete[] xa,ya,za;
//        gsl_spline2d_free(spline);
//        gsl_interp_accel_free(xacc);
//        gsl_interp_accel_free(yacc);
//    }
//
//};
//
//}

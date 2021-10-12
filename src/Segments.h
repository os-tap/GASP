#pragma once
#include "Params.h"
#include "Particle.h"

#include "Spline2d.h"

#include <vector>
#include <string>

#include <iostream>
#include <fstream>
#include <fmt/core.h>

#include <random>
#include <cmath>
#include <cassert>

#include <cstdio>
#include <utility>
#include <limits>

//#undef min
//#undef max
#define NOMINMAX 
//#include <tbb/task_scheduler_init.h>
#include <tbb/enumerable_thread_specific.h>
#undef NOMINMAX 
//#include <tbb/parallel_for.h>
//#include <tbb/blocked_range.h>
#include <tbb/concurrent_vector.h>
#include <execution>
#include <algorithm>


namespace ps {


    typedef std::vector <ps::Particle> ParticleList;
    typedef tbb::enumerable_thread_specific<int> CounterType;
    typedef tbb::enumerable_thread_specific<std::pair<std::uniform_real_distribution<double>, std::random_device>> GenType;

    class Segments {

        const Params* P;

        CounterType SageCounter{0};

        GenType ParticleGenerator;

        


    public:

        Segments(const Params&);
        void Update();

        ParticleList all_list;
        std::vector <Particle*> all_will_burn;

    private:
        tbb::concurrent_vector <size_t> will_burn_index;
        tbb::concurrent_vector <Particle> all_will_burn_concurrent;
        tbb::concurrent_vector <Particle> refilled;



    public:

        void Emit();
        void Refill();
        void FillParticles();
        void MoveParticles();
        void CrossParticles();
        void StepParticles();
        void ClearParticles();

        void CalcBurnRadius(int g = 1);
        void SplineGrid();
        void PlaceBurned();

        void LightsOut();
        void EraseParticles();

        void PrintGrid(int);
        void PrintCurvature(int);
        void PrintSwarm(int);
        void PrintCount(int n, int count);
        void PrintSVM(int n, int count);


    private:

        void CreateParticle(double x_cord, double z_cord);
        void MoveParticle(Particle& p);
        void StepParticle(Particle& p);
        void BurnParticle(Particle& p);
        void ParticleToSegment(Particle& p, size_t i);


    public:
        void Toggle_Fill();
    private:
        bool _fill_one = 0;

        void Fill_Sampling();

        void Fill_Grid();   
        void ResetFillGrid();
        std::vector <double> last_particles;
        int fill_grid_count;


        int delete_sage = 0;
        bool refill;





    public:

        int grid_count_x, grid_count_z, grid_count; 
        double grid_x_size, grid_z_size;
        double grid_min_size, grid_max_size;
        double grid_count_x_percent, grid_count_z_percent;
        int grid_particles_count;
        int grid_particles_min, grid_particles_max;
        double* gxa{nullptr}, * gya{ nullptr };
        Spline2d spline2d;



        struct SegPointBurn : public Point {
            coord_t r2;
            explicit SegPointBurn(Point& p, coord_t _r2) : Point(p.x, p.z), r2(_r2) {}
        };

        struct SegPointOk : public Point {
            size_t index;
            explicit SegPointOk(Point& p, size_t i) : Point(p.x, p.z), index(i) {}
            bool Cross(const SegPointBurn& p) const { return Point::Cross(p, p.r2); }
        };


        struct SegPoint {
            SegPoint(Point &p, double _r2, size_t i) : x(p.x), z(p.z), r2(_r2), index(i) {}
            SegPoint(Point &p, size_t i) : x(p.x), z(p.z), r2(0), index(i) {}
            bool Cross(SegPoint p) { return ((x - p.x) * (x - p.x) + (z - p.z) * (z - p.z)) <= p.r2; }
            bool CrossOk(SegPoint p) { return ((x - p.x) * (x - p.x) + (z - p.z) * (z - p.z)) <= r2; }
            float x, z, r2;
            size_t index;
        };
        struct FrontSegPoint : SegPoint {
            FrontSegPoint(Point& p, double _r2, size_t i) : SegPoint(p, _r2, i) {}
            unsigned count = 0;
        };

        struct Segment {
            int x, z;
            double br, br2;
            double curvature = 0, seg_start_x, seg_start_z;;
            int c_ok = 0, c_b = 0;
            int ok_size = 0, b_size = 0;
            tbb::concurrent_vector <SegPointOk> ok_list, b_list;
            tbb::concurrent_vector <SegPointBurn> burn_list;
            std::vector<int> burn_indexes;

            std::vector<FrontSegPoint> front_points;
            //std::vector<std::pair<int, int>> front_points;
        };

        std::vector<Segment> grid; 
        //tbb::concurrent_vector <Segment*> burn_segments;

        inline Segment& grids(size_t x, size_t z) {
            assert(z * grid_count_x + x < grid.size());
            return grid[z * grid_count_x + x];
        }

        int GetSegmentX(double x) const;
        int GetSegmentZ(double z) const;
        Segment& SegmentByPoint(double x, double z);
        void BurnSegmentByPoint(double x, double z);


        void SetSegmentsGrid(double);
        void UpdateSegments();
        void BurnSegment(Segment&);
        void DoSegment(Segment&); 
        void ClearSegments();



        std::vector<int> circle_squares;
        void CircleSquareSampling(double r);




        /**/
    public:

        void Density_Grid();

        /**/

    };


}; /* namespace ps */


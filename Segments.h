#pragma once
#include "Params.h"
#include "Particle.h"

#include <vector>
#include <string>

#include <iostream>
#include <fstream>
#include <fmt/core.h>

#include <random>
#include <cmath>
#include <cassert>


#define NOMINMAX 
#include <tbb/enumerable_thread_specific.h>

#include <tbb/concurrent_vector.h>
#include <pstl/execution>
#include <pstl/algorithm>

namespace ps {

    typedef std::vector <ps::Particle> ParticleList;
    typedef tbb::enumerable_thread_specific<std::pair<std::uniform_real_distribution<double>, std::random_device>> GenType;

    class Segments {

        const Params* P;
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

        void FillParticles();
        void Emit();
        void Refill();
        //void PlaceBurned();
        void MoveParticles();
        void CrossParticles();
        void StepParticles();
        void ClearParticles();

        void LightsOut();
        void EraseParticles();

        void PrintSwarm(int);
        void PrintCount(int n, int count);
        void PrintSVM(int n, int count);


    private:


        void CreateParticle(double x_cord, double z_cord, double p_speed);
        void CreateParticle(double x_cord, double z_cord);
        void MoveParticle(Particle& p);
        void StepParticle(Particle& p);
        void BurnParticle(Particle& p);
        void BurnParticle(Particle* p);
        void ParticleToSegment(Particle& p);
        void ParticleToSegment(Particle& p, size_t i);
        void ParticleToSegment(size_t i);
        void FrontPointToSegment(Point& p, size_t index);

        bool ParticleInBurnSegments(Particle* particle, int seg_x, int seg_z);


    public:
        void Toggle_Fill();
    private:
        bool _fill_one = 0;

        void Fill_Sampling();

        void Fill_Grid();   
        void ResetFillGrid();
        std::vector <double> last_particles;
        int fill_grid_count;





    public:
        bool refill;
        int grid_count_x, grid_count_z, grid_count; 
        double grid_x_size, grid_z_size;
        double grid_min_size, grid_max_size;
        double grid_count_x_percent, grid_count_z_percent;

        struct SegPoint {
            SegPoint(Point &p, double _r2, size_t i) : x(p.x), z(p.z), r2(_r2), index(i) {}
            SegPoint(Particle &p, size_t i) : x(p.x), z(p.z), r2(p.burn_radius_2), index(i) {}
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
            tbb::concurrent_vector <SegPoint> ok_list, burn_list, b_list;
            std::vector<int> burn_indexes;

            std::vector<FrontSegPoint> front_points;
            double seg_start_x, seg_start_z;// seg_end_x, , seg_end_z;
            //std::vector<std::pair<int, int>> front_points;
        };

        std::vector<Segment> grid; 
        //tbb::concurrent_vector <Segment*> burn_segments;

        Segment& grids(size_t x, size_t z) {
            assert(z * grid_count_x + x < grid.size(), "segments out of range");
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


        bool CheckSegmentNeighborsBurn(int seg_x, int seg_z);




        /**/
    public:
        int Line_Count();
        void Density_Grid();
        void Density_Radius();
        void Max_Radius();

        void CalcFrontlineRadius(std::vector <Point>& points);
        std::vector<double> front_crosses;
        void CalcFrontlineRadius2(std::vector <Point>& points);
        void RefractParticles();
        /**/

    };


} /* namespace ps */


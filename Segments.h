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

#include <tbb/concurrent_vector.h>
#include <pstl/execution>
#include <pstl/algorithm>

typedef std::vector <ps::Particle> ParticleList;

namespace ps {

    class Segments {

        const Params* P;

    public:

        Segments(const Params&);
        void Update();

        ParticleList all_list;
        std::vector <Particle*> all_will_burn;

    private:
        tbb::concurrent_vector <size_t> will_burn_index;
        tbb::concurrent_vector <Particle> all_will_burn_concurrent;



    public:

        void FillParticles();
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

        void CreateParticle(double x_cord, double y_cord, double z_cord, double p_speed);
        void MoveParticle(Particle& p);
        void StepParticle(Particle& p);
        void BurnParticle(Particle& p);
        void BurnParticle(Particle* p);
        void ParticleToSegment(Particle& p);
        void ParticleToSegment(Particle& p, size_t i);

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

        int grid_count_x, grid_count_y, grid_count_z, grid_count; 
        double grid_x_size, grid_y_size, grid_z_size;
        double grid_min_size, grid_max_size;
        double grid_count_x_percent, grid_count_y_percent, grid_count_z_percent;

        struct SegPoint {
            float x, y, z, r2;
            SegPoint(Particle &p, size_t i) : x(p.x), y(p.y), z(p.z), r2(p.burn_radius_2), index(i) { }
            bool Cross(SegPoint p) { 
                return ((x - p.x) * (x - p.x) + (y - p.y) * (y - p.y) + (z - p.z) * (z - p.z)) <= p.r2; 
            }
            size_t index;
        };

        struct Segment {
            //int x, y, z;
            tbb::concurrent_vector <SegPoint> ok_list, burn_list;
            std::vector<int> burn_indexes;
        };

        std::vector<Segment> grid; 
        //tbb::concurrent_vector <Segment*> burn_segments;

        Segment& grids(size_t x, size_t y, size_t z) {
            assert(z * grid_count_x * grid_count_y + y * grid_count_x + x < grid_count);
            return grid[z * grid_count_x * grid_count_y + y * grid_count_x + x];
        }

        int GetSegmentX(double x) const;
        int GetSegmentY(double y) const;
        int GetSegmentZ(double z) const;
        Segment& SegmentByPoint(double x, double y, double z);
        void BurnSegmentByPoint(double x, double z);


        void SetSegmentsGrid(double);
        void UpdateSegments();
        void BurnSegment(Segment&);
        void DoSegment(Segment&); 
        void ClearSegments();


        //bool CheckSegmentNeighborsBurn(int seg_x, int seg_z);




        /**/
    public:
        int Line_Count();
        void Density_Grid();
        void Density_Radius();
        void Max_Radius();

        void CalcFrontlineRadius(std::vector <Point>& points);
        void RefractParticles();
        /**/

    };


} /* namespace ps */


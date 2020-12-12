#pragma once
#include "Params.h"

#include "Particle.h"
#include "Frontline.h"

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
        tbb::concurrent_vector <Particle*> all_will_burn_concurrent;



    public:

        void FillParticles();
        void MoveParticles();
        void CrossParticles();
        void StepParticles();
        void ClearParticles();

        void LightsOut();
        void EraseParticles();

        void PrintSwarm(int);


    private:

        void CreateParticle(double x_cord, double z_cord, double p_speed);
        void MoveParticle(Particle& p);
        void StepParticle(Particle& p);
        void BurnParticle(Particle* p);
        void ParticleToSegment(Particle& p);

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

        int grid_count_x, grid_count_z, grid_count; 
        double grid_x_size, grid_z_size;
        double grid_min_size, grid_max_size;
        double grid_count_x_percent, grid_count_z_percent;

        struct Segment {
            int x, z;
            tbb::concurrent_vector <Particle*> ok_list, burn_list;
        };

        std::vector<Segment> grid; 
        tbb::concurrent_vector <Segment*> burn_segments;

        Segment& grids(size_t x, size_t z) {
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

        void CalcFrontlineRadius(std::vector <Frontline::front_line_point>& points);
        void RefractParticles();
        /**/

    };


} /* namespace ps */


#pragma once
#include "Params.h"

#include "Particle.h"
#include "Frontline.h"

//#include <list>
#include <vector>
#include <forward_list>

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
//#include <execution>
//#include <algorithm>
//#include <thread>
//#include <mutex>

typedef std::vector <ps::Particle> ParticleList;

namespace ps {

    class Segments {

        const Params* P;

    public:

        Segments(const Params&);
        void Update();


        ParticleList all_list;
        std::vector <Particle*> all_will_burn;
        tbb::concurrent_vector <Particle*> all_will_burn_concurrent;




    public:

        void Fill();
        void Toggle_Fill();

        void DoSegments(int beg, int end);
        void DoSegment(int i);
        void CrossParticles();
        void StepParticles();
        void MoveParticles();
        void ClearParticles();
        void ClearParticleList();

        void LightsOut();
        void Clear();

        void FinalLoop();


        void Print(int);


        int Line_Count();
        void Density_Grid();
        void Density_Radius();
        void Max_Radius();

        void CalcFrontlineRadius(std::vector <Frontline::front_line_point> &points);
        void RefractParticles();





        void Fill_Sampling();
        void Fill_Grid();
    private:

        bool _fill_one = 0;
        //void (Segments::* fill_func)() = 0;


        void ResetFillGrid();
        double* last_particles = 0;
        int fill_grid_count;

    private:



        void CreateParticle(double x_cord, double z_cord, double p_speed);
        void MoveParticle(Particle *p);
        void StepParticle(Particle *p);
        void BurnParticle(Particle *p);
        void ParticleToSegment(Particle *p);




    public:


        /*struct segment_mutex {
            std::mutex ok, burn;
        };*/

        struct Segment {
            bool has_burn = 0, refract = 0;
            int x,z;
//            std::vector<Particle*> ok_list, burn_list, will_burn;
            tbb::concurrent_vector <Particle*> ok_list, burn_list, will_burn;
//            std::mutex mtx_ok, mtx_burn;
//            segment_mutex *mtx = 0;
        };


        void DoSegment(Segment&);
        void DoSegment(Segment*);
        Segment* GetSeg(int x, int y);
        Segment* SegmentByPoint(double x, double z);
        Segment* GetSegment(double x, double z) const;
        void BurnSegment(Segment*);
        void UpdateSegments();


        int grid_count_x, grid_count_z, grid_count;


    private:

        std::vector<Segment> grids;
        //segment_mutex *mutex_mem = 0;


        Segment **grid = 0, *grid_mem = 0;

        double grid_count_x_percent, grid_count_z_percent;



        void SetSegmentsGrid(double);
        void ClearSegments();

        int GetSegmentX(double x) const;
        int GetSegmentZ(double z) const;

        bool ParticleInBurnSegment(Particle* particle, int seg_x, int seg_z);
        bool CheckSegmentBurn(int seg_x, int seg_z);

        tbb::concurrent_vector <Segment*> burn_segments;



    private:
        int particles_count = 0;
    public:
        int size() const {
            return particles_count;
        }

        ~Segments() {
            delete[]grid;
            delete[]grid_mem;
            delete[]last_particles;
        }




    };


} /* namespace ps */


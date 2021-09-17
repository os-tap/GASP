#pragma once

typedef float coord_t;


namespace ps {


    class Particle {

    public:


        Particle(coord_t x, coord_t z, coord_t speed, coord_t burn_radius);
        void SetBurnRadius(coord_t burn_radius);

        enum class State { OK, WARM, BURN, WAVE, SAGE, DIED };
        State state = State::OK;

        unsigned burn_counter = 0, warm_counter = 0, wave_counter = 0, sage_counter = 0;

        coord_t x, z, speed, x_speed = 0;
        coord_t burn_radius;
        coord_t burn_radius_2 = burn_radius * burn_radius;
        //unsigned seg_x, seg_z;

        //Particle(coord_t x, coord_t z, coord_t speed);


        //get coordinates
        coord_t _x()const;
        coord_t _z()const;

        //get state
        State getState()const;

        bool isBurn() const;
        bool isOk() const;
        void setBurn();

        void Step();
        void Move(const double dx, const double dz);
        coord_t Distance(const Particle &) const;
        bool Cross(const Particle &) const;
        bool CrossBurn(const Particle &) const;

        coord_t Distance(const Particle *) const;
        bool Cross(const Particle *) const;
        bool CrossBurn(const Particle *) const;

    };

} /* namespace ps */


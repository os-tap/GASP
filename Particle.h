#pragma once

typedef float coord_t;


namespace ps {


    class Particle {

    public:

        enum class State {
            OK = 0,
            WARM = 1,
            BURN = 2,
            WAVE = 3,
            SAGE = 4,
            DIED = 5
        };
        State state = State::OK;

        Particle(coord_t _x, coord_t _z, coord_t _speed, coord_t _burn_radius) :
            x(_x), z(_z), speed(_speed), burn_radius(_burn_radius) {};

        Particle(coord_t _x, coord_t _z, State _st) :
            x(_x), z(_z), state(_st) 
        {
            if (_st == Particle::State::BURN) burn_counter = 1;
            if (_st == Particle::State::SAGE) burn_counter = 2;
        }

        void SetBurnRadius(coord_t burn_radius);

        

        unsigned short burn_counter = 0, warm_counter = 0, wave_counter = 0, sage_counter = 0;

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

        const bool isBurn();
        const bool isOk();
        void setBurn();

        void Step();
        void Move();
        void Move(double dz, double dx);
        coord_t Distance(const Particle &) const;
        bool Cross(const Particle &) const;
        bool CrossBurn(const Particle &) const;

        coord_t Distance(const Particle *) const;
        bool Cross(const Particle *) const;
        bool CrossBurn(const Particle *) const;

    };

} /* namespace ps */


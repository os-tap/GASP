#pragma once

typedef float cord_t;


namespace ps {


    class Particle {

    public:

        cord_t x, y=0, z, speed;
        cord_t burn_radius = 0;
        cord_t burn_radius_2 = burn_radius * burn_radius;

        Particle(cord_t _x, cord_t _z, cord_t _speed, cord_t _burn_radius) :
            x(_x), y(0), z(_z), speed(_speed), burn_radius(_burn_radius) { };

        Particle(cord_t _x, cord_t _y, cord_t _z, cord_t _speed, cord_t _burn_radius) :
            x(_x), y(_y), z(_z), speed(_speed), burn_radius(_burn_radius) { };

        enum class State { OK, WARM, BURN, WAVE, SAGE, DIED };
        State state = State::OK;

        unsigned short burn_counter = 0, warm_counter = 0, wave_counter = 0, sage_counter = 0;

        //unsigned seg_x, seg_z;

        //Particle(cord_t x, cord_t z, cord_t speed);


        //get coordinates
        cord_t _x()const;
        cord_t _y()const;
        cord_t _z()const;

        //get state
        State getState()const;

        const bool isBurn();
        const bool isOk();
        void setBurn();


        void Move();
        cord_t Distance(const Particle &) const;
        bool Cross(const Particle &) const;
        bool CrossBurn(const Particle &) const;

        cord_t Distance(const Particle *) const;
        bool Cross(const Particle *) const;
        bool CrossBurn(const Particle *) const;

    };

} /* namespace ps */


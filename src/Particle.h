#pragma once
#include "Point.hpp"


typedef float coord_t;


namespace ps {


    class Particle : public Point {

    public:
        Particle(coord_t _x, coord_t _y, bool is_visible = true) : Point(_x,_y), visible(is_visible) {}

        enum class State { OK, WARM, BURN, WAVE, SAGE, DIED };
        State state = State::OK;

        unsigned burn_counter = 0;
        unsigned warm_counter = 0;

        bool visible;

        //get coordinates
        coord_t _x()const;
        coord_t _z()const;

        //get state
        State getState()const;

        bool isBurn() const;
        bool isOk() const;
        void setBurn();

        void Move(const double dx, const double dz);

    };

} /* namespace ps */


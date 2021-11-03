#pragma once
#include "Point.hpp"

typedef float coord_t;


namespace ps {


    class Particle : public Point {

    public:
        Particle(coord_t _x, coord_t _y) : Point(_x, _y) {}
        Particle(coord_t _x, coord_t _y, bool isVisible) : Point(_x, _y), visible(isVisible) {}

        enum class State { OK, WARM, BURN, WAVE, SAGE, DIED };
        State state = State::OK;

        unsigned burn_counter = 0, warm_counter = 0;
        bool visible = false;

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


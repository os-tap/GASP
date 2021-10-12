#include "Particle.h"

namespace ps {


    void Particle::setBurn() {
        state = State::WARM;
        burn_counter = 1;
    }


    bool Particle::isBurn() const {
        return state == State::BURN;
    }
    bool Particle::isOk() const {
        return state == State::OK;
    }



    void Particle::Move(const double dx, const double dz) {
        z += dz;// *(state != State::WAVE);
        x += dx;
    }


    coord_t Particle::_x() const {
        return x;
    }

    coord_t Particle::_z() const {
        return z;
    }

    Particle::State Particle::getState() const {
        return state;
    }

} /* namespace ps */

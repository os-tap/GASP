#include "Particle.h"

namespace ps {





    const bool Particle::isBurn() {
        return state == State::BURN;
    }
    const bool Particle::isOk() {
        return state == State::OK;
    }

    void Particle::setBurn() {
        state = State::WARM;
        burn_counter = 1;
    }


    void Particle::Move() {
        z+= speed * (state != State::WAVE);
    }


    cord_t Particle::_x() const {
        return x;
    }

    cord_t Particle::_y() const {
        return y;
    }

    cord_t Particle::_z() const {
        return z;
    }

    Particle::State Particle::getState() const {
        return state;
    }

    inline cord_t Particle::Distance(const Particle &p) const {
        return (x - p.x) * (x - p.x) + (y - p.y) * (y - p.y) + (z - p.z) * (z - p.z);
    }
    bool Particle::Cross(const Particle &p) const {
        return Distance(p) < p.burn_radius_2;
    }
    bool Particle::CrossBurn(const Particle &burn_particle) const {
        return Distance(burn_particle) < burn_particle.burn_radius_2;
    }

    inline cord_t Particle::Distance(const Particle* p) const {
        return (x - p->x) * (x - p->x) + (y - p->y) * (y - p->y) + (z - p->z) * (z - p->z);
    }
    bool Particle::Cross(const Particle* p) const {
        return Distance(p) < p->burn_radius_2;
    }
    bool Particle::CrossBurn(const Particle* p) const {
        return ((x - p->x) * (x - p->x) + (y - p->y) * (y - p->y) + (z - p->z) * (z - p->z)) < p->burn_radius_2;
    }

} /* namespace ps */

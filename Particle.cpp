#include "Particle.h"

namespace ps {


    /*Particle::Particle(coord_t _x, coord_t _z, coord_t _speed) :
        x(_x), z(_z), speed(_speed){}*/

    
    

    void Particle::SetBurnRadius(coord_t _burn_radius) {
        burn_radius = _burn_radius;
        burn_radius_2 = burn_radius * burn_radius;
    }



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
        z += speed;// *(state != State::WAVE);
        x += x_speed;
    }
    void Particle::Move(double dz, double dx) {
        z += dz;// *(state != State::WAVE);
        x += dx;
    }

    void Particle::Step()
    {
        //x_speed += 0.002;
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

    inline coord_t Particle::Distance(const Particle &p) const {
        return (x - p.x) * (x - p.x) + (z - p.z) * (z - p.z);
    }
    bool Particle::Cross(const Particle &p) const {
        return Distance(p) < p.burn_radius_2;
    }
    bool Particle::CrossBurn(const Particle &burn_particle) const {
        return Distance(burn_particle) < burn_particle.burn_radius_2;
    }

    inline coord_t Particle::Distance(const Particle* p) const {
        return (x - p->x) * (x - p->x) + (z - p->z) * (z - p->z);
    }
    bool Particle::Cross(const Particle* p) const {
        return Distance(p) < p->burn_radius_2;
    }
    bool Particle::CrossBurn(const Particle* p) const {
        return ((x - p->x) * (x - p->x) + (z - p->z) * (z - p->z)) < p->burn_radius_2;
    }

} /* namespace ps */

#include "Params.h"

namespace ps {
    class Field;
}

class Field
{
private:
    const ps::Params *P;

    double const_speed_x, const_speed_y;
    double base_speed;
public:
    Field(ps::Params&);
    double Vx(const double x, const double y) const;
    double Vy(const double x, const double y) const;

};

#pragma once

namespace ps
{
    typedef float coord_t;
    class Point
    {
    public:
        Point() {};
        Point(coord_t _x, coord_t _z) : x(_x), z(_z) {}
        coord_t x, z;
        inline coord_t Distance(const Point& p) const
        {
            return (x - p.x) * (x - p.x) + (z - p.z) * (z - p.z);
        }
        inline bool Cross(const Point& p, double r2) const
        {
            return Distance(p) < r2;
        }
    };

}
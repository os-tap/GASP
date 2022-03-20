#pragma once

namespace ps
{
    typedef float coord_t;
    class Point
    {
    public:
        Point(coord_t _x, coord_t _z) : x(_x), z(_z) {}
        coord_t x, z;
        inline static coord_t Distance(const Point &a, const Point &b)
        {
            return (a.x - b.x) * (a.x - b.x) + (a.z - b.z) * (a.z - b.z);
        }
        inline static coord_t Cross(const Point &a, const Point &b, double r2)
        {
            return Distance(a, b) <= r2;
        }
        inline coord_t Distance(const Point &p) const
        {
            return Distance(*this, p);
        }
        inline bool Cross(const Point &p, double r2) const
        {
            return Distance(p) < r2;
        }
    };

}
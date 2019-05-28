//---------------------------------------------------------------------------
#ifndef sphere3H
#define sphere3H
#include "rect3.h"
//---------------------------------------------------------------------------

class Sphere3 {
public:

    Sphere3 (double x=0, double y=0, double z=0, double rad=1):
        m_c(x,y,z), m_r(rad) {}

    Sphere3 (const Point3 &c, double rad):
        m_c(c), m_r(rad) {}

    Sphere3 (const Sphere3 &s):
        m_c(s.getCenter()), m_r(s.getRadius()) {}

    Sphere3 (const Rect3 &r):
        m_c(r.getMid()), m_r((m_c-r.getP1()).length()) {}

    inline const Point3 &getCenter(void) const { return m_c; }
    inline double getRadius(void) const { return m_r; }

    inline int isInside(Point3 p) const {
        p-=m_c;
        return p.length() < m_r;
    }

    inline int intersect(Point3 p, double r) const {
        p-=m_c;
        return (p.length() <= (m_r + r));
    }

    inline int intersect(const Sphere3 &s) const {
        return intersect(s.getCenter(), s.getRadius());
    }

    inline int intersect(const Rect3 &r) const {
	    return ((m_c.getX() + m_r) >= r.minX() &&
		    (m_c.getX() - m_r) <= r.maxX() &&
		    (m_c.getY() + m_r) >= r.minY() &&
		    (m_c.getY() - m_r) <= r.maxY() &&
		    (m_c.getZ() + m_r) >= r.minZ() &&
		    (m_c.getZ() - m_r) <= r.maxZ()); 
    }

protected:
    Point3 m_c;
    double m_r;
};
#endif

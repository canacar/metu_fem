//---------------------------------------------------------------------------
#ifdef __BORLANDC__
#pragma hdrstop
#endif

#include "rect3.h"
//---------------------------------------------------------------------------

void Rect3::init(void)
{
    register double t;
    
    if (m_p1.getX() > m_p2.getX()){
        t=m_p1.getX();
        m_p1.setX(m_p2.getX());
        m_p2.setX(t);
    }

    if (m_p1.getY() > m_p2.getY()){
        double t=m_p1.getY();
        m_p1.setY(m_p2.getY());
        m_p2.setY(t);
    }

    if (m_p1.getZ() > m_p2.getZ()){
        double t=m_p1.getZ();
        m_p1.setZ(m_p2.getZ());
        m_p2.setZ(t);
    }
}

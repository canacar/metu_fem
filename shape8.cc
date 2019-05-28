/* $Id: shape8.cc,v 1.2 2008/01/28 07:36:28 canacar Exp $ */
/* 
 * This file is part of the EMSI Tools Package developed at the
 * Brain Research Laboratory, Middle East Technical University
 * Department of Electrical and Electronics Engineering.
 *
 * Copyright (C) 2008 Zeynep Akalin Acar
 * Copyright (C) 2008 Can Erkin Acar
 * Copyright (C) 2008 Nevzat G. Gencer
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */

#include "shape8.h"
#include <math.h>

int FEShape8::m_surfnodes8[] = {
	 0, 1, 2, 3, 
	 0, 1, 4, 5,
	 1, 2, 5, 6,
	 2, 3, 6, 7,
	 0, 3, 4, 7,
	 4, 5, 6, 7
};


FEShape8::FEShape8(void):FEShape(8, 6, 4)
{
	setupVolGP();

	m_numsgp=9;
	m_sgp=new FEShape::GP3D[m_numsgp*m_numfaces];

	for (int s=0; s<m_numfaces; s++)
		setupSurfGP(s);

	m_surfnodes=m_surfnodes8;
}

void
FEShape8::setupVolGP()
{
	m_numvgp=27;
	m_vgp=new FEShape::GP3D[m_numvgp];

	int seq=0;
	for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++) {
			for (int k=0; k<3; k++) {
				m_vgp[seq].i=m_gpos[i];
				m_vgp[seq].j=m_gpos[j];
				m_vgp[seq].k=m_gpos[k];
				m_vgp[seq].wt=m_gval[i] *
					      m_gval[j] *
					      m_gval[k];
				seq++;
			}
		}
	}

}

void
FEShape8::setupSurfGP(int s)
{
	assert (s>=0 && s<m_numfaces);

	int seq=s*m_numsgp;

	for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++) {

			double x,y,z;

			switch (s) {
			case 0:
				z=-1;
				y=m_gpos[i];
				x=m_gpos[j];
				break;
			case 1:
				y=-1;
				z=m_gpos[i];
				x=m_gpos[j];
				break;
			case 2:
				x=1;
				z=m_gpos[i];
				y=m_gpos[j];
				break;
			case 3:
				y=1;
				z=m_gpos[i];
				x=m_gpos[j];
				break;
			case 4:
				x=-1;
				z=m_gpos[i];
				y=m_gpos[j];
				break;
			case 5:
				z=1;
				y=m_gpos[i];
				x=m_gpos[j];
				break;
			default:
				/* NOT POSSIBLE */
				assert(0);
				return;
			}

			m_sgp[seq].i=x;
			m_sgp[seq].j=y;
			m_sgp[seq].k=z;
			m_sgp[seq].wt=m_gval[i]*m_gval[j];

			seq++;
		}
	}
}


FEShape8::~FEShape8()
{
	if (m_vgp) delete[] m_vgp;
	if (m_sgp) delete[] m_sgp;
}

void
FEShape8::Shape(double x, double y, double z)
{
	double xm, xp, ym, yp, zm, zp;
 
	xm=(1-x)/2;
	ym=(1-y)/2;
	zm=(1-z)/2;
 
	xp=(1+x)/2;
	yp=(1+y)/2;
	zp=(1+z)/2;

	m_shape[0]=xm*ym*zm;
	m_shape[1]=xp*ym*zm;
	m_shape[2]=xp*yp*zm;
	m_shape[3]=xm*yp*zm;
	m_shape[4]=xm*ym*zp;
	m_shape[5]=xp*ym*zp;
	m_shape[6]=xp*yp*zp;
	m_shape[7]=xm*yp*zp;
 
	m_deriv[0][0]=-0.5*ym*zm;
	m_deriv[1][0]= 0.5*ym*zm;
	m_deriv[2][0]= 0.5*yp*zm;
	m_deriv[3][0]=-0.5*yp*zm;
	m_deriv[4][0]=-0.5*ym*zp;
	m_deriv[5][0]= 0.5*ym*zp;
	m_deriv[6][0]= 0.5*yp*zp;
	m_deriv[7][0]=-0.5*yp*zp;

	m_deriv[0][1]=xm*-0.5*zm;
	m_deriv[1][1]=xp*-0.5*zm;
	m_deriv[2][1]=xp* 0.5*zm;
	m_deriv[3][1]=xm* 0.5*zm;
	m_deriv[4][1]=xm*-0.5*zp;
	m_deriv[5][1]=xp*-0.5*zp;
	m_deriv[6][1]=xp* 0.5*zp;
	m_deriv[7][1]=xm* 0.5*zp;
 
	m_deriv[0][2]=xm*ym*-0.5;
	m_deriv[1][2]=xp*ym*-0.5;
	m_deriv[2][2]=xp*yp*-0.5;
	m_deriv[3][2]=xm*yp*-0.5;
	m_deriv[4][2]=xm*ym* 0.5;
	m_deriv[5][2]=xp*ym* 0.5;
	m_deriv[6][2]=xp*yp* 0.5;
	m_deriv[7][2]=xm*yp* 0.5;
}


/* shape function must already be calculated for required point
   (m_shape and m_deriv set) the function basically calculates
   a correct normal vector (vn*) given shape funcs. and surface */
 
int
FEShape8::SurfInt(const point *pts, int surf, 
		   double &vnx, double &vny, double &vnz)
{
	int it, jt;
	double J[3][3];

	switch(surf){
	case 0:
		it=1;
		jt=0;
		break;
	case 1:
		it=0;
		jt=2;
		break;
	case 2:
		it=1;
		jt=2;
		break;
	case 3:
		it=2;
		jt=0;
		break;
	case 4:
		it=2;
		jt=1;
		break;
	case 5:
		it=0;
		jt=1;
		break;
	default:
		/* NOT POSSIBLE */
		assert(0);
		return 1;
	}

	/* calculate only the jacobian, not inverse */
	if(Jacob(pts, NULL, J, NULL)) return 1;
 
	vnx=(J[it][1]*J[jt][2]-J[it][2]*J[jt][1])/*/ty*/;
	vny=(J[it][2]*J[jt][0]-J[it][0]*J[jt][2])/*/ty*/;
	vnz=(J[it][0]*J[jt][1]-J[it][1]*J[jt][0])/*/ty*/;
 
	/* we can calculate tempy through vn* more easily */
 
	return 0;
}

int
FEShape8::SurfIntR(const point *pts, const point sens, int surf,
		    double &bx, double &by, double &bz)
{
	int it, jt;
	double J[3][3];
	double vnx, vny, vnz;

	switch(surf){
	case 0:
		it=1;
		jt=0;
		break;
	case 1:
		it=0;
		jt=2;
		break;
	case 2:
		it=1;
		jt=2;
		break;
	case 3:
		it=2;
		jt=0;
		break;
	case 4:
		it=2;
		jt=1;
		break;
	case 5:
		it=0;
		jt=1;
		break;
	default:
		/* NOT POSSIBLE */
		assert(0);
		return 1;
	}

// calculate only the jacobian, not inverse 
	if(Jacob(pts, NULL, J, NULL)) return 1;
 
	vnx=(J[it][1]*J[jt][2]-J[it][2]*J[jt][1]); // /ty;
	vny=(J[it][2]*J[jt][0]-J[it][0]*J[jt][2]); // /ty;
	vnz=(J[it][0]*J[jt][1]-J[it][1]*J[jt][0]); // /ty;
 
// global coordinates of the sensor point 
	double rx=sens[0];
	double ry=sens[1];
	double rz=sens[2];
	for(int n=0; n<8; n++){
	  rx-=m_shape[n]*pts[n][0];
	  ry-=m_shape[n]*pts[n][1];
	  rz-=m_shape[n]*pts[n][2];
	}
// now we have (R/R^3) 
	double r3=sqrt(rx*rx+ry*ry+rz*rz);
	r3=r3*r3*r3;
	rx/=r3;
	ry/=r3;
	rz/=r3;
 
// finally obtain b* 
	bx=ry*vnz-vny*rz;
	by=vnx*rz-rx*vnz;
	bz=rx*vny-vnx*ry;
 
	return 0;
}

 

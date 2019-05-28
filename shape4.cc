/* $Id: shape4.cc,v 1.2 2008/11/24 04:02:41 canacar Exp $ */
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

#include "shape4.h"
#include <math.h>
#include <stdio.h>


/*
 *       y
 *       |
 *       +1
 *      /|\
 *      || \
 *      ||  \
 *     / |   \    /
 *<-s0|  |    \ s2
 *    |  |     \
 *    |  |      \
 *    |  |       \0
 *    | 3+--------+- x
 *   /  /        /
 *  |  /      --
 *  | /    --
 *  |/   /  s1
 *  +---     |
 * /2        v
 *z
 *
 */
int FEShape4::m_surfnodes4[] = {
	 3, 2, 1,
	 0, 2, 3,
	 0, 3, 1,
	 0, 1, 2,
};


FEShape4::FEShape4(void):FEShape(4, 4, 3)
{
	m_surfnodes = m_surfnodes4;

	setupVolGP();

	m_numsgp = 3;
	m_sgp = new FEShape::GP3D[m_numsgp * m_numfaces];

	for (int s = 0; s < m_numfaces; s++)
		setupSurfGP(s);
}


/* 
 * Use Gauss Points of order 3 (5 points), from: O C Zienkiewicz, The
 * Finite Element Method Vol 1, 2000 p.223"
 */
void
FEShape4::setupVolGP()
{
	m_numvgp= 5;
	m_vgp = new FEShape::GP3D[m_numvgp];

	m_vgp[0].i = 0.25;
	m_vgp[0].j = 0.25;
	m_vgp[0].k = 0.25;
	m_vgp[0].wt = -4.0/30.0;

	for (int i = 1; i <= 4; i++) {
		m_vgp[i].i = (i == 1) ? 0.5 : (1.0/6.0);
		m_vgp[i].j = (i == 2) ? 0.5 : (1.0/6.0);
		m_vgp[i].k = (i == 3) ? 0.5 : (1.0/6.0);
		m_vgp[i].wt = 9.0/120.0;
	}
}

void
FEShape4::setupSurfGP(int s)
{
	assert (s>=0 && s<m_numfaces);
	double i[3], j[3], k[3];
	int n;

	int seq = s * m_numsgp;

	/* we cheat and use quadratic GP which are at the mid point of
	 * each triangle, thus easy to compute from surface coordinates.
	 */
	
	for (n = 0; n < 3; n++)
		i[n] = j[n] = k[n] = 0;

	switch (s) {
	case 0:
		/* 0, 0, 0.5 */
		k[0] = 0.5;
		/* 0, 0.5, 0 */
		j[1] = 0.5;
		/* 0, 0.5, 0.5 */
		j[2] = k[2] = 0.5;
		break;
	case 1:
		/* 0, 0, 0.5 */
		k[0] = 0.5;
		/* 0.5, 0, 0 */
		i[1] = 0.5;
		/* 0.5, 0, 0.5 */
		i[2] = k[2] = 0.5;
		break;
	case 2:
		/* 0, 0.5, 0 */
		j[0] = 0.5;
		/* 0.5, 0, 0 */
		i[1] = 0.5;
		/* 0.5, 0.5, 0 */
		i[2] = j[2] = 0.5;
		break;
	case 3:
		/* 0.5, 0, 0.5 */
		i[0] = k[0] = 0.5;
		/* 0, 0.5, 0.5 */
		j[1] = k[1] = 0.5;
		/* 0.5, 0.5, 0 */
		i[2] = j[2] = 0.5;
		break;
	default:
		/* NOT POSSIBLE */
		assert(0);
		return;
	}

	/* for each GP */
	for (n = 0; n < 3; n++, seq++) {
		m_sgp[seq].i = i[n];
		m_sgp[seq].j = j[n];
		m_sgp[seq].k = k[n];
		m_sgp[seq].wt = 1.0/6.0;
	}
}


FEShape4::~FEShape4()
{
	if (m_vgp)
		delete[] m_vgp;
	if (m_sgp)
		delete[] m_sgp;
}

void
FEShape4::Shape(double x, double y, double z)
{
	m_shape[0] = x;
	m_shape[1] = y;
	m_shape[2] = z;
	m_shape[3] = 1 - x - y - z;
 
	m_deriv[0][0] = 1;
	m_deriv[1][0] = 0;
	m_deriv[2][0] = 0;
	m_deriv[3][0] = -1;

	m_deriv[0][1] = 0;
	m_deriv[1][1] = 1;
	m_deriv[2][1] = 0;
	m_deriv[3][1] = -1;

	m_deriv[0][2] = 0;
	m_deriv[1][2] = 0;
	m_deriv[2][2] = 1;
	m_deriv[3][2] = -1;
}


/* shape function must already be calculated for required point
   (m_shape and m_deriv set) the function basically calculates
   a correct normal vector (vn*) given shape funcs. and surface */
 
int
FEShape4::SurfInt(const point *pts, int surf, 
		   double &vnx, double &vny, double &vnz)
{
	double J[3][3];
	double v1[3], v2[3];

	/* calculate only the jacobian, not ainverse */
	if (Jacob(pts, NULL, J, NULL))
		return 1;
 
	switch(surf){
	case 0:
		for (int i = 0; i < 3; i++) {
			v1[i] = J[2][i];
			v2[i] = J[1][i];
		}
		break;
	case 1:
		for (int i = 0; i < 3; i++) {
			v1[i] = J[0][i];
			v2[i] = J[2][i];
		}
		break;
	case 2:
		for (int i = 0; i < 3; i++) {
			v1[i] = J[1][i];
			v2[i] = J[0][i];
		}
		break;
	case 3:
		for (int i = 0; i < 3; i++) {
			v1[i] = (J[1][i] - J[0][i]);
			v2[i] = (J[2][i] - J[0][i]);
		}
		break;
	default:
		/* NOT POSSIBLE */
		assert(0);
		return 1;
	}

	vnx=(v1[1]*v2[2] - v1[2]*v2[1]);
	vny=(v1[2]*v2[0] - v1[0]*v2[2]);
	vnz=(v1[0]*v2[1] - v1[1]*v2[0]);
 
	return 0;
}

int
FEShape4::SurfIntR(const point *pts, const point sens, int surf,
		    double &bx, double &by, double &bz)
{
	double J[3][3];
	double v1[3], v2[3];
	double vnx, vny, vnz;

	/* calculate only the jacobian, not inverse */
	if (Jacob(pts, NULL, J, NULL))
		return 1;

 
	switch(surf){
	case 0:
		for (int i = 0; i < 3; i++) {
			v1[i] = J[2][i];
			v2[i] = J[1][i];
		}
		break;
	case 1:
		for (int i = 0; i < 3; i++) {
			v1[i] = J[0][i];
			v2[i] = J[2][i];
		}
		break;
	case 2:
		for (int i = 0; i < 3; i++) {
			v1[i] = J[1][i];
			v2[i] = J[0][i];
		}
		break;
	case 3:
		for (int i = 0; i < 3; i++) {
			v1[i] = (J[1][i] - J[0][i]);
			v2[i] = (J[2][i] - J[0][i]);
		}
		break;
	default:
		/* NOT POSSIBLE */
		assert(0);
		return 1;
	}
 

	vnx=(v1[1]*v2[2] - v1[2]*v2[1])/2;
	vny=(v1[2]*v2[0] - v1[0]*v2[2])/2;
	vnz=(v1[0]*v2[1] - v1[1]*v2[0])/2;
 
// global coordinates of the sensor point 
	double rx = sens[0];
	double ry = sens[1];
	double rz = sens[2];

	for (int n = 0; n < 8; n++) {
		rx-=m_shape[n]*pts[n][0];
		ry-=m_shape[n]*pts[n][1];
		rz-=m_shape[n]*pts[n][2];
	}

// now we have (R/R^3) 
	double r3 = sqrt(rx*rx+ry*ry+rz*rz);
	r3 = r3*r3*r3;
	rx /= r3;
	ry /= r3;
	rz /= r3;
 
// finally obtain b* 
	bx = ry*vnz - vny*rz;
	by = vnx*rz - rx*vnz;
	bz = rx*vny - vnx*ry;

	return 0;
}

/* $Id$ */
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

#include "shape10.h"
#include <math.h>
#include <stdio.h>


/*
 *       y
 *       |
 *       +1
 *      /|\
 *      || \
 *      ||  \ 4
 *     / |   +    /
 *<-s0|  +8   \ s2
 *    |  |     \
 *   5+  |      \
 *    |  |  7    \0
 *    | 3+--+-----+- x
 *   /  /        /
 *  |  /    6 --
 *  | +9   -+
 *  |/   /  s1
 *  +---     |
 * /2        v
 *z
 *
 */
int FEShape10::m_surfnodes10[] = {
	3, 9, 2, 5, 1, 8,
	0, 6, 2, 9, 3, 7,
	0, 7, 3, 8, 1, 4,
	0, 4, 1, 5, 2, 6,
};

struct GP {
	double i, j, k;
	double wt;
};

/* 
 * Use Gauss Points of order 4 (6 points), from: Solin et.al.,
 * "Higher-Order Finite Element Methods", 2004, p.238"
 * 
 */
#define NUM_SGP 6

static struct GP sgp[NUM_SGP] = {
	{ 4.45948491e-01, 4.45948491e-01, 0, 4.46763179e-01},
	{ 4.45948491e-01, 1.08103018e-01, 0, 4.46763179e-01},
	{ 1.08103018e-01, 4.45948491e-01, 0, 4.46763179e-01},
	{ 9.15762135e-02, 9.15762135e-02, 0, 2.19903487e-01},
	{ 9.15762135e-02, 8.16847573e-01, 0, 2.19903487e-01},
	{ 8.16847573e-01, 9.15762135e-02, 0, 2.19903487e-01}
};

/* 
 * Use Gauss Points of order 4 (11 points), from: Solin et.al.,
 * "Higher-Order Finite Element Methods", 2004, p.247"
 * 
 */
/* This is the normalized table from the book [-1 1]
 * 
 * -5.00000000e-01 -5.00000000e-01 -5.00000000e-01 -1.05244444e-01
 * -8.57142857e-01 -8.57142857e-01 -8.57142857e-01 6.09777777e-02
 * -8.57142857e-01 -8.57142857e-01 5.71428571e-01 6.09777777e-02
 * -8.57142857e-01 5.71428571e-01 -8.57142857e-01 6.09777777e-02
 *  5.71428571e-01 -8.57142857e-01 -8.57142857e-01 6.09777777e-02
 * -2.01192848e-01 -2.01192848e-01 -7.98807152e-01 1.99111111e-01
 * -2.01192848e-01 -7.98807152e-01 -2.01192848e-01 1.99111111e-01
 * -7.98807152e-01 -2.01192848e-01 -2.01192848e-01 1.99111111e-01
 * -2.01192848e-01 -7.98807152e-01 -7.98807152e-01 1.99111111e-01
 * -7.98807152e-01 -2.01192848e-01 -7.98807152e-01 1.99111111e-01
 * -7.98807152e-01 -7.98807152e-01 -2.01192848e-01 1.99111111e-01
*/

#define NUM_VGP 11
static struct GP vgp[NUM_VGP] = {
	{ 2.50000000e-01, 2.50000000e-01, 2.50000000e-01, -1.05244444e-01},
	{ 7.14285714e-02, 7.14285714e-02, 7.14285714e-02,  6.09777777e-02},
	{ 7.14285714e-02, 7.14285714e-02, 7.85714286e-01,  6.09777777e-02},
	{ 7.14285714e-02, 7.85714286e-01, 7.14285714e-02,  6.09777777e-02},
	{ 7.85714286e-01, 7.14285714e-02, 7.14285714e-02,  6.09777777e-02},
	{ 3.99403576e-01, 3.99403576e-01, 1.00596424e-01,  1.99111111e-01},
	{ 3.99403576e-01, 1.00596424e-01, 3.99403576e-01,  1.99111111e-01},
	{ 1.00596424e-01, 3.99403576e-01, 3.99403576e-01,  1.99111111e-01},
	{ 3.99403576e-01, 1.00596424e-01, 1.00596424e-01,  1.99111111e-01},
	{ 1.00596424e-01, 3.99403576e-01, 1.00596424e-01,  1.99111111e-01},
	{ 1.00596424e-01, 1.00596424e-01, 3.99403576e-01,  1.99111111e-01}
};


FEShape10::FEShape10(void):FEShape(10, 4, 6)
{
	m_surfnodes = m_surfnodes10;

	setupVolGP();

	m_numsgp = NUM_SGP;
	m_sgp = new FEShape::GP3D[m_numsgp * m_numfaces];

	for (int s = 0; s < m_numfaces; s++)
		setupSurfGP(s);
}


void
FEShape10::setupVolGP()
{
	m_numvgp= NUM_VGP;
	m_vgp = new FEShape::GP3D[m_numvgp];
	for (int n = 0; n < NUM_VGP; n++) {
		m_vgp[n].i = vgp[n].i;
		m_vgp[n].j = vgp[n].j;
		m_vgp[n].k = vgp[n].k;
		m_vgp[n].wt = vgp[n].wt / 8;
	}
}

void
FEShape10::setupSurfGP(int s)
{
	assert (s>=0 && s<m_numfaces);
	int seq = s * m_numsgp;
	int n;

	for (n = 0; n < NUM_SGP; n++, seq++) {
		switch (s) {
		case 0:
			m_sgp[seq].i = 0;
			m_sgp[seq].j = sgp[n].i;
			m_sgp[seq].k = sgp[n].j;
			break;
		case 1:
			m_sgp[seq].i = sgp[n].i;
			m_sgp[seq].j = 0;
			m_sgp[seq].k = sgp[n].j;
			break;
		case 2:
			m_sgp[seq].i = sgp[n].i;
			m_sgp[seq].j = sgp[n].j;
			m_sgp[seq].k = 0;
			break;
		case 3:
			m_sgp[seq].i = sgp[n].i;
			m_sgp[seq].j = sgp[n].j;
			m_sgp[seq].k = 1 - (sgp[n].i + sgp[n].j);
			break;
		default:
			/* NOT POSSIBLE */
			assert(0);
			return;
		}
		m_sgp[seq].wt = sgp[n].wt / 4;
	}
}


FEShape10::~FEShape10()
{
	if (m_sgp)
		delete[] m_sgp;
	if (m_vgp)
		delete[] m_vgp;
}

void
FEShape10::Shape(double x, double y, double z)
{
	double x2 = x * x;
	double y2 = y * y;
	double z2 = z * z;

	double v = 1 - (x + y + z);
	double v2 = v * v; /* x2 + y2 + z2 + 2xy + 2xz + 2yz -2x -2y -2z + 1*/

	m_shape[0] = 2 * x2 - x;
	m_shape[1] = 2 * y2 - y;
	m_shape[2] = 2 * z2 - z;
	m_shape[3] = 2 * v2 - v;
	m_shape[4] = 4 * x * y;
	m_shape[5] = 4 * y * z;
	m_shape[6] = 4 * x * z;
	m_shape[7] = 4 * x * v;
 	m_shape[8] = 4 * y * v;
 	m_shape[9] = 4 * z * v;

	m_deriv[0][0] = 4 * x - 1;
	m_deriv[1][0] = 0;
	m_deriv[2][0] = 0;
	m_deriv[3][0] = 1 - 4 * v;
	m_deriv[4][0] = 4 * y;
	m_deriv[5][0] = 0;
	m_deriv[6][0] = 4 * z;
	m_deriv[7][0] = 4 * (v - x);
	m_deriv[8][0] = -4 * y;
	m_deriv[9][0] = -4 * z;

	m_deriv[0][1] = 0;
	m_deriv[1][1] = 4 * y - 1;
	m_deriv[2][1] = 0;
	m_deriv[3][1] = 1 - 4 * v;
	m_deriv[4][1] = 4 * x;
	m_deriv[5][1] = 4 * z;
	m_deriv[6][1] = 0;
	m_deriv[7][1] = -4 * x;
	m_deriv[8][1] = 4 * (v - y);
	m_deriv[9][1] = -4 * z;

	m_deriv[0][2] = 0;
	m_deriv[1][2] = 0;
	m_deriv[2][2] = 4 * z - 1;
	m_deriv[3][2] = 1 - 4 * v;
	m_deriv[4][2] = 0;
	m_deriv[5][2] = 4 * y;
	m_deriv[6][2] = 4 * x;
	m_deriv[7][2] = -4 * x;
	m_deriv[8][2] = - 4 * y;
	m_deriv[9][2] = 4 * (v - z);
}


/* shape function must already be calculated for required point
   (m_shape and m_deriv set) the function basically calculates
   a correct normal vector (vn*) given shape funcs. and surface */
 
int
FEShape10::SurfInt(const point *pts, int surf, 
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
FEShape10::SurfIntR(const point *pts, const point sens, int surf,
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

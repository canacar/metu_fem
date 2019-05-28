/* $Id: shape20.cc,v 1.2 2008/01/28 07:36:28 canacar Exp $ */
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

#include "shape20.h"
#include <math.h>

int FEShape20::m_surfnodes20[] = {
	 0, 1, 2, 3, 4, 5, 6, 7,
	 0, 1, 2, 8, 9,12,13,14,
	 2, 3, 4, 9,10,14,15,16,
	 4, 5, 6,10,11,16,17,18,
	 0, 6, 7, 8,11,12,18,19,
	12,13,14,15,16,17,18,19
};


FEShape20::FEShape20(void):FEShape(20, 6, 8)
{

	setupVolGP();

	m_numsgp=9;
	m_sgp=new FEShape::GP3D[m_numsgp*m_numfaces];

	for (int s=0; s<m_numfaces; s++)
		setupSurfGP(s);

	m_surfnodes=m_surfnodes20;
}

void
FEShape20::setupVolGP()
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
FEShape20::setupSurfGP(int s)
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


FEShape20::~FEShape20()
{
	if (m_vgp) delete[] m_vgp;
	if (m_sgp) delete[] m_sgp;
}

void
FEShape20::Shape(double x, double y, double z)
{
	double xy,xz,yz, xx,yy,zz;
	double xxy, xxz, yyx, yyz, zzx, zzy, xyz;
	double xxyz, yyxz, zzxy;
	double x2, y2, z2, xy2, xz2, yz2, xyz2;
 
	xy=x*y;
	xz=x*z;
	yz=y*z;
	xx=x*x;
	yy=y*y;
	zz=z*z;
	xxy=xx*y;
	xxz=xx*z;
	yyx=yy*x;
	yyz=yy*z;
	zzx=zz*x;
	zzy=zz*y;
	xyz=x*y*z;
	xxyz=xx*y*z;
	yyxz=yy*x*z;
	zzxy=zz*x*y;
	x2=x*2;
	xy2=xy*2;
	xz2=xz*2;
	xyz2=xyz*2;
	y2=y*2;
	yz2=yz*2;
	z2=z*2;

	m_shape[0]=(-2.0+x+y+z+xx+yy+zz-xxy-xxz-yyx-yyz-zzx-zzy-xyz+xxyz+yyxz+zzxy)*0.125;
	m_shape[1]=(1.0-y-z+yz-xx+xxy+xxz-xxyz)*0.25;
	m_shape[2]=(-2.0-x+y+z+xx+yy+zz-xxy-xxz+yyx-yyz+zzx-zzy+xyz+xxyz-yyxz-zzxy)*0.125;
	m_shape[3]=(1.0+x-z-xz-yy-yyx+yyz+yyxz)*0.25;
	m_shape[4]=(-2.0-x-y+z+xx+yy+zz+xxy-xxz+yyx-yyz+zzx+zzy-xyz-xxyz-yyxz+zzxy)*0.125;
	m_shape[5]=(1.0+y-z-yz-xx-xxy+xxz+xxyz)*0.25;
	m_shape[6]=(-2.0+x-y+z+xx+yy+zz+xxy-xxz-yyx-yyz-zzx+zzy+xyz-xxyz+yyxz-zzxy)*0.125;
	m_shape[7]=(1.0-x-z+xz-yy+yyx+yyz-yyxz)*0.25;
	m_shape[8]=(1.0-x-y+xy-zz+zzx+zzy-zzxy)*0.25;
	m_shape[9]=(1.0+x-y-xy-zz-zzx+zzy+zzxy)*0.25;
	m_shape[10]=(1.0+x+y+xy-zz-zzx-zzy-zzxy)*0.25;
	m_shape[11]=(1.0-x+y-xy-zz+zzx-zzy+zzxy)*0.25;
	m_shape[12]=(-2.0+x+y-z+xx+yy+zz-xxy+xxz-yyx+yyz-zzx-zzy+xyz-xxyz-yyxz+zzxy)*0.125;
	m_shape[13]=(1.0-y+z-yz-xx+xxy-xxz+xxyz)*0.25;
	m_shape[14]=(-2.0-x+y-z+xx+yy+zz-xxy+xxz+yyx+yyz+zzx-zzy-xyz-xxyz+yyxz-zzxy)*0.125;
	m_shape[15]=(1.0+x+z+xz-yy-yyx-yyz-yyxz)*0.25;
	m_shape[16]=(-2.0-x-y-z+xx+yy+zz+xxy+xxz+yyx+yyz+zzx+zzy+xyz+xxyz+yyxz+zzxy)*0.125;
	m_shape[17]=(1.0+y+z+yz-xx-xxy-xxz-xxyz)*0.25;
	m_shape[18]=(-2.0+x-y-z+xx+yy+zz+xxy+xxz-yyx+yyz-zzx+zzy-xyz+xxyz-yyxz-zzxy)*0.125;
	m_shape[19]=(1.0-x+z-xz-yy+yyx-yyz+yyxz)*0.25;
 
	m_deriv[0][0]=(1.0+x2-xy2-xz2-yy-zz-yz+xyz2+yyz+zzy)*0.125;
	m_deriv[0][1]=(1.0+y2-xx-xy2-yz2-zz-xz+xxz+xyz2+zzx)*0.125;
	m_deriv[0][2]=(1.0+z2-xx-yy-xz2-yz2-xy+xxy+yyx+xyz2)*0.125;
	m_deriv[1][0]=(-x2+xy2+xz2-xyz2)*0.25;
	m_deriv[1][1]=(-1.0+z+xx-xxz)*0.25;
	m_deriv[1][2]=(-1.0+y+xx-xxy)*0.25;
	m_deriv[2][0]=(-1.0+x2-xy2-xz2+yy+zz+yz+xyz2-yyz-zzy)*0.125;
	m_deriv[2][1]=(1.0+y2-xx+xy2-yz2-zz+xz+xxz-xyz2-zzx)*0.125;
	m_deriv[2][2]=(1.0+z2-xx-yy+xz2-yz2+xy+xxy-yyx-xyz2)*0.125;
	m_deriv[3][0]=(1.0-z-yy+yyz)*0.25;
	m_deriv[3][1]=(-y2-xy2+yz2+xyz2)*0.25;
	m_deriv[3][2]=(-1.0-x+yy+yyx)*0.25;
	m_deriv[4][0]=(-1.0+x2+xy2-xz2+yy+zz-yz-xyz2-yyz+zzy)*0.125;
	m_deriv[4][1]=(-1.0+y2+xx+xy2-yz2+zz-xz-xxz-xyz2+zzx)*0.125;
	m_deriv[4][2]=(1.0+z2-xx-yy+xz2+yz2-xy-xxy-yyx+xyz2)*0.125;
	m_deriv[5][0]=(-x2-xy2+xz2+xyz2)*0.25;
	m_deriv[5][1]=(1.0-z-xx+xxz)*0.25;
	m_deriv[5][2]=(-1.0-y+xx+xxy)*0.25;
	m_deriv[6][0]=(1.0+x2+xy2-xz2-yy-zz+yz-xyz2+yyz-zzy)*0.125;
	m_deriv[6][1]=(-1.0+y2+xx-xy2-yz2+zz+xz-xxz+xyz2-zzx)*0.125;
	m_deriv[6][2]=(1.0+z2-xx-yy-xz2+yz2+xy-xxy+yyx-xyz2)*0.125;
	m_deriv[7][0]=(-1.0+z+yy-yyz)*0.25;
	m_deriv[7][1]=(-y2+xy2+yz2-xyz2)*0.25;
	m_deriv[7][2]=(-1.0+x+yy-yyx)*0.25;
	m_deriv[8][0]=(-1.0+y+zz-zzy)*0.25;
	m_deriv[8][1]=(-1.0+x+zz-zzx)*0.25;
	m_deriv[8][2]=(-z2+xz2+yz2-xyz2)*0.25;
	m_deriv[9][0]=(1.0-y-zz+zzy)*0.25;
	m_deriv[9][1]=(-1.0-x+zz+zzx)*0.25;
	m_deriv[9][2]=(-z2-xz2+yz2+xyz2)*0.25;
	m_deriv[10][0]=(1.0+y-zz-zzy)*0.25;
	m_deriv[10][1]=(1.0+x-zz-zzx)*0.25;
	m_deriv[10][2]=(-z2-xz2-yz2-xyz2)*0.25;
	m_deriv[11][0]=(-1.0-y+zz+zzy)*0.25;
	m_deriv[11][1]=(1.0-x-zz+zzx)*0.25;
	m_deriv[11][2]=(-z2+xz2-yz2+xyz2)*0.25;
	m_deriv[12][0]=(1.0+x2-xy2+xz2-yy-zz+yz-xyz2-yyz+zzy)*0.125;
	m_deriv[12][1]=(1.0+y2-xx-xy2+yz2-zz+xz-xxz-xyz2+zzx)*0.125;
	m_deriv[12][2]=(-1.0+z2+xx+yy-xz2-yz2+xy-xxy-yyx+xyz2)*0.125;
	m_deriv[13][0]=(-x2+xy2-xz2+xyz2)*0.25;
	m_deriv[13][1]=(-1.0-z+xx+xxz)*0.25;
	m_deriv[13][2]=(1.0-y-xx+xxy)*0.25;
	m_deriv[14][0]=(-1.0+x2-xy2+xz2+yy+zz-yz-xyz2+yyz-zzy)*0.125;
	m_deriv[14][1]=(1.0+y2-xx+xy2+yz2-zz-xz-xxz+xyz2-zzx)*0.125;
	m_deriv[14][2]=(-1.0+z2+xx+yy+xz2-yz2-xy-xxy+yyx-xyz2)*0.125;
	m_deriv[15][0]=(1.0+z-yy-yyz)*0.25;
	m_deriv[15][1]=(-y2-xy2-yz2-xyz2)*0.25;
	m_deriv[15][2]=(1.0+x-yy-yyx)*0.25;
	m_deriv[16][0]=(-1.0+x2+xy2+xz2+yy+zz+yz+xyz2+yyz+zzy)*0.125;
	m_deriv[16][1]=(-1.0+y2+xx+xy2+yz2+zz+xz+xxz+xyz2+zzx)*0.125;
	m_deriv[16][2]=(-1.0+z2+xx+yy+xz2+yz2+xy+xxy+yyx+xyz2)*0.125;
	m_deriv[17][0]=(-x2-xy2-xz2-xyz2)*0.25;
	m_deriv[17][1]=(1.0+z-xx-xxz)*0.25;
	m_deriv[17][2]=(1.0+y-xx-xxy)*0.25;
	m_deriv[18][0]=(1.0+x2+xy2+xz2-yy-zz-yz+xyz2-yyz-zzy)*0.125;
	m_deriv[18][1]=(-1.0+y2+xx-xy2+yz2+zz-xz+xxz-xyz2-zzx)*0.125;
	m_deriv[18][2]=(-1.0+z2+xx+yy-xz2+yz2-xy+xxy-yyx-xyz2)*0.125;
	m_deriv[19][0]=(-1.0-z+yy+yyz)*0.25;
	m_deriv[19][1]=(-y2+xy2-yz2+xyz2)*0.25;
	m_deriv[19][2]=(1.0-x-yy+yyx)*0.25;
}


/* shape function must already be calculated for required point
   (m_shape and m_deriv set) the function basically calculates
   a correct normal vector (vn*) given shape funcs. and surface */
 
int
FEShape20::SurfInt(const point *pts, int surf, 
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

/* Not yet active! */

int
FEShape20::SurfIntR(const point *pts, const point sens, int surf,
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
	for(int n=0; n<20; n++){
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


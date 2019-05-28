/* $Id: shape.cc,v 1.2 2008/01/28 07:36:28 canacar Exp $ */
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
#include <math.h>
#include "define.h"
#include "shape.h"

// Gauss point constants
const double
FEShape::m_gpos[3]={0.7745966692, 0, -0.7745966692};
 
const double
FEShape::m_gval[3]={0.5555555556,0.888888889,0.5555555556};


FEShape::FEShape(int nnodes, int nfaces, int nnodeface):
	m_numnodes(nnodes), m_numfaces(nfaces), m_nodeface(nnodeface)
{
	assert (nnodes >0 && nfaces > 0 && nnodeface > 0);

	m_shape=new double[nnodes];
	m_deriv=new point[nnodes];
	m_cartp=new point[nnodes];
}


FEShape::~FEShape()
{
	if (m_shape) delete[] m_shape;
	if (m_deriv) delete[] m_deriv;
	if (m_cartp) delete[] m_cartp;
}

void
FEShape::Local2Global(const point *pts, double &x,
		      double &y, double &z)
{
	assert (pts);

	Shape(x,y,z);

	x=y=z=0;

	for (int n=0; n<m_numnodes; n++) {
		x+=m_shape[n]*pts[n][0];
		y+=m_shape[n]*pts[n][1];
		z+=m_shape[n]*pts[n][2];
	}
}

double
FEShape::interp(const double *fld, double x,
		double y, double z)
{
	assert (fld);

	Shape(x,y,z);

	double v=0;
	for (int n=0; n<m_numnodes; n++)
		v+=m_shape[n]*fld[n];

	return v;
}

int
FEShape::grad(const point *pts, const double *fld, point &gr, double x, double y, double z)
{
	double J[3][3], IJ[3][3];

	assert (fld); 

	gr[0] = gr[1] = gr[2] = 0;

	Shape(x, y, z);

	if(Jacob(pts, NULL, J, IJ)) return 1;
 
	/* calculate cartesian product (IJ*Deriv) = dNi/dxyz */

	for(int n=0; n<m_numnodes;n++){
		for(int id=0; id<3; id++){
			double tmp=0;
			for(int jd=0; jd<3; jd++)
				tmp+=IJ[id][jd] * m_deriv[n][jd];
			gr[id] += tmp * fld[n];
		}
	}

	return 0;
}

int
FEShape::jacob(const point *pts,  double x, double y, double z,
	       double J[3][3], double IJ[3][3])
{
	Shape(x, y, z);
	return Jacob(pts, NULL, J, IJ);
}

int
FEShape::Jacob(const point *pts,  double *detjb,
	       double jcb[3][3], double ijacob[3][3])
{	
	int i,j,n;
	const double *p;
	double dj;

	if(jcb==NULL || pts == NULL) return 1;
	
	/* clear array */
	memset(jcb, 0, sizeof(double)*3*3);
	
	for(n=0; n<m_numnodes; n++){  /* calculate jacobian */
		p=pts[n];
		for(j=0; j<3; j++){
			for(i=0; i<3; i++){
				jcb[i][j]+=m_deriv[n][i]*(*p);
      			}
			p++;
		}
	}
 
	dj=jcb[0][0]*(jcb[1][1]*jcb[2][2]-jcb[2][1]*jcb[1][2])
	   -jcb[1][0]*(jcb[0][1]*jcb[2][2]-jcb[0][2]*jcb[2][1])
	   +jcb[2][0]*(jcb[0][1]*jcb[1][2]-jcb[0][2]*jcb[1][1]);

	if(dj<=0) return 1;
 
	if(detjb) *detjb=dj;
 
	if(ijacob==NULL) return 0;

	/* invert jacobian matrix */
	ijacob[0][0]=(jcb[1][1]*jcb[2][2]-jcb[2][1]*jcb[1][2])/dj;
	ijacob[0][1]=(jcb[2][1]*jcb[0][2]-jcb[0][1]*jcb[2][2])/dj;
	ijacob[0][2]=(jcb[0][1]*jcb[1][2]-jcb[1][1]*jcb[0][2])/dj;
	ijacob[1][0]=(jcb[2][0]*jcb[1][2]-jcb[1][0]*jcb[2][2])/dj;
	ijacob[1][1]=(jcb[0][0]*jcb[2][2]-jcb[2][0]*jcb[0][2])/dj;
	ijacob[1][2]=(jcb[1][0]*jcb[0][2]-jcb[0][0]*jcb[1][2])/dj;
	ijacob[2][0]=(jcb[1][0]*jcb[2][1]-jcb[2][0]*jcb[1][1])/dj;
	ijacob[2][1]=(jcb[2][0]*jcb[0][1]-jcb[0][0]*jcb[2][1])/dj;
	ijacob[2][2]=(jcb[0][0]*jcb[1][1]-jcb[1][0]*jcb[0][1])/dj;

	return 0;
}

int
FEShape::Matrix(const point *pts, double sigma, DFSMatrix &mat, double *vol)
{
	int id,jd,n;
	double Detjb, tmp, volsig;
	double J[3][3], IJ[3][3];
	double *ci, *cj;

	assert (mat._size == m_numnodes); 

	/* clear matrices */
	mat.clear();
 
	if(vol) *vol=0;

	int ngp=numVolGP();

	for (int g=0; g<ngp; g++) {
		double gi, gj, gk, gwt;
		getVolGP(g, gi, gj, gk, gwt);
		Shape(gi, gj, gk);

		if(Jacob(pts, &Detjb, J, IJ)) return 1;
 
        	/* calculate cartesian product (IJ*Deriv) = dNi/dxyz */

		for(n=0; n<m_numnodes;n++){
			for(id=0; id<3; id++){
				tmp=0;
				for(jd=0; jd<3; jd++)
					tmp+=IJ[id][jd] * m_deriv[n][jd];
				m_cartp[n][id]=tmp;
			}
		}
		volsig = Detjb * gwt;
		if(vol) *vol+=volsig;
		volsig*=sigma;

		/* update matrix entries */
		for(id=0; id<m_numnodes; id++){
			cj=&(m_cartp[0][0]);
			for(jd=0; jd<m_numnodes; jd++){
				ci=&(m_cartp[id][0]);
				tmp=0;
				for(n=0; n<3; n++)
					tmp+=(*(ci++)) * (*(cj++));
				mat.add(id, jd, tmp*volsig);
			}
		}
	}

	return 0;
}

int
FEShape::Matrix(const point *pts, double *sarray, DFSMatrix &mat, double *vol)
{
	int id,jd,n;
	double Detjb, tmp, volsig;
	double J[3][3], IJ[3][3];
	double *ci, *cj;

	assert (mat._size == m_numnodes); 

	/* clear matrices */
	mat.clear();
 
	if(vol) *vol=0;

	int ngp=numVolGP();

	for (int g=0; g<ngp; g++) {
		double gi, gj, gk, gwt;
		getVolGP(g, gi, gj, gk, gwt);
		Shape(gi, gj, gk);

		if(Jacob(pts, &Detjb, J, IJ)) return 1;
 
        	/* calculate cartesian product (IJ*Deriv) = dNi/dxyz */

	        double sigma=0;
		for(n=0; n<m_numnodes;n++){
			for(id=0; id<3; id++){
				tmp=0;
				for(jd=0; jd<3; jd++)
					tmp+=IJ[id][jd] * m_deriv[n][jd];
				m_cartp[n][id]=tmp;
			}
			sigma+=m_shape[n]*sarray[n];
		}
		volsig = Detjb * gwt;
		if(vol) *vol+=volsig;
		volsig*=sigma;

		/* update matrix entries */
		for(id=0; id<m_numnodes; id++){
			cj=&(m_cartp[0][0]);
			for(jd=0; jd<m_numnodes; jd++){
				ci=&(m_cartp[id][0]);
				tmp=0;
				for(n=0; n<3; n++)
					tmp+=(*(ci++)) * (*(cj++));
				mat.add(id, jd, tmp*volsig);
			}
		}
	}

	return 0;
}


int
FEShape::Matrix(const point *pts, double (*getgp)(double,double,double),
		DFSMatrix &mat, double *vol)
{
	int id,jd,n;
	double Detjb, tmp, volsig;
	double J[3][3], IJ[3][3];
	double *ci, *cj;

	assert (mat._size == m_numnodes); 
	assert (getgp);

	/* clear matrices */
	mat.clear();
 
	if(vol) *vol=0;

	int ngp=numVolGP();

	for (int g=0; g<ngp; g++) {
		double gi, gj, gk, gwt;
		getVolGP(g, gi, gj, gk, gwt);
		Shape(gi, gj, gk);

		if(Jacob(pts, &Detjb, J, IJ)) return 1;
 
        	/* calculate cartesian product (IJ*Deriv) = dNi/dxyz */

		double gx=gi, gy=gj, gz=gk;
		
		Local2Global(pts, gx, gy, gz);
	        double sigma = getgp(gx, gy, gz);

		for(n=0; n<m_numnodes;n++){
			for(id=0; id<3; id++){
				tmp=0;
				for(jd=0; jd<3; jd++)
					tmp+=IJ[id][jd] * m_deriv[n][jd];
				m_cartp[n][id]=tmp;
			}
		}

		volsig = Detjb * gwt;
		if(vol) *vol+=volsig;
		volsig*=sigma;

		/* update matrix entries */
		for(id=0; id<m_numnodes; id++){
			cj=&(m_cartp[0][0]);
			for(jd=0; jd<m_numnodes; jd++){
				ci=&(m_cartp[id][0]);
				tmp=0;
				for(n=0; n<3; n++)
					tmp+=(*(ci++)) * (*(cj++));
				mat.add(id, jd, tmp*volsig);
			}
		}
	}

	return 0;
}

int
FEShape::SurfaceArea(const point *pts, int surf, double &sa)
{
	int ngp=numSurfGP();

	sa = 0;
	for (int g=0; g<ngp; g++) {
		double gi, gj, gk, gwt;
		double vnx, vny, vnz;
		getSurfGP(surf, g, gi, gj, gk, gwt);

		Shape(gi, gj, gk);
		if(SurfInt(pts, surf, vnx, vny, vnz))
			return 1;
		
		sa += sqrt(vnx*vnx + vny*vny + vnz*vnz) * gwt;
	}

	return 0;
}

int
FEShape::Volume(const point *pts, double &vol)
{
	double Detjb;
	double J[3][3];
 
	vol=0;
	int ngp=numVolGP();

        for (int g=0; g<ngp; g++) {
                double gi, gj, gk, gwt;
                getVolGP(g, gi, gj, gk, gwt);
                Shape(gi, gj, gk); 
		if(Jacob(pts, &Detjb, J, 0)) return 1;
		vol+=Detjb*gwt;
	}

	return 0;
}

int
FEShape::Bound(const point Ji, const point *pts, double *rhs)
{
	double Detjb;
	double J[3][3], IJ[3][3];
 
	int ngp=numVolGP();

	memset(rhs, 0, m_numnodes * sizeof(double));

        for (int g=0; g<ngp; g++) {
                double gi, gj, gk, gwt;
                getVolGP(g, gi, gj, gk, gwt);
                Shape(gi, gj, gk); 
		if(Jacob(pts, &Detjb, J, IJ)) return 1;
		double dv=Detjb*gwt;

		for(int n=0; n<m_numnodes;n++){
			for(int id=0; id<3; id++){
				double tmp=0;
				for(int jd=0; jd<3; jd++)
					tmp+=IJ[id][jd] * m_deriv[n][jd];
				m_cartp[n][id]=tmp;
			}
		}
		
		for(int n=0; n<m_numnodes;n++){
			rhs[n]-= dv * (m_cartp[n][0] * Ji[0] +
                                       m_cartp[n][1] * Ji[1] +
                                       m_cartp[n][2] * Ji[2]);
                }

	}

	return 0;
}

// calls a callback function at each gp to obtain the value of the
// vector field (which can be computed analytically) intended to be
// used with the magnetic vector potential
int
FEShape::BoundCB(const cbVec cbfunc, const point *pts, double *rhs, void *arg)
{
	double Detjb;
	double J[3][3], IJ[3][3];
 
	int ngp=numVolGP();

	memset(rhs, 0, m_numnodes * sizeof(double));

        for (int g=0; g<ngp; g++) {
                double gi, gj, gk, gwt;
                getVolGP(g, gi, gj, gk, gwt);
                Shape(gi, gj, gk); 
		if(Jacob(pts, &Detjb, J, IJ)) return 1;
		double dv=Detjb*gwt;

		for(int n=0; n<m_numnodes;n++){
			for(int id=0; id<3; id++){
				double tmp=0;
				for(int jd=0; jd<3; jd++)
					tmp+=IJ[id][jd] * m_deriv[n][jd];
				m_cartp[n][id]=tmp;
			}
		}

		point V;

		// compute global coord of gp
		double x, y, z;
		x=y=z=0;
		for (int n=0; n<m_numnodes; n++) {
			x+=m_shape[n]*pts[n][0];
			y+=m_shape[n]*pts[n][1];
			z+=m_shape[n]*pts[n][2];
		}

		// invoke callback
		if (cbfunc(x, y, z, arg, V))
			return 1;

		for(int n=0; n<m_numnodes;n++){
			rhs[n] -= dv * (m_cartp[n][0] * V[0] +
					m_cartp[n][1] * V[1] +
					m_cartp[n][2] * V[2]);
                }

	}

	return 0;
}


int
FEShape::BoundS(const point Ji, const point *pts, double *rhs)
{
	int n;
	for(n=0; n<m_numnodes; n++) rhs[n]=0;
 
	if(Ji[0] || Ji[1] || Ji[2])
		for(n=0; n<m_numfaces; n++) 
			if(BoundSurf(Ji, pts, n, rhs)) return 1;
	return 0;
}

int
FEShape::BoundSurf(const point Ji, const point *pts, int surf, double *rhs)
{
	double vnx, vny, vnz;
	double Jn;

	int ngp=numSurfGP();

	for (int g=0; g<ngp; g++) {
		double gi, gj, gk, gwt;
		getSurfGP(surf, g, gi, gj, gk, gwt);

		Shape(gi, gj, gk);
		if(SurfInt(pts, surf, vnx, vny, vnz))
			return 1;

		Jn=Ji[0]*vnx+Ji[1]*vny+Ji[2]*vnz;
		for(int n=0; n<m_nodeface; n++){
			int p=SurfNode(surf, n);
			rhs[p]-=m_shape[p]*Jn*gwt;
		}
	}

	return 0;
}

int
FEShape::BoundSCB(const cbVec cbfunc, const point *pts, double *rhs, void *arg)
{
	int n;
	for(n=0; n<m_numnodes; n++) rhs[n]=0;
 
	for(n=0; n<m_numfaces; n++) 
		if(BoundSurfCB(cbfunc, pts, n, rhs, arg)) return 1;

	return 0;
}

int
FEShape::BoundSurfCB(const cbVec cbfunc, const point *pts,
		     int surf, double *rhs, void *arg)
{
	double vnx, vny, vnz;
	double Jn;

	int ngp=numSurfGP();

	for (int g=0; g<ngp; g++) {
		double gi, gj, gk, gwt;
		getSurfGP(surf, g, gi, gj, gk, gwt);

		Shape(gi, gj, gk);
		if(SurfInt(pts, surf, vnx, vny, vnz))
			return 1;

		point V;

		// compute global coord of gp
		double x, y, z;
		x=y=z=0;
		for (int n=0; n<m_numnodes; n++) {
			x+=m_shape[n]*pts[n][0];
			y+=m_shape[n]*pts[n][1];
			z+=m_shape[n]*pts[n][2];
		}

		// invoke callback
		if (cbfunc(x, y, z, arg, V))
			return 1;		

		Jn=V[0]*vnx+V[1]*vny+V[2]*vnz;
		for(int n=0; n<m_nodeface; n++){
			int p=SurfNode(surf, n);
			rhs[p]-=m_shape[p]*Jn*gwt;
		}
	}

	return 0;
}


int
FEShape::BoundYan(const point D, const point Ji, const point *pts, double *rhs)
{
	double Detjb;
	double J[3][3], IJ[3][3];

	Shape(D[0], D[1], D[2]);

	if(Jacob(pts, &Detjb, J, IJ)) return 1;

	for(int n=0; n<m_numnodes;n++){
		for(int id=0; id<3; id++){
			double tmp=0;
			for(int jd=0; jd<3; jd++)
				tmp+=IJ[id][jd] * m_deriv[n][jd];
				m_cartp[n][id]=tmp;
		}
	}

	for(int n=0; n<m_numnodes; n++) {
		rhs[n]=- (m_cartp[n][0] * Ji[0] +
		          m_cartp[n][1] * Ji[1] +
		          m_cartp[n][2] * Ji[2]);
	}

	return 0;
}


int
FEShape::MagPri(const point *pts, const point sens, double mul, point *mat)
{
	int n;
	
	double r3,rx,ry,rz;
	double Detjb, volsig;
	double J[3][3];
	
	/* clear matrix */
	memset(mat, 0, 3 * sizeof(point));
	
	int ngp=numVolGP();

        for (int g=0; g<ngp; g++) {
                double gi, gj, gk, gwt;
                getVolGP(g, gi, gj, gk, gwt);
		/* for each gaussian point */
		Shape(gi, gj, gk);
		if(Jacob(pts, &Detjb, J, NULL))
			return 1;
		/* global coordinates of the gauss point */
		rx=sens[0];
		ry=sens[1];
		rz=sens[2];
		for(n=0; n<m_numnodes; n++){
			rx-=m_shape[n]*pts[n][0];
			ry-=m_shape[n]*pts[n][1];
			rz-=m_shape[n]*pts[n][2];
		}
		/* now we have (R/R^3) */
		
		r3=sqrt(rx*rx+ry*ry+rz*rz);
		r3=r3*r3*r3;
		
		volsig=mul*Detjb*gwt/r3;
		
		rx*=volsig;
		ry*=volsig;
		rz*=volsig;
		
		/* obtain the matrix */
//		mat[0][0]+=0;
		mat[0][1]+=rz;
		mat[0][2]-=ry;
		
		mat[1][0]-=rz;
//		mat[1][1]+=0;
		mat[1][2]+=rx;
		
		mat[2][0]+=ry;
		mat[2][1]-=rx;
//		mat[2][2]+=0;
	}
	
	return 0;
	
}

int
FEShape::MagSec(const point *pts, const point sens, double mul, point *mat)
{
	if (mat == NULL)
		return 1;

	memset(mat, 0, sizeof(point) * m_numnodes);
	
	for(int n=0; n<m_numfaces; n++)
		if(MagSecSurf(n, pts, sens, mul, mat)) return 1;
	
	return 0;
	
}

int
FEShape::MagSec(const point *pts, const point sens, double mul,
		double *mx, double *my, double *mz)
{
	if (mx == NULL || my == NULL || mz == NULL)
		return 1;

	memset(mx, 0, sizeof(double) * m_numnodes);
	memset(my, 0, sizeof(double) * m_numnodes);
	memset(mz, 0, sizeof(double) * m_numnodes);
	
	for(int n=0; n<m_numfaces; n++)
		if(MagSecSurf(n, pts, sens, mul, mx, my, mz)) return 1;
	
	return 0;
	
}

int
FEShape::MagSecSurf(int surf, const point *pts, const point sens,
		    double mul, point *mat)
{
	int n;
	double bx, by, bz;
	double ds;
/*  const char *p;*/
	
	int ngp=numSurfGP();
	
        for (int g=0; g<ngp; g++) {
                double gi, gj, gk, gwt;
                getSurfGP(surf, g, gi, gj, gk, gwt);
		Shape(gi, gj, gk);
      		if(SurfIntR(pts, sens, surf, bx, by, bz)) return 1;
      		ds=mul*gwt;
      		for(n=0; n<m_numnodes; n++){
        		mat[n][0]+=bx*m_shape[n]*ds;
        		mat[n][1]+=by*m_shape[n]*ds;
        		mat[n][2]+=bz*m_shape[n]*ds;
      		}
	}
	return 0;
}

int
FEShape::MagSecSurf(int surf, const point *pts, const point sens,
		    double mul, double *mx, double *my, double *mz)
{
	int n;
	double bx, by, bz;
	double ds;
/*  const char *p;*/
	
	int ngp=numSurfGP();
	
        for (int g=0; g<ngp; g++) {
                double gi, gj, gk, gwt;
                getSurfGP(surf, g, gi, gj, gk, gwt);
		Shape(gi, gj, gk);
      		if(SurfIntR(pts, sens, surf, bx, by, bz)) return 1;
      		ds=mul*gwt;
      		for(n=0; n<m_numnodes; n++){
        		mx[n]+=bx*m_shape[n]*ds;
        		my[n]+=by*m_shape[n]*ds;
        		mz[n]+=bz*m_shape[n]*ds;
      		}
	}
	return 0;
}



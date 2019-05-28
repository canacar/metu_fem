/* $Id: shape.h,v 1.2 2008/01/28 07:36:28 canacar Exp $ */
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
#ifndef _SHAPE_H_
#define _SHAPE_H_

#include "define.h"
#include "dmatrix.h"


// callback function prototype first three are the global coordinates
// for which the evalluation is required. The 4th argument is the
// user defined pointer passed to the BoundCB function, The point
// is returned by the final argument. The callback function should
// return 0 for success.

typedef int (*cbVec)(double , double , double, void *, point &);



// abstract base class for FEM element shapes
class FEShape {
public:
	FEShape(int nnodes, int nfaces, int nnodeface);
	virtual ~FEShape();

	inline int numNodes(void) const {return m_numnodes;}
	inline int numFaces(void) const {return m_numfaces;}
	inline int nodePerFace(void) const {return m_nodeface;}

	void Local2Global(const point *pts, double &x,
			  double &y, double &z);
	double interp(const double *fld, double x,
			  double y, double z);
	int grad(const point *pts, const double *fld, point &gr,
		 double x, double y, double z);
	int jacob(const point *pts, double x, double y, double z,
		  double J[3][3], double IJ[3][3]);
	int Matrix(const point *pts, double sigma,
		   DFSMatrix &mat, double *vol);
	int Matrix(const point *pts, double *sarray,
		   DFSMatrix &mat, double *vol);
	int Matrix(const point *pts, double (*getgp)(double,double,double),
		   DFSMatrix &mat, double *vol);
	int Volume(const point *pts, double &vol);
	int SurfaceArea(const point *pts, int surf, double &sa);

	int Bound(const point Ji, const point *pts, double *rhs);
	int BoundS(const point Ji, const point *pts, double *rhs);
	int BoundYan(const point D, const point Ji,
		     const point *pts, double *rhs);
	int BoundCB(const cbVec Vec, const point *pts,
		    double *rhs, void *arg);
	int BoundSCB(const cbVec Vec, const point *pts,
		    double *rhs, void *arg);
/*
	int BoundDipole(const point D, const point Rl, const point *pts,
			double *rhs);
	int BoundDipoleVol(const point D, double x, double y, double z,
			   double rhs);
	int BoundDipoleSurf(const point D, double x, double y, double z,
                      int it, int jt, const point *pts, double *rhs);
*/
	int MagPri(const point *pts, const point sens, double mul,
		   point *mat);
	int MagSec(const point *pts, const point sens, double mul,
		   point *mat);
	int MagSecSurf(int surf, const point *pts, const point sens,
		       double mul, point *mat);
	int MagSec(const point *pts, const point sens, double mul,
		   double *mx, double *my, double *mz);
	int MagSecSurf(int surf, const point *pts, const point sens,
		       double mul, double *mx, double *my, double *mz);

protected:

	struct GP3D {
		double i, j, k;
		double wt;
	};

	int Jacob(const point *pts,  double *detjb,
		  double jacob[3][3], double ijacob[3][3]);
	int BoundSurf(const point Ji, const point *pts,
		      int surf, double *rhs);

	int BoundSurfCB(const cbVec Vec, const point *pts,
		      int surf, double *rhs, void *arg);

	virtual void Shape(double x, double y, double z) = 0;

	virtual int SurfInt(const point *pts, int surf,
		    double &vnx, double &vny, double &vnz) = 0;

	virtual int SurfIntR(const point *pts, const point sens, int surf,
		     double &bx, double &by, double &bz) = 0;

	inline int numVolGP(void) {return m_numvgp;}
	inline int numSurfGP(void) {return m_numsgp;}

	inline void getVolGP(int gp, double &i, double &j,
			      double &k, double &wt) {
		assert (gp >=0 && gp < m_numvgp);
		i=m_vgp[gp].i; j=m_vgp[gp].j; k=m_vgp[gp].k;
		wt=m_vgp[gp].wt;
	}

	inline void getSurfGP(int surf, int gp, double &i,
			       double &j, double &k, double &wt) {
		assert (gp >=0 && gp < m_numsgp);
		assert (surf >=0 && surf <m_numfaces);
		int ind=surf*m_numsgp+gp;
		i=m_sgp[ind].i; j=m_sgp[ind].j; k=m_sgp[ind].k;
		wt=m_sgp[ind].wt;
	}

	inline int SurfNode(int surf, int seq) {
		assert (surf >=0 && surf <m_numfaces);
		assert (seq >=0 && seq < m_nodeface);
		return m_surfnodes[surf*m_nodeface+seq];
	}

// default, linear gauss points
	static const double m_gpos[3];
        static const double m_gval[3];

// volume and surface GP derived from above
	int m_numvgp, m_numsgp;
	GP3D *m_vgp, *m_sgp;

	int *m_surfnodes;

	double *m_shape;
	point *m_deriv, *m_cartp;

	int m_numnodes;
	int m_numfaces;
	int m_nodeface;
};

#endif


/* $Id: engine.h,v 1.2 2008/01/28 07:36:28 canacar Exp $ */
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
//---------------------------------------------------------------------------
#ifndef engineH
#define engineH
//---------------------------------------------------------------------------
#include <zlib.h>

#include "define.h"
#include "hptimer.h"
#include "petscksp.h"
#include "femmesh3.h"
#include "meshutil.h"

#include <list>

typedef list<DInfo> DIList;

class RHSVec {
 public:

	~RHSVec();

	inline int numNodes(void)
		{ return m_nnodes; }

	inline int numDipoles(void)
		{ return m_dipoles.size(); }

	inline int clear(void) {
		int ierr = VecSet(m_vec, 0);
		CHKERRQ(ierr);
		m_dipoles.clear();
		return 0;
	}

	inline void addDipole(const DInfo &di)
		{ m_dipoles.push_back(di); }

	inline DIList::iterator begin(void)
		{ return m_dipoles.begin(); }

	inline DIList::iterator end()
		{ return m_dipoles.end(); }

 protected:
	friend class FEngine;

	// only FEngine can construct
	RHSVec(MPI_Comm comm, int gsize, int lsize = PETSC_DECIDE);

	//  only for use by FEngine and descendants
	inline Vec *rhs(void)
		{ return &m_vec; }

	int m_nnodes;
	MPI_Comm m_comm;

	Vec m_vec;
	DIList m_dipoles;
};

class EShellMat {
public:
	EShellMat(MPI_Comm comm, int nnodes, int nloc, int mloc);
	~EShellMat();

	static int Mult(Mat, Vec, Vec);

	DFSMatrix *m_emat;
	int m_nnodes, m_err;
	int *m_idx;
	double *m_tmp;
	Mat m_mat;
};

class FEngine {
 public:
	FEngine (MPI_Comm comm, FEMesh &msh);
	virtual ~FEngine();
	
	int setupSolver(void);

	inline int clearRHS(void) {
		return m_rhs->clear();
	}

	inline int numNodes(void) { return m_nNodes; }
	inline int numElements(void) { return m_nElements; }

	inline int numProc(void) { return m_size; }
	inline int myRank(void) { return m_rank; }

//	inline int pMesh2Mat(int &p) { return m_mesh->pMesh2Mat(p); }

	inline double elemVolume(int elem) {
		double vol;
		m_mesh->elementVolume(elem, vol);
		return vol;
	}

	inline int mesh2Mat(int i) {
		assert (i >= 0 && i < m_nNodes);
		return m_ordering[i];
	}

	inline int mat2Mesh(int i) {
		assert (i >= 0 && i < m_nNodes);
		return m_revorder[i];
	}

	inline void mesh2Mat(int ni, int *iv) {
		assert(iv != NULL);
		for (int n = 0; n < ni; n++)
			iv[n] = mesh2Mat(iv[n]);
	}

	inline void mat2Mesh(int ni, int *iv) {
		assert(iv != NULL);
		for (int n = 0; n < ni; n++)
			iv[n] = mat2Mesh(iv[n]);
	}

	inline void mesh2Mat(int ni, const int *iv, int *dst) {
		assert(iv != NULL && dst != NULL);
		for (int n = 0; n < ni; n++)
			dst[n] = mesh2Mat(iv[n]);
	}

	inline void mat2Mesh(int ni, const int *iv, int *dst) {
		assert(iv != NULL && dst != NULL);
		for (int n = 0; n < ni; n++)
			dst[n] = mat2Mesh(iv[n]);
	}

	inline int addSource(const DInfo &di) {
		if (di.smodel == SM_J)
			return addSourceJ(di);
		else if (di.smodel == SM_YAN)
			return addSourceYan(di);
		else return -1;
	}

	int solvePot(void);
	int solvePot(Vec x);
	int solvePot(Vec *C, int numvec);

	int getVectors(int nvec, Vec *vec[]);
	int clearVectors(int nvec, Vec vec[]);

	double *solveMag(int numsens, const point *sens);
	double *solveMag(const Vec *Cmat, int numsens, const point *sens);

	// this the seperated version, allows use of different phi
	// in secondary computation
	double *solveMagPri(int numsens, const point *sens);
	double *solveMagSec(const Vec *Cmat, int numsens, Vec *phi);

	int findNode(double x, double y, double z);

	int invertSensCols(int numsens, const int *sens, Vec *vec[]);
	int invertSensCols(int numsens, const point *sens, Vec *vec[]);
	int calcMagSecMat(int numsens, const point *sens, Vec *vec[]);
	int calcMagSecMat(int numsens, const point *sens, Vec vec[]);

	int saveSensMatDip(const char *fnbase, const Vec *Ainv, int nAinv);
	int saveMagSensMatDip(const char *fnbase, const Vec *Cmat,
			      int ns, const point *sens);

	inline int saveRHS(const char *fname) {
		return saveVector(fname, *m_rhs->rhs());
	}

	inline int savePot(const char *fname) {
		return saveVector(fname, m_phi);
	}

	int saveVector(const char *fn, Vec vec);

	int saveMag(const char *fname, const double *mag, int numsens);
	FILE *createDipoleInfo(const char *, int);
	int writeDipoleInfo(FILE *fi, int num, const DInfo *dip);

	int computeMagLead (int numsens, const point *sloc,
			    const point *sdir, Vec *vec[], double rad);

 protected:
	void writeInfoHeader(FILE * fi, int nd);
	int fillMatrix(void);
	int findNodeLocal(double x, double y, double z);

	int addSourceJ(const DInfo &di);
	int addSourceYan(const DInfo &di);

	int addMagPriJ(DInfo &di, const point sens, double Bp[3]);
	int addMagSecondary(const point sens, const Vec phi, Vec *dm,
			    double Bs[3]);

	int sensMatrix(gzFile fd, const Vec *Ainv, int nAinv);
	int sensMatrixMag(gzFile fd, const Vec *Cmat,
			  int ns, const point *sens);

	int assignProcessor(const int *vpart, int *vfreq, int *vcount);
	int PartitionElements(void);
	int PartitionElementsSimple(void);
	int PartitionVertices(void);

	static int computeVecPotCB(double x, double y, double z,
				   void *arg, point &out);

 private:
	FEMesh *m_mesh;		// mesh
	RHSVec *m_rhs;		// RHS vector
	Vec	m_phi;		// Solution vector

	double *m_weights;	// Temporary storage for shape function weights

	// timers - XXX replace with petsc variants
	HPTimer *t_create, *t_fill, *t_sync, *t_fact;
	HPTimer *t_solve, *t_sinv, *t_diff;
	static HPTimer *t_cback;

	// petsc variables
	Mat         m_A;        // linear system matrix
	KSP         m_ksp;     // linear solver context

	int	   *m_ordering;	// ordering vector derived from AO
	int	   *m_revorder;	// reverse ordering vector derived from AO
	
	MPI_Comm m_comm;	// MPI communicator
	int m_rank, m_size;	// processor rank and # of processors

	// element and vector partitions
	// (vector partition is in Petsc Order)
	int m_vstart, m_vend, m_vsize;

	int m_nNodes, m_nElements;
	int m_zpmesh, m_zpmat;

	int m_infoseq;		// dipole info sequence number
};

#endif


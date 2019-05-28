/* $Id: engine.cc,v 1.5 2008/11/24 04:24:47 canacar Exp $ */
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

#include <math.h>
#ifdef _WIN32_
#include <conio.h>
#else
#include <signal.h>
#include <stdlib.h>
#endif

#include "engine.h"

#include "petscksp.h"
#include "petscao.h"
#include "femmesh3.h"
#include "hptimer.h"
#include "errno.h"

#ifndef _WIN32_
extern "C" {
#ifdef LINUX_ASSERT
	extern void __assert_fail (const char *ex, const char *file,
				   unsigned int line, const char *func)
	{
		mprintf("Assertion failed! %s:%d:%s: %s\n",
			file, line, func, ex);
		abort();
	}
#else  
	void __assert (const char *file, int line, const char *ex)
	{
		mprintf("Assertion failed! %s:%d %s\n", file, line, ex);
		abort();
	}
#endif
}
#endif


//extern int errno;


//---------------------------------------------------------------------------
RHSVec::RHSVec(MPI_Comm comm, int gsize, int lsize)
{
	m_comm = comm;
	m_nnodes = gsize;
	int ierr =  VecCreateMPI(m_comm, lsize, gsize, &m_vec);
	CHKERRCONTINUE(ierr);
	if (ierr)
		throw(ierr);
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
RHSVec::~RHSVec()
{
	int ierr = VecDestroy(m_vec);
	CHKERRCONTINUE(ierr);
	if (ierr)
		throw(ierr);
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// source info must be in m_rhs, solution must be in m_phi
// mag vector is reduced to processor 0 (others are incomplete)
double*
FEngine::solveMag(int numsens, const point *sens)
{
	int ierr;

	if (numsens < 1 || sens == NULL)
		return NULL;

	double *mag = new double[numsens * 6];

	double B[6];
	double *Br = mag;

	Vec *dm;
	ierr = VecDuplicateVecs(m_phi , 3, &dm);
	CHKERRCONTINUE(ierr);
	if (ierr)
		return NULL;

	for (int s = 0; s < numsens; s++) {
		memset(B, 0, sizeof(B));

		// primary field
		for (DIList::iterator i = m_rhs->begin();
		     i != m_rhs->end(); i++) {
			DInfo &di = *i;
			if (di.smodel == SM_J) {
				if (addMagPriJ(di, sens[s], B))
					mprintf("Error calculating PriMagJ\n");
			} else
				mprintf("Only J source model is supported\n");
		}

		// secondary field
		if (addMagSecondary(sens[s], *m_rhs->rhs(), dm, B + 3))
			mprintf("Error calculating sec mag field\n");

		// sync
		MPI_Reduce(B, Br, 6, MPI_DOUBLE, MPI_SUM, 0, m_comm);
		Br += 6;
	}

	ierr = VecDestroyVecs(dm, 3);
	CHKERRCONTINUE(ierr);
	if (ierr)
		return NULL;

	return mag;
}
//---------------------------------------------------------------------------
// source info and potential field must be in RHS Cmat is a Vec array
// containing 3 * numsens vectors, the result vector m_phi must
// contain the solution for m_rhs
// mag vector is reduced to processor 0 (others are incomplete)
double*
FEngine::solveMag(const Vec *Cmat, int numsens, const point *sens)
{
	if (Cmat == NULL || numsens < 1 || sens == NULL)
		return NULL;

	double *mag = new double[numsens * 6];

	const Vec *cmat = Cmat;

	double Bp[3];
	double *B = mag;

	for (int s = 0; s < numsens; s++) {
		memset(Bp, 0, 3 * sizeof(double));
		memset(B, 0, 6 * sizeof(double));
		
		// primary field
		for (DIList::iterator i = m_rhs->begin();
		     i != m_rhs->end(); i++) {
			DInfo &di = *i;
			if (di.smodel == SM_J) {
				if (addMagPriJ(di, sens[s], Bp))
					mprintf("Error calculating PriMagJ\n");
			} else
				mprintf("Only J source model is supported\n");
		}

		MPI_Reduce(Bp, B, 3, MPI_DOUBLE, MPI_SUM, 0, m_comm);

		// secondary field
		VecMDot(m_phi, 3, cmat, B+3);
		cmat += 3;
		B+=6;
	}

	mprintf("Done\n");
	return mag;
}
//---------------------------------------------------------------------------
// source info and potential field must be in RHS solves only
// the primary magnetic field containing 3 * numsens vectors
// mag vector is reduced to processor 0 (others are incomplete)
double*
FEngine::solveMagPri(int numsens, const point *sens)
{
	if (numsens < 1 || sens == NULL)
		return NULL;

	double *mag = new double[numsens * 3];

	double Bp[3];
	double *B = mag;

	for (int s = 0; s < numsens; s++) {
		memset(Bp, 0, 3 * sizeof(double));
		
		for (DIList::iterator i = m_rhs->begin();
		     i != m_rhs->end(); i++) {
			DInfo &di = *i;
			if (di.smodel == SM_J) {
				if (addMagPriJ(di, sens[s], Bp))
					mprintf("Error calculating PriMagJ\n");
			} else
				mprintf("Only J source model is supported\n");
		}

		MPI_Reduce(Bp, B, 3, MPI_DOUBLE, MPI_SUM, 0, m_comm);
		B += 3;
	}

	return mag;
}
//---------------------------------------------------------------------------
// Cmat is a Vec array containing 3 * numsens vectors,
// the result vector is in phi
// mag vector is reduced to processor 0 (others are incomplete)
double*
FEngine::solveMagSec(const Vec *Cmat, int numsens, Vec *phi)
{
	if (Cmat == NULL || numsens < 1 || phi == NULL)
		return NULL;

	int nvec = numsens * 3;

	double *mag = new double[nvec];

	memset(mag, 0, nvec * sizeof(double));
		
	VecMDot(*phi, nvec, Cmat, mag);

	return mag;
}
//---------------------------------------------------------------------------
// tmp MUST be an array containing three initialized vectors 
int
FEngine::addMagSecondary(const point sens, Vec phi, Vec *dm, double Bs[3])
{
	int ierr;
	
	if(sens==NULL || phi==NULL || Bs==NULL) return 1;
	
	Bs[0]=Bs[1]=Bs[2]=0;

	FEShape *shape=m_mesh->getShape();
	int nnodes=shape->numNodes();
	
	point parray[NUM_NODES];	// XXX
	double dx[nnodes];
	double dy[nnodes];
	double dz[nnodes];
	int idx[nnodes];
	
	int ne = numElements();;
	for (int n = 0; n < 3; n++) {
		ierr = VecSet(dm[n], 0);
		CHKERRQ(ierr);
	}

	for(int n=0; n<ne; n++){
		if (m_mesh->getElemClass(n) != m_rank)
			continue;
		
		mesh2Mat(nnodes, *m_mesh->getElem(n), idx);
		
		if(m_mesh->pArray(n, parray))
			return 1;

		if(shape->MagSec(parray, sens,
				 m_mesh->getSigmaE(n)*M_MU04PI, dx, dy, dz)){
			mprintf("Error calculating"
				"sec. mag. field entry %d!\n",n);
			return 1;
		}

		ierr = VecSetValues(dm[0], nnodes, idx, dx, ADD_VALUES);
		CHKERRQ(ierr);
		ierr = VecSetValues(dm[1], nnodes, idx, dy, ADD_VALUES);
		CHKERRQ(ierr);
		ierr = VecSetValues(dm[2], nnodes, idx, dz, ADD_VALUES);
		CHKERRQ(ierr);
	}

	for (int n = 0; n < 3; n++) {
		ierr = VecAssemblyBegin(dm[n]);
		CHKERRQ(ierr);
	}

	for (int n = 0; n < 3; n++) {
		double dot;

		ierr = VecAssemblyEnd(dm[n]);
		CHKERRQ(ierr);
		VecDot(m_phi, dm[n], &dot);
		Bs[n] += dot;
	}

	return 0;
}
//---------------------------------------------------------------------------
int
FEngine::getVectors(int nvec, Vec *vec[])
{
	if(nvec  < 1 || vec == NULL)
		SETERRQ(PETSC_ERR_ARG_OUTOFRANGE, "Invalid argument\n");

	int ierr = VecDuplicateVecs(m_phi, nvec, vec);
	CHKERRQ(ierr);

	return 0;
}
//---------------------------------------------------------------------------
int
FEngine::clearVectors(int nvec, Vec vec[])
{
	if(nvec  < 1 || vec == NULL)
		SETERRQ(PETSC_ERR_ARG_OUTOFRANGE, "Invalid argument\n");

	for (int n = 0; n < nvec; n++) {
		int ierr = VecSet(vec[n], 0);
		CHKERRQ(ierr);
	}

	return 0;
}
//---------------------------------------------------------------------------
int
FEngine::calcMagSecMat(int numsens, const point *sens, Vec *vec[])
{
	if(sens == NULL || numsens < 1)
		SETERRQ(PETSC_ERR_ARG_OUTOFRANGE, "Invalid argument\n");

	int ierr = VecDuplicateVecs(m_phi, 3 * numsens, vec);
	CHKERRQ(ierr);

	return calcMagSecMat(numsens, sens, *vec);
}
//---------------------------------------------------------------------------
int
FEngine::calcMagSecMat(int numsens, const point *sens, Vec *Cmat)
{
	int ierr;

	int nvec = 3 * numsens;
      
	if(sens == NULL || numsens < 1)
		SETERRQ(PETSC_ERR_ARG_OUTOFRANGE, "Invalid argument\n");

	mprintf("Computing C matrix with %d columns\n", nvec);

	FEShape *shape=m_mesh->getShape();
	int nnodes=shape->numNodes();
	
	point parray[NUM_NODES];

	int idx[nnodes];
	double dx[nnodes];
	double dy[nnodes];
	double dz[nnodes];

	int ne = numElements();

	int si = 0;
	for (int s = 0; s<numsens; s++) {
		for(int n=0; n<ne; n++){
			if (m_mesh->getElemClass(n) != m_rank)
				continue;
			
			mesh2Mat(nnodes, *m_mesh->getElem(n), idx);
			
			m_mesh->pArray(n, parray);

			if (shape->MagSec(parray, sens[s], m_mesh->getSigmaE(n)
					  * M_MU04PI, dx, dy, dz)){
				SETERRQ1(PETSC_ERR_LIB, "Error calculating"
					"sec. mag. field entry %d!\n",n);
			}

			ierr = VecSetValues(Cmat[si], nnodes, idx,
					    dx, ADD_VALUES);
			CHKERRQ(ierr);
			
			ierr = VecSetValues(Cmat[si+1], nnodes, idx,
					    dy, ADD_VALUES);
			CHKERRQ(ierr);

			ierr = VecSetValues(Cmat[si+2], nnodes, idx,
					    dz, ADD_VALUES);
			CHKERRQ(ierr);
		}
		// start assembly as soon as the vectors are complete
		for (int n = 0; n < 3; n++) {
			ierr = VecAssemblyBegin(Cmat[si++]);
			CHKERRQ(ierr);
		}
	}

	assert(si == nvec);

	// end assembly for all vectors
	for (int n = 0; n < nvec; n++) {
		ierr = VecAssemblyEnd(Cmat[n]);
		CHKERRQ(ierr);
	}

	return 0;
}
//---------------------------------------------------------------------------
/* calculate primary for a given dipole group (idx amd Ji) */
int
FEngine::addMagPriJ(DInfo &di, const point sens, double Bp[3])
{
	if (di.smodel != SM_J)
		return 1;

	if(Bp == NULL || sens == NULL)
		return 1;

	point parray[NUM_NODES];
	double Ck[3][3];

	int e = di.elem;
	if (m_mesh->getElemClass(e) != m_rank)
		return 0;

	FEShape *shape=m_mesh->getShape();
	
//	Bp[0]=Bp[1]=Bp[2]=0;

	double per=1/elemVolume(e);
	point J;
	
	J[0]=per*di.J[0];
	J[1]=per*di.J[1];
	J[2]=per*di.J[2];

	if(m_mesh->pArray(e, parray))
		return -1;
	
	// calculate primary magnetic field
	if(shape->MagPri(parray, sens, M_MU04PI, Ck))
		return -1;

	
	for(int i=0; i<3; i++)
		for(int j=0; j<3; j++)
			Bp[i]+=Ck[i][j]*J[j];

	return 0;
}
//---------------------------------------------------------------------------
int
FEngine::saveVector(const char *fn, const Vec vec)
{
	int ierr;

	if (fn == NULL)
		return 1;

	Vec vr;
	IS is1, is2;
	VecScatter vs;

	ierr = VecDuplicate(vec, &vr);
	CHKERRQ(ierr);
	ierr = ISCreateGeneral(m_comm, m_nNodes, m_ordering, &is1);
	CHKERRQ(ierr);
	ierr = ISCreateStride(m_comm, m_nNodes, 0, 1, &is2);
	CHKERRQ(ierr);

	ierr = VecScatterCreate(vec, is2, vr, is1, &vs);
	CHKERRQ(ierr);

	ierr = VecScatterBegin(vs, vec, vr, INSERT_VALUES, SCATTER_FORWARD);
	CHKERRQ(ierr);
	ierr = VecScatterEnd(vs, vec, vr, INSERT_VALUES, SCATTER_FORWARD);
	CHKERRQ(ierr);
	ierr = VecScatterDestroy(vs);
	CHKERRQ(ierr);

	PetscViewer viewer;

	ierr = PetscViewerASCIIOpen(m_comm, fn, &viewer);
	CHKERRQ(ierr);

	ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_COMMON);
	CHKERRQ(ierr);

	ierr = VecView(vr, viewer);
	CHKERRQ(ierr);

	PetscViewerDestroy(viewer);
	CHKERRQ(ierr);

	ierr = VecDestroy(vr);
	CHKERRQ(ierr);
	ierr = ISDestroy(is1);
	CHKERRQ(ierr);
	ierr = ISDestroy(is2);
	CHKERRQ(ierr);

	return 0;
}
//---------------------------------------------------------------------------
int
FEngine::fillMatrix(void)
{
	int ierr;

	FEShape *shape = m_mesh->getShape();
	int nnodes = shape->numNodes();
	DFSMatrix mat(nnodes);

	// XXX FEMesh::pArray _requires_ NUM_NODES sized array
	point parray[NUM_NODES];

	// other sigma types not supported yet
	double *sigma = m_mesh->getSigmaE();
	
	int idx[nnodes];

	int ne = m_mesh->numElements();

	for (int n = 0; n < ne; n++) {
		if (m_mesh->getElemClass(n) != m_rank)
			continue;
		
		mesh2Mat(nnodes, *m_mesh->getElem(n), idx);
		
		if (m_mesh->pArray(n, parray))
			return 1;

		if (shape->Matrix(parray, sigma[n], mat, NULL)) {
			PetscPrintf(PETSC_COMM_WORLD, "Error calculating "
				"matrix entry for element %d!\n" ,n);
			return 1;
		}

		// if we have zero node in idx, zero the corresponding
		// row and column in mat. The '1' is added in setupSolver
		int nz;
		for (nz = 0; nz <nnodes; nz++)
			if (idx[nz] == m_zpmat)
				break;

		if (nz < nnodes) {
			for (int i = 0; i < nnodes; i++) {
				mat.set(nz, i, 0);
				mat.set(i, nz, 0);
			}
		}


		// fill in petsc way
		ierr = MatSetValues(m_A, nnodes, idx, nnodes, idx,
				    mat._mat, ADD_VALUES);
		CHKERRQ(ierr);
	}

	// Assembled later

	return 0;
}
//---------------------------------------------------------------------------
int
FEngine::findNodeLocal(double x, double y, double z)
{
//	printf("[%d] searching for zero node in interval [%d %d)\n",
//	       m_rank, m_vstart, m_vend);
	for (int n = m_vstart; n < m_vend; n++) {
		const Point3 &pt = m_mesh->getNode(n);
		if (pt.getX() == x && pt.getY() == y && pt.getZ() == z) {
//			printf("[%d] found zero node %d\n", m_rank, n);
			return n;
		}
	}
//	printf("[%d] no zero node found.\n", m_rank);
	return -1;
}
//---------------------------------------------------------------------------
int
FEngine::setupSolver(void)
{
	int ierr;
	int nn = m_mesh->numNodes();

	t_create->start();
	// create matrix. let petsc decide the storage
	ierr = MatCreateMPIAIJ(m_comm, PETSC_DECIDE, PETSC_DECIDE, nn, nn,
			       81, PETSC_NULL, 80, PETSC_NULL, &m_A);
	CHKERRQ(ierr);

	MatSetOption (m_A, MAT_SYMMETRIC, PETSC_TRUE);
//	ierr = MatSetFromOptions(A);
	CHKERRQ(ierr);

	// get vector range from matrix (decided by PETSC)
	ierr = MatGetOwnershipRange(m_A, &m_vstart, &m_vend);
	CHKERRQ(ierr);
	m_vsize = m_vend - m_vstart;

	// on return, msh -> element_class contains assigned processor
	ierr = PartitionElementsSimple();
	CHKERRQ(ierr);

	ierr = PartitionVertices();
	CHKERRQ(ierr);

	t_create->stop();
		
	t_fill->start();

	// zero node
	m_zpmesh = findNode(0,0,0);
	if(m_zpmesh < 0){
		mprintf("No Zero node, using node 0\n");
		m_zpmesh = 0;
	}else{
		mprintf("Zero: %d\n", m_zpmesh);
	}
	m_zpmat = mesh2Mat(m_zpmesh);

	ierr = fillMatrix();
	CHKERRQ(ierr);

	if (m_rank == 0) {
		printf("[%d] zero node (mat) %d\n", m_rank, m_zpmat);
		ierr = MatSetValue(m_A, m_zpmat, m_zpmat, 1, ADD_VALUES);
		CHKERRQ(ierr);
	}

	// assemble matrix
	ierr = MatAssemblyBegin(m_A, MAT_FINAL_ASSEMBLY);
	CHKERRQ(ierr);
	ierr = MatAssemblyEnd(m_A, MAT_FINAL_ASSEMBLY);
	CHKERRQ(ierr);	
	
	t_fill->stop();

//#define VIEW_AMATRIX	
#ifdef VIEW_AMATRIX
	PetscViewer viewer;

	ierr = PetscViewerASCIIOpen(m_comm, "Amat", &viewer);
	CHKERRQ(ierr);

	ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_COMMON);
	CHKERRQ(ierr);

	ierr = MatView(m_A, viewer);
	CHKERRQ(ierr);

	PetscViewerDestroy(viewer);
	CHKERRQ(ierr);
#endif

	//Create linear solver context

	ierr = KSPCreate(m_comm, &m_ksp);
	CHKERRQ(ierr);
 
	// Set operators. Here the matrix that defines the linear system
	// also serves as the preconditioning matrix. And our A matrix
	// does NOT change between solutions

	ierr = KSPSetOperators(m_ksp, m_A, m_A, DIFFERENT_NONZERO_PATTERN);
	CHKERRQ(ierr);

	//  solver options (cg by default)

	ierr = KSPSetType (m_ksp, KSPCG);
	CHKERRQ(ierr);

	ierr = KSPCGSetType (m_ksp, KSP_CG_SYMMETRIC);
	CHKERRQ(ierr);

	//  preconditioner options
	PC pc;
	ierr = KSPGetPC(m_ksp, &pc);
	CHKERRQ(ierr);
 
//    Set runtime options, e.g.,
//    -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
//    These options will override those specified above as long as
//    KSPSetFromOptions() is called _after_ any other customization
//    routines.

	ierr = KSPSetFromOptions(m_ksp);
	CHKERRQ(ierr);

	return 0;
}
//---------------------------------------------------------------------------
int
FEngine::findNode(double x, double y, double z)
{
	int zp=findNodeLocal(x, y, z);
	int zv[2], zr[2];

	if(zp>=0){
		zv[0]=1;
		zv[1]=zp;
	} else
		zv[0] = zv[1] = 0;

	MPI_Allreduce(zv, zr, 2, MPI_INT, MPI_SUM, m_comm);

	if(zr[0]==0)
		return -1;

	zp= zr[1]/zr[0];

	return zp;
}
//---------------------------------------------------------------------------
FEngine::FEngine(MPI_Comm comm, FEMesh &msh)
{
	int ierr;

	m_mesh = &msh;
	m_comm = comm;
	m_nNodes = msh.numNodes();
	m_nElements = msh.numElements();

	t_create=HPTimerMgr::getTimer("Create Matrix");  
	t_fill=HPTimerMgr::getTimer("Fill Matrix");
	t_sync=HPTimerMgr::getTimer("Synchronize Matrix");
	t_fact=HPTimerMgr::getTimer("Factorize Matrix");

	t_solve=HPTimerMgr::getTimer("Solve Matrix");
	t_sinv=HPTimerMgr::getTimer("Invert Sensor Rows");
	t_diff=HPTimerMgr::getTimer("Compute Sens Column");

	if (t_cback == NULL)
		t_cback=HPTimerMgr::getTimer("Mag Vec Pot callback");

	m_zpmesh = m_zpmat = -1;
	m_rhs = new RHSVec(m_comm, m_nNodes);

	m_ordering = new int[m_nNodes];
	m_revorder = new int[m_nNodes];

	for (int n = 0; n < m_nNodes; n++)
		m_ordering[n] = n;

	for (int n = 0; n < m_nNodes; n++)
		m_revorder[n] = n;

	MPI_Comm_size(m_comm, &m_size);
        MPI_Comm_rank(m_comm, &m_rank);

	ierr = VecDuplicate(*m_rhs->rhs(), &m_phi);
	CHKERRCONTINUE(ierr);
	if (ierr)
		throw(ierr);
}
//---------------------------------------------------------------------------
FEngine::~FEngine()
{
	MatDestroy(m_A);
	KSPDestroy(m_ksp);
	delete m_ordering;
	delete m_revorder;
	delete m_rhs;
}
//---------------------------------------------------------------------------
int
FEngine::addSourceJ(const DInfo &di)
{
	int ierr;

	assert (di.smodel == SM_J);

	int elem=di.elem;

	if (elem < 0 || elem >= numElements())
		return 1;

	if (m_mesh->getElemClass(elem) != m_rank)
		return 0;

	double vol = elemVolume(elem);
	if (vol <= 0)
		return 1;

	printf("[%d] adding J dipole\n", m_rank);

	FEShape *shape=m_mesh->getShape();
	int nnodes=shape->numNodes();

	point parray[NUM_NODES];	// XXX
	double r[nnodes];
	int idx[nnodes];

	point j;

	// negative since shape->Bound also negates
	// XXX fix both
	j[0]=-di.J[0]/vol;
	j[1]=-di.J[1]/vol;
	j[2]=-di.J[2]/vol;

	if(m_mesh->pArray(elem, parray))
		return -1;

	if(shape->Bound(j, parray, r)) return -1;

	mesh2Mat(nnodes, *m_mesh->getElem(elem), idx);

	ierr = VecSetValues(*m_rhs->rhs(), nnodes, idx, r, ADD_VALUES);
	CHKERRQ(ierr);
	
	m_rhs->addDipole(di);

	printf("[%d] adding J dipole - done\n", m_rank);

	return 0;
}
//---------------------------------------------------------------------------
int
FEngine::addSourceYan(const DInfo &di)
{
	int ierr;

	assert (di.smodel == SM_YAN);

	int elem = di.elem;
	point d;
	d[0] = di.D[0];
	d[1] = di.D[1];
	d[2] = di.D[2];

//	mprintf("Dipole in elem: %d, (%g, %g, %g) direction: (%g, %g, %g)",
//		elem, d[0],d[1],d[2], di.J[0], di.J[1], di.J[2]);

	if (d[0] <-1 || d[0] > 1)
		return 1;
	if (d[1] <-1 || d[1] > 1)
		return 1;
	if (d[2] <-1 || d[2] > 1)
		return 1;
	
	if (elem < 0 || elem >= numElements())
		return 1;

	if (m_mesh->getElemClass(elem) != m_rank)
		return 0;

	FEShape *shape=m_mesh->getShape();
	int nnodes=shape->numNodes();

	point parray[NUM_NODES];	// XXX
	double r[nnodes];

	int idx[nnodes];
	mesh2Mat(nnodes, *m_mesh->getElem(elem), idx);

	if(m_mesh->pArray(elem, parray))
		return -1;

	// negative since shape->Bound also negates
	// XXX fix both
	point j;
	j[0]=-di.J[0];
	j[1]=-di.J[1];
	j[2]=-di.J[2];

	if(shape->BoundYan(d, j, parray, r))
		return -1;

	ierr = VecSetValues(*m_rhs->rhs(), nnodes, idx, r, ADD_VALUES);
	CHKERRQ(ierr);
	
	m_rhs->addDipole(di);

	return 0;
}
//---------------------------------------------------------------------------
// Assembly is done before solve, b in m_rhs, x in m_phi
int
FEngine::solvePot(void)
{
	int ierr, its;
	Vec *b = m_rhs->rhs();

	t_solve->start();

	ierr = VecAssemblyBegin(*b);
	CHKERRQ(ierr);
	ierr = VecAssemblyEnd(*b);
	CHKERRQ(ierr);

	ierr = KSPSolve(m_ksp, *b, m_phi);
	CHKERRQ(ierr);

	t_solve->stop();

	ierr = KSPGetIterationNumber(m_ksp, &its);
	CHKERRQ(ierr);

	printf("%d iterations\n", its);

#ifdef VERIFY_SOLUTION
	printf("[%d] Verify:\n", m_rank);

	Vec test;
	ierr = VecDuplicate(m_phi, &test);
	CHKERRQ(ierr);
	ierr = MatMult(m_A, m_phi, test);
	CHKERRQ(ierr);

	double nv1, nv2, nvd, nvx;

	ierr = VecNorm(*b, NORM_2, &nv1);
	CHKERRQ(ierr);
	ierr = VecNorm(test, NORM_2, &nv2);
	CHKERRQ(ierr);

	ierr = VecDot(*b, test, &nvx);
	CHKERRQ(ierr);

	PetscScalar mone = -1;
	ierr = VecAXPY(test, mone, *b);
	CHKERRQ(ierr);
	ierr = VecNorm(test, NORM_2, &nvd);
	CHKERRQ(ierr);

	printf("[%d] Norm: nv1:%g nv2:%g nvd:%g nvx:%g\n",
	       m_rank, nv1, nv2, nvd, sqrt(nvx));
#endif

	return 0;
}
//---------------------------------------------------------------------------
// Assembly is done before solve, b in m_rhs, x in m_phi
int
FEngine::solvePot(Vec x)
{
	int ierr, its;
	Vec *b = m_rhs->rhs();

	t_solve->start();

	ierr = VecAssemblyBegin(*b);
	CHKERRQ(ierr);
	ierr = VecAssemblyEnd(*b);
	CHKERRQ(ierr);

	ierr = KSPSolve(m_ksp, *b, x);
	CHKERRQ(ierr);

	t_solve->stop();

	ierr = KSPGetIterationNumber(m_ksp, &its);
	CHKERRQ(ierr);

	printf("%d iterations\n", its);

	return 0;
}
//---------------------------------------------------------------------------
int
FEngine::saveMag(const char *fname, const double *mag, int numsens)
{
	if (fname == NULL || mag == NULL)
		return 1;

	if (numsens < 0)
		return 1;

	if (m_rank != 0)
		return 0;

	FILE *f = fopen(fname, "w");
	if (f == NULL)
		return 1;

	fprintf(f, "%d\n", numsens);

	const double *m = mag;
	for (int n = 0; n < numsens; n++, m += 6) {
		fprintf(f, "%d %g %g %g %g %g %g\n", n+1,
			m[0], m[1], m[2], m[3], m[4], m[5]);
	}

	fclose(f);

	return 0;
}
//---------------------------------------------------------------------------
FILE *
FEngine::createDipoleInfo(const char *fn, int numdip)
{
	FILE *fi = NULL;
	if (m_rank == 0) {
		fi=fopen(fn, "wt");
		if (fi == NULL) {
			mprintf ("info file: %s\n", strerror(errno));
			return NULL;
		}
		writeInfoHeader(fi, numdip);
	}
	return fi;
}
//---------------------------------------------------------------------------
void
FEngine::writeInfoHeader(FILE * fi, int nd)
{
	m_infoseq = 0;
	if (fi != NULL && m_rank == 0) {
		fprintf(fi, "%s\n", m_mesh->getMeshFn());
		fprintf(fi, "%d ref point\n", m_zpmesh);
		fprintf(fi, "%d solutions\n", nd);
	}
}
//---------------------------------------------------------------------------
int
FEngine::writeDipoleInfo(FILE *fi, int num, const DInfo *dip)
{
	if (num < 1)
		SETERRQ(PETSC_ERR_ARG_OUTOFRANGE, "Invalid dipole count\n");

	if (m_rank != 0)
		return 0;

	if (fi == NULL)
		SETERRQ(PETSC_ERR_ARG_OUTOFRANGE, "invalid file pointer\n");
	
	fprintf(fi,"%d 0 %d dipoles\n", m_infoseq++, num);
	
	for (int n = 0; n < num; n++) {
		double x = dip[n].D[0];
		double y = dip[n].D[1];
		double z = dip[n].D[2];

		m_mesh->globalCoord(dip[n].elem, x, y, z);
		
		fprintf(fi,"%g %g %g %g %g %g\n", x, y, z,
				dip[n].J[0], dip[n].J[1], dip[n].J[2]);
	}

	return 0;
}
//---------------------------------------------------------------------------
// slow? (serial code)
int
FEngine::sensMatrix(gzFile fd, const Vec *Ainv, int nAinv)
{
	int ierr;
	if (m_rank == 0) {
		if (fd == NULL)
			SETERRQ(PETSC_ERR_ARG_OUTOFRANGE,
				"Invalid file pointer");
	}

	FEShape *shape = m_mesh->getShape();
	int nnodes=shape->numNodes();

	int ne = numElements();
	int nn = numNodes();

	point parray[NUM_NODES];
	DFSMatrix mat(nnodes);
	int idx[nnodes];

	// this is taken from VecMPIToSEQAll 
	// since all out parallel vectors have the same shape
	// most of the code can be combined
	Vec l_phi, l_Ainv;
	IS         is;
	VecScatter ctx;

	ierr = VecCreateSeq(PETSC_COMM_SELF, nn, &l_phi);
	CHKERRQ(ierr);
	ierr = VecDuplicate(l_phi, &l_Ainv);
	CHKERRQ(ierr);
	ierr = ISCreateStride(PETSC_COMM_SELF,nn, 0, 1, &is);
	CHKERRQ(ierr);
	ierr = VecScatterCreate(m_phi, is, l_phi, is, &ctx);
	CHKERRQ(ierr);

	// scatter Vector l_phi
	ierr = VecScatterBegin(ctx, m_phi, l_phi, INSERT_VALUES,
			       SCATTER_FORWARD);
	CHKERRQ(ierr);
	ierr = VecScatterEnd(ctx, m_phi, l_phi, INSERT_VALUES, 
			     SCATTER_FORWARD);
	CHKERRQ(ierr);

	// this is the vector for l_phi
	double *ax;
	ierr = VecGetArray(l_phi, &ax);
	CHKERRQ(ierr);

	// this is the output row of the sensitivity matrix
	double *smrow = new double[ne];
	double *rsmrow = new double[ne];
	
	for (int row = 0; row < nAinv; row++) {
		// scatter col'th column of Ainv
		ierr = VecScatterBegin(ctx, Ainv[row], l_Ainv, INSERT_VALUES,
				       SCATTER_FORWARD);
		CHKERRQ(ierr);
		ierr = VecScatterEnd(ctx, Ainv[row], l_Ainv, INSERT_VALUES, 
				     SCATTER_FORWARD);
		CHKERRQ(ierr);

		double *az;
		ierr = VecGetArray(l_Ainv, &az);
		CHKERRQ(ierr);

		memset(smrow, 0, ne * sizeof(double));

		for (int n = 0; n < ne; n++) {
			if (m_mesh->getElemClass(n) != m_rank)
				continue;			

			mesh2Mat(nnodes, *m_mesh->getElem(n), idx);
			m_mesh->pArray(n, parray);
		
			if(shape->Matrix(parray, 1, mat, 0)){
				SETERRQ1(PETSC_ERR_LIB,
					 "Error calculating entry %d!\n",n);
			}

			double *emat = mat._mat;
			double ay[nnodes];
			memset(ay, 0, nnodes * sizeof(double));

			for (int i = 0; i < nnodes; i++) {
				double val = ax[idx[i]];
				for (int j = 0; j < nnodes; j++)
					ay[j] += *(emat++) * val;
			}

			for (int i = 0; i < nnodes; i++)
				smrow[n] += ay[i] * az[idx[i]];
		}
		ierr = VecRestoreArray(l_Ainv, &az);
		CHKERRQ(ierr);

		MPI_Reduce(smrow, rsmrow, ne, MPI_DOUBLE, MPI_SUM, 0, m_comm);

		if (m_rank == 0) {
			gzprintf(fd, "%d", row);
			for (int n = 0; n < ne; n++)
				gzprintf(fd, " %g", rsmrow[n]);
			gzprintf(fd, "\n");
		}
	}

	ierr = VecRestoreArray(l_phi, &ax);
	CHKERRQ(ierr);

	ierr = VecScatterDestroy(ctx);
	CHKERRQ(ierr);
	ierr = ISDestroy(is);
	CHKERRQ(ierr);

	ierr = VecDestroy(l_phi);
	CHKERRQ(ierr);
	ierr = VecDestroy(l_Ainv);
	CHKERRQ(ierr);

	delete[] smrow;
	delete[] rsmrow;

	return 0;
}
//---------------------------------------------------------------------------
/* 
   for sensitivity matrix creation.
   sensors must be in MATRIX order
   Note instead of C * Ainv, it computes
   Ainv * Ct and returns the result as transposed matrix
   This holds since Ainv is symmetric
*/
int
FEngine::invertSensCols(int numsens, const int *sens, Vec *vec[])
{
	int ierr;

	if (sens == NULL || numsens <= 0 || numsens > m_nNodes)
		SETERRQ(PETSC_ERR_ARG_OUTOFRANGE, "Invalid arguments");

	ierr = VecDuplicateVecs(m_phi, numsens, vec);
	CHKERRQ(ierr);

	Vec *Cmat = *vec;

	for (int n = 0; n < numsens; n++) {
		int s = mesh2Mat(sens[n]);
		if (s >= m_vstart && s < m_vend) {
//			ierr = VecSetValue(Cmat[n], s, 1,INSERT_VALUES);
			PetscScalar one = 1;
			ierr = VecSetValues(Cmat[n], 1, &s,
					    &one, INSERT_VALUES);
			CHKERRQ(ierr);
		}
		ierr = VecAssemblyBegin(Cmat[n]);
		CHKERRQ(ierr);
	}

	for (int n = 0; n < numsens; n++) {
		ierr = VecAssemblyEnd(Cmat[n]);
		CHKERRQ(ierr);
	}

	ierr = solvePot(Cmat, numsens);
	CHKERRQ(ierr);

	return 0;
}
//---------------------------------------------------------------------------
// callback function that computes the magnetic vector potential at a
// given point. Arg is a pointer to a struct mdip that contains source
// information (magnetic dipole)

void
sph2rect(double &a1, double &a2, double &a3)
{
	double x = a1 * sin (a2) * cos(a3);
	double y = a1 * sin (a2) * sin(a3);
	double z = a1 * cos (a2);

	a1 = x;
	a2 = y;
	a3 = z;
}

void
rect2sph(double &a1, double &a2, double &a3)
{
	double R = sqrt(a1*a1 + a2*a2 + a3*a3);
	double theta = atan2(sqrt(a1*a1 + a2*a2), a3);
	double phi = atan2(a2, a1);

	a1 = R;
	a2 = theta;
	a3 = phi;
}

// rotate vector first around z axis by theta then around y axis by phi
void
align(double theta, double phi, double &x, double &y, double &z)
{
	double t = x;
	x = x  * cos(phi) - y * sin(phi);
	y = t * sin(phi) + y * cos(phi);

	t = z;
	z = z  * cos(theta) - x * sin(theta);
	x = t * sin(theta) + x * cos(theta);
}

// reverse prev. alignment
void
reverse(double theta, double phi, double &x, double &y, double &z)
{
	double t = z;
	z = z  * cos(-theta) - x * sin(-theta);
	x = t * sin(-theta) + x * cos(-theta);

	t = x;
	x = x  * cos(-phi) - y * sin(-phi);
	y = t * sin(-phi) + y * cos(-phi);
}

struct MDipole {
	Point3 c;	// center of the coil
	Point3 n;	// normal of the coil
	double r;	// radius of the coil
	double mul;	// multiplier
};

#define LOOP_DIV 100

HPTimer *FEngine::t_cback = NULL;

int
FEngine::computeVecPotCB(double x, double y, double z, void *arg, point &out)
{
	if (arg == NULL || out == NULL)
		return 1;

	t_cback->start();

	MDipole *mdip = (MDipole *) arg;

	// field point relative to coil
	Point3 r(x, y, z);
	r -= mdip->c;

	// spherical coordinates of the normal
	double nr, nth, nphi;
	mdip->n.getCoord(nr, nth, nphi);
	rect2sph(nr, nth, nphi);
	
	// rotate field point to align coil with z-axis
	align (-nth, -nphi, r.X(), r.Y(), r.Z());

	// align field point with y-axis
	rect2sph(r.X(), r.Y(), r.Z());
	double pa = r.Z() - M_PI_2;
	r.Z() = M_PI_2;
	sph2rect(r.X(), r.Y(), r.Z());
	
	double b = mdip->r;
	double dl = mdip->mul * b * 2 * M_PI / LOOP_DIV;

	double ax = 0;

	for (int n = 0; n < LOOP_DIV; n++) {
		double ang = 2 * M_PI * n / LOOP_DIV;
		double cosa = cos(ang);
		double sina = sin(ang);
		Point3 pc(cosa * b, sina * b, 0);
		pc -= r;
		double R = pc.length();
		ax += sina / R;
	}

	ax *= dl;

	out[0] = ax * cos(pa);
	out[1] = ax * sin(pa);
	out[2] = 0;

	reverse(-nth, -nphi, out[0], out[1], out[2]);

	t_cback->stop();
	return 0;
}
//---------------------------------------------------------------------------
/* 
   Computes the reciprocal field for a given sensor set
   sloc gives the location of sensors. sdir is the direction
   returns an array of numsens vectors of longth num_nodes each
*/
int
FEngine::computeMagLead (int numsens, const point *sloc,
			 const point *sdir, Vec *vec[], double rad)
{
	int ierr;
	MDipole mdip;

	if (sloc == NULL || sdir == NULL || numsens <= 0 || numsens > m_nNodes)
		SETERRQ(PETSC_ERR_ARG_OUTOFRANGE, "Invalid arguments");

	ierr = VecDuplicateVecs(m_phi, numsens, vec);
	CHKERRQ(ierr);
		
	FEShape *shape=m_mesh->getShape();
	int nnodes=shape->numNodes();

	Vec *Cmat = *vec;
	point parray[NUM_NODES];	// XXX
	double r[nnodes];
	int idx[nnodes];

	mdip.r = rad;

	printf("Coil radius: %g\n", mdip.r);

	for (int n = 0; n < numsens; n++) {
		// set mdip parameters
		mdip.c = sloc[n];
		mdip.n = sdir[n];
		mdip.n.normalize();
 
		for (int e = 0; e < numElements(); e++) {
			if (m_mesh->getElemClass(e) != m_rank)
				continue;

			// XXX set mdip mul - is it correct?
			mdip.mul = m_mesh->getSigmaE(e) * M_MU04PI;

			m_mesh->pArray(e, parray);

			if(shape->BoundSCB(computeVecPotCB, parray, r, &mdip))
				SETERRQ(1, "Failed to compute rhs ");
			
			mesh2Mat(nnodes, *m_mesh->getElem(e), idx);

			for (int i = 0; i < nnodes; i++) {
				if (idx[i] == m_zpmat)
					r[i] = 0;
			}

			ierr = VecSetValues(Cmat[n], nnodes,
					    idx, r, ADD_VALUES);
			CHKERRQ(ierr);
		}

		ierr = VecAssemblyBegin(Cmat[n]);
		CHKERRQ(ierr);
		ierr = VecAssemblyEnd(Cmat[n]);
		CHKERRQ(ierr);
	}

	// solve for all
	ierr = solvePot(Cmat, numsens);
	CHKERRQ(ierr);

	return 0;
}
//---------------------------------------------------------------------------
/*
  Compute and return Ainv * C, computing Ainv one at a time, and only
  for the columns of C, replacing C. Assumes C is available on every
  processor. The results are broadcast to each processor.
  Uses m_phi as temporary (XXX: need a different tmp vector?)
*/
int
FEngine::solvePot(Vec *C, int nvec)
{
	int ierr, its = 0;

	if (C == NULL)
		SETERRQ(PETSC_ERR_ARG_OUTOFRANGE, "Invalid Vector array");

	mprintf("Solving for %d columns\n", nvec);

	for (int n = 0; n < nvec; n++) {
		t_sinv->start();
		ierr = VecSwap(C[n], m_phi);
		CHKERRQ(ierr);
		ierr = KSPSolve(m_ksp, m_phi, C[n]);
		CHKERRQ(ierr);
		t_sinv->stop();
		ierr = KSPGetIterationNumber(m_ksp, &its);
		CHKERRQ(ierr);
		mprintf("%d: %d iterations, OK\n", n, its);
	}

	return 0;
}
//---------------------------------------------------------------------------
// requires dipole soultion available in RHS
// XXX save it as binary ?

int
FEngine::saveSensMatDip(const char *fnbase, const Vec *Ainv, int nAinv)
{
	int ierr;

	if (fnbase == 0)
		SETERRQ(PETSC_ERR_ARG_OUTOFRANGE, "invalid base filename");

	int len = strlen(fnbase);

#define LEN_INC 16
	if (len <0 || len + LEN_INC >= PATH_MAX)
		SETERRQ(PETSC_ERR_ARG_OUTOFRANGE, "filename too long\n");

	len +=LEN_INC;

	char buf[len];

	snprintf(buf, len, "%s.dat.gz", fnbase);

	gzFile fd = 0;
	if (m_rank == 0) {
		fd = gzopen(buf, "wt9");
		if (fd == NULL)
			SETERRQ1(PETSC_ERR_LIB,
				 "Failed to open output file %s\n", buf);
	}

	t_diff->start();
	mprintf("Calculating sensitivity matrix rows\n");

	ierr = sensMatrix(fd, Ainv, nAinv);
	CHKERRQ(ierr);

	t_diff->stop();

	if (m_rank == 0) {
		gzclose(fd);
	}

	return 0;
}
//---------------------------------------------------------------------------
// input: Cmat is (C * Ainv)t

int
FEngine::sensMatrixMag(gzFile fd, const Vec *Cmat, int ns, const point *sens)
{
	int ierr;

	FEShape *shape = m_mesh->getShape();
	int nnodes=shape->numNodes();

	point parray[NUM_NODES];
	DFSMatrix mat(nnodes);

	int nn = numNodes();
	int ne = numElements();
	int idx[nnodes];

	double dx[nnodes];
	double dy[nnodes];
	double dz[nnodes];
	
	// this is taken from VecMPIToSEQAll 
	// since all out parallel vectors have the same shape
	// most of the code can be combined
	Vec l_phi, *l_sens;
	IS         is;
	VecScatter ctx;

	ierr = VecCreateSeq(PETSC_COMM_SELF, nn, &l_phi);
	CHKERRQ(ierr);

	ierr = VecDuplicateVecs(l_phi, 3, &l_sens);
	CHKERRQ(ierr);

	ierr = ISCreateStride(PETSC_COMM_SELF,nn, 0, 1, &is);
	CHKERRQ(ierr);
	ierr = VecScatterCreate(m_phi, is, l_phi, is, &ctx);
	CHKERRQ(ierr);

	// scatter Vector l_phi to all processors
	ierr = VecScatterBegin(ctx, m_phi, l_phi, INSERT_VALUES,
			       SCATTER_FORWARD);
	CHKERRQ(ierr);
	ierr = VecScatterEnd(ctx, m_phi, l_phi, INSERT_VALUES, 
			     SCATTER_FORWARD);
	CHKERRQ(ierr);

	// this is the vector for l_phi
	double *ax;
	ierr = VecGetArray(l_phi, &ax);
	CHKERRQ(ierr);

	// this is the output row of the sensitivity matrix
	double *smrow = new double[3*ne];
	double *rsmrow = new double[3*ne];

	for (int s = 0; s < ns; s++) {
		int row = 3*s;
		double **az;

		// scatter 3 rows of C to l_sens corresponding arrays
		for (int i = 0; i < 3; i++) {
			ierr = VecScatterBegin(ctx, Cmat[row + i], l_sens[i],
					       INSERT_VALUES, SCATTER_FORWARD);
			CHKERRQ(ierr);
			ierr = VecScatterEnd(ctx, Cmat[row + i], l_sens[i],
					     INSERT_VALUES, SCATTER_FORWARD);
			CHKERRQ(ierr);
		}

		ierr = VecGetArrays(l_sens, 3,  &az);
		CHKERRQ(ierr);
		memset(smrow, 0, ne * 3 * sizeof(double));

		for (int n = 0; n < ne; n++) {
			if (m_mesh->getElemClass(n) != m_rank)
				continue;			

			mesh2Mat(nnodes, *m_mesh->getElem(n), idx);
			m_mesh->pArray(n, parray);
		
			if(shape->Matrix(parray, 1, mat, 0)){
				SETERRQ1(PETSC_ERR_LIB,
					 "Error calculating entry %d!\n",n);
			}

			double *emat = mat._mat;
			double ay[nnodes];
			memset(ay, 0, nnodes * sizeof(double));

			for (int i = 0; i < nnodes; i++) {
				double val = ax[idx[i]];
				for (int j = 0; j < nnodes; j++)
					ay[j] += *(emat++) * val;
			}

			double *sm = smrow;
			for (int j = 0; j < 3; j++) {
				for (int i = 0; i < nnodes; i++)
					sm[n] += ay[i] * az[j][idx[i]];
				sm+=ne;
			}

			// now dC/dsig * phi
			if(shape->MagSec(parray, sens[s], M_MU04PI,
					 dx, dy, dz))
				return 1;

			sm = smrow;
			for (int i = 0; i < nnodes; i++)
				sm[n] -= dx[i] * az[0][idx[i]];

			sm += ne;
			for (int i = 0; i < nnodes; i++)
				sm[n] -= dy[i] * az[1][idx[i]];

			sm += ne;
			for (int i = 0; i < nnodes; i++)
				sm[n] -= dz[i] * az[2][idx[i]];
		}

		ierr = VecRestoreArrays(l_sens, 3, &az);
		CHKERRQ(ierr);

		MPI_Reduce(smrow, rsmrow, 3 * ne,
			   MPI_DOUBLE, MPI_SUM, 0, m_comm);

		
		if (m_rank == 0) {
			double *r = rsmrow;
			for (int i = 0; i < 3; i++) {
				gzprintf(fd, "%d", row + i);
				for (int n = 0; n < ne; n++)
					gzprintf(fd, " %g", *(r++));
				gzprintf(fd, "\n");
			}
		}
	}

	ierr = VecRestoreArray(l_phi, &ax);
	CHKERRQ(ierr);

	ierr = VecScatterDestroy(ctx);
	CHKERRQ(ierr);
	ierr = ISDestroy(is);
	CHKERRQ(ierr);
	
	ierr = VecDestroy(l_phi);
	CHKERRQ(ierr);

	ierr = VecDestroyVecs(l_sens, 3);
	CHKERRQ(ierr);

	delete[] smrow;
	delete[] rsmrow;
	
	return 0;
}
//---------------------------------------------------------------------------
int
FEngine::saveMagSensMatDip(const char *fnbase, const Vec *Cmat,
			   int ns, const point *sens)
{
	int ierr;
	if (fnbase == 0)
		return 1;

	int len = strlen(fnbase);

#define LEN_INC 16
	if (len <0 || len + LEN_INC >= PATH_MAX)
		return 1;

	len +=LEN_INC;

	char buf[len];
	
	snprintf(buf, len, "%s.dat.gz", fnbase);

	gzFile fd = 0;
	if (m_rank == 0) {
		fd = gzopen(buf, "wt9");
		if (fd == NULL)
			SETERRQ1(PETSC_ERR_LIB,
				 "Failed to open output file %s\n", buf);
	}
	t_diff->start();

	mprintf("Calculating B sens. columns\n");

	ierr = sensMatrixMag(fd, Cmat, ns, sens);
	CHKERRQ(ierr);

	t_diff->stop();

	if (m_rank == 0)
		gzclose(fd);

	return 0;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
int
FEngine::assignProcessor(const int *vpart, int *vfreq, int *vcount)
{
	// try to assign on the wrt frequency order
	// tries nodes with same frequency in rank order
	while (1) {
		int mx = vfreq[0];
		int mid = 0;
		for (int n = 1; n < m_size; n++) {
			if (vfreq[n] > mx) {
				mx = vfreq[n];
				mid = n;
			}
		}
		vfreq[mid] = -1;

		// tried all processors and could not fit!
		if (mx < 0)
			return -1;

		// try to fit to mid
		if (vcount[mid] < vpart[mid]) {
			vcount[mid]++;
			return mid;
		}
	}
}
//---------------------------------------------------------------------------
int
FEngine::PartitionVertices(void)
{
	int nn = numNodes();
	int *vpart, *vcount, *vfreq;
	int ierr;

	// allocate partition, counter and frequency vectors
	ierr = PetscMalloc(m_size * sizeof(int), &vpart);
	CHKERRQ(ierr);
	ierr = PetscMalloc(m_size * sizeof(int), &vcount);
	CHKERRQ(ierr);
	ierr = PetscMalloc(m_size * sizeof(int), &vfreq);
	CHKERRQ(ierr);
	
	assert(m_vsize > 0);

	// use vcount as temporary source vector
	for (int i=0; i<m_size; i++)
		vcount[i] = (m_rank == i) ? m_vsize : 0;

	MPI_Allreduce(vcount, vpart, m_size,
		      MPI_INT, MPI_SUM, m_comm);

	// allocate ordering vectors
	int *aov, *pov;
	ierr = PetscMalloc(m_vsize * sizeof(int), &aov);
	CHKERRQ(ierr);
	ierr = PetscMalloc(m_vsize * sizeof(int), &pov);
	CHKERRQ(ierr);

	for (int n = 0; n < m_size; n++)
		vcount[n] = 0;

	int oidx = 0;
	for (int n = 0; n < nn; n++) {
		int elem[MAX_ELEM_NBR];
		int sz = m_mesh->getNeighborElem(n, elem, MAX_ELEM_NBR);

		memset (vfreq, 0, sizeof(int) * m_size);

		// replace vector with class (processor id)
		// and compute frequency of each processor
		for (int m = 0; m < sz; m++) {
			elem[m] = m_mesh->getElemClass(elem[m]);
			assert (elem[m] >= 0 && elem[m] < m_size);
			vfreq[elem[m]]++;
		}

		// assign vertex depending on frequency, partition size
		// and number of vertices in partitions
		int proc = assignProcessor(vpart, vfreq, vcount);

		if (proc < 0) {
			printf("Failed to assign node %d\n", n);
			return 1;
		}

		// if local, update index vectors
		if (proc == m_rank) {
			pov[oidx] = m_vstart + oidx;
			aov[oidx] = n;
			oidx ++;
		}
	}

	assert(oidx == m_vsize);

	ierr = PetscFree(vpart);
	CHKERRQ(ierr);
	ierr = PetscFree(vcount);
	CHKERRQ(ierr);
	ierr = PetscFree(vfreq);
	CHKERRQ(ierr);
	vpart = vcount = vfreq = NULL;

	// create application ordering

// #define DO_ORDERING
#ifdef DO_ORDERING
	AO ao;

	ierr = AOCreateBasic(m_comm, m_vsize, aov, pov , &ao);
	CHKERRQ(ierr);

	// convert ordering vector
	ierr = AOApplicationToPetsc(ao, m_nNodes, m_ordering);
	CHKERRQ(ierr);

	// convert ordering vector
	ierr = AOPetscToApplication(ao, m_nNodes, m_revorder);
	CHKERRQ(ierr);

	// cleanup
	ierr = AODestroy(ao);
	CHKERRQ(ierr);

#endif
	ierr = PetscFree(aov);
	CHKERRQ(ierr);
	ierr = PetscFree(pov);
	CHKERRQ(ierr);

#ifdef TEST_ORDERING
	// test ordering
	printf("[%d] Testing ordering\n", m_rank);
	for (int n = 0; n < m_nNodes; n++) {
		if (m_ordering[n] < 0 || m_ordering[n] >= m_nNodes ||
		    m_revorder[n] < 0 || m_revorder[n] >= m_nNodes) {
			printf("invalid vector!\n");
			return 1;
		}
		if (n != m_ordering[m_revorder[n]] ||
		    n != m_revorder[m_ordering[n]]) {
			printf("mapping failed\n");
			return 1;
		}
	}
	
	
	Vec vtest;
	
	ierr = VecCreateMPI(m_comm, PETSC_DECIDE, m_nNodes, &vtest);
	CHKERRQ(ierr);

	for (int i = m_vstart; i < m_vend; i++) {
		ierr = VecSetValue(vtest, i, m_ordering[i], INSERT_VALUES);
		CHKERRQ(ierr);
	}
	
	ierr = VecAssemblyBegin(vtest);
	CHKERRQ(ierr);
	ierr = VecAssemblyEnd(vtest);
	CHKERRQ(ierr);
	
	ierr = saveVector("test.vec", vtest);
	CHKERRQ(ierr);
	
	char buf[64];
	sprintf(buf, "order%d.vec", m_rank);
	FILE *f = fopen(buf, "w");
	for (int n = 0; n <m_nNodes; n++)
		fprintf(f, "%d %d %d\n", n, m_ordering[n], m_revorder[n]);
	fclose(f);

	printf("[%d] OK\n", m_rank);
#endif

	return 0;
}
//---------------------------------------------------------------------------
int
FEngine::PartitionElements(void)
{
	int ne, ierr;
	int *elempart;

	Mat		Adj;
	MatPartitioning part;
	IS              isnewproc, istotal; 

	ne = numElements();

	ierr = PetscMalloc(m_size * sizeof(int), &elempart);
	CHKERRQ(ierr);
	
//	PetscPrintf(m_comm,
	printf(
		"%d processors, this is processor %d\n", m_size, m_rank);
	
	for (int i=0; i<m_size; i++)
		elempart[i] = ne/m_size + ((ne % m_size) > i);

	int Estart, Eend, Esize;
	Estart = 0;
	Eend = elempart[0];
	Esize = elempart[m_rank];

	for (int n = 1; n <= m_rank; n++) {
		Estart = Eend;
		Eend += elempart[n];
	}

	ierr = PetscFree(elempart);
	CHKERRQ(ierr);
	elempart = NULL;

//	PetscPrintf(m_comm,
	printf(
		"[%d] Using %d elements from %d to %d\n",
		m_rank, Esize, Estart, Eend);

	int ntotal = 0;

	// count total number of neighbors
	for (int n = Estart; n < Eend; n++) {
		int nnbr = m_mesh->getElemNeighbors(n, NULL, 0);
		if (nnbr < 0) {
			printf("[%d] invalid neighbor count elem %d!\n",
			       m_rank, n);
			return 1;
		}
		// do not include self
		ntotal += nnbr - 1;
	}
	int *ia, *ja;

	ierr = PetscMalloc((Esize + 1) * sizeof(int), &ia);
	CHKERRQ(ierr);
	ierr = PetscMalloc(ntotal * sizeof(int), &ja);
	CHKERRQ(ierr);

	int iidx, jidx;
	iidx = 0;
	jidx = 0;
	for (int n = Estart; n < Eend; n++) {
		ia[iidx++] = jidx;
		int nnbr;
		const int *nbrs = m_mesh->getElemNeighbors(n, nnbr);
		if (nbrs == NULL) {
			printf("[%d] invalid neighbor array elem %d!\n",
			       m_rank, n);
			return 1;
		}
		// do not include self
		for (int m = 0; m < nnbr; m++) {
			if (nbrs[m] == n)
				continue;
			ja[jidx++] = nbrs[m];
		}
	}
	ia[iidx] = jidx;

	if (iidx != Esize || jidx != ntotal) 
		printf("[%d] index mismatch! %d, %d != %d, %d\n",
		       m_rank, iidx, jidx, Esize, ntotal);

	// from dm/ao/examples/tutorials/ex2.c
	// Create the adjacency graph matrix
	ierr = MatCreateMPIAdj(m_comm, Esize, ne, ia, ja,
			       PETSC_NULL, &Adj);
	CHKERRQ(ierr);

	// Create the partioning object
	ierr = MatPartitioningCreate(m_comm, &part);
	CHKERRQ(ierr);
	ierr = MatPartitioningSetAdjacency(part, Adj);
	CHKERRQ(ierr);
	ierr = MatPartitioningSetFromOptions(part);
	CHKERRQ(ierr);

	// isnewproc - indicates for each local element
	//             the new processor it is assigned to
	ierr = MatPartitioningApply(part, &isnewproc);
	CHKERRQ(ierr);
	ierr = MatPartitioningDestroy(part);
	CHKERRQ(ierr);
	
	// isall - gather all IS vectors XXX this wastes space but
	//         much cleaner since we already hold ALL the elements anyway
	//	   store the processor information in the element 'class'.
	ierr = ISAllGather(isnewproc, &istotal);
	CHKERRQ(ierr);

	const PetscInt *idx;
	ierr = ISGetIndices(istotal, &idx);
	CHKERRQ(ierr);

	int elocal = 0;
	for (int n=0; n<ne; n++) {
		m_mesh->setElemClass(n, idx[n]);
		if (idx[n] == m_rank)
			++elocal;
	}
	printf("[%d] %d local elements\n", m_rank, elocal);

	ierr = ISRestoreIndices(istotal, &idx);
	CHKERRQ(ierr);

	ierr = ISDestroy(istotal);
	CHKERRQ(ierr);
	ierr = ISDestroy(isnewproc);
	CHKERRQ(ierr);
	ierr = MatDestroy(Adj);
	CHKERRQ(ierr);
	return 0;
}
//---------------------------------------------------------------------------
int
FEngine::PartitionElementsSimple(void)
{
	int ne, ierr;
	int *elempart;

	ne = m_mesh->numElements();

	ierr = PetscMalloc(m_size * sizeof(int), &elempart);
	CHKERRQ(ierr);
	
//	PetscPrintf(m_comm,
	printf(
		    "%d processors, this is processor %d\n", m_size, m_rank);
	
	for (int i=0; i<m_size; i++)
		elempart[i] = ne/m_size + ((ne % m_size) > i);

	int Estart, Eend, Esize;
	Estart = 0;
	Eend = elempart[0];
	Esize = elempart[m_rank];

	for (int n = 1; n <= m_rank; n++) {
		Estart = Eend;
		Eend += elempart[n];
	}

	for (int n = 1; n <= m_rank; n++) {
		Estart = Eend;
		Eend += elempart[n];
	}

//	PetscPrintf(m_comm,
	printf(
		"[%d] Using %d elements from %d to %d\n",
		m_rank, Esize, Estart, Eend);

	for (int n = 1; n < m_size; n++)
		elempart[n] += elempart[n-1];

	int elocal = 0;
	int n = 0;
	for (int r = 0; r < m_size; r++) {
		for (;n < elempart[r]; n++) {
			m_mesh->setElemClass(n, r);
			if (r == m_rank)
				++elocal;
		}
	}

	PetscFree(elempart);
	elempart = NULL;

	printf("[%d] %d local elements\n", m_rank, elocal);

	return 0;
}
//---------------------------------------------------------------------------
// No longer used (very inefficient) but left as an example
//---------------------------------------------------------------------------
EShellMat::EShellMat(MPI_Comm comm, int nnodes, int mloc, int nloc):
	m_nnodes(nnodes)
{
	m_idx = new int[nnodes];
	m_tmp = new double[nnodes];
	m_emat = new DFSMatrix(nnodes);
	m_err = MatCreateShell(comm, mloc, nloc, PETSC_DETERMINE,
			       PETSC_DETERMINE, this, &m_mat);
	if (m_err)
		return;
	m_err = MatShellSetOperation(m_mat, MATOP_MULT,(void (*) ()) Mult);
}//---------------------------------------------------------------------------
EShellMat::~EShellMat()
{
	MatDestroy(m_mat);
	delete[] m_idx;
	delete[] m_tmp;
	delete m_emat;
}
//---------------------------------------------------------------------------
// Computes the matrix-vector product, y = Ax. 
int
EShellMat::Mult(Mat mat, Vec x, Vec y)
{
	int ierr, nn, xmin, xmax; 
	int *idx;
	double *ax, *ay, *emat;
	EShellMat *me;

	ierr = MatShellGetContext(mat, (void **) &me);
	CHKERRQ(ierr);

	if (me == NULL)
		SETERRQ(PETSC_ERR_ARG_OUTOFRANGE, "Invalid class pointer");

	ierr = VecGetArray(x, &ax);
	CHKERRQ(ierr);

	ierr = VecGetOwnershipRange(x, &xmin, &xmax);
	CHKERRQ(ierr);

	idx  = me->m_idx;
	ay   = me->m_tmp;
	nn   = me->m_nnodes;

#ifdef MULT_SANITY_CHECK
	if (idx == NULL)
		SETERRQ(PETSC_ERR_ARG_CORRUPT, "Invalid index pointer\n");

	if (ay == NULL)
		SETERRQ(PETSC_ERR_ARG_CORRUPT, "Invalid temp pointer\n");

	if (ay == NULL)
		SETERRQ(PETSC_ERR_ARG_CORRUPT, "Invalid vector pointer\n");

	if (nn != 20 && nn != 8)
		SETERRQ(PETSC_ERR_ARG_CORRUPT, "Invalid number of nodes\n");

	if (me->m_emat == NULL)
		SETERRQ(PETSC_ERR_ARG_CORRUPT, "Invalid Matrix\n");
#endif

	emat = me->m_emat->_mat;

#ifdef MULT_SANITY_CHECK
	if (emat == NULL)
		SETERRQ(PETSC_ERR_ARG_CORRUPT, "Invalid array\n");
#endif

 	memset(ay, 0, sizeof(double) * nn);

 	for (int n = 0; n < nn; n++) {
 		if (idx[n] < xmin || idx[n] >= xmax) {
			emat+=nn;
 			continue;
		}
 		double val = ax[idx[n] - xmin];
 		for (int m = 0; m < nn; m++)
 			ay[m] += *(emat++) * val;
 	}

	ierr = VecRestoreArray(x, &ax);
	CHKERRQ(ierr);

	ierr = VecSet(y, 0);
	CHKERRQ(ierr);

 	ierr = VecSetValues(y, nn, idx, ay, ADD_VALUES);
 	CHKERRQ(ierr);

 	ierr = VecAssemblyBegin(y);
 	CHKERRQ(ierr);
 	ierr = VecAssemblyEnd(y);
	CHKERRQ(ierr);

	return 0;
}

//--------------------------------------------------------------------------

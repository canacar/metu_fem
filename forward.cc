/* $Id: forward.cc,v 1.3 2008/11/24 03:49:28 canacar Exp $ */
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
/*
  This code is based on petsc example 'sles/examples/ex2.c'
*/

#include <stdio.h>
#include <stdlib.h>

#include "engine.h"


/* global variables */
int num_sens=0;
int num_bsens = 0;
int num_dipoles=0;

/* only one of sensors and psensors will be used depending on type */
int *sensors = NULL;
point *psensors = NULL;

DInfo *dipoles = NULL;
point *bsensors = NULL;

// read magnetic sensor data -- temporary - should be part of meshutil
point *
ReadMagSensors(char *fn, int &nsens)
{
	FILE *fs=fopen(fn,"rt");
	if (fs == NULL) return 0;
	
	nsens=0;
	point *sens = 0;

	if (MeshUtil::readListHdr(fs, &nsens))
		goto exit;

	if (nsens == 0)
		goto exit;

	sens=new point[nsens];
	if (MeshUtil::readPointList(fs, sens, 1, nsens)) {
		delete[] sens;
		sens = 0;
	}
 exit:
	fclose(fs);
	return sens;
}

#undef __FUNCT__
#define __FUNCT__ "forward"

int
forward(FEngine &engine)
{
	int ierr;

	char fname[PATH_MAX];
	FILE *fi = NULL;

	Vec *Cmat = NULL;

	ierr = engine.setupSolver();
	CHKERRQ(ierr);

	try {
		fi = engine.createDipoleInfo("xinfo.pot", num_dipoles);
	} catch(...) {
		SETERRQ(PETSC_ERR_FILE_OPEN, "Failed to open xInfo file!\n");
	}

	if (bsensors != NULL) {
		ierr = engine.calcMagSecMat(num_bsens, bsensors, &Cmat);
		CHKERRQ(ierr);
	}
	

	for (int n = 0; n< num_dipoles; n++) {
		engine.clearRHS();

		mprintf("Dipole %d\n", n);

		ierr = engine.writeDipoleInfo(fi, 1, &dipoles[n]);
		CHKERRQ(ierr);
		ierr = engine.addSource(dipoles[n]);
		CHKERRQ(ierr);

		ierr = engine.solvePot();
		CHKERRQ(ierr);

		double *mag = NULL;

		if (bsensors != NULL) {
			mag = engine.solveMag(Cmat, num_bsens, bsensors);
			if (mag == NULL)
				SETERRQ(1, "Failed to solve magnetic field\n");
		}
		
		snprintf(fname, sizeof(fname), "x%03d.pot", n);

		ierr = engine.savePot(fname);
		CHKERRQ(ierr);

		if (mag != NULL) {
			snprintf(fname, sizeof(fname), "x%03d.mag", n);
			ierr = engine.saveMag(fname, mag, num_bsens);
			CHKERRQ(ierr);
		}
	}

	if (fi)
		fclose(fi);

	if (bsensors != NULL) {
		ierr = VecDestroyVecs(Cmat, num_bsens);
		CHKERRQ(ierr);
	}
	
	return 0;
}


#undef __FUNCT__
#define __FUNCT__ "forward2"
// this is called if 3 * (num_bsens + 1) > num_dipoles
// in this case storing phi for each dipole instead of
// computing C and computing sensors one by one is better

int
forward2(FEngine &engine)
{
	int ierr;

	char fname[PATH_MAX];
	FILE *fi = NULL;

	Vec *Phimat = NULL;
	Vec *Cmat = NULL;

	if (num_dipoles < 1 || num_bsens < 1 || bsensors == 0)
		SETERRQ(1, "invalid arguments\n");

	ierr = engine.setupSolver();
	CHKERRQ(ierr);

	try {
		fi = engine.createDipoleInfo("xinfo.pot", num_dipoles);
	} catch(...) {
		SETERRQ(PETSC_ERR_FILE_OPEN, "Failed to open xInfo file!\n");
	}

	ierr = engine.getVectors(num_dipoles, &Phimat);
	CHKERRQ(ierr);

	ierr = engine.getVectors(3, &Cmat);
	CHKERRQ(ierr);

	double **dmag = new double*[num_dipoles];

	for (int n = 0; n < num_dipoles; n++)
		dmag[n] = new double[num_bsens * 6];

	for (int n = 0; n< num_dipoles; n++) {
		engine.clearRHS();

		mprintf("Dipole %d\n", n);

		ierr = engine.writeDipoleInfo(fi, 1, &dipoles[n]);
		CHKERRQ(ierr);
		ierr = engine.addSource(dipoles[n]);
		CHKERRQ(ierr);

		// write result to Phimat
		ierr = engine.solvePot(Phimat[n]);
		CHKERRQ(ierr);
		snprintf(fname, sizeof(fname), "x%03d.pot", n);

		ierr = engine.saveVector(fname, Phimat[n]);
		CHKERRQ(ierr);

		// solve primary magnetic field
		double *mag = engine.solveMagPri(num_bsens, bsensors);
		if (mag == NULL)
			SETERRQ(1, "Failed to solve magnetic field\n");

		double *mtot = dmag[n];
		double *mpri = mag;

		for (int s = 0; s < num_bsens; s++) {
			memcpy(mtot, mpri, 3 * sizeof(double));
			mtot += 6;
			mpri += 3;
		}
		delete[] mag;
	}

	for (int s = 0; s < num_bsens; s++) {
		mprintf("Sensor %d\n", s);

		// compute for 1 sensor (in place)
		ierr = engine.calcMagSecMat(1, &bsensors[s], Cmat);
		CHKERRQ(ierr);
		
		for (int d = 0; d < num_dipoles; d++) {
			double *mag = NULL;
			mag = engine.solveMagSec(Cmat, 1, &Phimat[d]);
			if (mag == NULL)
				SETERRQ(1, "Failed to solve magnetic field\n");
			double *dst = dmag[d];
			memcpy(dst+(6*s)+3, mag, 3 * sizeof(double));
			delete[] mag;
		}
		engine.clearVectors(3, Cmat);
	}

	for (int n = 0; n < num_dipoles; n++) {
		snprintf(fname, sizeof(fname), "x%03d.mag", n);
		ierr = engine.saveMag(fname, dmag[n], num_bsens);
		CHKERRQ(ierr);
		delete[] dmag[n];
		dmag[n] = 0;
	}

	delete[] dmag;

	if (fi)
		fclose(fi);

	ierr = VecDestroyVecs(Cmat, 3);
	CHKERRQ(ierr);

	ierr = VecDestroyVecs(Phimat, num_dipoles);
	CHKERRQ(ierr);

	return 0;
}

#undef __FUNCT__
#define __FUNCT__ "sens"

int
sens(FEngine &engine)
{
	int ierr;

	char fname[PATH_MAX];
	FILE *fi = NULL;

	if ((sensors == NULL && psensors == NULL) || num_sens < 1)
		SETERRQ(1, "no sensors\n");

	ierr = engine.setupSolver();
	CHKERRQ(ierr);

	Vec *Ainv = NULL;

	if (sensors != NULL)
		ierr = engine.invertSensCols(num_sens, sensors, &Ainv);
	else
		ierr = engine.invertSensCols(num_sens, psensors, &Ainv);

	CHKERRQ(ierr);

	assert(Ainv != NULL);

	try {
		fi = engine.createDipoleInfo("xinfo.sens", num_dipoles);
	} catch(...) {
		SETERRQ(PETSC_ERR_FILE_OPEN, "Failed to open xInfo file!\n");
	}

	for (int n = 0; n < num_dipoles; n++) {
		ierr = engine.clearRHS();
		CHKERRQ(ierr);

		mprintf("Dipole %d\n", n);

		ierr = engine.writeDipoleInfo(fi, 1, &dipoles[n]);
		CHKERRQ(ierr);

		ierr = engine.addSource(dipoles[n]);
		CHKERRQ(ierr);

		mprintf("Solving ...\n");

		ierr = engine.solvePot();
		CHKERRQ(ierr);
		
		snprintf(fname, sizeof(fname), "x%03d.sens", n);

		ierr = engine.savePot(fname);
		CHKERRQ(ierr);

		snprintf(fname, sizeof(fname), "sn%02d", n);

		ierr = engine.saveSensMatDip(fname,  Ainv, num_sens);
		CHKERRQ(ierr);
	}

	ierr = VecDestroyVecs(Ainv, num_sens);
	CHKERRQ(ierr);
	
	if (fi)
		fclose(fi);

	return 0;
}

#undef __FUNCT__
#define __FUNCT__ "sensB"

int
sensB(FEngine &engine)
{
	int ierr;

	char fname[PATH_MAX];
	FILE *fi = NULL;

	ierr = engine.setupSolver();
	CHKERRQ(ierr);

	if (bsensors == NULL || num_bsens < 1)
		SETERRQ(1, "No magnetic field sensors\n");

	Vec *Cmat;
	ierr = engine.calcMagSecMat(num_bsens, bsensors, &Cmat);
	CHKERRQ(ierr);

	if (Cmat == NULL)
		SETERRQ(1, "Failed to calculate matrix!\n");
	
	// compute nonzero columns

	ierr = engine.solvePot(Cmat, num_bsens * 3);
	CHKERRQ(ierr);

	try {
		fi = engine.createDipoleInfo("xinfo.bsp", num_dipoles);
	} catch(...) {
		SETERRQ(1, "Failed to open xInfo file!\n");
	}
	
	for (int n = 0; n < num_dipoles; n++) {
		ierr = engine.clearRHS();
		CHKERRQ(ierr);

		mprintf("Dipole %d\n", n);

		ierr = engine.writeDipoleInfo(fi, 1, &dipoles[n]);
		CHKERRQ(ierr);

		ierr = engine.addSource(dipoles[n]);
		CHKERRQ(ierr);

		mprintf("Solving ...\n");

		ierr = engine.solvePot();
		CHKERRQ(ierr);
		
		double *mag = engine.solveMag(Cmat, num_bsens, bsensors);
		if (mag == NULL)
			SETERRQ(1, "failed to solve magnetic field");

		snprintf(fname, sizeof(fname), "x%03d.bsp", n);
		ierr = engine.savePot(fname);
		CHKERRQ(ierr);

		if (mag != NULL) {
			snprintf(fname, sizeof(fname), "x%03d.bsm", n);
			ierr = engine.saveMag(fname, mag, num_bsens);
			CHKERRQ(ierr);
		}

		snprintf(fname, sizeof(fname), "bsn%02d", n);

		ierr = engine.saveMagSensMatDip(fname, Cmat,
						num_bsens, bsensors);
		CHKERRQ(ierr);
	}
	
	if (fi)
		fclose(fi);

	ierr = VecDestroyVecs(Cmat, num_bsens);
	CHKERRQ(ierr);

	return 0;
}


#undef __FUNCT__
#define __FUNCT__ "rfPot"
int
rfPot(FEngine &engine)
{
	int ierr;

	if ((sensors == NULL && psensors == NULL) || num_sens < 1)
		SETERRQ(1, "no sensors\n");
	
	ierr = engine.setupSolver();
	CHKERRQ(ierr);

	Vec *Ainv = NULL;
	if (sensors != NULL)
		ierr = engine.invertSensCols(num_sens, sensors, &Ainv);
	else
		ierr = engine.invertSensCols(num_sens, psensors, &Ainv);

	CHKERRQ(ierr);

	static char buf[PATH_MAX+1];

	for (int n = 0; n < num_sens; n++) {
		snprintf(buf, PATH_MAX, "rf%03d.pot", n);
		engine.saveVector(buf, Ainv[n]);
	}

	ierr = VecDestroyVecs(Ainv, num_sens);
	CHKERRQ(ierr);

	return 0;
}


#undef __FUNCT__
#define __FUNCT__ "rfMag"
int
rfMag(FEngine &engine, double rad)
{
	int ierr;

	Vec *Cmat = NULL;
	
	if (bsensors == NULL || num_bsens < 1)
		SETERRQ(1, "No magnetic field sensors\n");

	ierr = engine.setupSolver();
	CHKERRQ(ierr);

	// use radial direction for sensor
	ierr = engine.computeMagLead (num_bsens, bsensors,
				      bsensors, &Cmat, rad);
	CHKERRQ(ierr);

	static char buf[PATH_MAX+1];
	for (int n = 0; n < num_bsens; n++) {
		snprintf(buf, PATH_MAX, "rf%03d.mag", n);
		engine.saveVector(buf, Cmat[n]);
	}

	ierr = VecDestroyVecs(Cmat, num_bsens);
	CHKERRQ(ierr);

	return 0;
}


#undef __FUNCT__
#define __FUNCT__ "my_exit"

void
my_exit(int ret)
{
	PetscPrintf(PETSC_COMM_WORLD, "Exiting with code %d\n", ret);
	PetscFinalize();
	exit(ret);
}

static char help[] = \
"Solves EMSI forward problem in parallel with PetSc KSP Solver.\n\
Input parameters include:\n\
  -m <base_fn>   : Mesh filename\n\
  -d <dipole_fn> : Dipole filename(optional)\n\
  -s <sensor_fn> : Potential sensor filename(optional)\n\
  -b <bsens_fn>  : Magnetic sensor filename(optional)\n\
  -r <radius>    : Coil radius (default 1cm)\n\
  -sens          : Solve sensitivity to potentials\n\
  -bsens         : Solve magnetic field sensitivity\n\
  -rfpot         : Compute reciprocal field for potential sensors\n\
  -rfmag         : Compute reciprocal field for magbetic sensors\n\n";

#define MAX_BASE 1024
#define MAX_COND 100

#undef __FUNCT__
#define __FUNCT__ "main"
int
main(int argc, char **argv)
{
	vector<double> cmap;
	char buf[MAX_BASE];
	char *mpath, *cond[MAX_COND];
	char *dfn;
	char *mfn;
	char *sfn;

	int ierr;

	PetscTruth t, sflag, bsflag, rfpflag, rfmflag;
	PetscInt cmax;

	PetscInitialize(&argc, &argv, NULL, help);

	ierr = PetscOptionsGetString(PETSC_NULL, "-m", buf, MAX_BASE, &t);
	CHKERRQ(ierr);

	if (!t) {
		PetscPrintf(PETSC_COMM_WORLD, "%s", help);		
		my_exit(1);
	}

	if ((mpath = strdup(buf)) == NULL) {
		SETERRQ (PETSC_ERR_MEM, "strdup mpath");
	}

	cmax = MAX_COND; /* maximum 100 distinct conductivity labels supported */
	ierr = PetscOptionsGetStringArray(PETSC_NULL, "-c", cond, &cmax, &t);
	CHKERRQ(ierr);

	if (t) {
		/* XXX TODO Fill cmap */
		for (int i = 0; i < cmax; i++) {
			char *p = cond[i];
			char *ep = NULL;
			unsigned long ix = strtoul(p, &ep, 10);
			if (p == ep || ix <= 0 || ix > MAX_COND)
				SETERRQ (PETSC_ERR_MEM, "invalid conductivity index");
			if (*ep != '=')
				SETERRQ (PETSC_ERR_MEM, "expected =");
			p = ep + 1;

			double v = strtod(p, &ep);
			if (p == ep || *ep != '\0' || v < 0)
				SETERRQ (PETSC_ERR_MEM, "invalid conductivity value");
			if (cmap.size() < ix)
				cmap.resize(ix, -1);
			cmap[ix - 1] = v;
		}
	}

	for (unsigned int i = 0; i < cmap.size(); i++)
		PetscPrintf(PETSC_COMM_WORLD, "%d = %g\n", i, cmap[i]);

	/* Initialize Mesh */
	FEMesh *msh=new FEMesh();
	if (msh->loadMesh(mpath, &cmap)) {
		PetscPrintf(PETSC_COMM_WORLD, "Failed to load mesh!\n");
		my_exit(1);
	}

	for (PetscInt n = 0; n < cmax; n++) {
		ierr = PetscFree(cond[n]);
		CHKERRQ(ierr);
	}

	/* Print Mesh info */
	PetscPrintf(PETSC_COMM_WORLD, "Mesh initialized\n");
	PetscPrintf(PETSC_COMM_WORLD, "  %d nodes, %d elements\n",
		    msh->numNodes(), msh->numElements());
	FEShape *shape=msh->getShape();
	PetscPrintf(PETSC_COMM_WORLD,
		    "  Each element has %d nodes, %d faces\n",
		    shape->numNodes(), shape->numFaces());
	
	// Load Dipoles
	ierr = PetscOptionsGetString(PETSC_NULL, "-d", buf, MAX_BASE, &t);
	CHKERRQ(ierr);
	
	if (t) {
		if ((dfn = strdup(buf)) == NULL) {
			SETERRQ (PETSC_ERR_MEM, "strdup mpath");
		}
		mprintf("Reading dipoles from %s\n", dfn);
		dipoles = MeshUtil::readDipoles(dfn, num_dipoles,
						msh->numElements());
		if (num_dipoles)
			mprintf("%d Dipoles loaded\n", num_dipoles);
	}
	
	// Load MagSens
	ierr = PetscOptionsGetString(PETSC_NULL, "-b", buf, MAX_BASE, &t);
	CHKERRQ(ierr);
	
	if (t) {
		if ((mfn = strdup(buf)) == NULL) {
			SETERRQ (PETSC_ERR_MEM, "strdup mpath");
		}
		bsensors = ReadMagSensors(mfn, num_bsens);
		if (bsensors == NULL || num_bsens < 1) {
			mprintf("Failed to load magnetic field sensor file\n");
			bsensors = NULL;
			num_bsens = 0;
		} else
			mprintf("%d sensors loaded\n", num_bsens);
	}

	/* Load sensors */
	ierr = PetscOptionsGetString(PETSC_NULL, "-s", buf, MAX_BASE, &t);
	CHKERRQ(ierr);

	if (t) {
		if ((sfn = strdup(buf)) == NULL) {
			SETERRQ (PETSC_ERR_MEM, "strdup sfn");
		}
		
		mprintf("Reading sensors from %s\n", sfn);
		sensors = MeshUtil::readSensors(sfn, num_sens,
						msh->numNodes());
		if (num_sens)
			mprintf("%d Sensor nodes loaded\n", num_sens);
		else {
			sensors = NULL;
			psensors = ReadMagSensors(sfn, num_sens);
			if (num_sens)
				mprintf("%d Sensor coords loaded\n", num_sens);
		}
	}
	
	FEngine *eng = new FEngine(PETSC_COMM_WORLD, *msh);
	
	
	ierr = PetscOptionsGetTruth(PETSC_NULL, "-sens", &sflag, &t);
	CHKERRQ(ierr);
	if (t == PETSC_FALSE)
		sflag = PETSC_FALSE;
	
	ierr = PetscOptionsGetTruth(PETSC_NULL, "-bsens", &bsflag, &t);
	CHKERRQ(ierr);
	if (t == PETSC_FALSE)
		bsflag = PETSC_FALSE;

	ierr = PetscOptionsGetTruth(PETSC_NULL, "-rfpot", &rfpflag, &t);
	CHKERRQ(ierr);
	if (t == PETSC_FALSE)
		rfpflag = PETSC_FALSE;

	ierr = PetscOptionsGetTruth(PETSC_NULL, "-rfmag", &rfmflag, &t);
	CHKERRQ(ierr);
	if (t == PETSC_FALSE)
		rfmflag = PETSC_FALSE;
		
	// Coil Radius
	PetscReal radius;
	ierr = PetscOptionsGetReal(PETSC_NULL, "-r", &radius, &t);
	CHKERRQ(ierr);
	if (t == PETSC_FALSE || radius <= 0)
		radius = 0.01;

	printf("Coil radius %g\n", (double) radius);

	try {
		if (sflag) {
			ierr = sens(*eng);
			CHKERRQ(ierr);
		}
		
		if (bsflag) {
			ierr = sensB(*eng);
			CHKERRQ(ierr);
		}

		if (rfpflag) {
			ierr = rfPot(*eng);
			CHKERRQ(ierr);
		}

		if (rfmflag) {
			ierr = rfMag(*eng, radius);
			CHKERRQ(ierr);
		}
		
		if (!sflag && !bsflag && !rfpflag && !rfmflag) {
			if (num_dipoles == 0) {
				printf("Nothing to do (no dipoles)!\n");
				ierr = 0;
			} else if (3 * (num_bsens + 1) > num_dipoles) {
				// do not store C, store phi instead
				ierr = forward2(*eng);
			} else {
				ierr = forward(*eng);
			}
			CHKERRQ(ierr);
		}

	} catch (int ier) {
		CHKERRQ(ier);
	} catch(...) {
		SETERRQ(1, "Unknown exception\n");
	}
	
	ierr = PetscFinalize();
	CHKERRQ(ierr);
	
	HPTimerMgr::report();
	
	if(msh) delete msh;
	
	return 0;
}

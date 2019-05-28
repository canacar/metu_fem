/* $Id: meshutil.cc,v 1.3 2008/01/28 07:41:13 canacar Exp $ */
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
/*! \file meshutil.cpp
    \brief Mesh Utility functions.
    Contains MeshUtil class implementation.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "meshutil.h"

//---------------------------------------------------------------------------
char MeshUtil::m_buf[MMGR_BUFSIZE];

/*! reads in an integer list.
    Used internally by loadIntList to load a list of intgers of the form
    [index] [value]
    \param f mesh file descriptor
    \param lst pre-allocated list
    \param start starting index of list
    \param size number of items to load (size of lst)
    \return 0 if success
    \sa loadIntList()
*/
int MeshUtil::readIntList(FILE *f, int *lst, int start, int size){
  int n,i,v;

  for(n=0; n<size; n++){
    if(fgets(m_buf, MMGR_BUFSIZE, f)==NULL) return 1;
    if(sscanf(m_buf, "%d %d",&i, &v)!=2) return 1;
    if(i!=(n+start)) return 1;
    lst[n]=v;
  }

  return 0;
}

/*! reads in a list of doubles.
    Used internally by scanMesh to load a list of doubles of the form
    [index] [value]
    \param f mesh file descriptor
    \param lst pre-allocated list
    \param start starting index of list
    \param size number of items to load (size of lst)
    \return 0 if success
    \sa scanMesh()
*/
int MeshUtil::readDoubleList(FILE *f, double *lst, int start, int size){
  int n,i;
  double v;

  for(n=0; n<size; n++){
    if(fgets(m_buf, MMGR_BUFSIZE, f)==NULL) return 1;
    if(sscanf(m_buf, "%d %lf",&i, &v)!=2) return 1;
    if(i!=(n+start)) return 1;
    lst[n]=v;
  }

  return 0;
}

/*! reads in a list of points.
    Used to load a list of doubles of the form
    [index] [x] [y] [z]
    \param f mesh file descriptor
    \param pts pre-allocated list of points
    \param start starting index of list
    \param size number of items to load (size of lst)
    \return 0 if success
*/
int MeshUtil::readPointList(FILE *f, point *pts, int start, int size){
  int n,i;
  double x,y,z;

  for(n=0; n<size; n++){
    if(fgets(m_buf, MMGR_BUFSIZE, f)==NULL) return 1;
    if(sscanf(m_buf, "%d %lg %lg %lg",&i, &x,&y,&z)!=4) return 1;
    if(i!=(n+start)) return 1;
    pts[n][0]=x;
    pts[n][1]=y;
    pts[n][2]=z;
  }

  return 0;
}

/*! reads the list header.
    the list header contains the size followed by an optional description
    [size] [description (optional)]
    \param f mesh file descriptor
    \param size pointer to the variable that will hold the size of list
    \return 0 if success
*/
int MeshUtil::readListHdr(FILE *f, int *size){
  if(fgets(m_buf, MMGR_BUFSIZE, f)==NULL) return 1;     /* # of elements */
  if(sscanf(m_buf, "%d", size)!=1) return 1;
  if(*size<=0) return 1;
  return 0;
}

/*! reads the element list header.
    the list header contains the size info followed by an optional description
    [size] [node/elem] [description (optional)]
    \param f mesh file descriptor
    \param size pointer to the variable that will hold the size of list
    \param nen  number of nodes in element
    \return 0 if success
*/
int MeshUtil::readElemListHdr(FILE *f, int *size, int *nen){
  if(fgets(m_buf, MMGR_BUFSIZE, f)==NULL) return 1;     /* # of elements */
  int n=sscanf(m_buf, "%d %d",size, nen);
  if (n == 1) *nen=20;
  else if (n!=2) return 1;
  if (*nen <1) return 1;
  if(*size<=0) return 1;
  return 0;
}

/*! Allocates and Loads a list of integers.
    \param f mesh file descriptor
    \param size pointer to the variable that will hold the size of list
    \return allocated list, 0 if error
*/
int *MeshUtil::loadIntList(FILE *f, int *size){
  int sz, *l;

  if(readListHdr(f, &sz)) return NULL;
  l=new int[sz];
  if(l==NULL) return NULL;
  if(readIntList(f,l,1,sz)){
    delete[] l;
    return NULL;
  }
  *size=sz;
  return l;
}

/*! skips a given number of list items.
    The indices are verified but type is not important
    \param f mesh file descriptor
    \param start starting index of list
    \param size number of items to skip
    \return 0 if success
*/
int MeshUtil::skipLineItems(FILE *f, int start, int size){
  int n,i;

  for(n=0; n<size; n++){
    if(fgets(m_buf, MMGR_BUFSIZE, f)==NULL) return 1;
    if(sscanf(m_buf, "%d",&i)!=1) return 1;
    if(i!=(n+start)) return 1;
  }

  return 0;
}

/*! skips a given number of node items.
    The indices are verified but type is not important
    \param f mesh file descriptor
    \param start starting index of list
    \param size number of items to skip
    \return 0 if success
*/
int MeshUtil::skipNodeItems(FILE *f, int start, int size, int en){
  int n,m,i;

  for(n=0; n<size; n++){
    if(fscanf(f, "%d",&i)!=1) return 1;
    if(i!=(n+start)) return 1;
    for(m=0; m<en; m++)
      if(fscanf(f,"%d",&i)!=1) return 1;    // read individual nodes why?
    fgets(m_buf,MMGR_BUFSIZE,f);    // discard to end of line
  }
  return 0;
}

/*! reads in an list of elements.
    Used internally by loadIntList to load a list of elements of the form
    [index] [index1] ... [index20]
    \param f mesh file descriptor
    \param nds pre-allocated list of nodes
    \param start starting index of list
    \param size number of items to load (size of lst)
    \return 0 if success
*/
int MeshUtil::readElemList(FILE *f, NodeArray &nds, int start,
			   int size, int en, int nn){
  int n,i;

  for(n=0; n<size; n++){
    if(fscanf(f, "%d",&i)!=1) return 1;
    if(i!=(n+start)) return 1;

    for(int k=0; k<en; k++){
      if(fscanf(f,"%d",&i)!=1) return 1;
      i--;
      if(i<0 || i>=nn) return 1;
      nds[n].node(k)=i;
    }
    fgets(m_buf,MMGR_BUFSIZE,f);    // discard to end of line
  }

  return 0;
}

/*! reads in an list of dipoles.
    [index] [elem] [dx] [dy] [dz] [Jx] [Jy] [Jz]
    \param f mesh file descriptor
    \param nds pre-allocated list of nodes
    \param start starting index of list
    \param size number of items to load (size of lst)
    \param number of elements in mesh (for input validation)
    \return 0 if success
*/
int MeshUtil::readDipoleList(FILE *f, DInfo *dip, int start,
			   int size, int ne){
  int n,i,nf;
  int pdg = 0;

  for(n=0; n<size; n++){
    int el, sm, dg;
    double dx, dy, dz;
    double jx, jy, jz;

    fgets(m_buf, MMGR_BUFSIZE,f);
    nf = sscanf(m_buf, "%d %d %d %lg %lg %lg %lg %lg %lg %d",
               &i, &sm, &el, &dx, &dy, &dz, &jx, &jy, &jz, &dg);

    if (nf == 9)
      dg = n;
    else if (nf != 10)
      return 1;

    if(i!=(n+start)) return 1;
    if (el <0 || el >=ne) return 1;
    if (dx < -1 || dx > 1) return 1;
    if (dy < -1 || dy > 1) return 1;
    if (dz < -1 || dz > 1) return 1;

    /* make sure dipole group is increasing */
    if (n == 0)
      pdg = dg;
    else if (dg != pdg) {
        pdg ++;
      if (pdg != dg)
        return 1;
    }

    if (sm != SM_J && sm != SM_YAN)
      return 1;

    dip[n].elem=el;
    dip[n].smodel=(SourceModelID) sm;
    dip[n].dgroup = dg;
    dip[n].D[0]=dx;
    dip[n].D[1]=dy;
    dip[n].D[2]=dz;
    dip[n].J[0]=jx;
    dip[n].J[1]=jy;
    dip[n].J[2]=jz;
  }

  return 0;
}

/*! reads in an list of shell definitions for concentric spheres model.
    [index] [radius] [sigma]
    radius must be > 0 and increasing.
    \param f mesh file descriptor
    \param nds pre-allocated list of nodes
    \param start starting index of list
    \param size number of items to load (size of lst)
    \return 0 if success
*/
int MeshUtil::readShellList(FILE *f, Shell *sh, int start, int size){
	int n,i;
	double prev = 0;
	for(n=0; n<size; n++){
		double rad, sig;
		
		fgets(m_buf, MMGR_BUFSIZE,f);
		if(sscanf(m_buf, "%d %lg %lg", &i, &rad, &sig) != 3)
			return 1;
		if(i!=(n+start)) return 1;
		if (rad  <= prev) return 1;
		if (sig < 0) return 1;
		
		sh[n].sig=sig;
		sh[n].radius=rad;
		prev = rad;
	}
	
	return 0;
}

/*! checks the signature of a mesh file
    The mesh files start with
    0 [signature]
    \param f mesh file descriptor
    \param sign id the text that will be compared with the signatire
    \return 0 if success
*/
int MeshUtil::checkHeader(FILE *f, const char *sign){
  if(fgets(m_buf,MMGR_BUFSIZE,f)==NULL) return 1;
  if(m_buf[0]!='0') return 1;
  if(strncmp(&m_buf[2],sign,strlen(sign))) return 1;
  if(fgets(m_buf,MMGR_BUFSIZE,f)==NULL) return 1;
  return 0;
}

int *MeshUtil::readSensors(char *fn, int &num, int nn)
{
  int *sens=0;

  if (fn == NULL) return 0;

  FILE *f=fopen(fn,"r");

  if (f == NULL) {
    mprintf("Failed to open sensor file!\n");
    return 0;
  }

  if(checkHeader(f, "Parallel Mesh Sensor Nodes")) goto err;

  if(readListHdr(f, &num)) goto err;
  if(num <= 0) goto err;

  sens = new int[num];

  if (readIntList(f, sens, 1, num)) goto err;

  for (int n=0; n<num; n++)
    if (sens[n] < 0 || sens[n] >=nn) goto err;

  fclose(f);
  return sens;

 err:
  fclose(f);
  if (sens) delete[] sens;
  return 0;

}

DInfo *MeshUtil::readDipoles(char *fn, int &num, int ne)
{
  DInfo *dip = 0;

  if (fn == NULL) return 0;

  FILE *f=fopen(fn,"r");
  if (f == NULL) {
    mprintf("Failed to open dipole file %s!\n", fn);
    return 0;
  }

  if(checkHeader(f, "Parallel Mesh Dipole Information")) goto err;

  if(readListHdr(f, &num)) goto err;
  if(num <= 0) goto err;

  dip = new DInfo[num];

  if (readDipoleList(f, dip, 1, num, ne)) goto err;

  fclose(f);
  return dip;

 err:
  num = 0;
  fclose(f);
  if (dip) delete[] dip;
  return 0;
}

Shell *MeshUtil::readShells(char *fn, int &num)
{
	Shell *sh = 0;
	
	if (fn == NULL) return 0;
	
	FILE *f=fopen(fn,"r");
	if (f == NULL) {
		mprintf("Failed to open dipole file!\n");
		return 0;
	}

	// skip first line: mesh name
	fgets(m_buf, MMGR_BUFSIZE,f);

	if(readListHdr(f, &num)) goto err;
	if(num <= 0) goto err;
	
	sh = new Shell[num];
	
	if (readShellList(f, sh, 1, num)) goto err;
	
	fclose(f);
	return sh;
	
 err:
	fclose(f);
	if (sh) delete[] sh;
	return NULL;
}

/* $Id: femmesh3.cc,v 1.7 2008/11/24 04:15:52 canacar Exp $ */
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
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <deque>
#include <assert.h>
#include <math.h>
#include "femmesh3.h"
#include "shape20.h"
#include "shape10.h"
#include "shape8.h"
#include "shape4.h"

//---------------------------------------------------------------------------
using namespace std;

FaceList FEMesh::m_facemap20[]={{12, 14, 16, 18},
				 { 6,  0, 12, 18},
				 { 0,  2, 14, 12},
				 { 2,  4, 16, 14},
				 { 6,  4,  2,  0},
				 {18, 16,  4,  6}};
/*
  Needed later
FEMesh::m_face8map20[]={{12, 13, 14, 15, 16, 17, 18, 19},
                        { 6,  7,  0,  8, 12, 19, 18, 11},
			{ 0,  1,  2,  9, 14, 13, 12,  8},
			{ 2,  3,  4, 10, 16, 15, 14,  9},
			{ 6,  7,  0,  1,  2,  3,  4,  5},
			{18, 17, 16, 10,  4,  5,  6, 11}};
*/

FaceList FEMesh::m_facemap8[]={{4, 5, 6, 7},
				{3, 0, 4, 7},
				{0, 1, 5, 4},
				{1, 2, 6, 5},
				{3, 2, 1, 0},
				{7, 6, 2, 3}};


FaceList FEMesh::m_facemap10[]={{3, 2, 1, 0xff},
				{0, 2, 3, 0xff},
				{0, 3, 1, 0xff},
				{0, 1, 2, 0xff}};

FaceList FEMesh::m_facemap4[]={{3, 2, 1, 0xff},
				{0, 2, 3, 0xff},
				{0, 3, 1, 0xff},
				{0, 1, 2, 0xff}};


int FEMesh::m_nnbrmap20[][3]={{  1,  7,  8},  // 0
                              {  0,  2, -1},  // 1
                              {  1,  3,  9},  // 2
                              {  2,  4, -1},  // 3
                              {  3,  5, 10},  // 4
                              {  4,  6, -1},  // 5
                              {  5,  7, 11},  // 6
                              {  0,  6, -1},  // 7
                              {  0, 12, -1},  // 8
                              {  2, 14, -1},  // 9
                              {  4, 16, -1},  // 10
                              {  6, 18, -1},  // 11
                              {  8, 13, 19},  // 12
                              { 12, 14, -1},  // 13
                              {  9, 13, 15},  // 14
                              { 14, 16, -1},  // 15
                              { 10, 15, 17},  // 16
                              { 16, 18, -1},  // 17
                              { 11, 17, 19},  // 18
                              { 12, 18, -1}}; // 19

int FEMesh::m_nnbrmap8[][3]={{1, 3, 4},  // 0
                             {0, 2, 5},  // 1
                             {1, 3, 6},  // 2
                             {0, 2, 7},  // 3
                             {0, 5, 7},  // 4
                             {1, 4, 6},  // 5
                             {2, 5, 7},  // 6
                             {3, 4, 6}}; // 7

int FEMesh::m_nnbrmap10[][3]={{4, 6, 7},  // 0
			      {4, 5, 8},  // 1
			      {5, 6, 9},  // 2
			      {7, 8, 9},  // 3
			      {0, 1, -1},  // 4
			      {1, 2, -1},  // 5
			      {0, 2, -1},  // 6
			      {0, 3, -1},  // 7
			      {1, 3, -1},  // 8
			      {2, 3, -1}}; // 9

int FEMesh::m_nnbrmap4[][3]={{1, 2, 3},  // 0
                             {0, 2, 3},  // 1
                             {0, 1, 3},  // 2
                             {0, 1, 2}}; // 3

int FEMesh::m_edge20[][2]={{ 0, 2}, { 2, 4}, { 4, 6}, { 6, 0},
			   { 0,12}, { 2,14}, { 4,16}, { 6,18},
			   {12,14}, {14,16}, {16,18}, {18,12}};

int FEMesh::m_edge8[][2]={{0,1}, {1,2}, {2,3}, {3,0},
			  {0,4}, {1,5}, {2,6}, {3,7},
			  {4,5}, {5,6}, {6,7}, {7,4}};

int FEMesh::m_edge10[][2]={{0,1}, {0,2}, {0,3},
			  {1,2}, {1,3}, {2,3}};

int FEMesh::m_edge4[][2]={{0,1}, {0,2}, {0,3},
			  {1,2}, {1,3}, {2,3}};

#define MAX_CORNER_NODE 8
// defaults are for 20 noded, loadMesh modifies it for other elements
int FEMesh::m_corner[]={12,14,16,18,0,2,4,6};

//---------------------------------------------------------------------------
FENeighbor::FENeighbor()
{
    for(int n=0; n<MAX_FACE_ELEM; n++){
        nelem[n]=-1;
        nface[n]=0xff;
    }
}
//---------------------------------------------------------------------------
FENeighbor::~FENeighbor()
{
}
//---------------------------------------------------------------------------
FENeighbor::FENeighbor(const FENeighbor &nbr)
{
    memcpy(nelem, nbr.nelem, MAX_FACE_ELEM*sizeof(int));
    memcpy(nface, nbr.nface, MAX_FACE_ELEM*sizeof(unsigned char));
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
FENodeInfo &FENodeInfo::operator=(const FENodeInfo &i)
{
    flags = i.flags;
    nbrs = i.nbrs;
    return *this;
}
//---------------------------------------------------------------------------
FENodeInfo::FENodeInfo(const FENodeInfo &i) : flags(i.flags), nbrs(i.nbrs)
{
    flags = i.flags;
}
//---------------------------------------------------------------------------
void FENodeInfo::addNeighbor(int nbr)
{
    for (size_t n=0; n < nbrs.size(); n++)
        if (nbrs[n] == nbr) return;
    nbrs.push_back(nbr);
}
//---------------------------------------------------------------------------
void FENodeInfo::delNeighbor(int nbr)
{
    size_t n;
    for (n=0; n < nbrs.size(); n++)
        if (nbrs[n] == nbr) break;
    assert (n < nbrs.size());  // not found
    nbrs[n] = nbrs.back();
    nbrs.pop_back();
}
//---------------------------------------------------------------------------
int
FEMesh::parseSigmaLine(const char *buf, int idx, vector<double> *cmap, double *sig, int *cid)
{
	const char *sp;
	char *ep;
	double val;

	sp = buf;
	val = strtod(sp, &ep);
	if (ep == sp)
		return 1;
	if (!isspace(*ep))
		return 1;
	if (idx >= 0 && ((int)floor(val) != idx))
		return 1;
	for (sp = ep; *sp && isspace(*sp); sp++);

	if (*sp == 0)
		return 1;

	if (*sp == 'C' || *sp == 'c') {
		unsigned long ci = strtol(++sp, &ep, 10);
		if (cmap == NULL || sp == ep) {
			fprintf(stderr, "invalid label %s [%s]\n", buf, sp);
			return 1;
		}
		if (ci <= 0 || ci > cmap->size()) {
			fprintf(stderr, "invalid index %ld [1 -- %ld]\n", (long)ci, (long)cmap->size());
			return 1;
		}
		val = (*cmap)[ci - 1];
		if (cid != NULL)
			*cid = ci;
	} else {
		if (cid != NULL)
			*cid = -1;
		val = strtod(sp, &ep);
		if (sp == ep) {
			fprintf(stderr, "Invalid conductivity: %s\n", sp);
			return 1;
		}
	}

	if (sig != NULL)
		*sig = val;

	return 0;
}
//---------------------------------------------------------------------------
#define MESHBUFSIZE 4096
int FEMesh::loadMesh(const char *name, vector<double> *cmap)
{
    static char buf[MESHBUFSIZE];
    int nn, ne, ns;

    FILE *f=fopen(name,"rt");
    if (f==NULL) return 1;

    if(fgets(buf, MESHBUFSIZE, f)==NULL) return 1;     // # of nodes
    if(sscanf(buf, "%d", &nn)!=1) return 1;
    if(nn<0) return 1;

    if (nn == 0) {      // 0 if it is a parallel mesh
        if(fgets(buf, MESHBUFSIZE, f)==NULL) return 1;     // another zero
        if(sscanf(buf, "%d", &nn)!=1) return 1;
        if(nn!=0) return 1;

        if(fgets(buf, MESHBUFSIZE, f)==NULL) return 1;     // # of nodes
        if(sscanf(buf, "%d", &nn)!=1) return 1;
        if(nn<=0) return 1;
    }

    if (m_meshfn != NULL)
	    free(m_meshfn);	// allocated by strdup

    m_meshfn = NULL;

    // create node info
    int start=m_nnodes;

    // reserve space for nodes
    reserveSpace(0, nn);

fprintf(stderr, "loading %d nodes\n", nn);
    for(int n=0; n<nn; n++){
        int i;
        double it, x,y,z;

        if(fgets(buf, MESHBUFSIZE, f)==NULL) return 1;
        if(sscanf(buf, "%lg %lg %lg %lg",&it, &x,&y,&z)!=4) return 1;
        i=(int) floor(it);
        if(it!=it) return 1;
        if(i!=(n+1)) return 1;
        Point3 nd(x,y,z);
        if(addNode(nd)) return 1;
    }

    if(fgets(buf, MESHBUFSIZE, f)==NULL) return 1;     // # of elements
    int r=sscanf(buf, "%d %d", &ne, &ns);
    if (r == 1) ns = 20;
    else if (r !=2) return 1;
    else printf("%d nodes/elem\n",ns);

    if (ns !=20 && ns !=8 && ns != 10 && ns != 4) return 1;
    if(ne<=0) return 1;
    m_nnelem=ns;

    m_meshfn = strdup(name);

    if (m_nnelem == 4) {
        m_shape = new FEShape4();
	// correct m_corner
        for (int n=0; n < MAX_CORNER_NODE; n++)
                m_corner[n] = (n < m_nnelem) ? n : -1;
	m_nnface = 3;
	m_nfelem = 4;
    }else if (m_nnelem == 10) {
        m_shape = new FEShape10();
	// correct m_corner
        for (int n=0; n < MAX_CORNER_NODE; n++)
                m_corner[n] = (n < 4) ? n : -1;
	m_nnface = 3;
	m_nfelem = 4;
    } else if (m_nnelem==8) {
        m_shape=new FEShape8();
	// correct m_corner
        for (int n=0; n<MAX_CORNER_NODE; n++)
                m_corner[n]=n;
	m_nnface = 4;
	m_nfelem = 6;
    } else {
        m_shape=new FEShape20();
	m_nnface = 4;
	m_nfelem = 6;
    }

    // reserve space for elements
    reserveSpace(ne, 0);
fprintf(stderr, "loading %d elements\n", ne);
    for(int n=0; n<ne; n++){
        int i;
	int err = 1;
	do {
            RElem el;
           if(fscanf(f, "%d", &i)!=1) break;
            if(i!=(n+1)) break;
            for(int m=0; m<ns; m++){
                if(fscanf(f,"%d",&i)!=1) break;
                i--;
                if(i<0 || i>=nn) break;
                el.elem[m]=i+start;
            }
            if (addElem(el)) {
	   	    fprintf(stderr, "Failed to add Element %d!\n", n);
		    return 1;
	    }
            fgets(buf,MESHBUFSIZE,f);    // discard to end of line
	    err = 0;
	} while (0);
	if (err) {
		fprintf(stderr, "Error parsing line %d\n", n + 1);
	}
    }

    int nsig, cls;
   
fprintf(stderr, "loading %d sigma\n", ne);
    if(fgets(buf, MESHBUFSIZE, f)==NULL) return 0;     // elem sigma, optional
    if(sscanf(buf, "%d", &nsig)!=1) return 0;

    if(nsig == ne) {

        m_sigelem=new double[ne];
	if (cmap)
		m_sigelemcls = new int[ne];
	else
		m_sigelemcls = NULL;
    
        for(int n=0; n<ne; n++){
            double sig;
            if(fgets(buf, MESHBUFSIZE, f)==NULL) {
		fprintf(stderr, "failed to read sigma line for element %d\n", n + 1);
		return 1;
	    }
            if (parseSigmaLine(buf, n + 1, cmap, &sig, &cls)) {
		fprintf(stderr, "failed to parse sigma line for element %d\n", n + 1);
                return 1;
	    }
            m_sigelem[n]=sig;
            if (cmap)
		m_sigelemcls[n]=cls;
        }

        if(fgets(buf, MESHBUFSIZE, f)==NULL) return 0;     
        if(sscanf(buf, "%d", &nsig)!=1) return 0;
    }

    if(nsig == nn) {                                // node sigma, optional

        m_signode=new double[nn];
	if (cmap)
		m_signodecls = new int[nn];
	else
		m_signodecls = NULL;
    
        for(int n=0; n<nn; n++){
            double sig;
            if(fgets(buf, MESHBUFSIZE, f)==NULL) return 1;
            if (parseSigmaLine(buf, n + 1, cmap, &sig, &cls))
                return 1;
            m_signode[n]=sig;
            if (cmap)
		m_sigelemcls[n]=cls;
        }
    }

    return 0;
}
//---------------------------------------------------------------------------
int FEMesh::initMesh(int nnelem)
{
    m_meshfn = NULL;
    m_nnodes = 0;
    m_nelem = 0;

    if (nnelem !=20 && nnelem !=8 && nnelem != 4 && nnelem != 10)
	    return 1;

    m_nnelem = nnelem;

    if (m_nnelem == 4) {
        m_shape = new FEShape4();
	// correct m_corner
        for (int n=0; n < MAX_CORNER_NODE; n++)
                m_corner[n] = (n < m_nnelem) ? n : -1;
	m_nnface = 3;
	m_nfelem = 4;
    } else if (m_nnelem == 10) {
        m_shape = new FEShape10();
	// correct m_corner
        for (int n=0; n < MAX_CORNER_NODE; n++)
                m_corner[n] = (n < m_nnelem) ? n : -1;
	m_nnface = 3;
	m_nfelem = 4;
    } else if (m_nnelem==8) {
        m_shape=new FEShape8();
	// correct m_corner
        for (int n=0; n<MAX_CORNER_NODE; n++)
                m_corner[n]=n;
	m_nnface = 4;
	m_nfelem = 6;
    } else {
        m_shape=new FEShape20();
	m_nnface = 4;
	m_nfelem = 6;
    }

    return 0;
}
//---------------------------------------------------------------------------
FEMesh::FEMesh()
{
    m_nnodes=0;
    m_nelem=0;
    m_signode=0;
    m_sigelem=0;
    m_shape=0;
    m_nnface = 0;
    m_nfelem = 0;
    m_meshfn = NULL;
    m_scache = NULL;
}
//---------------------------------------------------------------------------
FEMesh::~FEMesh(void)
{
    if (m_sigelem) delete[] m_sigelem;
    if (m_signode) delete[] m_signode;
    if (m_shape) delete m_shape;
    if (m_meshfn != NULL)
	    free(m_meshfn);	// allocated by strdup
    if (m_scache != NULL)
	    delete m_scache;	// allocated by strdup
}
//---------------------------------------------------------------------------
int FEMesh::findNeighbor(int &ne, int &nf)
{
    assert(ne>=0 && ne<m_nelem);
    assert(nf>=0 && nf<m_nfelem);
    FENeighbor &fn=m_nbrs[ne];
    int e=fn.nelem[nf];
    int f=fn.nface[nf];
    if(e<0) return 1;
    ne=e;
    nf=f;
    return 0;
}
//---------------------------------------------------------------------------
int FEMesh::getFaceNodes(int elem, int face, int *nodes) const
{
    if(elem<0 || elem>=m_nelem)
        assert(elem>=0 && elem<m_nelem);

    assert(face>=0 && face<m_nfelem);
    assert(nodes);
    const node *nd=getElem(elem);
    if (m_nnelem == 20) {
        for(int n=0; n<m_nnface; n++)
            nodes[n]=(*nd)[m_facemap20[face][n]];
    } else if (m_nnelem == 8) {
        for(int n=0; n<m_nnface; n++)
            nodes[n]=(*nd)[m_facemap8[face][n]];
    } else if (m_nnelem == 10) {
        for(int n=0; n<m_nnface; n++)
            nodes[n]=(*nd)[m_facemap10[face][n]];
    } else {
        for(int n=0; n<m_nnface; n++)
            nodes[n]=(*nd)[m_facemap4[face][n]];
    }

    return 0;
}
//---------------------------------------------------------------------------
void FEMesh::addNodeNeighbors(int el, const node *nds)
{
    assert(nds);
    assert(el>=0 && el < m_nelem);
    // add element as neighbor of nodes
    for (int n=0; n < m_nnelem; n++) {
        int nd=(*nds)[n];
        assert(nd>=0 && nd<m_nnodes);
        m_ninfo[nd].addNeighbor(el);
    }
}
//---------------------------------------------------------------------------
void FEMesh::delNodeNeighbors(int el, const node *nds)
{
    assert(nds);
    assert(el>=0 && el < m_nelem);
    // remove element from neighbor list of nodes
    for (int n=0; n < m_nnelem; n++) {
        int nd=(*nds)[n];
        assert(nd>=0 && nd<m_nnodes);
        m_ninfo[nd].delNeighbor(el);
    }
}
//---------------------------------------------------------------------------
// check if elements are neighbors, if so find matching faces
void FEMesh::checkNeighbor(int e1, const node* n1, int e2, const node *n2)
{
    assert(n1);
    assert(n2);
    assert(e1>=0 && e1<m_nelem);
    assert(e2>=0 && e2<m_nelem);

    static int common[MAX_NODE_FACE];
    static int cn1[MAX_NODE_FACE];  // common nodes in e1
    static int cn2[MAX_NODE_FACE];  // common nodes in e2

    int ncom=0;

    for(int i1=0; i1<MAX_CORNER_NODE; i1++){
        int nd1 = m_corner[i1];
	if (nd1 < 0) continue;
        nd1=(*n1)[nd1];
        int i2;
        for(i2=0; i2<MAX_CORNER_NODE; i2++){
	    int nd2 = m_corner[i2];
	    if (nd2 < 0) continue;
            nd2=(*n2)[nd2];
            if(nd1==nd2) break;
        }
        if (i2==MAX_CORNER_NODE) continue;
	assert(ncom < m_nnface);

        common[ncom]=nd1;
        cn1[ncom]=i1;
        cn2[ncom]=i2;
        ncom++;
    }
    if (ncom!=m_nnface)
	    return; // not face neighbor

    int f1=locateFace(cn1);
    assert(f1>=0 && f1<m_nfelem);

    int f2=locateFace(cn2);
    assert(f2>=0 && f2<m_nfelem);

    // check needed ???
    assert(m_nbrs[e1].nelem[f1]==-1);
    m_nbrs[e1].nelem[f1]=e2;
    m_nbrs[e1].nface[f1]=(unsigned char)f2;
    assert(m_nbrs[e2].nelem[f2]==-1);
    m_nbrs[e2].nelem[f2]=e1;
    m_nbrs[e2].nface[f2]=(unsigned char)f1;
}
//---------------------------------------------------------------------------
int FEMesh::locateFace(int *nodes)
{
    assert(nodes);
    int f;
    for(f=0; f<m_nfelem; f++){
        int n;
        for(n=0; n<m_nnface; n++){
            int nd=m_corner[nodes[n]];
            int i;
            for(i=0; i<m_nnface; i++)
                if (m_nnelem == 20) {
                        if (m_facemap20[f][i]==nd) break;
                } else if (m_nnelem == 8) {
                        if (m_facemap8[f][i]==nd) break;
                } else if (m_nnelem == 10) {
                        if (m_facemap10[f][i]==nd) break;
                } else {
                        if (m_facemap4[f][i]==nd) break;
		}
            if (i == m_nnface) break;
        }
        if (n == m_nnface) return f;
    }
    return -1;
}
//---------------------------------------------------------------------------
void FEMesh::elemCoord(int elem, int idx, double &x, double &y, double &z)
{
    assert(elem>=0 && elem<m_nelem);
    assert(idx>=0 && idx<m_nnelem);

    node *nd=getElem(elem);
    const Point3 &pt=getNode((*nd)[idx]);
    x=pt.getX();
    y=pt.getY();
    z=pt.getZ();
}
//---------------------------------------------------------------------------
node *FEMesh::getElem(int elem)
{
    assert(elem>=0 && elem<m_nelem);
    return &m_elements[elem].elem;
}
//---------------------------------------------------------------------------
const node *FEMesh::getElem(int elem) const
{
    assert(elem>=0 && elem<m_nelem);
    return &m_elements[elem].elem;
}
//---------------------------------------------------------------------------
Point3 &FEMesh::getNodeRef(int nd)
{
    assert(nd>=0 && nd<m_nnodes);
    return m_nodes[nd];
}
//---------------------------------------------------------------------------
const Point3 &FEMesh::getNode(int nd) const
{
    assert(nd>=0 && nd<m_nnodes);
    return m_nodes[nd];
}
//---------------------------------------------------------------------------
void FEMesh::globalCoord(int elem, double &x, double &y, double &z) const
{
    assert(elem>=0 && elem<m_nelem);
    assert(x>=-1 && x<=1);
    assert(y>=-1 && y<=1);
    assert(z>=-1 && z<=1);

    static point pa[NUM_NODES];
    pArray(elem,pa);

    m_shape->Local2Global(pa, x, y, z);
}
//---------------------------------------------------------------------------
int
FEMesh::localCoord(int elem, double &x, double &y, double &z,
		   double *err, int *ni) const
{
    assert(elem>=0 && elem<m_nelem);

    point pa[NUM_NODES];
    pArray(elem,pa);

    double i=0, j=0, k=0;
    int iter = 0;
    double J[3][3];	// Jacobian
    double IJ[3][3];	// Jacobian

    double px, py, pz;
    double dx, dy, dz;
    px = py = pz = 0;

    double alpha = 0.8;
    double derr, perr;

#define MAX_ITER 1000
#define PE_TRESH 1e-12
#define DE_TRESH 1e-25

    if (m_nnelem == 4 || m_nnelem == 10)
	    i = j = k = 0.25;

    while (iter++ < MAX_ITER) {
	    // present position
	    dx = i; dy = j; dz = k;
	    // compute present global coordinate
	    m_shape->Local2Global(pa, dx, dy, dz);

	    derr = fabs(px - dx) + fabs(py - dy) + fabs(pz - dz);
	    
	    // save previous location
	    px = dx; py = dy; pz = dz;
	    // -dX
	    dx -= x;
	    dy -= y;
	    dz -= z;

	    perr = fabs(dx) + fabs(dy) + fabs(dz);

	    if (perr < PE_TRESH || derr < DE_TRESH)
		    break;

	    // compute jacobian
	    if (m_shape->jacob(pa, i, j, k, J, IJ)) {
		    fprintf(stderr, "Error computing Jacobian!\n");
		    break;
	    }

	    // correct I
	    i -= alpha * (dx * IJ[0][0] + dy * IJ[1][0] + dz * IJ[2][0]);
	    j -= alpha * (dx * IJ[0][1] + dy * IJ[1][1] + dz * IJ[2][1]);
	    k -= alpha * (dx * IJ[0][2] + dy * IJ[1][2] + dz * IJ[2][2]);
    }

    int ret = 0;

    if (i > 1 || j > 1 || k > 1)
	    ret++;
    else if (m_nnelem == 4 || m_nnelem == 10) {
	    if (i < 0 || j < 0 || k < 0 || (i + j + k) > 1)
		    ret++;
    } else if (i < -1 || j < -1 || k < -1)
	    ret ++;

    x = i;
    y = j;
    z = k;

    if (ni)
	    *ni = iter;
    if (err)
	    *err = perr;

    return ret;
}

//---------------------------------------------------------------------------
void
FEMesh::limitCoord(double &x, double &y, double &z)
{
	if (x > 1)
		x = 1;
	if (y > 1)
		y = 1;
	if (z > 1)
		z = 1;

	if (m_nnelem == 4 || m_nnelem == 10) {
		if (x < 0)
			x = 0;
		if (y < 0)
			y = 0;
		if (z < 0)
			z = 0;
		if ((x + y + z) > 1)
			z = 1 - x - y;
	} else {
		if (x < -1)
			x = -1;
		if (y < -1)
			y = -1;
		if (z < -1)
			z = -1;
	}
}
//---------------------------------------------------------------------------
// return the element and local coordinate for a given global coord
// brute force for now, need better search algorithm
int
FEMesh::localElem(double &x, double &y, double &z)
{
	int n, e0, ni;
	double err, err0 = -1;
	point pa[NUM_NODES];

	Point3 p(x, y, z);

	if (m_scache == NULL) {
		Point3 pl0, pl1;
		getLimits(pl0, pl1);
		m_scache = new SCache(pl0, pl1 - pl0);
		m_scache->resize(numElements());

		for (n = 0; n < numElements(); n++) {
			Point3 c;
			double r;
			elemBoundSph(n, c, r);
			r *= 1.01;
			m_scache->addSphere(n, Sphere3(c, r));
		}
	}

	SCache::spset_t eset;

	// test sphere centered on the point
	Sphere3 se(p, 0.001);

	m_scache->intersect(se, eset);

	SCache::spset_t::iterator it;
	for (it = eset.begin(); it != eset.end(); it++) {
		int e = *it;

//	for (e = 0; e < numElements(); e++) {
		p.getCoord(x, y, z);
		if (!localCoord(e, x, y, z))
			return e;

		limitCoord(x, y, z);
		pArray(e, pa);
		m_shape->Local2Global(pa, x, y, z);

		Point3 p1(x, y, z);
		p1 -= p;
		err =  fabs(p1.X()) + fabs(p1.Y()) + fabs(p1.Z());
		if (err < err0 || err0 == -1) {
			err0 = err;
			e0 = e;
		}
	}

	p.getCoord(x, y, z);
	localCoord(e0, x, y, z, &err, &ni);
	limitCoord(x, y, z);

	if (err0 < 1e-3)
		return e0;

	pArray(e0, pa);
//	m_shape->Local2Global(pa, x, y, z);

	fprintf(stderr, "e: %d, err0: %g, err:%g, ni:%d, p0:<%g, %g, %g> p1:<%g, %g, %g>\n",
		e0, err0, err, ni, p.X(), p.Y(), p.Z(), x, y, z);

	return -1;
}
//---------------------------------------------------------------------------
double FEMesh::interpField(int elem, double x, double y, double z, const double *fld)
{
    assert(elem>=0 && elem<m_nelem);
    assert(x>=-1 && x<=1);
    assert(y>=-1 && y<=1);
    assert(z>=-1 && z<=1);

    static double va[NUM_NODES];
    vArray(elem, va, fld);

    return  m_shape->interp(va, x, y, z);
}
//---------------------------------------------------------------------------
Point3 FEMesh::gradField(int elem, double x, double y, double z, const double *fld)
{
    assert(elem>=0 && elem<m_nelem);
    assert(x>=-1 && x<=1);
    assert(y>=-1 && y<=1);
    assert(z>=-1 && z<=1);

    static double fa[NUM_NODES];
    static point pa[NUM_NODES];
    point gr;

    pArray(elem, pa);
    vArray(elem, fa, fld);

    m_shape->grad(pa, fa, gr, x, y, z);

    return Point3(gr);
}
//---------------------------------------------------------------------------
int FEMesh::pArray(int elem, point *pa) const
{
    assert(pa);
    assert(elem>=0 && elem<m_nelem);
    const node *nd=getElem(elem);
    
    int n;
    for(n=0; n<m_nnelem; n++){
        const Point3 &pt=getNode((*nd)[n]);
        (pa[n])[0]=pt.getX();
        (pa[n])[1]=pt.getY();
        (pa[n])[2]=pt.getZ();
    }
    for(; n<NUM_NODES; n++){    // zero the rest
        (pa[n])[0]=0;
        (pa[n])[1]=0;
        (pa[n])[2]=0;
    }
    return 0;
}
//---------------------------------------------------------------------------
int FEMesh::vArray(int elem, double *arr, const double *field)
{
    assert(arr);
    assert(field);
    assert(elem>=0 && elem<m_nelem);
    node *nd=getElem(elem);
    
    int n;
    for(n=0; n<m_nnelem; n++)
	arr[n]=field[(*nd)[n]];

    for(; n<NUM_NODES; n++)    // zero the rest
	arr[n]=0;
    
    return 0;
}
//---------------------------------------------------------------------------
#define RESERVE_EMIN 256 // size: 5*4096
#define RESERVE_NMIN 1024 // size: 3*4096
void FEMesh::reserveSpace(int nelem, int nnodes)
{
    if(nelem>0){
        if(m_nbrs.capacity() < (unsigned)(m_nelem+nelem)) {
            if(nelem<RESERVE_EMIN) nelem=RESERVE_EMIN;
            m_nbrs.reserve(m_nelem+nelem);
            m_elements.reserve(m_nelem+nelem);
        }
    }

    if(nnodes>0){
        if(m_nodes.capacity() < (unsigned)(m_nnodes+nnodes)) {
            if(nnodes<RESERVE_NMIN) nnodes=RESERVE_NMIN;
            m_nodes.reserve(m_nnodes+nnodes);
            m_ninfo.resize(m_nnodes+nnodes);
        }
    }
}
//---------------------------------------------------------------------------
const int *
FEMesh::getElemNeighbors(int el, int &size)
{
	static int nbrs[MAX_ELEM_NBR];

	size = 0;

	if (el < 0 || el >= m_nelem)
		return NULL;

	const node *nd = getElem(el);
	assert(nd);
	int nnbr=0;
	
	for(int n = 0; n < m_nnelem; n++){
		if(mergeNeighborElem((*nd)[n], nbrs, nnbr, MAX_ELEM_NBR))
			return NULL;
	}

	size = nnbr;

	return nbrs;

}
//---------------------------------------------------------------------------
int
FEMesh::getElemNeighbors(int el, int *elem, int size)
{
	static int nbrs[MAX_ELEM_NBR];

	if (el < 0 || el >= m_nelem)
		return -1;

	const node *nd = getElem(el);
	assert(nd);
	int nnbr=0;
	
	for(int n = 0; n < m_nnelem; n++){
		if(mergeNeighborElem((*nd)[n], nbrs, nnbr, MAX_ELEM_NBR))
			return -1;
	}

	if (elem == NULL || size < 0)
		return nnbr;

	for (int n = 0; n < size; n++) 
		elem[n] = nbrs[n];

	return nnbr;
}
//---------------------------------------------------------------------------
int FEMesh::checkNeighbors(int el, const node *nd)
{
    static int nbrs[MAX_ELEM_NBR];
    assert(nd);
    int nnbr=0;

    for(int n=0; n < m_nnelem; n++){
        if(mergeNeighborElem((*nd)[n], nbrs, nnbr, MAX_ELEM_NBR)) {
		fprintf(stderr, "Merge neighbor failed for node %d\n", n);
            return 1;
        }
    }

    for(int n=0; n<nnbr; n++){
        int e2=nbrs[n];
        if(e2==el) continue;
        node *n2=getElem(e2);
        checkNeighbor(el,nd,e2,n2);
    }
    return 0;
}
//---------------------------------------------------------------------------
int FEMesh::addElem(const RElem &el)
{
    for (int n=0; n<m_nnelem; n++)
        if(el.elem[n]<0 || el.elem[n]>=m_nnodes)
	    return 1;

    int e1=m_nelem++;
    const node *n1=&el.elem;

    m_nbrs.resize(m_nelem);
    m_elements.push_back(el);

    addNodeNeighbors(e1,n1);

    return checkNeighbors(e1,n1);
}
//---------------------------------------------------------------------------
int FEMesh::addNode(const Point3 &nd)
{
    m_nodes.push_back(nd);
    m_ninfo.push_back(FENodeInfo());
    m_nnodes++;
    return 0;
}
//---------------------------------------------------------------------------
int FEMesh::replaceElem(int e1, RElem &el)
{
    for(int n=0; n<m_nnelem; n++)
        if(el.elem[n]<0 || el.elem[n]>=m_nnodes) return 1;

    node *n1=getElem(e1);
    FENeighbor &nb1=m_nbrs[e1];

    // clear neighbors
    for(int f1=0; f1<m_nfelem; f1++){
        int e2=nb1.nelem[f1];
        if(e2<0) continue;
        int f2=nb1.nface[f1];
        assert(f2>=0 && f2<m_nfelem);
        FENeighbor &nb2=m_nbrs[e2];
        assert(nb2.nelem[f2]==e1);
        nb2.nelem[f2]=-1;
        nb2.nface[f2]=-1;
        nb1.nelem[f1]=-1;
        nb1.nface[f1]=-1;
    }
    delNodeNeighbors(e1,n1);

    // set new data
    for(int n=0; n<m_nnelem; n++)
        (*n1)[n]=el.elem[n];

    // create neighbors
    checkNeighbors(e1, n1);
    addNodeNeighbors(e1,n1);

    return 0;
}
//---------------------------------------------------------------------------
int FEMesh::deleteElem(int el)
{
    if(el<0 || el>=m_nelem) return 1;

    int last=m_nelem-1;

    node *n1=getElem(last);
    FENeighbor &nb1=m_nbrs[last];

    // clear neighbors of last node
    for(int f1=0; f1<m_nfelem; f1++){
        int e2=nb1.nelem[f1];
        if(e2<0) continue;
        int f2=nb1.nface[f1];
        assert(f2>=0 && f2<m_nfelem);
        FENeighbor &nb2=m_nbrs[e2];
        assert(nb2.nelem[f2]==last);
        nb2.nelem[f2]=-1;
        nb2.nface[f2]=-1;
        nb1.nelem[f1]=-1;
        nb1.nface[f1]=-1;
    }
    delNodeNeighbors(last,n1);

    m_nelem=last;

    if(el!=last){
        RElem re;
        for(int n=0; n<m_nnelem; n++)
            re.elem[n]=(*n1)[n];

        if(replaceElem(el,re)) return 1;
    }

    m_nbrs.pop_back();
    m_elements.pop_back();

    return 0;
}
//---------------------------------------------------------------------------
void FEMesh::saveElements(FILE *f)
{
    assert(f);
    for(int i=0; i<m_nelem; i++){
        node *nd=getElem(i);
        assert(nd);
        fprintf(f,"%d",i+1);
        for(int n=0; n<m_nnelem; n++)
            fprintf(f," %d",(*nd)[n]+1);
        fprintf(f,"\n");
    }
}
//---------------------------------------------------------------------------
void FEMesh::saveNodes(FILE *f)
{
    assert(f);
    for(int i=0; i<m_nnodes; i++){
        const Point3 &pt=getNode(i);
        fprintf(f,"%d %g %g %g\n",i+1,pt.getX(), pt.getY(), pt.getZ());
    }
}
//---------------------------------------------------------------------------
void FEMesh::saveSigmaE(FILE *f, double *sigma)
{
    assert(f && sigma);
    for(int i=0; i<m_nelem; i++)
        fprintf(f,"%d %g\n",i+1, sigma[i]);
}
//---------------------------------------------------------------------------
void FEMesh::saveSigmaN(FILE *f, double *sigma)
{
    assert(f && sigma);
    for(int i=0; i<m_nnodes; i++)
        fprintf(f,"%d %g\n",i+1, sigma[i]);
}
//---------------------------------------------------------------------------
#define BUF_LEN 1024
int FEMesh::save(const char *fn, double *esig, double *nsig)
{
    static char buf[BUF_LEN+1];
    if(fn==0) return 1;

    snprintf(buf,BUF_LEN,"%s.con",fn);
    buf[BUF_LEN]=0;

    FILE *f=fopen(buf,"wt");
    if(f==0) return 1;
    saveElements(f);
    fclose(f);

    snprintf(buf,BUF_LEN,"%s.cor",fn);
    buf[BUF_LEN]=0;

    f=fopen(buf,"wt");
    if(f==0) return 1;
    saveNodes(f);
    fclose(f);

    if (esig) {
        snprintf(buf,BUF_LEN,"%s.esg",fn);
        buf[BUF_LEN]=0;

        f=fopen(buf,"wt");
        if(f==0) return 1;
        saveSigmaE(f, esig);
        fclose(f);
    }

    if (nsig) {
        snprintf(buf,BUF_LEN,"%s.nsg",fn);
        buf[BUF_LEN]=0;

        f=fopen(buf,"wt");
        if(f==0) return 1;
        saveSigmaN(f, nsig);
        fclose(f);
    }

    if (esig == 0 && nsig == 0) return 0;

    snprintf(buf,BUF_LEN,"%s.msh",fn);
    buf[BUF_LEN]=0;

    f=fopen(buf,"wt");
    if(f==0) return 1;
    fprintf(f,"%d\n", m_nnodes);
    saveNodes(f);
    fprintf(f,"%d %d\n", m_nelem, m_nnelem);
    saveElements(f);
    if (esig) {
        fprintf(f,"%d\n", m_nelem);
        saveSigmaE(f, esig);
    }
    if (nsig) {
        fprintf(f,"%d\n", m_nnodes);
        saveSigmaN(f, nsig);
    }

    fclose(f);
    return 0;
}
//---------------------------------------------------------------------------
int FEMesh::offsetMesh(double dx, double dy, double dz)
{
    for(int n=0; n<m_nnodes; n++)
        moveNodeRel(n, dx, dy, dz);
    return 0;
}
//---------------------------------------------------------------------------
int FEMesh::scaleMesh(double sx, double sy, double sz)
{
    for(int n=0; n<m_nnodes; n++)
        scaleNode(n, sx, sy, sz);
    return 0;
}
//---------------------------------------------------------------------------
void FEMesh::moveNodeRel(int n, double dx, double dy, double dz)
{
    Point3 &pt=getNodeRef(n);
    pt+=Point3(dx,dy,dz);
//    setNode(n,pt);
}
//---------------------------------------------------------------------------
void FEMesh::scaleNode(int n, double sx, double sy, double sz)
{
    Point3 &pt=getNodeRef(n);
    pt.setX(pt.getX()*sx);
    pt.setY(pt.getY()*sy);
    pt.setZ(pt.getZ()*sz);
//    setNode(n,pt);
}
//---------------------------------------------------------------------------
void FEMesh::setNode(int nd, const Point3 &pt)
{
    assert(nd>=0 && nd<m_nnodes);
    Point3 &p=m_nodes[nd];
    p=pt;
    markNeighborNodes(nd);
}
//---------------------------------------------------------------------------
void FEMesh::setNode(int nd, double x, double y, double z)
{
    assert(nd>=0 && nd<m_nnodes);
    Point3 &p=m_nodes[nd];
    p=Point3(x,y,z);
    markNeighborNodes(nd);
}
//---------------------------------------------------------------------------
int FEMesh::markBoundaryNodes(void)
{
    int count = 0;

    for(int e=0; e<m_nelem; e++){
        FENeighbor &fn=m_nbrs[e];
        node *nds=getElem(e);
        for(int f=0; f<m_nfelem; f++){
            if(fn.nelem[f]>=0) continue;
            FaceList *fl;
            if (m_nnelem == 20)
	        fl=&(m_facemap20[f]);
            else if (m_nnelem == 8)
		fl=&(m_facemap8[f]);
            else if (m_nnelem == 10)
		fl=&(m_facemap10[f]);
            else
		fl=&(m_facemap4[f]);

            for(int n=0; n<m_nnface; n++){
                int nd=(*fl)[n];
                assert(nd>=0 && nd<m_nnelem);
                nd=(*nds)[nd];
                if(isNodeBound(nd)) continue;
                markNode(nd,FEM_NODE_BOUND);
                count++;
            }
        }
    }
    return count;
}
//---------------------------------------------------------------------------
int FEMesh::markClassBoundary(int cls)
{
    int count = 0;

    for(int e=0; e<m_nelem; e++){
	if (getElemClass(e) != cls) continue;

        FENeighbor &fn=m_nbrs[e];
        node *nds=getElem(e);
        for(int f=0; f<m_nfelem; f++){
            if(fn.nelem[f] >= 0)
		if (getElemClass(fn.nelem[f]) == cls)
			continue;
            FaceList *fl;
            if (m_nnelem == 20)
		fl=&(m_facemap20[f]);
            else if (m_nnelem == 8)
		fl=&(m_facemap8[f]);
            else if (m_nnelem == 10)
		fl=&(m_facemap10[f]);
            else
		fl=&(m_facemap4[f]);

            for(int n=0; n<m_nnface; n++){
                int nd=(*fl)[n];
                assert(nd>=0 && nd<m_nnelem);
                nd=(*nds)[nd];
//                if(isNodeBound(nd)) continue;
                markNode(nd,FEM_NODE_EBCLS);
                count++;
            }
        }
    }
    return count;
}
//---------------------------------------------------------------------------
int FEMesh::getNeighborNodes(int nd, int *nodes, int size)
{
    assert(nodes);
    assert(nd>=0 && nd<m_nnodes);

    int nnbr=0;
    FENodeInfo &fi=m_ninfo[nd];

    for (int n=0; n<fi.numNeighbors(); n++) {
        int el=fi.getNeighbor(n);
        node *nds=getElem(el);
        int i, j;
        for (i=0; i<m_nnelem; i++)
            if ((*nds)[i]==nd) break;
        assert(i < m_nnelem);
        for (int m=0; m<3; m++) {
            int nn;

            if(m_nnelem == 20)
		nn=m_nnbrmap20[i][m];
            else if(m_nnelem == 8)
		nn=m_nnbrmap8[i][m];
            else if(m_nnelem == 10)
		nn=m_nnbrmap10[i][m];
            else
		nn=m_nnbrmap4[i][m];

            if (nn < 0) continue;
            nn=(*nds)[nn];  // map to node id
            for (j=0; j<nnbr; j++)
                if (nodes[j]==nn) break;
            if (j < nnbr) continue;
            if (nnbr == size) return -1;  // full
            nodes[nnbr++]=nn;
        }
    }
    return nnbr;
}
//---------------------------------------------------------------------------
int FEMesh::markNeighborNodes(int nd)
{
    assert(nd>=0 && nd<m_nnodes);

    int nnbr=0;
    FENodeInfo &fi=m_ninfo[nd];

    for (int n=0; n<fi.numNeighbors(); n++) {
        int el=fi.getNeighbor(n);
        node *nds=getElem(el);
        int i;
        for (i=0; i<m_nnelem; i++)
            if ((*nds)[i]==nd) break;
        assert(i < m_nnelem);
        for (int m=0; m<3; m++) {
            int nn;
            if (m_nnelem == 20)
		nn=m_nnbrmap20[i][m];
            else if (m_nnelem == 8)
		nn=m_nnbrmap8[i][m];
            else if (m_nnelem == 10)
		nn=m_nnbrmap10[i][m];
            else
		nn=m_nnbrmap4[i][m];

            if (nn < 0) continue;
            nn=(*nds)[nn];  // map to node id
            markNode(nn, FEM_NODE_NBCHG);
        }
    }
    return nnbr;
}
//---------------------------------------------------------------------------
int FEMesh::getNeighborElem(int nd, int *elem, int size)
{
    assert(nd>=0 && nd<m_nnodes);

    FENodeInfo &fi=m_ninfo[nd];
    int nnbr=fi.numNeighbors();

    if (elem == NULL) 
	    return nnbr;

    for (int n=0; n<nnbr; n++) {
        if(n >= size) return -1; // full
        elem[n]=fi.getNeighbor(n);
    }
    return nnbr;
}
//---------------------------------------------------------------------------
int FEMesh::mergeNeighborElem(int nd, int *elem, int &size, int maxsize)
{
    assert(elem);
    assert(nd>=0 && nd<m_nnodes);

    FENodeInfo &fi=m_ninfo[nd];
    int nnbr=fi.numNeighbors();

    for (int n=0; n<nnbr; n++) {
        int nb=fi.getNeighbor(n);
        int m;
        for (m=0; m<size; m++)
            if (elem[m]==nb) break;
        if (m!=size) continue;
        if (size >= maxsize) return 1;
        elem[size++]=nb;
    }
    return 0;
}
//---------------------------------------------------------------------------
void
FEMesh::getLimits(double &xi, double &yi, double &zi,
		  double &xf, double &yf, double &zf) const
{
	const Point3 &p0=getNode(0);
	xi=xf=p0.getX();
	yi=yf=p0.getY();
	zi=zf=p0.getZ();
	
	for (int n=1; n<m_nnodes; n++) {
		const Point3 &p=getNode(n);
		if (xi > p.getX()) xi=p.getX();
		if (yi > p.getY()) yi=p.getY();
		if (zi > p.getZ()) zi=p.getZ();
		if (xf < p.getX()) xf=p.getX();
		if (yf < p.getY()) yf=p.getY();
		if (zf < p.getZ()) zf=p.getZ();
	}
}
//---------------------------------------------------------------------------
int FEMesh::replaceNode(int n1, int n2)
{
    static int elem[MAX_ELEM_NBR];

    assert(n1>=0 && n1<m_nnodes);
    assert(n2>=0 && n2<m_nnodes);
    assert(m_ninfo[n1].numNeighbors() == 0);

    int nnbr=m_ninfo[n2].numNeighbors();
    assert(nnbr <= MAX_ELEM_NBR);

    for (int n=0; n<nnbr; n++)
        elem[n]=m_ninfo[n2].getNeighbor(n);

    for (int n=0; n<nnbr; n++) {
        RElem re;
        node *nds=getElem(elem[n]);
        for(int m=0; m<m_nnelem; m++) {
            int e=(*nds)[m];
            if (e == n2) e=n1;
            re.elem[m]=e;
        }
        if(replaceElem(elem[n],re)) return 1;
    }

    setNode(n1, getNode(n2));

    return 0;

}
//---------------------------------------------------------------------------
int FEMesh::discardNodes(void)
{
    int start=0;
    int end=m_nnodes-1;

    while(1){
        while(end > start) {
            if (m_ninfo[end].numNeighbors()) break;
            end--;
        }
        while(end > start) {
            if (m_ninfo[start].numNeighbors() == 0) break;
            start++;
        }
        if (start >= end) break;

        if (replaceNode(start, end))
            return 1;
    }
    m_nnodes=end+1;
    return 0;
}
//---------------------------------------------------------------------------
int FEMesh::elemBoundSph(int el, Point3 &c, double &rad) const
{
	if (el < 0 || el > m_nelem) return 1;

	static double *radcache = NULL;

	if (radcache == NULL) {
		radcache = new double[m_nelem];
		for (int n = 0; n < m_nelem; n++)
			radcache[n] = -1;
	}

	// set center
	double x,y,z;
	if (m_nnelem == 4 || m_nnelem == 10)
		x = y = z = 0.25;
	else
		x = y = z = 0;

	globalCoord(el, x, y, z);
	c.setX(x);
	c.setY(y);
	c.setZ(z);

	rad=radcache[el];
	if (rad >= 0)
		return 0;

	int start = (m_nnelem == 4 || m_nnelem == 10) ? 0 : -1;

	for (int i = start; i <= 1; i++) {
		for (int j = start; j <= 1; j++) {
			for (int k = start; k <= 1; k++) {
				if (m_nnelem == 4 && (i + j + k) > 1)
					continue;
				if (m_nnelem == 10 && (i + j + k) > 1)
					continue;
				x=(double)i;
				y=(double)j;
				z=(double)k;
				globalCoord(el,x,y,z);
				Point3 pt(x,y,z);
				pt-=c;
				if (rad < pt.length())
					rad = pt.length();
			}
		}
	}
	radcache[el] = rad;
	
	return 0;
}
//---------------------------------------------------------------------------
void FEMesh::getEdge(int elem, int edge, Point3 &p1, Point3 &p2)
{
	assert(elem>=0 && elem<m_nelem);
	assert(edge>=0 && edge <NUM_EDGE_ELEM);
	node *nd=getElem(elem);
	
	int n1, n2;

	switch(m_nnelem) {
	case 20:
		n1=m_edge20[edge][0];
		n2=m_edge20[edge][1];
		break;
	case 10:
		n1=m_edge10[edge][0];
		n2=m_edge10[edge][1];
		break;
	case 8:
		n1=m_edge8[edge][0];
		n2=m_edge8[edge][1];
		break;
	case 4:
		n1=m_edge4[edge][0];
		n2=m_edge4[edge][1];
		break;
	default:
		assert(0);
	}

	p1 = getNode((*nd)[n1]);
	p2 = getNode((*nd)[n2]);
}
//---------------------------------------------------------------------------
void FEMesh::getEdge(int elem, int edge, int &n1, int &n2)
{
	assert(elem>=0 && elem<m_nelem);
	assert(edge>=0 && edge <numEdgeElem());
	node *nd=getElem(elem);

	switch(m_nnelem) {
	case 20:
		n1=m_edge20[edge][0];
		n2=m_edge20[edge][1];
		break;
	case 10:
		n1=m_edge10[edge][0];
		n2=m_edge10[edge][1];
		break;
	case 8:
		n1=m_edge8[edge][0];
		n2=m_edge8[edge][1];
		break;
	case 4:
		n1=m_edge4[edge][0];
		n2=m_edge4[edge][1];
		break;
	default:
		assert(0);
	}

	n1 = (*nd)[n1];
	n2 = (*nd)[n2];

}
//---------------------------------------------------------------------------
int
FEMesh::elementVolume(int elem, double &vol) const
{
	assert(elem>=0 && elem<m_nelem);

	static point pa[NUM_NODES];
	pArray(elem,pa);

	return m_shape->Volume(pa, vol);
}
//---------------------------------------------------------------------------
int
FEMesh::faceSurfaceArea(int elem, int face, double &sa) const
{
	// map face numbers to surface numbers in shape.cc
	static char face2surf[]={5, 4, 1, 2, 0, 3};

	assert(elem>=0 && elem<m_nelem);
	assert(face>=0 && face<m_nfelem);

	static point pa[NUM_NODES];
	pArray(elem,pa);

	if (m_nnelem == 4 || m_nnelem == 10)
		return m_shape->SurfaceArea(pa, face, sa);
	else
		return m_shape->SurfaceArea(pa, face2surf[face], sa);
}
//---------------------------------------------------------------------------

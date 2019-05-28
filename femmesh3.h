/* $Id: femmesh3.h,v 1.4 2008/11/24 04:12:13 canacar Exp $ */
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
#ifndef femmesh3H
#define femmesh3H
#include <stdio.h>
#include "define.h"
#include "point3.h"
#include "shape.h"
#include "scache.h"
#include <vector>
#include <deque>
using namespace std;
//---------------------------------------------------------------------------

// wrappers for difficult-to-use node/point typedefs
struct RElem{
    int cls;
    node elem;
};

struct FENeighbor{
    FENeighbor();
    FENeighbor(const FENeighbor &nbr);
    virtual ~FENeighbor();
    int nelem[MAX_FACE_ELEM];
    char nface[MAX_FACE_ELEM];
};

// XXX need something more efficient!
#define MAX_ELEM_NBR 8000
// node flags
#define FEM_NODE_BOUND 1        // node fixed to a boundary
#define FEM_NODE_NBCHG 2        // change in neighbor/environment
#define FEM_NODE_ERROR 4        // error occured
#define FEM_NODE_EBCLS 8	// Element class boundary

struct FENodeInfo{
    FENodeInfo(int f=0): flags(f) {}
    FENodeInfo(const FENodeInfo &i);

    FENodeInfo &operator=(const FENodeInfo &i);

    inline int numNeighbors(void) const {return nbrs.size();}
    inline int getNeighbor(int n) const {return nbrs[n]; }
    void addNeighbor(int nbr);
    void delNeighbor(int nbr);

    unsigned flags;
    deque<int> nbrs;
};

class FEMesh{
public:
    FEMesh();
    virtual ~FEMesh(void);

    inline int numNodes(void) const {return m_nnodes;}
    inline int numElements(void) const {return m_nelem;}
    inline int numNodeElem(void) const {return m_nnelem;}
    inline int numFaceElem(void) const {return m_nfelem;}
    inline int numNodeFace(void) const {return m_nnface;}

    inline int numEdgeElem(void) const {
	    return ((m_nnelem == 4 || m_nnelem == 10) ? 6 : 12);
    }

    inline unsigned nodeFlags(int n) const {return m_ninfo[n].flags;}

    inline int getElemClass(int e) const {
	assert (e>=0 && e<m_nelem);
	return m_elements[e].cls;
    }
    inline void setElemClass(int e, int cls) {
	assert (e>=0 && e<m_nelem);
	m_elements[e].cls=cls;
    }

    inline void setNodeFlags(int n, unsigned f) {m_ninfo[n].flags = f;}
    inline void clearNodeFlag(int n, unsigned f) {m_ninfo[n].flags &= ~f;}
    inline void markNode(int n, unsigned f) {m_ninfo[n].flags |= f;}

    inline bool isNodeChange(int n) {return m_ninfo[n].flags & FEM_NODE_NBCHG;}
    inline bool isNodeBound(int n) {return m_ninfo[n].flags & FEM_NODE_BOUND;}
    inline bool isClassBound(int n) {return m_ninfo[n].flags & FEM_NODE_EBCLS;}

    inline FEShape *getShape(void) {return m_shape;}

    inline const char *getMeshFn(void) {return m_meshfn;}

    int findNeighbor(int &ne, int &nf);
    int getNeighborNodes(int nd, int *nodes, int size);
    int getNeighborElem(int nd, int *elem, int size);
    int getElemNeighbors(int el, int *elem, int size);
    const int *getElemNeighbors(int el, int &size);
    int getFaceNodes(int elem, int face, int *nodes) const;
    void globalCoord(int elem, double &x, double &y, double &z) const;
    int localCoord(int elem, double &x, double &y, double &z,
		   double *err = NULL, int *ni = NULL) const;
    int localElem(double &x, double &y, double &z);
    double interpField(int elem, double x, double y, double z,
		       const double *fld);
    Point3 gradField(int elem, double x, double y, double z,
		     const double *fld);
    void elemCoord(int elem, int i, double &x, double &y, double &z);

    int pArray(int elem, point *pa) const;
    int vArray(int elem, double *arr, const double *field);

    const node *getElem(int elem) const;
    node *getElem(int elem);
    const Point3 &getNode(int node) const;
    void getEdge(int elem, int edge, Point3 &p1, Point3 &p2);
    void getEdge(int elem, int edge, int &n1, int &n2);

    void getLimits(double &x0, double &y0, double &z0,
                   double &x1, double &y1, double &z1) const;
    inline void getLimits(Point3 &p0, Point3 &p1) const
        {getLimits(p0.X(), p0.Y(), p0.Z(), p1.X(), p1.Y(), p1.Z());}

    void setNode(int node, const Point3 &pt);
    void setNode(int node, double x, double y, double z);
    void moveNodeRel(int n, double dx, double dy, double dz);
    void scaleNode(int n, double sx, double sy, double sz);
    int offsetMesh(double dx, double dy, double dz);
    inline void offsetMesh(const Point3 &pt) 
	{offsetMesh(pt.getX(), pt.getY(), pt.getZ());}
    int scaleMesh(double dx, double dy, double dz);
    int markBoundaryNodes();
    int markClassBoundary(int cls);

    int addElem(const RElem &el);
    int deleteElem(int el);
    int addNode(const Point3 &nd);
    int replaceElem(int elem, RElem &el);
    void reserveSpace(int nelem, int nnodes);

    int discardNodes(void);

    int loadMesh(const char *name, vector<double> *cmap = NULL);
    int initMesh(int nnelem);

    inline double getSigmaE(int e) {
	assert (e>=0 && e<m_nelem);
	    return m_sigelem[e];
    }
    inline double *getSigmaE(void) {return m_sigelem;}
    inline double *getSigmaN(void) {return m_signode;}
    int save(const char *fn, double *esig=0, double *nsig=0);

//    findNodeNeighbors(void);
    int elemBoundSph(int el, Point3 &c, double &rad) const;
    
    int elementVolume(int el, double &vol) const;
    int faceSurfaceArea(int elem, int face, double &sa) const;

protected:

    static FaceList m_facemap20[];
    static int m_nnbrmap20[][3];
    static int m_edge20[][2];
    static int m_corner[];

    static FaceList m_facemap8[];
    static int m_nnbrmap8[][3];
    static int m_edge8[][2];

    static FaceList m_facemap10[];
    static int m_nnbrmap10[][3];
    static int m_edge10[][2];

    static FaceList m_facemap4[];
    static int m_nnbrmap4[][3];
    static int m_edge4[][2];

    void checkNeighbor(int e1, const node* n1, int e2, const node *n2);
    void addNodeNeighbors(int el, const node *nds);
    void delNodeNeighbors(int el, const node *nds);
    int locateFace(int *nodes);

    Point3 &getNodeRef(int node);

    void saveNodes(FILE *f);
    void saveElements(FILE *f);
    void saveSigmaE(FILE *f, double *sigma);
    void saveSigmaN(FILE *f, double *sigma);

    int mergeNeighborElem(int nd, int *elem, int &size, int maxsize);
    int checkNeighbors(int el, const node *nd);

    int replaceNode(int n1, int n2);

    int markNeighborNodes(int nd);
    void limitCoord(double &x, double &y, double &z);

    int parseSigmaLine(const char *buf, int idx, vector<double> *cmap, double *sig);

    vector<FENeighbor> m_nbrs; // neighbor info
    vector<Point3> m_nodes;
    vector<RElem> m_elements;
    deque<FENodeInfo> m_ninfo;

    SCache *m_scache;
    FEShape *m_shape;
    int m_nnodes;
    int m_nelem;
    int m_nnelem;
    int m_nnface;
    int m_nfelem;
    double *m_signode;
    double *m_sigelem;
    char   *m_meshfn;
};

#endif

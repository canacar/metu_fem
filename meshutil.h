/* $Id: meshutil.h,v 1.3 2008/01/28 07:41:13 canacar Exp $ */
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
/*! \file meshutil.h
    \brief Mesh Utility class header file.
    Contains Mesh utility class decleration.
*/
#ifndef meshUtilH
#define meshUtilH
//---------------------------------------------------------------------------
#include <stdio.h>
#include "define.h"
#define MMGR_BUFSIZE 1024

class Node {
public:
	Node(){}
	virtual ~Node();
	virtual int &node(int i) = 0 ;
};

class Node8:public Node {
public:
	Node8():Node() {}
	virtual ~Node8();
	virtual int &node(int i){
		assert (i>=0 && i<8);
		return _node[i];
	}
private:
	int _node[8];
};

class Node20:public Node {
public:
	Node20():Node() {}
	virtual ~Node20();
	virtual int &node(int i){
		assert (i>=0 && i<20);
		return _node[i];
	}
private:
	int _node[20];
};

class NodeArray {
public:
	NodeArray(int sz):_size(sz) {}
	virtual ~NodeArray();
	inline int size(void) {return _size;}
	virtual Node &operator[](int i) = 0;

private:
	int _size;
};

class Node8Array:public NodeArray {
public:
        Node8Array(int sz):NodeArray(sz){
                assert(sz>0);
                _nodes=new Node8[sz];
        }
        virtual ~Node8Array();
        virtual Node &operator[](int i){
                assert(i>=0 && i<size());
                return _nodes[i];
        }
private:
        Node8 *_nodes;
};

class Node20Array:public NodeArray {
public:
	Node20Array(int sz):NodeArray(sz){
		assert(sz>0);
		_nodes=new Node20[sz];
	}
	virtual ~Node20Array();
	virtual Node &operator[](int i){
		assert(i>=0 && i<size());
		return _nodes[i];
	}
private:
	Node20 *_nodes;
};


enum SourceModelID { SM_J, SM_YAN };

struct DInfo {
	point J; // dipole strength
	point D; // local coords
	SourceModelID smodel; // source model
	int elem; // element
	int dgroup; // dipole group
};

struct Shell {
	double radius;
	double sig;
};

class MeshUtil {
 public:
	static int readListHdr(FILE *f, int *size);
	static int readElemListHdr(FILE *f, int *size, int *nen);
	static int readIntList(FILE *f, int *lst, int start, int size);
	static int readDoubleList(FILE *f, double *lst, int start, int size);
	static int readPointList(FILE *f, point *pts, int start, int size);
	static int readElemList(FILE *f, NodeArray &nds, int start, 
				int size, int nn, int en);
	static int *loadIntList(FILE *f, int *size);
	static int skipLineItems(FILE *f, int start, int size);
	static int skipNodeItems(FILE *f, int start, int size, int en);
	static int checkHeader(FILE *f, const char *sign);
	static int readDipoleList(FILE *f, DInfo *dip, int start,
				  int size, int ne);
	static int readShellList(FILE *f, Shell *sh, int start, int size);

	static int *readSensors(char *fn, int &num, int nn);
	static DInfo *readDipoles(char *fn, int &num, int ne);
	static Shell *readShells(char *fn, int &num);

 private:
	static char m_buf[MMGR_BUFSIZE];
};

#endif

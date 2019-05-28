/* $Id: dmatrix.h,v 1.2 2008/01/28 07:36:28 canacar Exp $ */
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
#ifndef dmatrixH
#define dmatrixH

#include "define.h"
#include <string.h>

//---------------------------------------------------------------------------
/* DFSMatrix: A FULL square matrix of sz entries (0 to sz-1)
 * access elements directly from the _vector variable
 */

struct DFSMatrix{
  DFSMatrix(int sz);
  ~DFSMatrix(){if(_mat) delete[] _mat;}

  void mem(int &alloc, int &used, int &bytes);
  inline double get(int i, int j){
	assert(i>=0 && i<_size && j>=0 && j<_size);
	return _mat[i*_size+j];
  }
  inline void set(int i, int j, double val) {
	assert(i>=0 && i<_size && j>=0 && j<_size);
	_mat[i*_size+j]=val;
  }
  inline void add(int i, int j, double val) {
	assert(i>=0 && i<_size && j>=0 && j<_size);
	_mat[i*_size+j]+=val;
  }
  inline void clear(void){memset(_mat,0,_size*_size*sizeof(double));}
  double *_mat;
  int _size;
};

#endif

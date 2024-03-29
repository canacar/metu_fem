/* $Id: dmatrix.cc,v 1.2 2008/01/28 07:36:28 canacar Exp $ */
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
#ifdef _WIN32_
#include <mem.h>
#endif

#include "dmatrix.h"

DFSMatrix::DFSMatrix(int sz){
  assert(sz>=0);
  if(sz){
	_mat=new double[sz*sz];
	_size=sz;
	clear();
  }
}

void DFSMatrix::mem(int &alloc, int &used, int &bytes){
  int len=_size*_size;
  alloc+=_size;
  for(int n=0; n<len; n++)
    if(_mat[n]) used++;
  bytes+=(len*sizeof(double))+sizeof(DFSMatrix);
}

//---------------------------------------------------------------------------

/* $Id: shape20.h,v 1.2 2008/01/28 07:36:28 canacar Exp $ */
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

#ifndef _SHAPE20_H_
#define _SHAPE20_H_

#include "shape.h"

// abstract base class for FEM element shapes
class FEShape20:public FEShape {
public:
	FEShape20(void);
	virtual ~FEShape20();

protected:

	virtual void Shape(double x, double y, double z);

	virtual int SurfInt(const point *pts, int surf,
		    double &vnx, double &vny, double &vnz);

	virtual int SurfIntR(const point *pts, const point sens, int surf,
                    double &bx, double &by, double &bz);
	void setupVolGP(void);
	void setupSurfGP(int s);

private:
	static int m_surfnodes20[];
};

#endif


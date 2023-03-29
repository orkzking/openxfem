#
# Makefile - To compile OpenXFEM++ with MinGW GCC Compiler
# Copyright (C) 2011  Biplab Kumar Modak (bkmodak@gmail.com)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

BUILD	= debug
CC  	= gcc
CPP 	= g++

ifeq ($(BUILD), release)
	CPPFLAGS	= -Wall -O2
	LDFLAGS 	= -s -static
else
	CPPFLAGS	= -Wall -g
	LDFLAGS 	= -static
endif

SOURCES = 	3disoelement.cpp \
			abssigneddistance.cpp\
			ansplate.cpp \
			auxiliaryfield.cpp \
			bimatcrackasym.cpp \
			boundary.cpp \
			circle.cpp \
			clock.cpp \
			column.cpp \
			constant.cpp \
			constantstiffness.cpp \
			crackgrowthdirectionlaw.cpp \
			crackgrowthincrementlaw.cpp \
			crackinterior.cpp \
			crackjunction.cpp \
			cracktip.cpp \
			deadwght.cpp \
			delaunay.cpp \
			diagmtrx.cpp \
			dictionr.cpp \
			discontinuousfunction.cpp \
			dof.cpp \
			domain.cpp \
			elasticmaterial.cpp \
			element.cpp \
			enrichmentdetector.cpp \
			enrichmentfunction.cpp \
			enrichmentitem.cpp \
			fei2dquadlin.cpp \
			fei2dtrilin.cpp \
			fei2dtriqua.cpp \
			fei3dlineartetra.cpp \
			femcmpnn.cpp \
			fixedincrement.cpp \
			flormtrx.cpp \
			flotarry.cpp \
			freader.cpp \
			freestor.cpp \
			gausspnt.cpp \
			geometrydescription.cpp \
			geometryentity.cpp \
			geometry_2D.cpp \
			geometry_base.cpp \
			hole.cpp \
			homogeneouscrackasym.cpp \
			initial.cpp \
			intarray.cpp \
			integrationrule.cpp \
			junctionfunction.cpp \
			levelsetdescription.cpp \
			lhs.cpp \
			linsyst.cpp \
			list.cpp \
			load.cpp \
			loadtime.cpp \
			main.cpp \
			material.cpp \
			materialinterface.cpp \
			mathfem.cpp \
			mathutil.cpp \
			matrix.cpp \
			matrixops.cpp \
			maxhoopstress.cpp \
			mitc4.cpp \
			modifiedhomocrackasymp.cpp \
			newmark.cpp \
			newtonraphson.cpp \
			nlsolver.cpp \
			node.cpp \
			nodload.cpp \
			peak.cpp \
			piecewis.cpp \
			piecewiselinear.cpp \
			planelast.cpp \
			planeproblem.cpp \
			plateiso.cpp \
			plateiso4.cpp \
			polymtrx.cpp \
			polynoxy.cpp \
			quad_u.cpp \
			samplingcriterion.cpp \
			skyline.cpp \
			splitgaussquadrature.cpp \
			standarddescriptio.cpp \
			standardquadrature.cpp \
			stdafx.cpp \
			stressarray.cpp \
			string.cpp \
			tetra4.cpp \
			timestep.cpp \
			timinteg.cpp \
			tri6.cpp \
			tri_u.cpp \
			vertex.cpp \
			voidenrichfunction.cpp \
			vonmisesmaterial.cpp \
			vonmisesmaterial_h.cpp

OBJS = $(SOURCES:.cpp=.o)

all:	openxfem++.exe

openxfem++.exe : $(OBJS)
	$(CPP) $(LDFLAGS) $(OBJS) -o $@

.c.o:
	$(CPP) -c $(CPPFLAGS) $< -o $@

clean:
	del *.o *.exe
/*******************************************************************************
* CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
* version 0.1                                                                  *
* Copyright (C) 2009-2012, IGG Team, LSIIT, University of Strasbourg           *
*                                                                              *
* This library is free software; you can redistribute it and/or modify it      *
* under the terms of the GNU Lesser General Public License as published by the *
* Free Software Foundation; either version 2.1 of the License, or (at your     *
* option) any later version.                                                   *
*                                                                              *
* This library is distributed in the hope that it will be useful, but WITHOUT  *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License  *
* for more details.                                                            *
*                                                                              *
* You should have received a copy of the GNU Lesser General Public License     *
* along with this library; if not, write to the Free Software Foundation,      *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.           *
*                                                                              *
* Web site: http://cgogn.unistra.fr/                                           *
* Contact information: cgogn@unistra.fr                                        *
*                                                                              *
*******************************************************************************/

#ifndef __TRAVERSOR_CELL_H__
#define __TRAVERSOR_CELL_H__

#include "Topology/generic/dart.h"
#include "Topology/generic/dartmarker.h"
#include "Topology/generic/cellmarker.h"
#include "Topology/generic/traversorGen.h"

namespace CGoGN
{

/// Macro for simple syntax traversal (can be nested)
/// parameters:
///   - cell type
///   - name of iterator variable (automatically declared)
///   - map type
///   - map used
///
/// Ex: FOR_ALL_CELLS(VERTEX, di, myMap) { ... }
///
#define FOR_ALL_CELLS(ORBIT, NAME_ITER, MAP_TYPE, MAP_PARAM) TraversorCell<MAP_TYPE, (ORBIT)>  NAME_ITER_TraversalMacroLocalLoop(MAP_PARAM); \
	for (Dart NAME_ITER = NAME_ITER_TraversalMacroLocalLoop.begin(); NAME_ITER != NAME_ITER_TraversalMacroLocalLoop.end(); NAME_ITER = NAME_ITER_TraversalMacroLocalLoop.next())


template <typename MAP, unsigned int ORBIT>
class TraversorCell //: public Traversor<MAP>
{
private:
	const MAP& m ;

	const AttributeContainer* cont ;
	unsigned int qCurrent ;

	DartMarker* dmark ;
	CellMarker<ORBIT>* cmark ;
	const AttributeMultiVector<Dart>* quickTraversal ;

	Dart current ;
	bool firstTraversal ;

public:
	TraversorCell(const MAP& map, bool forceDartMarker = false, unsigned int thread = 0) ;

	~TraversorCell() ;

	inline Dart begin() ;

	inline Dart end() ;

	inline Dart next() ;

	inline void skip(Dart d);
} ;


template <unsigned int ORBIT>
class TraversorCell<GenericMap, ORBIT> //: public Traversor<GenericMap>
{
private:
	const GenericMap& m ;

	const AttributeContainer* cont ;
	unsigned int qCurrent ;

	DartMarker* dmark ;
	CellMarker<ORBIT>* cmark ;
	const AttributeMultiVector<Dart>* quickTraversal ;

	Dart current ;
	bool firstTraversal ;

public:
	TraversorCell(const GenericMap& map, bool forceDartMarker = false, unsigned int thread = 0) ;

	~TraversorCell() ;

	inline Dart begin() ;

	inline Dart end() ;

	inline Dart next() ;

	inline void skip(Dart d);
} ;


template <typename MAP>
class TraversorV : public TraversorCell<MAP, VERTEX>
{
public:
	TraversorV(const MAP& m, unsigned int thread = 0) : TraversorCell<MAP, VERTEX>(m, false, thread)
	{}
};

template <typename MAP>
class TraversorE : public TraversorCell<MAP, EDGE>
{
public:
	TraversorE(const MAP& m, unsigned int thread = 0) : TraversorCell<MAP, EDGE>(m, false, thread)
	{}
};

template <typename MAP>
class TraversorF : public TraversorCell<MAP, FACE>
{
public:
	TraversorF(const MAP& m, unsigned int thread = 0) : TraversorCell<MAP, FACE>(m, false, thread)
	{}
};

template <typename MAP>
class TraversorW : public TraversorCell<MAP, VOLUME>
{
public:
	TraversorW(const MAP& m, unsigned int thread = 0) : TraversorCell<MAP, VOLUME>(m, false, thread)
	{}
};

} // namespace CGoGN

#include "Topology/generic/traversorCell.hpp"

#endif

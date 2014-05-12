/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 INRIA, USTL, UJF, CNRS, MGH                    *
*                                                                             *
* This library is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This library is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this library; if not, write to the Free Software Foundation,     *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.          *
*******************************************************************************
*                              SOFA :: Framework                              *
*                                                                             *
* Authors: The SOFA Team (see Authors.txt)                                    *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/

#include <sofa/core/loader/VoxelLoader.h>

#include <sofa/core/ObjectFactory.h>
#include <sofa/helper/io/Image.h>
#include <sofa/helper/io/ImageRAW.h>
#include <iostream>
#include <string>
#include <map>
#include <algorithm>

namespace sofa
{

namespace core
{

namespace loader
{

using namespace sofa::defaulttype;
using namespace sofa::core::loader;

VoxelLoader::VoxelLoader()
    :BaseLoader()
    ,positions(initData(&positions,"position","Coordinates of the nodes loaded"))
    ,hexahedra(initData(&hexahedra,"hexahedra","Hexahedra loaded"))
{
}

VoxelLoader::~VoxelLoader()
{
}



void VoxelLoader::addHexahedron(helper::vector< helper::fixed_array<size_type,8> >* pHexahedra,
        size_type p0, size_type p1, size_type p2, size_type p3,
        size_type p4, size_type p5, size_type p6, size_type p7)
{
    addHexahedron(pHexahedra, helper::fixed_array <size_type,8>(p0, p1, p2, p3, p4, p5, p6, p7));
}

void VoxelLoader::addHexahedron(helper::vector< helper::fixed_array<size_type,8> >* pHexahedra, const helper::fixed_array<size_type,8> &p)
{
    pHexahedra->push_back(p);
}

} // namespace loader

} // namespace core

} // namespace sofa


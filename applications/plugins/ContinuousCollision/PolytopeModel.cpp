/******************************************************************************
 *       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
 *                (c) 2006-2011 MGH, INRIA, USTL, UJF, CNRS                    *
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
 *                               SOFA :: Plugins                               *
 *                                                                             *
 * Authors: Ricardo Ortiz <ricardo.ortiz@kitware.com>                          *
 *                                                                             *
 ******************************************************************************/

#define CCD_PLUGIN_POLYTOPE_MODEL

#include "initContinuousCollision.h"
#include "PolytopeModel.hpp"
#include <sofa/core/ObjectFactory.h>

namespace sofa
{

namespace component
{

namespace collision
{

SOFA_DECL_CLASS(Polytope)
SOFA_DECL_CLASS(PolytopeModel)

using namespace sofa::defaulttype;

int PolytopeModelClass = core::RegisterObject("Collision model representing a k-DOP")
#ifndef SOFA_FLOAT
                         .add< TPolytopeModel<defaulttype::Vec3dTypes> >()
#endif
#ifndef SOFA_DOUBLE
                         .add< TPolytopeModel<defaulttype::Vec3fTypes> >()
#endif
                         .addAlias("DopModel")
                         .addAlias("PolytopeModel")
                         ;

#ifndef SOFA_FLOAT
template class SOFA_CONTINUOUS_COLLISION_API TPolytopeModel<defaulttype::Vec3dTypes>;
#endif
#ifndef SOFA_DOUBLE
template class SOFA_CONTINUOUS_COLLISION_API TPolytopeModel<defaulttype::Vec3fTypes>;
#endif
} // namespace collision

} // namespace component

} // namespace sofa


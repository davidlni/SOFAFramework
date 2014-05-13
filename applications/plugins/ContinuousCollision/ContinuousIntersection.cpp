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
 *                               SOFA :: Modules                               *
 *                                                                             *
 * Authors: The SOFA Team and external contributors (see Authors.txt)          *
 *                                                                             *
 * Contact information: contact@sofa-framework.org                             *
 ******************************************************************************/
#include <sofa/helper/system/config.h>
#include <sofa/helper/FnDispatcher.inl>
#include "ContinuousIntersection.h"
#include <sofa/core/ObjectFactory.h>
#include <sofa/core/collision/Intersection.inl>
#include <sofa/helper/proximity.h>
#include <iostream>
#include <algorithm>

namespace sofa
{

namespace component
{

namespace collision
{

using namespace sofa::defaulttype;
using namespace sofa::core::collision;

SOFA_DECL_CLASS(ContinuousIntersection)

int ContinuousIntersectionClass = core::RegisterObject("TODO-ContinuousIntersectionClass")
                                .add< ContinuousIntersection >()
                                ;


ContinuousIntersection::ContinuousIntersection()
{
    intersectors.add<CubeModel,       CubeModel,         ContinuousIntersection> (this);

    intersectors.add<SphereModel,     SphereModel,       ContinuousIntersection> (this);

    intersectors.add<CapsuleModel,CapsuleModel, ContinuousIntersection> (this);
    intersectors.add<CapsuleModel,SphereModel, ContinuousIntersection> (this);

    intersectors.add<OBBModel,OBBModel,ContinuousIntersection>(this);
    intersectors.add<SphereModel,OBBModel, ContinuousIntersection> (this);
    intersectors.add<CapsuleModel,OBBModel,ContinuousIntersection>(this);

    intersectors.add<RigidSphereModel,RigidSphereModel,ContinuousIntersection>(this);
    intersectors.add<SphereModel,RigidSphereModel, ContinuousIntersection> (this);
    intersectors.add<CapsuleModel,RigidSphereModel,ContinuousIntersection>(this);
    intersectors.add<RigidSphereModel,OBBModel,ContinuousIntersection>(this);

    intersectors.add<CapsuleModel,RigidCapsuleModel, ContinuousIntersection> (this);
    intersectors.add<RigidCapsuleModel,RigidCapsuleModel, ContinuousIntersection> (this);
    intersectors.add<RigidCapsuleModel,SphereModel, ContinuousIntersection> (this);
    intersectors.add<RigidCapsuleModel,OBBModel,ContinuousIntersection>(this);
    intersectors.add<RigidCapsuleModel,RigidSphereModel,ContinuousIntersection>(this);


    //IntersectorFactory::getInstance()->addIntersectors(this);
}

/// Return the intersector class handling the given pair of collision models, or NULL if not supported.
ElementIntersector* ContinuousIntersection::findIntersector(core::CollisionModel* object1, core::CollisionModel* object2, bool& swapModels)
{
    return intersectors.get(object1, object2, swapModels);
}


} // namespace collision

} // namespace component

namespace core
{
namespace collision
{
template class SOFA_BASE_COLLISION_API IntersectorFactory<component::collision::ContinuousIntersection>;
}
}

} // namespace sofa


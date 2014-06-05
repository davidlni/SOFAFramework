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
#ifndef SOFA_COMPONENT_COLLISION_CONTINUOUSCONTACTMAPPER_INL
#define SOFA_COMPONENT_COLLISION_CONTINUOUSCONTACTMAPPER_INL

#include "ContinuousContactMapper.h"
#include <sofa/simulation/common/Node.h>
#include <sofa/simulation/common/DeleteVisitor.h>

namespace sofa
{

namespace component
{

namespace collision
{

template < class TCollisionModel, class DataTypes >
void ContinuousContactMapper<TCollisionModel,DataTypes>::cleanup()
{
    if (mapping!=NULL)
    {
        simulation::Node* parent = dynamic_cast<simulation::Node*>(model->getContext());
        if (parent!=NULL)
        {
            simulation::Node::SPtr child = dynamic_cast<simulation::Node*>(mapping->getContext());
            child->detachFromGraph();
            child->execute<simulation::DeleteVisitor>(sofa::core::ExecParams::defaultInstance());
            child.reset();
            mapping.reset();
        }
    }
}

template < class TCollisionModel, class DataTypes >
typename ContinuousContactMapper<TCollisionModel,DataTypes>::MMechanicalState* ContinuousContactMapper<TCollisionModel,DataTypes>::createMapping(const char* name)
{
    if (model==NULL) return NULL;
    simulation::Node* parent = dynamic_cast<simulation::Node*>(model->getContext());
    if (parent==NULL)
    {
        std::cerr << "ERROR: ContinuousContactMapper only works for scenegraph scenes.\n";
        return NULL;
    }
    simulation::Node::SPtr child = parent->createChild(name);
    typename MMechanicalObject::SPtr mstate = sofa::core::objectmodel::New<MMechanicalObject>();
    child->addObject(mstate);
    mstate->useMask.setValue(true);
    mapper = sofa::core::objectmodel::New<mapping::BarycentricMapperMeshTopology<InDataTypes,DataTypes> >(model->getMeshTopology(), (topology::PointSetTopologyContainer*)NULL, &model->getMechanicalState()->forceMask, &mstate->forceMask);
    mapping =  sofa::core::objectmodel::New<MMapping>(model->getMechanicalState(), mstate.get(), mapper);
    child->addObject(mapping);
    return mstate.get();
}

} // namespace collision

} // namespace component

} // namespace sofa

#endif

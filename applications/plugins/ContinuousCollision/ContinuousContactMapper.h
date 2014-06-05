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
#ifndef SOFA_COMPONENT_COLLISION_CONTINUOUSCONTACTMAPPER_H
#define SOFA_COMPONENT_COLLISION_CONTINUOUSCONTACTMAPPER_H

#include <sofa/helper/system/config.h>
#include <sofa/component/mapping/BarycentricMapping.h>
#include <sofa/component/container/MechanicalObject.h>
#include <sofa/component/collision/BaseContactMapper.h>
#include "RTriangleModel.h"

#include <iostream>


namespace sofa
{

namespace component
{

namespace collision
{

using namespace sofa::defaulttype;

/// Base class for all mappers using ContinuousMapping
template < class TCollisionModel, class DataTypes >
class ContinuousContactMapper : public BaseContactMapper<DataTypes>
{
public:
    typedef typename DataTypes::Real Real;
    typedef typename DataTypes::Coord Coord;
    typedef TCollisionModel MCollisionModel;
    typedef typename MCollisionModel::InDataTypes InDataTypes;
    typedef typename MCollisionModel::Topology InTopology;
    typedef core::behavior::MechanicalState< InDataTypes> InMechanicalState;
    typedef core::behavior::MechanicalState<  DataTypes> MMechanicalState;
    typedef component::container::MechanicalObject<DataTypes> MMechanicalObject;
    typedef mapping::BarycentricMapping< InDataTypes, DataTypes > MMapping;
    typedef mapping::TopologyBarycentricMapper<InDataTypes, DataTypes> MMapper;
    MCollisionModel* model;
    typename MMapping::SPtr mapping;
    typename MMapper::SPtr mapper;

    ContinuousContactMapper()
        : model(NULL), mapping(NULL), mapper(NULL)
    {
    }

    void setCollisionModel(MCollisionModel* model)
    {
        this->model = model;
    }

    void cleanup();

    MMechanicalState* createMapping(const char* name="contactPoints");

    void resize(int size)
    {
        if (mapping!=NULL)
        {
            mapper->clear();
            mapping->getMechTo()[0]->resize(size);
        }
    }

    void update()
    {
        if (mapping!=NULL)
        {
            core::BaseMapping* map = mapping.get();
            map->apply(core::MechanicalParams::defaultInstance(), core::VecCoordId::position(), core::ConstVecCoordId::position());
            map->applyJ(core::MechanicalParams::defaultInstance(), core::VecDerivId::velocity(), core::ConstVecDerivId::velocity());
        }
    }

    void updateXfree()
    {
        if (mapping!=NULL)
        {
            core::BaseMapping* map = mapping.get();
            map->apply(core::MechanicalParams::defaultInstance(), core::VecCoordId::freePosition(), core::ConstVecCoordId::freePosition());
            map->applyJ(core::MechanicalParams::defaultInstance(), core::VecDerivId::freeVelocity(), core::ConstVecDerivId::freeVelocity());
        }
    }

    void updateX0()
    {
        if (mapping!=NULL)
        {
            core::BaseMapping* map = mapping.get();
            map->apply(core::MechanicalParams::defaultInstance(), core::VecCoordId::restPosition(), core::ConstVecCoordId::restPosition());
        }
    }
};

/// Mapper for RTriangleModel
template<class DataTypes>
class ContactMapper<RTriangleModel, DataTypes> : public ContinuousContactMapper<RTriangleModel, DataTypes>
{
public:
    typedef typename DataTypes::Real Real;
    typedef typename DataTypes::Coord Coord;
    int addPoint(const Coord& P, int index, Real&)
    {
        int nbt = this->model->getMeshTopology()->getNbTriangles();
        if (index < nbt)
            return this->mapper->createPointInTriangle(P, index, this->model->getMechanicalState()->getX());
        else
        {
            int qindex = (index - nbt)/2;
            int nbq = this->model->getMeshTopology()->getNbQuads();
            if (qindex < nbq)
                return this->mapper->createPointInQuad(P, qindex, this->model->getMechanicalState()->getX());
            else
            {
                std::cerr << "ContactMapper<TriangleMeshModel>: ERROR invalid contact element index "<<index<<" on a topology with "<<nbt<<" triangles and "<<nbq<<" quads."<<std::endl;
                std::cerr << "model="<<this->model->getName()<<" size="<<this->model->getSize()<<std::endl;
                return -1;
            }
        }
    }
    int addPointB(const Coord& P, int index, Real& /*r*/, const Vector3& baryP)
    {

        int nbt = this->model->getMeshTopology()->getNbTriangles();
        if (index < nbt)
            return this->mapper->addPointInTriangle(index, baryP.ptr());
        else
        {
            // TODO: barycentric coordinates usage for quads
            int qindex = (index - nbt)/2;
            int nbq = this->model->getMeshTopology()->getNbQuads();
            if (qindex < nbq)
                return this->mapper->createPointInQuad(P, qindex, this->model->getMechanicalState()->getX());
            else
            {
                std::cerr << "ContactMapper<TriangleMeshModel>: ERROR invalid contact element index "<<index<<" on a topology with "<<nbt<<" triangles and "<<nbq<<" quads."<<std::endl;
                std::cerr << "model="<<this->model->getName()<<" size="<<this->model->getSize()<<std::endl;
                return -1;
            }
        }
    }

    inline int addPointB(const Coord& P, int index, Real& r ){return addPoint(P,index,r);}

};

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(CCD_PLUGIN_CONTINUOUS_CONTACT_MAPPER)
extern template class SOFA_CONTINUOUS_COLLISION_API ContactMapper<RTriangleModel>;
#endif

} // namespace collision

} // namespace component

} // namespace sofa

#endif

#ifndef SOFA_COMPONENT_COLLISION_POLYTOPEMODEL_H
#define SOFA_COMPONENT_COLLISION_POLYTOPEMODEL_H
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

#include <sofa/core/CollisionModel.h>
#include <sofa/component/container/MechanicalObject.h>
#include <sofa/defaulttype/Vec3Types.h>

#include "DiscreteOrientedPolytope.h"

namespace sofa
{

namespace component
{

namespace collision
{

using namespace sofa::defaulttype;

template<typename TDataTypes, size_t NumberOfPlanes>
class TPolytopeModel;

template<typename TDataTypes, size_t K>
class Polytope : public core::TCollisionElementIterator<PolytopeModel<TDataTypes,K> >
{
public:
  Polytope(PolytopeModel<TDataTypes,K>* model=NULL, int index=0);

    explicit Polytope(const core::CollisionElementIterator& i);

    const DiscreteOrientedPolytope<Vec3Types::Real,K>& getDOP() const;
    const std::pair<Polytope<TDataTypes,K>,Polytope<TDataTypes,K> >& subcells() const;
};

template<typename TDataTypes, size_t NumberOfPlanes = 18>
class SOFA_BASE_COLLISION_API TPolytopeModel : public core::CollisionModel
{
public:
    enum { K = NumberOfPlanes };
    SOFA_CLASS(SOFA_TEMPLATE2(TPolytopeModel,TDataTypes,K),sofa::core::CollisionModel);

public:
    typedef core::CollisionElementIterator ChildIterator;
    typedef TDataTypes DataTypes;
    typedef typename DataTypes::Real Real;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename DOPType::DistanceArrayType DistanceArrayType;
    typedef Polytope<TDataTypes,K> Element;
    typedef DiscreteOrientedPolytope<Real,K> DOPType;
    friend class Element;

    struct PolytopeData : public DOPType
    {
        std::pair<Element,Element> subcells;
        std::pair<core::CollisionElementIterator,core::CollisionElementIterator> children; ///< Note that children is only meaningfull if subcells in empty
    };

    struct PolytopeSortPredicate
    {
        int axis;
        PolytopeSortPredicate(int _axis) : axis(_axis) {}
        bool operator()(const PolytopeData& c1,const PolytopeData& c2) const
        {
            return c1.GetCenter(axis) < c2.GetCenter(axis);
        }
    };

protected:
    TPolytopeModel();

public:
    virtual void resize(int size);

    void setParentOf(int childIndex, const Vector3& p);
    void setLeafCube(int cubeIndex, int childIndex);
    void setLeafCube(int cubeIndex, std::pair<core::CollisionElementIterator,core::CollisionElementIterator> children, const Vector3& p);


    unsigned int getNumberCells() {
        return elems.size();
    }

    void getBoundingTree ( sofa::helper::vector< std::pair< Vector3, Vector3> > &bounding )
    {
        bounding.resize(elems.size());
        for (unsigned int index=0; index<elems.size(); index++)
        {
            bounding[index] = std::make_pair( elems[index].minBBox, elems[index].maxBBox);
        }
    }

    int getLeafIndex(int index) const
    {
        return elems[index].children.first.getIndex();
    }

    int getLeafEndIndex(int index) const
    {
        return elems[index].children.second.getIndex();
    }

    const PolytopeData & getPolytopeData(int index) const {
        return elems[index];
    }

    // -- CollisionModel interface

    /**
      *Here we make up the hierarchy (a tree) of bounding boxes which contain final CollisionElements like Spheres or Triangles.
      *The leafs of the tree contain final CollisionElements. This hierarchy is made up from the top to the bottom, i.e., we begin
      *to compute a bounding box containing all CollisionElements, then we divide this big bounding box into two boxes.
      *These new two boxes inherit from the root box and have depth 1. Then we can do the same operation for the new boxes.
      *The division is done only if the box contains more than 4 final CollisionElements and if the depth doesn't exceed
      *the max depth. The division is made along an axis. This axis corresponds to the biggest dimension of the current bounding box.
      *Note : a bounding box is a k-DOP here.
      */
    virtual void computeBoundingTree(int maxDepth=0);

    virtual std::pair<core::CollisionElementIterator,core::CollisionElementIterator> getInternalChildren(int index) const;

    virtual std::pair<core::CollisionElementIterator,core::CollisionElementIterator> getExternalChildren(int index) const;

    virtual bool isLeaf( int index ) const;

    void draw(const core::visual::VisualParams*,int index);

    void draw(const core::visual::VisualParams* vparams);

    int addPolytope(Element subcellsBegin, Element subcellsEnd);
    void updatePolytope(PolytopeData &element);
    void updatePolytope(int index);
    void updatePolytopes();

protected:
    sofa::helper::vector<PolytopeData> elems;
    sofa::helper::vector<int> parentOf; ///< Given the index of a child leaf element, store the index of the parent cube

};

template<typename TDataTypes, size_t K>
inline Polytope<TDataTypes,K>::Polytope(PolytopeModel<K>* model, int index)
    : core::TCollisionElementIterator<PolytopeModel<K> >(model, index)
{}

template<typename TDataTypes, size_t K>
inline Polytope<TDataTypes,K>::Polytope(const core::CollisionElementIterator& i)
    : core::TCollisionElementIterator<PolytopeModel<K> >(static_cast<PolytopeModel<K>*>(i.getCollisionModel()), i.getIndex())
{}


template<typename TDataTypes, size_t K>
inline const DiscreteOrientedPolytope<Vec3Types::Real,K>& Polytope<TDataTypes,K>::getDOP() const
{
    return this->model->elems[index];
}

template<typename TDataTypes, size_t K>
inline const std::pair<Polytope<TDataTypes,K>,Polytope<TDataTypes,K> >& Polytope<TDataTypes,K>::subcells() const
{
    return this->model->elems[index].subcells;
}

} // namespace collision

} // namespace component

} // namespace sofa

#endif

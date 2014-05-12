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

template<typename TDataTypes, size_t K = 18>
class TPolytope : public core::TCollisionElementIterator<TPolytopeModel<TDataTypes,K> >
{
public:
    TPolytope(TPolytopeModel<TDataTypes,K>* model=NULL, size_t index=0);

    explicit TPolytope(const core::CollisionElementIterator& i);

    const std::pair<TPolytope<TDataTypes,K>,TPolytope<TDataTypes,K> >& subcells() const;
    const typename TPolytopeModel<TDataTypes,K>::PolytopeData& getElement() const;
};

template<typename TDataTypes, size_t NumberOfPlanes = 18>
class TPolytopeModel : public core::CollisionModel
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

    typedef DiscreteOrientedPolytope<Real,K> DOPType;
    typedef typename DOPType::DistanceArrayType DistanceArrayType;
    typedef TPolytope<DataTypes,K> Element;
    friend class TPolytope<DataTypes,K>;

    struct PolytopeData : public DOPType
    {
        std::pair<Element,Element> subcells;
        std::pair<ChildIterator,ChildIterator> children;
    };

    struct PolytopeSortPredicate
    {
        PolytopeSortPredicate(const size_t &_axis, const Real &_center) : axis(_axis), center(_center) {}
        bool operator()(const PolytopeData& c1) const
        {
            return c1.GetCenter(axis) < center;
        }
        const size_t &axis;
        const Real &center;
    };

    typedef typename sofa::helper::vector<PolytopeData> ElementListType;

protected:
    TPolytopeModel();

public:
    virtual void resize(size_t size);

    void setParentOf(size_t childIndex, const Vector3& p);

    // For continuous collision
    void enlarge(size_t childIndex, const Real &distance);
    void setLeafPolytope(size_t polytopeIndex, size_t childIndex);
    void setLeafPolytope(size_t polytopeIndex, std::pair<ChildIterator,ChildIterator> children, const Vector3& p);

    size_t getNumberCells() {
        return polytopes.size();
    }

    void getBoundingTree ( sofa::helper::vector< PolytopeData > &bounding )
    {
        bounding = polytopes;
    }

    size_t getLeafIndex(size_t index) const
    {
        return polytopes[index].children.first.getIndex();
    }

    size_t getLeafEndIndex(size_t index) const
    {
        return polytopes[index].children.second.getIndex();
    }

    const PolytopeData & getPolytopeData(size_t index) const {
        return polytopes[index];
    }

    PolytopeData & getPolytopeData(size_t index)  {
        return polytopes[index];
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
    virtual void computeBoundingTree(size_t maxDepth=0);

    virtual std::pair<ChildIterator,ChildIterator> getInternalChildren(size_t index) const;
    virtual std::pair<ChildIterator,ChildIterator> getExternalChildren(size_t index) const;

    virtual bool isLeaf( size_t index ) const;

    void draw(const core::visual::VisualParams*,size_t index);

    void draw(const core::visual::VisualParams* vparams);

    size_t addPolytope(Element subcellsBegin, Element subcellsEnd);
    void updatePolytope(PolytopeData &element);
    void updatePolytope(size_t index);
    void updatePolytopes();

protected:
    ElementListType polytopes;
    sofa::helper::vector<size_t> parentOf; ///< Given the index of a child leaf element, store the index of the parent cube

};

template<typename TDataTypes, size_t K>
inline TPolytope<TDataTypes,K>::TPolytope(TPolytopeModel<TDataTypes,K>* model, size_t index)
    : core::TCollisionElementIterator<TPolytopeModel<TDataTypes,K> >(model, index)
{}

template<typename TDataTypes, size_t K>
inline TPolytope<TDataTypes,K>::TPolytope(const core::CollisionElementIterator& i)
    : core::TCollisionElementIterator<TPolytopeModel<TDataTypes,K> >(static_cast<TPolytopeModel<TDataTypes,K>*>(i.getCollisionModel()), i.getIndex())
{}


template<typename TDataTypes, size_t K>
inline const typename TPolytopeModel<TDataTypes,K>::PolytopeData& TPolytope<TDataTypes,K>::getElement() const
{
    return this->model->polytopes[this->index];
}

template<typename TDataTypes, size_t K>
inline const std::pair<TPolytope<TDataTypes,K>,TPolytope<TDataTypes,K> >& TPolytope<TDataTypes,K>::subcells() const
{
    return this->getElement().subcells;
}

typedef TPolytopeModel<Vec3Types> PolytopeModel;
typedef TPolytope<Vec3Types> Polytope;

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(CCD_PLUGIN_POLYTOPE_MODEL)
#ifndef SOFA_FLOAT
extern template class SOFA_CONTINUOUS_COLLISION_API TPolytopeModel<defaulttype::Vec3dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_CONTINUOUS_COLLISION_API TPolytopeModel<defaulttype::Vec3fTypes>;
#endif
#endif

} // namespace collision

} // namespace component

} // namespace sofa

#endif

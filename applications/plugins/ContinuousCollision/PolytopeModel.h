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

    const typename TPolytopeModel<TDataTypes,K>::PolytopeData& getElement() const;
    const std::pair<TPolytope<TDataTypes,K>,TPolytope<TDataTypes,K> >& subcells() const;
    bool overlaps(const TPolytope<TDataTypes,K> &other) const;
    const typename TPolytopeModel<TDataTypes,K>::DOPType &getBox() const;
};

template<typename TDataTypes, size_t NumberOfPlanes = 18>
class TPolytopeModel : public core::CollisionModel
{
public:
    enum { K = NumberOfPlanes };
    SOFA_CLASS(SOFA_TEMPLATE2(TPolytopeModel,TDataTypes,K),sofa::core::CollisionModel);
    typedef core::topology::BaseMeshTopology::index_type index_type;

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
      PolytopeSortPredicate(const index_type &_axis) : axis(_axis) {}
        bool operator()(const PolytopeData& c1, const PolytopeData& c2) const
        {
            return c1.GetCenter(axis) < c2.GetCenter(axis);
        }
        const index_type &axis;
    };

    typedef typename sofa::helper::vector<PolytopeData> ElementListType;

protected:
    TPolytopeModel();

public:
    virtual void resize(size_t size);

    void cleanPolytopes();

    void setParentOf(const index_type &childIndex, const Vector3& p);

    // For continuous collision
    void enlarge(const index_type &childIndex, const Real &distance);
    void setLeafPolytope(const index_type &polytopeIndex, const index_type &childIndex);
    void setLeafPolytope(const index_type &polytopeIndex, const std::pair<ChildIterator,ChildIterator> &children, const Vector3& p);

    index_type getNumberCells() {
        return polytopes.size();
    }

    void getBoundingTree ( sofa::helper::vector< PolytopeData > &bounding )
    {
        bounding = polytopes;
    }

    index_type getLeafIndex(const index_type &index) const
    {
        return polytopes[index].children.first.getIndex();
    }

    index_type getLeafEndIndex(const index_type &index) const
    {
        return polytopes[index].children.second.getIndex();
    }

    const PolytopeData & getPolytopeData(const index_type &index) const {
        return polytopes[index];
    }

    PolytopeData & getPolytopeData(const index_type &index)  {
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

    size_t addPolytope(const Element &subcellsBegin, const Element &subcellsEnd);
    void updatePolytope(PolytopeData &element);
    void updatePolytope(const index_type &index);
    void updatePolytopes();

protected:
    ElementListType polytopes;
    sofa::helper::vector<index_type> parentOf; ///< Given the index of a child leaf element, store the index of the parent cube

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

template<typename TDataTypes, size_t K>
inline const typename TPolytopeModel<TDataTypes,K>::DOPType& TPolytope<TDataTypes,K>::getBox() const
{
  return this->model->polytopes[this->index];
}

template<typename TDataTypes, size_t K>
inline bool TPolytope<TDataTypes,K>::overlaps(const TPolytope<TDataTypes,K> &other) const
{
  return this->getElement().overlaps(other.getElement());
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

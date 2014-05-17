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
#ifndef SOFA_COMPONENT_COLLISION_TRIANGLEMODEL_H
#define SOFA_COMPONENT_COLLISION_TRIANGLEMODEL_H

#include <sofa/core/CollisionModel.h>
#include <sofa/component/collision/LocalMinDistanceFilter.h>
#include <sofa/component/container/MechanicalObject.h>
#include <sofa/component/topology/TopologyData.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/component/collision/PointModel.h>
#include <CGALPlugin/DecimateMesh.h>
#include <map>
#include "DiscreteOrientedPolytope.h"

namespace sofa
{

namespace component
{

namespace collision
{

using namespace sofa::defaulttype;

template<class DataTypes>
class TRTriangleModel;

class TriangleLocalMinDistanceFilter;

template<class TDataTypes>
class TREdge : public core::TCollisionElementIterator< TRTriangleModel<TDataTypes> >
{
public:
  typedef sofa::core::topology::BaseMeshTopology::index_type index_type;
  typedef std::pair<index_type,index_type> index_pair_type;
  typedef TDataTypes DataTypes;
  typedef typename DataTypes::Coord Coord;
  typedef typename DataTypes::Deriv Deriv;
  typedef TREdgeModel<DataTypes> ParentModel;
  typedef typename DataTypes::Real Real;

  TREdge(ParentModel* model, index_type index);
  TREdge() {}
  explicit TREdge(const core::CollisionElementIterator& i);

  const Coord& p1() const;
  const Coord& p2() const;
  const Coord& p(index_type i) const;

  const index_type &p1Index() const;
  const index_type &p2Index() const;
  const index_type &Index(int i) const;

  const Coord& p1Free() const;
  const Coord& p2Free() const;

  const index_type &e1Index() const;

  const Coord& operator[](int i) const;

  index_type covertices(const TREdge &other, helper::vector<index_pair_type> &verticesPairs);
  bool covertices(const TREdge &other);
  bool coedge(const TREdge &other);

  const Deriv& v1() const;
  const Deriv& v2() const;
  const Deriv& v(int i) const;

  const Deriv& n() const;
  Deriv& n();

  /// Return true if the element stores a free position vector
  bool hasFreePosition() const;

  int flags() const;

  TREdge& shape() {
    return *this;
  }
  const TREdge& shape() const {
    return *this;
  }

  Coord interpX(defaulttype::Vec<2,Real> bary) const
  {
    return (p1()*(1-bary[0]-bary[1])) + (p2()*bary[0]) + (p3()*bary[1]);
  }

};

template<class TDataTypes>
class TRVertex : public core::TCollisionElementIterator< TRTriangleModel<TDataTypes> >
{
public:
  typedef sofa::core::topology::BaseMeshTopology::index_type index_type;
  typedef std::pair<index_type,index_type> index_pair_type;
  typedef TDataTypes DataTypes;
  typedef typename DataTypes::Coord Coord;
  typedef typename DataTypes::Deriv Deriv;
  typedef TRTriangleModel<DataTypes> ParentModel;
  typedef typename DataTypes::Real Real;

  TRVertex(ParentModel* model, index_type index);
  TRVertex() {}
  explicit TRVertex(const core::CollisionElementIterator& i);

  const Coord& p() const;
  const index_type &Index() const;
  const Coord& pFree() const;
  const Coord& operator[](int i) const;

  index_type covertices(const TRVertex &other, helper::vector<index_pair_type> &vertexPair);
  bool covertices(const TRVertex &other);
  bool coedge(const TRVertex &other);

  const Deriv& v() const;

  /// Return true if the element stores a free position vector
  bool hasFreePosition() const;

  int flags() const;

  TRVertex& shape() {
    return *this;
  }
  const TRVertex& shape() const {
    return *this;
  }

};

template<class TDataTypes>
class TRTriangle : public core::TCollisionElementIterator< TRTriangleModel<TDataTypes> >
{
public:
    typedef sofa::core::topology::BaseMeshTopology::index_type index_type;
    typedef std::pair<index_type,index_type> index_pair_type;
    typedef TDataTypes DataTypes;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef TRTriangleModel<DataTypes> ParentModel;
    typedef typename DataTypes::Real Real;

    TRTriangle(ParentModel* model, index_type index);
    TRTriangle() {}
    explicit TRTriangle(const core::CollisionElementIterator& i);

    const Coord& p1() const;
    const Coord& p2() const;
    const Coord& p3() const;
    const Coord& p(index_type i) const;

    const index_type &p1Index() const;
    const index_type &p2Index() const;
    const index_type &p3Index() const;
    const index_type &Index(int i) const;

    const Coord& p1Free() const;
    const Coord& p2Free() const;
    const Coord& p3Free() const;

    const index_type &e1Index() const;
    const index_type &e2Index() const;
    const index_type &e3Index() const;

    const Coord& operator[](int i) const;

    index_type covertices(const TRTriangle &other, helper::vector<index_pair_type> &verticesPairs);
    bool covertices(const TRTriangle &other);
    bool coedge(const TRTriangle &other);

    const Deriv& v1() const;
    const Deriv& v2() const;
    const Deriv& v3() const;
    const Deriv& v(int i) const;

    const Deriv& n() const;
    Deriv& n();

    /// Return true if the element stores a free position vector
    bool hasFreePosition() const;

    int flags() const;

    TRTriangle& shape() {
        return *this;
    }
    const TRTriangle& shape() const {
        return *this;
    }

    Coord interpX(defaulttype::Vec<2,Real> bary) const
    {
        return (p1()*(1-bary[0]-bary[1])) + (p2()*bary[0]) + (p3()*bary[1]);
    }
};

template<class TDataTypes>
class TRTriangleModel : public core::CollisionModel
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(TRTriangleModel, TDataTypes), core::CollisionModel);

    typedef sofa::core::topology::BaseMeshTopology::index_type index_type;
    typedef std::pair<index_type,index_type> index_pair_type;
    typedef TDataTypes DataTypes;
    typedef DataTypes InDataTypes;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef TRTriangle<DataTypes> TriangleElement;
    typedef TREdge<DataTypes> EdgeElement;
    typedef TRVertex<DataTypes> VertexElement;
//     typedef TRTriangle<DataTypes> Element;
    friend class TRTriangle<DataTypes>;

    enum TriangleFlag
    {
        FLAG_P1  = 1<<0, ///< Point 1  is attached to this triangle
        FLAG_P2  = 1<<1, ///< Point 2  is attached to this triangle
        FLAG_P3  = 1<<2, ///< Point 3  is attached to this triangle
        FLAG_E23 = 1<<3, ///< Edge 2-3 is attached to this triangle
        FLAG_E31 = 1<<4, ///< Edge 3-1 is attached to this triangle
        FLAG_E12 = 1<<5, ///< Edge 1-2 is attached to this triangle
        FLAG_BE23 = 1<<6, ///< Edge 2-3 is attached to this triangle and is a boundary
        FLAG_BE31 = 1<<7, ///< Edge 3-1 is attached to this triangle and is a boundary
        FLAG_BE12 = 1<<8, ///< Edge 1-2 is attached to this triangle and is a boundary
        FLAG_POINTS  = FLAG_P1|FLAG_P2|FLAG_P3,
        FLAG_EDGES   = FLAG_E12|FLAG_E23|FLAG_E31,
        FLAG_BEDGES  = FLAG_BE12|FLAG_BE23|FLAG_BE31,
    };

    enum { NBARY = 2 };

    Data<bool> bothSide; // to activate collision on both side of the triangle model

    struct EdgeFeature : public sofa::core::topology::BaseMeshTopology::Edge
    {
        typedef sofa::core::topology::BaseMeshTopology::index_type index_type;
    private:
        index_type FacesIds[2];

        void set(const index_type &id0, const index_type &id1)
        {
            if (id0 > id1)
            {
                this->elems[0] = id0;
                this->elems[1] = id1;
            }
            else
            {
                this->elems[1] = id0;
                this->elems[0] = id1;
            }
        }

    public:
        EdgeFeature()
        {
            this->elems[0] = this->elems[1] = std::numeric_limits<index_type>::max();
            this->FacesIds[0] = this->FacesIds[1] = std::numeric_limits<index_type>::max();
        }

        EdgeFeature(const index_type &id0, const index_type &id1, const index_type &fid)
            : sofa::core::topology::BaseMeshTopology::Edge(id0,id1)
        {
            this->FacesIds[0] = fid;
            this->FacesIds[1] = std::numeric_limits<index_type>::max();
        }

        EdgeFeature(const index_type &id0, const index_type &id1)
            : sofa::core::topology::BaseMeshTopology::Edge(id0,id1)
        {
            this->FacesIds[0] = 0;
            this->FacesIds[1] = std::numeric_limits<index_type>::max();
        }

        void SetFaceId(const index_type &id)
        {
            this->FacesIds[1] = id;
        }

        index_type &VertexId(const index_type &i)
        {
            return this->elems[i];
        }

        index_type const &VertexId(const index_type &i) const
        {
            return this->elems[i];
        }

        index_type &FaceId(const index_type &i)
        {
            return this->FacesIds[i];
        }

        index_type const &FaceId(const index_type &i) const
        {
            return this->FacesIds[i];
        }

        bool operator == (const EdgeFeature &other) const
        {
            return (this->elems[0] == other.elems[0] && this->elems[1] == other.elems[1]);
        }

        bool operator < (const EdgeFeature &other) const
        {
            if (this->elems[0] == other.elems[0])
                return this->elems[1] < other.elems[1];
            else
                return this->elems[0] < other.elems[0];
        }

        bool Covertex(const EdgeFeature &other)
        {
            for(size_t i = 0; i < 2; ++i)
                for(size_t j = 0; j < 2; ++j)
                    if(this->elems[i] == other.VertexId(j))
                        return true;

            return false;
        }
    };

    struct TriangleEdgesFeature : public sofa::core::topology::BaseMeshTopology::EdgesInTriangle
    {
        typedef sofa::core::topology::BaseMeshTopology::index_type index_type;

    public:
        TriangleEdgesFeature()
        {
            this->elems[0] = this->elems[1] = this->elems[2] = std::numeric_limits<index_type>::max();;
        }

        TriangleEdgesFeature(const index_type &id0, const index_type &id1, const index_type &id2)
            : sofa::core::topology::BaseMeshTopology::EdgesInTriangle(id0,id1,id2)
        {
        }

        void set(const index_type &id0, const index_type &id1, const index_type &id2)
        {
            this->elems[0] = id0;
            this->elems[1] = id1;
            this->elems[2] = id2;
        }

        index_type const &EdgeId(const index_type &i) const {
            return this->elems[i];
        }

        index_type const &operator[](const index_type &i) const {
            return this->elems[i];
        }

        index_type &EdgeId(const index_type &i) {
            return this->elems[i];
        }
    };

    struct VertexTriangleFeatures : public helper::vector<index_type> {};

    struct AdjacentPair
    {
        typedef sofa::core::topology::BaseMeshTopology::index_type index_type;
    private:
        index_type Id[2];
        char State[2];
        char Status;

    public:
        AdjacentPair(const index_type &id1, const index_type &id2, const char &state1, const char &state2, const char &status)
        {
            this->Id[0] = id1;
            this->Id[1] = id2;
            this->State[0] = state1;
            this->State[1] = state2;
            this->Status = status;
        }

        void GetParameters(index_type &id1, index_type &id2, char &state1, char &state2)
        {
          id1 = this->Id[0];
          id2 = this->Id[1];
          state1 = this->State[0];
          state2 = this->State[1];
        }

        void GetParameters(index_type &id1, index_type &id2, char &state1, char &state2, char &status)
        {
            this->GetParameters(id1,id2,state1,state2);
            status = this->Status;
        }

        bool operator == (const AdjacentPair &other) const
        {
            return  (this->Id[0] == other.Id[0] && this->Id[1] == other.Id[1]);
        }

        bool operator = (const AdjacentPair &other) const
        {
          this->Id[0] == other.Id[0];
          this->Id[1] == other.Id[1];
          this->State[0] == other.State[0];
          this->State[1] == other.State[1];
          this->Status == other.Status;

        }

        bool operator < (const AdjacentPair &other) const
        {
            if (this->Id[0] == other.Id[0])
                return this->Id[1] < other.Id[1];
            else
                return this->Id[0] < other.Id[0];
        }
    };

    struct NonAdjacentPair
    {
        typedef sofa::core::topology::BaseMeshTopology::index_type index_type;
    private:
        index_type Id[2];

    public:
        NonAdjacentPairs(const index_type &id1, const index_type &id2)
        {
            if (id1 > id2) {
                Id[0] = id1;
                Id[1] = id2;
            } else {
                Id[0] = id2;
                Id[1] = id1;
            }
        }

        void
        GetParameters(index_type &id1, index_type &id2)
        {
            id1 = Id[0];
            id2 = Id[1];
        }

        bool operator < (const NonAdjacentPairs &other) {
            if (Id[0] == other.Id[0])
                return Id[1] < other.Id[1];
            else
                return Id[0] < other.Id[0];
        }
    };

    struct FeaturePair : public helper::fixed_array<2,index_type>
    {
      EdgeFeaturePair(const index_type &i, const index_type &j) : helper::fixed_array<2,index_type>(i,j){}

      void GetParameters(index_type &i, index_type &j)
      {
        i = this->elems[0];
        j = this->elems[1];
      }

      bool operator == (const EdgeFeaturePair &other) const
      {
        return  (this->elems[0] == other.elems[0] && this->elems[1] == other.elems[1]);
      }

      bool operator < (const EdgeFeaturePair &other) const
      {
        if (this->elems[0] == other.elems[0])
          return this->elems[1] < other.elems[1];
        else
          return this->elems[0] < other.elems[0];
      }
    };


protected:
    VecDeriv normals;

    const sofa::core::topology::BaseMeshTopology::SeqTriangles* triangles;
    sofa::core::topology::BaseMeshTopology::SeqTriangles mytriangles;
    helper::vector<EdgeFeature> edgeFeatures;
    helper::vector<TriangleEdgesFeature> triangleEdgeFeatures;
    helper::vector<VertexTriangleFeatures> vertexTriangleFeatures;
    helper::vector<DiscreteOrientedPolytope<Real> > edgeBoxes;
    helper::vector<DiscreteOrientedPolytope<Real> > vertexBoxes;
    helper::vector<DiscreteOrientedPolytope<Real> > faceBoxes;
    helper::vector<AdjacentPair> adjacentPairs[2];
    helper::vector<NonAdjacentPair> nonAdjacentPairs;
    helper::set<FeaturePair> edgeEdgeFeaturePairs;
    helper::set<FeaturePair> vertexFaceFeaturePairs;

    bool needsUpdate;
    virtual void updateFromTopology(double dt);
    virtual void updateFlags(int ntri=-1);
    virtual void updateNormals();
    int getTriangleFlags(int i);

    core::behavior::MechanicalState<DataTypes>* mstate;
    Data<bool> computeNormals;
    int meshRevision;

    sofa::core::topology::BaseMeshTopology* _topology;

    PointModel* mpoints;

    TriangleLocalMinDistanceFilter *m_lmdFilter;

protected:

    TRTriangleModel();
public:
    virtual void init();

    // -- CollisionModel interface

    virtual void resize(int size);

    virtual void computeContinuousBoundingTree(double dt, size_t maxDepth=0);
    void updateFeatureBoxes(double dt);
    void bufferAdjacentLists();
    index_type covertexFace(const index_pair_type &triangleIndices,
                            index_pair_type &vertexPair);
    char getStatus1(const index_pair_type &triangleIndices,
                    index_pair_type &vertexPair);
    char getStatus2(const index_pair_type &triangleIndices,
                    index_pair_type &vertexPair);
    void getFeature1(const index_pair_type &triangleIndices,
                      const index_pair_type &vertexPair);
    void getFeature2(const index_pair_type &triangleIndices,
                      const index_pair_type &vertexPair);

    void setOrphans();
    void testOrphans();
    bool testOrphansEdgeToEdge(const index_type &i1, const index_type &i2);
    bool testOrphansVertexToFace(const index_type &i1, const index_type &i2);
    void insertEdgeToEdgeFeature(const index_type &i1, const index_type &i2);
    void insertVertexToFaceFeature(const index_type &i1, const index_type &i2);
    float intersectVertexFace(const double &dt, const index_type &v, const index_type &f);
    float intersectVertexFace(const double &dt, const index_type &v, const index_type &f1, const index_type &f2);

    float intersectEdgeEdge(const double &dt, const index_type &e1, const index_type &e2);
    float intersectEdgeEdge(const double &dt, const index_type &e1, const index_type &e2, const index_type &f1, const index_type &f2);
    bool testVertexFace(const double &dt, const index_type &e, const index_type &e);
    bool testEdgeEdge(const double &dt, const index_type &e, const index_type &e);
    bool testCoplanarity(const Vector3 &p0, const Vector3 &p1, const Vector3 &p2,
                         const Vector3 &p3, const Vector3 &p4, const Vector3 &p5,
                         const Vector3 &p6, const Vector3 &p7);
    float vertexFaceIntersection(const double &dt, const index_type &v, const index_type &f);
    float edgeEdgeIntersection(const double &dt, const index_type &e1, const index_type &e2);
    void computeCubicCoefficients(const Vector3 &a0, const Vector3 &ad, const Vector3 &b0, const Vector3 &bd,
                                  const Vector3 &c0, const Vector3 &cd, const Vector3 &p0, const Vector3 &pd,
                                  Real &a, Real &b, Real &c, Real &d);
    void computeCubicCoefficientsVertexFace(const Vector3 &a0, const Vector3 &ad, const Vector3 &b0, const Vector3 &bd,
                                       const Vector3 &c0, const Vector3 &cd, const Vector3 &p0, const Vector3 &pd,
                                       Real &a, Real &b, Real &c, Real &d);
    void computeCubicCoefficientsEdgeEdge(const Vector3 &a0, const Vector3 &ad, const Vector3 &b0, const Vector3 &bd,
                                     const Vector3 &c0, const Vector3 &cd, const Vector3 &d0, const Vector3 &dd,
                                     Real &a, Real &b, Real &c, Real &d);
    void draw(const core::visual::VisualParams*,int index);

    void draw(const core::visual::VisualParams* vparams);

    virtual bool canCollideWithElement(index_type index, CollisionModel* model2, index_type index2);

    virtual void handleTopologyChange();

    core::behavior::MechanicalState<DataTypes>* getMechanicalState() {
        return mstate;
    }
    const core::behavior::MechanicalState<DataTypes>* getMechanicalState() const {
        return mstate;
    }

    const VecCoord& getX() const {
        return *(getMechanicalState()->getX());
    }
    const sofa::core::topology::BaseMeshTopology::SeqTriangles& getTriangles() const {
        return *triangles;
    }
    const VecDeriv& getNormals() const {
        return normals;
    }

    TriangleLocalMinDistanceFilter *getFilter() const;

    void setFilter(TriangleLocalMinDistanceFilter * /*lmdFilter*/);

    Deriv velocity(int index) const;

    /// Pre-construction check method called by ObjectFactory.
    /// Check that DataTypes matches the MechanicalState.
    template<class T>
    static bool canCreate(T*& obj, core::objectmodel::BaseContext* context, core::objectmodel::BaseObjectDescription* arg)
    {
        if (dynamic_cast<core::behavior::MechanicalState<DataTypes>*>(context->getMechanicalState()) == NULL)
            return false;
        return BaseObject::canCreate(obj, context, arg);
    }

    virtual std::string getTemplateName() const
    {
        return templateName(this);
    }

    static std::string templateName(const TRTriangleModel<DataTypes>* = NULL)
    {
        return DataTypes::Name();
    }
};

template<class DataTypes>
inline TRTriangle<DataTypes>::TRTriangle(ParentModel* model,
                                         typename TRTriangle<DataTypes>::index_type index)
    : core::TCollisionElementIterator<ParentModel>(model, index)
{}

template<class DataTypes>
inline TRTriangle<DataTypes>::TRTriangle(const core::CollisionElementIterator& i)
    : core::TCollisionElementIterator<ParentModel>(static_cast<ParentModel*>(i.getCollisionModel()), i.getIndex())
{}

template<class DataTypes>
inline TRTriangle<DataTypes>::TRTriangle(ParentModel* model,
                                         typename TRTriangle<DataTypes>::index_type index)
    : core::TCollisionElementIterator<ParentModel>(model, index)
{}

template<class DataTypes>
inline const typename DataTypes::Coord& TRTriangle<DataTypes>::p1() const {
    return (*this->model->mstate->getX())[(*(this->model->triangles))[this->index][0]];
}
template<class DataTypes>
inline const typename DataTypes::Coord& TRTriangle<DataTypes>::p2() const {
    return (*this->model->mstate->getX())[(*(this->model->triangles))[this->index][1]];
}
template<class DataTypes>
inline const typename DataTypes::Coord& TRTriangle<DataTypes>::p3() const {
    return (*this->model->mstate->getX())[(*(this->model->triangles))[this->index][2]];
}
template<class DataTypes>
inline const typename DataTypes::Coord& TRTriangle<DataTypes>::p(size_t i) const {
    return (*this->model->mstate->getX())[(*(this->model->triangles))[this->index][i]];
}

template<class DataTypes>
inline const bool TRTriangle<DataTypes>::covertices(const TRTriangle<DataTypes> &other)
{
    if(p1Index() == other.p1Index())
        return true;

    if(p1Index() == other.p2Index())
        return true;

    if(p1Index() == other.p3Index())
        return true;

    if(p2Index() == other.p1Index())
        return true;

    if(p2Index() == other.p2Index())
        return true;

    if(p2Index() == other.p3Index())
        return true;

    if(p3Index() == other.p1Index())
        return true;

    if(p3Index() == other.p2Index())
        return true;

    if(p3Index() == other.p3Index())
        return true;

    return false;
}

template<class DataTypes>
inline const bool TRTriangle<DataTypes>::coedge(const TRTriangle<DataTypes> &other)
{
    int numberOfCommonVertices = 0;
    if(p1Index() == other.p1Index())
        numberOfCommonVertices++;

    if(p1Index() == other.p2Index())
        numberOfCommonVertices++;

    if(p1Index() == other.p3Index())
        numberOfCommonVertices++;

    if(p2Index() == other.p1Index())
        numberOfCommonVertices++;

    if(p2Index() == other.p2Index())
        numberOfCommonVertices++;

    if(p2Index() == other.p3Index())
        numberOfCommonVertices++;

    if(p3Index() == other.p1Index())
        numberOfCommonVertices++;

    if(p3Index() == other.p2Index())
        numberOfCommonVertices++;

    if(p3Index() == other.p3Index())
        numberOfCommonVertices++;

    return numberOfCommonVertices > 1;
}

template<class DataTypes>
inline const typename TRTriangle<DataTypes>::index_type TRTriangle<DataTypes>::covertices(const TRTriangle<DataTypes> &other,
                                                                                          typename TRTriangle<DataTypes>::index_pair_type &vertexPair)
{
    helper::vector<index_pair_type> keeps;
    if(p1Index() == other.p1Index())
        keeps.push_back(std::make_pair(0,0));

    if(p1Index() == other.p2Index())
        keeps.push_back(std::make_pair(0,1));

    if(p1Index() == other.p3Index())
        keeps.push_back(std::make_pair(0,2));

    if(p2Index() == other.p1Index())
        keeps.push_back(std::make_pair(1,0));

    if(p2Index() == other.p2Index())
        keeps.push_back(std::make_pair(1,1));

    if(p2Index() == other.p3Index())
        keeps.push_back(std::make_pair(1,2));

    if(p3Index() == other.p1Index())
        keeps.push_back(std::make_pair(2,0));

    if(p3Index() == other.p2Index())
        keeps.push_back(std::make_pair(2,1));

    if(p3Index() == other.p3Index())
        keeps.push_back(std::make_pair(2,2));

    size_t numberOfCommonVertices = keeps.size();
    assert(numberOfCommonVertices <= 3);

    vertexPair.first = vertexPair.second = 0;

    if(numberOfCommonVertices == 1)
    {
        vertexPair.first = keeps[0].first;
        vertexPair.second = keeps[0].second;
    }
    else if(numberOfCommonVertices == 2)
    {
        if(keeps[0].first == 0 && keeps[1].first == 0)
            vertexPair.first = 0;
        if(keeps[0].first == 1 && keeps[1].first == 1)
            vertexPair.first = 1;
        if(keeps[0].first == 2 && keeps[1].first == 2)
            vertexPair.first = 2;
        if(keeps[0].second == 0 && keeps[1].second == 0)
            vertexPair.second = 0;
        if(keeps[0].second == 1 && keeps[1].second == 1)
            vertexPair.second = 1;
        if(keeps[0].second == 2 && keeps[1].second == 2)
            vertexPair.second = 2;
    }
    else if (numberOfCommonVertices == 3)
        return numberOfCommonVertices-1;

}

template<class DataTypes>
inline const typename DataTypes::Coord& TRTriangle<DataTypes>::operator[](size_t i) const {
    return (*this->model->mstate->getX())[(*(this->model->triangles))[this->index][i]];
}

template<class DataTypes>
inline const typename DataTypes::Coord& TRTriangle<DataTypes>::p1Free() const {
    return (*this->model->mstate->getXfree())[(*(this->model->triangles))[this->index][0]];
}
template<class DataTypes>
inline const typename DataTypes::Coord& TRTriangle<DataTypes>::p2Free() const {
    return (*this->model->mstate->getXfree())[(*(this->model->triangles))[this->index][1]];
}
template<class DataTypes>
inline const typename DataTypes::Coord& TRTriangle<DataTypes>::p3Free() const {
    return (*this->model->mstate->getXfree())[(*(this->model->triangles))[this->index][2]];
}

template<class DataTypes>
inline const size_t &TRTriangle<DataTypes>::p1Index() const {
    return (*(this->model->triangles))[this->index][0];
}
template<class DataTypes>
inline const size_t &TRTriangle<DataTypes>::p2Index() const {
    return (*(this->model->triangles))[this->index][1];
}
template<class DataTypes>
inline const size_t &TRTriangle<DataTypes>::p3Index() const {
    return (*(this->model->triangles))[this->index][2];
}
template<class DataTypes>
inline const size_t &TRTriangle<DataTypes>::Index(int i) const {
  return (*(this->model->triangles))[this->index][i];
}

template<class DataTypes>
inline const size_t &TRTriangle<DataTypes>::e1Index() const {
    return this->model->triangleEdgeFeatures[this->index][0];
}
template<class DataTypes>
inline const size_t &TRTriangle<DataTypes>::e2Index() const {
    return this->model->triangleEdgeFeatures[this->index][1];
}
template<class DataTypes>
inline const size_t &TRTriangle<DataTypes>::e3Index() const {
    return this->model->triangleEdgeFeatures[this->index][2];
}

template<class DataTypes>
inline const typename DataTypes::Deriv& TRTriangle<DataTypes>::v1() const {
    return (*this->model->mstate->getV())[(*(this->model->triangles))[this->index][0]];
}
template<class DataTypes>
inline const typename DataTypes::Deriv& TRTriangle<DataTypes>::v2() const {
    return (*this->model->mstate->getV())[(*(this->model->triangles))[this->index][1]];
}
template<class DataTypes>
inline const typename DataTypes::Deriv& TRTriangle<DataTypes>::v3() const {
    return (*this->model->mstate->getV())[(*(this->model->triangles))[this->index][2]];
}
template<class DataTypes>
inline const typename DataTypes::Deriv& TRTriangle<DataTypes>::v(size_t i) const {
    return (*this->model->mstate->getV())[(*(this->model->triangles))[this->index][i]];
}

template<class DataTypes>
inline const typename DataTypes::Deriv& TRTriangle<DataTypes>::n() const {
    return this->model->normals[this->index];
}
template<class DataTypes>
inline       typename DataTypes::Deriv& TRTriangle<DataTypes>::n()       {
    return this->model->normals[this->index];
}

template<class DataTypes>
inline size_t TRTriangle<DataTypes>::flags() const {
    return this->model->getTriangleFlags(this->index);
}

template<class DataTypes>
inline bool TRTriangle<DataTypes>::hasFreePosition() const {
    return this->model->mstate->read(core::ConstVecCoordId::freePosition())->isSet();
}

template<class DataTypes>
inline typename DataTypes::Deriv TRTriangleModel<DataTypes>::velocity(size_t index) const {
    return ((*mstate->getV())[(*(triangles))[index][0]] + (*mstate->getV())[(*(triangles))[index][1]] +
            (*mstate->getV())[(*(triangles))[index][2]])/((Real)(3.0));
}
typedef TRTriangleModel<Vec3Types> RTriangleModel;
typedef TRTriangle<Vec3Types> RTriangle;

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_BUILD_MESH_COLLISION)
#ifndef SOFA_FLOAT
extern template class SOFA_MESH_COLLISION_API TRTriangleModel<defaulttype::Vec3dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_MESH_COLLISION_API TRTriangleModel<defaulttype::Vec3fTypes>;
#endif
#endif

} // namespace collision

} // namespace component

} // namespace sofa

#endif

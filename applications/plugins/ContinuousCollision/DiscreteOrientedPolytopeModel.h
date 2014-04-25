#ifndef SOFA_COMPONENT_COLLISION_DISCRETEORIENTEDPOLYTOPEMODEL_H
#define SOFA_COMPONENT_COLLISION_DISCRETEORIENTEDPOLYTOPEMODEL_H

#include <sofa/core/CollisionModel.h>
#include <sofa/component/container/MechanicalObject.h>
#include <sofa/defaulttype/Vec3Types.h>
#include<limits>
#include<array>

namespace sofa
{

namespace component
{

namespace collision
{


using namespace sofa::defaulttype;

template<int k>
class DiscreteOrientedPolytopeModel;

template<int k>
class Polytope : public core::TCollisionElementIterator<DiscreteOrientedPolytopeModel<k> >
{

public:
  Polytope(DiscreteOrientedPolytopeModel<k>* model=NULL, int index=0);

  explicit Polytope(const core::CollisionElementIterator& i);

  const Vector3& minVect() const;

  const Vector3& maxVect() const;

  const std::pair<DiscreteOrientedPolytopeModel,DiscreteOrientedPolytopeModel>& subcells() const;
};

template<int _k>
class SOFA_BASE_COLLISION_API DiscreteOrientedPolytopeModel : public core::CollisionModel
{
public:
  SOFA_CLASS(SOFA_TEMPLATE(DiscreteOrientedPolytopeModel,_k),sofa::core::CollisionModel);

  enum { K = _k };

  void inline setDistances6(Vector3 &p)
  {
    dist[0] = p[0];
    dist[1] = p[1];
    dist[2] = p[2];
  }

  template<>
  void inline setDistances18(Vector3 &p)
  {
    dist[3] = p[0] + p[1];
    dist[4] = p[0] + p[2];
    dist[5] = p[1] + p[2];
    dist[6] = p[0] - p[1];
    dist[7] = p[0] - p[2];
    dist[7] = p[1] - p[2];
  }

  template<>
  void inline setDistances<6>(Vector3 &p)
  {
    this->setDistances<5>(p);
    dist[8] = p[1] - p[2];
  }

  template<>
  void inline setDistances<9>(Vector3 &p)
  {
    this->setDistances<6>(p);
    dist[9] = dist[3] - p[2];
    dist[10] = dist[6] + p[2];
    dist[11] = dist[5] - p[0];
  }

  template<>
  void inline setDistances<16>(Vector3 &p)
  {
    this->setDistances<0>(p);
    this->setDistances<5>(p);
  }

  template<>
  void inline setDistances<18>(Vector3 &p)
  {
    this->setDistances<0>(p);
    this->setDistances<6>(p);
  }

  template<>
  void inline setDistances<24>(Vector3 &p)
  {
    this->setDistances<0>(p);
    this->setDistances<9>(p);
  }

  DiscreteOrientedPolytope()
  {
    T real_max = std::numeric_limits<T>::max();
    for(size_t i = 0; i < K / 2; ++i)
    {
      const size_t half = K/2;
      dist[i] = real_max;
      dist[i + half] = -real_max;
    }
  }
  DiscreteOrientedPolytope(const Vector3& v)
  {
    this->setDistances<(K - 6) / 2>(v);
  }

  DiscreteOrientedPolytope(const Vector3& v,const Vector3& v)
  {

  }

private:
  std::array<float,K> dist;
};

} // namespace collision

} // namespace component

} // namespace sofa

#endif

#ifndef SOFA_COMPONENT_COLLISION_DISCRETEORIENTEDPOLYTOPEMODEL_H
#define SOFA_COMPONENT_COLLISION_DISCRETEORIENTEDPOLYTOPEMODEL_H

// STL includes
#include<limits>

// SOFA includes
#include <sofa/core/CollisionModel.h>
#include <sofa/component/container/MechanicalObject.h>
#include <sofa/defaulttype/Vec3Types.h>

// Local includes
#include"NormalProjections.hpp"

namespace sofa
{

namespace component
{

namespace collision
{
template<typename T, int _K>
struct DiscreteOrientedPolytope
{
public:
    enum { K = _K, KHalf = _K/2 };

    DiscreteOrientedPolytope()
    {
        this->clean();
    }

    void DiscreteOrientedPolytope(const defaulttype::Vector3 &p)
    {
        NormalProjections<K>::ComputeAll(p,this->Distance);
    }

    void DiscreteOrientedPolytope(const defaulttype::Vector3 &p, const defaulttype::Vector3 &q)
    {
        NormalProjections<K>::ComputeAll(p,q,this->Distance);
    }

    bool overlaps(const DiscreteOrientedPolytope<T,K> &other)
    {
        for(size_t i = 0; i < KHalf; ++i)
        {
            if(this->Distance[i] > other.Distance[i] ) return false;
            if(this->Distance[i+KHalf] < other.Distance[i+KHalf]) return false;
        }
        return true;
    }

    bool overlaps(const DiscreteOrientedPolytope<T,K> &a, DiscreteOrientedPolytope<T,K> &b)
    {
        if(!this->overlaps(a))
            return false;
        for(size_t i = 0; i < KHalf; ++i)
        {
            b.Distance[i] = std::max(this->Distance[i],a.Distance[i]);
            b.Distance[i+KHalf] = std::min(this->Distance[i+KHalf],a.Distance[i+KHalf]);
        }
        return true;
    }

    bool inside(const defaulttype::Vector3 &p)
    {
        std::array<T,KHalf> d = {0};
        NormalProjections<K>::ComputeAll(p,d);
        for(size_t i = 0; i < KHalf; ++i)
            if(this->Distance[i] > d[i] || this->Distance[i+KHalf] < d[i])
                return false;
        return true;
    }

    DiscreteOrientedPolytope<T,K> &operator+=(const defaulttype::Vector3 &p)
    {
        std::array<T,KHalf> d = {0};
        NormalProjections<K>::ComputeAll(p,d);
        for(size_t i = 0; i < KHalf; ++i)
        {
            this->Distance[i] = std::min(d[i],this->Distance[i]);
            this->Distance[i+KHalf] = std::min(d[i],this->Distance[i+KHalf]);
        }
        return *this;
    }

    DiscreteOrientedPolytope<T,K> &operator+=(const DiscreteOrientedPolytope<T,K> &other)
    {
        for(size_t i = 0; i < KHalf; ++i)
        {
            this->Distance[i] = std::min(this->Distance[i],other.Distance[i]);
            this->Distance[i+KHalf] = std::max(this->Distance[i+KHalf],other.Distance[i+KHalf]);
        }
        return *this;
    }

    DiscreteOrientedPolytope<T,K> operator+(const DiscreteOrientedPolytope<T,K> &other)
    {
        DiscreteOrientedPolytope<T,K> result(*this);
        return (result += other);
    }

    T length(size_t i)
    {
        return this->Distance[i+KHalf] - this->Distance[i];
    }

    T width()  const
    {
        return this->Distance[KHalf] - this->Distance[0];
    }

    T height() const
    {
        return this->Distance[KHalf+1] - this->Distance[1];
    }

    T depth()  const
    {
        return this->Distance[KHalf+2] - this->Distance[2];
    }

    T volume() const
    {
        return width()*height()*depth();
    }

    defaulttype::Vector3 center() const {
        return defaulttype::Vector3(this->Distance[0]+this->Distance[KHalf], this->Distance[1]+this->Distance[KHalf+1], this->Distance[2]+this->Distance[KHalf+2])*T(0.5);
    }

    void clean()
    {
        for(size_t i = 0; i < KHalf; ++i )
        {
            this->Distance[i] = std::numeric_limits<T>::infinity();
            this->Distance[i+KHalf] = -this->Distance[i];
        }
    }

private:
    std::array<T,K> Distance;
};


} // namespace collision

} // namespace component

} // namespace sofa

#endif

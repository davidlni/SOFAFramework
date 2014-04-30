#ifndef SOFA_COMPONENT_COLLISION_DISCRETEORIENTEDPOLYTOPEMODEL_H
#define SOFA_COMPONENT_COLLISION_DISCRETEORIENTEDPOLYTOPEMODEL_H
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
// STL includes
#include<limits>

// SOFA includes
#include <sofa/core/CollisionModel.h>
#include <sofa/component/container/MechanicalObject.h>
#include <sofa/defaulttype/Vec3Types.h>

// Local includes
#include"NormalProjections.h"

namespace sofa
{

namespace component
{

namespace collision
{

template<typename T, size_t NumberOfPlanes = 6>
class DiscreteOrientedPolytope
{
public:
    enum { K = NumberOfPlanes, KHalf = NumberOfPlanes/2 };
    typedef std::array<T,K> DistanceArrayType;
    typedef std::array<T,KHalf> HalfDistanceArrayType;


public:

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

    friend bool overlaps(const DiscreteOrientedPolytope<T,K> &other)
    {
        for(size_t i = 0; i < KHalf; ++i)
        {
            if( this->Distance[i] > other.Distance[i] || this->Distance[i+KHalf] < other.Distance[i+KHalf] )
                return false;
        }
        return true;
    }

    friend bool overlaps(const DiscreteOrientedPolytope<T,K> &a, DiscreteOrientedPolytope<T,K> &out)
    {
        if(!this->overlaps(a))
            return false;

        for(size_t i = 0; i < KHalf; ++i)
        {
            out.Distance[i] = std::max(this->Distance[i],a.Distance[i]);
            out.Distance[i+KHalf] = std::min(this->Distance[i+KHalf],a.Distance[i+KHalf]);
        }
        return true;
    }

    bool inside(const defaulttype::Vector3 &p)
    {
        HalfDistanceArrayType d = {0};
        NormalProjections<K>::ComputeAll(p,d);
        for(size_t i = 0; i < KHalf; ++i)
            if(this->Distance[i] > d[i] || this->Distance[i+KHalf] < d[i])
                return false;
        return true;
    }

    DiscreteOrientedPolytope<T,K> &operator+=(const defaulttype::Vector3 &p)
    {
        HalfDistanceArrayType d = {0};
        NormalProjections<K>::ComputeAll(p,d);
        for(size_t i = 0; i < KHalf; ++i)
        {
            this->Distance[i] = std::min(d[i],this->Distance[i]);
            this->Distance[i+KHalf] = std::min(d[i],this->Distance[i+KHalf]);
        }
        return *this;
    }

    friend DiscreteOrientedPolytope<T,K> &operator+=(const DiscreteOrientedPolytope<T,K> &other)
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

    void operator=(const DiscreteOrientedPolytope<T,K> &other)
    {
        this->Distance = other.Distance;
    }

    T length(size_t i) const
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

    int GetSplitAxis()
    {
        T w = width();
        T h = height();
        T d = depth();
        if(w >= h && w >= d)
            return 0;
        else if(h >= w && h >= d)
            return 1;
        return 2;
    }

    defaulttype::Vector3 GetCenter() const {
        return defaulttype::Vector3(this->Distance[0]+this->Distance[KHalf],
                                    this->Distance[1]+this->Distance[KHalf+1],
                                    this->Distance[2]+this->Distance[KHalf+2])*T(0.5);
    }

    T GetCenter(size_t i) const {
      return (this->Distance[i]+this->Distance[i+KHalf])*T(0.5);
    }

    defaulttype::Vector3 GetLength() const {
        return defaulttype::Vector3(this->Distance[KHalf]-this->Distance[0],
                                    this->Distance[KHalf+1]-this->Distance[1],
                                    this->Distance[KHalf+2]-this->Distance[2]);
    }

    void clean()
    {
        for(size_t i = 0; i < KHalf; ++i )
        {
            this->Distance[i] = std::numeric_limits<T>::infinity();
            this->Distance[i+KHalf] = -this->Distance[i];
        }
    }

    inline const T& GetDistance(size_t i) const
    {
        return this->Distance[i];
    }

    inline T& GetDistance(size_t i)
    {
        return this->Distance[i];
    }

    inline DistanceArrayType& GetDistance()
    {
        return this->Distance;
    }

protected:
    DistanceArrayType Distance;
};


} // namespace collision

} // namespace component

} // namespace sofa

#endif

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

template<typename T, size_t NumberOfPlanes = 18>
class DiscreteOrientedPolytope
{
public:
    enum { K = NumberOfPlanes, KHalf = NumberOfPlanes/2 };
    typedef std::array<T,K> DistanceArrayType;
    typedef std::array<T,KHalf> HalfDistanceArrayType;

public:

    /**
     * @brief Constructor
     *  Builds a default dop covering the entire 3D space
     *
     */
    DiscreteOrientedPolytope()
    {
      for(size_t i = 0; i < KHalf; ++i )
      {
        this->Distance[i] = std::numeric_limits<T>::max();
        this->Distance[i+KHalf] = -this->Distance[i];
      }
    }

    DiscreteOrientedPolytope(const defaulttype::Vector3 &p)
    {
        NormalProjections<K>::ComputeAll(p,this->Distance);
    }

    DiscreteOrientedPolytope(const defaulttype::Vector3 &p, const defaulttype::Vector3 &q)
    {
        NormalProjections<K>::ComputeAll(p,q,this->Distance);
    }

    bool Overlaps(const DiscreteOrientedPolytope<T,K> &other) const
    {
        for(size_t i = 0; i < KHalf; ++i)
        {
            if( this->Distance[i] > other.Distance[i] || this->Distance[i+KHalf] < other.Distance[i+KHalf] )
                return false;
        }
        return true;
    }

    bool Overlaps(const DiscreteOrientedPolytope<T,K> &a, DiscreteOrientedPolytope<T,K> &out) const
    {
        if(!this->Overlaps(a))
            return false;

        for(size_t i = 0; i < KHalf; ++i)
        {
            out.Distance[i] = std::max(this->Distance[i],a.Distance[i]);
            out.Distance[i+KHalf] = std::min(this->Distance[i+KHalf],a.Distance[i+KHalf]);
        }
        return true;
    }

    bool Inside(const defaulttype::Vector3 &p) const
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
            this->Distance[i+KHalf] = std::max(d[i],this->Distance[i+KHalf]);
        }
        return *this;
    }

    DiscreteOrientedPolytope<T,K> &operator+=(const T &distance)
    {
      for(size_t i = 0; i < KHalf; ++i)
      {
        this->Distance[i] -= distance;
        this->Distance[i+KHalf] += distance;
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

    inline DiscreteOrientedPolytope<T,K> operator+(const DiscreteOrientedPolytope<T,K> &other) const
    {
        DiscreteOrientedPolytope<T,K> result(*this);
        return (result += other);
    }

    inline DiscreteOrientedPolytope<T,K> operator+(const T &other) const
    {
      DiscreteOrientedPolytope<T,K> result(*this);
      return (result += other);
    }

    inline DiscreteOrientedPolytope<T,K> operator+(const defaulttype::Vector3 &p) const
    {
      DiscreteOrientedPolytope<T,K> result(*this);
      return (result += p);
    }

    inline DiscreteOrientedPolytope<T,K> &operator*=(const defaulttype::Vector3 &p)
    {
        DiscreteOrientedPolytope<T,K> dop(p);
        for(size_t i = 0; i < K; ++i)
            this->Distance[i] += dop.GetDistance(i);
        return *this;
    }

    inline DiscreteOrientedPolytope<T,K> operator*(const defaulttype::Vector3 &p)
    {
        DiscreteOrientedPolytope<T,K> dop(p);
        return dop*=(*this);
    }

    inline bool operator%(const DiscreteOrientedPolytope<T,K> &other)
    {
        return this->Overlaps(other);
    }

    inline void operator=(const DiscreteOrientedPolytope<T,K> &other)
    {
        this->Distance = other.Distance;
    }

    inline T GetLength(size_t i) const
    {
        return this->Distance[i+KHalf] - this->Distance[i];
    }

    inline T GetWidth() const
    {
        return this->Distance[KHalf] - this->Distance[0];
    }

    inline T GetHeight() const
    {
        return this->Distance[KHalf+1] - this->Distance[1];
    }

    inline T GetDepth() const
    {
        return this->Distance[KHalf+2] - this->Distance[2];
    }

    inline T GetVolume() const
    {
        return this->GetWidth()*this->GetHeight()*this->GetDepth();
    }

    inline size_t GetSplitAxis() const
    {
        T w = this->GetWidth();
        T h = this->GetHeight();
        T d = this->GetDepth();
        if(w >= h && w >= d)
            return 0;
        else if(h >= w && h >= d)
            return 1;
        return 2;
    }

    inline defaulttype::Vector3 GetCenter() const {
        return defaulttype::Vector3(this->Distance[0]+this->Distance[KHalf],
                                    this->Distance[1]+this->Distance[KHalf+1],
                                    this->Distance[2]+this->Distance[KHalf+2])*T(0.5);
    }

    inline T GetCenter(const size_t &i) const {
        return (this->Distance[i]+this->Distance[i+KHalf])*T(0.5);
    }

    inline defaulttype::Vector3 GetLength() const {
        return defaulttype::Vector3(this->Distance[KHalf]-this->Distance[0],
                                    this->Distance[KHalf+1]-this->Distance[1],
                                    this->Distance[KHalf+2]-this->Distance[2]);
    }

    /**
     * @brief ...
     *
     * @return void
     */
    void Clean()
    {
        for(size_t i = 0; i < KHalf; ++i )
        {
            this->Distance[i] = std::numeric_limits<T>::max();
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

    inline DistanceArrayType& Enlarge(const T &distance)
    {
      return this->Distance;
    }
    
    inline defaulttype::Vector3 GetBoundingBoxMin()
    {
      return defaulttype::Vector3(this->Distance[0],this->Distance[1],this->Distance[2]);
    }
    
    inline defaulttype::Vector3 GetBoundingBoxMax()
    {
      return defaulttype::Vector3(this->Distance[KHalf],this->Distance[KHalf+1],this->Distance[KHalf+2]);
    }

protected:
    DistanceArrayType Distance;
};


} // namespace collision

} // namespace component

} // namespace sofa

#endif

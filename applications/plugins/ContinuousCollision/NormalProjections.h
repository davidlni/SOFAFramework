#ifndef NORMAL_PROJECTIONS_H
#define NORMAL_PROJECTIONS_H
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

#include<sofa/defaulttype/Vec3Types.h>
#include<array>
#include<algorithm>
#include<type_traits>
///
/// These specialized clsses compute projections of points in 3D
/// space onto normals of the k-DOP bounding planes. The computation
/// is done in a very specific order that the denpends on the value of k.
/// The value of k must be either 6,12,18 or 24
///
/// For K = 6, this is the commonly used AABB:
/// (-1,0,0) and (1,0,0)  -> indices 0 and 4
/// (0,-1,0) and (0,1,0)  -> indices 1 and 5
/// (0,0,-1) and (0,0,1)  -> indices 2 and 6
/// For K = 12, the planes are 6 AABB planes and 6 diagonal planes that cut off some space of the edges:
/// (-1,0,0) and (1,0,0)  -> indices 0 and 6
/// (0,-1,0) and (0,1,0)  -> indices 1 and 7
/// (0,0,-1) and (0,0,1)  -> indices 2 and 8
/// (-1,-1,0) and (1,1,0) -> indices 3 and 9
/// (-1,0,-1) and (1,0,1) -> indices 4 and 10
/// (0,-1,-1) and (0,1,1) -> indices 5 and 11
/// For K = 18, the planes are 12-DOP planes and 6 aditional diagonal planes that cut off some space of the edges:
/// (-1,0,0) and (1,0,0)  -> indices 0 and 9
/// (0,-1,0) and (0,1,0)  -> indices 1 and 10
/// (0,0,-1) and (0,0,1)  -> indices 2 and 11
/// (-1,-1,0) and (1,1,0) -> indices 3 and 12
/// (-1,0,-1) and (1,0,1) -> indices 4 and 13
/// (0,-1,-1) and (0,1,1) -> indices 5 and 14
/// (-1,1,0) and (1,-1,0) -> indices 6 and 15
/// (-1,0,1) and (1,0,-1) -> indices 7 and 16
/// (0,-1,1) and (0,1,-1) -> indices 8 and 17
/// For K = 24, the planes are 18 planes and 6 additional diagonal planes that cut off some space of the edges:
/// (-1,0,0) and (1,0,0)  -> indices 0 and 12
/// (0,-1,0) and (0,1,0)  -> indices 1 and 13
/// (0,0,-1) and (0,0,1)  -> indices 2 and 14
/// (-1,-1,0) and (1,1,0) -> indices 3 and 15
/// (-1,0,-1) and (1,0,1) -> indices 4 and 16
/// (0,-1,-1) and (0,1,1) -> indices 5 and 17
/// (-1,1,0) and (1,-1,0) -> indices 6 and 18
/// (-1,0,1) and (1,0,-1) -> indices 7 and 19
/// (0,-1,1) and (0,1,-1) -> indices 8 and 20
/// (-1, -1, 1) and (1, 1, -1) --> indices 9 and 21
/// (-1, 1, -1) and (1, -1, 1) --> indices 10 and 22
/// (1, -1, -1) and (-1, 1, 1) --> indices 11 and 23
///
namespace sofa
{

namespace component
{

namespace collision
{

template<int K>
struct NormalProjections;

template<>
struct NormalProjections<6>
{
    template<int K, typename T>
    inline static void Compute(const defaulttype::Vector3 &p, std::array<T,K/2> &d)
    {
        static_assert(K >= 6,"Projection requested is can't be computed.");
        for(size_t i = 0; i < 3; ++i)
            d[i] = p[i];
    }

    template<int K, typename T>
    inline static void Compute(const defaulttype::Vector3 &p, std::array<T,K> &d)
    {
        static_assert(K >= 6,"Projection requested is can't be computed.");
        for(size_t i = 0; i < 3; ++i)
            d[i] = d[i+K/2] = p[i];
    }

    template<int K, typename T>
    inline static void Compute(const defaulttype::Vector3 &p, const defaulttype::Vector3 &q, std::array<T,K> &d)
    {
        static_assert(K >= 6,"Projection requested is can't be computed.");
        std::pair<T,T> minmax;
        for(size_t i = 0; i < 3; ++i)
        {
            minmax = std::minmax(p[i],q[i]);
            d[i] = minmax.first;
            d[i+K/2] = minmax.second;
        }
    }

    template<typename T>
    inline static void ComputeAll(const defaulttype::Vector3 &p, std::array<T,K> &d)
    {
        Compute<6>(p,d);
    }

    template<typename T>
    inline static void ComputeAll(const defaulttype::Vector3 &p, std::array<T,K/2> &d)
    {
        Compute<6>(p,d);
    }

};

template<>
struct NormalProjections<12>
{
    template<int K, typename T>
    inline static void Compute(const defaulttype::Vector3 &p, std::array<T,K/2> &d)
    {
        static_assert(K >= 12,"Projection requested is can't be computed.");
        for(size_t i = 0, idx = 3; i < 2; ++i)
            for(size_t j = i+1; j < 3; ++j, ++idx)
                d[idx] = p[i]+p[j];
    }

    template<int K, typename T>
    inline static void Compute(const defaulttype::Vector3 &p, std::array<T,K> &d)
    {
        static_assert(K >= 12,"Projection requested is can't be computed.");
        for(size_t i = 0, idx = 3; i < 2; ++i)
            for(size_t j = i+1; j < 3; ++j, ++idx)
                d[idx] = d[idx+K/2] = p[i]+p[j];
    }

    template<int K, typename T>
    inline static void Compute(const defaulttype::Vector3 &p, const defaulttype::Vector3 &q, std::array<T,K> &d)
    {
        static_assert(K >= 12,"Projection requested is can't be computed.");
        std::pair<T,T> minmax;
        for(size_t i = 0, idx = 3; i < 2; ++i)
            for(size_t j = i+1; j < 3; ++j, ++idx)
            {
                minmax = std::minmax(p[i]+p[j],q[i]+q[j]);
                d[idx] = minmax.first;
                d[idx+K/2] = minmax.second;
            }
    }

    template<typename T>
    inline static void ComputeAll(const defaulttype::Vector3 &p, std::array<T,6> &d)
    {
        NormalProjections<6>::Compute<12>(p,d);
        Compute<12>(p,d);
    }

    template<typename T>
    inline static void ComputeAll(const defaulttype::Vector3 &p, std::array<T,12> &d)
    {
        NormalProjections<6>::Compute<12>(p,d);
        Compute<12>(p,d);
    }

    template<typename T>
    inline static void ComputeAll(const defaulttype::Vector3 &p, const defaulttype::Vector3 &q, std::array<T,12> &d)
    {
        NormalProjections<6>::Compute<12>(p,q,d);
        Compute<12>(p,q,d);
    }

};

template<>
struct NormalProjections<18>
{
    template<int K, typename T>
    inline static void Compute(const defaulttype::Vector3 &p, std::array<T,K/2> &d)
    {
        static_assert(K >= 18,"Projection requested is can't be computed.");
        for(size_t i = 0, idx = 6; i < 2; ++i)
            for(size_t j = i+1; j < 3; ++j, ++idx)
                d[idx] = p[i]-p[j];
    }

    template<int K, typename T>
    inline static void Compute(const defaulttype::Vector3 &p, std::array<T,K> &d)
    {
        static_assert(K >= 18,"Projection requested is can't be computed.");
        for(size_t i = 0, idx = 6; i < 2; ++i)
            for(size_t j = i+1; j < 3; ++j, ++idx)
                d[idx] = d[idx+K/2] = p[i]-p[j];
    }

    template<int K, typename T>
    inline static void Compute(const defaulttype::Vector3 &p, const defaulttype::Vector3 &q, std::array<T,K> &d)
    {
        static_assert(K >= 18,"Projection requested is can't be computed.");
        std::pair<T,T> minmax;
        for(size_t i = 0, idx = 6; i < 2; ++i)
            for(size_t j = i+1; j < 3; ++j, ++idx)
            {
                minmax = std::minmax(p[i]-p[j],q[i]-q[j]);
                d[idx] = minmax.first;
                d[idx+K/2] = minmax.second;
            }
    }

    template<typename T>
    inline static void ComputeAll(const defaulttype::Vector3 &p, std::array<T,9> &d)
    {
        NormalProjections<6>::Compute<18>(p,d);
        NormalProjections<12>::Compute<18>(p,d);
        Compute<18>(p,d);
    }

    template<typename T>
    inline static void ComputeAll(const defaulttype::Vector3 &p, std::array<T,18> &d)
    {
        NormalProjections<6>::Compute<18>(p,d);
        NormalProjections<12>::Compute<18>(p,d);
        Compute<18>(p,d);
    }

    template<typename T>
    inline static void ComputeAll(const defaulttype::Vector3 &p, const defaulttype::Vector3 &q, std::array<T,18> &d)
    {
        NormalProjections<6>::Compute<18>(p,q,d);
        NormalProjections<12>::Compute<18>(p,q,d);
        Compute<18>(p,q,d);
    }

};

template<>
struct NormalProjections<24>
{
    template<int K, typename T>
    inline static void Compute(const defaulttype::Vector3 &p, std::array<T,K/2> &d)
    {
        static_assert(K >= 24,"Projection requested is can't be computed.");
        for(size_t i = 0, idx = 9; i < 2; ++i)
            for(size_t j = i+1; j < 3; ++j, ++idx)
                d[idx] = p[i]+p[j]-p[3-(i+j)];
    }

    template<int K, typename T>
    inline static void Compute(const defaulttype::Vector3 &p, std::array<T,K> &d)
    {
        static_assert(K >= 24,"Projection requested is can't be computed.");
        for(size_t i = 0, idx = 9; i < 2; ++i)
            for(size_t j = i+1; j < 3; ++j, ++idx)
                d[idx] = d[idx+K/2] = p[i]+p[j]-p[3-(i+j)];
    }

    template<int K, typename T>
    inline static void Compute(const defaulttype::Vector3 &p, const defaulttype::Vector3 &q, std::array<T,K> &d)
    {
        static_assert(K >= 24,"Projection requested is can't be computed.");
        std::pair<T,T> minmax;
        for(size_t i = 0, idx = 9; i < 2; ++i)
            for(size_t j = i+1; j < 3; ++j)
            {
                minmax = std::minmax(p[i]+p[j]-p[3-(i+j)],q[i]+q[j]-q[3-(i+j)]);
                d[idx] = minmax.first;
                d[idx+K/2] = minmax.second;
            }
    }

    template<typename T>
    inline static void ComputeAll(const defaulttype::Vector3 &p,  std::array<T,12> &d)
    {
        NormalProjections<6>::Compute<24>(p,d);
        NormalProjections<12>::Compute<24>(p,d);
        NormalProjections<18>::Compute<24>(p,d);
        Compute<24>(p,d);
    }

    template<typename T>
    inline static void ComputeAll(const defaulttype::Vector3 &p,  std::array<T,24> &d)
    {
        NormalProjections<6>::Compute<24>(p,d);
        NormalProjections<12>::Compute<24>(p,d);
        NormalProjections<18>::Compute<24>(p,d);
        Compute<24>(p,d);
    }

    template<typename T>
    inline static void ComputeAll(const defaulttype::Vector3 &p, const defaulttype::Vector3 &q, std::array<T,24> &d)
    {
        NormalProjections<6>::Compute<24>(p,q,d);
        NormalProjections<12>::Compute<24>(p,q,d);
        NormalProjections<18>::Compute<24>(p,q,d);
        Compute<24>(p,q,d);
    }

};

} // namespace collision

} // namespace component

} // namespace sofa

#endif

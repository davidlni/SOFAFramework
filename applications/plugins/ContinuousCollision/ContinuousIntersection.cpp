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
#include <iostream>
#include <algorithm>

#include "initContinuousCollision.h"
#include <sofa/core/collision/Intersection.inl>
#include <sofa/core/ObjectFactory.h>
#include <sofa/helper/FnDispatcher.inl>
#include <sofa/helper/proximity.h>
#include <sofa/helper/system/config.h>

#include "ContinuousIntersection.h"
#include "PolytopeModel.h"

namespace sofa
{

namespace component
{

namespace collision
{

using namespace sofa::defaulttype;
using namespace sofa::core::collision;

SOFA_DECL_CLASS(ContinuousIntersection)

int ContinuousIntersectionClass = core::RegisterObject("ContinuousIntersectionClass")
                                  .add< ContinuousIntersection >()
                                  ;


struct NonAdjacentPair
{
private:
    RTriangle elementPair[2];

public:
    NonAdjacentPair(const RTriangle &t1, const RTriangle &t2)
    {
        if (t1.getIndex() > t2.getIndex())
        {
            elementPair[0] = t1;
            elementPair[1] = t2;
        }
        else
        {
            elementPair[0] = t2;
            elementPair[1] = t1;
        }
    }

    void GetParameters(RTriangle &t1, RTriangle &t2)
    {
        t1 = elementPair[0];
        t2 = elementPair[1];
    }

    bool operator < (const NonAdjacentPair &other) {
        if (elementPair[0].getIndex() == other.elementPair[0].getIndex())
            return elementPair[1].getIndex() < other.elementPair[1].getIndex();
        else
            return elementPair[0].getIndex() < other.elementPair[0].getIndex();
    }

    const RTriangle &getElement1() const
    {
      return this->elementPair[0];
    }

    const RTriangle &getElement2() const
    {
      return this->elementPair[1];
    }
};

class ContinuousIntersection::Impl
{
public:
    typedef Vector3::value_type Real;
    typedef RTriangle::index_type index_type;
    typedef RTriangleModel::FeaturePair FeaturePair;
    void processNonAdjacentList(const double &dt, OutputVector*);
    void processOrphans(OutputVector*);
    bool testFeatures(const double &dt, const RTriangle &face1, const RTriangle &face2, OutputVector *output);
    float intersectVertexFace(const double &dt, const RVertex &v, const RTriangle &f, OutputVector *output);
    float intersectVertexFace(const double &dt, const RVertex &v, const RTriangle &f1, const RTriangle &f2, OutputVector *output);

    float intersectEdgeEdge(const double &dt, const REdge &e1, const REdge &e2, OutputVector *output);
    float intersectEdgeEdge(const double &dt, const REdge &e1, const REdge &e2, const RTriangle &f1, const RTriangle &f2, OutputVector *output);
    bool testVertexFace(const double &dt, const RVertex &vertex, const RTriangle &face);
    bool testCoplanarity(const Vector3 &p0, const Vector3 &p1, const Vector3 &p2, const Vector3 &p3,
                         const Vector3 &p4, const Vector3 &p5, const Vector3 &p6, const Vector3 &p7 );
    float vertexFaceIntersection(const double &dt, const RVertex &vertex, const RTriangle &face, OutputVector *output);
    void computeCubicCoefficientsVertexFace(
        const Vector3 &a0, const Vector3 &ad, const Vector3 &b0, const Vector3 &bd,
        const Vector3 &c0, const Vector3 &cd, const Vector3 &p0, const Vector3 &pd,
        helper::fixed_array<Real,4> &coeff);
    void computeCubicCoefficientsEdgeEdge(
        const Vector3 &a0, const Vector3 &ad, const Vector3 &b0, const Vector3 &bd,
        const Vector3 &c0, const Vector3 &cd, const Vector3 &d0, const Vector3 &dd,
        helper::fixed_array<Real,4> &coeff);
    void computeCubicCoefficients(
        const Vector3 &a, const Vector3 &b, const Vector3 &c, const Vector3 &d,
        const Vector3 &e, const Vector3 &f, helper::fixed_array<Real,4> &coeff);
    bool solveCubicWithIntervalNewton(Real &l, Real &r, bool bVF, NewtonCheckData &data, helper::fixed_array<Real,4> &coeff, Vector3 &point, Vector3 &normal);
    bool insideTriangle(const Vector3 &a, const Vector3 &b, const Vector3 &c, const Vector3 &p, Vector3 &point, Vector3 &n);
    bool edgeEdgeIntersect(const Vector3 &p1, const Vector3 &p2, const Vector3 &p3, const Vector3 &p4, Vector3 &point, Vector3 &n);
    bool testEdgeEdge(const double &dt, const REdge &e1, const REdge &e2);
    float edgeEdgeIntersection(const double &dt, const REdge &edge1, const REdge &edge2, OutputVector *output);
    void testOrphans(const double &dt, RTriangleModel* model, OutputVector *output);

    helper::vector<NonAdjacentPair> nonAdjacentPairs;
};

float ContinuousIntersection::Impl::intersectEdgeEdge(const double &dt, const REdge &e1, const REdge &e2, OutputVector *output)
{
    PolytopeModel::DOPType box1(e1.p1(),e1.p2());
    box1+=PolytopeModel::DOPType(e1.p1()+e1.v1()*dt,e1.p2()+e1.v2()*dt);
    PolytopeModel::DOPType box2(e2.p1(),e2.p2());
    box1+=PolytopeModel::DOPType(e2.p1()+e2.v1()*dt,e2.p2()+e2.v2()*dt);

//     if(e1.getBox() != box1)
//     {
//         std::cout << "Edge-Edge" << std::endl;
//         std::cout << e1.getBox();
//         std::cout << box1;
//         std::cout << std::endl;
//     }
//     if(e2.getBox() != box2)
//     {
//         std::cout << "Edge-Edge" << std::endl;
//         std::cout << e2.getBox();
//         std::cout << box2;
//         std::cout << std::endl;
//     }

    if (!box1.Overlaps(box2))
//     if (!e1.getBox().Overlaps(e2.getBox()))
        return -1.f;

    if (this->testEdgeEdge(dt, e1, e2))
        return this->edgeEdgeIntersection(dt,e1,e2,output);

    return -1.f;
}

void ContinuousIntersection::Impl::testOrphans(const double &dt, RTriangleModel* model, OutputVector *output)
{
    for (typename helper::set<FeaturePair>::iterator i=model->getEdgeEdgePairs().begin(), end=model->getEdgeEdgePairs().end(); i!=end; ++i)
    {
        index_type e1, e2;
        i->GetParameters(e1,e2);
        this->intersectEdgeEdge(dt,REdge(model,e1),REdge(model,e2),output);
    }
    for (typename helper::set<FeaturePair>::iterator i=model->getVertexFacePairs().begin(), end=model->getVertexFacePairs().end(); i!=end; ++i)
    {
        index_type v, f;
        i->GetParameters(f,v);
        this->intersectVertexFace(dt,RVertex(model,v),RTriangle(model,f),output);
    }
}

float ContinuousIntersection::Impl::edgeEdgeIntersection(const double &dt, const REdge &edge1, const REdge &edge2, OutputVector *output)
{
    /* Default value returned if no collision occurs */
    float collisionTime = -1.0f;
    Real epsilon = std::numeric_limits<Real>::epsilon();

    /* diff. vectors for linear interpolation: x1-x0 == dt*v0 */
    Vector3 dd = edge2.v2()*dt;
    Vector3 ad = edge1.v1()*dt;
    Vector3 bd = edge1.v2()*dt;
    Vector3 cd = edge2.v1()*dt;

    /*
    * Compute scalar coefficients by evaluating dot and cross-products.
    */
    helper::fixed_array<Real,4> coeffs; /* cubic polynomial coefficients */
    this->computeCubicCoefficientsEdgeEdge(edge1.p1(), ad, edge1.p2(), bd, edge2.p1(), cd, edge2.p2(), dd, coeffs);

    if (std::fabs(coeffs[0]) <= epsilon &&
            std::fabs(coeffs[1]) <= epsilon &&
            std::fabs(coeffs[2]) <= epsilon &&
            std::fabs(coeffs[3]) <= epsilon)
        return collisionTime;

    NewtonCheckData data;
    data.a0 = edge1.p1(), data.b0 = edge1.p2();
    data.c0 = edge2.p1(), data.p0 = edge2.p2();
    data.ad = ad, data.bd = bd;
    data.cd = cd, data.pd = dd;

    /*
    * iteratively solve the cubic (scalar) equation and test for validity of the solution.
    */
    Real l = 0;
    Real r = 1;
    Vector3 point,n;
    if (this->solveCubicWithIntervalNewton(l, r, false, data, coeffs, point,n))
    {
        collisionTime = (l+r)*0.5f;
        output->push_back(DetectionOutput());
        DetectionOutput &detection = output->back();

        detection.point[0] = point;
        detection.point[1] = point;
        detection.normal = -n;
        detection.value = 0;
        detection.elem.first = edge1.getTriangle(0);
        detection.elem.second = edge2.getTriangle(0);
        detection.id = edge1.getIndex();
        detection.deltaT = collisionTime;

        std::cout << "EdgeEdge - Normal = " << -n << std::endl;

    }

    return collisionTime;
}

bool ContinuousIntersection::Impl::testEdgeEdge(const double &dt, const REdge &e1, const REdge &e2)
{
    RVertex pe1(e1.model,e1.p1Index());
    RVertex pe2(e1.model,e1.p2Index());
    RVertex pe3(e2.model,e2.p1Index());
    RVertex pe4(e2.model,e2.p2Index());

    return this->testCoplanarity(pe1.p(), pe2.p(), pe3.p(), pe4.p(),
                                 pe1.p()+pe1.v()*dt, pe2.p()+pe2.v()*dt, pe3.p()+pe3.v()*dt, pe4.p()+pe4.v()*dt);
}

float ContinuousIntersection::Impl::intersectEdgeEdge(const double &dt, const REdge &e1, const REdge &e2, const RTriangle &f1, const RTriangle &f2, OutputVector *output)
{
  PolytopeModel::DOPType box1(e1.p1(),e1.p2());
  box1+=PolytopeModel::DOPType(e1.p1()+e1.v1()*dt,e1.p2()+e1.v2()*dt);
  PolytopeModel::DOPType box2(e2.p1(),e2.p2());
  box1+=PolytopeModel::DOPType(e2.p1()+e2.v1()*dt,e2.p2()+e2.v2()*dt);

  if (!box1.Overlaps(box2))
    return -1.f;

    REdge e[2] = {e1,e2};
    RTriangle f[2] = {f1,f2};

    if (e1.getIndex() > e2.getIndex())
    {
        e[0] = e2;
        e[1] = e1;
        f[0] = f2;
        f[1] = f1;
    }

    for (int i=0; i<2; ++i)
    {
        RTriangle adjacentTriangle1 = e[0].getTriangle(i);
        if(adjacentTriangle1.valid())
        for (int j=0; j<2; ++j)
        {
            RTriangle adjacentTriangle2 = e[1].getTriangle(j);
            if (!adjacentTriangle2.valid())
                continue;

            if (!adjacentTriangle1.covertices(adjacentTriangle2))
            {
              if (adjacentTriangle1 == f[0] && adjacentTriangle2 == f[1])
                {
                    if (this->testEdgeEdge(dt, e1, e2))
                        return this->edgeEdgeIntersection(dt,e1,e2,output);
                    return -1.f;
                }
                return -1.f;
            }
        }
    }

    return -1.f;
}

bool ContinuousIntersection::Impl::edgeEdgeIntersect(const Vector3 &p1, const Vector3 &p2, const Vector3 &p3, const Vector3 &p4, Vector3 &point, Vector3 &n)
{
    Real epsilon = std::numeric_limits<Real>::epsilon();

    Vector3 p21 = p2 - p1;
    Vector3 p43 = p4 - p3;
    if(p21.norm2() < epsilon || p43.norm2() < epsilon)
//     if (fabs(p43[0])  < epsilon &&
//             fabs(p43[1])  < epsilon &&
//             fabs(p43[2])  < epsilon)
        return false;

//     if(p21.norm2() < epsilon)
//     if (fabs(p21[0])  < epsilon &&
//             fabs(p21[1])  < epsilon &&
//             fabs(p21[2])  < epsilon)
//         return false;

    Vector3 p13 = p1 - p3;

    Real d4321 = p43*p21;
    Real d4343 = p43.norm2();
    Real d2121 = p21.norm2();
    Real denom = d2121 * d4343 - d4321 * d4321;

    if (fabs(denom) < epsilon)
        return false;

    Real d1321 = p13*p21;
    Real d1343 = p13*p43;
    Real numer = d1343 * d4321 - d1321 * d4343;

    Real mua = numer / denom;
    if (mua < 0 || mua > 1)
        return false;

    Real mub = (d1343 + d4321 * mua) / d4343;
    if (mub < 0 || mub > 1)
        return false;

    // Compute intersection points and normal
    point = p1 + p21*mua;
    Vector3 pb = p3 + p43*mub;
    n = p21.cross(p43);
    n.normalize();
    return true;
}

bool ContinuousIntersection::Impl::insideTriangle(const Vector3 &a, const Vector3 &b, const Vector3 &c, const Vector3 &p, Vector3 &point, Vector3 &n)
{
    Vector3 da = a - p;
    Vector3 db = b - p;
    Vector3 dc = c - p;

    n = (b-a).cross(c-a);

    Real wa = db.cross(dc)*n;
    Real wb = dc.cross(da)*n;
    Real wc = da.cross(db)*n;

    if (wa < Real(0) || wb < Real(0) || wc < Real(0))
        return false;

    //Compute barycentric coordinates
    Real area2 = n.norm2();
    assert(area2 > 0);
    Real invArea2 = 1.0/area2;
    wa *= invArea2, wb *= invArea2, wc *= invArea2;

    point = Vector3(wa, wb, wc);

    n.normalize();
    return true;
}

void ContinuousIntersection::Impl::computeCubicCoefficientsVertexFace(
    const Vector3 &a0, const Vector3 &ad, const Vector3 &b0, const Vector3 &bd,
    const Vector3 &c0, const Vector3 &cd, const Vector3 &p0, const Vector3 &pd,
    helper::fixed_array<Real,4> &coeff)
{
    this->computeCubicCoefficients(bd-ad,cd-ad,pd-ad,c0-a0,b0-a0,p0-a0,coeff);
}

void ContinuousIntersection::Impl::computeCubicCoefficientsEdgeEdge(
    const Vector3 &a0, const Vector3 &ad, const Vector3 &b0, const Vector3 &bd,
    const Vector3 &c0, const Vector3 &cd, const Vector3 &d0, const Vector3 &dd,
    helper::fixed_array<Real,4> &coeff)
{
    this->computeCubicCoefficients(bd-ad,dd-cd,cd-ad,d0-c0,b0-a0,c0-a0,coeff);
}

void ContinuousIntersection::Impl::computeCubicCoefficients(
    const Vector3 &a, const Vector3 &b, const Vector3 &c, const Vector3 &d,
    const Vector3 &e, const Vector3 &f, helper::fixed_array<Real,4> &coeff)
{
    Vector3 v1 = a.cross(b);
    Vector3 v2 = a.cross(d);
    Vector3 v3 = e.cross(b);
    Vector3 v4 = e.cross(d);

    coeff[0] = f*v4;
    coeff[1] = c*v4 + f*(v2 + v3);
    coeff[2] = f*v1 + c*(v2 + v3);
    coeff[3] = c*v1;
}

bool ContinuousIntersection::Impl::solveCubicWithIntervalNewton(Real &l, Real &r, bool bVF, NewtonCheckData &data, helper::fixed_array<Real,4> &coeff, Vector3 &point, Vector3 &n)
{
    Real v2[2] = {l*l,r*r};
    Real v[2] = {l,r};

    unsigned char min3, min2, min1, max3, max2, max1;

    min3=*((unsigned char*)&coeff[3]+7) >> 7;
    max3=min3^1;
    min2=*((unsigned char*)&coeff[2]+7) >> 7;
    max2=min2^1;
    min1=*((unsigned char*)&coeff[1]+7) >> 7;
    max1=min1^1;

    // Compute bounds for the cubic
    Real bounds[2] = { coeff[3]*v2[min3]*v[min3]+coeff[2]*v2[min2]+coeff[1]*v[min1]+coeff[0],
                       coeff[3]*v2[max3]*v[max3]+coeff[2]*v2[max2]+coeff[1]*v[max1]+coeff[0]
                     };

    if (bounds[0]>0 || bounds[1]<0) return false;

    // Starting at the middle of the interval
    Real m = 0.5*(r+l);

    // Compute bounds for the derivative
    Real dbounds[2] = { 3.0*coeff[3]*v2[min3]+2.0*coeff[2]*v[min2]+coeff[1],
                        3.0*coeff[3]*v2[max3]+2.0*coeff[2]*v[max2]+coeff[1]
                      };

    if (dbounds[0] > 0 || dbounds[1] < 0) // we can use Newton
    {
        Real m2 = m*m;
        Real fm = coeff[3]*m2*m+coeff[2]*m2+coeff[1]*m+coeff[0];
        Real nl = m;
        Real nu = m;
        if (fm > 0)
        {
            nl -= fm*(1.0/dbounds[0]);
            nu -= fm*(1.0/dbounds[1]);
        }
        else
        {
            nu -= fm*(1.0/dbounds[0]);
            nl -= fm*(1.0/dbounds[1]);
        }

        // Intersect with [l,r]
        if (nu<l || nl>r)
            return false; // outside the interval

        if (nl>l)
        {
            if (nu<r)
            {
                l = nl;
                r = nu;
                m = 0.5*(l+r);
            }
            else
            {
                l = nl;
                m = 0.5*(l+r);
            }
        }
        else
        {
            if (nu<r)
            {
                r = nu;
                m = 0.5*(l+r);
            }
        }
    }

    // Check root validity
    if ((r-l)< Real(1e-10))
    {
        if (bVF)
            return this->insideTriangle(data.ad*r + data.a0, data.bd*r + data.b0, data.cd*r + data.c0, data.pd*r + data.p0,point,n);
        else
            return this->edgeEdgeIntersect(data.ad*r + data.a0, data.bd*r + data.b0, data.cd*r + data.c0, data.pd*r + data.p0,point,n);
    }

    Real rBkUp = r;
    r = m;
    if (this->solveCubicWithIntervalNewton(l,r,bVF,data,coeff,point,n))
        return true;
    l = m, r = rBkUp;
    return this->solveCubicWithIntervalNewton(l,r,bVF,data,coeff,point,n);
}

float ContinuousIntersection::Impl::vertexFaceIntersection(const double &dt, const RVertex &vertex, const RTriangle &face, OutputVector *output)
{
    /* Default value returned if no collision occurs */
    float collisionTime = -1.0f;
    Real epsilon = std::numeric_limits<Real>::epsilon();

    /* diff. vectors for linear interpolation: x1-x0 == dt*v0 */
    Vector3 ad = face.v1()*dt;
    Vector3 bd = face.v2()*dt;
    Vector3 cd = face.v3()*dt;
    Vector3 qd = vertex.v()*dt;

    /*
     * Compute scalar coefficients by evaluating dot and cross-products.
     */
    helper::fixed_array<Real,4> coeffs; /* cubic polynomial coefficients */
    this->computeCubicCoefficientsVertexFace(face.p1(), ad, face.p2(), bd, face.p3(), cd, vertex.p(), qd, coeffs);

    if (std::fabs(coeffs[0]) <= epsilon &&
            std::fabs(coeffs[1]) <= epsilon &&
            std::fabs(coeffs[2]) <= epsilon &&
            std::fabs(coeffs[3]) <= epsilon)
        return collisionTime;

    NewtonCheckData data;
    data.a0 = face.p1(), data.b0 = face.p2();
    data.c0 = face.p3(), data.p0 = vertex.p();
    data.ad = ad, data.bd = bd;
    data.cd = cd, data.pd = qd;

    /*
     * iteratively solve the cubic (scalar) equation and test for validity of the solution.
     */
    Real l = 0;
    Real r = 1;
    Vector3 barycenter, n;
    if (this->solveCubicWithIntervalNewton(l, r, true, data, coeffs,barycenter,n))
    {
        collisionTime = (l+r)*0.5;
        Vector3 qi = qd*collisionTime+vertex.p();
        output->push_back(DetectionOutput());
        DetectionOutput &detection = output->back();

        detection.point[0] = qi;
        detection.point[1] = qi;
        detection.normal = -n;
        detection.value = 0;
        detection.elem.first = face;
        detection.elem.second = RTriangle(vertex.model, vertex.getTriangles()[0]);
        detection.id = face.getIndex();
        detection.deltaT = collisionTime;
    }

    return collisionTime;
}

bool ContinuousIntersection::Impl::testVertexFace(const double &dt, const RVertex &vertex, const RTriangle &face)
{
//   std::cout << face.p1() << std::endl;
//   std::cout << face.p2() << std::endl;
//   std::cout << face.p3() << std::endl;
//   std::cout << vertex.p() << std::endl;
//   std::cout << face.p1()+face.v1()*dt << std::endl;
//   std::cout << face.p2()+face.v2()*dt << std::endl;
//   std::cout << face.p3()+face.v3()*dt << std::endl;
//   std::cout << vertex.p()+vertex.v()*dt << std::endl;
    return this->testCoplanarity(face.p1(), face.p2(), face.p3(), vertex.p(),face.p1()+face.v1()*dt,
                                 face.p2()+face.v2()*dt, face.p3()+face.v3()*dt,vertex.p()+vertex.v()*dt);
}

bool ContinuousIntersection::Impl::testCoplanarity(const Vector3 &p0, const Vector3 &p1, const Vector3 &p2, const Vector3 &p3,
        const Vector3 &p4, const Vector3 &p5, const Vector3 &p6, const Vector3 &p7 )
{
    Vector3 n0 = (p1-p0).cross(p2-p1);
    Vector3 n1 = (p5-p4).cross(p6-p5);
    Vector3 delta = ((p5-p1)-(p4-p0)).cross((p6-p2)-(p5-p1));
    Vector3 nX = (n0+n1-delta)*0.5;

    Vector3 pa0 = p3-p0;
    Vector3 pa1 = p7-p4;

    Real A = n0*pa0;
    Real B = n1*pa1;
    Real C = nX*pa0;
    Real D = nX*pa1;
    Real E = n1*pa0;
    Real F = n0*pa1;

    if (A > 0 && B > 0 && (2*C+F) > 0 && (2*D+E) > 0)
        return false;

    if (A < 0 && B < 0 && (2*C+F) < 0 && (2*D+E) < 0)
        return false;

    return true;
}

float ContinuousIntersection::Impl::intersectVertexFace(const double &dt, const RVertex &vertex, const RTriangle &face, OutputVector *output)
{
    PolytopeModel::DOPType box1(vertex.p(),vertex.p()+vertex.v()*dt);
    PolytopeModel::DOPType box2(face.p1(),face.p2(),face.p3());
    box2 += PolytopeModel::DOPType(face.p1()+face.v1()*dt,face.p2()+face.v1()*dt,face.p3()+face.v1()*dt);
    if (!box1.Overlaps(box2))
//   if (!vertex.getBox().Overlaps(face.getBox()))
    return -1.f;

  if (this->testVertexFace(dt,vertex,face))
    return this->vertexFaceIntersection(dt,vertex,face,output);

  return -1.f;
}

float ContinuousIntersection::Impl::intersectVertexFace(const double &dt, const RVertex &vertex, const RTriangle &face1, const RTriangle &face2, OutputVector *output)
{
    PolytopeModel::DOPType box1(vertex.p(),vertex.p()+vertex.v()*dt);
    PolytopeModel::DOPType box2(face1.p1(),face1.p2(),face1.p3());
    box2 += PolytopeModel::DOPType(face1.p1()+face1.v1()*dt,face1.p2()+face1.v1()*dt,face1.p3()+face1.v1()*dt);

//     if(vertex.getBox() != box1)
//     {
//         std::cout << "Vertex-Face" << std::endl;
//         std::cout << vertex.getBox();
//         std::cout << box1;
//         std::cout << std::endl;
//     }
//     if(face1.getBox() != box2)
//     {
//         std::cout << "Vertex-Face" << std::endl;
//         std::cout << face1.getBox();
//         std::cout << box2;
//         std::cout << std::endl;
//     }
//     std::cout.precision(16);
//     std::cout << vertex.p() << std::endl;
//     std::cout << vertex.p()+vertex.v()*dt<< std::endl;
//     std::cout << face1.p1() << std::endl;
//     std::cout << face1.p2() << std::endl;
//     std::cout << face1.p3() << std::endl;

    if (!box1.Overlaps(box2))
        return -1.f;

    typename helper::vector<RVertex::index_type>::const_iterator i, end = vertex.getTriangles().end();
    for (i = vertex.getTriangles().begin(); i != end; ++i)
    {
        RTriangle t(vertex.model,*i);
        if (!face1.covertices(t))
        {
            if (t == face2)
            {
                if (this->testVertexFace(dt,vertex,face1))
                    return this->vertexFaceIntersection(dt,vertex,face1,output);
                return -1.f;
            }
            return -1.f;
        }
    }
    return -1.f;
}

void ContinuousIntersection::Impl::processNonAdjacentList(const double &dt, OutputVector* output)
{
    helper::vector<NonAdjacentPair>::iterator i = this->nonAdjacentPairs.begin();
    for(helper::vector<NonAdjacentPair>::iterator end = this->nonAdjacentPairs.end(); i != end; ++i)
        this->testFeatures(dt,i->getElement1(),i->getElement2(),output);
}

bool ContinuousIntersection::Impl::testFeatures(const double &dt, const RTriangle &face1, const RTriangle &face2, OutputVector *output)
{
    // 6 VF test
    RVertex v1(face1.model,face1.p1Index());
    RVertex v2(face1.model,face1.p2Index());
    RVertex v3(face1.model,face1.p3Index());
    RVertex v4(face2.model,face2.p1Index());
    RVertex v5(face2.model,face2.p2Index());
    RVertex v6(face2.model,face2.p3Index());
//     std::cout << "Testing vertex " << v1.getIndex() << " against face " << face2.model->getTriangles()[face2.getIndex()] << std::endl;
//     std::cout << "output.size() = " << output->size() << std::endl;
//     std::cout << "Testing vertex " << v2.getIndex() << " against face " << face2.model->getTriangles()[face2.getIndex()] << std::endl;
//     std::cout << "output.size() = " << output->size() << std::endl;
//     std::cout << "Testing vertex " << v3.getIndex() << " against face " << face2.model->getTriangles()[face2.getIndex()] << std::endl;
//     std::cout << "output.size() = " << output->size() << std::endl;
    this->intersectVertexFace(dt,v1,face2,face1,output);
    this->intersectVertexFace(dt,v2,face2,face1,output);
    this->intersectVertexFace(dt,v3,face2,face1,output);
        std::cout << "output.size() = " << output->size() << std::endl;

//     std::cout << "Testing vertex " << v4.getIndex() << " against face " << face1.model->getTriangles()[face1.getIndex()] << std::endl;
//     std::cout << "output.size() = " << output->size() << std::endl;
//     std::cout << "Testing vertex " << v5.getIndex() << " against face " << face1.model->getTriangles()[face1.getIndex()] << std::endl;
//     std::cout << "output.size() = " << output->size() << std::endl;
//     std::cout << "Testing vertex " << v6.getIndex() << " against face " << face1.model->getTriangles()[face1.getIndex()] << std::endl;
//     std::cout << "output.size() = " << output->size() << std::endl;
    this->intersectVertexFace(dt,v4,face1,face2,output);
    this->intersectVertexFace(dt,v5,face1,face2,output);
    this->intersectVertexFace(dt,v6,face1,face2,output);
        std::cout << "output.size() = " << output->size() << std::endl;

    // 9 EE test
    REdge e(face1.model,face1.e1Index()), e1(face2.model,face2.e1Index()), e2(face2.model,face2.e2Index()), e3(face2.model,face2.e3Index());
//     std::cout << "Testing edge [" << face1.model->getEdges()[e.getIndex()] << "] against edge [" << face2.model->getEdges()[e1.getIndex()] << "]"<< std::endl;
//     std::cout << "output.size() = " << output->size() << std::endl;
//     std::cout << "Testing edge [" << face1.model->getEdges()[e.getIndex()] << "] against edge [" << face2.model->getEdges()[e2.getIndex()] << "]"<< std::endl;
//     std::cout << "output.size() = " << output->size() << std::endl;
//     std::cout << "Testing edge [" << face1.model->getEdges()[e.getIndex()] << "] against edge [" << face2.model->getEdges()[e3.getIndex()] << "]"<< std::endl;
//     std::cout << "output.size() = " << output->size() << std::endl;
    this->intersectEdgeEdge(dt,e,e1,face1,face2,output);
    this->intersectEdgeEdge(dt,e,e2,face1,face2,output);
    this->intersectEdgeEdge(dt,e,e3,face1,face2,output);
        std::cout << "output.size() = " << output->size() << std::endl;

    e = REdge(face1.model,face1.e2Index());
//     std::cout << "Testing edge [" << face1.model->getEdges()[e.getIndex()] << "] against edge [" << face2.model->getEdges()[e1.getIndex()] << "]"<< std::endl;
//     std::cout << "output.size() = " << output->size() << std::endl;
//     std::cout << "Testing edge [" << face1.model->getEdges()[e.getIndex()] << "] against edge [" << face2.model->getEdges()[e2.getIndex()] << "]"<< std::endl;
//     std::cout << "output.size() = " << output->size() << std::endl;
//     std::cout << "Testing edge [" << face1.model->getEdges()[e.getIndex()] << "] against edge [" << face2.model->getEdges()[e3.getIndex()] << "]"<< std::endl;
//     std::cout << "output.size() = " << output->size() << std::endl;
    this->intersectEdgeEdge(dt,e,e1,face1,face2,output);
    this->intersectEdgeEdge(dt,e,e2,face1,face2,output);
    this->intersectEdgeEdge(dt,e,e3,face1,face2,output);
        std::cout << "output.size() = " << output->size() << std::endl;

    e = REdge(face1.model,face1.e3Index());
//     std::cout << "Testing edge [" << face1.model->getEdges()[e.getIndex()] << "] against edge [" << face2.model->getEdges()[e1.getIndex()] << "]"<< std::endl;
//     std::cout << "output.size() = " << output->size() << std::endl;
//     std::cout << "Testing edge [" << face1.model->getEdges()[e.getIndex()] << "] against edge [" << face2.model->getEdges()[e2.getIndex()] << "]"<< std::endl;
//     std::cout << "output.size() = " << output->size() << std::endl;
//     std::cout << "Testing edge [" << face1.model->getEdges()[e.getIndex()] << "] against edge [" << face2.model->getEdges()[e3.getIndex()] << "]"<< std::endl;
//     std::cout << "output.size() = " << output->size() << std::endl;
    this->intersectEdgeEdge(dt,e,e1,face1,face2,output);
    this->intersectEdgeEdge(dt,e,e2,face1,face2,output);
    this->intersectEdgeEdge(dt,e,e3,face1,face2,output);
        std::cout << "output.size() = " << output->size() << std::endl;
    return true;
}

ContinuousIntersection::ContinuousIntersection()
{
    this->intersectors.add<PolytopeModel,PolytopeModel,ContinuousIntersection> (this);
    this->intersectors.add<RTriangleModel,RTriangleModel,ContinuousIntersection> (this);
    this->pimpl = new Impl;
}

ContinuousIntersection::~ContinuousIntersection()
{
    delete this->pimpl;
}

int ContinuousIntersection::computeIntersection(const RTriangle& t1, const RTriangle& t2, OutputVector*)
{
    if(!t1.covertices(t2))
    {
      PolytopeModel::DOPType box1(t1.p1(),t1.p2(),t1.p3());
      box1+=PolytopeModel::DOPType(t1.p1()+t1.v1(),t1.p2()+t1.v2(),t1.p3()+t1.v3());
      PolytopeModel::DOPType box2(t2.p1(),t2.p2(),t2.p3());
      box2+=PolytopeModel::DOPType(t2.p1()+t2.v1(),t2.p2()+t2.v2(),t2.p3()+t2.v3());
//       if(t1.getBox().Overlaps(t2.getBox()))
      if(box1.Overlaps(box2))
            this->pimpl->nonAdjacentPairs.push_back(NonAdjacentPair(t1,t2));
    }
    return 0;
}

// int ContinuousIntersection::computeIntersection(const Polytope& t1, const RTriangle& t2, OutputVector*)
// {
//     if(this->testIntersection(t1,t2))
// //       this->pimpl->nonAdjacentPairs.push_back(NonAdjacentPair(t1,t2));
//     return 0;
// }

void ContinuousIntersection::beginBroadPhase()
{
    this->pimpl->nonAdjacentPairs.clear();
}

int ContinuousIntersection::beginIntersection(sofa::core::CollisionModel* model1, sofa::core::CollisionModel* model2, OutputVector* )
{
    if(model1->getLast() == model2->getLast())
    {
      static_cast<RTriangleModel*>(model1->getLast())->bufferAdjacentLists();
      static_cast<RTriangleModel*>(model1->getLast())->setOrphans();
    }
    return 0;
}

int ContinuousIntersection::endIntersection(sofa::core::CollisionModel* model1, sofa::core::CollisionModel* model2, OutputVector* contacts)
{
  double dt = this->getContext()->getDt();
  this->pimpl->processNonAdjacentList(dt,contacts);
//   sout << "contacts.size()" << contacts->size() << std::endl;
  if(model1->getLast() == model2->getLast())
  {
    this->pimpl->testOrphans(dt,static_cast<RTriangleModel*>(model1->getLast()),contacts);
  }
  return 0;
}

/// Return the intersector class handling the given pair of collision models, or NULL if not supported.
ElementIntersector* ContinuousIntersection::findIntersector(core::CollisionModel* object1, core::CollisionModel* object2, bool& swapModels)
{
    return intersectors.get(object1, object2, swapModels);
}


} // namespace collision

} // namespace component

namespace core
{
namespace collision
{
template class SOFA_CONTINUOUS_COLLISION_API IntersectorFactory<component::collision::ContinuousIntersection>;
}
}

} // namespace sofa


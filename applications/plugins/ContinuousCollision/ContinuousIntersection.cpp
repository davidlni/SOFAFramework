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

#include <sofa/core/collision/Intersection.inl>
#include <sofa/core/ObjectFactory.h>
#include <sofa/helper/FnDispatcher.inl>
#include <sofa/helper/proximity.h>
#include <sofa/helper/system/config.h>

#include "ContinuousIntersection.h"
#include "PolytopeModel.h"
#include "RTriangleModel.h"

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
  core::CollisionElementIterator elementPair[2];

public:
  NonAdjacentPair(const RTriangle &t1, const RTriangle &t2)
    {
        if (t1 > t2) {
            elementPair[0] = t1;
            elementPair[1] = t2;
        } else {
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
};

class ContinuousIntersection::Impl
{
public:
  typedef Vector3::value_type Real;
  Impl(ContinuousIntersection* external);
  void processNonAdjacentList(OutputVector*);
  void processOrphans(OutputVector*);
  bool testFeatures(const double &dt, const RTriangle &face1, const RTriangle &face2, OutputVector *output);
  float intersectVertexFace(const double &dt, const RVertex &v, const RTriangle &f, OutputVector *output);
  float intersectVertexFace(const double &dt, const RVertex &v, const RTriangle &f1, const RTriangle &f2, OutputVector *output);

  float intersectEdgeEdge(const double &dt, const REdge &e1, const REdge &e2, OutputVector *output);
  float intersectEdgeEdge(const double &dt, const REdge &e1, const REdge &e2, const RTriangle &f1, const RTriangle &f2, OutputVector *output);
  bool testVertexFace(const double &dt, const RVertex &vertex, const RTriangle &face);
  bool testCoplanarity(const Vector3 &p0, const Vector3 &p1, const Vector3 &p2, const Vector3 &p3,
        const Vector3 &p4, const Vector3 &p5, const Vector3 &p6, const Vector3 &p7 );
  float vertexFaceIntersection(const double &dt, const VertexElement &vertex, const Element &face, OutputVector *output);
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
  bool solveCubicWithIntervalNewton(Real &l, Real &r, bool bVF, NewtonCheckData &data, helper::fixed_array<Real,4> &coeff, Vector3 &point);
  bool insideTriangle(const Vector3 &a, const Vector3 &b, const Vector3 &c, const Vector3 &p, Vector3 &n);
  bool lineLineIntersect(const Vector3 &p1, const Vector3 &p2, const Vector3 &p3, const Vector3 &p4, Vector3 &n);
  bool testEdgeEdge(const double &dt, const REdge &e1, const REdge &e2);
  float edgeEdgeIntersection(const double &dt, const REdge &edge1, const REdge &edge2, OutputVector *output);
  void testOrphans(const double &dt, OutputVector *output);
private:
  helper::vector<NonAdjacentPair> nonAdjacentPairs;
};

float ContinuousIntersection::Impl::intersectVertexFace(const double &dt, const RVertex &vertex, const RTriangle &face, OutputVector *output)
{
    if (!vertex.getBox().Overlaps(face.getBox()))
        return false;

    if (this->testVertexFace(dt,vertex,face))
        return this->vertexFaceIntersection(dt,vertex,face,output);

    return -1.f;
}

float ContinuousIntersection::Impl::intersectEdgeEdge(const double &dt, const REdge &e1, const REdge &e2, OutputVector *output)
{
    if (!e1.getBox().Overlaps(e2.getBox()))
        return -1.f;

    if (this->testEdgeEdge(dt, e1, e2))
        return this->edgeEdgeIntersection(dt,e1,e2,output);

    return -1.f;
}

void ContinuousIntersection::Impl::testOrphans(const double &dt, RTriangleModel* model, OutputVector *output)
{
    typedef RTriangleModel::FeaturePair FeaturePair;
    typedef RTriangleModel::index_type index_type;
    for (typename helper::set<FeaturePair>::iterator i=model->getEdgeEdgePairs().begin(), end=model->getEdgeEdgePairs().end(); i!=end; ++i)
    {
        index_type e1, e2;
        i->GetParameters(e1,e2);
        this->intersectEdgeEdge(dt,REdge(model,e1),REdge(model,e2),output);
    }
    for (typename helper::set<FeaturePair>::iterator i=model->getVertexFacePairs().begin(), end=model->getVertexFacePairs().end(); i!=end; ++i)
    {
        index_type v, f;
        i->GetParameters(f, v);
        this->intersectVertexFace(dt,RVertex(model,v), RTriangle(model,f),output);
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
    if (solveCubicWithIntervalNewton(l, r, false, data, coeffs, point,n))
    {
        collisionTime = (l+r)*0.5f;
	output->push_back(DetectionOutput());
	DetectionOutput &detection = output->back();
	
	detection.point[0] = point;
	detection.point[1] = point;    
	detection.normal = -n;
	detection.value = 0;
	detection.elem.first = edge1;
	detection.elem.second = edge2;
	// WARNING: Does this index need to be unique?
	detection.id = edge1.getIndex();
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

float ContinuousIntersection::Impl::intersectEdgeEdge(const double &dt, const REdge &e1, const REdge &e2, const RTriangle &f1, const RTriangle &f2)
{
  if (!e1.getBox().Overlaps(e2.getBox()))
        return -1.f;

    REdge e[2] = {e1,e2};
    RTriangle f[2] = {f1,f2};

    if (e1.getIndex() > e2.getIndex())
    {
        e[0] = e2, e[1] = e1;
        f[0] = f2, f[1] = f1;
    }

    for (int i=0; i<2; i++)
        for (int j=0; j<2; j++)
        {
            index_type ff1 = e[0].FaceId(i);
            index_type ff2 = e[1].FaceId(j);

            if (ff1 == -1 || ff2 == -1)
                continue;

            RTriangle adjacentFace1(this,ff1);
            RTriangle adjacentFace2(this,ff2);
            if (!adjacentFace1.covertices(adjacentFace2))
            {
              if (adjacentFace1 == f[0] && adjacentFace2 == f[1])
                {
                    if (this->testEdgeEdge(dt, e1, e2))
                        return this->edgeEdgeIntersection(dt,e1,e2);
                    return -1.f;
                }
                return -1.f;
            }
        }

    return -1.f;
}

bool ContinuousIntersection::Impl::lineLineIntersect(const Vector3 &p1, const Vector3 &p2, const Vector3 &p3, const Vector3 &p4, Vector3 &point, Vector3 &n)
{
    Vector3 p13,p43,p21;
    Real d1343,d4321,d1321,d4343,d2121;
    Real numer,denom;
    Real epsilon = std::numeric_limits<Real>::epsilon();

    p13 = p1 - p3;
    p43 = p4 - p3;
    if (fabs(p43[0])  < epsilon &&
            fabs(p43[1])  < epsilon &&
            fabs(p43[2])  < epsilon)
        return false;

    p21 = p2 - p1;
    if (fabs(p21[0])  < epsilon &&
            fabs(p21[1])  < epsilon &&
            fabs(p21[2])  < epsilon)
        return false;

    d1343 = p13*p43;
    d4321 = p43*p21;
    d1321 = p13*p21;
    d4343 = p43*p43;
    d2121 = p21*p21;

    denom = d2121 * d4343 - d4321 * d4321;
    if (fabs(denom) < epsilon)
        return false;
    numer = d1343 * d4321 - d1321 * d4343;

    Real mua = numer / denom;
    if (mua < 0 || mua > 1)
        return false;

    Real mub = (d1343 + d4321 * mua) / d4343;
    if (mub < 0 || mub > 1)
        return false;

    // Do not return these
    point = p1 + p21*mua;
    Vector3 pb = p3 + p43*mub;
    n = point-pb;
    return true;
}


bool ContinuousIntersection::Impl::insideTriangle(const Vector3 &a, const Vector3 &b, const Vector3 &c, const Vector3 &p, Vector3 &point, Vector3 &n)
{
    Vector3 da, db, dc;
    Real wa, wb, wc;

    Vector3 ba = b-a;
    Vector3 ca = c-a;
    n = ba.cross(ca);

    da = a - p, db = b - p, dc = c - p;
    if ((wa = db.cross(dc)*n) < 0.0f) return false;
    if ((wb = dc.cross(da)*n) < 0.0f) return false;
    if ((wc = da.cross(db)*n) < 0.0f) return false;

    //Compute barycentric coordinates
    Real area2 = n.norm2();
    wa /= area2, wb /= area2, wc /= area2;
    
    point = Vector3(wa, wb, wc);

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

bool ContinuousIntersection::Impl::solveCubicWithIntervalNewton(Real &l, Real &r, bool bVF, NewtonCheckData &data, helper::fixed_array<Real,4> &coeff, Vector3 &point)
{
    Real v2[2] = {l*l,r*r};
    Real v[2] = {l,r};
    Real rBkUp;

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

    if (bounds[0]<0) return false;
    if (bounds[1]>0) return false;

    // Starting starts at the middle of the interval
    Real m = 0.5*(r+l);

    // Compute bounds for the derivative
    Real dbounds[2] = { 3.0*coeff[3]*v2[min3]+2.0*coeff[2]*v[min2]+coeff[1],
                        3.0*coeff[3]*v2[max3]+2.0*coeff[2]*v[max2]+coeff[1]
                      };

    if ((dbounds[0] > 0) || (dbounds[1] < 0)) // we can use Newton
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
            return this->insideTriangle(data.ad*r + data.a0, data.bd*r + data.b0, data.cd*r + data.c0, data.pd*r + data.p0,point);
        else
            return this->lineLineIntersect(data.ad*r + data.a0, data.bd*r + data.b0, data.cd*r + data.c0, data.pd*r + data.p0,point);
    }
    rBkUp = r, r = m;
    if (this->solveCubicWithIntervalNewton(l,r, bVF, data, coeff, point)) return true;
    l = m, r = rBkUp;
    return (this->solveCubicWithIntervalNewton(l,r, bVF, data, coeff, point));
}

float ContinuousIntersection::Impl::vertexFaceIntersection(const double &dt, const RVertex &vertex, const RTriangle &face, OutputVector *output)
{
    /* Default value returned if no collision occurs */
    float collisionTime = -1.0f;
    Real epsilon = std::numeric_limits<Real>::epsilon();

    /* diff. vectors for linear interpolation: x1-x0 == dt*v0 */
    Vector3 qd = vertex.v()*dt;
    Vector3 ad = face.v1()*dt;
    Vector3 bd = face.v2()*dt;
    Vector3 cd = face.v3()*dt;

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
    Vertex barycenter, n;
    if (this->solveCubicWithIntervalNewton(l, r, true, data, coeffs,barycenter,n))
    {
        collisionTime = (l+r)*0.5f;
	Vector3 qi = qd*collisionTime+vertex.p();
	output->push_back(DetectionOutput());
	DetectionOutput &detection = output->back();
	
	detection.point[0] = qi;
	detection.point[1] = qi;    
	detection.normal = -n;
	detection.value = 0;
	detection.elem.first = vertex;
	detection.elem.second = face;
	// WARNING: Does this index need to be unique?
	detection.id = face.getIndex();
    }

    return collisionTime;
}

bool ContinuousIntersection::Impl::testVertexFace(const double &dt, const RVertex &vertex, const RTriangle &face)
{
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

float ContinuousIntersection::Impl::intersectVertexFace(const double &dt, const RVertex &vertex, const RTriangle &face1, const RTriangle &face2, OutputVector *output)
{
  if (!vertex.getBox().Overlaps(face1.getBox()))
    return -1.f;

  typename helper::vector<RVertex::index_type>::const_iterator i, end = vertex.getFaces().end();
  for (i = vertex.getFaces().begin(); i != end; ++i)
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
  double dt;
  helper::vector<NonAdjacentPair>::iterator i = this->nonAdjacentPairs.begin();
  for(helper::vector<NonAdjacentPair>::iterator end = this->nonAdjacentPairs.end(); i != end; ++i)
  {
    RTriangle t1, t2;
    i->getParameters(t1,t2);
    this->testFeatures(dt,t1,t2,output);
  }
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
  this->intersectVertexFace(dt,v1,face2,face1,output);
  this->intersectVertexFace(dt,v2,face2,face1,output);
  this->intersectVertexFace(dt,v3,face2,face1,output);

  this->intersectVertexFace(dt,v4,face1,face2,output);
  this->intersectVertexFace(dt,v5,face1,face2,output);
  this->intersectVertexFace(dt,v6,face1,face2,output);

  // 9 EE test
  REdge e(face1.model,face1.e1Index()), e1(face2.model,face2.e1Index()), e2(face2.model,face2.e2Index()), e3(face2.model,face2.e3Index());
  this->intersectEdgeEdge(dt,e,e1,face1,face2,output);
  this->intersectEdgeEdge(dt,e,e2,face1,face2,output);
  this->intersectEdgeEdge(dt,e,e3,face1,face2,output);

  e = REdge(face1.model,face1.e2Index());
  this->intersectEdgeEdge(dt,e,e1,face1,face2,output);
  this->intersectEdgeEdge(dt,e,e2,face1,face2,output);
  this->intersectEdgeEdge(dt,e,e3,face1,face2,output);

  e = REdge(face1.model,face1.e3Index());
  this->intersectEdgeEdge(dt,e,e1,face1,face2,output);
  this->intersectEdgeEdge(dt,e,e2,face1,face2,output);
  this->intersectEdgeEdge(dt,e,e3,face1,face2,output);
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

template <>
int ContinuousIntersection::computeIntersection(const RTriangle& t1, const RTriangle& t2, OutputVector*)
{
    if(!t1.covertices(t2))
    {
        if(t1.getBox().Overlaps(t2.getBox()))
          this->pimpl->nonAdjacentPairs.push_back(NonAdjacentPair(t1,t2));
    }
    return 0;
}

template <>
bool ContinuousIntersection::testIntersection(Polytope& p1,Polytope& p2)
{
    return p1.overlaps(p2);
}

int ContinuousIntersection::beginIntersection(sofa::core::CollisionModel* model1, sofa::core::CollisionModel* model2, OutputVector* contacts)
{
    return 0;
}

int ContinuousIntersection::endIntersection(sofa::core::CollisionModel* model1, sofa::core::CollisionModel* model2, OutputVector* contacts)
{
    this->pimpl->processNonAdjacentList(dt,contacts);
    this->pimpl->testOrphans(dt,model1,contacts);
    this->pimpl->testOrphans(dt,model2,contacts);
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
template class SOFA_BASE_COLLISION_API IntersectorFactory<component::collision::ContinuousIntersection>;
}
}

} // namespace sofa


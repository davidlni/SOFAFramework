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
//#ifndef SOFA_COMPONENT_COLLISION_TRIANGLEMODEL_INL
//#define SOFA_COMPONENT_COLLISION_TRIANGLEMODEL_INL

#include <vector>
#include <iostream>
#include <cmath>

#include "RTriangleModel.h"
#include "PolytopeModel.h"
#include <sofa/core/visual/VisualParams.h>
#include <sofa/component/collision/TriangleLocalMinDistanceFilter.h>
#include <sofa/component/collision/CubeModel.h>
#include <sofa/component/collision/Triangle.h>
#include <sofa/component/topology/TopologyData.inl>
#include <sofa/simulation/common/Node.h>
#include <sofa/component/topology/RegularGridTopology.h>
#include <sofa/core/CollisionElement.h>
#include <sofa/helper/gl/template.h>

#include <sofa/core/topology/TopologyChange.h>


namespace sofa
{

namespace component
{

namespace collision
{

template<class DataTypes>
TRTriangleModel<DataTypes>::TRTriangleModel()
    : bothSide(initData(&bothSide, false, "bothSide", "activate collision on both side of the triangle model") )
    , mstate(NULL)
    , computeNormals(initData(&computeNormals, true, "computeNormals", "set to false to disable computation of triangles normal"))
    , meshRevision(-1)
    , m_lmdFilter(NULL)
{
    triangles = &mytriangles;
    enum_type = TRIANGLE_TYPE;
}

template<class DataTypes>
void TRTriangleModel<DataTypes>::resize(int size)
{
    this->core::CollisionModel::resize(size);
    normals.resize(size);
}

template<class DataTypes>
void TRTriangleModel<DataTypes>::init()
{
    _topology = this->getContext()->getMeshTopology();

    this->CollisionModel::init();
    mstate = dynamic_cast< core::behavior::MechanicalState<DataTypes>* > (this->getContext()->getMechanicalState());

    this->getContext()->get(mpoints);

    if (mstate==NULL)
    {
        serr << "TriangleModel requires a Vec3 Mechanical Model" << sendl;
        return;
    }

    if (!_topology)
    {
        serr << "TriangleModel requires a BaseMeshTopology" << sendl;
        return;
    }

    simulation::Node* node = dynamic_cast< simulation::Node* >(this->getContext());
    if (node != 0)
    {
        m_lmdFilter = node->getNodeObject< TriangleLocalMinDistanceFilter >();
    }

    //sout << "INFO_print : Col - init TRIANGLE " << sendl;
    sout << "TriangleModel: initially "<<_topology->getNbTriangles()<<" triangles." << sendl;
    triangles = &_topology->getTriangles();
    resize(_topology->getNbTriangles());

    updateFromTopology(this->getContext()->getDt());
    updateNormals();
}

template<class DataTypes>
void TRTriangleModel<DataTypes>::updateNormals()
{
    for (int i=0; i<size; i++)
    {
        Element t(this,i);
        const Vector3& pt1 = t.p1();
        const Vector3& pt2 = t.p2();
        const Vector3& pt3 = t.p3();

        t.n() = cross(pt2-pt1,pt3-pt1);
        t.n().normalize();
        //sout << i << " " << t.n() << sendl;
    }
}

template<class DataTypes>
void TRTriangleModel<DataTypes>::updateFromTopology(const double &dt)
{
    //    needsUpdate = false;
    const unsigned npoints = mstate->getX()->size();
    const unsigned ntris = _topology->getNbTriangles();
    const unsigned nquads = _topology->getNbQuads();
    const unsigned newsize = ntris+2*nquads;

    int revision = _topology->getRevision();
    if (newsize==(unsigned)size && revision == meshRevision)
        return;
    meshRevision = revision;

    needsUpdate=true;

    resize(newsize);

    if (newsize == ntris)
    {
        // no need to copy the triangle indices
        triangles = & _topology->getTriangles();
    }
    else
    {
        helper::vector<EdgeFeature> edges;
        triangles = &mytriangles;
        mytriangles.resize(newsize);
        int index = 0;
        for (unsigned i=0; i<ntris; i++)
        {
            topology::BaseMeshTopology::Triangle idx = _topology->getTriangle(i);
            if (idx[0] >= npoints || idx[1] >= npoints || idx[2] >= npoints)
            {
                serr << "ERROR: Out of range index in triangle "<<i<<": "<<idx[0]<<" "<<idx[1]<<" "<<idx[2]<<" ( total points="<<npoints<<")"<<sendl;
                if (idx[0] >= npoints) idx[0] = npoints-1;
                if (idx[1] >= npoints) idx[1] = npoints-1;
                if (idx[2] >= npoints) idx[2] = npoints-1;
            }
            mytriangles[index] = idx;

            edges.push_back(EdgeFeature(idx[0],idx[1],i));
            edges.push_back(EdgeFeature(idx[1],idx[2],i));
            edges.push_back(EdgeFeature(idx[2],idx[0],i));

            ++index;
        }
        for (unsigned i=0; i<nquads; i++)
        {
            topology::BaseMeshTopology::Quad idx = _topology->getQuad(i);
            if (idx[0] >= npoints || idx[1] >= npoints || idx[2] >= npoints || idx[3] >= npoints)
            {
                serr << "ERROR: Out of range index in quad "<<i<<": "<<idx[0]<<" "<<idx[1]<<" "<<idx[2]<<" "<<idx[3]<<" ( total points="<<npoints<<")"<<sendl;
                if (idx[0] >= npoints) idx[0] = npoints-1;
                if (idx[1] >= npoints) idx[1] = npoints-1;
                if (idx[2] >= npoints) idx[2] = npoints-1;
                if (idx[3] >= npoints) idx[3] = npoints-1;
            }
            mytriangles[index][0] = idx[1];
            mytriangles[index][1] = idx[2];
            mytriangles[index][2] = idx[0];
            ++index;
            mytriangles[index][0] = idx[3];
            mytriangles[index][1] = idx[0];
            mytriangles[index][2] = idx[2];
            ++index;

            edges.push_back(EdgeFeature(idx[1],idx[2],i));
            edges.push_back(EdgeFeature(idx[2],idx[0],i));
            edges.push_back(EdgeFeature(idx[0],idx[1],i));

            edges.push_back(EdgeFeature(idx[3],idx[0],i));
            edges.push_back(EdgeFeature(idx[0],idx[2],i));
            edges.push_back(EdgeFeature(idx[2],idx[3],i));
        }

        std::sort(edges.begin(),edges.end());
        // Eliminate duplicated edges
        for(typename helper::vector<EdgeFeature>::iterator i = edges.begin(), end = edges.end(); i != end; ++i)
        {
            if(!this->edgeFeatures.empty() && *i == edgeFeatures.back())
                this->edgeFeatures.back().SetFaceId(i->FaceId(0));
            else
                this->edgeFeatures.push_back(*i);
        }

        typename helper::vector<EdgeFeature>::iterator begin = this->edgeFeatures.begin();
        typename helper::vector<EdgeFeature>::iterator end = this->edgeFeatures.end();

        this->vertexTriangleFeatures.resize(npoints);
        for (size_t i = 0; i < newsize; ++i)
        {
            Element t(this,i);
            size_t idx0 = t.p1Index();
            size_t idx1 = t.p2Index();
            size_t idx2 = t.p3Index();
            typename helper::vector<EdgeFeature>::iterator mid0 = std::lower_bound(begin,end,EdgeFeature(idx0,idx1));
            typename helper::vector<EdgeFeature>::iterator mid1 = std::lower_bound(begin,end,EdgeFeature(idx1,idx2));
            typename helper::vector<EdgeFeature>::iterator mid2 = std::lower_bound(begin,end,EdgeFeature(idx2,idx0));

            this->triangleEdgeFeatures.push_back(TriangleEdgesFeature());
            this->triangleEdgeFeatures.back().set(std::distance(begin,mid0),
                                                  std::distance(begin,mid1),
                                                  std::distance(begin,mid2));
            // populate vertex-face map
            this->vertexTriangleFeatures[idx0].push_back(i);
            this->vertexTriangleFeatures[idx1].push_back(i);
            this->vertexTriangleFeatures[idx2].push_back(i);
        }

        this->vertexBoxes.resize(npoints);
        this->edgeBoxes.resize(this->edgeFeatures.size());
//         this->faceBoxes.resize(newsize);
        this->updateFeatureBoxes(dt);

    }
    updateFlags();
    updateNormals();
}

template<class DataTypes>
void TRTriangleModel<DataTypes>::updateFlags(int /*ntri*/)
{

}

template<class DataTypes>
void TRTriangleModel<DataTypes>::handleTopologyChange(const double &dt)
{
    //bool debug_mode = false;
    updateFromTopology(dt);
    if (triangles != &mytriangles)
    {
        // We use the same triangle array as the topology -> only resize and recompute flags
        std::list<const sofa::core::topology::TopologyChange *>::const_iterator itBegin=_topology->beginChange();
        std::list<const sofa::core::topology::TopologyChange *>::const_iterator itEnd=_topology->endChange();

        while( itBegin != itEnd )
        {
            core::topology::TopologyChangeType changeType = (*itBegin)->getChangeType();

            switch( changeType )
            {
            case core::topology::ENDING_EVENT:
            {
                sout << "TriangleModel: now "<<_topology->getNbTriangles()<<" triangles." << sendl;
                resize(_topology->getNbTriangles());
                needsUpdate=true;
                updateFlags();

                //                 updateNormals();
                break;
            }

            default:
                break;
            }
            ++itBegin;
        }
        return;
    }
}

template<class DataTypes>
void TRTriangleModel<DataTypes>::draw(const core::visual::VisualParams* vparams ,int index)
{
    Element t(this,index);
    vparams->drawTool()->setPolygonMode(0,vparams->displayFlags().getShowWireFrame());
    vparams->drawTool()->setLightingEnabled(true);
    vparams->drawTool()->drawTriangle( t.p1(), t.p2(), t.p3(), t.n() );
    vparams->drawTool()->setLightingEnabled(false);
}

template<class DataTypes>
void TRTriangleModel<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
    if (vparams->displayFlags().getShowCollisionModels())
    {
        if (bothSide.getValue() || vparams->displayFlags().getShowWireFrame())
            vparams->drawTool()->setPolygonMode(0,vparams->displayFlags().getShowWireFrame());
        else
        {
            vparams->drawTool()->setPolygonMode(2,true);
            vparams->drawTool()->setPolygonMode(1,false);
        }

        std::vector< Vector3 > points;
        std::vector< Vec<3,int> > indices;
        std::vector< Vector3 > normals;
        int index=0;
        for (int i=0; i<size; i++)
        {
            Element t(this,i);
            normals.push_back(t.n());
            points.push_back(t.p1());
            points.push_back(t.p2());
            points.push_back(t.p3());
            indices.push_back(Vec<3,int>(index,index+1,index+2));
            index+=3;
        }

        vparams->drawTool()->setLightingEnabled(true);
        vparams->drawTool()->drawTriangles(points, indices, normals, Vec<4,float>(getColor4f()));
        vparams->drawTool()->setLightingEnabled(false);
        vparams->drawTool()->setPolygonMode(0,false);


        if (vparams->displayFlags().getShowNormals())
        {
            std::vector< Vector3 > points;
            for (int i=0; i<size; i++)
            {
                Element t(this,i);
                points.push_back((t.p1()+t.p2()+t.p3())/3.0);
                points.push_back(points.back()+t.n());
            }

            vparams->drawTool()->drawLines(points, 1, Vec<4,float>(1,1,1,1));

        }
    }
    if (getPrevious()!=NULL && vparams->displayFlags().getShowBoundingCollisionModels())
        getPrevious()->draw(vparams);
}

template<class DataTypes>
bool TRTriangleModel<DataTypes>::canCollideWithElement(size_t index, CollisionModel* model2, size_t index2)
{
    if (!this->bSelfCollision.getValue()) return true; // we need to perform this verification process only for the selfcollision case.
    if (this->getContext() != model2->getContext()) return true;

    Element t(this,index);
    if (model2 == mpoints)
    {
        // if point belong to the triangle, return false
        if ( index2==t.p1Index() || index2==t.p2Index() || index2==t.p3Index())
            return false;
    }

    //// TODO : case with auto-collis with segment and auto-collis with itself

    return true;

}

template<class DataTypes>
void TRTriangleModel<DataTypes>::computeBoundingTree(size_t maxDepth)
{
    PolytopeModel* cubeModel = createPrevious<PolytopeModel>();

    if (needsUpdate && !cubeModel->empty()) cubeModel->resize(0);
    if (!isMoving() && !cubeModel->empty() && !needsUpdate) return; // No need to recompute BBox if immobile

    needsUpdate=false;
    const VecCoord& x = *this->mstate->getX();

    const bool calcNormals = computeNormals.getValue();

    cubeModel->resize(this->size);  // size = number of triangles
    if (!this->empty())
    {
        const SReal distance = this->proximity.getValue();
        for (size_t i = 0; i < this->size; ++i)
        {
            Element t(this,i);
            size_t idx0 = t.p1Index();
            size_t idx1 = t.p2Index();
            size_t idx2 = t.p3Index();

            const Vector3& pt1 = x[idx0];
            const Vector3& pt2 = x[idx1];
            const Vector3& pt3 = x[idx2];

            if (calcNormals)
            {
                // Also recompute normal vector
                t.n() = (pt2-pt1).cross(pt3-pt1);
                t.n().normalize();
            }
            cubeModel->setParentOf(i, pt1); // define the bounding box of the current triangle
            cubeModel->setParentOf(i, pt2);
            cubeModel->setParentOf(i, pt3);
            cubeModel->enlarge(i,distance);

        }
        cubeModel->computeBoundingTree(maxDepth);

    }


    if (m_lmdFilter != 0)
    {
        m_lmdFilter->invalidate();
    }
}

template<class DataTypes>
void TRTriangleModel<DataTypes>::updateFeatureBoxes(const double &dt)
{
    const VecCoord& x = *this->mstate->getX();
    const VecDeriv& v = *this->mstate->getV();
    for(size_t i = 0, end = this->vertexBoxes.size(); i < end; ++i)
    {
        DiscreteOrientedPolytope<Real> box(x[i]+v[i]*dt);
        this->vertexBoxes[i] = box + x[i];
    }

    for(size_t i = 0, end = this->edgeFeatures.size(); i < end; ++i)
    {
        const size_t &id0 = this->edgeFeatures[i].VertexId(0);
        const size_t &id1 = this->edgeFeatures[i].VertexId(1);

        this->edgeBoxes[i] = this->vertexBoxes[id0] + this->vertexBoxes[id1];
    }

    // Get boxes from BVH hierarchy
    core::CollisionModel *boxes = this->getPrevious();
    if(boxes)
      static_cast<PolytopeModel*>(boxes)->updatePolytopes();
//     for(size_t i = 0; i < this->size; ++i)
//     {
//         Element t(this,i);
//         const size_t &id0 = t.p1Index();
//         const size_t &id1 = this->triangleEdgeFeatures[i].EdgeId(1);
//         this->faceBoxes[i] = this->vertexBoxes[id0];
//         this->faceBoxes[i] += this->edgeBoxes[id1];
//     }
}


template<class DataTypes>
void TRTriangleModel<DataTypes>::bufferAdjacentLists()
{
    adjacentPairs[1].clear();
    for(size_t i = 0, end = this->edgeFeatures.size(); i < end; ++i)
    {
        index_pair_type trianglePair(this->edgeFeatures[i].FaceId(0),this->edgeFeatures[i].FaceId(1));
        if (trianglePair.first == -1 || trianglePair.second == -1) continue;
        if (trianglePair.first < trianglePair.second)
        {   //Swap
            index_type tmp = trianglePair.first;
            trianglePair.first = trianglePair.second;
            trianglePair.second = tmp;
        }

        index_pair_type vertexPair;
        index_type cov = this->covertexFace(trianglePair,vertexPair);
        assert(cov == 2);

        char status = this->getStatus2(trianglePair,vertexPair);

        adjacentPairs[1].push_back(AdjacentPair(trianglePair.first, trianglePair.second, vertexPair.first, vertexPair.second, status));
    }

    helper::set<AdjacentPair> adjacentPairSet;
    for (size_t i = 0, end = this->vertexTriangleFeatures.size(); i < end; ++i)
    {
        for (typename VertexTriangleFeatures::iterator it = this->vertexTriangleFeatures[i].begin(), it_end = this->vertexTriangleFeatures[i].end(); it != it_end; ++it)
        {
          for (typename VertexTriangleFeatures::iterator it2 = it; it2 != it_end; ++it2)
            {
                if (it == it2) continue;

                index_pair_type trianglePair(*it,*it2);
                if (trianglePair.first < trianglePair.second)
                {   //Swap
                    index_type tmp = trianglePair.first;
                    trianglePair.first = trianglePair.second;
                    trianglePair.second = tmp;
                }

                index_pair_type vertexPair;
                index_type cov = this->covertexFace(trianglePair,vertexPair);
                if (cov == 2) continue;

                char status = this->getStatus1(trianglePair,vertexPair);

                adjacentPairSet.insert(AdjacentPair(trianglePair.first, trianglePair.second, vertexPair.first, vertexPair.second, status));
            }
        }
    }

    adjacentPairs[0].clear();
//     adjacentPairs[0].resize(adjacentPairSet.size());
    std::copy(adjacentPairSet.begin(),adjacentPairSet.end(),adjacentPairs[0].begin());
}

template<class DataTypes>
char  TRTriangleModel<DataTypes>::getStatus1(const typename TRTriangleModel<DataTypes>::index_pair_type &triangleIndices,
        const typename TRTriangleModel<DataTypes>::index_pair_type &vertexPair)
{
    index_type edge0 = this->triangleEdgeFeatures[triangleIndices.first].EdgeId(vertexPair.first);
    index_type face0 = this->edgeFeatures[edge0].FaceId(0);
    if(face0 == triangleIndices.first)
        face0 = this->edgeFeatures[edge0].FaceId(1);

    index_type edge1 = this->triangleEdgeFeatures[triangleIndices.second].EdgeId(vertexPair.second);
    index_type e11 = this->triangleEdgeFeatures[triangleIndices.second].EdgeId((vertexPair.second+2)%3);

    if (this->edgeFeatures[edge0].Covertex(this->edgeFeatures[edge1]))
    {   // Swap
        index_type tmp = edge1;
        edge1 = e11;
        e11 = tmp;
    }

    index_type face1 = this->edgeFeatures[e11].FaceId(0);
    if(face1 == triangleIndices.second)
        face1 = this->edgeFeatures[e11].FaceId(1);

    if (face0 != -1)
        if (face1 != -1)
            return 3; // just skip it
        else
            return 2; //group 2 only
    else if (face1 != -1)
        return 1; //group 1 only
    else
        return 0; //both group
}

template<class DataTypes>
char TRTriangleModel<DataTypes>::getStatus2(const typename TRTriangleModel<DataTypes>::index_pair_type &triangleIndices,
        const typename TRTriangleModel<DataTypes>::index_pair_type &vertexPair)
{
    index_type edge = this->triangleEdgeFeatures[triangleIndices.first].EdgeId((vertexPair.first+1)%3);
    index_type face0 = this->edgeFeatures[edge].FaceId(0);
    if(face0 == triangleIndices.first)
        face0 = this->edgeFeatures[edge].FaceId(1);

    edge = this->triangleEdgeFeatures[triangleIndices.second].EdgeId((vertexPair.second+1)%3);
    index_type face1 = this->edgeFeatures[edge].FaceId(0);
    if(face1 == triangleIndices.second)
        face1 = this->edgeFeatures[edge].FaceId(1);

    if (face0 != -1)
        if (face1 != -1)
            return 3; //just skip it ...
        else
            return 2; //group 2 only
    else if (face1 != -1)
        return 1; //group 1 only
    else
        return 0; //both group
}

template<class DataTypes>
typename TRTriangleModel<DataTypes>::index_type
TRTriangleModel<DataTypes>::covertexFace(const typename TRTriangleModel<DataTypes>::index_pair_type &triangleIndices,
        typename TRTriangleModel<DataTypes>::index_pair_type &vertexPair)
{
    Element face[2] = { Element(this,triangleIndices.first),
                        Element(this,triangleIndices.second)
                      };
    return face[0].covertices(face[1],vertexPair);
}

template<class DataTypes>
void TRTriangleModel<DataTypes>::setOrphans()
{
    std::pair<index_type,index_type> trianglePair, vertexPair;

    for (typename helper::vector<AdjacentPair>::iterator i=adjacentPairs[0].begin(); i!=adjacentPairs[0].end(); ++i)
    {
        i->GetParameters(trianglePair,vertexPair);
        this->getFeature1(trianglePair,vertexPair);
    }

    for (typename helper::vector<AdjacentPair>::iterator i=adjacentPairs[1].begin(); i!=adjacentPairs[1].end(); ++i)
    {
        i->GetParameters(trianglePair,vertexPair);
        this->getFeature1(trianglePair,vertexPair);
    }

}

// template<class DataTypes>
// void TRTriangleModel<DataTypes>::testOrphans(const double &dt)
// {
//     for (typename helper::set<FeaturePair>::iterator i=this->edgeEdgeFeaturePairs.begin(), end=this->edgeEdgeFeaturePairs.end(); i!=end; ++i)
//     {
//         index_type e1, e2;
//         i->GetParameters(e1,e2);
//         this->intersectEdgeEdge(dt,EdgeElement(this,e1),EdgeElement(this,e2));
//     }
//     for (typename helper::set<FeaturePair>::iterator i=this->vertexFaceFeaturePairs.begin(), end=this->vertexFaceFeaturePairs.end(); i!=end; ++i)
//     {
//         index_type v, f;
//         i->GetParameters(f, v);
//         this->intersectVertexFace(dt,VertexElement(this,v), Element(this,f));
//     }
// }

template<class DataTypes>
bool TRTriangleModel<DataTypes>::testOrphansEdgeToEdge(const index_type &i1, const index_type &i2)
{
    for (index_type i=0; i<2; i++)
    {
        index_type f1 = this->edgeFeatures[i1].FaceId(i);
        if (f1 == -1) continue;

        for (index_type j=0; j<2; j++)
        {
            index_type f2 = this->edgeFeatures[i2].FaceId(j);
            if (f2 == -1) continue;

            if (!Element(this,f1).covertices(Element(this,f2)))
                return false;
        }
    }

    return true;
}

template<class DataTypes>
bool TRTriangleModel<DataTypes>::testOrphansVertexToFace(const index_type &v, const index_type &f)
{
    for (typename VertexTriangleFeatures::iterator i=vertexTriangleFeatures[v].begin(),
      end = vertexTriangleFeatures[v].end(); i != end; i++)
    {
        index_type f2 = *i;
        if (!Element(this,f2).covertices(Element(this,f)))
            return false;
    }
    return true;
}

template<class DataTypes>
void TRTriangleModel<DataTypes>::getFeature1(const index_pair_type &triangleIndices, const index_pair_type &indexPair)
{
    // 2 VF test
    Element t = Element(this,triangleIndices.second);
    for (int i=0; i<3; i++)
    {
        if (i == indexPair.second) continue; // skip one
        this->insertVertexToFaceFeature(triangleIndices.first, t.Index(i));
    }

    // 2 EE test
    index_type e1 = this->triangleEdgeFeatures[triangleIndices.second].EdgeId((indexPair.second+1)%3);
    index_type e0 = this->triangleEdgeFeatures[triangleIndices.first].EdgeId((indexPair.first)%3);
    this->insertEdgeToEdgeFeature(e0, e1);
    e0 = this->triangleEdgeFeatures[triangleIndices.first].EdgeId((indexPair.first+2)%3);
    this->insertEdgeToEdgeFeature(e0, e1);

    // 2 VF test
    t = Element(this,triangleIndices.first);
    for (int i=0; i<3; i++)
    {
        if (i == indexPair.first) continue; // skip one
        this->insertVertexToFaceFeature(triangleIndices.second, t.Index(i));
    }

    // 3 EE test
    e0 = this->triangleEdgeFeatures[triangleIndices.first].EdgeId((indexPair.first+1)%3);
    for (int i=0; i<3; i++)
    {
      index_type e1 = this->triangleEdgeFeatures[triangleIndices.second].EdgeId(i);
        this->insertEdgeToEdgeFeature(e0, e1);
    }
}

template<class DataTypes>
void TRTriangleModel<DataTypes>::getFeature2(const index_pair_type &triangleIndices, const index_pair_type &indexPair)
{
    index_type e0 = this->triangleEdgeFeatures[triangleIndices.first].EdgeId(indexPair.first);
    index_type e00 = this->triangleEdgeFeatures[triangleIndices.first].EdgeId((indexPair.first+2)%3);
    index_type e1 = this->triangleEdgeFeatures[triangleIndices.second].EdgeId(indexPair.second);
    index_type e11 = this->triangleEdgeFeatures[triangleIndices.second].EdgeId((indexPair.second+2)%3);

    if (this->edgeFeatures[e0].Covertex(this->edgeFeatures[e1]))
    {   // Swap
        index_type tmp = e1;
        e1 = e11;
        e11 = tmp;
    }

    Element face = Element(this,triangleIndices.first);
    index_type fid = triangleIndices.second;
    index_type vid = face.Index(indexPair.first);
    this->insertVertexToFaceFeature(fid, vid);
    this->insertEdgeToEdgeFeature(e0, e1);

    // 2 VF test
    face = Element(this,triangleIndices.second);
    fid = triangleIndices.first;
    vid = face.Index(indexPair.second);
    this->insertVertexToFaceFeature(fid, vid);
    this->insertEdgeToEdgeFeature(e00, e11);

}

template<class DataTypes>
void TRTriangleModel<DataTypes>::insertEdgeToEdgeFeature(const index_type &e1, const index_type &e2)
{
    if (this->testOrphansEdgeToEdge(e1, e2))
    {
        this->edgeEdgeFeaturePairs.insert(FeaturePair(e1, e2));
    }
}

template<class DataTypes>
void TRTriangleModel<DataTypes>::insertVertexToFaceFeature(const index_type &f, const index_type &v)
{
    if (this->testOrphansVertexToFace(f, v))
    {
        this->vertexFaceFeaturePairs.insert(FeaturePair(f, v));
    }
}
// 
// template<class DataTypes>
// float TRTriangleModel<DataTypes>::intersectVertexFace(const double &dt, const VertexElement &vertex, const Element &f1, const Element &f2)
// {
//   if (!Polytope(static_cast<PolytopeModel*>(this->getPrevious()),f1.getIndex()).getElement().Overlaps(this->vertexBoxes[vertex.getIndex()]))
//         return -1.f;
//     typename VertexTriangleFeatures::const_iterator i, end = this->vertexTriangleFeatures[vertex.getIndex()].end();
//     for (i = this->vertexTriangleFeatures[vertex.getIndex()].begin(); i != end; ++i)
//     {
//         index_type fid = *i;
//         if (!f1.covertices(Element(this,fid)))
//         {
//             if (fid == f2.getIndex())
//             {
//                 if (this->testVertexFace(dt,vertex,f1))
//                     return this->vertexFaceIntersection(dt,vertex,f1);
//                 return -1.f;
//             }
//             return -1.f;
//         }
//     }
//     return -1.f;
// }

// template<class DataTypes>
// float TRTriangleModel<DataTypes>::intersectVertexFace(const double &dt, const VertexElement &vertex, const Element &face)
// {
//     if (!Polytope(static_cast<PolytopeModel*>(this->getPrevious()),face.getIndex()).getElement().Overlaps(this->vertexBoxes[vertex.getIndex()]))
//         return false;
// 
//     if (this->testVertexFace(dt,vertex,face))
//         return this->vertexFaceIntersection(dt,vertex,face);
// 
//     return -1.f;
// }

// template<class DataTypes>
// float TRTriangleModel<DataTypes>::intersectEdgeEdge(const double &dt, const EdgeElement &e1, const EdgeElement &e2)
// {
//     if (!e1.getBox().Overlaps(e2.getBox()))
//         return -1.f;
// 
//     if (this->testEdgeEdge(dt, e1, e2))
//         return this->edgeEdgeIntersection(dt,e1,e2);
// 
//     return -1.f;
// }

// template<class DataTypes>
// float TRTriangleModel<DataTypes>::intersectEdgeEdge(const double &dt, const EdgeElement &e1, const EdgeElement &e2, const Element &f1, const Element &f2)
// {
//   if (!this->edgeBoxes[e1.getIndex()].Overlaps(this->edgeBoxes[e2.getIndex()]))
//         return -1.f;
// 
//     EdgeElement e[2] = {e1,e2};
//     Element f[2] = {f1,f2};
// 
//     if (e1.getIndex() > e2.getIndex())
//     {
//         e[0] = e2, e[1] = e1;
//         f[0] = f2, f[1] = f1;
//     }
// 
//     for (int i=0; i<2; i++)
//         for (int j=0; j<2; j++)
//         {
//             index_type ff1 = e[0].FaceId(i);
//             index_type ff2 = e[1].FaceId(j);
// 
//             if (ff1 == -1 || ff2 == -1)
//                 continue;
// 
//             Element adjacentFace1(this,ff1);
//             Element adjacentFace2(this,ff2);
//             if (!adjacentFace1.covertices(adjacentFace2))
//             {
//               if (adjacentFace1 == f[0] && adjacentFace2 == f[1])
//                 {
//                     if (this->testEdgeEdge(dt, e1, e2))
//                         return this->edgeEdgeIntersection(dt,e1,e2);
//                     return -1.f;
//                 }
//                 return -1.f;
//             }
//         }
// 
//     return -1.f;
// }
// 
// template<class DataTypes>
// float TRTriangleModel<DataTypes>::vertexFaceIntersection(const double &dt, const VertexElement &vertex, const Element &face)
// {
//     /* Default value returned if no collision occurs */
//     float collisionTime = -1.0f;
// 
//     /* diff. vectors for linear interpolation: x1-x0 == dt*v0 */
//     Vector3 qd = vertex.v()*dt;
//     Vector3 ad = face.v1()*dt;
//     Vector3 bd = face.v2()*dt;
//     Vector3 cd = face.v3()*dt;
// 
//     /*
//      * Compute scalar coefficients by evaluating dot and cross-products.
//      */
//     helper::fixed_array<Real,4> coeffs; /* cubic polynomial coefficients */
//     this->computeCubicCoefficientsVertexFace(face.p1(), ad, face.p2(), bd, face.p3(), cd, vertex.p(), qd, coeffs);
// 
//     if (std::fabs(coeffs[0]) <= std::numeric_limits<Real>::epsilon() &&
//       std::fabs(coeffs[1]) <= std::numeric_limits<Real>::epsilon() &&
//       std::fabs(coeffs[2]) <= std::numeric_limits<Real>::epsilon() &&
//       std::fabs(coeffs[3]) <= std::numeric_limits<Real>::epsilon())
//         return collisionTime;
// 
//     NewtonCheckData data;
//     data.a0 = face.p1(), data.b0 = face.p2();
//     data.c0 = face.p3(), data.p0 = vertex.p();
//     data.ad = ad, data.bd = bd;
//     data.cd = cd, data.pd = qd;
// 
//     /*
//      * iteratively solve the cubic (scalar) equation and test for validity of the solution.
//      */
//     Real l = 0;
//     Real r = 1;
// 
//     if (this->solveCubicWithIntervalNewton(l, r, true, data, coeffs))
//     {
//         collisionTime = (l+r)*0.5f;
//     }
// 
//     return collisionTime;
// }
// 
// template<class DataTypes>
// float TRTriangleModel<DataTypes>::edgeEdgeIntersection(const double &dt, const EdgeElement &edge1, const EdgeElement &edge2)
// {
//     /* Default value returned if no collision occurs */
//     float collisionTime = -1.0f;
// 
//     /* diff. vectors for linear interpolation: x1-x0 == dt*v0 */
//     Vector3 dd = edge2.v2()*dt;
//     Vector3 ad = edge1.v1()*dt;
//     Vector3 bd = edge1.v2()*dt;
//     Vector3 cd = edge2.v1()*dt;
// 
//     /*
//     * Compute scalar coefficients by evaluating dot and cross-products.
//     */
//     helper::fixed_array<Real,4> coeffs; /* cubic polynomial coefficients */
//     this->computeCubicCoefficientsEdgeEdge(edge1.p1(), ad, edge1.p2(), bd, edge2.p1(), cd, edge2.p2(), dd, coeffs);
// 
//     if (std::fabs(coeffs[0]) <= std::numeric_limits<Real>::epsilon() &&
//       std::fabs(coeffs[1]) <= std::numeric_limits<Real>::epsilon() &&
//       std::fabs(coeffs[2]) <= std::numeric_limits<Real>::epsilon() &&
//       std::fabs(coeffs[3]) <= std::numeric_limits<Real>::epsilon())
//         return collisionTime;
// 
//     NewtonCheckData data;
//     data.a0 = edge1.p1(), data.b0 = edge1.p2();
//     data.c0 = edge2.p1(), data.p0 = edge2.p2();
//     data.ad = ad, data.bd = bd;
//     data.cd = cd, data.pd = dd;
// 
//     /*
//     * iteratively solve the cubic (scalar) equation and test for validity of the solution.
//     */
//     Real l = 0;
//     Real r = 1;
// 
//     if (solveCubicWithIntervalNewton(l, r, false, data, coeffs))
//     {
//         collisionTime = (l+r)*0.5f;
//     }
// 
//     return collisionTime;
// }
// 
// template<class DataTypes>
// void TRTriangleModel<DataTypes>::computeCubicCoefficientsVertexFace(
//     const Vector3 &a0, const Vector3 &ad, const Vector3 &b0, const Vector3 &bd,
//     const Vector3 &c0, const Vector3 &cd, const Vector3 &p0, const Vector3 &pd,
//     helper::fixed_array<Real,4> &coeff)
// {
//     this->computeCubicCoefficients(bd-ad,cd-ad,pd-ad,c0-a0,b0-a0,p0-a0,coeff);
// }
// 
// template<class DataTypes>
// void TRTriangleModel<DataTypes>::computeCubicCoefficientsEdgeEdge(
//     const Vector3 &a0, const Vector3 &ad, const Vector3 &b0, const Vector3 &bd,
//     const Vector3 &c0, const Vector3 &cd, const Vector3 &d0, const Vector3 &dd,
//     helper::fixed_array<Real,4> &coeff)
// {
//     this->computeCubicCoefficients(bd-ad,dd-cd,cd-ad,d0-c0,b0-a0,c0-a0,coeff);
// }
// 
// template<class DataTypes>
// void TRTriangleModel<DataTypes>::computeCubicCoefficients(
//     const Vector3 &a, const Vector3 &b, const Vector3 &c, const Vector3 &d,
//     const Vector3 &e, const Vector3 &f, helper::fixed_array<Real,4> &coeff)
// {
//     Vector3 v1 = a.cross(b);
//     Vector3 v2 = a.cross(d);
//     Vector3 v3 = e.cross(b);
//     Vector3 v4 = e.cross(d);
// 
//     coeff[0] = f*v4;
//     coeff[1] = c*v4 + f*(v2 + v3);
//     coeff[2] = f*v1 + c*(v2 + v3);
//     coeff[3] = c*v1;
// }
// 
// template<class DataTypes>
// bool TRTriangleModel<DataTypes>::solveCubicWithIntervalNewton(Real &l, Real &r, bool bVF, NewtonCheckData &data, helper::fixed_array<Real,4> &coeff)
// {
//     Real v2[2] = {l*l,r*r};
//     Real v[2] = {l,r};
//     Real rBkUp;
// 
//     unsigned char min3, min2, min1, max3, max2, max1;
// 
//     min3=*((unsigned char*)&coeff[3]+7) >> 7;
//     max3=min3^1;
//     min2=*((unsigned char*)&coeff[2]+7) >> 7;
//     max2=min2^1;
//     min1=*((unsigned char*)&coeff[1]+7) >> 7;
//     max1=min1^1;
// 
//     // Compute bounds for the cubic
//     Real bounds[2] = { coeff[3]*v2[min3]*v[min3]+coeff[2]*v2[min2]+coeff[1]*v[min1]+coeff[0],
//                        coeff[3]*v2[max3]*v[max3]+coeff[2]*v2[max2]+coeff[1]*v[max1]+coeff[0]
//                      };
// 
//     if (bounds[0]<0) return false;
//     if (bounds[1]>0) return false;
// 
//     // Starting starts at the middle of the interval
//     Real m = 0.5*(r+l);
// 
//     // Compute bounds for the derivative
//     Real dbounds[2] = { 3.0*coeff[3]*v2[min3]+2.0*coeff[2]*v[min2]+coeff[1],
//                         3.0*coeff[3]*v2[max3]+2.0*coeff[2]*v[max2]+coeff[1]
//                       };
// 
//     if ((dbounds[0] > 0) || (dbounds[1] < 0)) // we can use Newton
//     {
//         Real m2 = m*m;
//         Real fm = coeff[3]*m2*m+coeff[2]*m2+coeff[1]*m+coeff[0];
//         Real nl = m;
//         Real nu = m;
//         if (fm > 0)
//         {
//             nl -= fm*(1.0/dbounds[0]);
//             nu -= fm*(1.0/dbounds[1]);
//         }
//         else
//         {
//             nu -= fm*(1.0/dbounds[0]);
//             nl -= fm*(1.0/dbounds[1]);
//         }
// 
//         // Intersect with [l,r]
//         if (nu<l || nl>r)
//           return false; // outside the interval
// 
//         if (nl>l)
//         {
//             if (nu<r)
//             {
//                 l = nl;
//                 r = nu;
//                 m = 0.5*(l+r);
//             }
//             else
//             {
//                 l = nl;
//                 m = 0.5*(l+r);
//             }
//         }
//         else
//         {
//             if (nu<r)
//             {
//                 r = nu;
//                 m = 0.5*(l+r);
//             }
//         }
//     }
// 
//     // Check root validity
//     if ((r-l)< Real(1e-10))
//     {
//         if (bVF)
//             return this->insideTriangle(data.ad*r + data.a0, data.bd*r + data.b0, data.cd*r + data.c0, data.pd*r + data.p0);
//         else
//             return this->lineLineIntersect(data.ad*r + data.a0, data.bd*r + data.b0, data.cd*r + data.c0, data.pd*r + data.p0);
//     }
//     rBkUp = r, r = m;
//     if (this->solveCubicWithIntervalNewton(l,r, bVF, data, coeff)) return true;
//     l = m, r = rBkUp;
//     return (this->solveCubicWithIntervalNewton(l,r, bVF, data, coeff));
// }
// 
// /*
//  * Ordinary inside-triangle test for p. The triangle normal is computed from the vertices.
//  */
// template<class DataTypes>
// bool TRTriangleModel<DataTypes>::insideTriangle(const Vector3 &a, const Vector3 &b, const Vector3 &c, const Vector3 &p)
// {
//     Vector3 n, da, db, dc;
//     Real wa, wb, wc;
// 
//     Vector3 ba = b-a;
//     Vector3 ca = c-a;
//     n = ba.cross(ca);
// 
//     da = a - p, db = b - p, dc = c - p;
//     if ((wa = db.cross(dc)*n) < 0.0f) return false;
//     if ((wb = dc.cross(da)*n) < 0.0f) return false;
//     if ((wc = da.cross(db)*n) < 0.0f) return false;
// 
//     //Compute barycentric coordinates
//     Real area2 = n.norm2();
//     wa /= area2, wb /= area2, wc /= area2;
// 
//     // Does not returns it
//     Vector3 baryc = Vector3(wa, wb, wc);
// 
//     return true;
// }
// 
// template<class DataTypes>
// bool TRTriangleModel<DataTypes>::insideLineSegment(const Vector3 &a, const Vector3 &b, const Vector3 &c)
// {
//     return (a-b)*(a-c)<=0;
// }
// 
// template<class DataTypes>
// bool TRTriangleModel<DataTypes>::testFeatures(const double &dt, const Element &face1, const Element &face2)
// {
//     // 6 VF test
//     this->intersectVertexFace(dt,VertexElement(this,face1.p1Index()), face2, face1);
//     this->intersectVertexFace(dt,VertexElement(this,face1.p2Index()), face2, face1);
//     this->intersectVertexFace(dt,VertexElement(this,face1.p3Index()), face2, face1);
// 
//     this->intersectVertexFace(dt,VertexElement(this,face2.p1Index()), face1, face2);
//     this->intersectVertexFace(dt,VertexElement(this,face2.p2Index()), face1, face2);
//     this->intersectVertexFace(dt,VertexElement(this,face2.p3Index()), face1, face2);
// 
//     // 9 EE test
//     EdgeElement e(this,face1.e1Index()), e1(this,face2.e1Index()), e2(this,face2.e2Index()), e3(this,face2.e3Index());
//     this->intersectEdgeEdge(dt,e, e1, face1, face2);
//     this->intersectEdgeEdge(dt,e, e2, face1, face2);
//     this->intersectEdgeEdge(dt,e, e3, face1, face2);
// 
//     e = EdgeElement(this,face1.e2Index());
//     this->intersectEdgeEdge(dt,e, e1, face1, face2);
//     this->intersectEdgeEdge(dt,e, e2, face1, face2);
//     this->intersectEdgeEdge(dt,e, e3, face1, face2);
// 
//     e = EdgeElement(this,face1.e3Index());
//     this->intersectEdgeEdge(dt,e, e1, face1, face2);
//     this->intersectEdgeEdge(dt,e, e2, face1, face2);
//     this->intersectEdgeEdge(dt,e, e3, face1, face2);
//     return true;
// }
// 
// /*
//  Calculate the line segment PaPb that is the shortest route between
//  two lines P1P2 and P3P4. Calculate also the values of mua and mub where
//  Pa = P1 + mua (P2 - P1)
//  Pb = P3 + mub (P4 - P3)
//  Return FALSE if no solution exists.
//  */
// template<class DataTypes>
// bool TRTriangleModel<DataTypes>::lineLineIntersect(
//     const Vector3 &p1, const Vector3 &p2, const Vector3 &p3, const Vector3 &p4)
// {
//     Vector3 p13,p43,p21;
//     Real d1343,d4321,d1321,d4343,d2121;
//     Real numer,denom;
// 
//     p13 = p1 - p3;
//     p43 = p4 - p3;
//     if (fabs(p43[0])  < std::numeric_limits<Real>::epsilon() &&
//             fabs(p43[1])  < std::numeric_limits<Real>::epsilon() &&
//             fabs(p43[2])  < std::numeric_limits<Real>::epsilon())
//         return false;
// 
//     p21 = p2 - p1;
//     if (fabs(p21[0])  < std::numeric_limits<Real>::epsilon() &&
//             fabs(p21[1])  < std::numeric_limits<Real>::epsilon() &&
//             fabs(p21[2])  < std::numeric_limits<Real>::epsilon())
//         return false;
// 
//     d1343 = p13*p43;
//     d4321 = p43*p21;
//     d1321 = p13*p21;
//     d4343 = p43*p43;
//     d2121 = p21*p21;
// 
//     denom = d2121 * d4343 - d4321 * d4321;
//     if (fabs(denom) < std::numeric_limits<Real>::epsilon())
//         return false;
//     numer = d1343 * d4321 - d1321 * d4343;
// 
//     Real mua = numer / denom;
//     if (mua < 0 || mua > 1)
//         return false;
// 
//     Real mub = (d1343 + d4321 * mua) / d4343;
//     if (mub < 0 || mub > 1)
//         return false;
// 
//     // Do not return these
//     Vector3 pa, pb;
//     pa = p1 + p21*mua;
//     pb = p3 + p43*mub;
//     return true;
// }
// 
// template<class DataTypes>
// bool TRTriangleModel<DataTypes>::testEdgeEdge(const double &dt, const EdgeElement &e1, const EdgeElement &e2)
// {
//     VertexElement pe1(this,e1.p1Index());
//     VertexElement pe2(this,e1.p2Index());
//     VertexElement pe3(this,e2.p1Index());
//     VertexElement pe4(this,e2.p2Index());
// 
//     return this->testCoplanarity(pe1.p(), pe2.p(), pe3.p(), pe4.p(),
//                                  pe1.p()+pe1.v()*dt, pe2.p()+pe2.v()*dt, pe3.p()+pe3.v()*dt, pe4.p()+pe4.v()*dt);
// }
// 
// template<class DataTypes>
// bool TRTriangleModel<DataTypes>::testVertexFace(const double &dt, const VertexElement &vertex, const Element &face)
// {
//     return this->testCoplanarity(face.p1(), face.p2(), face.p3(), vertex.p(),face.p1()+face.v1()*dt,
//                                  face.p2()+face.v2()*dt, face.p3()+face.v3()*dt,vertex.p()+vertex.v()*dt);
// }
// 
// template<class DataTypes>
// bool TRTriangleModel<DataTypes>::testCoplanarity(const Vector3 &p0, const Vector3 &p1, const Vector3 &p2, const Vector3 &p3,
//         const Vector3 &p4, const Vector3 &p5, const Vector3 &p6, const Vector3 &p7
//                                                 )
// {
//     Vector3 n0 = (p1-p0).cross(p2-p1);
//     Vector3 n1 = (p5-p4).cross(p6-p5);
//     Vector3 delta = ((p5-p1)-(p4-p0)).cross((p6-p2)-(p5-p1));
//     Vector3 nX = (n0+n1-delta)*0.5;
// 
//     Vector3 pa0 = p3-p0;
//     Vector3 pa1 = p7-p4;
// 
//     Real A = n0*pa0;
//     Real B = n1*pa1;
//     Real C = nX*pa0;
//     Real D = nX*pa1;
//     Real E = n1*pa0;
//     Real F = n0*pa1;
// 
//     if (A > 0 && B > 0 && (2*C+F) > 0 && (2*D+E) > 0)
//         return false;
// 
//     if (A < 0 && B < 0 && (2*C+F) < 0 && (2*D+E) < 0)
//         return false;
// 
//     return true;
// }

template<class DataTypes>
void TRTriangleModel<DataTypes>::computeContinuousBoundingTree(double dt, size_t maxDepth)
{
    PolytopeModel* cubeModel = this->createPrevious<PolytopeModel>();

    if (this->needsUpdate) cubeModel->resize(0);
    if (!this->isMoving() && !cubeModel->empty() && !this->needsUpdate) return; // No need to recompute BBox if immobile

    this->needsUpdate = false;
    const VecCoord& x = *this->mstate->getX();
    const VecDeriv& v = *this->mstate->getV();

    const bool calcNormals = this->computeNormals.getValue();
    cubeModel->resize(this->size);
    if (!this->empty())
    {
        const SReal distance = this->proximity.getValue();
        for (int i = 0; i < this->size; i++)
        {
            Element t(this,i);
            const Vector3& pt1 = x[t.p1Index()];
            const Vector3& pt2 = x[t.p2Index()];
            const Vector3& pt3 = x[t.p3Index()];

            if (calcNormals)
            {
                // Also recompute normal vector
                t.n() = (pt2-pt1).cross(pt3-pt1);
                t.n().normalize();
            }
            cubeModel->setParentOf(i,pt1);
            cubeModel->setParentOf(i,pt2);
            cubeModel->setParentOf(i,pt3);

            cubeModel->setParentOf(i,pt1 + v[t.p1Index()]*dt);
            cubeModel->setParentOf(i,pt2 + v[t.p2Index()]*dt);
            cubeModel->setParentOf(i,pt3 + v[t.p3Index()]*dt);

            cubeModel->enlarge(i,distance);
        }
        cubeModel->computeBoundingTree(maxDepth);
    }
}

template<class DataTypes>
TriangleLocalMinDistanceFilter *TRTriangleModel<DataTypes>::getFilter() const
{
    return m_lmdFilter;
}

template<class DataTypes>
void TRTriangleModel<DataTypes>::setFilter(TriangleLocalMinDistanceFilter *lmdFilter)
{
    m_lmdFilter = lmdFilter;
}

template<class DataTypes>
int TRTriangleModel<DataTypes>::getTriangleFlags(int i)
{
    int f = 0;
    sofa::core::topology::BaseMeshTopology::Triangle t = (*triangles)[i];

    if (i < _topology->getNbTriangles())
    {
        for (unsigned int j=0; j<3; ++j)
        {
            const sofa::core::topology::BaseMeshTopology::TrianglesAroundVertex& tav = _topology->getTrianglesAroundVertex(t[j]);
            if (tav[0] == (sofa::core::topology::BaseMeshTopology::TriangleID)i)
                f |= (FLAG_P1 << j);
        }

        const sofa::core::topology::BaseMeshTopology::EdgesInTriangle& e = _topology->getEdgesInTriangle(i);

        for (unsigned int j=0; j<3; ++j)
        {
            const sofa::core::topology::BaseMeshTopology::TrianglesAroundEdge& tae = _topology->getTrianglesAroundEdge(e[j]);
            if (tae[0] == (sofa::core::topology::BaseMeshTopology::TriangleID)i)
                f |= (FLAG_E23 << j);
            if (tae.size() == 1)
                f |= (FLAG_BE23 << j);
        }
    }
    else
    {
        /// \TODO flags for quads
    }
    return f;
}

} // namespace collision

} // namespace component

} // namespace sofa

//#endif

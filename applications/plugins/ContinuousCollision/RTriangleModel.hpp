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
#include <vector>
#include <sofa/helper/gl/template.h>
#include <iostream>

#include <sofa/core/topology/TopologyChange.h>

#include <sofa/simulation/common/Simulation.h>

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

    updateFromTopology();
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
void TRTriangleModel<DataTypes>::updateFromTopology()
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
        for(helper::vector<EdgeFeature>::iterator i = edges.begin(), end = edges.end(); i != end; ++i)
        {
            if(!this->edgeFeatures.empty() && *i == edgeFeatures.back())
                this->edgeFeatures.back().SetFaceId(i->FaceId(0));
            else
                this->edgeFeatures.push_back(*i);
        }

        helper::vector<EdgeFeature>::iterator begin = this->edgeFeatures.begin();
        helper::vector<EdgeFeature>::iterator end = this->edgeFeatures.end();

        this->vertexTriangleFeatures.resize(npoints);
        for (size_t i = 0; i < newsize; ++i)
        {
            Element t(this,i);
            size_t idx0 = t.p1Index();
            size_t idx1 = t.p2Index();
            size_t idx2 = t.p3Index();
            helper::vector<EdgeFeature>::iterator mid0 = std::lower_bound(first,last,EdgeFeature(idx0,idx1));
            helper::vector<EdgeFeature>::iterator mid1 = std::lower_bound(first,last,EdgeFeature(idx1,idx2));
            helper::vector<EdgeFeature>::iterator mid2 = std::lower_bound(first,last,EdgeFeature(idx2,idx0));

            this->triangleEdgeFeatures.push_back(TriangleEdgesFeature());
            this->triangleEdgeFeatures.back().set(std::distance(first,mid0),
                                                  std::distance(first,mid1),
                                                  std::distance(first,mid2));
            // populate vertex-face map
            this->vertexTriangleFeatures[idx0].push_back(i);
            this->vertexTriangleFeatures[idx1].push_back(i);
            this->vertexTriangleFeatures[idx2].push_back(i);
        }

        this->vertexBoxes.resize(npoints);
        this->edgeBoxes.resize(this->edgeFeatures.size());
        this->faceBoxes.resize(newsize);
        this->updateFeatureBoxes();

    }
    updateFlags();
    updateNormals();
}

template<class DataTypes>
void TRTriangleModel<DataTypes>::updateFlags(int /*ntri*/)
{

}

template<class DataTypes>
void TRTriangleModel<DataTypes>::handleTopologyChange()
{
    //bool debug_mode = false;
    updateFromTopology();
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
        const SReal distance = (SReal)this->proximity.getValue();
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
                t.n() = (pt2-pt1).cross(pt3-pt1).normalize();
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
void TRTriangleModel<DataTypes>::updateFeatureBoxes(double dt)
{
    const VecCoord& x = *this->mstate->getX();
    const VecDeriv& v = *this->mstate->getV();

    for(size_t i = 0, end = this->vertexBoxes.size(); i < end; ++i)
    {
        DiscreteOrientedPolytope<Real> box(x[i]+dt*v[i]);
        this->vertexBoxes[i] = box + x[i];
    }

    for(size_t i = 0, end = this->edgeFeatures.size(); i < end; ++i)
    {
        const size_t &id0 = this->edgeFeatures[i].VertexId(0);
        const size_t &id1 = this->edgeFeatures[i].VertexId(1);

        this->edgeBoxes[i] = this->vertexBoxes[id0] + this->vertexBoxes[id1];
    }

    for(size_t i = 0; i < this->size; ++i)
    {
        Element t(this,i);
        size_t id0 = t.p1Index();
        size_t id1 = this->triangleEdgeFeatures[i].VertexId(1);
        this->faceBoxes[i] = this->vertexBoxes[id0] + edgeBoxes[id1];
    }
}


template<class DataTypes>
void TRTriangleModel<DataTypes>::bufferAdjacentLists()
{
    adjacentPairList[1].clear();
    for(size_t i = 0, end = this->edgeFeatures.size(); i < end; ++i)
    {
        std::pair<size_t,size_t> trianglePair(this->edgeFeatures[i].FaceId(0),this->edgeFeatures[i].FaceId(1));
        if (trianglePair.first == -1 || trianglePair.second == -1) continue;
        if (trianglePair.first < trianglePair.second)
        {   //Swap
            size_t tmp = trianglePair.first;
            trianglePair.first = trianglePair.second;
            trianglePair.second = tmp;
        }

        std::pair<size_t,size_t> vertexPair;
        size_t cov = this->covertexFace(trianglePair,vertexPair);
        assert(cov == 2);

        char status = this->getStatus2(trianglePair,vertexPair);
        //if (status == 3) continue;

        adjacentPairList[1].push_back(AdjacentPairs(trianglePair.first, trianglePair.second, vertexPair.first, vertexPair.second, status));
    }

    adjacentPairList[2].clear();
    for (size_t i=0; i<_num_vtx; i++) {
      for (id_list::iterator it1=_vtx_fids[i].begin(); it1!=_vtx_fids[i].end(); it1++) {
        for (id_list::iterator it2=it1; it2!=_vtx_fids[i].end(); it2++) {
          if (it1 == it2) continue;

          size_t id1 = *it1;
          size_t id2 = *it2;
          if (id1 < id2) swapI(id1, id2);

          size_t st1 = 0, st2 = 0;
          size_t cov = Covertex_F(id1, id2, st1, st2);
          if (cov == 2) continue;

          char status = get_status_1(id1, id2, st1, st2);
          //if (status == 3) continue;

          adj_1_list.push_back(adjacent_pair(id1, id2, st1, st2, status));
        }
      }
    }
}

template<class DataTypes>
char  TRTriangleModel<DataTypes>::getStatus1(const std::pair<size_t,size_t> &triangleIndices,
        std::pair<size_t,size_t> &vertexPair)
{
    size_t edge0 = this->triangleEdgeFeatures[triangleIndices.first].EdgeId(vertexPair.first);
    size_t face0 = this->edgeFeatures[edge0].FaceId(0);
    if(face0 == triangleIndices.first)
        face0 = this->edgeFeatures[edge0].FaceId(1);

    size_t edge1 = this->triangleEdgeFeatures[triangleIndices.second].EdgeId(vertexPair.second);
    size_t e11 = this->triangleEdgeFeatures[triangleIndices.second].EdgeId((vertexPair.second+2)%3);

    if (this->edgeFeatures[edge0].Covertex(this->edgeFeatures[edge1]))
    {   // Swap
        size_t tmp = edge1;
        edge1 = e11;
        e11 = tmp;
    }

    size_t face1 = this->edgeFeatures[e11].FaceId(0);
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
char TRTriangleModel<DataTypes>::getStatus2(const std::pair<size_t,size_t> &triangleIndices,
        std::pair<size_t,size_t> &vertexPair)
{
    size_t edge = this->triangleEdgeFeatures[triangleIndices.first].EdgeId((vertexPair.first+1)%3);
    size_t face0 = this->edgeFeatures[edge].FaceId(0);
    if(face0 == triangleIndices.first)
        face0 = this->edgeFeatures[edge].FaceId(1);

    edge = this->triangleEdgeFeatures[triangleIndices.second].EdgeId((vertexPair.second+1)%3);
    size_t face1 = this->edgeFeatures[edge].FaceId(0);
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
size_t TRTriangleModel<DataTypes>::covertexFace(const std::pair<size_t,size_t> &triangleIndices,
        std::pair<size_t,size_t> &vertexPair)
{
    Element face[2] = { Element(this,triangleIndices.first),
                        Element(this,triangleIndices.second)
                      };
    return face[0].covertex(face[1],vertexPair);
}

template<class DataTypes>
void TRTriangleModel<DataTypes>::computeContinuousBoundingTree(double dt, size_t maxDepth)
{
    PolytopeModel* cubeModel = createPrevious<PolytopeModel>();

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
            const Vector3 pt1v = pt1 + v[t.p1Index()]*dt;
            const Vector3 pt2v = pt2 + v[t.p2Index()]*dt;
            const Vector3 pt3v = pt3 + v[t.p3Index()]*dt;

            if (calcNormals)
            {
                // Also recompute normal vector
                t.n() = (pt2-pt1).cross(pt3-pt1).normalize();
            }
            cubeModel->setParentOf(i,pt1);
            cubeModel->setParentOf(i,pt2);
            cubeModel->setParentOf(i,pt3);

            cubeModel->setParentOf(i,pt1v);
            cubeModel->setParentOf(i,pt2v);
            cubeModel->setParentOf(i,pt3v);

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

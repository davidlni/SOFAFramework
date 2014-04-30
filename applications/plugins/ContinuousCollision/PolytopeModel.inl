#ifndef SOFA_COMPONENT_COLLISION_POLYTOPEMODEL_INL
#define SOFA_COMPONENT_COLLISION_POLYTOPEMODEL_INL
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
#include "PolytopeModel.h"

#include <sofa/core/visual/VisualParams.h>
#include <sofa/simulation/common/Simulation.h>
#include <sofa/core/ObjectFactory.h>
#include <algorithm>
#include <math.h>
#include <sofa/helper/system/gl.h>

namespace sofa
{

namespace component
{

namespace collision
{

template<typename TDataTypes, size_t K>
TPolytopeModel<TDataTypes,K>::TPolytopeModel()
{
    enum_type = KDOP_TYPE;
}

template<typename TDataTypes, size_t K>
void TPolytopeModel<TDataTypes,K>::resize(int s)
{
    int previousSize = this->size;
    if (s == previousSize) return;
    // reset parent
    CollisionModel* parent = getPrevious();
    while(parent != NULL)
    {
        parent->resize(0);
        parent = parent->getPrevious();
    }
    this->core::CollisionModel::resize(s);
    this->elems.resize(s);
    this->parentOf.resize(s);
    // set additional indices
    for (int i=previousSize; i<s; ++i)
    {
        this->elems[i].children.first=core::CollisionElementIterator(getNext(), i);
        this->elems[i].children.second=core::CollisionElementIterator(getNext(), i+1);
        this->parentOf[i] = i;
    }
}

template<typename TDataTypes, size_t K>
void TPolytopeModel<TDataTypes,K>::setParentOf(int childIndex, const Vector3& p)
{
    int i = parentOf[childIndex];
    elems[i] += p;
}

template<typename TDataTypes, size_t K>
void TPolytopeModel<TDataTypes,K>::setLeafPolytope(int cubeIndex, int childIndex)
{
    parentOf[childIndex] = cubeIndex;
    this->elems[cubeIndex].children.first=core::CollisionElementIterator(getNext(), childIndex);
    this->elems[cubeIndex].children.second=core::CollisionElementIterator(getNext(), childIndex+1);
}

template<typename TDataTypes, size_t K>
void TPolytopeModel<TDataTypes,K>::setLeafPolytope(int cubeIndex, std::pair<core::CollisionElementIterator,core::CollisionElementIterator> children, const Vector3& p)
{
    elems[cubeIndex] += p;
    elems[cubeIndex].children = children;
}

template<typename TDataTypes, size_t K>
int TPolytopeModel<TDataTypes,K>::addPolytope(Polytope<TDataTypes,K> subcellsBegin, Polytope<TDataTypes,K> subcellsEnd)
{
    int i = this->size;
    elems.push_back(PolytopeData());
    PolytopeData &element = elems.back();
    element.subcells.first = subcellsBegin;
    element.subcells.second = subcellsEnd;
    element.children.first = core::CollisionElementIterator();
    element.children.second = core::CollisionElementIterator();
    this->core::CollisionModel::resize(elems.size());
    updatePolytope(element);
    return i;
}

template<typename TDataTypes, size_t K>
void TPolytopeModel<TDataTypes,K>::updatePolytope(PolytopeData &element)
{
    const std::pair<Element,Element>& subcells = element.subcells;
    if (subcells.first != subcells.second)
    {
        Element c = subcells.first;
        while(c != subcells.second)
        {
            element += static_cast<const PolytopeData&>(c.getDOP());
            ++c;
        }
    }
}

template<typename TDataTypes, size_t K>
void TPolytopeModel<TDataTypes,K>::updatePolytope(int index)
{
    this->updatePolytope(elems[index]);
}

template<typename TDataTypes, size_t K>
void TPolytopeModel<TDataTypes,K>::updatePolytopes()
{
    for (int i=0; i<size; i++)
        updatePolytope(i);
}

template<typename TDataTypes, size_t K>
std::pair<core::CollisionElementIterator,core::CollisionElementIterator> TPolytopeModel<TDataTypes,K>::getInternalChildren(int index) const
{
    return elems[index].subcells;
}

template<typename TDataTypes, size_t K>
std::pair<core::CollisionElementIterator,core::CollisionElementIterator> TPolytopeModel<TDataTypes,K>::getExternalChildren(int index) const
{
    return elems[index].children;
}

template<typename TDataTypes, size_t K>
bool TPolytopeModel<TDataTypes,K>::isLeaf( int index ) const
{
    return elems[index].children.first.valid();
}

template<typename TDataTypes, size_t K>
void TPolytopeModel<TDataTypes,K>::computeBoundingTree(int maxDepth)
{
    if(maxDepth <= 0)
        return;

    //sout << ">PolytopeModel::computeBoundingTree("<<maxDepth<<")"<<sendl;
    std::list<TPolytopeModel<TDataTypes,K>*> levels;
    levels.push_front(createPrevious<PolytopeModel>());
    for (int i=0; i<maxDepth; i++)
        levels.push_front(levels.front()->createPrevious<TPolytopeModel<TDataTypes,K> >());

    TPolytopeModel<TDataTypes,K>* root = levels.front();
    //if (isStatic() && root->getPrevious() == NULL && !root->empty()) return; // No need to recompute BBox if immobile

    if (root->empty() || root->getPrevious() != NULL)
    {
        // Tree must be reconstructed
        // First remove extra levels
        while(root->getPrevious()!=NULL)
        {
            core::CollisionModel::SPtr m = root->getPrevious();
            root->setPrevious(m->getPrevious());
            if (m->getMaster()) m->getMaster()->removeSlave(m);
            m.reset();
        }

        // Then clear all existing levels
        for (std::list<TPolytopeModel<TDataTypes,K>*>::iterator it = levels.begin(); it != levels.end(); ++it)
            (*it)->resize(0);

        // Then build root cell
        //sout << "PolytopeModel: add root cube"<<sendl;
        root->addPolytope(Element(this,0),Element(this,size));

        // Construct tree by splitting cells along their biggest dimension
        std::list<TPolytopeModel<TDataTypes,K>*>::iterator it = levels.begin();
        TPolytopeModel<TDataTypes,K>* level = *it;
        ++it;
        int lvl = 0;
        while(it != levels.end())
        {
            //sout << "PolytopeModel: split level "<<lvl<<sendl;
            TPolytopeModel<TDataTypes,K>* clevel = *it;
            clevel->elems.reserve(level->size*2);
            for(Element cell = Element(level->begin()); cell != level->end() ; ++cell)
            {
                const std::pair<Element,Element>& subcells = cell.subcells();
                int ncells = subcells.second.getIndex() - subcells.first.getIndex();
                //sout << "PolytopeModel: level "<<lvl<<" cell "<<cell.getIndex()<<": current subcells "<<subcells.first.getIndex() << " - "<<subcells.second.getIndex()<<sendl;
                if (ncells > 4)
                {
                    // Only split cells with more than 4 childs
                    // Find the biggest dimension
                    int splitAxis = cell.getDOP().GetSplitAxis();
                    int middle = subcells.first.getIndex()+(ncells+1)/2;

                    // Separate cells on each side of the median cell
                    PolytopeSortPredicate sortPredicate(splitAxis);
                    std::sort(this->elems.begin()+subcells.first.getIndex(),this->elems.begin()+subcells.second.getIndex(), sortPredicate);

                    // Create the two new subcells
                    Element cmiddle(this, middle);
                    int c1 = clevel->addPolytope(subcells.first, cmiddle);
                    int c2 = clevel->addPolytope(cmiddle, subcells.second);
                    //sout << "L"<<lvl<<" cell "<<cell.getIndex()<<" split along "<<(splitAxis==0?'X':splitAxis==1?'Y':'Z')<<" in cell "<<c1<<" size "<<middle-subcells.first.getIndex()<<" and cell "<<c2<<" size "<<subcells.second.getIndex()-middle<<"."<<sendl;
                    level->elems[cell.getIndex()].subcells.first = Element(clevel,c1);
                    level->elems[cell.getIndex()].subcells.second = Element(clevel,c2+1);
                }
            }
            ++it;
            level = clevel;
            ++lvl;
        }
        if (!parentOf.empty())
        {
            // Finally update parentOf to reflect new cell order
            for (int i=0; i<size; i++)
                parentOf[elems[i].children.first.getIndex()] = i;
        }
    }
    else
    {
        // Simply update the existing tree, starting from the bottom
        int lvl = 0;
        for (std::list<TPolytopeModel<TDataTypes,K>*>::reverse_iterator it = levels.rbegin(); it != levels.rend(); it++)
        {
            //sout << "PolytopeModel: update level "<<lvl<<sendl;
            (*it)->updatePolytopes();
            ++lvl;
        }
    }
    //sout << "<PolytopeModel::computeBoundingTree("<<maxDepth<<")"<<sendl;
}


template<typename TDataTypes, size_t K>
void TPolytopeModel<TDataTypes,K>::draw(const core::visual::VisualParams* , int index)
{
#ifndef SOFA_NO_OPENGL
    const Vector3& vmin = elems[index].minBBox;
    const Vector3& vmax = elems[index].maxBBox;

    glBegin(GL_LINES);
    {
        glVertex3d(vmin[0], vmin[1], vmin[2]);
        glVertex3d(vmin[0], vmin[1], vmax[2]);
        glVertex3d(vmin[0], vmax[1], vmin[2]);
        glVertex3d(vmin[0], vmax[1], vmax[2]);
        glVertex3d(vmax[0], vmin[1], vmin[2]);
        glVertex3d(vmax[0], vmin[1], vmax[2]);
        glVertex3d(vmax[0], vmax[1], vmin[2]);
        glVertex3d(vmax[0], vmax[1], vmax[2]);

        glVertex3d(vmin[0], vmin[1], vmin[2]);
        glVertex3d(vmin[0], vmax[1], vmin[2]);
        glVertex3d(vmin[0], vmin[1], vmax[2]);
        glVertex3d(vmin[0], vmax[1], vmax[2]);
        glVertex3d(vmax[0], vmin[1], vmin[2]);
        glVertex3d(vmax[0], vmax[1], vmin[2]);
        glVertex3d(vmax[0], vmin[1], vmax[2]);
        glVertex3d(vmax[0], vmax[1], vmax[2]);

        glVertex3d(vmin[0], vmin[1], vmin[2]);
        glVertex3d(vmax[0], vmin[1], vmin[2]);
        glVertex3d(vmin[0], vmax[1], vmin[2]);
        glVertex3d(vmax[0], vmax[1], vmin[2]);
        glVertex3d(vmin[0], vmin[1], vmax[2]);
        glVertex3d(vmax[0], vmin[1], vmax[2]);
        glVertex3d(vmin[0], vmax[1], vmax[2]);
        glVertex3d(vmax[0], vmax[1], vmax[2]);
    }
    glEnd();
#endif /* SOFA_NO_OPENGL */
}

template<typename TDataTypes, size_t K>
void TPolytopeModel<TDataTypes,K>::draw(const core::visual::VisualParams* vparams)
{
    if (!isActive() || !((getNext()==NULL)?vparams->displayFlags().getShowCollisionModels():vparams->displayFlags().getShowBoundingCollisionModels())) return;

    int level=0;
    CollisionModel* m = getPrevious();
    float color = 1.0f;
    while (m!=NULL)
    {
        m = m->getPrevious();
        ++level;
        color *= 0.5f;
    }
    Vec<4,float> c;
    if (isSimulated())
        c=Vec<4,float>(1.0f, 1.0f, 1.0f, color);
    else
        c=Vec<4,float>(1.0f, 1.0f, 1.0f, color);

    std::vector< Vector3 > points;
    for (int i=0; i<size; i++)
    {
        const Vector3& vmin = elems[i].minBBox;
        const Vector3& vmax = elems[i].maxBBox;

        points.push_back(Vector3(vmin[0], vmin[1], vmin[2]));
        points.push_back(Vector3(vmin[0], vmin[1], vmax[2]));
        points.push_back(Vector3(vmin[0], vmax[1], vmin[2]));
        points.push_back(Vector3(vmin[0], vmax[1], vmax[2]));
        points.push_back(Vector3(vmax[0], vmin[1], vmin[2]));
        points.push_back(Vector3(vmax[0], vmin[1], vmax[2]));
        points.push_back(Vector3(vmax[0], vmax[1], vmin[2]));
        points.push_back(Vector3(vmax[0], vmax[1], vmax[2]));

        points.push_back(Vector3(vmin[0], vmin[1], vmin[2]));
        points.push_back(Vector3(vmin[0], vmax[1], vmin[2]));
        points.push_back(Vector3(vmin[0], vmin[1], vmax[2]));
        points.push_back(Vector3(vmin[0], vmax[1], vmax[2]));
        points.push_back(Vector3(vmax[0], vmin[1], vmin[2]));
        points.push_back(Vector3(vmax[0], vmax[1], vmin[2]));
        points.push_back(Vector3(vmax[0], vmin[1], vmax[2]));
        points.push_back(Vector3(vmax[0], vmax[1], vmax[2]));

        points.push_back(Vector3(vmin[0], vmin[1], vmin[2]));
        points.push_back(Vector3(vmax[0], vmin[1], vmin[2]));
        points.push_back(Vector3(vmin[0], vmax[1], vmin[2]));
        points.push_back(Vector3(vmax[0], vmax[1], vmin[2]));
        points.push_back(Vector3(vmin[0], vmin[1], vmax[2]));
        points.push_back(Vector3(vmax[0], vmin[1], vmax[2]));
        points.push_back(Vector3(vmin[0], vmax[1], vmax[2]));
        points.push_back(Vector3(vmax[0], vmax[1], vmax[2]));
    }

    vparams->drawTool()->drawLines(points, 1, Vec<4,float>(c));


    if (getPrevious()!=NULL)
        getPrevious()->draw(vparams);
}

} // namespace collision

} // namespace component

} // namespace sofa

#endif

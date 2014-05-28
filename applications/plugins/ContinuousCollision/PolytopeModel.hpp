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

#include <stack>

#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/ObjectFactory.h>

#include "PolytopeModel.h"

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
void TPolytopeModel<TDataTypes,K>::resize(size_t s)
{
    if (s == this->getSize()) return;
    // reset parent
    CollisionModel* parent = getPrevious();
    while(parent != NULL)
    {
        parent->resize(0);
        parent = parent->getPrevious();
    }
    this->polytopes.resize(s);
    this->parentOf.resize(s);
    // set additional indices
    for (size_t i = this->getSize(); i < s; ++i)
    {
        this->polytopes[i].children.first=core::CollisionElementIterator(getNext(), i);
        this->polytopes[i].children.second=core::CollisionElementIterator(getNext(), i+1);
        this->parentOf[i] = i;
    }
    this->getSize() = s;
}

template<typename TDataTypes, size_t K>
inline void TPolytopeModel<TDataTypes,K>::enlarge(const index_type &childIndex, const Real& p)
{
    index_type i = this->parentOf[childIndex];
    this->polytopes[i] += p;
}

template<typename TDataTypes, size_t K>
inline void TPolytopeModel<TDataTypes,K>::setParentOf(const index_type &childIndex, const Vector3& p)
{
  index_type i = this->parentOf[childIndex];
  this->polytopes[i] += p;
}

template<typename TDataTypes, size_t K>
void TPolytopeModel<TDataTypes,K>::setLeafPolytope(const index_type &polytopeIndex,
                                                   const index_type &childIndex)
{
    this->parentOf[childIndex] = polytopeIndex;
    this->polytopes[polytopeIndex].children.first=core::CollisionElementIterator(getNext(), childIndex);
    this->polytopes[polytopeIndex].children.second=core::CollisionElementIterator(getNext(), childIndex+1);
}

template<typename TDataTypes, size_t K>
void TPolytopeModel<TDataTypes,K>::setLeafPolytope(const index_type &polytopeIndex,
                                                   const std::pair<core::CollisionElementIterator,core::CollisionElementIterator> &children,
                                                   const Vector3& p)
{
    this->polytopes[polytopeIndex] += p;
    this->polytopes[polytopeIndex].children = children;
}

template<typename TDataTypes, size_t K>
size_t TPolytopeModel<TDataTypes,K>::addPolytope(const TPolytope<TDataTypes,K> &subcellsBegin,
                                                 const TPolytope<TDataTypes,K> &subcellsEnd)
{
    size_t i = this->getSize();
    this->polytopes.push_back(PolytopeData());
    PolytopeData &element = this->polytopes.back();
    element.subcells.first = subcellsBegin;
    element.subcells.second = subcellsEnd;
    element.children.first = core::CollisionElementIterator();
    element.children.second = core::CollisionElementIterator();
    this->setSize(this->polytopes.size());
    this->updatePolytope(element);
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
            element += c.getElement();
            ++c;
        }
    }
}

template<typename TDataTypes, size_t K>
void TPolytopeModel<TDataTypes,K>::updatePolytope(const index_type &index)
{
    this->updatePolytope(this->polytopes[index]);
}

template<typename TDataTypes, size_t K>
void TPolytopeModel<TDataTypes,K>::updatePolytopes()
{
    for (size_t i=0, end = this->getSize(); i < end; i++)
        this->updatePolytope(i);
}

template<typename TDataTypes, size_t K>
std::pair<core::CollisionElementIterator,core::CollisionElementIterator>
TPolytopeModel<TDataTypes,K>::getInternalChildren(size_t index) const
{
    return this->polytopes[index].subcells;
}

template<typename TDataTypes, size_t K>
std::pair<core::CollisionElementIterator,core::CollisionElementIterator>
TPolytopeModel<TDataTypes,K>::getExternalChildren(size_t index) const
{
    return this->polytopes[index].children;
}

template<typename TDataTypes, size_t K>
bool TPolytopeModel<TDataTypes,K>::isLeaf( size_t index ) const
{
    return this->polytopes[index].children.first.valid();
}

template<typename TDataTypes, size_t K>
void TPolytopeModel<TDataTypes,K>::computeBoundingTree(size_t maxDepth)
{
    if(maxDepth <= 0)
        return;

    helper::vector<MyType*> levelStack;
    levelStack.push_back(this->createPrevious<MyType>());
    for (index_type i = 0; i < maxDepth; i++)
      levelStack.push_back(levelStack.back()->createPrevious<MyType>());

    MyType* root = levelStack.back();
    if (root->empty() || root->getPrevious())
    {
        // Remove the root from stack
        levelStack.pop_back();

        // Tree must be reconstructed
        // First remove extra levels
        while(root->getPrevious())
        {
            core::CollisionModel::SPtr m = root->getPrevious();
            root->setPrevious(m->getPrevious());
            if (m->getMaster()) m->getMaster()->removeSlave(m);
            m.reset();
        }

        // Build root cell
        // Root consist of one polytope at the 0 level with all the elements
        root->resize(0);
        root->addPolytope(Element(this,0),Element(this,this->getSize()));

        // Root is at level 0
        MyType* parentLevel = root;

        // Construct tree by splitting cells along their biggest dimension
        while(!levelStack.empty())
        {
            MyType* currentLevel = levelStack.back();
            levelStack.pop_back();
            currentLevel->resize(0);
            currentLevel->polytopes.reserve(parentLevel->getSize()*2);
            for(Element cell = Element(parentLevel->begin()); cell != parentLevel->end() ; ++cell)
            {
                const std::pair<Element,Element>& subcells = cell.subcells();
                index_type ncells = subcells.second.getIndex() - subcells.first.getIndex();
                // Only split cells with more than 4 childs
                if (ncells > 4)
                {
                    // Find the biggest dimension
                    index_type splitAxis = cell.getElement().GetSplitAxis();
                    Real center = cell.getElement().GetCenter(splitAxis);

                    // Separate cells on each side of the median cell
                    typename ElementListType::iterator start, mid, finish;
                    start = this->polytopes.begin()+subcells.first.getIndex();
                    finish = this->polytopes.begin()+subcells.second.getIndex();
                    mid = std::partition(start,finish,PolytopeSortPredicate(splitAxis,center));

                    index_type middleIndex = std::distance(start,mid);

                    // Create the two new polytops in current level containing the splited indices
                    Element cmiddle(this, middleIndex);
                    index_type c1 = currentLevel->addPolytope(subcells.first, cmiddle);
                    index_type c2 = currentLevel->addPolytope(cmiddle, subcells.second);

                    // Update subcells of the parent level with the two newly created cells
                    parentLevel->polytopes[cell.getIndex()].subcells.first = Element(currentLevel,c1);
                    parentLevel->polytopes[cell.getIndex()].subcells.second = Element(currentLevel,c2+1);
                }
            }
            parentLevel = currentLevel;
        }
        if (!parentOf.empty())
        {
            // Finally update parentOf to reflect new cell order
            for (size_t i=0, end = this->getSize(); i < end; i++)
                parentOf[polytopes[i].children.first.getIndex()] = i;
        }
    }
    else
    {
        // Update the existing tree, starting from the bottom
        for (size_t i = 0; i < levelStack.size(); ++i)
        {
            sout << "PolytopeModel: update level "<<i<<sendl;
            levelStack[i]->updatePolytopes();
        }
    }
    //sout << "<PolytopeModel::computeBoundingTree("<<maxDepth<<")"<<sendl;
}

template<typename TDataTypes, size_t K>
void TPolytopeModel<TDataTypes,K>::draw(const core::visual::VisualParams* , size_t /*index*/)
{
// #ifndef SOFA_NO_OPENGL
//     const Vector3& vmin = polytopes[index].minBBox;
//     const Vector3& vmax = polytopes[index].maxBBox;
//
//     glBegin(GL_LINES);
//     {
//         glVertex3d(vmin[0], vmin[1], vmin[2]);
//         glVertex3d(vmin[0], vmin[1], vmax[2]);
//         glVertex3d(vmin[0], vmax[1], vmin[2]);
//         glVertex3d(vmin[0], vmax[1], vmax[2]);
//         glVertex3d(vmax[0], vmin[1], vmin[2]);
//         glVertex3d(vmax[0], vmin[1], vmax[2]);
//         glVertex3d(vmax[0], vmax[1], vmin[2]);
//         glVertex3d(vmax[0], vmax[1], vmax[2]);
//
//         glVertex3d(vmin[0], vmin[1], vmin[2]);
//         glVertex3d(vmin[0], vmax[1], vmin[2]);
//         glVertex3d(vmin[0], vmin[1], vmax[2]);
//         glVertex3d(vmin[0], vmax[1], vmax[2]);
//         glVertex3d(vmax[0], vmin[1], vmin[2]);
//         glVertex3d(vmax[0], vmax[1], vmin[2]);
//         glVertex3d(vmax[0], vmin[1], vmax[2]);
//         glVertex3d(vmax[0], vmax[1], vmax[2]);
//
//         glVertex3d(vmin[0], vmin[1], vmin[2]);
//         glVertex3d(vmax[0], vmin[1], vmin[2]);
//         glVertex3d(vmin[0], vmax[1], vmin[2]);
//         glVertex3d(vmax[0], vmax[1], vmin[2]);
//         glVertex3d(vmin[0], vmin[1], vmax[2]);
//         glVertex3d(vmax[0], vmin[1], vmax[2]);
//         glVertex3d(vmin[0], vmax[1], vmax[2]);
//         glVertex3d(vmax[0], vmax[1], vmax[2]);
//     }
//     glEnd();
// #endif /* SOFA_NO_OPENGL */
}

template<typename TDataTypes, size_t K>
void TPolytopeModel<TDataTypes,K>::draw(const core::visual::VisualParams* /*vparams*/)
{
//     if (!isActive() || !((getNext()==NULL)?vparams->displayFlags().getShowCollisionModels():vparams->displayFlags().getShowBoundingCollisionModels())) return;
//
//     size_t level=0;
//     CollisionModel* m = getPrevious();
//     float color = 1.0f;
//     while (m!=NULL)
//     {
//         m = m->getPrevious();
//         ++level;
//         color *= 0.5f;
//     }
//     Vec<4,float> c;
//     if (isSimulated())
//         c=Vec<4,float>(1.0f, 1.0f, 1.0f, color);
//     else
//         c=Vec<4,float>(1.0f, 1.0f, 1.0f, color);
//
//     std::vector< Vector3 > points;
//     for (size_t i=0; i<size; i++)
//     {
//         const Vector3& vmin = polytopes[i].minBBox;
//         const Vector3& vmax = polytopes[i].maxBBox;
//
//         points.push_back(Vector3(vmin[0], vmin[1], vmin[2]));
//         points.push_back(Vector3(vmin[0], vmin[1], vmax[2]));
//         points.push_back(Vector3(vmin[0], vmax[1], vmin[2]));
//         points.push_back(Vector3(vmin[0], vmax[1], vmax[2]));
//         points.push_back(Vector3(vmax[0], vmin[1], vmin[2]));
//         points.push_back(Vector3(vmax[0], vmin[1], vmax[2]));
//         points.push_back(Vector3(vmax[0], vmax[1], vmin[2]));
//         points.push_back(Vector3(vmax[0], vmax[1], vmax[2]));
//
//         points.push_back(Vector3(vmin[0], vmin[1], vmin[2]));
//         points.push_back(Vector3(vmin[0], vmax[1], vmin[2]));
//         points.push_back(Vector3(vmin[0], vmin[1], vmax[2]));
//         points.push_back(Vector3(vmin[0], vmax[1], vmax[2]));
//         points.push_back(Vector3(vmax[0], vmin[1], vmin[2]));
//         points.push_back(Vector3(vmax[0], vmax[1], vmin[2]));
//         points.push_back(Vector3(vmax[0], vmin[1], vmax[2]));
//         points.push_back(Vector3(vmax[0], vmax[1], vmax[2]));
//
//         points.push_back(Vector3(vmin[0], vmin[1], vmin[2]));
//         points.push_back(Vector3(vmax[0], vmin[1], vmin[2]));
//         points.push_back(Vector3(vmin[0], vmax[1], vmin[2]));
//         points.push_back(Vector3(vmax[0], vmax[1], vmin[2]));
//         points.push_back(Vector3(vmin[0], vmin[1], vmax[2]));
//         points.push_back(Vector3(vmax[0], vmin[1], vmax[2]));
//         points.push_back(Vector3(vmin[0], vmax[1], vmax[2]));
//         points.push_back(Vector3(vmax[0], vmax[1], vmax[2]));
//     }
//
//     vparams->drawTool()->drawLines(points, 1, c);
//
//     if (getPrevious()!=NULL)
//         getPrevious()->draw(vparams);
}

} // namespace collision

} // namespace component

} // namespace sofa

#endif

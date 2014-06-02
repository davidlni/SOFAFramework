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

#include "initContinuousCollision.h"
#include "ContinuousDetection.h"
#include "RTriangleModel.h"
#include <sofa/core/visual/VisualParams.h>
#include <sofa/component/collision/CapsuleModel.h>
#include <sofa/component/collision/Sphere.h>
#include <sofa/component/collision/Triangle.h>
#include <sofa/component/collision/Line.h>
#include <sofa/component/collision/Point.h>
#include <sofa/helper/FnDispatcher.h>
#include <sofa/core/ObjectFactory.h>
#include <map>
#include <queue>
#include <stack>

#include <sofa/helper/system/gl.h>
#include <sofa/helper/system/glut.h>

namespace sofa
{

namespace component
{

namespace collision
{

using namespace defaulttype;
using namespace helper;
using namespace collision;

SOFA_DECL_CLASS(ContinuousDetection)

int ContinuousDetectionClass = core::RegisterObject("Continuous Collision Detection")
                               .add< ContinuousDetection >()
                               ;

using namespace core::objectmodel;


class ContinuousDetection::Impl
{
public:
  void bufferAdjacentList(RTriangleModel *cm);
  void setOrphans(RTriangleModel *cm);
};

inline void ContinuousDetection::Impl::bufferAdjacentList(RTriangleModel* cm)
{
  cm->bufferAdjacentLists();
}

inline void ContinuousDetection::Impl::setOrphans(RTriangleModel* cm)
{
  cm->setOrphans();
}

ContinuousDetection::ContinuousDetection()
    : bDraw(initData(&bDraw, false, "draw", "enable/disable display of results"))
{
  this->pimpl = new Impl;
}

ContinuousDetection::~ContinuousDetection()
{
  delete this->pimpl;
}

void ContinuousDetection::init()
{
    reinit();
}

void ContinuousDetection::reinit()
{

}

void ContinuousDetection::addCollisionModel(core::CollisionModel *cm)
{
    //sout<<"--------- add Collision Model : "<<cm->getLast()->getName()<<" -------"<<sendl;
    if (cm->empty())
        return;


    if (cm->isSimulated() && cm->getLast()->canCollideWith(cm->getLast()))
    {
        // self collision
        //sout << "Test broad phase Self "<<cm->getLast()->getName()<<sendl;
        bool swapModels = false;
        core::collision::ElementIntersector* intersector = intersectionMethod->findIntersector(cm, cm, swapModels);
        if (intersector != NULL)
            if (intersector->canIntersect(cm->begin(), cm->begin()))
            {
                //sout << "Broad phase Self "<<cm->getLast()->getName()<<sendl;
                cmPairs.push_back(std::make_pair(cm, cm));
            }
    }

    for (helper::vector<core::CollisionModel*>::iterator it = collisionModels.begin(); it != collisionModels.end(); ++it)
    {
        core::CollisionModel* cm2 = *it;

        if (!cm->isSimulated() && !cm2->isSimulated())
        {
            continue;
        }

        if (!keepCollisionBetween(cm->getLast(), cm2->getLast()))
            continue;

        bool swapModels = false;
        core::collision::ElementIntersector* intersector = intersectionMethod->findIntersector(cm, cm2, swapModels);
        if (intersector == NULL)
            continue;

        core::CollisionModel* cm1 = (swapModels?cm2:cm);
        cm2 = (swapModels?cm:cm2);

        // Here we assume a single root element is present in both models
        if (intersector->canIntersect(cm1->begin(), cm2->begin()))
        {
            sout << "Broad phase "<<cm1->getLast()->getName()<<" - "<<cm2->getLast()->getName()<<sendl;
            cmPairs.push_back(std::make_pair(cm1, cm2));
        }
    }

    collisionModels.push_back(cm);
}

bool ContinuousDetection::keepCollisionBetween(core::CollisionModel *cm1, core::CollisionModel *cm2)
{
    if (!cm1->canCollideWith(cm2) || !cm2->canCollideWith(cm1))
    {
        return false;
    }

    return true;
}

class MirrorIntersector : public core::collision::ElementIntersector
{
public:
  core::collision::ElementIntersector* intersector;

  /// Test if 2 elements can collide. Note that this can be conservative (i.e. return true even when no collision is present)
  virtual bool canIntersect(core::CollisionElementIterator elem1, core::CollisionElementIterator elem2)
  {
    return intersector->canIntersect(elem2, elem1);
  }

  /// Begin intersection tests between two collision models. Return the number of contacts written in the contacts vector.
  /// If the given contacts vector is NULL, then this method should allocate it.
  virtual int beginIntersect(core::CollisionModel* model1, core::CollisionModel* model2, core::collision::DetectionOutputVector*& contacts)
  {
    return intersector->beginIntersect(model2, model1, contacts);
  }

  /// Compute the intersection between 2 elements. Return the number of contacts written in the contacts vector.
  virtual int intersect(core::CollisionElementIterator elem1, core::CollisionElementIterator elem2, core::collision::DetectionOutputVector* contacts)
  {
    return intersector->intersect(elem2, elem1, contacts);
  }

  /// End intersection tests between two collision models. Return the number of contacts written in the contacts vector.
  virtual int endIntersect(core::CollisionModel* model1, core::CollisionModel* model2, core::collision::DetectionOutputVector* contacts)
  {
    return intersector->endIntersect(model2, model1, contacts);
  }

  virtual std::string name() const
  {
    return intersector->name() + std::string("<SWAP>");
  }

};

void ContinuousDetection::addCollisionPair(const std::pair<core::CollisionModel*, core::CollisionModel*>& cmPair)
{
    typedef std::pair<core::CollisionElementIterator,core::CollisionElementIterator> ElementIteratorPair;
    typedef std::pair< ElementIteratorPair, ElementIteratorPair > TestPair;

    core::CollisionModel *cm1 = cmPair.first; //->getNext();
    core::CollisionModel *cm2 = cmPair.second; //->getNext();

    if (!cm1->isSimulated() && !cm2->isSimulated())
        return;

    if (cm1->empty() || cm2->empty())
        return;


    // Get the leaf CollisionModels
    core::CollisionModel *finalcm1 = cm1->getLast();
    core::CollisionModel *finalcm2 = cm2->getLast();

    //sout << "Final phase "<<gettypename(typeid(*finalcm1))<<" - "<<gettypename(typeid(*finalcm2))<<sendl;
    // Find the interrsection method for the leaf CollisionModels
    bool swapModels = false;
    core::collision::ElementIntersector* finalintersector = intersectionMethod->findIntersector(finalcm1, finalcm2, swapModels);
    if (finalintersector == NULL)
        return;

    if (swapModels)
    {
        core::CollisionModel* tmp;
        tmp = cm1;
        cm1 = cm2;
        cm2 = tmp;
        tmp = finalcm1;
        finalcm1 = finalcm2;
        finalcm2 = tmp;
    }

    const bool self = (finalcm1->getContext() == finalcm2->getContext());
    if (self)
    {
       sout << "SELF: Final intersector " << finalintersector->name() << " for "<<finalcm1->getName()<<" - "<<finalcm2->getName()<<sendl;
    }

    core::collision::DetectionOutputVector*& outputs = this->getDetectionOutputs(finalcm1, finalcm2);

    finalintersector->beginIntersect(finalcm1, finalcm2, outputs);//creates outputs if null

    if (finalcm1 == cm1 || finalcm2 == cm2)
    {
        // The last model also contains the root element -> it does not only contains the final level of the tree
        finalcm1 = NULL;
        finalcm2 = NULL;
        finalintersector = NULL;
    }

    std::queue< TestPair > externalCells;

    ElementIteratorPair internalChildren1 = cm1->begin().getInternalChildren();
    ElementIteratorPair internalChildren2 = cm2->begin().getInternalChildren();
    ElementIteratorPair externalChildren1 = cm1->begin().getExternalChildren();
    ElementIteratorPair externalChildren2 = cm2->begin().getExternalChildren();
    if (internalChildren1.first != internalChildren1.second)
    {
        if (internalChildren2.first != internalChildren2.second)
            externalCells.push(std::make_pair(internalChildren1,internalChildren2));
        if (externalChildren2.first != externalChildren2.second)
            externalCells.push(std::make_pair(internalChildren1,externalChildren2));
    }

    if (externalChildren1.first != externalChildren1.second)
    {
        if (internalChildren2.first != internalChildren2.second)
            externalCells.push(std::make_pair(externalChildren1,internalChildren2));
        if (externalChildren2.first != externalChildren2.second)
            externalCells.push(std::make_pair(externalChildren1,externalChildren2));
    }

    core::collision::ElementIntersector* intersector = NULL;
    MirrorIntersector mirror;

    // Force later init of intersector
    cm1 = NULL;
    cm2 = NULL;

    while (!externalCells.empty())
    {
        TestPair root = externalCells.front();
        externalCells.pop();

        // If the CollisionElements do not belong to cm1 and cm2, update cm1 and cm2
        if (cm1 != root.first.first.getCollisionModel() || cm2 != root.second.first.getCollisionModel())
        {
            cm1 = root.first.first.getCollisionModel();
            cm2 = root.second.first.getCollisionModel();
            if (!cm1 || !cm2) continue;
            intersector = intersectionMethod->findIntersector(cm1, cm2, swapModels);

            if (intersector == NULL)
            {
                sout << "ContinuousDetection: Error finding intersector " << intersectionMethod->getName() << " for "<<cm1->getClassName()<<" - "<<cm2->getClassName()<<sendl;
            }
            //else sout << "ContinuousDetection: intersector " << intersector->name() << " for " << intersectionMethod->getName() << " for "<<gettypename(typeid(*cm1))<<" - "<<gettypename(typeid(*cm2))<<sendl;
            if (swapModels)
            {
                mirror.intersector = intersector;
                intersector = &mirror;
            }
        }
        if (intersector == NULL)
            continue;
        std::stack< TestPair > internalCells;
        internalCells.push(root);

        while (!internalCells.empty())
        {
            TestPair current = internalCells.top();
            internalCells.pop();

            core::CollisionElementIterator begin1 = current.first.first;
            core::CollisionElementIterator end1 = current.first.second;
            core::CollisionElementIterator begin2 = current.second.first;
            core::CollisionElementIterator end2 = current.second.second;

            if (begin1.getCollisionModel() == finalcm1 && begin2.getCollisionModel() == finalcm2)
            {
                // Final collision pairs
                for (core::CollisionElementIterator it1 = begin1; it1 != end1; ++it1)
                {
                    for (core::CollisionElementIterator it2 = begin2; it2 != end2; ++it2)
                    {
                        if (!self && it1.canCollideWith(it2))
                        {
                            intersector->intersect(it1,it2,outputs);
                        }
                    }
                }
            }
            else
            {
                for (core::CollisionElementIterator it1 = begin1; it1 != end1; ++it1)
                {
                    for (core::CollisionElementIterator it2 = begin2; it2 != end2; ++it2)
                    {
                        if (intersector->canIntersect(it1,it2))
                        {
                            // Need to test recursively
                            // Note that an element cannot have both internal and external children
                            TestPair newInternalTests(it1.getInternalChildren(),it2.getInternalChildren());
                            TestPair newExternalTests(it1.getExternalChildren(),it2.getExternalChildren());
                            if (newInternalTests.first.first != newInternalTests.first.second)
                            {
                                if (newInternalTests.second.first != newInternalTests.second.second)
                                {
                                    internalCells.push(newInternalTests);
                                }
                                else
                                {
                                    newInternalTests.second.first = it2;
                                    newInternalTests.second.second = it2;
                                    ++newInternalTests.second.second;
                                    internalCells.push(newInternalTests);
                                }
                            }
                            else
                            {
                                if (newInternalTests.second.first != newInternalTests.second.second)
                                {
                                    newInternalTests.first.first = it1;
                                    newInternalTests.first.second = it1;
                                    ++newInternalTests.first.second;
                                    internalCells.push(newInternalTests);
                                }
                                else
                                {
                                    // If no internal children exist then we are at the leafs
                                    // Test leafs
                                    if (newExternalTests.first.first != newExternalTests.first.second)
                                    {
                                        if (newExternalTests.second.first != newExternalTests.second.second)
                                        {
                                            if (newExternalTests.first.first.getCollisionModel() == finalcm1 && newExternalTests.second.first.getCollisionModel() == finalcm2)
                                            {
                                                core::CollisionElementIterator begin1 = newExternalTests.first.first;
                                                core::CollisionElementIterator end1 = newExternalTests.first.second;
                                                core::CollisionElementIterator begin2 = newExternalTests.second.first;
                                                core::CollisionElementIterator end2 = newExternalTests.second.second;
                                                for (core::CollisionElementIterator it1 = begin1; it1 != end1; ++it1)
                                                {
                                                    for (core::CollisionElementIterator it2 = begin2; it2 != end2; ++it2)
                                                    {
                                                        // Final collision pair
                                                        if (!self && it1.canCollideWith(it2))
                                                            finalintersector->intersect(it1,it2,outputs);
                                                    }
                                                }
                                            }
                                            else
                                                externalCells.push(newExternalTests);
                                        }
                                        else
                                        {
                                            // only first element has external children
                                            // test them against the second element
                                            newExternalTests.second.first = it2;
                                            newExternalTests.second.second = it2;
                                            ++newExternalTests.second.second;
                                            externalCells.push(std::make_pair(newExternalTests.first, newInternalTests.second));
                                        }
                                    }
                                    else if (newExternalTests.second.first != newExternalTests.second.second)
                                    {
                                        // only first element has external children
                                        // test them against the first element
                                        newExternalTests.first.first = it1;
                                        newExternalTests.first.second = it1;
                                        ++newExternalTests.first.second;
                                        externalCells.push(std::make_pair(newExternalTests.first, newExternalTests.second));
                                    }
                                    else
                                    {
                                        // No child -> final collision pair
                                        if (!self && it1.canCollideWith(it2))
                                          finalintersector->intersect(it1,it2, outputs);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }


    finalintersector->endIntersect(finalcm1, finalcm2, outputs);

    //sout << "Narrow phase "<<cm1->getLast()->getName()<<"("<<gettypename(typeid(*cm1->getLast()))<<") - "<<cm2->getLast()->getName()<<"("<<gettypename(typeid(*cm2->getLast()))<<"): "<<elemPairs.size()-size0<<" contacts."<<sendl;
}

void ContinuousDetection::draw(const core::visual::VisualParams* )
{
    if (!bDraw.getValue()) return;
    /*
     *        glDisable(GL_LIGHTING);
     *        glColor3f(1.0, 0.0, 1.0);
     *        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
     *        glLineWidth(3);
     *        glPointSize(5);
     *
     *        for (DetectionOutputMap::iterator it = outputsMap.begin(); it!=outputsMap.end(); it++)
     *        {
     *            core::collision::DetectionOutputVector& outputs = it->second;
     *            for (core::collision::DetectionOutputVector::iterator it2 = outputs.begin(); it2!=outputs.end(); it2++)
     *            {
     *                it2->elem.first.draw(vparams);
     *                it2->elem.second.draw(vparams);
    }
    }
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glLineWidth(1);
    glPointSize(1);
    */
}

} // namespace collision

} // namespace component

} // namespace sofa


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

#include "CompliantAttachPerformer.h"
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/BaseMapping.h>
#include <sofa/component/collision/MouseInteractor.h>
#include <sofa/simulation/common/Simulation.h>
#include <iostream>
using std::cerr;
using std::endl;

#include "compliance/UniformCompliance.h"
#include <sofa/simulation/common/InitVisitor.h>

#include <sofa/component/mapping/DistanceFromTargetMapping.inl>


namespace sofa
{

namespace component
{

namespace collision
{


template <class DataTypes>
CompliantAttachPerformer<DataTypes>::CompliantAttachPerformer(BaseMouseInteractor *i):TInteractionPerformer<DataTypes>(i)
    , mapper(0)
    , mouseState(0)
{
    //    cerr<<"CompliantAttachPerformer<DataTypes>::CompliantAttachPerformer()" << endl;
    this->interactor->setMouseAttached(false);
    _compliance = 1e-3;
    _isCompliance = true;
}


template <class DataTypes>
CompliantAttachPerformer<DataTypes>::~CompliantAttachPerformer()
{
    //    cerr<<"CompliantAttachPerformer<DataTypes>::~CompliantAttachPerformer()" << endl;
    clear();
}

template <class DataTypes>
void CompliantAttachPerformer<DataTypes>::clear()
{
    pickedNode->removeChild(interactionNode);
    interactionNode.reset();  // delete the subgraph if no other reference to it

    if (mapper)
    {
        mapper->cleanup();
        delete mapper; mapper=NULL;
    }

    this->interactor->setDistanceFromMouse(0);
    this->interactor->setMouseAttached(false);
}


template <class DataTypes>
void CompliantAttachPerformer<DataTypes>::start()
{
    typedef sofa::component::collision::BaseContactMapper< DataTypes >        MouseContactMapper;


    if (interactionNode)  // previous interaction still holding
    {
        cerr<<"CompliantAttachPerformer<DataTypes>::start(), releasing previous interaction" << endl;
        clear();            // release it
        return;
    }
    //    cerr<<"CompliantAttachPerformer<DataTypes>::start()" << endl;


    //--------- Picked object
    BodyPicked picked=this->interactor->getBodyPicked();
    if (!picked.body && !picked.mstate) return;

    core::behavior::MechanicalState<DataTypes>* pickedState=dynamic_cast< core::behavior::MechanicalState<DataTypes>*  >(picked.mstate);
    pickedNode = dynamic_cast<simulation::Node*> (pickedState->getContext());
    assert(pickedNode);


    // get a picked state and a particle index
    core::behavior::MechanicalState<DataTypes>* mstateCollision=NULL;
    if (picked.body)  // picked a collision element: create a MechanicalState mapped under the collision model, to store the picked point - NOT TESTED YET !!!!
    {
//        std::cerr<<"mapped\n";

        mapper = MouseContactMapper::Create(picked.body);
        if (!mapper)
        {
            this->interactor->serr << "Problem with Mouse Mapper creation : " << this->interactor->sendl;
            return;
        }
        std::string name = "contactMouse";
        mstateCollision = mapper->createMapping(name.c_str());
        mapper->resize(1);

        const typename DataTypes::Coord pointPicked=picked.point;
        const int idx=picked.indexCollisionElement;
        typename DataTypes::Real r=0.0;

        pickedParticleIndex = mapper->addPoint(pointPicked, idx, r);
        mapper->update();

        // copy the tags of the collision model to the mapped state
        if (mstateCollision->getContext() != picked.body->getContext())
        {

            simulation::Node *mappedNode=(simulation::Node *) mstateCollision->getContext();
            simulation::Node *mainNode=(simulation::Node *) picked.body->getContext();
            core::behavior::BaseMechanicalState *mainDof=dynamic_cast<core::behavior::BaseMechanicalState *>(mainNode->getMechanicalState());
            const core::objectmodel::TagSet &tags=mainDof->getTags();
            for (core::objectmodel::TagSet::const_iterator it=tags.begin(); it!=tags.end(); ++it)
            {
                mstateCollision->addTag(*it);
                mappedNode->mechanicalMapping->addTag(*it);
            }
            mstateCollision->setName("AttachedPoint");
            mappedNode->mechanicalMapping->setName("MouseMapping");
        }
    }
    else // picked an existing particle
    {
//        std::cerr<<"Already\n";

        mstateCollision = dynamic_cast< core::behavior::MechanicalState<DataTypes>*  >(picked.mstate);
        pickedParticleIndex = picked.indexCollisionElement;
        //        cerr<<"CompliantAttachPerformer<DataTypes>::attach, pickedParticleIndex = " << pickedParticleIndex << endl;
        if (!mstateCollision)
        {
            this->interactor->serr << "incompatible MState during Mouse Interaction " << this->interactor->sendl;
            return;
        }
    }



    //-------- Mouse manipulator
    mouseMapping = this->interactor->core::objectmodel::BaseObject::template searchUp<sofa::core::BaseMapping>();
    this->mouseState = dynamic_cast<Point3dState*>(this->interactor->getMouseContainer());
//    typename Point3dState::ReadVecCoord xmouse = mouseState->readPositions();
//    typename Point3dState::Coord pointOnRay = mouseState->readPositions()[0];

    // set target point to closest point on the ray
    double distanceFromMouse=picked.rayLength;
    Ray ray = this->interactor->getMouseRayModel()->getRay(0);
    defaulttype::Vector3 pointOnRay = ray.origin() + ray.direction()*distanceFromMouse;
//    ray.setOrigin(pointOnRay);
    this->interactor->setMouseAttached(true);
    this->interactor->setDistanceFromMouse(distanceFromMouse);

    mouseState->writePositions()[0] = pointOnRay;

    //---------- Set up the interaction

    // look for existing interactions
    static const std::string distanceMappingName="InteractionDistanceFromTargetMapping_createdByCompliantAttachPerformer";

    interactionNode = pickedNode->createChild("InteractionDistanceNode");

    typedef component::container::MechanicalObject<DataTypes1> MechanicalObject1;
    typename MechanicalObject1::SPtr extensions = New<MechanicalObject1>();
    interactionNode->addObject(extensions);
    extensions->setName("extensionValues");

    distanceMapping = New<DistanceFromTargetMapping31>();
    distanceMapping->setModels(mstateCollision,extensions.get());
    interactionNode->addObject( distanceMapping );
    distanceMapping->setName(distanceMappingName.c_str());
    distanceMapping->createTarget(/*picked.indexCollisionElement*/ pickedParticleIndex, pointOnRay, /*(picked.point-pointOnRay).norm()*/ 0);
    distanceMapping->_arrowSize = _arrowSize;
    distanceMapping->_color = _color;

    //       cerr<<"CompliantAttachPerformer<DataTypes>::start(), create target of " << picked.indexCollisionElement << " at " <<  (*mstateCollision->getX())[picked.indexCollisionElement] << " to " << pointOnRay << ", distance = " << (picked.point-pointOnRay).norm() << endl;

    typedef forcefield::UniformCompliance<DataTypes1> UniformCompliance1;
    typename UniformCompliance1::SPtr compliance = New<UniformCompliance1>();
    compliance->setName("pickCompliance");
    compliance->compliance.setValue(_compliance);
    compliance->isCompliance.setValue(_isCompliance);
    interactionNode->addObject(compliance);
    compliance->rayleighStiffness.setValue(0.1);

    interactionNode->execute<simulation::InitVisitor>(sofa::core::ExecParams::defaultInstance());


}

template <class DataTypes>
void CompliantAttachPerformer<DataTypes>::execute()
{
    // update the distance mapping using the target position
    typename Point3dState::ReadVecCoord xmouse = mouseState->readPositions();
    distanceMapping->updateTarget(pickedParticleIndex,xmouse[0]);

    mouseMapping->apply(core::MechanicalParams::defaultInstance());
    mouseMapping->applyJ(core::MechanicalParams::defaultInstance());

    this->interactor->setMouseAttached(true);

//        cerr<<"CompliantAttachPerformer<DataTypes>::execute(), mouse position = " << xmouse[0] << endl;
}


template <class DataTypes>
void CompliantAttachPerformer<DataTypes>::configure(configurationsetting::MouseButtonSetting* setting)
{
    //Get back parameters from the MouseButtonSetting Component
    configurationsetting::CompliantAttachButtonSetting* s = dynamic_cast<configurationsetting::CompliantAttachButtonSetting*>(setting);
    if (s)
    {
        _compliance = s->compliance.getValue();
        _isCompliance = s->isCompliance.getValue();
        _arrowSize = s->arrowSize.getValue();
        _color = s->color.getValue();
    }
}


}
}
}

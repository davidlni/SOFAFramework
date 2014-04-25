
#include <sofa/helper/system/SetDirectory.h>

//Including Simulation
#include <sofa/simulation/common/Simulation.h>

//Including components for collision detection
#include <sofa/component/collision/DefaultPipeline.h>
#include <sofa/component/collision/DefaultContactManager.h>
#include <sofa/component/collision/DefaultCollisionGroupManager.h>

#include "../initContinuousCollision.h"
#include "../ContinuousDetection.h"
#include "../ContinuousIntersection.h"
#include "../RTriangleModel.h"


using namespace sofa::helper;
using sofa::helper::vector;
using namespace sofa::simulation;
using namespace sofa::core::objectmodel;
using namespace sofa::component;


Node::SPtr CreateRootWithCollisionPipeline(const std::string& responseType = "default")
{

  Node::SPtr root = sofa::simulation::getSimulation()->createNewGraph("root");

  //Components for collision management
  //------------------------------------
  //--> adding collision pipeline
  collision::DefaultPipeline::SPtr collisionPipeline = sofa::core::objectmodel::New<collision::DefaultPipeline>();
  collisionPipeline->setName("Collision Pipeline");
  root->addObject(collisionPipeline);

  //--> adding collision detection system
  collision::ContinuousDetection::SPtr detection = sofa::core::objectmodel::New<collision::ContinuousDetection>();
  detection->setName("Detection");
  root->addObject(detection);

  //--> adding component to detection intersection of elements
  collision::ContinuousIntersection::SPtr detectionProximity = sofa::core::objectmodel::New<collision::ContinuousIntersection>();
  detectionProximity->setName("Proximity");
  detectionProximity->setAlarmDistance(0.3);   //warning distance
  detectionProximity->setContactDistance(0.2); //min distance before setting a spring to create a repulsion
  root->addObject(detectionProximity);

  //--> adding contact manager
  collision::DefaultContactManager::SPtr contactManager = sofa::core::objectmodel::New<collision::DefaultContactManager>();
  contactManager->setName("Contact Manager");
  contactManager->setDefaultResponseType(responseType);
  root->addObject(contactManager);

  //--> adding component to handle groups of collision.
  collision::DefaultCollisionGroupManager::SPtr collisionGroupManager = sofa::core::objectmodel::New<collision::DefaultCollisionGroupManager>();
  collisionGroupManager->setName("Collision Group Manager");
  root->addObject(collisionGroupManager);

  return root;
}

int main(int ac, char **av)
{

}

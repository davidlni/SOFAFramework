
#include <sofa/helper/system/SetDirectory.h>

//Including Simulation
#include <sofa/simulation/common/Simulation.h>

#include <sofa/core/VecId.h>
//Including components for collision detection
#include <sofa/component/contextobject/Gravity.h>
#include <sofa/component/contextobject/CoordinateSystem.h>
#include <sofa/component/odesolver/EulerImplicitSolver.h>
#include <sofa/component/linearsolver/CGLinearSolver.h>

#include <sofa/component/collision/DefaultPipeline.h>
#include <sofa/component/collision/DefaultContactManager.h>
#include <sofa/component/collision/DefaultCollisionGroupManager.h>
#include <sofa/component/loader/MeshObjLoader.h>
#include <sofa/component/loader/MeshSTLLoader.h>
#include <sofa/helper/ArgumentParser.h>
#include <sofa/component/typedef/Sofa_typedef.h>
#include <sofa/helper/system/FileRepository.h>
#include <sofa/helper/system/glut.h>

#include <sofa/component/topology/RegularGridTopology.h>
#include <sofa/component/topology/MeshTopology.h>
#include <sofa/component/mapping/SubsetMapping.h>
#include <sofa/component/misc/STLExporter.h>

#include <sofa/simulation/tree/GNode.h>
#include <sofa/simulation/tree/TreeSimulation.h>
#include <sofa/gui/GUIManager.h>
#include <sofa/gui/Main.h>

#include "../initContinuousCollision.h"
#include "../ContinuousDetection.h"
#include "../ContinuousIntersection.h"
#include "../RTriangleModel.h"


using sofa::component::odesolver::EulerImplicitSolver;
typedef sofa::component::linearsolver::CGLinearSolver<sofa::component::linearsolver::GraphScatteredMatrix, sofa::component::linearsolver::GraphScatteredVector> CGLinearSolver;

sofa::simulation::Node::SPtr CreateRootWithCollisionPipeline(const std::string& responseType = "default")
{

  sofa::simulation::Node::SPtr root = sofa::simulation::getSimulation()->createNewGraph("root");

  //Components for collision management
  //------------------------------------
  //--> adding collision pipeline
  sofa::component::collision::DefaultPipeline::SPtr collisionPipeline = sofa::core::objectmodel::New<sofa::component::collision::DefaultPipeline>();
  collisionPipeline->setName("Collision Pipeline");
  root->addObject(collisionPipeline);

  //--> adding collision detection system
  sofa::component::collision::ContinuousDetection::SPtr detection = sofa::core::objectmodel::New<sofa::component::collision::ContinuousDetection>();
  detection->setName("Detection");
  root->addObject(detection);

  //--> adding component to detection intersection of elements
  sofa::component::collision::ContinuousIntersection::SPtr detectionProximity = sofa::core::objectmodel::New<sofa::component::collision::ContinuousIntersection>();
  detectionProximity->setName("Proximity");
  detectionProximity->setAlarmDistance(0.3);   //warning distance
  detectionProximity->setContactDistance(0.2); //min distance before setting a spring to create a repulsion
  root->addObject(detectionProximity);

  //--> adding contact manager
  sofa::component::collision::DefaultContactManager::SPtr contactManager = sofa::core::objectmodel::New<sofa::component::collision::DefaultContactManager>();
  contactManager->setName("Contact Manager");
  contactManager->setDefaultResponseType(responseType);
  root->addObject(contactManager);

  //--> adding component to handle groups of collision.
  sofa::component::collision::DefaultCollisionGroupManager::SPtr collisionGroupManager = sofa::core::objectmodel::New<sofa::component::collision::DefaultCollisionGroupManager>();
  collisionGroupManager->setName("Collision Group Manager");
  root->addObject(collisionGroupManager);

  return root;
}

MechanicalObject3::SPtr createFEM(sofa::simulation::Node::SPtr node)
{
  sofa::simulation::Node::SPtr FEMNode = node->createChild("FEMNode");

  // Tetrahedron degrees of freedom
  sofa::component::topology::RegularGridTopology::SPtr grid = sofa::core::objectmodel::New<sofa::component::topology::RegularGridTopology>();
  grid->setNumVertices(2,2,2);
  grid->setPos(-3.5,3.5,-3.5,3.5,-3.5,3.5);
  FEMNode->addObject(grid);

  MechanicalObject3::SPtr DOF = sofa::core::objectmodel::New<MechanicalObject3>();
  FEMNode->addObject(DOF);

  UniformMass3::SPtr mass = sofa::core::objectmodel::New<UniformMass3>();
  mass->mass.setValue( 0.25 );
  FEMNode->addObject(mass);

  TetrahedronFEMForceField3::SPtr tetraFEMFF = sofa::core::objectmodel::New<TetrahedronFEMForceField3>();
  tetraFEMFF->setName("FEM");
  tetraFEMFF->setComputeGlobalMatrix(false);
  tetraFEMFF->setMethod("large");
  tetraFEMFF->setPoissonRatio(0.3);
  tetraFEMFF->setYoungModulus(50);
  FEMNode->addObject(tetraFEMFF);

  return DOF;
}

MechanicalObject3::SPtr createOneTetraFEM(sofa::simulation::Node::SPtr node)
{
  sofa::simulation::Node::SPtr FEMNode = node->createChild("FEMNode");

  // Tetrahedron degrees of freedom
  MechanicalObject3::SPtr DOF = sofa::core::objectmodel::New<MechanicalObject3>();
  const Deriv3 rotation(30,0,24);
  FEMNode->addObject(DOF);
  DOF->resize(4);
  DOF->setName("VolumeDOF");
  DOF->setRotation(rotation[0],rotation[1],rotation[2]);
  //get write access to the position vector of mechanical object DOF
  sofa::helper::WriteAccessor<Data<VecCoord3> > x = *DOF->write(sofa::core::VecId::position());
  x[0] = Coord3(0,20,0);
  x[1] = Coord3(10,10,0);
  x[2] = Coord3(-5,10,10);
  x[3] = Coord3(-5,10,-10);
//   DOF->showObject.setValue(true);
//   DOF->showObjectScale.setValue(10.);

  UniformMass3::SPtr mass = sofa::core::objectmodel::New<UniformMass3>();
  mass->mass.setValue( 0.25 );
  FEMNode->addObject(mass);

  // Tetrahedron topology
  sofa::component::topology::MeshTopology::SPtr topology = sofa::core::objectmodel::New<sofa::component::topology::MeshTopology>();
  topology->setName("VolumeTopology");
  topology->addTetra(0,1,2,3);
  FEMNode->addObject( topology );

  TetrahedronFEMForceField3::SPtr tetraFEMFF = sofa::core::objectmodel::New<TetrahedronFEMForceField3>();
  tetraFEMFF->setName("FEM");
  tetraFEMFF->setComputeGlobalMatrix(false);
  tetraFEMFF->setMethod("large");
  tetraFEMFF->setPoissonRatio(0.45);
  tetraFEMFF->setYoungModulus(1);
  FEMNode->addObject(tetraFEMFF);

  return DOF;
}

MechanicalObject3::SPtr createFloor(sofa::simulation::Node::SPtr node)
{
  sofa::simulation::Node::SPtr FloorNode = node->createChild("FloorNode");

  const Deriv3 translation(0,-28,0);
  const Deriv3 rotation(0,0,24);
  sofa::component::loader::MeshSTLLoader::SPtr loader_surf = sofa::core::objectmodel::New<sofa::component::loader::MeshSTLLoader>();
  loader_surf->setName("loader");
  loader_surf->setFilename(sofa::helper::system::DataRepository.getFile("mesh/floor_tri.stl"));
  loader_surf->load();
  FloorNode->addObject(loader_surf);

  sofa::component::topology::MeshTopology::SPtr floorMesh= sofa::core::objectmodel::New<sofa::component::topology::MeshTopology>();
  floorMesh->setSrc("", loader_surf.get());
  FloorNode->addObject(floorMesh);

  MechanicalObject3::SPtr DOF = sofa::core::objectmodel::New<MechanicalObject3>();
  DOF->setName("Floor Object");
//   DOF->setTranslation(translation[0],translation[1],translation[2]);
//   DOF->setRotation(rotation[0],rotation[1],rotation[2]);
  FloorNode->addObject(DOF);
  //
  sofa::component::collision::RTriangleModel::SPtr triangle = sofa::core::objectmodel::New<sofa::component::collision::RTriangleModel>();
  triangle->setName("FloorTriangleCollision");
  triangle->setSelfCollision(true);
  triangle->setSimulated(false);
  triangle->setMoving(false);
  FloorNode->addObject(triangle);

  sofa::component::misc::STLExporter::SPtr exporter = sofa::core::objectmodel::New<sofa::component::misc::STLExporter>();
  exporter->stlFilename.setValue("floor");
  exporter->exportEveryNbSteps.setValue(1);
  exporter->m_fileFormat.setValue(0);
//   FloorNode->addObject(exporter);

  return DOF;
}

MechanicalObject3::SPtr createOneTetCollision(sofa::simulation::Node::SPtr node)
{
  sofa::simulation::Node::SPtr collisionNode = node->createChild("collisionNode");
  // Tetrahedron degrees of freedom

  const Deriv3 rotation(30,0,24);
  MechanicalObject3::SPtr DOF = sofa::core::objectmodel::New<MechanicalObject3>();
  collisionNode->addObject(DOF);
  DOF->resize(4);
  DOF->setName("SurfaceDOF");
  DOF->setRotation(rotation[0],rotation[1],rotation[2]);

  //get write access to the position vector of mechanical object DOF
  sofa::helper::WriteAccessor<Data<VecCoord3> > x = *DOF->write(sofa::core::VecId::position());
  x[0] = Coord3(0,20,0);
  x[1] = Coord3(10,10,0);
  x[2] = Coord3(-5,10,10);
  x[3] = Coord3(-5,10,-10);
//   DOF->showObject.setValue(true);
//   DOF->showObjectScale.setValue(10.);

  // Tetrahedron topology
  sofa::component::topology::MeshTopology::SPtr topology = sofa::core::objectmodel::New<sofa::component::topology::MeshTopology>();
  topology->setName("SurfaceTopology");
  topology->seqPoints.setParent(&DOF->x);
  topology->addTriangle(0,1,2);
  topology->addTriangle(1,2,3);
  topology->addTriangle(2,3,0);
  topology->addTriangle(3,0,1);
  collisionNode->addObject( topology );

  //
  sofa::component::collision::RTriangleModel::SPtr triangle = sofa::core::objectmodel::New<sofa::component::collision::RTriangleModel>();
  triangle->setName("TriangleCollision");
  triangle->setSelfCollision(true);
  collisionNode->addObject(triangle);

  sofa::component::misc::STLExporter::SPtr exporter = sofa::core::objectmodel::New<sofa::component::misc::STLExporter>();
//   collisionNode->addObject(exporter);
  exporter->stlFilename.setValue("tet");
  exporter->exportEveryNbSteps.setValue(1);
  exporter->m_fileFormat.setValue(0);
  return DOF;
}

MechanicalObject3::SPtr createCollision(sofa::simulation::Node::SPtr node)
{
  sofa::simulation::Node::SPtr collisionNode = node->createChild("collisionNode");

  sofa::component::loader::MeshObjLoader::SPtr loader_surf = sofa::core::objectmodel::New<sofa::component::loader::MeshObjLoader>();
  loader_surf->setName("loader");
  loader_surf->setFilename(sofa::helper::system::DataRepository.getFile("mesh/smCube125.obj"));
  loader_surf->load();
  collisionNode->addObject(loader_surf);

  sofa::component::topology::MeshTopology::SPtr cubeMesh= sofa::core::objectmodel::New<sofa::component::topology::MeshTopology>();
  cubeMesh->setSrc("", loader_surf.get());
  collisionNode->addObject(cubeMesh);

  MechanicalObject3::SPtr DOF = sofa::core::objectmodel::New<MechanicalObject3>();
  DOF->setName("Collision Object ");
  collisionNode->addObject(DOF);
  //
  sofa::component::collision::RTriangleModel::SPtr triangle = sofa::core::objectmodel::New<sofa::component::collision::RTriangleModel>();
  triangle->setName("TriangleCollision");
  triangle->setSelfCollision(true);
  collisionNode->addObject(triangle);

  return DOF;
}

int main(int ac, char **av)
{
  glutInit(&ac,av);
  sofa::helper::parse("This is a SOFA application. Here are the command line arguments")(ac,av);
  sofa::simulation::setSimulation(new sofa::simulation::tree::TreeSimulation());

  sofa::gui::initMain();
  sofa::gui::GUIManager::Init(av[0]);

  // The graph root node
  sofa::simulation::Node::SPtr root = CreateRootWithCollisionPipeline();
  root->setGravity(sofa::defaulttype::Vector3(0,-10,0) );
  // One solver for all the graph

  EulerImplicitSolver::SPtr solver = sofa::core::objectmodel::New<EulerImplicitSolver>();
  solver->setName("solver");
  solver->f_printLog.setValue(false);
  root->addObject(solver);

  CGLinearSolver::SPtr linearSolver = sofa::core::objectmodel::New<CGLinearSolver>();
  linearSolver->setName("linearSolver");
  root->addObject(linearSolver);

  MechanicalObject3::SPtr femDOF = createOneTetraFEM(root);
  MechanicalObject3::SPtr collisionDOF = createOneTetCollision(root);
//
  BarycentricMapping3_to_3::SPtr collisionMapping = sofa::core::objectmodel::New<BarycentricMapping3_to_3>();
  collisionMapping->setModels(femDOF.get(),collisionDOF.get());
  collisionDOF->getContext()->addObject(collisionMapping);

  createFloor(root);

  root->setAnimate(false);

  sofa::simulation::getSimulation()->init(root.get());

  //=======================================
  // Run the main loop
  sofa::gui::GUIManager::MainLoop(root);

}

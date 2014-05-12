/******************************************************************************
 *       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
 *                (c) 2006-2011 INRIA, USTL, UJF, CNRS, MGH                    *
 *                                                                             *
 * This program is free software; you can redistribute it and/or modify it     *
 * under the terms of the GNU General Public License as published by the Free  *
 * Software Foundation; either version 2 of the License, or (at your option)   *
 * any later version.                                                          *
 *                                                                             *
 * This program is distributed in the hope that it will be useful, but WITHOUT *
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for    *
 * more details.                                                               *
 *                                                                             *
 * You should have received a copy of the GNU General Public License along     *
 * with this program; if not, write to the Free Software Foundation, Inc., 51  *
 * Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.                   *
 *******************************************************************************
 *                            SOFA :: Applications                             *
 *                                                                             *
 * Authors: The SOFA Team and external contributors (see Authors.txt)          *
 *                                                                             *
 * Contact information: contact@sofa-framework.org                             *
 ******************************************************************************/

#include <plugins/SofaTest/Sofa_test.h>
#include <sofa/helper/system/SetDirectory.h>
#include <sofa/helper/system/FileRepository.h>
#include <sofa/component/init.h>
#include <sofa/core/ExecParams.h>

//Including Simulation
#include <sofa/simulation/common/Simulation.h>
#include <sofa/simulation/graph/DAGSimulation.h>
#include <sofa/simulation/common/Node.h>

// Including component
#include <sofa/component/container/MechanicalObject.h>

namespace sofa {

  using namespace component;
  using namespace defaulttype;

  // Polytope test: Test BVH creation using the polytope model

  template <typename _DataTypes>
  struct PolytopeModel_test : public Sofa_test<typename _DataTypes::Real>
  {
    typedef _DataTypes DataTypes;
    typedef typename DataTypes::CPos CPos;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef container::MechanicalObject<DataTypes> MechanicalObject;

    /// Root of the scene graph
    simulation::Node::SPtr root;

    /// Simulation
    simulation::Simulation* simulation;

    /// Create the context for the scene
    void SetUp()
    {
      // Init simulation
      sofa::component::init();
      sofa::simulation::setSimulation(simulation = new sofa::simulation::graph::DAGSimulation());

      root = sofa::simulation::getSimulation()->createNewGraph("root");
    }

    /// Load the scene to test SmallCorotationalPatchTest.scn
    void loadScene(std::string sceneName)
    {
      // Load the scene from the xml file
      std::string fileName = std::string(CCD_TEST_SCENES_DIR) + "/" + sceneName;
      root = sofa::core::objectmodel::SPtr_dynamic_cast<sofa::simulation::Node>( sofa::simulation::getSimulation()->load(fileName.c_str()));
    }

    /// Unload the scene
    void TearDown()
    {
      if (root!=NULL)
        sofa::simulation::getSimulation()->unload(root);
    }

  };

  // Define the list of DataTypes to instantiate
  using testing::Types;
  typedef testing::Types<Vec3Types> DataTypes; // the types to instantiate.

  // Test suite for all the instantiations
  TYPED_TEST_CASE(PolytopeModel_test, DataTypes);


} // namespace sofa

/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 INRIA, USTL, UJF, CNRS, MGH                    *
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



// File automatically generated by "generateTypedef"


#ifndef SOFA_TYPEDEF_Mass_double_H
#define SOFA_TYPEDEF_Mass_double_H

//Default files containing the declaration of the vector type
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/defaulttype/Mat.h>




#include <BaseMechanics/DiagonalMass.h>
#include <SimpleFem/forcefield/HexahedralFEMForceFieldAndMass.h>
#include <NonUniformFem/forcefield/HexahedronCompositeFEMForceFieldAndMass.h>
#include <SimpleFem/forcefield/HexahedronFEMForceFieldAndMass.h>
#include <MiscForceField/MatrixMass.h>
#include <MiscForceField/MeshMatrixMass.h>
#include <NonUniformFem/forcefield/NonUniformHexahedralFEMForceFieldAndMass.h>
#include <NonUniformFem/forcefield/NonUniformHexahedronFEMForceFieldAndMass.h>
#include <BaseMechanics/UniformMass.h>



//---------------------------------------------------------------------------------------------
//Typedef for DiagonalMass
typedef sofa::component::mass::DiagonalMass<sofa::defaulttype::StdRigidTypes<3, double>, sofa::defaulttype::RigidMass<3, double> > DiagonalMassRigid3d;
typedef sofa::component::mass::DiagonalMass<sofa::defaulttype::StdRigidTypes<2, double>, sofa::defaulttype::RigidMass<2, double> > DiagonalMassRigid2d;
typedef sofa::component::mass::DiagonalMass<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<1, double>, sofa::defaulttype::Vec<1, double>, double>, double> DiagonalMass1d;
typedef sofa::component::mass::DiagonalMass<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<2, double>, sofa::defaulttype::Vec<2, double>, double>, double> DiagonalMass2d;
typedef sofa::component::mass::DiagonalMass<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double>, double> DiagonalMass3d;



//---------------------------------------------------------------------------------------------
//Typedef for HexahedralFEMForceFieldAndMass
typedef sofa::component::forcefield::HexahedralFEMForceFieldAndMass<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double> > HexahedralFEMForceFieldAndMass3d;



//---------------------------------------------------------------------------------------------
//Typedef for HexahedronCompositeFEMForceFieldAndMass
typedef sofa::component::forcefield::HexahedronCompositeFEMForceFieldAndMass<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double> > HexahedronCompositeFEMForceFieldAndMass3d;



//---------------------------------------------------------------------------------------------
//Typedef for HexahedronFEMForceFieldAndMass
typedef sofa::component::forcefield::HexahedronFEMForceFieldAndMass<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double> > HexahedronFEMForceFieldAndMass3d;



//---------------------------------------------------------------------------------------------
//Typedef for MatrixMass
typedef sofa::component::mass::MatrixMass<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<1, double>, sofa::defaulttype::Vec<1, double>, double>, sofa::defaulttype::Mat<1, 1, double> > MatrixMass1d;
typedef sofa::component::mass::MatrixMass<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<2, double>, sofa::defaulttype::Vec<2, double>, double>, sofa::defaulttype::Mat<2, 2, double> > MatrixMass2d;
typedef sofa::component::mass::MatrixMass<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double>, sofa::defaulttype::Mat<3, 3, double> > MatrixMass3d;



//---------------------------------------------------------------------------------------------
//Typedef for MeshMatrixMass
typedef sofa::component::mass::MeshMatrixMass<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<1, double>, sofa::defaulttype::Vec<1, double>, double>, double> MeshMatrixMass1d;
typedef sofa::component::mass::MeshMatrixMass<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<2, double>, sofa::defaulttype::Vec<2, double>, double>, double> MeshMatrixMass2d;
typedef sofa::component::mass::MeshMatrixMass<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double>, double> MeshMatrixMass3d;



//---------------------------------------------------------------------------------------------
//Typedef for NonUniformHexahedralFEMForceFieldAndMass
typedef sofa::component::forcefield::NonUniformHexahedralFEMForceFieldAndMass<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double> > NonUniformHexahedralFEMForceFieldAndMass3d;



//---------------------------------------------------------------------------------------------
//Typedef for NonUniformHexahedronFEMForceFieldAndMass
typedef sofa::component::forcefield::NonUniformHexahedronFEMForceFieldAndMass<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double> > NonUniformHexahedronFEMForceFieldAndMass3d;



//---------------------------------------------------------------------------------------------
//Typedef for UniformMass
typedef sofa::component::mass::UniformMass<sofa::defaulttype::StdRigidTypes<3, double>, sofa::defaulttype::RigidMass<3, double> > UniformMassRigid3d;
typedef sofa::component::mass::UniformMass<sofa::defaulttype::StdRigidTypes<2, double>, sofa::defaulttype::RigidMass<2, double> > UniformMassRigid2d;
typedef sofa::component::mass::UniformMass<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<1, double>, sofa::defaulttype::Vec<1, double>, double>, double> UniformMass1d;
typedef sofa::component::mass::UniformMass<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<2, double>, sofa::defaulttype::Vec<2, double>, double>, double> UniformMass2d;
typedef sofa::component::mass::UniformMass<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double>, double> UniformMass3d;
typedef sofa::component::mass::UniformMass<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<6, double>, sofa::defaulttype::Vec<6, double>, double>, double> UniformMass6d;





#ifndef SOFA_FLOAT
typedef DiagonalMassRigid3d DiagonalMassRigid3;
typedef DiagonalMassRigid2d DiagonalMassRigid2;
typedef DiagonalMass1d DiagonalMass1;
typedef DiagonalMass2d DiagonalMass2;
typedef DiagonalMass3d DiagonalMass3;
typedef HexahedralFEMForceFieldAndMass3d HexahedralFEMForceFieldAndMass3;
typedef HexahedronCompositeFEMForceFieldAndMass3d HexahedronCompositeFEMForceFieldAndMass3;
typedef HexahedronFEMForceFieldAndMass3d HexahedronFEMForceFieldAndMass3;
typedef MatrixMass1d MatrixMass1;
typedef MatrixMass2d MatrixMass2;
typedef MatrixMass3d MatrixMass3;
typedef MeshMatrixMass1d MeshMatrixMass1;
typedef MeshMatrixMass2d MeshMatrixMass2;
typedef MeshMatrixMass3d MeshMatrixMass3;
typedef NonUniformHexahedralFEMForceFieldAndMass3d NonUniformHexahedralFEMForceFieldAndMass3;
typedef NonUniformHexahedronFEMForceFieldAndMass3d NonUniformHexahedronFEMForceFieldAndMass3;
typedef UniformMassRigid3d UniformMassRigid3;
typedef UniformMassRigid2d UniformMassRigid2;
typedef UniformMass1d UniformMass1;
typedef UniformMass2d UniformMass2;
typedef UniformMass3d UniformMass3;
typedef UniformMass6d UniformMass6;
#endif

#endif

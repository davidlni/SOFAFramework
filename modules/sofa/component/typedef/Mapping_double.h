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


#ifndef SOFA_TYPEDEF_Mapping_double_H
#define SOFA_TYPEDEF_Mapping_double_H

//Default files containing the declaration of the vector type
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/defaulttype/Mat.h>


//Default files needed to create a Mapping
#include <sofa/core/State.h>
#include <sofa/core/Mapping.h>


#include <sofa/component/mapping/ArticulatedSystemMapping.h>
#include <sofa/component/mapping/BarycentricMapping.h>
#include <sofa/component/mapping/BeamLinearMapping.h>
#include <sofa/component/mapping/CenterOfMassMapping.h>
#include <sofa/component/mapping/CenterOfMassMulti2Mapping.h>
#include <sofa/component/mapping/CenterOfMassMultiMapping.h>
#include <sofa/component/mapping/CenterPointMechanicalMapping.h>
#include <sofa/component/mapping/CurveMapping.h>
#include <sofa/component/mapping/ExternalInterpolationMapping.h>
#include <sofa/component/mapping/HexahedronCompositeFEMMapping.h>
#include <sofa/component/mapping/IdentityMapping.h>
#include <sofa/component/mapping/ImplicitSurfaceMapping.h>
#include <sofa/component/mapping/LineSetSkinningMapping.h>
#include <sofa/component/mapping/Mesh2PointMechanicalMapping.h>
#include <sofa/component/mapping/RigidMapping.h>
#include <sofa/component/mapping/RigidRigidMapping.h>
#include <sofa/component/mapping/SPHFluidSurfaceMapping.h>
#include <sofa/component/mapping/SimpleTesselatedTetraMechanicalMapping.h>
#include <sofa/component/mapping/SkinningMapping.h>
#include <sofa/component/mapping/SubsetMapping.h>
#include <sofa/component/mapping/SubsetMultiMapping.h>
#include <sofa/component/mapping/TubularMapping.h>



//---------------------------------------------------------------------------------------------
//Typedef for ArticulatedSystemMapping
typedef sofa::component::mapping::ArticulatedSystemMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<1, double>, sofa::defaulttype::Vec<1, double>, double>, sofa::defaulttype::StdRigidTypes<3, double>, sofa::defaulttype::StdRigidTypes<3, double> > ArticulatedSystemMapping1d_Rigid3d_to_Rigid3d;



//---------------------------------------------------------------------------------------------
//Typedef for BarycentricMapping
typedef sofa::component::mapping::BarycentricMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double>, sofa::defaulttype::ExtVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > BarycentricMapping3d_to_Ext3f;
typedef sofa::component::mapping::BarycentricMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double> > BarycentricMapping3d_to_3d;
typedef sofa::component::mapping::BarycentricMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double>, sofa::defaulttype::StdRigidTypes<3, double> > BarycentricMapping3d_to_Rigid3d;



//---------------------------------------------------------------------------------------------
//Typedef for BeamLinearMapping
typedef sofa::component::mapping::BeamLinearMapping<sofa::defaulttype::StdRigidTypes<3, double>, sofa::defaulttype::ExtVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > BeamLinearMappingRigid3d_to_Ext3f;
typedef sofa::component::mapping::BeamLinearMapping<sofa::defaulttype::StdRigidTypes<3, double>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double> > BeamLinearMappingRigid3d_to_3d;



//---------------------------------------------------------------------------------------------
//Typedef for CenterOfMassMapping
typedef sofa::component::mapping::CenterOfMassMapping<sofa::defaulttype::StdRigidTypes<3, double>, sofa::defaulttype::ExtVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double> > CenterOfMassMappingRigid3d_to_Ext3d;
typedef sofa::component::mapping::CenterOfMassMapping<sofa::defaulttype::StdRigidTypes<3, double>, sofa::defaulttype::ExtVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > CenterOfMassMappingRigid3d_to_Ext3f;
typedef sofa::component::mapping::CenterOfMassMapping<sofa::defaulttype::StdRigidTypes<3, double>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double> > CenterOfMassMappingRigid3d_to_3d;
typedef sofa::component::mapping::CenterOfMassMapping<sofa::defaulttype::StdRigidTypes<2, double>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<2, double>, sofa::defaulttype::Vec<2, double>, double> > CenterOfMassMappingRigid2d_to_2d;



//---------------------------------------------------------------------------------------------
//Typedef for CenterOfMassMulti2Mapping
typedef sofa::component::mapping::CenterOfMassMulti2Mapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double>, sofa::defaulttype::StdRigidTypes<3, double>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double> > CenterOfMassMulti2Mapping3d_Rigid3d_to_3d;



//---------------------------------------------------------------------------------------------
//Typedef for CenterOfMassMultiMapping
typedef sofa::component::mapping::CenterOfMassMultiMapping<sofa::defaulttype::StdRigidTypes<3, double>, sofa::defaulttype::StdRigidTypes<3, double> > CenterOfMassMultiMappingRigid3d_to_Rigid3d;
typedef sofa::component::mapping::CenterOfMassMultiMapping<sofa::defaulttype::StdRigidTypes<3, double>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double> > CenterOfMassMultiMappingRigid3d_to_3d;
typedef sofa::component::mapping::CenterOfMassMultiMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double> > CenterOfMassMultiMapping3d_to_3d;



//---------------------------------------------------------------------------------------------
//Typedef for CenterPointMechanicalMapping
typedef sofa::component::mapping::CenterPointMechanicalMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double>, sofa::defaulttype::ExtVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > CenterPointMechanicalMapping3d_to_Ext3f;
typedef sofa::component::mapping::CenterPointMechanicalMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double> > CenterPointMechanicalMapping3d_to_3d;



//---------------------------------------------------------------------------------------------
//Typedef for CurveMapping
typedef sofa::component::mapping::CurveMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double>, sofa::defaulttype::StdRigidTypes<3, double> > CurveMapping3d_to_Rigid3d;



//---------------------------------------------------------------------------------------------
//Typedef for ExternalInterpolationMapping
typedef sofa::component::mapping::ExternalInterpolationMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<1, double>, sofa::defaulttype::Vec<1, double>, double>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<1, double>, sofa::defaulttype::Vec<1, double>, double> > ExternalInterpolationMapping1d_to_1d;
typedef sofa::component::mapping::ExternalInterpolationMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<2, double>, sofa::defaulttype::Vec<2, double>, double>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<2, double>, sofa::defaulttype::Vec<2, double>, double> > ExternalInterpolationMapping2d_to_2d;
typedef sofa::component::mapping::ExternalInterpolationMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double>, sofa::defaulttype::ExtVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > ExternalInterpolationMapping3d_to_Ext3f;
typedef sofa::component::mapping::ExternalInterpolationMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double> > ExternalInterpolationMapping3d_to_3d;



//---------------------------------------------------------------------------------------------
//Typedef for HexahedronCompositeFEMMapping
typedef sofa::component::mapping::HexahedronCompositeFEMMapping<sofa::core::Mapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double>, sofa::defaulttype::ExtVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > > HexahedronCompositeFEMMapping3d_to_Ext3f;
typedef sofa::component::mapping::HexahedronCompositeFEMMapping<sofa::core::Mapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double> > > HexahedronCompositeFEMMapping3d_to_3d;



//---------------------------------------------------------------------------------------------
//Typedef for IdentityMapping
typedef sofa::component::mapping::IdentityMapping<sofa::defaulttype::StdRigidTypes<3, double>, sofa::defaulttype::ExtVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > IdentityMappingRigid3d_to_Ext3f;
typedef sofa::component::mapping::IdentityMapping<sofa::defaulttype::StdRigidTypes<3, double>, sofa::defaulttype::StdRigidTypes<3, double> > IdentityMappingRigid3d_to_Rigid3d;
typedef sofa::component::mapping::IdentityMapping<sofa::defaulttype::StdRigidTypes<3, double>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double> > IdentityMappingRigid3d_to_3d;
typedef sofa::component::mapping::IdentityMapping<sofa::defaulttype::StdRigidTypes<2, double>, sofa::defaulttype::StdRigidTypes<2, double> > IdentityMappingRigid2d_to_Rigid2d;
typedef sofa::component::mapping::IdentityMapping<sofa::defaulttype::StdRigidTypes<2, double>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<2, double>, sofa::defaulttype::Vec<2, double>, double> > IdentityMappingRigid2d_to_2d;
typedef sofa::component::mapping::IdentityMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<1, double>, sofa::defaulttype::Vec<1, double>, double>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<1, double>, sofa::defaulttype::Vec<1, double>, double> > IdentityMapping1d_to_1d;
typedef sofa::component::mapping::IdentityMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<2, double>, sofa::defaulttype::Vec<2, double>, double>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<2, double>, sofa::defaulttype::Vec<2, double>, double> > IdentityMapping2d_to_2d;
typedef sofa::component::mapping::IdentityMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double>, sofa::defaulttype::ExtVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > IdentityMapping3d_to_Ext3f;
typedef sofa::component::mapping::IdentityMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double> > IdentityMapping3d_to_3d;
typedef sofa::component::mapping::IdentityMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<6, double>, sofa::defaulttype::Vec<6, double>, double>, sofa::defaulttype::ExtVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > IdentityMapping6d_to_Ext3f;
typedef sofa::component::mapping::IdentityMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<6, double>, sofa::defaulttype::Vec<6, double>, double>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double> > IdentityMapping6d_to_3d;
typedef sofa::component::mapping::IdentityMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<6, double>, sofa::defaulttype::Vec<6, double>, double>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<6, double>, sofa::defaulttype::Vec<6, double>, double> > IdentityMapping6d_to_6d;



//---------------------------------------------------------------------------------------------
//Typedef for ImplicitSurfaceMapping
typedef sofa::component::mapping::ImplicitSurfaceMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double>, sofa::defaulttype::ExtVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > ImplicitSurfaceMapping3d_to_Ext3f;
typedef sofa::component::mapping::ImplicitSurfaceMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double> > ImplicitSurfaceMapping3d_to_3d;



//---------------------------------------------------------------------------------------------
//Typedef for LineSetSkinningMapping
typedef sofa::component::mapping::LineSetSkinningMapping<sofa::defaulttype::StdRigidTypes<3, double>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double> > LineSetSkinningMappingRigid3d_to_3d;



//---------------------------------------------------------------------------------------------
//Typedef for Mesh2PointMechanicalMapping
typedef sofa::component::mapping::Mesh2PointMechanicalMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double>, sofa::defaulttype::ExtVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > Mesh2PointMechanicalMapping3d_to_Ext3f;
typedef sofa::component::mapping::Mesh2PointMechanicalMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double> > Mesh2PointMechanicalMapping3d_to_3d;



//---------------------------------------------------------------------------------------------
//Typedef for RigidMapping
typedef sofa::component::mapping::RigidMapping<sofa::defaulttype::StdRigidTypes<3, double>, sofa::defaulttype::ExtVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > RigidMappingRigid3d_to_Ext3f;
typedef sofa::component::mapping::RigidMapping<sofa::defaulttype::StdRigidTypes<3, double>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double> > RigidMappingRigid3d_to_3d;
typedef sofa::component::mapping::RigidMapping<sofa::defaulttype::StdRigidTypes<2, double>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<2, double>, sofa::defaulttype::Vec<2, double>, double> > RigidMappingRigid2d_to_2d;



//---------------------------------------------------------------------------------------------
//Typedef for RigidRigidMapping
typedef sofa::component::mapping::RigidRigidMapping<sofa::defaulttype::StdRigidTypes<3, double>, sofa::defaulttype::StdRigidTypes<3, double> > RigidRigidMappingRigid3d_to_Rigid3d;



//---------------------------------------------------------------------------------------------
//Typedef for SPHFluidSurfaceMapping
typedef sofa::component::mapping::SPHFluidSurfaceMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double>, sofa::defaulttype::ExtVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > SPHFluidSurfaceMapping3d_to_Ext3f;
typedef sofa::component::mapping::SPHFluidSurfaceMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double> > SPHFluidSurfaceMapping3d_to_3d;



//---------------------------------------------------------------------------------------------
//Typedef for SimpleTesselatedTetraMechanicalMapping
typedef sofa::component::mapping::SimpleTesselatedTetraMechanicalMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double>, sofa::defaulttype::ExtVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > SimpleTesselatedTetraMechanicalMapping3d_to_Ext3f;
typedef sofa::component::mapping::SimpleTesselatedTetraMechanicalMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double> > SimpleTesselatedTetraMechanicalMapping3d_to_3d;



//---------------------------------------------------------------------------------------------
//Typedef for SkinningMapping
typedef sofa::component::mapping::SkinningMapping<sofa::defaulttype::StdRigidTypes<3, double>, sofa::defaulttype::ExtVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > SkinningMappingRigid3d_to_Ext3f;
typedef sofa::component::mapping::SkinningMapping<sofa::defaulttype::StdRigidTypes<3, double>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double> > SkinningMappingRigid3d_to_3d;



//---------------------------------------------------------------------------------------------
//Typedef for SubsetMapping
typedef sofa::component::mapping::SubsetMapping<sofa::defaulttype::StdRigidTypes<3, double>, sofa::defaulttype::StdRigidTypes<3, double> > SubsetMappingRigid3d_to_Rigid3d;
typedef sofa::component::mapping::SubsetMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<1, double>, sofa::defaulttype::Vec<1, double>, double>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<1, double>, sofa::defaulttype::Vec<1, double>, double> > SubsetMapping1d_to_1d;
typedef sofa::component::mapping::SubsetMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double>, sofa::defaulttype::ExtVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > SubsetMapping3d_to_Ext3f;
typedef sofa::component::mapping::SubsetMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double> > SubsetMapping3d_to_3d;



//---------------------------------------------------------------------------------------------
//Typedef for SubsetMultiMapping
typedef sofa::component::mapping::SubsetMultiMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double> > SubsetMultiMapping3d_to_3d;



//---------------------------------------------------------------------------------------------
//Typedef for TubularMapping
typedef sofa::component::mapping::TubularMapping<sofa::defaulttype::StdRigidTypes<3, double>, sofa::defaulttype::ExtVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > TubularMappingRigid3d_to_Ext3f;
typedef sofa::component::mapping::TubularMapping<sofa::defaulttype::StdRigidTypes<3, double>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double> > TubularMappingRigid3d_to_3d;





#ifndef SOFA_FLOAT
typedef ArticulatedSystemMapping1d_Rigid3d_to_Rigid3d ArticulatedSystemMapping1_Rigid3_to_Rigid3;
typedef BarycentricMapping3d_to_Ext3f BarycentricMapping3_to_Ext3;
typedef BarycentricMapping3d_to_3d BarycentricMapping3_to_3;
typedef BarycentricMapping3d_to_Rigid3d BarycentricMapping3_to_Rigid3;
typedef BeamLinearMappingRigid3d_to_Ext3f BeamLinearMappingRigid3_to_Ext3;
typedef BeamLinearMappingRigid3d_to_3d BeamLinearMappingRigid3_to_3;
typedef CenterOfMassMappingRigid3d_to_Ext3d CenterOfMassMappingRigid3_to_Ext3d;
typedef CenterOfMassMappingRigid3d_to_Ext3f CenterOfMassMappingRigid3_to_Ext3;
typedef CenterOfMassMappingRigid3d_to_3d CenterOfMassMappingRigid3_to_3;
typedef CenterOfMassMappingRigid2d_to_2d CenterOfMassMappingRigid2_to_2;
typedef CenterOfMassMulti2Mapping3d_Rigid3d_to_3d CenterOfMassMulti2Mapping3_Rigid3_to_3;
typedef CenterOfMassMultiMappingRigid3d_to_Rigid3d CenterOfMassMultiMappingRigid3_to_Rigid3;
typedef CenterOfMassMultiMappingRigid3d_to_3d CenterOfMassMultiMappingRigid3_to_3;
typedef CenterOfMassMultiMapping3d_to_3d CenterOfMassMultiMapping3_to_3;
typedef CenterPointMechanicalMapping3d_to_Ext3f CenterPointMechanicalMapping3_to_Ext3;
typedef CenterPointMechanicalMapping3d_to_3d CenterPointMechanicalMapping3_to_3;
typedef CurveMapping3d_to_Rigid3d CurveMapping3_to_Rigid3;
typedef ExternalInterpolationMapping1d_to_1d ExternalInterpolationMapping1_to_1;
typedef ExternalInterpolationMapping2d_to_2d ExternalInterpolationMapping2_to_2;
typedef ExternalInterpolationMapping3d_to_Ext3f ExternalInterpolationMapping3_to_Ext3;
typedef ExternalInterpolationMapping3d_to_3d ExternalInterpolationMapping3_to_3;
typedef HexahedronCompositeFEMMapping3d_to_Ext3f HexahedronCompositeFEMMapping3_to_Ext3;
typedef HexahedronCompositeFEMMapping3d_to_3d HexahedronCompositeFEMMapping3_to_3;
typedef IdentityMappingRigid3d_to_Ext3f IdentityMappingRigid3_to_Ext3;
typedef IdentityMappingRigid3d_to_Rigid3d IdentityMappingRigid3_to_Rigid3;
typedef IdentityMappingRigid3d_to_3d IdentityMappingRigid3_to_3;
typedef IdentityMappingRigid2d_to_Rigid2d IdentityMappingRigid2_to_Rigid2;
typedef IdentityMappingRigid2d_to_2d IdentityMappingRigid2_to_2;
typedef IdentityMapping1d_to_1d IdentityMapping1_to_1;
typedef IdentityMapping2d_to_2d IdentityMapping2_to_2;
typedef IdentityMapping3d_to_Ext3f IdentityMapping3_to_Ext3;
typedef IdentityMapping3d_to_3d IdentityMapping3_to_3;
typedef IdentityMapping6d_to_Ext3f IdentityMapping6_to_Ext3;
typedef IdentityMapping6d_to_3d IdentityMapping6_to_3;
typedef IdentityMapping6d_to_6d IdentityMapping6_to_6;
typedef ImplicitSurfaceMapping3d_to_Ext3f ImplicitSurfaceMapping3_to_Ext3;
typedef ImplicitSurfaceMapping3d_to_3d ImplicitSurfaceMapping3_to_3;
typedef LineSetSkinningMappingRigid3d_to_3d LineSetSkinningMappingRigid3_to_3;
typedef Mesh2PointMechanicalMapping3d_to_Ext3f Mesh2PointMechanicalMapping3_to_Ext3;
typedef Mesh2PointMechanicalMapping3d_to_3d Mesh2PointMechanicalMapping3_to_3;
typedef RigidMappingRigid3d_to_Ext3f RigidMappingRigid3_to_Ext3;
typedef RigidMappingRigid3d_to_3d RigidMappingRigid3_to_3;
typedef RigidMappingRigid2d_to_2d RigidMappingRigid2_to_2;
typedef RigidRigidMappingRigid3d_to_Rigid3d RigidRigidMappingRigid3_to_Rigid3;
typedef SPHFluidSurfaceMapping3d_to_Ext3f SPHFluidSurfaceMapping3_to_Ext3;
typedef SPHFluidSurfaceMapping3d_to_3d SPHFluidSurfaceMapping3_to_3;
typedef SimpleTesselatedTetraMechanicalMapping3d_to_Ext3f SimpleTesselatedTetraMechanicalMapping3_to_Ext3;
typedef SimpleTesselatedTetraMechanicalMapping3d_to_3d SimpleTesselatedTetraMechanicalMapping3_to_3;
typedef SkinningMappingRigid3d_to_Ext3f SkinningMappingRigid3_to_Ext3;
typedef SkinningMappingRigid3d_to_3d SkinningMappingRigid3_to_3;
typedef SubsetMappingRigid3d_to_Rigid3d SubsetMappingRigid3_to_Rigid3;
typedef SubsetMapping1d_to_1d SubsetMapping1_to_1;
typedef SubsetMapping3d_to_Ext3f SubsetMapping3_to_Ext3;
typedef SubsetMapping3d_to_3d SubsetMapping3_to_3;
typedef SubsetMultiMapping3d_to_3d SubsetMultiMapping3_to_3;
typedef TubularMappingRigid3d_to_Ext3f TubularMappingRigid3_to_Ext3;
typedef TubularMappingRigid3d_to_3d TubularMappingRigid3_to_3;
#endif

#endif

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


#ifndef SOFA_TYPEDEF_Mapping_float_H
#define SOFA_TYPEDEF_Mapping_float_H

//Default files containing the declaration of the vector type
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/defaulttype/Mat.h>


//Default files needed to create a Mapping
#include <sofa/core/State.h>
#include <sofa/core/Mapping.h>


#include <Rigid/ArticulatedSystemMapping.h>
#include <BaseMechanics/BarycentricMapping.h>
#include <MiscMapping/BeamLinearMapping.h>
#include <MiscMapping/CenterOfMassMapping.h>
#include <MiscMapping/CenterOfMassMulti2Mapping.h>
#include <MiscMapping/CenterOfMassMultiMapping.h>
#include <MiscMapping/CenterPointMechanicalMapping.h>
#include <MiscMapping/CurveMapping.h>
#include <MiscMapping/ExternalInterpolationMapping.h>
#include <NonUniformFem/mapping/HexahedronCompositeFEMMapping.h>
#include <BaseMechanics/IdentityMapping.h>
#include <VolumetricData/mapping/ImplicitSurfaceMapping.h>
#include <Rigid/LineSetSkinningMapping.h>
#include <TopologyMapping/mapping/Mesh2PointMechanicalMapping.h>
#include <Rigid/RigidMapping.h>
#include <Rigid/RigidRigidMapping.h>
#include <SphFluid/SPHFluidSurfaceMapping.h>
#include <TopologyMapping/mapping/SimpleTesselatedTetraMechanicalMapping.h>
#include <Rigid/SkinningMapping.h>
#include <BaseMechanics/SubsetMapping.h>
#include <MiscMapping/SubsetMultiMapping.h>
#include <MiscMapping/TubularMapping.h>



//---------------------------------------------------------------------------------------------
//Typedef for ArticulatedSystemMapping
typedef sofa::component::mapping::ArticulatedSystemMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<1, float>, sofa::defaulttype::Vec<1, float>, float>, sofa::defaulttype::StdRigidTypes<3, float>, sofa::defaulttype::StdRigidTypes<3, float> > ArticulatedSystemMapping1f_Rigid3f_to_Rigid3f;



//---------------------------------------------------------------------------------------------
//Typedef for BarycentricMapping
typedef sofa::component::mapping::BarycentricMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float>, sofa::defaulttype::ExtVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > BarycentricMapping3f_to_Ext3f;
typedef sofa::component::mapping::BarycentricMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > BarycentricMapping3f_to_3f;
typedef sofa::component::mapping::BarycentricMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float>, sofa::defaulttype::StdRigidTypes<3, float> > BarycentricMapping3f_to_Rigid3f;



//---------------------------------------------------------------------------------------------
//Typedef for BeamLinearMapping
typedef sofa::component::mapping::BeamLinearMapping<sofa::defaulttype::StdRigidTypes<3, float>, sofa::defaulttype::ExtVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > BeamLinearMappingRigid3f_to_Ext3f;
typedef sofa::component::mapping::BeamLinearMapping<sofa::defaulttype::StdRigidTypes<3, float>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > BeamLinearMappingRigid3f_to_3f;



//---------------------------------------------------------------------------------------------
//Typedef for CenterOfMassMapping
typedef sofa::component::mapping::CenterOfMassMapping<sofa::defaulttype::StdRigidTypes<2, float>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<2, float>, sofa::defaulttype::Vec<2, float>, float> > CenterOfMassMappingRigid2f_to_2f;
typedef sofa::component::mapping::CenterOfMassMapping<sofa::defaulttype::StdRigidTypes<3, float>, sofa::defaulttype::ExtVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > CenterOfMassMappingRigid3f_to_Ext3f;
typedef sofa::component::mapping::CenterOfMassMapping<sofa::defaulttype::StdRigidTypes<3, float>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > CenterOfMassMappingRigid3f_to_3f;



//---------------------------------------------------------------------------------------------
//Typedef for CenterOfMassMulti2Mapping
typedef sofa::component::mapping::CenterOfMassMulti2Mapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float>, sofa::defaulttype::StdRigidTypes<3, float>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > CenterOfMassMulti2Mapping3f_Rigid3f_to_3f;



//---------------------------------------------------------------------------------------------
//Typedef for CenterOfMassMultiMapping
typedef sofa::component::mapping::CenterOfMassMultiMapping<sofa::defaulttype::StdRigidTypes<3, float>, sofa::defaulttype::StdRigidTypes<3, float> > CenterOfMassMultiMappingRigid3f_to_Rigid3f;
typedef sofa::component::mapping::CenterOfMassMultiMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > CenterOfMassMultiMapping3f_to_3f;



//---------------------------------------------------------------------------------------------
//Typedef for CenterPointMechanicalMapping
typedef sofa::component::mapping::CenterPointMechanicalMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float>, sofa::defaulttype::ExtVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > CenterPointMechanicalMapping3f_to_Ext3f;
typedef sofa::component::mapping::CenterPointMechanicalMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > CenterPointMechanicalMapping3f_to_3f;



//---------------------------------------------------------------------------------------------
//Typedef for CurveMapping
typedef sofa::component::mapping::CurveMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float>, sofa::defaulttype::StdRigidTypes<3, float> > CurveMapping3f_to_Rigid3f;



//---------------------------------------------------------------------------------------------
//Typedef for ExternalInterpolationMapping
typedef sofa::component::mapping::ExternalInterpolationMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<1, float>, sofa::defaulttype::Vec<1, float>, float>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<1, float>, sofa::defaulttype::Vec<1, float>, float> > ExternalInterpolationMapping1f_to_1f;
typedef sofa::component::mapping::ExternalInterpolationMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<2, float>, sofa::defaulttype::Vec<2, float>, float>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<2, float>, sofa::defaulttype::Vec<2, float>, float> > ExternalInterpolationMapping2f_to_2f;
typedef sofa::component::mapping::ExternalInterpolationMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float>, sofa::defaulttype::ExtVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > ExternalInterpolationMapping3f_to_Ext3f;
typedef sofa::component::mapping::ExternalInterpolationMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > ExternalInterpolationMapping3f_to_3f;



//---------------------------------------------------------------------------------------------
//Typedef for HexahedronCompositeFEMMapping
typedef sofa::component::mapping::HexahedronCompositeFEMMapping<sofa::core::Mapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float>, sofa::defaulttype::ExtVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > > HexahedronCompositeFEMMapping3f_to_Ext3f;
typedef sofa::component::mapping::HexahedronCompositeFEMMapping<sofa::core::Mapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > > HexahedronCompositeFEMMapping3f_to_3f;



//---------------------------------------------------------------------------------------------
//Typedef for IdentityMapping
typedef sofa::component::mapping::IdentityMapping<sofa::defaulttype::StdRigidTypes<2, float>, sofa::defaulttype::StdRigidTypes<2, float> > IdentityMappingRigid2f_to_Rigid2f;
typedef sofa::component::mapping::IdentityMapping<sofa::defaulttype::StdRigidTypes<2, float>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<2, float>, sofa::defaulttype::Vec<2, float>, float> > IdentityMappingRigid2f_to_2f;
typedef sofa::component::mapping::IdentityMapping<sofa::defaulttype::StdRigidTypes<3, float>, sofa::defaulttype::ExtVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > IdentityMappingRigid3f_to_Ext3f;
typedef sofa::component::mapping::IdentityMapping<sofa::defaulttype::StdRigidTypes<3, float>, sofa::defaulttype::StdRigidTypes<3, float> > IdentityMappingRigid3f_to_Rigid3f;
typedef sofa::component::mapping::IdentityMapping<sofa::defaulttype::StdRigidTypes<3, float>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > IdentityMappingRigid3f_to_3f;
typedef sofa::component::mapping::IdentityMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<1, float>, sofa::defaulttype::Vec<1, float>, float>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<1, float>, sofa::defaulttype::Vec<1, float>, float> > IdentityMapping1f_to_1f;
typedef sofa::component::mapping::IdentityMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<2, float>, sofa::defaulttype::Vec<2, float>, float>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<2, float>, sofa::defaulttype::Vec<2, float>, float> > IdentityMapping2f_to_2f;
typedef sofa::component::mapping::IdentityMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float>, sofa::defaulttype::ExtVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > IdentityMapping3f_to_Ext3f;
typedef sofa::component::mapping::IdentityMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > IdentityMapping3f_to_3f;
typedef sofa::component::mapping::IdentityMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<6, float>, sofa::defaulttype::Vec<6, float>, float>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<6, float>, sofa::defaulttype::Vec<6, float>, float> > IdentityMapping6f_to_6f;



//---------------------------------------------------------------------------------------------
//Typedef for ImplicitSurfaceMapping
typedef sofa::component::mapping::ImplicitSurfaceMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float>, sofa::defaulttype::ExtVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > ImplicitSurfaceMapping3f_to_Ext3f;
typedef sofa::component::mapping::ImplicitSurfaceMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > ImplicitSurfaceMapping3f_to_3f;



//---------------------------------------------------------------------------------------------
//Typedef for LineSetSkinningMapping
typedef sofa::component::mapping::LineSetSkinningMapping<sofa::defaulttype::StdRigidTypes<3, float>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > LineSetSkinningMappingRigid3f_to_3f;



//---------------------------------------------------------------------------------------------
//Typedef for Mesh2PointMechanicalMapping
typedef sofa::component::mapping::Mesh2PointMechanicalMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float>, sofa::defaulttype::ExtVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > Mesh2PointMechanicalMapping3f_to_Ext3f;
typedef sofa::component::mapping::Mesh2PointMechanicalMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > Mesh2PointMechanicalMapping3f_to_3f;



//---------------------------------------------------------------------------------------------
//Typedef for RigidMapping
typedef sofa::component::mapping::RigidMapping<sofa::defaulttype::StdRigidTypes<2, float>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<2, float>, sofa::defaulttype::Vec<2, float>, float> > RigidMappingRigid2f_to_2f;
typedef sofa::component::mapping::RigidMapping<sofa::defaulttype::StdRigidTypes<3, float>, sofa::defaulttype::ExtVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > RigidMappingRigid3f_to_Ext3f;
typedef sofa::component::mapping::RigidMapping<sofa::defaulttype::StdRigidTypes<3, float>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > RigidMappingRigid3f_to_3f;



//---------------------------------------------------------------------------------------------
//Typedef for RigidRigidMapping
typedef sofa::component::mapping::RigidRigidMapping<sofa::defaulttype::StdRigidTypes<3, float>, sofa::defaulttype::StdRigidTypes<3, float> > RigidRigidMappingRigid3f_to_Rigid3f;



//---------------------------------------------------------------------------------------------
//Typedef for SPHFluidSurfaceMapping
typedef sofa::component::mapping::SPHFluidSurfaceMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float>, sofa::defaulttype::ExtVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > SPHFluidSurfaceMapping3f_to_Ext3f;
typedef sofa::component::mapping::SPHFluidSurfaceMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > SPHFluidSurfaceMapping3f_to_3f;



//---------------------------------------------------------------------------------------------
//Typedef for SimpleTesselatedTetraMechanicalMapping
typedef sofa::component::mapping::SimpleTesselatedTetraMechanicalMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float>, sofa::defaulttype::ExtVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > SimpleTesselatedTetraMechanicalMapping3f_to_Ext3f;
typedef sofa::component::mapping::SimpleTesselatedTetraMechanicalMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > SimpleTesselatedTetraMechanicalMapping3f_to_3f;



//---------------------------------------------------------------------------------------------
//Typedef for SkinningMapping
typedef sofa::component::mapping::SkinningMapping<sofa::defaulttype::StdRigidTypes<3, float>, sofa::defaulttype::ExtVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > SkinningMappingRigid3f_to_Ext3f;
typedef sofa::component::mapping::SkinningMapping<sofa::defaulttype::StdRigidTypes<3, float>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > SkinningMappingRigid3f_to_3f;



//---------------------------------------------------------------------------------------------
//Typedef for SubsetMapping
typedef sofa::component::mapping::SubsetMapping<sofa::defaulttype::StdRigidTypes<3, float>, sofa::defaulttype::StdRigidTypes<3, float> > SubsetMappingRigid3f_to_Rigid3f;
typedef sofa::component::mapping::SubsetMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<1, float>, sofa::defaulttype::Vec<1, float>, float>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<1, float>, sofa::defaulttype::Vec<1, float>, float> > SubsetMapping1f_to_1f;
typedef sofa::component::mapping::SubsetMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float>, sofa::defaulttype::ExtVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > SubsetMapping3f_to_Ext3f;
typedef sofa::component::mapping::SubsetMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > SubsetMapping3f_to_3f;



//---------------------------------------------------------------------------------------------
//Typedef for SubsetMultiMapping
typedef sofa::component::mapping::SubsetMultiMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > SubsetMultiMapping3f_to_3f;



//---------------------------------------------------------------------------------------------
//Typedef for TubularMapping
typedef sofa::component::mapping::TubularMapping<sofa::defaulttype::StdRigidTypes<3, float>, sofa::defaulttype::ExtVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > TubularMappingRigid3f_to_Ext3f;
typedef sofa::component::mapping::TubularMapping<sofa::defaulttype::StdRigidTypes<3, float>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > TubularMappingRigid3f_to_3f;





#ifdef SOFA_FLOAT
typedef ArticulatedSystemMapping1f_Rigid3f_to_Rigid3f ArticulatedSystemMapping1_Rigid3_to_Rigid3;
typedef BarycentricMapping3f_to_Ext3f BarycentricMapping3_to_Ext3;
typedef BarycentricMapping3f_to_3f BarycentricMapping3_to_3;
typedef BarycentricMapping3f_to_Rigid3f BarycentricMapping3_to_Rigid3;
typedef BeamLinearMappingRigid3f_to_Ext3f BeamLinearMappingRigid3_to_Ext3;
typedef BeamLinearMappingRigid3f_to_3f BeamLinearMappingRigid3_to_3;
typedef CenterOfMassMappingRigid2f_to_2f CenterOfMassMappingRigid2_to_2;
typedef CenterOfMassMappingRigid3f_to_Ext3f CenterOfMassMappingRigid3_to_Ext3;
typedef CenterOfMassMappingRigid3f_to_3f CenterOfMassMappingRigid3_to_3;
typedef CenterOfMassMulti2Mapping3f_Rigid3f_to_3f CenterOfMassMulti2Mapping3_Rigid3_to_3;
typedef CenterOfMassMultiMappingRigid3f_to_Rigid3f CenterOfMassMultiMappingRigid3_to_Rigid3;
typedef CenterOfMassMultiMapping3f_to_3f CenterOfMassMultiMapping3_to_3;
typedef CenterPointMechanicalMapping3f_to_Ext3f CenterPointMechanicalMapping3_to_Ext3;
typedef CenterPointMechanicalMapping3f_to_3f CenterPointMechanicalMapping3_to_3;
typedef CurveMapping3f_to_Rigid3f CurveMapping3_to_Rigid3;
typedef ExternalInterpolationMapping1f_to_1f ExternalInterpolationMapping1_to_1;
typedef ExternalInterpolationMapping2f_to_2f ExternalInterpolationMapping2_to_2;
typedef ExternalInterpolationMapping3f_to_Ext3f ExternalInterpolationMapping3_to_Ext3;
typedef ExternalInterpolationMapping3f_to_3f ExternalInterpolationMapping3_to_3;
typedef HexahedronCompositeFEMMapping3f_to_Ext3f HexahedronCompositeFEMMapping3_to_Ext3;
typedef HexahedronCompositeFEMMapping3f_to_3f HexahedronCompositeFEMMapping3_to_3;
typedef IdentityMappingRigid2f_to_Rigid2f IdentityMappingRigid2_to_Rigid2;
typedef IdentityMappingRigid2f_to_2f IdentityMappingRigid2_to_2;
typedef IdentityMappingRigid3f_to_Ext3f IdentityMappingRigid3_to_Ext3;
typedef IdentityMappingRigid3f_to_Rigid3f IdentityMappingRigid3_to_Rigid3;
typedef IdentityMappingRigid3f_to_3f IdentityMappingRigid3_to_3;
typedef IdentityMapping1f_to_1f IdentityMapping1_to_1;
typedef IdentityMapping2f_to_2f IdentityMapping2_to_2;
typedef IdentityMapping3f_to_Ext3f IdentityMapping3_to_Ext3;
typedef IdentityMapping3f_to_3f IdentityMapping3_to_3;
typedef IdentityMapping6f_to_6f IdentityMapping6_to_6;
typedef ImplicitSurfaceMapping3f_to_Ext3f ImplicitSurfaceMapping3_to_Ext3;
typedef ImplicitSurfaceMapping3f_to_3f ImplicitSurfaceMapping3_to_3;
typedef LineSetSkinningMappingRigid3f_to_3f LineSetSkinningMappingRigid3_to_3;
typedef Mesh2PointMechanicalMapping3f_to_Ext3f Mesh2PointMechanicalMapping3_to_Ext3;
typedef Mesh2PointMechanicalMapping3f_to_3f Mesh2PointMechanicalMapping3_to_3;
typedef RigidMappingRigid2f_to_2f RigidMappingRigid2_to_2;
typedef RigidMappingRigid3f_to_Ext3f RigidMappingRigid3_to_Ext3;
typedef RigidMappingRigid3f_to_3f RigidMappingRigid3_to_3;
typedef RigidRigidMappingRigid3f_to_Rigid3f RigidRigidMappingRigid3_to_Rigid3;
typedef SPHFluidSurfaceMapping3f_to_Ext3f SPHFluidSurfaceMapping3_to_Ext3;
typedef SPHFluidSurfaceMapping3f_to_3f SPHFluidSurfaceMapping3_to_3;
typedef SimpleTesselatedTetraMechanicalMapping3f_to_Ext3f SimpleTesselatedTetraMechanicalMapping3_to_Ext3;
typedef SimpleTesselatedTetraMechanicalMapping3f_to_3f SimpleTesselatedTetraMechanicalMapping3_to_3;
typedef SkinningMappingRigid3f_to_Ext3f SkinningMappingRigid3_to_Ext3;
typedef SkinningMappingRigid3f_to_3f SkinningMappingRigid3_to_3;
typedef SubsetMappingRigid3f_to_Rigid3f SubsetMappingRigid3_to_Rigid3;
typedef SubsetMapping1f_to_1f SubsetMapping1_to_1;
typedef SubsetMapping3f_to_Ext3f SubsetMapping3_to_Ext3;
typedef SubsetMapping3f_to_3f SubsetMapping3_to_3;
typedef SubsetMultiMapping3f_to_3f SubsetMultiMapping3_to_3;
typedef TubularMappingRigid3f_to_Ext3f TubularMappingRigid3_to_Ext3;
typedef TubularMappingRigid3f_to_3f TubularMappingRigid3_to_3;
#endif

#endif

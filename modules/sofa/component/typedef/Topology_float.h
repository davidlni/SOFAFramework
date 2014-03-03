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


#ifndef SOFA_TYPEDEF_Topology_float_H
#define SOFA_TYPEDEF_Topology_float_H

//Default files containing the declaration of the vector type
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/defaulttype/Mat.h>




#include <VolumetricData/mapping/ImplicitSurfaceMapping.h>
#include <SphFluid/SPHFluidSurfaceMapping.h>



//---------------------------------------------------------------------------------------------
//Typedef for ImplicitSurfaceMapping
typedef sofa::component::mapping::ImplicitSurfaceMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float>, sofa::defaulttype::ExtVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > ImplicitSurfaceMapping3f_Ext3f;
typedef sofa::component::mapping::ImplicitSurfaceMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > ImplicitSurfaceMapping3f_3f;



//---------------------------------------------------------------------------------------------
//Typedef for SPHFluidSurfaceMapping
typedef sofa::component::mapping::SPHFluidSurfaceMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float>, sofa::defaulttype::ExtVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > SPHFluidSurfaceMapping3f_Ext3f;
typedef sofa::component::mapping::SPHFluidSurfaceMapping<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float>, sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, float>, sofa::defaulttype::Vec<3, float>, float> > SPHFluidSurfaceMapping3f_3f;





#ifdef SOFA_FLOAT
typedef ImplicitSurfaceMapping3f_Ext3f ImplicitSurfaceMapping3_Ext3;
typedef ImplicitSurfaceMapping3f_3f ImplicitSurfaceMapping3_3;
typedef SPHFluidSurfaceMapping3f_Ext3f SPHFluidSurfaceMapping3_Ext3;
typedef SPHFluidSurfaceMapping3f_3f SPHFluidSurfaceMapping3_3;
#endif

#endif

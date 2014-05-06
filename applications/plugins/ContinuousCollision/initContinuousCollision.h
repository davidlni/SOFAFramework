#ifndef INITContinuousCollision_H
#define INITContinuousCollision_H
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

#include <sofa/helper/system/config.h>

#ifdef SOFA_BUILD_PLUGIN_CONTINUOUS_COLLISION
#define SOFA_CONTINUOUS_COLLISION_API SOFA_EXPORT_DYNAMIC_LIBRARY
#else
#define SOFA_CONTINUOUS_COLLISION_API  SOFA_IMPORT_DYNAMIC_LIBRARY
#endif

/** \mainpage
  This is a the starting page of the plugin documentation, defined in file initContinuousCollision.h
  */

#endif // INITContinuousCollision_H

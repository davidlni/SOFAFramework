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
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#ifndef FLEXIBLE_AffineTYPES_H
#define FLEXIBLE_AffineTYPES_H

#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/defaulttype/Mat.h>
#include <sofa/helper/vector.h>
#include <sofa/helper/rmath.h>
#include <sofa/helper/decompose.h>
#ifdef SOFA_SMP
#include <sofa/defaulttype/SharedTypes.h>
#endif /* SOFA_SMP */

#include "DeformableFrameMass.h"

namespace sofa
{

namespace defaulttype
{

using std::endl;
using helper::vector;

/** DOF types associated with deformable frames. Each deformable frame generates an affine displacement field, with 12 independent degrees of freedom.
 */
template<int _spatial_dimensions, typename _Real>
class StdAffineTypes
{
public:
    static const unsigned int spatial_dimensions = _spatial_dimensions;  ///< Number of dimensions the frame is moving in, typically 3
    static const unsigned int VSize = spatial_dimensions +  spatial_dimensions * spatial_dimensions;  // number of entries
    enum { coord_total_size = VSize };
    enum { deriv_total_size = VSize };
    typedef _Real Real;
    typedef vector<Real> VecReal;

    // ------------    Types and methods defined for easier data access
    typedef Vec<spatial_dimensions, Real> SpatialCoord;                   ///< Position or velocity of a point
    typedef Mat<spatial_dimensions,spatial_dimensions, Real> Frame;       ///< Matrix representing a frame

    typedef SpatialCoord CPos;
    typedef SpatialCoord DPos;

    class Coord : public Vec<VSize,Real>
    {
        typedef Vec<VSize,Real> MyVec;

    public:

        enum { spatial_dimensions = _spatial_dimensions }; // different from Vec::spatial_dimensions == 12

        Coord() { clear(); }
        Coord( const Vec<VSize,Real>& d):MyVec(d) {}
        Coord( const SpatialCoord& c, const Frame& a) { getCenter()=c; getAffine()=a;}
        void clear()  { MyVec::clear(); for(unsigned int i=0; i<_spatial_dimensions; ++i) getAffine()[i][i]=(Real)1.0; } // init affine part to identity

        typedef Real value_type;

        /// point
        SpatialCoord& getCenter() { return *reinterpret_cast<SpatialCoord*>(&this->elems[0]); }
        const SpatialCoord& getCenter() const { return *reinterpret_cast<const SpatialCoord*>(&this->elems[0]); }

        /// local frame
        Frame& getAffine() { return *reinterpret_cast<Frame*>(&this->elems[_spatial_dimensions]); }
        const Frame& getAffine() const { return *reinterpret_cast<const Frame*>(&this->elems[_spatial_dimensions]); }


        /// write to an output stream
        inline friend std::ostream& operator << ( std::ostream& out, const Coord& c )
        {
            out<<c.getCenter()<<" "<<c.getAffine();
            return out;
        }
        /// read from an input stream
        inline friend std::istream& operator >> ( std::istream& in, Coord& c )
        {
            in>>c.getCenter()>>c.getAffine();
            return in;
        }

        /// Write the OpenGL transformation matrix
        void writeOpenGlMatrix ( float m[16] ) const
        {
            BOOST_STATIC_ASSERT(_spatial_dimensions == 3);
            m[0] = (float)getAffine()(0,0);
            m[4] = (float)getAffine()(0,1);
            m[8] = (float)getAffine()(0,2);
            m[1] = (float)getAffine()(1,0);
            m[5] = (float)getAffine()(1,1);
            m[9] = (float)getAffine()(1,2);
            m[2] = (float)getAffine()(2,0);
            m[6] = (float)getAffine()(2,1);
            m[10] = (float)getAffine()(2,2);
            m[3] = 0;
            m[7] = 0;
            m[11] = 0;
            m[12] = ( float ) getCenter()[0];
            m[13] = ( float ) getCenter()[1];
            m[14] = ( float ) getCenter()[2];
            m[15] = 1;
        }

        /// Project a point from the child frame to the parent frame
        SpatialCoord pointToParent( const SpatialCoord& v ) const
        {
            return getAffine()*v + getCenter();
        }

        /// Project a point from the parent frame to the child frame
        SpatialCoord pointToChild( const SpatialCoord& v ) const
        {
            Frame affineInv;
            affineInv.invert( getAffine() );
            return affineInv * (v-getCenter());
        }

        /// Project a point from the child frame to the parent frame
        Coord pointToParent( const Coord& v ) const
        {
            return Coord( getAffine()*v.getCenter() + getCenter(), getAffine()*v.getAffine() );
        }

        /// Project a point from the parent frame to the child frame
        Coord pointToChild( const Coord& v ) const
        {
            Frame affineInv;
            affineInv.invert( getAffine() );
            return Coord( affineInv * (v.getCenter()-getCenter()), affineInv*v.getAffine() );
        }

        /// project to a rigid displacement
        void setRigid()
        {
            Frame rotation;
            helper::Decompose<Real>::polarDecomposition( getAffine(), rotation );
            getAffine() = rotation;
        }
    };

    typedef vector<Coord> VecCoord;

    static const char* Name();


    static Coord interpolate ( const helper::vector< Coord > & ancestors, const helper::vector< Real > & coefs )
    {
        assert ( ancestors.size() == coefs.size() );
        Coord c;
        for ( unsigned int i = 0; i < ancestors.size(); i++ ) c += ancestors[i] * coefs[i];  // Position and deformation gradient linear interpolation.
        return c;
    }

    static Coord inverse( const Coord& c )
    {
        Frame m;
#ifdef DEBUG
        bool invertible = invertMatrix(m,c.getAffine());
        assert(invertible);
#else
        invertMatrix(m,c.getAffine());
#endif
        return Coord( -(m*c.getCenter()),m );
    }


    class Deriv : public Vec<VSize,Real>
    {
        typedef Vec<VSize,Real> MyVec;
    public:

        enum { spatial_dimensions = _spatial_dimensions }; // different from Vec::spatial_dimensions == 12

        Deriv() { MyVec::clear(); }
        Deriv( const Vec<VSize,Real>& d):MyVec(d) {}
        Deriv( const SpatialCoord& c, const Frame& a) { getVCenter()=c; getVAffine()=a;}

        //static const unsigned int total_size = VSize;
        typedef Real value_type;

        /// point
        SpatialCoord& getVCenter() { return *reinterpret_cast<SpatialCoord*>(&this->elems[0]); }
        const SpatialCoord& getVCenter() const { return *reinterpret_cast<const SpatialCoord*>(&this->elems[0]); }

        /// local frame
        Frame& getVAffine() { return *reinterpret_cast<Frame*>(&this->elems[spatial_dimensions]); }
        const Frame& getVAffine() const { return *reinterpret_cast<const Frame*>(&this->elems[spatial_dimensions]); }

        /// project to a rigid motion
        void setRigid(const Coord& c)
        {
            // Compute velocity tensor W = Adot.Ainv
            Frame Ainv;  invertMatrix(Ainv,c.getAffine());
            Frame W = getVAffine() * Ainv;

            // make it skew-symmetric
            for(unsigned i=0; i<spatial_dimensions; i++) W[i][i] = 0.0;
            for(unsigned i=0; i<spatial_dimensions; i++)
            {
                for(unsigned j=i+1; j<spatial_dimensions; j++)
                {
                    W[i][j] = (W[i][j] - W[j][i]) *(Real)0.5;
                    W[j][i] = - W[i][j];
                }
            }

            // retrieve global velocity : Rdot = W.R
            Frame R;
            helper::Decompose<Real>::polarDecomposition( c.getAffine() , R );
            getVAffine() = W*R;
        }

        template< int N, class Real2 > // N <= VSize
        void operator+=( const Vec<N,Real2>& p ) { for(int i=0;i<N;++i) this->elems[i] += (Real)p[i]; }
        template< int N, class Real2 > // N <= VSize
        void operator=( const Vec<N,Real2>& p ) { for(int i=0;i<N;++i) this->elems[i] = (Real)p[i]; }

    };

    typedef vector<Deriv> VecDeriv;
    typedef MapMapSparseMatrix<Deriv> MatrixDeriv;

    static Deriv interpolate ( const helper::vector< Deriv > & ancestors, const helper::vector< Real > & coefs )
    {
        assert ( ancestors.size() == coefs.size() );
        Deriv c;
        for ( unsigned int i = 0; i < ancestors.size(); i++ )     c += ancestors[i] * coefs[i];
        return c;
    }



    /** @name Conversions
              * Convert to/from points in space
             */
    //@{
    static const SpatialCoord& getCPos(const Coord& c) { return c.getCenter(); }
    static void setCPos(Coord& c, const SpatialCoord& v) { c.getCenter() = v; }

    static const SpatialCoord& getDPos(const Deriv& d) { return d.getVCenter(); }
    static void setDPos(Deriv& d, const SpatialCoord& v) { d.getVCenter() = v; }

    template<typename T>
    static void set ( Deriv& c, T x, T y, T z )
    {
        c.clear();
        c.getVCenter()[0] = ( Real ) x;
        c.getVCenter()[1] = ( Real ) y;
        c.getVCenter()[2] = ( Real ) z;
    }

    template<typename T>
    static void get ( T& x, T& y, T& z, const Deriv& c )
    {
        x = ( T ) c.getVCenter() [0];
        y = ( T ) c.getVCenter() [1];
        z = ( T ) c.getVCenter() [2];
    }

    template<typename T>
    static void add ( Deriv& c, T x, T y, T z )
    {
        c.getVCenter() [0] += ( Real ) x;
        c.getVCenter() [1] += ( Real ) y;
        c.getVCenter() [2] += ( Real ) z;
    }

    template<typename T>
    static void set ( Coord& c, T x, T y, T z )
    {
        c.clear();
        c.getCenter()[0] = ( Real ) x;
        c.getCenter() [1] = ( Real ) y;
        c.getCenter() [2] = ( Real ) z;
    }

    template<typename T>
    static void get ( T& x, T& y, T& z, const Coord& c )
    {
        x = ( T ) c.getCenter() [0];
        y = ( T ) c.getCenter() [1];
        z = ( T ) c.getCenter() [2];
    }

    template<typename T>
    static void add ( Coord& c, T x, T y, T z )
    {
        c.getCenter() [0] += ( Real ) x;
        c.getCenter() [1] += ( Real ) y;
        c.getCenter() [2] += ( Real ) z;
    }
    //@}

};



#ifndef SOFA_FLOAT
typedef StdAffineTypes<3, double> Affine3dTypes;

// Specialization of the defaulttype::DataTypeInfo type traits template
template<> struct DataTypeInfo< sofa::defaulttype::Affine3dTypes::Coord > : public FixedArrayTypeInfo< sofa::defaulttype::Affine3dTypes::Coord, sofa::defaulttype::Affine3dTypes::Coord::total_size >
{
    static std::string name() { std::ostringstream o; o << "AffineCoord<" << sofa::defaulttype::Affine3dTypes::Coord::total_size << "," << DataTypeName<sofa::defaulttype::Affine3dTypes::Real>::name() << ">"; return o.str(); }
};
template<> struct DataTypeInfo< sofa::defaulttype::Affine3dTypes::Deriv > : public FixedArrayTypeInfo< sofa::defaulttype::Affine3dTypes::Deriv, sofa::defaulttype::Affine3dTypes::Deriv::total_size >
{
    static std::string name() { std::ostringstream o; o << "AffineDeriv<" << sofa::defaulttype::Affine3dTypes::Deriv::total_size << "," << DataTypeName<sofa::defaulttype::Affine3dTypes::Real>::name() << ">"; return o.str(); }
};
#endif
#ifndef SOFA_DOUBLE
typedef StdAffineTypes<3, float> Affine3fTypes;

// Specialization of the defaulttype::DataTypeInfo type traits template
template<> struct DataTypeInfo< sofa::defaulttype::Affine3fTypes::Coord > : public FixedArrayTypeInfo< sofa::defaulttype::Affine3fTypes::Coord, sofa::defaulttype::Affine3fTypes::Coord::total_size >
{
    static std::string name() { std::ostringstream o; o << "AffineCoord<" << sofa::defaulttype::Affine3fTypes::Coord::total_size << "," << DataTypeName<sofa::defaulttype::Affine3fTypes::Real>::name() << ">"; return o.str(); }
};
template<> struct DataTypeInfo< sofa::defaulttype::Affine3fTypes::Deriv > : public FixedArrayTypeInfo< sofa::defaulttype::Affine3fTypes::Deriv, sofa::defaulttype::Affine3fTypes::Deriv::total_size >
{
    static std::string name() { std::ostringstream o; o << "AffineDeriv<" << sofa::defaulttype::Affine3fTypes::Deriv::total_size << "," << DataTypeName<sofa::defaulttype::Affine3fTypes::Real>::name() << ">"; return o.str(); }
};
#endif

/// Note: Many scenes use Affine as template for 3D double-precision rigid type. Changing it to Affine3d would break backward compatibility.
#ifdef SOFA_FLOAT
template<> inline const char* Affine3fTypes::Name() { return "Affine"; }
#else
template<> inline const char* Affine3dTypes::Name() { return "Affine"; }
#ifndef SOFA_DOUBLE
template<> inline const char* Affine3fTypes::Name() { return "Affine3f"; }
#endif
#endif

#ifdef SOFA_FLOAT
typedef Affine3fTypes Affine3Types;
#else
typedef Affine3dTypes Affine3Types;
#endif
//typedef Affine3Types AffineTypes;







// The next line hides all those methods from the doxygen documentation
/// \cond TEMPLATE_OVERRIDES


#ifndef SOFA_FLOAT
template<> struct DataTypeName< defaulttype::Affine3dTypes::Coord > { static const char* name() { return "Affine3dTypes::Coord"; } };
#endif
#ifndef SOFA_DOUBLE
template<> struct DataTypeName< defaulttype::Affine3fTypes::Coord > { static const char* name() { return "Affine3fTypes::Coord"; } };
#endif


/// \endcond




// ====================================================================
// AffineMass


#ifndef SOFA_FLOAT
typedef DeformableFrameMass<3, StdAffineTypes<3,double>::deriv_total_size, double> Affine3dMass;
#endif
#ifndef SOFA_DOUBLE
typedef DeformableFrameMass<3, StdAffineTypes<3,float>::deriv_total_size, float> Affine3fMass;
#endif

#ifdef SOFA_FLOAT
typedef Affine3fMass Affine3Mass;
#else
typedef Affine3dMass Affine3Mass;
#endif



// The next line hides all those methods from the doxygen documentation
/// \cond TEMPLATE_OVERRIDES

#ifndef SOFA_FLOAT
template<> struct DataTypeName< defaulttype::Affine3dMass > { static const char* name() { return "Affine3dMass"; } };
#endif
#ifndef SOFA_DOUBLE
template<> struct DataTypeName< defaulttype::Affine3fMass > { static const char* name() { return "Affine3fMass"; } };
#endif

/// \endcond



} // namespace defaulttype





} // namespace sofa



#endif

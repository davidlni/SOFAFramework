
#ifndef IMAGE_BranchingImage_H
#define IMAGE_BranchingImage_H



#include <sofa/defaulttype/Vec.h>
#include <sofa/helper/vector.h>
#include <map>

#include "ImageTypes.h"
#include "Containers.h"

namespace sofa
{

namespace defaulttype
{

using helper::vector;
using helper::NoPreallocationVector;

/// type identifier, must be unique
static const int IMAGELABEL_BRANCHINGIMAGE = 1;







/// index indentifying a single connectionVoxel = index1D + SuperimposedIndex
struct BranchingImageVoxelIndex
{
    BranchingImageVoxelIndex() : index1d(-1), offset(-1) {}

    BranchingImageVoxelIndex( unsigned index1d, unsigned offset ) : index1d(index1d), offset(offset) {}

    BranchingImageVoxelIndex( const BranchingImageVoxelIndex& other )
    {
        index1d = other.index1d;
        offset = other.offset;
    }

    unsigned index1d; // in BranchingImage3D
    unsigned offset; // in SuperimposedVoxels vector

    void operator=( const BranchingImageVoxelIndex& other )
    {
        index1d = other.index1d;
        offset = other.offset;
    }

    bool operator==( const BranchingImageVoxelIndex& other ) const
    {
        return index1d==other.index1d && offset==other.offset;
    }

    bool operator<( const BranchingImageVoxelIndex& other ) const
    {
        return index1d<other.index1d || offset<other.offset;
    }

}; // struct BranchingImageVoxelIndex


/// offsets in each direction x/y/z around a voxel
struct BranchingImageNeighbourOffset
{
    typedef enum { FACE=1, EDGE=2, CORNER=3, ONPLACE=0, NOTCLOSE=4 } ConnectionType; // warning: do not change enum values they are significant

    /// default constructor
    /// @warning no initialization
    BranchingImageNeighbourOffset() {}

    BranchingImageNeighbourOffset( const BranchingImageNeighbourOffset& other )
    {
        *this = other;
    }

    BranchingImageNeighbourOffset( int xOffset, int yOffset, int zOffset )
    {
        set( xOffset, yOffset, zOffset );
    }

    inline void operator=( const BranchingImageNeighbourOffset& other )
    {
        memcpy( offset, other.offset, 3*sizeof(int) );
    }

    inline void set( int xOffset, int yOffset, int zOffset )
    {
        offset[0] = xOffset;
        offset[1] = yOffset;
        offset[2] = zOffset;
    }

    /// \returns the opposite direction of a given direction  left->right,  right->left
    inline BranchingImageNeighbourOffset opposite() const
    {
        BranchingImageNeighbourOffset op( *this );
        for( unsigned int i=0 ; i<3 ; ++i )
            if( op[i] ) op[i] = -op[i];
        return op;
    }

    /// dir represent the direction  0->x, 1->y, 2->z
    inline int& operator[]( unsigned dir ) { assert( dir<=2 ); return offset[dir]; }
    inline const int& operator[]( unsigned dir ) const { assert( dir<=2 ); return offset[dir]; }

    bool operator==( const BranchingImageNeighbourOffset& other ) const
    {
        return offset[0]==other.offset[0] && offset[1]==other.offset[1] && offset[2]==other.offset[2];
    }

    bool operator!=( const BranchingImageNeighbourOffset& other ) const
    {
        return offset[0]!=other.offset[0] || offset[1]!=other.offset[1] || offset[2]!=other.offset[2];
    }

    ConnectionType connectionType() const
    {
        for( unsigned int i=0 ; i<3 ; ++i )
            if( offset[i] < -1 || offset[i] > 1 ) return NOTCLOSE;
        return ConnectionType( abs(offset[0])+abs(offset[1])+abs(offset[2]) );
    }

protected:

    int offset[3]; ///< in each direction x,y,z, the relative offset of the neighour

}; // struct BranchingImageNeighbourOffset



typedef enum { CONNECTIVITY_6=6, CONNECTIVITY_26=26 } BranchingImageConnectivity;



/// A BranchingImage is an array (size of t) of vectors (one vector per pixel (x,y,z)).
/// Each pixel corresponds to a SuperimposedVoxels, alias a vector of ConnectionVoxel.
/// a ConnectionVoxel stores a value for each channels + its neighbours indices
/// Nesme, Kry, Jeřábková, Faure, "Preserving Topology and Elasticity for Embedded Deformable Models", Siggraph09
template<typename _T>
class BranchingImage
{

public:

    /// stored type
    typedef _T T;

    /// type identifier, must be unique
    static const int label = IMAGELABEL_BRANCHINGIMAGE;


    typedef BranchingImageVoxelIndex VoxelIndex;
    typedef BranchingImageNeighbourOffset NeighbourOffset;


    /// a list of neighbours is stored as a vector of VoxelIndex
    typedef NoPreallocationVector<VoxelIndex> Neighbours;


    /// a ConnectionVoxel stores a value for each channels + its neighbour indices /*+ its 1D index in the image*/
    /// NB: a ConnectionVoxel does not know its spectrum size (nb of channels) to save memory
    class ConnectionVoxel
    {

    public:

        /// default constructor = no allocation
        ConnectionVoxel() : value(0), neighbours() {}
        /// with allocation constructor
        ConnectionVoxel( size_t spectrum ) : neighbours() { if( !spectrum ) { value = 0; return; } value = new T[spectrum]; }
        ~ConnectionVoxel() { if( value ) delete [] value; }

        /// copy
        /// @warning neighbourhood is copied but impossible to check its validity
        void clone( const ConnectionVoxel& cv, unsigned spectrum )
        {
            equal( cv, spectrum ); // copy value

            neighbours = cv.neighbours;

            //index = cv.index;
        }

        /// copy only the topology and initialize value size acoording to spectrum
        template<typename T2>
        void cloneTopology( const typename BranchingImage<T2>::ConnectionVoxel& cv, unsigned spectrum, const T defaultValue=(T)0 )
        {
            if( value ) delete [] value;

            if( !spectrum || !cv.value ) { value=0; return; }

            value = new T[spectrum];
            for( unsigned i=0 ; i<spectrum ; ++i ) value[i]=defaultValue;

            neighbours = cv.neighbours;

            //index = cv.index;
        }


        /// alloc or realloc without keeping existing data and without initialization
        void resize( size_t newSize )
        {
            if( value ) delete [] value;
            if( !newSize ) { value=0; return; }
            value = new T[newSize];
            // could do a realloc keeping existing data
            // should it clear neighbourhood?
        }

        /// computes a norm over all channels
        double magnitude( unsigned spectrum, const int magnitude_type=2 ) const
        {
            double res = 0;
            switch (magnitude_type) {
            case -1 : {
                for( unsigned i=0 ; i<spectrum ; ++i ) { const double val = (double)abs(value[i]); if (val>res) res = val; }
            } break;
            case 1 : {
                for( unsigned i=0 ; i<spectrum ; ++i ) res += (double)abs(value[i]);
            } break;
            default : {
                for( unsigned i=0 ; i<spectrum ; ++i ) res += (double)(value[i]*value[i]);
                res = (double)sqrt(res);
            }
            }
            return res;
        }

        /// @returns the min channel value
        T min( unsigned spectrum ) const
        {
            if( !value ) return 0;
            T m = value[0];
            for( unsigned i=1 ; i<spectrum ; ++i ) if( value[i]<m ) m=value[i];
            return m;
        }

        /// @returns the max channel value
        T max( unsigned spectrum ) const
        {
            if( !value ) return 0;
            T m = value[0];
            for( unsigned i=1 ; i<spectrum ; ++i ) if( value[i]>m ) m=value[i];
            return m;
        }

        /// @returns true iff all channels are 0
        bool empty( unsigned spectrum ) const
        {
            if( !value ) return true;
            for( unsigned i=0 ; i<spectrum ; ++i ) if( value[i] ) return false;
            return true;
        }

        //unsigned index; ///< the 1D position in the BranchingImage3D // @TODO is it really necessary

        T* value; ///< value of the voxel for each channel (value is the size of the C dimension of the ConnectionImage)

        Neighbours neighbours;


        /// accessor
        /// @warning index must be less than the spectrum
        const T& operator [] (size_t index) const
        {
            return value[ index ];
        }
        T& operator [] (size_t index)
        {
            return value[ index ];
        }

        /// ==
        bool isEqual( const ConnectionVoxel& other, unsigned spectrum ) const
        {
            for( unsigned i=0 ; i<spectrum ; ++i )
                if( value[i] != other.value[i] ) return false;
            return true;
        }

        /// =
        void equal( const ConnectionVoxel& other, unsigned spectrum )
        {
            if( value ) delete [] value;

            if( !spectrum || !other.value ) { value=0; return; }

            value = new T[spectrum];
            memcpy( value, other.value, spectrum*sizeof(T) );
        }

        /// +=
        void addEqual( const ConnectionVoxel& other, unsigned spectrum )
        {
            for( unsigned i=0 ; i<spectrum ; ++i )
                value[i] += other.value[i];
        }

        /// /=
        void divideEqual( SReal d, unsigned spectrum )
        {
            for( unsigned i=0 ; i<spectrum ; ++i )
                value[i] /= d;
        }


        /// add the given voxel as a neigbour
        /// @warning it is doing only one way (this has to be added as a neighbour of n)
        /// if testUnicity==true, the neighbour is added only if it is not already there
        void addNeighbour( const VoxelIndex& neighbour, bool testUnicity = true )
        {
            if( !testUnicity || !isNeighbour(neighbour) ) neighbours.push_back( neighbour );
        }

        /// is the given voxel a neighbour of this voxel?
        bool isNeighbour( const VoxelIndex& neighbour ) const
        {
            return neighbours.isPresent( neighbour );
        }


    private:

        // no pure copy constructor (spectrum is needed to be able to copy the ConnectionVoxel)
        ConnectionVoxel( const ConnectionVoxel& ) { assert(false); }
        void operator=( const ConnectionVoxel& ) { assert(false); }
        bool operator==( const ConnectionVoxel& ) const { assert(false); return false; }

    }; // class ConnectionVoxel




    /// An array of ConnectionVoxel
    class SuperimposedVoxels : public NoPreallocationVector< ConnectionVoxel >
    {

    public:

        typedef NoPreallocationVector< ConnectionVoxel > Inherited;

        SuperimposedVoxels() : Inherited() {}

        /// add a superimposed voxel
        void push_back( const ConnectionVoxel& v, unsigned spectrum )
        {
            size_t newSize = this->_size+1;

            if( !this->_array )
            {
                this->_size = newSize;
                this->_array = new ConnectionVoxel[this->_size];
            }
            else
            {
                ConnectionVoxel* oldArray = this->_array;
                this->_array = new ConnectionVoxel[newSize];
                for( size_t i=0 ; i<this->_size ; ++i )
                    this->_array[i].clone( oldArray[i], spectrum );
                this->_size = newSize;
                delete [] oldArray;
            }

            this->last().clone( v, spectrum );
        }

        /// copy superimposed voxels
        /// @warning about voxel connectivity
        void clone( const SuperimposedVoxels& other, unsigned spectrum )
        {
            if( other.empty() ) Inherited::clear();
            else
            {
                Inherited::resize( other._size );
                for( unsigned i=0 ; i<this->_size ; ++i )
                {
                    this->_array[i].clone( other._array[i], spectrum );
                }
            }
        }

        /// copy only the topology and initialize value size acoording to spectrum
        template<typename T2>
        void cloneTopology( const typename BranchingImage<T2>::SuperimposedVoxels& other, unsigned spectrum , const T defaultValue=(T)0 )
        {
            if( other.empty() ) Inherited::clear();
            else
            {
                Inherited::resize( other.size());
                for( unsigned i=0 ; i<this->_size ; ++i )
                {
                    this->_array[i].template cloneTopology<T2>( other[i], spectrum, defaultValue );
                }
            }
        }

        /// equivalent to ==
        bool isEqual( const SuperimposedVoxels& other, unsigned spectrum ) const
        {
            if( this->_size != other._size ) return false;
            for( unsigned i=0 ; i<this->_size ; ++i )
                if( !this->_array[i].isEqual( other._array[i], spectrum ) ) return false;
            return true;
        }

        /// convert to a unique voxel
        /// conversionType : 0->first voxel, 1->average, 2->nb superimposed voxels, 3->sum
        template<class T2>
        void toFlatVoxel( T2& v, unsigned conversionType, unsigned channel ) const
        {
            if( this->empty() ) return;

            switch( conversionType )
            {
            case 1: // average
                assert( this->_array[0][channel] <= std::numeric_limits<T2>::max() );
                v = (T2)this->_array[0][channel];
                for( unsigned i=1 ; i<this->_size ; ++i )
                {
                    assert( v <= std::numeric_limits<T2>::max()-this->_array[i][channel] );
                    v += (T2)this->_array[i][channel];
                }
                v = (T2)( v / (float)this->_size );
                break;
            case 0: // first
                assert( this->_array[0][channel] <= std::numeric_limits<T2>::max() );
                v = (T2)this->_array[0][channel];
                break;
            case 3: // sum
                assert( typeid(T)==typeid(bool) || this->_array[0][channel] <= std::numeric_limits<T2>::max() );
                v = (T2)this->_array[0][channel];
                for( unsigned i=1 ; i<this->_size ; ++i )
                {
                    assert( typeid(T)==typeid(bool) || v <= std::numeric_limits<T2>::max()-this->_array[i][channel] );
                    v += (T2)this->_array[i][channel];
                }
                break;
            case 2: // count
            default:
                assert( typeid(T)==typeid(bool) || this->_size <= (unsigned)std::numeric_limits<T2>::max() );
                v = (T2)this->_size;
                break;
            }
        }

        void resize( size_t newSize, unsigned spectrum )
        {
            Inherited::resize( newSize );
            for( unsigned i=0 ; i<this->_size ; ++i )
                this->_array[i].resize( spectrum );
        }

        void fill( const ConnectionVoxel& v, unsigned spectrum )
        {
            for( unsigned i=0 ; i<this->_size ; ++i )
                this->_array[i].clone( v, spectrum );
        }


        /// all needed operators +, +=, etc. can be overloaded here


    private :

        /// impossible to copy a ConnectedVoxel without the spectrum size
        void push_back( const ConnectionVoxel& ) { assert(false); }
        SuperimposedVoxels( const SuperimposedVoxels& ) : Inherited() { assert(false); }
        void operator=( const SuperimposedVoxels& ) { assert(false); }
        bool operator==( const SuperimposedVoxels& ) const { assert(false); return false; }
        void resize( size_t ) { assert(false); }
        void resizeAndKeep( size_t ) { assert(false); }
        void fill( const ConnectionVoxel& ) { assert(false); }

    }; // class SuperimposedVoxels



    /// a BranchingImage is a dense image with a vector of SuperimposedVoxels at each pixel
    class BranchingImage3D : public NoPreallocationVector<SuperimposedVoxels>
    {
    public:

        typedef NoPreallocationVector<SuperimposedVoxels> Inherited;

        BranchingImage3D() : Inherited() {}

        /// copy
        void clone( const BranchingImage3D& other, unsigned spectrum )
        {
            this->resize( other._size );
            for( unsigned i=0 ; i<this->_size ; ++i )
            {
                this->_array[i].clone( other._array[i], spectrum );
            }
        }

        /// copy only the topology and initialize value size acoording to spectrum
        template<typename T2>
        void cloneTopology( const typename BranchingImage<T2>::BranchingImage3D& other, unsigned spectrum, const T defaultValue=(T)0)
        {
            this->resize( other.size() );
            for( unsigned i=0 ; i<this->_size ; ++i )
            {
                this->_array[i].template cloneTopology<T2>( other[i], spectrum, defaultValue );
            }
        }

        /// equivalent to ==
        bool isEqual( const BranchingImage3D& other, unsigned spectrum ) const
        {
            if( this->_size != other._size ) return false;
            for( unsigned i=0 ; i<this->_size ; ++i )
                if( !this->_array[i].isEqual( other._array[i], spectrum ) ) return false;
            return true;
        }

        /// \returns the offset of the given ConnectionVoxel in its SuperImposedVoxels vector
        int getOffset( unsigned index1d, const ConnectionVoxel& v ) const
        {
            return this->_array[index1d].getOffset( &v );
        }

        const ConnectionVoxel& operator()( const VoxelIndex& index ) const { return (*this)[index.index1d][index.offset]; }
              ConnectionVoxel& operator()( const VoxelIndex& index )       { return (*this)[index.index1d][index.offset]; }

    private:

        /// impossible to copy a ConnectedVoxel without the spectrum size
        BranchingImage3D( const BranchingImage3D& ) : Inherited() { assert(false); }
        void operator=( const BranchingImage3D& ) { assert(false); }
        bool operator==( const BranchingImage3D& ) const { assert(false); return false; }
        void push_back( const SuperimposedVoxels& ) { assert(false); }
        void resizeAndKeep( size_t ) { assert(false); }
        void fill( const SuperimposedVoxels& ) { assert(false); }

    }; // class BranchingImage3D




    /////////////////////////////////////////
    ////// BRANCHING IMAGE
    /////////////////////////////////////////


    /// the 5 dimension labels of an image ( x, y, z, spectrum=nb channels , time )
    typedef enum{ DIMENSION_X=0, DIMENSION_Y, DIMENSION_Z, DIMENSION_S /* spectrum = nb channels*/, DIMENSION_T /*4th dimension = time*/, NB_DimensionLabel } DimensionLabel;
    /// the 5 dimensions of an image ( x, y, z, spectrum=nb channels , time )
    typedef Vec<NB_DimensionLabel,unsigned int> Dimension; // [x,y,z,s,t]
    typedef Dimension imCoord; // to have a common api with Image

protected:

    /// @todo maybe this could be implicit?
    struct BranchingImageDimension
    {
        BranchingImageDimension() : dimension(), sliceSize(0), imageSize(0) {}
        void set( const Dimension& d )
        {
            dimension = d;
            sliceSize = dimension[DIMENSION_X] * dimension[DIMENSION_Y];
            imageSize = sliceSize * dimension[DIMENSION_Z];
        }
        void clear() { dimension.clear(); sliceSize = 0; imageSize = 0; }
        Dimension dimension; ///< the image dimensions [x,y,z,s,t]
        unsigned sliceSize; ///< (x,y) slice size (stored only for performance)
        unsigned imageSize; ///< (x,y,z) image size (stored only for performance)
    };

    BranchingImageDimension *_branchingImageDimension;

public:

    BranchingImage3D* imgList; ///< array of BranchingImage over time t
    bool managingMemory;



    static const char* Name();

    ///constructors/destructors
    BranchingImage() : imgList(0), managingMemory(true)
    {
        _branchingImageDimension = new BranchingImageDimension();
    }

    ~BranchingImage()
    {
        if( managingMemory )
        {
            if( imgList ) delete [] imgList;
            /*if( _branchingImageDimension )*/ delete _branchingImageDimension;
        }
    }


    /// copy constructor
    /// @warning by default it uses shared memory (copy data iff !shared)
    BranchingImage(const BranchingImage<T>& img, bool shared=true)
    {
        if( shared )
        {
            _branchingImageDimension = img._branchingImageDimension;
            imgList = img.imgList;
            managingMemory = false;
        }
        else
        {
            _branchingImageDimension = new BranchingImageDimension();
            imgList = 0;
            managingMemory = true;
            *this = img; // copy data
        }
    }

    /// clone (data is copied)
    BranchingImage<T>& operator=(const BranchingImage<T>& im)
    {
        // allocate & copy everything
        setDimensions( im.getDimension() );

        for( unsigned t=0 ; t<getDimension()[DIMENSION_T] ; ++t )
        {
            imgList[t].clone( im.imgList[t], getDimension()[DIMENSION_S] );
        }

        return *this;
    }


    /// conversion from flat image to connection image
    BranchingImage( const Image<T>& img, BranchingImageConnectivity connectivity )
    {
        this->fromImage( img, connectivity );
    }

    /// conversion from flat image to connection image
    template<class T2>
    void fromImage( const Image<T2>& im, BranchingImageConnectivity connectivity )
    {
        setDimensions( im.getDimensions() );

        for( unsigned t=0 ; t<getDimension()[DIMENSION_T] ; ++t )
        {
            BranchingImage3D& imt = imgList[t];
            const CImg<T2>& cimgt = im.getCImg(t);
            unsigned index1D = 0;
            cimg_forXYZ(cimgt,x,y,z)
            {
                CImg<long double> vect=cimgt.get_vector_at(x,y,z);
                if( vect.magnitude(1) != 0 )
                {
                    //                    assert( index1D == index3Dto1D(x,y,z) );
                    {
                        ConnectionVoxel v( getDimension()[DIMENSION_S] );
                        for( unsigned c = 0 ; c<getDimension()[DIMENSION_S] ; ++c )
                            v[c] = cimgt(x,y,z,c);
                        //                        v.index = index1D;
                        v.neighbours.clear();
                        imt[index1D].push_back( v, getDimension()[DIMENSION_S] );
                    }

                    // neighbours

                    if( connectivity == CONNECTIVITY_6 )
                    {
                        // face
                        if( x>0 && !imt[index1D-1].empty() ) { imt[index1D][0].addNeighbour( VoxelIndex( index1D-1, 0 ), false ); imt[index1D-1][0].addNeighbour( VoxelIndex( index1D, 0 ), false ); }
                        if( y>0 && !imt[index1D-getDimension()[DIMENSION_X]].empty() ) { imt[index1D][0].addNeighbour( VoxelIndex( index1D-getDimension()[DIMENSION_X], 0 ), false ); imt[index1D-getDimension()[DIMENSION_X]][0].addNeighbour( VoxelIndex( index1D, 0 ), false ); }
                        if( z>0 && !imt[index1D-getSliceSize()].empty() ) { imt[index1D][0].addNeighbour( VoxelIndex( index1D-getSliceSize(), 0 ), false ); imt[index1D-getSliceSize()][0].addNeighbour( VoxelIndex( index1D, 0 ), false ); }
                    }
                    else
                    {
                        // neighbours is all directions
                        // TODO could be improved by testing only 3/8 face-, 9/12 edge- and 4/8 corner-neighbours (and so not testing unicity while adding the neighbour)
                        for( int gx = -1 ; gx <= 1 ; ++gx )
                        {
                            if( x+gx<0 || x+gx>(int)getDimension()[DIMENSION_X]-1 ) continue;
                            for( int gy = -1 ; gy <= 1 ; ++gy )
                            {
                                if( y+gy<0 || y+gy>(int)getDimension()[DIMENSION_Y]-1 ) continue;
                                for( int gz = -1 ; gz <= 1 ; ++gz )
                                {
                                    if( z+gz<0 || z+gz>(int)getDimension()[DIMENSION_Z]-1 ) continue;
                                    if( !gx && !gy && !gz ) continue; // do not test with itself

                                    const NeighbourOffset no( gx, gy, gz );

                                    const unsigned neighbourIndex = getNeighbourIndex( no, index1D );

                                    if( !imt[neighbourIndex].empty() )
                                    {
                                        imt[index1D][0].addNeighbour( VoxelIndex( neighbourIndex, 0 ), true );
                                        imt[neighbourIndex][0].addNeighbour( VoxelIndex( index1D, 0 ), true );
                                    }
                                }
                            }
                        }
                    }
                }
                ++index1D;
            }
        }
    }


    /// conversion to a flat image
    /// conversionType : 0->first voxel, 1->average, 2->nb superimposed voxels, 3->sum
    template<class T2>
    void toImage( Image<T2>& img, unsigned conversionType ) const
    {
        img.clear();
        typename Image<T2>::imCoord dim = getDimension();
        img.setDimensions( dim );
        for( unsigned t=0 ; t<getDimension()[DIMENSION_T] ; ++t )
        {
            const BranchingImage3D& bimt = imgList[t];
            CImg<T2>& cimgt = img.getCImg(t);

            cimgt.fill((T2)0);

            unsigned index1D = 0;
            cimg_forXYZ(cimgt,x,y,z)
            {
                cimg_forC( cimgt, c )
                        bimt[index1D].toFlatVoxel( cimgt(x,y,z,c), conversionType, c );

                ++index1D;
            }
        }
    }


    /// delete everything, free memory
    void clear()
    {
        if( managingMemory )
        {
            if( imgList )
            {
                delete [] imgList;
                imgList = 0;
            }
            _branchingImageDimension->clear();
        }
        else
        {
            imgList = 0;
            _branchingImageDimension = new BranchingImageDimension();
            managingMemory = true;
        }
    }


    /// check if image coordinates are inside bounds
    template<class t>
    inline bool isInside( t x, t y, t z ) const
    {
        if(isEmpty()) return false;
        if(x<0) return false;
        if(y<0) return false;
        if(z<0) return false;
        if(x>=(t)getDimension()[DIMENSION_X]) return false;
        if(y>=(t)getDimension()[DIMENSION_Y]) return false;
        if(z>=(t)getDimension()[DIMENSION_Z]) return false;
        return true;
    }

    /// compute the map key in BranchingImage from the pixel position
    inline unsigned index3Dto1D( unsigned x, unsigned y, unsigned z ) const
    {
        return ( z * getDimension()[DIMENSION_Y]  + y ) * getDimension()[DIMENSION_X] + x;
    }

    /// compute the pixel position from the map key in BranchingImage
    inline void index1Dto3D( unsigned key, unsigned& x, unsigned& y, unsigned& z ) const
    {
        //        x = key % getDimension()[DIMENSION_X];
        //        y = ( key / getDimension()[DIMENSION_X] ) % getDimension()[DIMENSION_Y];
        //        z = key / sliceSize;
        y = key / getDimension()[DIMENSION_X];
        x = key - y * getDimension()[DIMENSION_X];
        z = y / getDimension()[DIMENSION_Y];
        y = y - z * getDimension()[DIMENSION_Y];
    }

    /// \returns the index of the neighbour given by its offset (index in the BranchingImage3D)
    inline unsigned getNeighbourIndex( const NeighbourOffset& d, unsigned index1D ) const
    {
        return index1D+d[0]+d[1]*getDimension()[DIMENSION_X]+d[2]*getSliceSize();
    }


    /// \returns the list of connected neighbours of a given voxel
    const Neighbours& getNeighbours(const VoxelIndex& index, const unsigned t=0) const
    {
        const ConnectionVoxel& voxel = this->imgList[t][index.index1d][index.offset];
        return voxel.neighbours;
    }

    /// \returns the list of connected neighbours of a given voxel and the corresponding euclidean distances (image transform supposed to be linear)
    // TO DO: bias the distance using a lut (= multiply distance by lut(v1,v2))
    template<typename real>
    const Neighbours& getNeighboursAndDistances(std::vector< real > &dist, const VoxelIndex& index, const sofa::defaulttype::Vec<3,real>& voxelsize, const unsigned t=0, const bool bias = false) const
    {
        const ConnectionVoxel& voxel = this->imgList[t][index.index1d][index.offset];
        const Neighbours& neighbours = voxel.neighbours;

        dist.resize( neighbours.size() );
        real b1 = bias? (real)operator()(index,0,t) : 1.0 ,b2;

        for( unsigned n = 0 ; n < neighbours.size() ; ++n )
        {
            NeighbourOffset dir = getDirection( index.index1d, neighbours[n].index1d );
            b2 = bias? (real)operator()(neighbours[n],0,t): 1.0;
            real d = Vec3d( abs(dir[0])*voxelsize[0], abs(dir[1])*voxelsize[1], abs(dir[2])*voxelsize[2] ).norm();
            if(d==0) d= (voxelsize[0]+voxelsize[1]+voxelsize[2])/6.; // enforce a certain distance (half voxel size) for superimposed voxels
            dist[n] = d * 1.0/sofa::helper::rmin(b1,b2);
        }
        return voxel.neighbours;
    }

    /// \returns the number of superimposed voxels
    /// @warning validity of indices and time not checked
    unsigned int Nb_superimposed(const unsigned& off1D, const unsigned t=0) const { return this->imgList[t][off1D].size(); }

    /// \returns image value at a given voxel index, time and channel
    /// @warning validity of indices, channel and time not checked
    inline const T& operator()(const unsigned& off1D, const unsigned& v, const unsigned c=0, const unsigned t=0) const  { return this->imgList[t][off1D][v].value[c]; }
    inline       T& operator()(const unsigned& off1D, const unsigned& v, const unsigned c=0, const unsigned t=0)        { return this->imgList[t][off1D][v].value[c]; }
    inline const T& operator()(const VoxelIndex& index, const unsigned c=0, const unsigned t=0)                 const   { return this->imgList[t][index.index1d][index.offset].value[c]; }
    inline       T& operator()(const VoxelIndex& index, const unsigned c=0, const unsigned t=0)                         { return this->imgList[t][index.index1d][index.offset].value[c]; }

    /// \returns the offset between two neighbour voxels
    /// example: returning (-1,0,0) means neighbourIndex is at the LEFT position of index
    inline NeighbourOffset getDirection( unsigned index1D, unsigned neighbourIndex1D ) const
    {
        unsigned X0,Y0,Z0;
        index1Dto3D( index1D, X0,Y0,Z0 );

        for( int x=-1 ; x<=1 ; ++x ) if((int)X0+x>=0) if((int)X0+x<(int)getDimension()[DIMENSION_X])
        for( int y=-1 ; y<=1 ; ++y ) if((int)Y0+y>=0) if((int)Y0+y<(int)getDimension()[DIMENSION_Y])
        for( int z=-1 ; z<=1 ; ++z ) if((int)Z0+z>=0) if((int)Z0+z<(int)getDimension()[DIMENSION_Z])
            if( neighbourIndex1D == index1D+x+y*getDimension()[DIMENSION_X]+z*getSliceSize() ) { return NeighbourOffset(x,y,z); } // connected neighbours

        // not connected neighbours
        // TODO @todo important
        return NeighbourOffset(0,0,0);
    }

    bool isEmpty() const { if( imgList ) return false; else return true;}

    /// \returns the 5 image dimensions (x,y,z,s,t)
    inline const Dimension& getDimension() const
    {
        return _branchingImageDimension->dimension;
    }
    inline const Dimension& getDimensions() const
    {
        return _branchingImageDimension->dimension;
    }
    inline unsigned getSliceSize() const
    {
        return _branchingImageDimension->sliceSize;
    }

    inline unsigned getImageSize() const
    {
        return _branchingImageDimension->imageSize;
    }

    /// resizing
    /// @warning data is deleted
    void setDimensions( const Dimension& newDimension )
    {
        clear();

        for( unsigned i=0 ; i<NB_DimensionLabel ; ++i ) if( !newDimension[i] ) { return; }

        _branchingImageDimension->set( newDimension );
        imgList = new BranchingImage3D[getDimension()[DIMENSION_T]];
        for( unsigned t=0 ; t<getDimension()[DIMENSION_T] ; ++t )
            imgList[t].resize( getImageSize() );
    }

    /// copy only the topology and initialize values based on available DIMENSION_T and DIMENSION_S
    template<typename T2>
    void cloneTopology( const BranchingImage<T2>& other, const T defaultValue=(T)0)
    {
        for( unsigned t=0 ; t<getDimension()[DIMENSION_T] ; ++t )
            imgList[t].template cloneTopology<T2>(other.imgList[t],getDimension()[DIMENSION_S],defaultValue);
    }

    /// read dimensions
    inline friend std::istream& operator >> ( std::istream& in, BranchingImage<T>& im )
    {
        Dimension dim;
        in >> dim;
        im.setDimensions( dim );
        return in;
    }

    /// write dimensions
    friend std::ostream& operator << ( std::ostream& out, const BranchingImage<T>& im )
    {
        out << im.getDimension();
        return out;
    }

    /// comparison
    bool operator==( const BranchingImage<T>& other ) const
    {
        for( unsigned t=0 ; t<getDimension()[DIMENSION_T] ; ++t )
            if( !imgList[t].isEqual( other.imgList[t], getDimension()[DIMENSION_S] ) ) return false;
        return true;
    }

    /// comparison
    bool operator!=( const BranchingImage<T>& other ) const
    {
        return !(*this==other);
    }

    /// count the nb of sumperimposed voxels
    unsigned count() const
    {
        unsigned total = 0;
        for( unsigned t=0 ; t<getDimension()[DIMENSION_T] ; ++t )
        {
            int index1d = 0;
            for( unsigned z=0 ; z<getDimension()[DIMENSION_Z] ; ++z )
                for( unsigned y=0 ; y<getDimension()[DIMENSION_Y] ; ++y )
                    for( unsigned x=0 ; x<getDimension()[DIMENSION_X] ; ++x )
                    {
                        total += imgList[t][index1d].size();
                        ++index1d;
                    }
        }
        return total;
    }


    /// sum every values (every channels)
    T sum() const
    {
        T total = 0;
        for( unsigned t=0 ; t<getDimension()[DIMENSION_T] ; ++t )
        {
            int index1d = 0;
            for( unsigned z=0 ; z<getDimension()[DIMENSION_Z] ; ++z )
                for( unsigned y=0 ; y<getDimension()[DIMENSION_Y] ; ++y )
                    for( unsigned x=0 ; x<getDimension()[DIMENSION_X] ; ++x )
                    {
                        for( unsigned v=0 ; v<imgList[t][index1d].size() ; ++v )
                            for( unsigned s=0 ; s<getDimension()[DIMENSION_S] ; ++s )
                                total += imgList[t][index1d][v][s];
                        ++index1d;
                    }
        }
        return total;
    }

    /// \returns an approximative size in bytes, useful for debugging
    size_t approximativeSizeInBytes() const
    {
        size_t total = getDimension()[DIMENSION_T]*(getImageSize()+1)*( sizeof(unsigned) + sizeof(void*) ); // superimposed voxel vectors + BranchingImage3D vector

        for( unsigned t=0 ; t<getDimension()[DIMENSION_T] ; ++t ) // per BranchingImage3D
        {
            const BranchingImage3D& imt = imgList[t];
            for( unsigned int index=0 ; index<imt.size() ; ++index ) // per SumperimposedVoxels
            {
                total += imt[index].size() * ( /*sizeof(unsigned)*/ /*index*/ +
                                               sizeof(void*) /* channel vector*/ +
                                               getDimension()[DIMENSION_S]*sizeof(T) /*channel entries*/
                                               );

                for( unsigned v=0 ; v<imt[index].size() ; ++v ) // per ConnnectedVoxel
                {
                    total += imt[index][v].neighbours.size() * 2 * sizeof(unsigned); //neighbours
                }
            }
        }
        return total;
    }

    /// check neighbourhood validity
    int isNeighbourhoodValid() const
    {
        for( unsigned t=0 ; t<getDimension()[DIMENSION_T] ; ++t )
        {
            const BranchingImage3D& imt = imgList[t];
            unsigned index1d = 0;
            for( unsigned z = 0 ; z < getDimension()[DIMENSION_Z] ; ++z )
            {
                for( unsigned y = 0 ; y < getDimension()[DIMENSION_Y] ; ++y )
                {
                    for( unsigned x = 0 ; x < getDimension()[DIMENSION_X] ; ++x )
                    {
                        const SuperimposedVoxels& voxels = imt[index1d];
                        for( unsigned v = 0 ; v < voxels.size() ; ++v )
                        {
                            const Neighbours& neighbours = voxels[v].neighbours;
                            //if( neighbours.empty() ) return 3;

                            for( unsigned n = 0 ; n < neighbours.size() ; ++n )
                            {
                                unsigned neighbourIndex = neighbours[n].index1d;
                                unsigned neighbourOffset = neighbours[n].offset;
                                if( neighbourOffset >= imt[neighbourIndex].size() )  return 1; // there is nobody where there should be the neighbour
                                if( !imt[neighbourIndex][neighbourOffset].isNeighbour( VoxelIndex( index1d, imt.getOffset(index1d,voxels[v]) ) ) ) return 2; // complementary neighbour is no inserted
                            }
                        }
                        ++ index1d;
                    }
                }
            }
        }
        return 0;
    }

    /// \returns 0 iff flatImg==*this otherwise returns an error code
    template<class T2>
    int isEqual( const Image<T2>& flatImg, BranchingImageConnectivity connectivity, bool valueTest = true, bool neighbourTest = false ) const
    {
        cimglist_for(flatImg.getCImgList(),l)
        {
              const CImg<T2>& cimgl = flatImg.getCImg(l);
              const BranchingImage3D& iml = this->imgList[l];
              unsigned index1d = -1;
              cimg_forXYZ(cimgl,x,y,z)
              {
                  ++index1d;

                  const SuperimposedVoxels& voxels = iml[index1d];

                  if( voxels.empty() ) //the pixel x,y,z is not present in the branching image
                  {
                      if ( cimgl.get_vector_at(x,y,z).magnitude(1)!=0 ) return 1; // if the pixel is present in the flat image, there is a pb
                      else continue; // no pixel -> nothing to compare, test the next pixel
                  }

                  if( voxels.size()>1 ) return 2; // the branching image has been built from a flat image, so there should be no superimposed voxels

                  if( valueTest )
                  {
                      for( unsigned c=0 ; c<flatImg.getDimensions()[3] ; ++c ) // for all channels
                          if( (T)cimgl(x,y,z,c) != voxels[0][c] ) return 3; // check that the value is the same
                  }

                  if( neighbourTest )
                  {
                      if( connectivity == CONNECTIVITY_6 )
                      {
                          // test neighbourhood connections
                          if( x>0 && ( ( cimgl.get_vector_at(x-1,y,z).magnitude(1)==0 ) == ( voxels[0].isNeighbour( VoxelIndex(index1d-1,0) ) ) ) ) return 4;
                          if( (unsigned)x<flatImg.getDimensions()[0]-1 && ( ( cimgl.get_vector_at(x+1,y,z).magnitude(1)==0 ) == ( voxels[0].isNeighbour( VoxelIndex(index1d+1,0) ) ) ) ) return 4;
                          if( y>0 && ( ( cimgl.get_vector_at(x,y-1,z).magnitude(1)==0 ) == ( voxels[0].isNeighbour( VoxelIndex(index1d-getDimension()[DIMENSION_X],0) ) ) ) ) return 4;
                          if( (unsigned)y<flatImg.getDimensions()[1]-1 && ( ( cimgl.get_vector_at(x,y+1,z).magnitude(1)==0 ) == ( voxels[0].isNeighbour( VoxelIndex(index1d+getDimension()[DIMENSION_X],0) ) ) ) ) return 4;
                          if( z>0 && ( ( cimgl.get_vector_at(x,y,z-1).magnitude(1)==0 ) == ( voxels[0].isNeighbour( VoxelIndex(index1d-getSliceSize(),0) ) ) ) ) return 4;
                          if( (unsigned)z<flatImg.getDimensions()[2]-1 && ( ( cimgl.get_vector_at(x,y,z+1).magnitude(1)==0 ) == ( voxels[0].isNeighbour( VoxelIndex(index1d+getSliceSize(),0) ) ) ) ) return 4;
                      }
                      else // CONNECTIVITY_26
                      {
                          for( int gx = -1 ; gx <= 1 ; ++gx )
                          {
                              if( x+gx<0 || x+gx>(int)getDimension()[DIMENSION_X]-1 ) continue;
                              for( int gy = -1 ; gy <= 1 ; ++gy )
                              {
                                  if( y+gy<0 || y+gy>(int)getDimension()[DIMENSION_Y]-1 ) continue;
                                  for( int gz = -1 ; gz <= 1 ; ++gz )
                                  {
                                      if( z+gz<0 || z+gz>(int)getDimension()[DIMENSION_Z]-1 ) continue;
                                      if( !gx && !gy && !gz ) continue; // do not test with itself

                                      const NeighbourOffset no( gx, gy, gz );

                                      const unsigned neighbourIndex = getNeighbourIndex( no, index1d );

                                      if( ( cimgl.get_vector_at(x+gx,y+gy,z+gz).magnitude(1)==0 ) == ( voxels[0].isNeighbour( VoxelIndex(neighbourIndex,0) ) ) ) return 4;
                                  }
                              }
                          }
                      }
                  }
              }
        }
        return 0;
    }




    template<typename F>
    bool save( const char *const headerFilename, const F *const scale=0, const F *const translation=0, const F *const affine=0, F offsetT=0, F scaleT=0, bool isPerspective=false ) const
    {
        if( !getDimension()[DIMENSION_T] ) return false;

        std::ofstream fileStream (headerFilename, std::ofstream::out);
        if (!fileStream.is_open())	{	std::cout << "Can not open " << headerFilename << std::endl;	return false; }

        fileStream << "ObjectType = BranchingImage" << std::endl;

        unsigned int nbdims=(getDimension()[DIMENSION_T]==1)?3:4; //  for 2-d, we still need z scale dimension

        fileStream << "NDims = " << nbdims << std::endl;

        fileStream << "ElementNumberOfChannels = " << getDimension()[DIMENSION_S] << std::endl;

        fileStream << "DimSize = "; for(unsigned int i=0;i<nbdims;i++) fileStream << getDimension()[i] << " "; fileStream << std::endl;

        fileStream << "ElementType = ";
        if(!strcmp(cimg::type<T>::string(),"char")) fileStream << "MET_CHAR" << std::endl;
        else if(!strcmp(cimg::type<T>::string(),"double")) fileStream << "MET_DOUBLE" << std::endl;
        else if(!strcmp(cimg::type<T>::string(),"float")) fileStream << "MET_FLOAT" << std::endl;
        else if(!strcmp(cimg::type<T>::string(),"int")) fileStream << "MET_INT" << std::endl;
        else if(!strcmp(cimg::type<T>::string(),"long")) fileStream << "MET_LONG" << std::endl;
        else if(!strcmp(cimg::type<T>::string(),"short")) fileStream << "MET_SHORT" << std::endl;
        else if(!strcmp(cimg::type<T>::string(),"unsigned char")) fileStream << "MET_UCHAR" << std::endl;
        else if(!strcmp(cimg::type<T>::string(),"unsigned int")) fileStream << "MET_UINT" << std::endl;
        else if(!strcmp(cimg::type<T>::string(),"unsigned long")) fileStream << "MET_ULONG" << std::endl;
        else if(!strcmp(cimg::type<T>::string(),"unsigned short")) fileStream << "MET_USHORT" << std::endl;
        else if(!strcmp(cimg::type<T>::string(),"bool")) fileStream << "MET_BOOL" << std::endl;
        else fileStream << "MET_UNKNOWN" << std::endl;

        if(scale) { fileStream << "ElementSpacing = "; for(unsigned int i=0;i<3;i++) fileStream << scale[i] << " "; if(nbdims==4) fileStream << scaleT; fileStream << std::endl; }

        if(translation) { fileStream << "Position = "; for(unsigned int i=0;i<3;i++) fileStream << translation[i] << " "; if(nbdims==4) fileStream << offsetT; fileStream << std::endl; }

        if(affine) { fileStream << "Orientation = "; for(unsigned int i=0;i<9;i++) fileStream << affine[i] << " "; fileStream << std::endl; }

        fileStream << "isPerpective = " << isPerspective << std::endl;

        std::string imageFilename(headerFilename); imageFilename.replace(imageFilename.find_last_of('.')+1,imageFilename.size(),"bia");
        fileStream << "ElementDataFile = " << imageFilename.c_str() << std::endl;
        fileStream.close();


        fileStream.open( imageFilename.c_str(), std::ofstream::out );
        if( !fileStream.is_open() ) { std::cout << "Can not open " << imageFilename << std::endl; return false; }

        // each line corresponds to a superimposed voxel list
        // each superimposed voxel is separated by #
        // the end of list is given by }
        // a superimposed has a double for each channels and a list of neighbours stored as 2 unsigned (index, offset)

        for( unsigned t=0 ; t<getDimension()[DIMENSION_T] ; ++t )
        {
            unsigned index1d = 0;
            for( unsigned x=0 ; x<getDimension()[DIMENSION_X] ; ++x )
            for( unsigned y=0 ; y<getDimension()[DIMENSION_Y] ; ++y )
            for( unsigned z=0 ; z<getDimension()[DIMENSION_Z] ; ++z )
            {
                for( unsigned v=0 ; v<imgList[t][index1d].size() ; ++v )
                {
                    for( unsigned s=0 ; s<getDimension()[DIMENSION_S] ; ++s )
                        fileStream << (double)imgList[t][index1d][v][s] << " ";

                    for( unsigned n=0 ; n<imgList[t][index1d][v].neighbours.size() ; ++n )
                    {
                        const VoxelIndex& vi = imgList[t][index1d][v].neighbours[n];
                        fileStream << vi.index1d << " " << vi.offset << " ";
                    }

                    if( v!=imgList[t][index1d].size()-1 ) fileStream << "#";
                }
                fileStream << "}";
                fileStream << std::endl;
                ++index1d;
            }
        }

        fileStream.close();
        return true;
    }


    template<typename F>
    bool load( const char *const headerFilename, F *const scale=0, F *const translation=0, F *const affine=0, F *const offsetT=0, F *const scaleT=0, bool *const isPerspective=0 )
    {
        std::ifstream fileStream(headerFilename, std::ifstream::in);
        if (!fileStream.is_open())	{ std::cout << "Can not open " << headerFilename << std::endl; clear(); return false; }

        std::string str,str2,imageFilename;
        unsigned int nbchannels=1,nbdims=4,dim[] = {1,1,1,1}; // 3 spatial dimas + time
        std::string inputType(cimg::type<T>::string());
        while(!fileStream.eof())
        {
            fileStream >> str;

            if(!str.compare("ObjectType"))
            {
                fileStream >> str2; // '='
                fileStream >> str2;
                if(str2.compare("BranchingImage")) { std::cout << "BranchingImage::load: not a BranchingImage ObjectType "<<std::endl; clear(); return false;}
            }
            else if(!str.compare("ElementDataFile"))
            {
                fileStream >> str2; // '='
                fileStream >> imageFilename;
            }
            else if(!str.compare("NDims"))
            {
                fileStream >> str2;  // '='
                fileStream >> nbdims;
                if(nbdims>4) { std::cout << "BranchingImage::load: dimensions > 4 not supported  "<<std::endl; clear(); return false;}
            }
            else if(!str.compare("ElementNumberOfChannels"))
            {
                fileStream >> str2;  // '='
                fileStream >> nbchannels;
            }
            else if(!str.compare("DimSize") || !str.compare("dimensions") || !str.compare("dim"))
            {
                fileStream >> str2;  // '='
                for(unsigned int i=0;i<nbdims;i++) fileStream >> dim[i];
            }
            else if(!str.compare("ElementSpacing") || !str.compare("spacing") || !str.compare("scale3d") || !str.compare("voxelSize"))
            {
                fileStream >> str2; // '='
                double val[4];
                for(unsigned int i=0;i<nbdims;i++) fileStream >> val[i];
                if(scale) for(unsigned int i=0;i<3;i++) if(i<nbdims) scale[i] = (F)val[i];
                if(scaleT) if(nbdims>3) *scaleT = (F)val[3];
           }
            else if(!str.compare("Position") || !str.compare("Offset") || !str.compare("translation") || !str.compare("origin"))
            {
                fileStream >> str2; // '='
                double val[4];
                for(unsigned int i=0;i<nbdims;i++) fileStream >> val[i];
                if(translation) for(unsigned int i=0;i<3;i++) if(i<nbdims) translation[i] = (F)val[i];
                if(offsetT) if(nbdims>3) *offsetT = (F)val[3];
            }
            else if(!str.compare("Orientation"))
            {
                fileStream >> str2; // '='
                double val[4*4];
                for(unsigned int i=0;i<nbdims*nbdims;i++) fileStream >> val[i];
                if(affine) { for(unsigned int i=0;i<3;i++) if(i<nbdims) for(unsigned int j=0;j<3;j++) if(j<nbdims) affine[i*3+j] = (F)val[i*nbdims+j]; }
                // to do: handle "CenterOfRotation" Tag
            }
            else if(!str.compare("isPerpective")) { fileStream >> str2; bool val; fileStream >> val; if(isPerspective) *isPerspective=val; }
            else if(!str.compare("ElementType") || !str.compare("voxelType"))  // not used (should be known in advance for template)
            {
                fileStream >> str2; // '='
                fileStream >> str2;

                if(!str2.compare("MET_CHAR"))           inputType=std::string("char");
                else if(!str2.compare("MET_DOUBLE"))    inputType=std::string("double");
                else if(!str2.compare("MET_FLOAT"))     inputType=std::string("float");
                else if(!str2.compare("MET_INT"))       inputType=std::string("int");
                else if(!str2.compare("MET_LONG"))      inputType=std::string("long");
                else if(!str2.compare("MET_SHORT"))     inputType=std::string("short");
                else if(!str2.compare("MET_UCHAR"))     inputType=std::string("unsigned char");
                else if(!str2.compare("MET_UINT"))      inputType=std::string("unsigned int");
                else if(!str2.compare("MET_ULONG"))     inputType=std::string("unsigned long");
                else if(!str2.compare("MET_USHORT"))    inputType=std::string("unsigned short");
                else if(!str2.compare("MET_BOOL"))      inputType=std::string("bool");

                if(inputType!=std::string(cimg::type<T>::string()))
                {
                    std::cout<<"BranchingImage::load: BranchingImage type ( "<< str2 <<" ) is different from "<< cimg::type<T>::string() <<" ) - no conversion available"<<std::endl;
                    clear();
                    return false;
                }
            }
        }
        fileStream.close();

        if(!imageFilename.size()) // no specified file name -> replace .mhd by .bia
        {
            imageFilename = std::string(headerFilename);
            imageFilename .replace(imageFilename.find_last_of('.')+1,imageFilename.size(),"bia");
        }
        else // add path to the specified file name
        {
            std::string tmp(headerFilename);
            std::size_t pos=tmp.find_last_of('/');
            if(pos==std::string::npos) pos=tmp.find_last_of('\\');
            if(pos!=std::string::npos) {tmp.erase(pos+1); imageFilename.insert(0,tmp);}
        }

        fileStream.open( imageFilename.c_str(), std::ifstream::in );
        if( !fileStream.is_open() ) { std::cout << "Can not open " << imageFilename << std::endl; clear(); return false; }

        setDimensions( Dimension( dim[0], dim[1], dim[2], nbchannels, dim[3] ) );

        for( unsigned t=0 ; t<getDimension()[DIMENSION_T] ; ++t )
        {
            unsigned index1d = 0;
            for( unsigned x=0 ; x<getDimension()[DIMENSION_X] ; ++x )
            for( unsigned y=0 ; y<getDimension()[DIMENSION_Y] ; ++y )
            for( unsigned z=0 ; z<getDimension()[DIMENSION_Z] ; ++z )
            {

                char symbol;
                fileStream >> symbol;
                while( symbol != '}' ) // end of superimposed voxels
                {
                    if( symbol != '#' ) fileStream.seekg( -1, std::ifstream::cur ); // unread symbol

                    ConnectionVoxel cv( getDimension()[DIMENSION_S] );

                    for( unsigned s=0 ; s<getDimension()[DIMENSION_S] ; ++s )
                    {
                        double c;
                        fileStream >> c;
                        cv[s] = (T)c;

//                        std::cerr<<(T)c<<",";
                    }
//                    std::cerr<<"-";
                    fileStream >> symbol;
                    while( symbol != '#' && symbol != '}' ) // end of neighbours
                    {
                        fileStream.seekg( -1, std::ifstream::cur ); // unread symbol

                        unsigned index, offset;
                        fileStream >> index;
                        fileStream >> offset;
                        cv.addNeighbour( VoxelIndex( index, offset ) );

//                        std::cerr<<index<<"/"<<offset<<",";

                        fileStream >> symbol;
                    }

                    imgList[t][index1d].push_back( cv, getDimension()[DIMENSION_S] );

//                    std::cerr<<"#";
                }

//                std::cerr<<std::endl;

                ++index1d;
            }
        }

        fileStream.close();

        return true;
    }


};




/// macros for image loops

#define bimg_for1(bound,i) for (unsigned i = 0; i<bound; ++i)
#define bimg_forC(img,c) bimg_for1(img.getDimension()[img.DIMENSION_S],c)
#define bimg_forX(img,x) bimg_for1(img.getDimension()[img.DIMENSION_X],x)
#define bimg_forY(img,y) bimg_for1(img.getDimension()[img.DIMENSION_Y],y)
#define bimg_forZ(img,z) bimg_for1(img.getDimension()[img.DIMENSION_Z],z)
#define bimg_forT(img,t) bimg_for1(img.getDimension()[img.DIMENSION_T],t)
#define bimg_foroff1D(img,off1D)  bimg_for1(img.getImageSize(),off1D)
#define bimg_forXY(img,x,y) bimg_forY(img,y) bimg_forX(img,x)
#define bimg_forXZ(img,x,z) bimg_forZ(img,z) bimg_forX(img,x)
#define bimg_forYZ(img,y,z) bimg_forZ(img,z) bimg_forY(img,y)
#define bimg_forXT(img,x,t) bimg_forT(img,t) bimg_forX(img,x)
#define bimg_forYT(img,y,t) bimg_forT(img,t) bimg_forY(img,y)
#define bimg_forZT(img,z,t) bimg_forT(img,t) bimg_forZ(img,z)
#define bimg_forXYZ(img,x,y,z) bimg_forZ(img,z) bimg_forXY(img,x,y)
#define bimg_forXYT(img,x,y,t) bimg_forT(img,t) bimg_forXY(img,x,y)
#define bimg_forXZT(img,x,z,t) bimg_forT(img,t) bimg_forXZ(img,x,z)
#define bimg_forYZT(img,y,z,t) bimg_fort(img,t) bimg_forYZ(img,y,z)
#define bimg_forXYZT(img,x,y,z,t) bimg_forT(img,t) bimg_forXYZ(img,x,y,z)
#define bimg_forVXYZT(img,v,x,y,z,t) bimg_forT(img,t) bimg_forXYZ(img,x,y,z) for( unsigned v=0 ; v<img.imgList[t][img.index3Dto1D(x,y,z)].size() ; ++v )
#define bimg_forCVXYZT(img,c,v,x,y,z,t) bimg_forT(img,t) bimg_forXYZ(img,x,y,z) for( unsigned v=0 ; v<img.imgList[t][img.index3Dto1D(x,y,z)].size() ; ++v ) bimg_forC(img,c)
#define bimg_forVoffT(img,v,off1D,t) bimg_forT(img,t) bimg_foroff1D(img,off1D) for( unsigned v=0 ; v<img.imgList[t][off1D].size() ; ++v )
#define bimg_forCVoffT(img,c,v,off1D,t) bimg_forT(img,t) bimg_foroff1D(img,off1D) for( unsigned v=0 ; v<img.imgList[t][off1D].size() ; ++v )  bimg_forC(img,c)


typedef BranchingImage<char> BranchingImageC;
typedef BranchingImage<unsigned char> BranchingImageUC;
typedef BranchingImage<int> BranchingImageI;
typedef BranchingImage<unsigned int> BranchingImageUI;
typedef BranchingImage<short> BranchingImageS;
typedef BranchingImage<unsigned short> BranchingImageUS;
typedef BranchingImage<long> BranchingImageL;
typedef BranchingImage<unsigned long> BranchingImageUL;
typedef BranchingImage<float> BranchingImageF;
typedef BranchingImage<double> BranchingImageD;
typedef BranchingImage<bool> BranchingImageB;

template<> inline const char* BranchingImageC::Name() { return "BranchingImageC"; }
template<> inline const char* BranchingImageUC::Name() { return "BranchingImageUC"; }
template<> inline const char* BranchingImageI::Name() { return "BranchingImageI"; }
template<> inline const char* BranchingImageUI::Name() { return "BranchingImageUI"; }
template<> inline const char* BranchingImageS::Name() { return "BranchingImageS"; }
template<> inline const char* BranchingImageUS::Name() { return "BranchingImageUS"; }
template<> inline const char* BranchingImageL::Name() { return "BranchingImageL"; }
template<> inline const char* BranchingImageUL::Name() { return "BranchingImageUL"; }
template<> inline const char* BranchingImageF::Name() { return "BranchingImageF"; }
template<> inline const char* BranchingImageD::Name() { return "BranchingImageD"; }
template<> inline const char* BranchingImageB::Name() { return "BranchingImageB"; }

// The next line hides all those methods from the doxygen documentation
/// \cond TEMPLATE_OVERRIDES

template<> struct DataTypeName< defaulttype::BranchingImageC > { static const char* name() { return "BranchingImageC"; } };
template<> struct DataTypeName< defaulttype::BranchingImageUC > { static const char* name() { return "BranchingImageUC"; } };
template<> struct DataTypeName< defaulttype::BranchingImageI > { static const char* name() { return "BranchingImageI"; } };
template<> struct DataTypeName< defaulttype::BranchingImageUI > { static const char* name() { return "BranchingImageUI"; } };
template<> struct DataTypeName< defaulttype::BranchingImageS > { static const char* name() { return "BranchingImageS"; } };
template<> struct DataTypeName< defaulttype::BranchingImageUS > { static const char* name() { return "BranchingImageUS"; } };
template<> struct DataTypeName< defaulttype::BranchingImageL > { static const char* name() { return "BranchingImageL"; } };
template<> struct DataTypeName< defaulttype::BranchingImageUL > { static const char* name() { return "BranchingImageUL"; } };
template<> struct DataTypeName< defaulttype::BranchingImageF > { static const char* name() { return "BranchingImageF"; } };
template<> struct DataTypeName< defaulttype::BranchingImageD > { static const char* name() { return "BranchingImageD"; } };
template<> struct DataTypeName< defaulttype::BranchingImageB > { static const char* name() { return "BranchingImageB"; } };

/// \endcond


} // namespace defaulttype


} // namespace sofa


#endif // IMAGE_BranchingImage_H

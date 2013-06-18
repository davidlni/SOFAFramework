#ifndef BASEINTTOOL_H
#define BASEINTTOOL_H


#include <sofa/core/collision/Intersection.h>
#include <sofa/helper/FnDispatcher.h>
#include <sofa/component/collision/CapsuleModel.h>
#include <sofa/component/collision/SphereModel.h>
#include <sofa/component/collision/OBBModel.h>
#include <sofa/component/collision/IntrCapsuleOBB.h>
#include <sofa/component/collision/CapsuleIntTool.h>
#include <sofa/component/collision/OBBIntTool.h>
#include <cmath>

namespace sofa
{
namespace component
{
namespace collision
{

class SOFA_BASE_COLLISION_API BaseIntTool : public CapsuleIntTool,public OBBIntTool{
public:
    typedef sofa::helper::vector<sofa::core::collision::DetectionOutput> OutputVector;

    static bool testIntersection(Cube& ,Cube&,double alarmDist);

    template <class DataTypes1,class DataTypes2>
    static bool testIntersection(TSphere<DataTypes1>& sph1, TSphere<DataTypes2>& sph2,double alarmDist)
    {
        typename TSphere<DataTypes1>::Real r = sph1.r() + sph2.r() + alarmDist;
        return ( sph1.center() - sph2.center() ).norm2() <= r*r;
    }


    template <class DataTypes1,class DataTypes2>
    static int computeIntersection(TSphere<DataTypes1>& sph1, TSphere<DataTypes2>& sph2,double alarmDist,double contactDist,OutputVector* contacts);


    inline static int computeIntersection(Capsule &c1, Capsule &c2, double alarmDist, double contactDist, OutputVector *contacts){
        return CapsuleIntTool::computeIntersection(c1,c2,alarmDist,contactDist,contacts);
    }

    template <class DataType>
    inline static int computeIntersection(Capsule &cap, TSphere<DataType> &sph, double alarmDist, double contactDist, OutputVector *contacts){
        return CapsuleIntTool::computeIntersection(cap,sph,alarmDist,contactDist,contacts);
    }

    inline static int computeIntersection(Capsule &cap, OBB & obb, double alarmDist, double contactDist, OutputVector *contacts){
        return CapsuleIntTool::computeIntersection(cap,obb,alarmDist,contactDist,contacts);
    }

    inline static int computeIntersection(OBB &obb0, OBB &obb1, double alarmDist, double contactDist, OutputVector *contacts){
        return OBBIntTool::computeIntersection(obb0,obb1,alarmDist,contactDist,contacts);
    }

    template <class DataType>
    inline static int computeIntersection(TSphere<DataType> &sph, OBB &obb, double alarmDist, double contactDist, OutputVector *contacts){
        return OBBIntTool::computeIntersection(sph,obb,alarmDist,contactDist,contacts);
    }
};

template <class DataTypes1,class DataTypes2>
int BaseIntTool::computeIntersection(TSphere<DataTypes1>& sph1, TSphere<DataTypes2>& sph2, double alarmDist, double contactDist, OutputVector * contacts)
{
    double r = sph1.r() + sph2.r();
    double myAlarmDist = alarmDist + r;
    Vector3 dist = sph2.center() - sph1.center();
    SReal norm2 = dist.norm2();

    if (norm2 >= myAlarmDist*myAlarmDist)
        return 0;

    contacts->resize(contacts->size()+1);
    DetectionOutput *detection = &*(contacts->end()-1);
    detection->normal = dist;
    double distSph1Sph2 = helper::rsqrt(norm2);
    detection->normal /= distSph1Sph2;
    detection->point[0] = sph1.center() + detection->normal * sph1.r();
    detection->point[1] = sph2.center() - detection->normal * sph2.r();

    detection->value = distSph1Sph2 - r - contactDist;
    detection->elem.first = sph1;
    detection->elem.second = sph2;
    detection->id = (sph1.getCollisionModel()->getSize() > sph2.getCollisionModel()->getSize()) ? sph1.getIndex() : sph2.getIndex();

    return 1;
}

}
}
}
#endif // BASEINTTOOL_H

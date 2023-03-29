/*
    imager.h

    Copyright (C) 2013 by Don Cross  -  http://cosinekitty.com/raytrace

    This software is provided 'as-is', without any express or implied
    warranty. In no event will the author be held liable for any damages
    arising from the use of this software.

    Permission is granted to anyone to use this software for any purpose,
    including commercial applications, and to alter it and redistribute it
    freely, subject to the following restrictions:

    1. The origin of this software must not be misrepresented; you must not
       claim that you wrote the original software. If you use this software
       in a product, an acknowledgment in the product documentation would be
       appreciated but is not required.

    2. Altered source versions must be plainly marked as such, and must not be
       misrepresented as being the original software.

    3. This notice may not be removed or altered from any source
       distribution.
*/

#ifndef __DDC_IMAGER_H
#define __DDC_IMAGER_H

#include <vector>
#include <cmath>
#include "algebra.h"

namespace Imager
{
    const double Pi = 3.141592653589793238462643383279502884;

    // EPSILON is a tolerance value for floating point roundoff error.
    // It is used in many calculations where we want to err
    // on a certain side of a threshold, such as determining
    // whether or not a point is inside a solid or not,
    // or whether a point is at least a minimum distance
    // away from another point.
    const double EPSILON = 1.0e-6;

    inline double RadiansFromDegrees(double degrees)
    {
        return degrees * (Pi / 180.0);
    }

    //------------------------------------------------------------------------
    // Forward declarations
    class SolidObject;
    class ImageBuffer;

    //------------------------------------------------------------------------

    // An exception thrown by imager code when a fatal error occurs.
    class ImagerException
    {
    public:
        explicit ImagerException(const char *_message)
            : message(_message)
        {
        }

        const char *GetMessage() const { return message; }

    private:
        const char * const message;
    };

    //------------------------------------------------------------------------

    // An exception thrown when multiple intersections lie at the
    // same distance from the vantage point.  SaveImage catches
    // these and marks such pixels as ambiguous.  It performs a second
    // pass later that averages the color values of surrounding
    // non-ambiguous pixels.
    class AmbiguousIntersectionException
    {
    };

    //------------------------------------------------------------------------

    class Vector
    {
    public:
        double x;
        double y;
        double z;

        // Default constructor: create a vector whose
        // x, y, z components are all zero.
        Vector()
            : x(0.0)
            , y(0.0)
            , z(0.0)
        {
        }

        // This constructor initializes a vector
        // to any desired component values.
        Vector(double _x, double _y, double _z)
            : x(_x)
            , y(_y)
            , z(_z)
        {
        }

        // Returns the square of the magnitude of this vector.
        // This is more efficient than computing the magnitude itself,
        // and is just as good for comparing two vectors to see which
        // is longer or shorter.
        const double MagnitudeSquared() const
        {
            return (x*x) + (y*y) + (z*z);
        }

        const double Magnitude() const
        {
            return sqrt(MagnitudeSquared());
        }

        const Vector UnitVector() const
        {
            const double mag = Magnitude();
            return Vector(x/mag, y/mag, z/mag);
        }

        Vector& operator *= (const double factor)
        {
            x *= factor;
            y *= factor;
            z *= factor;
            return *this;
        }

        Vector& operator += (const Vector& other)
        {
            x += other.x;
            y += other.y;
            z += other.z;
            return *this;
        }
    };


    //------------------------------------------------------------------------

    inline Vector operator + (const Vector &a, const Vector &b)
    {
        return Vector(a.x + b.x, a.y + b.y, a.z + b.z);
    }

    inline Vector operator - (const Vector &a, const Vector &b)
    {
        return Vector(a.x - b.x, a.y - b.y, a.z - b.z);
    }

    inline Vector operator - (const Vector& a)
    {
        return Vector(-a.x, -a.y, -a.z);
    }

    inline double DotProduct (const Vector& a, const Vector& b)
    {
        return (a.x*b.x) + (a.y*b.y) + (a.z*b.z);
    }

    inline Vector CrossProduct (const Vector& a, const Vector& b)
    {
        return Vector(
            (a.y * b.z) - (a.z * b.y),
            (a.z * b.x) - (a.x * b.z),
            (a.x * b.y) - (a.y * b.x));
    }

    inline Vector operator * (double s, const Vector& v)
    {
        return Vector(s*v.x, s*v.y, s*v.z);
    }

    inline Vector operator / (const Vector& v, double s)
    {
        return Vector(v.x/s, v.y/s, v.z/s);
    }

    //------------------------------------------------------------------------

    struct Color
    {
        double  red;
        double  green;
        double  blue;

        Color(double _red, double _green, double _blue, double _luminosity = 1.0)
            : red  (_luminosity * _red)
            , green(_luminosity * _green)
            , blue (_luminosity * _blue)
        {
        }

        Color()
            : red(0.0)
            , green(0.0)
            , blue(0.0)
        {
        }

        Color& operator += (const Color& other)
        {
            red   += other.red;
            green += other.green;
            blue  += other.blue;
            return *this;
        }

        Color& operator *= (const Color& other)
        {
            red   *= other.red;
            green *= other.green;
            blue  *= other.blue;
            return *this;
        }

        Color& operator *= (double factor)
        {
            red   *= factor;
            green *= factor;
            blue  *= factor;
            return *this;
        }

        Color& operator /= (double denom)
        {
            red   /= denom;
            green /= denom;
            blue  /= denom;
            return *this;
        }

        void Validate() const
        {
            if ((red < 0.0) || (green < 0.0) || (blue < 0.0))
            {
                throw ImagerException("Negative color values not allowed.");
            }
        }
    };

    inline Color operator * (const Color& aColor, const Color& bColor)
    {
        return Color(
            aColor.red   * bColor.red,
            aColor.green * bColor.green,
            aColor.blue  * bColor.blue);
    }

    inline Color operator * (double scalar, const Color &color)
    {
        return Color(
            scalar * color.red,
            scalar * color.green,
            scalar * color.blue);
    }

    inline Color operator + (const Color& a, const Color& b)
    {
        return Color(
            a.red   + b.red,
            a.green + b.green,
            a.blue  + b.blue);
    }

    //------------------------------------------------------------------------
    // struct Intersection provides information about a ray intersecting
    // with a point on the surface of a SolidObject.

    struct Intersection
    {
        // The square of the distance from the
        // vantage point to the intersection point.
        double distanceSquared;

        // The location of the intersection point.
        Vector point;

        // The unit vector perpendicular to the
        // surface at the intersection point.
        Vector surfaceNormal;

        // A pointer to the solid object that the ray
        // intersected with.
        const SolidObject* solid;

        // An optional tag for classes derived from SolidObject to cache
        // arbitrary information about surface optics.  Most classes can
        // safely leave this pointer as nullptr, its default value.
        const void* context;

        // An optional tag used for debugging.
        // Anything that finds an intersection may elect to make tag point
        // at a static string to help the programmer figure out, for example,
        // which of multiple surfaces was involved.  This is just a char*
        // instead of std::string to minimize overhead by eliminating dynamic
        // memory allocation.
        const char* tag;

        // This constructor initializes to deterministic values
        // in case some code forgets to initialize something.
        Intersection()
            : distanceSquared(1.0e+20)  // larger than any reasonable value
            , point()
            , surfaceNormal()
            , solid(nullptr)
            , context(nullptr)
            , tag(nullptr)
        {
        }
    };

    typedef std::vector<Intersection> IntersectionList;

    int PickClosestIntersection(
        const IntersectionList& list,
        Intersection& intersection);

    //------------------------------------------------------------------------

    class Taggable       // helps debugging; allows caller to assign names to things
    {
    public:
        Taggable(std::string _tag = "")
            : tag(_tag)
        {
        }

        void SetTag(std::string _tag)
        {
            tag = _tag;
        }

        std::string GetTag() const
        {
            return tag;
        }

    private:
        std::string tag;
    };

    //------------------------------------------------------------------------

    class SolidObject: public Taggable
    {
    public:
        SolidObject(const Vector& _center = Vector(), bool _isFullyEnclosed = true)
            : center(_center)
            , isFullyEnclosed(_isFullyEnclosed)
        {
        }

        virtual ~SolidObject()
        {
        }

        // Appends to 'intersectionList' all the
        // intersections found starting at the specified vantage
        // point in the direction of the direction vector.
        virtual void AppendAllIntersections(
            const Vector& vantage,
            const Vector& direction,
            IntersectionList& intersectionList) const = 0;

        // Searches for any intersections with this solid from the
        // vantage point in the given direction.  If none are found, the
        // function returns 0 and the 'intersection' parameter is left
        // unchanged.  Otherwise, returns the positive number of
        // intersections that lie at minimal distance from the vantage point
        // in that direction.  Usually this number will be 1 (a unique
        // intersection is closer than all the others) but it can be greater
        // if multiple intersections are equally close (e.g. the ray hitting
        // exactly at the corner of a cube could cause this function to
        // return 3).  If this function returns a value greater than zero,
        // it means the 'intersection' parameter has been filled in with the
        // closest intersection (or one of the equally closest intersections).
        int FindClosestIntersection(
            const Vector& vantage,
            const Vector& direction,
            Intersection &intersection) const
        {
            cachedIntersectionList.clear();
            AppendAllIntersections(vantage, direction, cachedIntersectionList);
            return PickClosestIntersection(cachedIntersectionList, intersection);
        }

        // Returns true if the given point is inside this solid object.
        // This is a default implementation that counts intersections
        // that enter or exit the solid in a given direction from the point.
        // Derived classes can often implement a more efficient algorithm
        // to override this default algorithm.
        virtual bool Contains(const Vector& point) const;

        // Returns the optical properties (reflection and refraction)
        // at a given point on the surface of this solid.
        // By default, the optical properties are the same everywhere,
        // but a derived class may override this behavior to create
        // patterns of different colors or gloss.
        // It is recommended to keep constant refractive index
        // throughout the solid, or the results may look weird.
        virtual Color SurfaceOptics(
            const Vector& surfacePoint,
            const void *context) const
        {
            return uniformColor;
        }

        // The following three member functions rotate this
        // object counterclockwise around a line parallel
        // to the x, y, or z axis, as seen from the positive
        // axis direction.
        virtual SolidObject& RotateX(double angleInDegrees) = 0;
        virtual SolidObject& RotateY(double angleInDegrees) = 0;
        virtual SolidObject& RotateZ(double angleInDegrees) = 0;

        // Moves the entire solid object by the delta values dx, dy, dz.
        // Derived classes that override this method must chain to it
        // in order to translate the center of rotation.
        virtual SolidObject& Translate(double dx, double dy, double dz)
        {
            center.x += dx;
            center.y += dy;
            center.z += dz;
            return *this;
        }

        // Moves the center of the solid object to
        // the new location (cx, cy, cz).
        SolidObject& Move(double cx, double cy, double cz)
        {
            Translate(cx - center.x, cy - center.y, cz - center.z);
            return *this;
        }

        // Moves the center of the solid object to the
        // location specified by the position vector newCenter.
        SolidObject& Move(const Vector& newCenter)
        {
            Move(newCenter.x, newCenter.y, newCenter.z);
            return *this;
        }

        const Vector& Center() const { return center; }

        void SetFullMatte(const Color& matteColor)
        {
            uniformColor = matteColor;
        }

    protected:
        const Color& GetUniformOptics() const
        {
            return uniformColor;
        }

    private:
        Vector center;  // The point in space about which this object rotates.

        // By default, a solid object has uniform optical properties
        // across its entire surface.  Unless a derived class
        // overrides the virtual member function SurfaceOptics(),
        // the member variable uniformOptics holds these optical
        // properties.
        Color uniformColor;

        // A flag that indicates whether the Contains() method
        // should try to determine whether a point is inside this
        // solid.  If true, containment calculations proceed;
        // if false, Contains() always returns false.
        // Many derived classes will override the Contains() method
        // and therefore make this flag irrelevant.
        const bool isFullyEnclosed;

        // The following members are an optimization to minimize
        // the overhead and fragmentation caused by repeated
        // memory allocations creating and destroying
        // std::vector contents.
        mutable IntersectionList cachedIntersectionList;
        mutable IntersectionList enclosureList;
    };

    //------------------------------------------------------------------------

    // This class encapsulates the notion of a binary operator
    // that operates on two SolidObjects.  Both SolidObjects
    // must support the Contains() method, or an exception
    // will occur during rendering.
    class SolidObject_BinaryOperator: public SolidObject
    {
    public:
        // The parameters '_left' and '_right' must be dynamically
        // allocated using operator new. This class will own
        // responsibility for deleting them when it is itself deleted.
        SolidObject_BinaryOperator(
            const Vector& _center,
            SolidObject* _left,
            SolidObject* _right)
                : SolidObject(_center)
                , left(_left)
                , right(_right)
        {
        }

        virtual ~SolidObject_BinaryOperator()
        {
            delete left;
            left = nullptr;

            delete right;
            right = nullptr;
        }

        // All rotations and translations are applied
        // to the two nested solids in tandem.

        // The following three member functions rotate this
        // object counterclockwise around a line parallel
        // to the x, y, or z axis, as seen from the positive
        // axis direction.
        virtual SolidObject& RotateX(double angleInDegrees);
        virtual SolidObject& RotateY(double angleInDegrees);
        virtual SolidObject& RotateZ(double angleInDegrees);

        virtual SolidObject& Translate(double dx, double dy, double dz);

    protected:
        SolidObject& Left()  const { return *left;  }
        SolidObject& Right() const { return *right; }

        void NestedRotateX(
            SolidObject &nested,
            double angleInDegrees,
            double a,
            double b);

        void NestedRotateY(
            SolidObject &nested,
            double angleInDegrees,
            double a,
            double b);

        void NestedRotateZ(
            SolidObject &nested,
            double angleInDegrees,
            double a,
            double b);

        // The following list is for caching and filtering
        // intersections with the left and right nested solids.
        // It is mutable to allow modification from const methods.
        mutable IntersectionList tempIntersectionList;

    private:
        SolidObject* left;
        SolidObject* right;
    };

    //------------------------------------------------------------------------

    class SetUnion: public SolidObject_BinaryOperator
    {
    public:
        SetUnion(const Vector& _center, SolidObject* _left, SolidObject* _right)
            : SolidObject_BinaryOperator(_center, _left, _right)
        {
            SetTag("SetUnion");
        }

        virtual void AppendAllIntersections(
            const Vector& vantage,
            const Vector& direction,
            IntersectionList& intersectionList) const
        {
            // Find all intersections with the left solid.
            Left().AppendAllIntersections(vantage, direction, intersectionList);

            // Append all intersections with the right solid.
            Right().AppendAllIntersections(vantage, direction, intersectionList);
        }

        virtual bool Contains(const Vector& point) const
        {
            // A point is inside the set union if
            // it is in either of the nested solids.
            return Left().Contains(point) || Right().Contains(point);
        }
    };

    //------------------------------------------------------------------------

    class SetIntersection: public SolidObject_BinaryOperator
    {
    public:
        SetIntersection(
            const Vector& _center,
            SolidObject* _left,
            SolidObject* _right)
                : SolidObject_BinaryOperator(_center, _left, _right)
        {
            SetTag("SetIntersection");
        }

        virtual void AppendAllIntersections(
            const Vector& vantage,
            const Vector& direction,
            IntersectionList& intersectionList) const;

        virtual bool Contains(const Vector& point) const
        {
            // A point is inside the set intersection if
            // it is inside both of the nested solids.
            return Left().Contains(point) && Right().Contains(point);
        }

    private:
        void AppendOverlappingIntersections(
            const Vector& vantage,
            const Vector& direction,
            const SolidObject& aSolid,
            const SolidObject& bSolid,
            IntersectionList& intersectionList) const;

        bool HasOverlappingIntersection(
            const Vector& vantage,
            const Vector& direction,
            const SolidObject& aSolid,
            const SolidObject& bSolid) const;
    };

    //------------------------------------------------------------------------

    // This derived abstract class is specialized for objects (like torus)
    // that are easy to define in terms of a fixed orientation and position
    // in space, but for which generalized rotation makes the algebra
    // annoyingly difficult. Instead, we allow defining the object in terms
    // of a new coordinate system <r,s,t> and translate locations and rays
    // from <x,y,z> camera coordinates into <r,s,t> object coordinates.
    class SolidObject_Reorientable: public SolidObject
    {
    public:
        explicit SolidObject_Reorientable(const Vector& _center = Vector())
            : SolidObject(_center)
            , rDir(1.0, 0.0, 0.0)
            , sDir(0.0, 1.0, 0.0)
            , tDir(0.0, 0.0, 1.0)
            , xDir(1.0, 0.0, 0.0)
            , yDir(0.0, 1.0, 0.0)
            , zDir(0.0, 0.0, 1.0)
        {
        }

        // Fills in 'intersectionList' with a list of all the
        // intersections found starting at the specified
        // vantage point in the specified direction.
        // Any pre-existing content in 'intersectionList'
        // is discarded first.
        // Returns the number of intersections found,
        // which will have the same value as intersectionList.size().
        virtual void AppendAllIntersections(
            const Vector& vantage,
            const Vector& direction,
            IntersectionList& intersectionList) const;

        virtual SolidObject& RotateX(double angleInDegrees);
        virtual SolidObject& RotateY(double angleInDegrees);
        virtual SolidObject& RotateZ(double angleInDegrees);

        virtual bool Contains(const Vector& point) const
        {
            return ObjectSpace_Contains(ObjectPointFromCameraPoint(point));
        }

        virtual Color SurfaceOptics(
            const Vector& surfacePoint,
            const void *context) const
        {
            return ObjectSpace_SurfaceOptics(
                ObjectPointFromCameraPoint(surfacePoint),
                context);
        }

    protected:
        // The following method is called by AppendAllIntersections,
        // but with 'vantage' and 'direction' vectors transformed
        // from <x,y,z> camera space into <r,s,t> object space.
        // Intersection objects are returned in terms of object coordinates,
        // and they are automatically translated back into camera
        // coordinates by the caller.
        virtual void ObjectSpace_AppendAllIntersections(
            const Vector& vantage,
            const Vector& direction,
            IntersectionList& intersectionList) const = 0;

        // Returns true if the specified point in object space
        // is on or inside the solid object.
        // Actually, well-behaved derived classes should provide
        // a tolerance for points slightly outside the object's
        // boundaries and return true then also.
        // This tolerance handles small floating point rounding
        // errors that may cause a point that is supposed to be
        // considered part of the solid to be incorrectly excluded.
        virtual bool ObjectSpace_Contains(const Vector& point) const = 0;

        virtual Color ObjectSpace_SurfaceOptics(
            const Vector& surfacePoint,
            const void *context) const
        {
            return GetUniformOptics();
        }

        Vector ObjectDirFromCameraDir(const Vector& cameraDir) const
        {
            return Vector(
                DotProduct(cameraDir,rDir),
                DotProduct(cameraDir,sDir),
                DotProduct(cameraDir,tDir));
        }

        Vector ObjectPointFromCameraPoint(const Vector &cameraPoint) const
        {
            return ObjectDirFromCameraDir(cameraPoint - Center());
        }

        Vector CameraDirFromObjectDir(const Vector& objectDir) const
        {
            return Vector(
                DotProduct(objectDir,xDir),
                DotProduct(objectDir,yDir),
                DotProduct(objectDir,zDir));
        }

        Vector CameraPointFromObjectPoint(const Vector& objectPoint) const
        {
            return Center() + CameraDirFromObjectDir(objectPoint);
        }

        void UpdateInverseRotation()
        {
            // See the following Wikipedia articles to explain why
            // the inverse of a rotation matrix is just its transpose.
            // http://en.wikipedia.org/wiki/Rotation_matrix
            // http://en.wikipedia.org/wiki/Orthogonal_matrix

            xDir = Vector(rDir.x, sDir.x, tDir.x);
            yDir = Vector(rDir.y, sDir.y, tDir.y);
            zDir = Vector(rDir.z, sDir.z, tDir.z);
        }

    private:
        // The members rDir, sDir, tDir are unit vectors in the direction of
        // the <r,s,t> object axes, each expressed in <x,y,z> camera space.
        // For any point P = <Px,Py,Pz> in camera coordinates, we can
        // determine object-relative coordinates as dot products
        // <(P-C).rDir,(P-C).sDir,(P-C).tDir>,
        // where C = the center of the object as returned by method Center().
        // Another way to look at this is that (rDir, sDir, tDir) taken
        // together are really just a 3*3 rotation matrix.
        Vector  rDir;
        Vector  sDir;
        Vector  tDir;

        // The members xDir, yDir, zDir are unit vectors in the direction
        // of the <x,y,z> camera axes, each expressed in <r,s,t> object space.
        // These are maintained in tandem with rDir, sDir, tDir as various
        // rotations take place.  Taken together, they form an inverse
        // rotation matrix, so (xDir,yDir,zDir) as a 3*3 matrix
        // is calculated as the transpose of the 3*3 matrix (rDir,sDir,tDir).
        // Because an object is never rotated during the rendering of a given
        // frame, it is more efficient to have both matrices pre-calculated.
        Vector  xDir;
        Vector  yDir;
        Vector  zDir;
    };

    //------------------------------------------------------------------------

    // A thin ring is a zero-thickness circular disc with an optional
    // disc-shaped hole in the center.
    class ThinRing: public SolidObject_Reorientable
    {
    public:
        ThinRing(double _innerRadius, double _outerRadius)
            : SolidObject_Reorientable()
            , r1(_innerRadius)
            , r2(_outerRadius)
        {
            SetTag("ThinRing");
        }

    protected:
        virtual void ObjectSpace_AppendAllIntersections(
            const Vector& vantage,
            const Vector& direction,
            IntersectionList& intersectionList) const;

        virtual bool ObjectSpace_Contains(const Vector& point) const
        {
            if (fabs(point.z) <= EPSILON)
            {
                const double magSquared = point.x*point.x + point.y*point.y;
                return
                    (r1*r1 <= EPSILON + magSquared) &&
                    (magSquared <= EPSILON + r2*r2);
            }

            return false;
        }

    private:
        double  r1;     // The radius of the hole at the center of the ring.
        double  r2;     // The outer radius of the ring.

        // A temporary intersection list, cached inside this object
        // to avoid repeated memory allocations.
        // Marked mutable to allow const functions to cache
        // lists whose memory may be reused as needed.
        mutable IntersectionList tempIntersectionList;
    };

    //------------------------------------------------------------------------

    // A thin disc is a zero-thickness disc.
    // It is implemented as a thin ring with a zero-radius hole.
    class ThinDisc: public ThinRing
    {
    public:
        ThinDisc(double _radius)
            : ThinRing(0.0, _radius)
        {
            SetTag("ThinDisc");
        }
    };

    //------------------------------------------------------------------------

    // A sphere-like object, only with different dimensions allowed in
    // the x, y, and z directions.
    class Spheroid: public SolidObject_Reorientable
    {
    public:
        Spheroid(double _a, double _b, double _c)
            : SolidObject_Reorientable()
            , a(_a)
            , b(_b)
            , c(_c)
            , a2(_a * _a)
            , b2(_b * _b)
            , c2(_c * _c)
        {
            SetTag("Spheroid");
        }

    protected:
        virtual void ObjectSpace_AppendAllIntersections(
            const Vector& vantage,
            const Vector& direction,
            IntersectionList& intersectionList) const;

        virtual bool ObjectSpace_Contains(const Vector& point) const
        {
            const double xr = point.x / a;
            const double yr = point.y / b;
            const double zr = point.z / c;
            return xr*xr + yr*yr + zr*zr <= 1.0 + EPSILON;
        }

    private:
        const double  a;      // radius along the x-axis
        const double  b;      // radius along the y-axis
        const double  c;      // radius along the z-axis

        const double  a2;     // a*a, cached for efficiency
        const double  b2;     // b*b, cached for efficiency
        const double  c2;     // c*c, cached for efficiency
    };

    //------------------------------------------------------------------------

    // A sphere that is more efficient than Spheroid with equal dimensions.
    class Sphere: public SolidObject
    {
    public:
        Sphere(const Vector& _center, double _radius)
            : SolidObject(_center)
            , radius(_radius)
        {
            SetTag("Sphere");   // tag for debugging
        }

        virtual void AppendAllIntersections(
            const Vector& vantage,
            const Vector& direction,
            IntersectionList& intersectionList) const;

        virtual bool Contains(const Vector& point) const
        {
            // Add a little bit to the actual radius to be more tolerant
            // of rounding errors that would incorrectly exclude a
            // point that should be inside the sphere.
            const double r = radius + EPSILON;

            // A point is inside the sphere if the square of its distance
            // from the center is within the square of the radius.
            return (point - Center()).MagnitudeSquared() <= (r * r);
        }

        // The nice thing about a sphere is that rotating
        // it has no effect on its appearance!
        virtual SolidObject& RotateX(double angleInDegrees) { return *this; }
        virtual SolidObject& RotateY(double angleInDegrees) { return *this; }
        virtual SolidObject& RotateZ(double angleInDegrees) { return *this; }

    private:
        double  radius;
    };

    //------------------------------------------------------------------------

    // For now, all light sources are single points with an inherent color.
    // Luminosity of the light source can be changed by multiplying
    // color.red, color.green, color.blue all by a constant value.
    struct LightSource: public Taggable
    {
        LightSource(const Vector& _location, const Color& _color, std::string _tag = "")
            : Taggable(_tag)
            , location(_location)
            , color(_color)
        {
        }

        Vector  location;
        Color   color;
    };

    //------------------------------------------------------------------------

    class Aimer        // base class for arbitrary vector aiming logic
    {
    public:
        virtual Vector Aim(const Vector& raw) const = 0;
    };

    //------------------------------------------------------------------------

    // The Scene object renders a collection of SolidObjects and
    // LightSources that illuminate them.
    // SolidObjects are added one by one using the method AddSolidObject.
    // Likewise, LightSources are added using AddLightSource.
    class Scene
    {
    public:
        explicit Scene(const Color& _backgroundColor = Color())
            : backgroundColor(_backgroundColor)
            , aimer(nullptr)
        {
        }

        virtual ~Scene()
        {
            ClearSolidObjectList();
        }

        void SetAimer(Aimer *_aimer)
        {
            aimer = _aimer;
        }

        // Caller must allocate solidObject via operator new.
        // This class will then own the responsibility of deleting it.
        SolidObject& AddSolidObject(SolidObject* solidObject)
        {
            solidObjectList.push_back(solidObject);
            return *solidObject;
        }

        void AddLightSource(const LightSource &lightSource)
        {
            lightSourceList.push_back(lightSource);
        }

        // Renders an image of the current scene, with the camera
        // at <0, 0, 0> and looking into the +z axis, with the +y axis upward.
        // Writes the image to the specified PNG file, which should have a
        // ".png" extension.
        // The resulting image will have pixel dimensions pixelsWide wide
        // by pixelsHigh high.
        // The zoom factor specifies magnification level: use 1.0
        // to start with, and try larger/smaller values to
        // increase/decrease magnification.
        // antiAliasFactor specifies what multiplier to use
        // for oversampling.  Note that this causes run time and memory usage
        // to increase O(N^2), so it is best to use a value between 1
        // (fastest but most "jaggy") to 4 (16 times slower but results
        // in much smoother images).
        void SaveImage(
            const char *outPngFileName,
            size_t pixelsWide,
            size_t pixelsHigh,
            double zoom,
            size_t antiAliasFactor) const;

    private:
        void ClearSolidObjectList();

        int FindClosestIntersection(
            const Vector& vantage,
            const Vector& direction,
            Intersection& intersection) const;

        bool HasClearLineOfSight(
            const Vector& point1,
            const Vector& point2) const;

        Color TraceRay(
            const Vector& vantage,
            const Vector& direction,
            Color rayIntensity,
            int recursionDepth) const;

        Color CalculateLighting(
            const Intersection& intersection,
            const Vector& direction,
            Color rayIntensity,
            int recursionDepth) const;

        Color CalculateMatte(const Intersection& intersection) const;

        void ResolveAmbiguousPixel(ImageBuffer& buffer, size_t i, size_t j) const;

        // Convert a floating point color component value,
        // based on the maximum component value,
        // to a byte RGB value in the range 0x00 to 0xff.
        static unsigned char ConvertPixelValue(
            double colorComponent,
            double maxColorValue)
        {
            int pixelValue =
                static_cast<int> (255.0 * colorComponent / maxColorValue);

            // Clamp to the allowed range of values 0..255.
            if (pixelValue < 0)
            {
                pixelValue = 0;
            }
            else if (pixelValue > 255)
            {
                pixelValue = 255;
            }

            return static_cast<unsigned char>(pixelValue);
        }

        // The color to use for pixels where no solid
        // object intersection was found.
        Color backgroundColor;

        // Define some list types used by member variables below.
        typedef std::vector<SolidObject*> SolidObjectList;
        typedef std::vector<LightSource> LightSourceList;

        // Define types needed to hold a list of pixel coordinates.
        struct PixelCoordinates
        {
            size_t i;
            size_t j;

            PixelCoordinates(size_t _i, size_t _j)
                : i(_i)
                , j(_j)
            {
            }
        };
        typedef std::vector<PixelCoordinates> PixelList;

        // A list of all the solid objects in the scene.
        SolidObjectList solidObjectList;

        // A list of all the point light sources in the scene.
        LightSourceList lightSourceList;

        // Help performance by avoiding constant construction/destruction
        // of intersection lists.
        mutable IntersectionList cachedIntersectionList;
        Aimer *aimer;
    };

    //------------------------------------------------------------------------

    // The information available for any pixel in an ImageBuffer
    struct PixelData
    {
        Color   color;
        bool    isAmbiguous;

        PixelData()
            : color()
            , isAmbiguous(false)
        {
        }
    };

    //------------------------------------------------------------------------
    // Holds an image in memory as it is being rendered.
    // Once calculated, the image in the buffer can be translated
    // into a graphics format like PNG.
    class ImageBuffer
    {
    public:
        ImageBuffer (
            size_t _pixelsWide,
            size_t _pixelsHigh,
            const Color &backgroundColor)
                : pixelsWide(_pixelsWide)
                , pixelsHigh(_pixelsHigh)
                , numPixels(_pixelsWide * _pixelsHigh)
        {
            array = new PixelData[numPixels];
        }

        virtual ~ImageBuffer()
        {
            delete[] array;
            array = nullptr;
            pixelsWide = pixelsHigh = numPixels = 0;
        }

        // Returns a read/write reference to the pixel data for the
        // specified column (i) and row (j).
        // Throws an exception if the coordinates are out of bounds.
        PixelData& Pixel(size_t i, size_t j) const
        {
            if ((i < pixelsWide) && (j < pixelsHigh))
                return array[(j * pixelsWide) + i];

            throw ImagerException("Pixel coordinate(s) out of bounds");
        }

        size_t GetPixelsWide() const
        {
            return pixelsWide;
        }

        size_t GetPixelsHigh() const
        {
            return pixelsHigh;
        }

        // Finds the maximum red, green, or blue value in the image.
        // Used for automatically scaling the image brightness.
        double MaxColorValue() const
        {
            double max = 0.0;
            for (size_t i=0; i < numPixels; ++i)
            {
                array[i].color.Validate();
                if (array[i].color.red > max)
                    max = array[i].color.red;
                if (array[i].color.green > max)
                    max = array[i].color.green;
                if (array[i].color.blue > max)
                    max = array[i].color.blue;
            }
            if (max == 0.0)
            {
                // Safety feature: the image is solid black anyway,
                // so there is no point trying to scale it.
                // If we did, we would end up dividing by zero.
                max = 1.0;
            }
            return max;
        }

    private:
        size_t  pixelsWide;     // the width of the image in pixels (columns).
        size_t  pixelsHigh;     // the height of the image in pixels (rows).
        size_t  numPixels;      // the total number of pixels.
        PixelData*  array;      // flattened array [pixelsWide * pixelsHigh].
    };

    // Output operators (print helpful debug information).
    std::ostream& operator<< (std::ostream&, const Color&);
    std::ostream& operator<< (std::ostream&, const Vector&);
    std::ostream& operator<< (std::ostream&, const Intersection&);
    void Indent(std::ostream&, int depth);
}

#endif // __DDC_IMAGER_H

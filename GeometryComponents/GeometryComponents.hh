#ifndef GEOMETRYCOMPONENTS_HH
#define GEOMETRYCOMPONENTS_HH
/*
 * @brief:  this file gives the class declaration of geometry componnets
 * which will be used in PET simulations and reconstruction.
 * @auther: chengaoyu2013@gmail.com
 * @date:   2017/05/24
*/
/* Note:the unit of length is mm. */

#include<cstddef>
#include<cmath>
#include<limits>
#include<iostream>
#include<algorithm>
#include<memory>
#include<vector>

//#include<math.h>
//#include<stdio.h>


class SphericalPoint3D;

/*
 * Class declaration:
 * this class set the grid size
 * e.g. an image, a Block with nX*nY*nZ size.
*/
class GridSize
{
public:
    GridSize();
    GridSize(size_t X,size_t Y, size_t Z);
    virtual ~GridSize();
    //set the size in 3 coordinates of the grid.
    virtual void setnX(size_t nX);
    virtual void setnY(size_t nY);
    virtual void setnZ(size_t nZ);
    virtual void setGridSize(size_t nX, size_t nY, size_t nZ);
    virtual void setn(size_t num, int index);

    //get the grid size of 3 coordinates.
    virtual size_t getnX() const;
    virtual size_t getnY() const;
    virtual size_t getnZ() const;

    //get nth index
    virtual size_t getn(int index) const;
    virtual bool CheckSameSize(const GridSize& gs) const;

private:
    size_t nx;
    size_t ny;
    size_t nz;
};



/*
 *this class set a point in 3D space.
 *the values (px,py,pz) are set to be float to reduce the computing cost,
//and the type "float" is accurate enough.
*/
class Point3D{
public:
    Point3D();
    Point3D(const Point3D& pt);
    Point3D(float px,float py, float pz);
    //transform a GridSize object and create a Point3D object;
    Point3D(const GridSize &gs);
    virtual ~Point3D();
    // set the three coordinates of the point.
    virtual void setPoint(float px, float py, float pz);
    virtual void setpx(float px);
    virtual void setpy(float py);
    virtual void setpz(float pz);
    // getter
    virtual float getpx() const;
    virtual float getpy() const;
    virtual float getpz() const;
    // get the components by index.(for loops)
    virtual float getp(int index) const;
    virtual void  setp(float value, int index);
    // calculate the distance between two points.
    virtual float calculateDist(const Point3D& pt1) const;

    Point3D operator+(const Point3D& pt1);
    Point3D operator-(const Point3D& pt1);
    Point3D operator*(const Point3D& pt1);
    Point3D operator/(const Point3D& pt);
    Point3D operator/(float ref);
    bool operator==(const Point3D& pt) const;
    //Point3D operator=(const Point3D& pt1);
    //transform a point in rectangular coordinate system to a point in spherical coordinate system.
    virtual SphericalPoint3D PT2SPT() const;
    virtual float getNorm2() const;

private:
    float px;
    float py;
    float pz;
};

/*
 * Class declaration:
 * this class set a point in spherical coordinate system
*/
class SphericalPoint3D{
public:
    SphericalPoint3D();
    SphericalPoint3D(float radius, float theta, float phi);
    virtual ~SphericalPoint3D();
    //set the three coordinate of the point.
    virtual void setSPoint(float radius,float theta, float phi);
    virtual void setRadius(float radius);
    virtual void setTheta(float theta);
    virtual void setPhi(float phi);
    //getter
    virtual float getRadius() const;
    virtual float getTheta() const;
    virtual float getPhi() const;
    // calculate the distance between two points
    virtual float calculateDist(const SphericalPoint3D& spt1) const;
    //const SphericalPoint3D& getSPT();
    //transform a point in spherical coordinate system to a point in rectangular coordinate system.
    virtual Point3D SPT2PT() const;

private:
    float radius;
    float theta;
    float phi;
};


/*
 * Class declaration:
 * set a Ray emitted from a point to another.
*/
typedef Point3D direction;

class Ray{
public:
    Ray();
    Ray(const Point3D& pt1, const Point3D& pt2);
    virtual ~Ray();
    //
    virtual void setstartP(const Point3D& pt);
    virtual void setendP(const Point3D& pt);
    virtual void setTs(float mint, float maxt);
    virtual void setRay(const Point3D& pt1 ,const Point3D pt2);

    // get the two end points.
    virtual const Point3D& getstartP() const;
    virtual const Point3D& getendP() const;
    // get the min_t and max_t
    virtual float getMint() const;
    virtual float getMaxt() const;

    // get the direction of this ray.
    virtual const direction& getDir() const;
    // get the length of this Ray.
    virtual float getLen() const;
    // update the dir and distance when the ends are changed.
    virtual void update();


private:
    Point3D startPt;
    Point3D endPt;

    direction dir;   // the unit vector from startPt to endPt.(direction cosines)
    float     len;   // length between the two points.
    float     min_t; // the 1st initial point(t=0), reset by Intersect... routines.
    float     max_t; // the 2nd initial point(t = len), reset by Intersect... routines.
};



/*
 * Class declaration:
 * this class models a Block in the RingPET.
*/
typedef Point3D BlockSize;
typedef Point3D RegionSize;
//typedef GridSize BlockGrid;
//typedef Point3D RegionSize;

class Block{
public:
    Block();
    virtual ~Block();
    //setter and getter
    virtual void setBlockGrid(const GridSize &blockgrid);
    virtual void setBlockSize(const Point3D &size);
    virtual void setBlockCenter(const Point3D &center);

    virtual const GridSize& getBlockGrid() const;
    virtual const BlockSize& getBlockSize() const;
    virtual const BlockSize& getBlockCenter() const;

    // locate a mesh in the block and return its position(a 3D point)
    virtual Point3D LocatePoint(size_t mi) const;            //when given an index
    virtual Point3D LocatePoint(const GridSize& iMesh) const;   //when given the index of a mesh
    virtual Point3D LocatePoint(const Point3D& pt) const;       //when given the position of a point in space

    virtual GridSize LocateIndex(const Point3D& pt) const;

    //calculate the size of a single mesh in the grid.
    virtual Point3D calculateInterval() const;
    //calculate the center of the first mesh in the block.
    virtual Point3D calculateOffset() const;
    //calculate the number of meshes
    virtual size_t getTotalMeshes() const;

    //get the bound of this block.
    virtual Point3D getTopBound() const;
    virtual Point3D getBottomBound() const;

private:
    BlockSize  bs;     //the region size of block
    GridSize   bg;    //the number of meshes of the block
    Point3D    bc;    //the center postion of block


};

/*
 * Class declaration:
 * this class defines a 3D grid in the space.
*/
typedef GridSize MeshIndex;

class Grid3D:public Block
{
public:
    Grid3D();
    virtual ~Grid3D();
    virtual int initialization();

    //the all the meshes to one value.
    virtual void setAllGridValue(float value);
    virtual void addAllGridValue(float value);
    virtual void minusAllGridValue(float value);
    virtual void multiAllGridValue(float value);
    virtual void divideAllGridValue(float value);

    virtual void setGridValue(const MeshIndex& mi, float value);
    virtual float getGridValue(const MeshIndex& mi) const;

    virtual void setGridValue(size_t mi, float value);
    virtual float getGridValue(size_t mi) const;

    //override to update the size of value vector.
    virtual void setBlockGrid(const GridSize& blockgrid) override;

    // compare the  blockgrid of two Grid3D object.
    virtual bool SizeCheck(const Grid3D& gd);

    // 3D filtering of the grid.
    virtual void conv3(const Grid3D & gd);

    // arithmetic operation of two grids.
    Grid3D  operator+(const Grid3D& gd);
    Grid3D  operator-(const Grid3D& gd);
    Grid3D  operator*(const Grid3D& gd);
    Grid3D  operator/(const Grid3D& gd);

    //calculate the exponent of the mesh values.
    virtual void power(float gd);

	//compute the sqrt of the mesh values
	virtual void sqrt();
    //get the absolution of the mesh values
    virtual void abs();
    //find the minimum and maxmium value of the 'GDValue'
    virtual void findMinMax(float& min, float& max) const;
    //sum the GDvalue in the grid.
    virtual float sumGDvalue() const;

private:
    //the value in each mesh
    std::vector<std::vector<std::vector<float> > > GDValue;
};

/*
 * Class declaration:
 * this class defines an valid event recorded by detector
*/

//the event need to be added and the Ray is used
//typedef Ray Event;

/*
 * Class declaration:
 * this class set a ray cast.
 * reset by SetupRayCast in backprojection.
*/
class RayCast{
public:
    //friend class Backprojection;
    RayCast();
    virtual ~RayCast();
    //
    virtual void SetupRayCast(const Block &imgbox, const Ray& ray);
    virtual int  SetupRayCastComponent(const Block& imgbox, const Ray& ray, int comp);

    //the following members are declared as public for convenience
    //but this is unsafe!!
    //must make them private later
    int	boffs;            // buffer offset for voxel
    float nextT[3];       // distance to next crossing
    float deltaT[3];      // distance between successive crossings
    int   deltaBuf[3];       // increment in the buffer offset for a crossing
    int   inBuf[3];     // the number of crossings which take us out of the imaging volume

};

/*
 * Class declaration:
 * this class define the rotation of a Point along the axis
*/

enum class Axis{
    x = 1,
    y = 2,
    z = 3
};

class Rotation{
public:
    Rotation();
    virtual ~Rotation();
    virtual void setMatrix(float radian, Axis dimension);
    virtual Point3D rotate(const Point3D& pt) const;
    virtual void reset();
    //virtual Rotation operator=(const Rotation& r) const;
private:
    float rotation_matrix[3][3];
    // rotation_matrix
    // [0,0][0,1][0,2]
    // [1,0][1,1][1,2]
    // [2,0][2,1][2,2]
};

#endif // GEOMETRYCOMPONENTS_HH

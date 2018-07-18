/*
 * @brief:
 * @auther: chengaoyu2013@gmail.com
 * @date:   2017/05/24
*/
#include "../GeometryComponents.hh"

/* class--- Block
 * function definitions
*/

Block::Block():bs(Point3D(1,1,1)),bg(GridSize(10,10,10)),bc(Point3D(0,0,0)){
    //default initialization
    /*
     * the blocksize is set to be positive
     */
}

Block::~Block(){
    //do nothing!
}

void Block::setBlockGrid(const GridSize &blockgrid){
    //the grid number of the block should be positive!
    if(blockgrid.getnX()>0&&blockgrid.getnY()>0&&blockgrid.getnZ()>0)
        this->bg = blockgrid;
    else{
        std::cout<<"Invalid input! The block grid should be positive! Use default setting!"<<std::endl;
        this->bg = GridSize(10,10,10);
    }
}

void Block::setBlockSize(const Point3D& size){
    //the size of the block should be positive!
    if(size.getpx()>0&&size.getpy()>0&&size.getpz()>0)
        this->bs = size;
    else{
        std::cout<<"Invalid input! The block size should be positive! Use default setting!"<<std::endl;
        this->bs = Point3D(1,1,1);
    }
}

void Block::setBlockCenter(const Point3D &center ){
    this->bc = center;
}

const GridSize& Block::getBlockGrid() const{
    return this->bg;
}

const Point3D& Block::getBlockSize() const{
    return this->bs;
}

const Point3D& Block::getBlockCenter() const{
    return this->bc;
}

Point3D Block::LocatePoint(size_t mi) const{
    int plane = this->getBlockGrid().getnZ()*this->getBlockGrid().getnY();
    int line = this->getBlockGrid().getnZ();
    size_t ix =(size_t) (mi/plane);
    size_t iy =(size_t)((mi-ix*plane)/line);
    size_t iz =(size_t)((mi%plane)%line);
    GridSize iMesh(ix,iy,iz);
    Point3D pt = LocatePoint(iMesh);
    return  pt;
}

Point3D Block::LocatePoint(const GridSize& iMesh) const{

    //calculate the position of the first mesh(0,0,0)
    Point3D firstmesh = this->calculateOffset();
    //transform the mesh index into float for calculating
    Point3D iMeshfloat(iMesh);
    //check if the index is in the range of the block.
    //because size_t is unsigned, so they are alway positive.
    if(     iMesh.getnX()>this->getBlockGrid().getnX()-1||
            iMesh.getnY()>this->getBlockGrid().getnY()-1||
            iMesh.getnZ()>this->getBlockGrid().getnZ()-1){
        std::cout<<"mesh index out of range! "<<std::endl;
        // if invalid input,return the center position of the block.
        return this->getBottomBound();
    }

    Point3D interval  = this->calculateInterval();
    return Point3D(firstmesh+iMeshfloat*interval);
}

Point3D Block::LocatePoint(const Point3D& pt) const{
    //calculate the position of the first mesh(0,0,0)
    Point3D firstmesh = this->calculateOffset();
    Point3D bbound = this->getBottomBound();
    Point3D meshpt(pt.getpx(),pt.getpy(),pt.getpz());
    Point3D dist = meshpt-bbound;
    Point3D interval  = this->calculateInterval();
    Point3D meshindex = (dist/interval);
    Point3D iMeshfloat(std::floor(meshindex.getpx()),std::floor(meshindex.getpy()),std::floor(meshindex.getpz()));
    //check if the index is in the range of the block.
    if(((int)iMeshfloat.getpx())<0||
            ((int)iMeshfloat.getpy())<0||
            ((int)iMeshfloat.getpz())<0||
            ((int)iMeshfloat.getpx())>((int)(this->getBlockGrid().getnX()-1))||
            ((int)iMeshfloat.getpy())>((int)(this->getBlockGrid().getnY()-1))||
            ((int)iMeshfloat.getpz())>((int)(this->getBlockGrid().getnZ()-1))){
        //std::cout<<"point out of range! "<<std::endl;
        // if invalid input,return the center position of the block.
        return bbound;
    }
    return Point3D(firstmesh+iMeshfloat*interval);
}

GridSize Block::LocateIndex(const Point3D &pt) const{
    //calculate the position of the first mesh(0,0,0)
    Point3D thept(pt);
    Point3D firstmesh = this->calculateOffset();
    Point3D diff = thept - firstmesh;
    Point3D interval =this->calculateInterval();
    Point3D indiff = diff/interval;
    GridSize ptIndex((size_t)std::round(indiff.getpx()),(size_t)std::round(indiff.getpy()),(size_t)std::round(indiff.getpz()));
    return ptIndex;
}

Point3D Block::calculateInterval() const{
    //obtain the block size
    Point3D thegs = this->getBlockGrid();
    //obtain the size of a mesh
    Point3D thebs = this->getBlockSize();
    //calculate the size of a mesh
    Point3D interval = Point3D(thebs/thegs);
    return interval;
}

//calculate the first mesh center in the block
Point3D Block::calculateOffset() const{

    Point3D interval = this->calculateInterval();
    //obtain the center of the block
    Point3D thebc = this->getBlockCenter();
    //obtain the size of a mesh
    Point3D thebs = this->getBlockSize();

    Point3D mesh = Point3D(thebc+interval/2-thebs/2);
    return mesh;
}

Point3D Block::getBottomBound() const{
    Point3D bbound;
    Point3D thebc(this->getBlockCenter());
    Point3D thebs(this->getBlockSize());
    bbound = thebc-(thebs/2);
    return bbound;
}

Point3D Block::getTopBound() const{
    Point3D tbound;
    Point3D thebc(this->getBlockCenter());
    Point3D thebs(this->getBlockSize());
    tbound = thebc+(thebs/2);
    return tbound;
}

size_t Block::getTotalMeshes() const{
    size_t X = this->getBlockGrid().getnX();
    size_t Y = this->getBlockGrid().getnY();
    size_t Z = this->getBlockGrid().getnZ();
    return X*Y*Z;
}

/*
 * @brief:  Define a 3 dimensional grid.
 * @auther: chengaoyu2013@gmail.com
 * @date:   2017/05/24
*/
#include "../GeometryComponents.hh"

GridSize::GridSize()
{
    //default initialization
    this->setGridSize(10,10,10);
}

GridSize::GridSize(size_t X, size_t Y, size_t Z)
{
    setGridSize(X,Y,Z);
}

GridSize::~GridSize(){
    //do nothing
}

void GridSize::setnX(size_t nX){
    this->nx = nX;
}

void GridSize::setnY(size_t nY){
    this->ny = nY;
}

void GridSize::setnZ(size_t nZ){
    this->nz = nZ;
}

void GridSize::setGridSize(size_t nX, size_t nY, size_t nZ){
    this->setnX(nX);
    this->setnY(nY);
    this->setnZ(nZ);
}

void GridSize::setn(size_t num, int index){
    switch (index) {
    case 0:
        this->setnX(num);
        break;
    case 1:
        this->setnY(num);
        break;
    case 2:
        this->setnZ(num);
        break;
    default:
        break;
    }
}

size_t GridSize::getnX() const{
    return nx;
}

size_t GridSize::getnY() const{
    return ny;
}

size_t GridSize::getnZ() const{
    return nz;
}

size_t GridSize::getn(int index) const{
    switch (index) {
    case 0:
        return this->getnX();
    case 1:
        return this->getnY();
    case 2:
        return this->getnZ();
    default:
        std::cout<<"GridSize::getn(): "<<index<<" invalid value! "<<std::endl;
        std::exit(-1);
    }
}

bool GridSize::CheckSameSize(const GridSize &gs) const{
    if (this->getnX()==gs.getnX()&&this->getnY() ==gs.getnY() &&this->getnZ() == gs.getnZ())
        return true;
    else
        return false;
}

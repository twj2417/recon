/*
 * @brief:  definite Grid3D in "Geometry.hh"
 * @auther: chengaoyu2013@gmail.com
 * @date:   2017/05/24
*/
#include "../GeometryComponents.hh"
/* class--- Grid3D
 * function definitions
*/

Grid3D::Grid3D(){
    this->initialization();
}

Grid3D::~Grid3D(){

}
int Grid3D::initialization(){
    GridSize block = this->getBlockGrid();
    size_t nx = block.getnX();
    size_t ny = block.getnY();
    size_t nz = block.getnZ();
    //initialize the grid with zeros.
    GDValue.clear();
    GDValue.resize(nx);
    for(auto& column : GDValue){
        column.resize(ny);
        for(auto& layer : column){
            layer.resize(nz);
        }
    }
    return this->getTotalMeshes();
}

void Grid3D::setGridValue(const MeshIndex &mi, float value){
    /*
     * bound check to be added!
     */
    size_t ix = mi.getnX();
    size_t iy = mi.getnY();
    size_t iz = mi.getnZ();
    GDValue[ix][iy][iz] =value;
}

void Grid3D::setGridValue(size_t mi, float value){
    /*
     * bound check to be added!
     */
    int plane = this->getBlockGrid().getnX()*this->getBlockGrid().getnY();
    int line = this->getBlockGrid().getnX();
    size_t iz =(size_t) (mi/plane);
    size_t iy =(size_t)((mi-iz*plane)/line);
    size_t ix =(size_t)((mi%plane)%line);
    GDValue[ix][iy][iz] =value;
}

void Grid3D::setAllGridValue(float value){
    size_t nX = this->getBlockGrid().getnX();
    size_t nY = this->getBlockGrid().getnY();
    size_t nZ = this->getBlockGrid().getnZ();
    for(size_t iX = 0;iX<nX;iX++){
        for(size_t iY = 0;iY<nY;iY++){
            for(size_t iZ = 0; iZ<nZ;iZ++){
                GridSize loc(iX,iY,iZ);
                this->setGridValue(loc,value);
            }
        }
    }
}

void Grid3D::addAllGridValue(float value){
    size_t nX = this->getBlockGrid().getnX();
    size_t nY = this->getBlockGrid().getnY();
    size_t nZ = this->getBlockGrid().getnZ();
    for(size_t iX = 0;iX<nX;iX++){
        for(size_t iY = 0;iY<nY;iY++){
            for(size_t iZ = 0; iZ<nZ;iZ++){
                GridSize loc(iX,iY,iZ);
                this->setGridValue(loc,this->getGridValue(loc)+value);
            }
        }
    }
}

void Grid3D::minusAllGridValue(float value){
    addAllGridValue(-value);
}

void Grid3D::divideAllGridValue(float value){
    if(value == 0){
        std::cout<<"devideAllGridValue():zero dividor! "<<std::endl;
        std::exit(-1);
    }
    size_t nX = this->getBlockGrid().getnX();
    size_t nY = this->getBlockGrid().getnY();
    size_t nZ = this->getBlockGrid().getnZ();
    for(size_t iX = 0;iX<nX;iX++){
        for(size_t iY = 0;iY<nY;iY++){
            for(size_t iZ = 0; iZ<nZ;iZ++){
                GridSize loc(iX,iY,iZ);
                this->setGridValue(loc,this->getGridValue(loc)/value);
            }
        }
    }
}
void Grid3D::multiAllGridValue(float value){
    size_t nX = this->getBlockGrid().getnX();
    size_t nY = this->getBlockGrid().getnY();
    size_t nZ = this->getBlockGrid().getnZ();
    for(size_t iX = 0;iX<nX;iX++){
        for(size_t iY = 0;iY<nY;iY++){
            for(size_t iZ = 0; iZ<nZ;iZ++){
                GridSize loc(iX,iY,iZ);
                this->setGridValue(loc,this->getGridValue(loc)*value);
            }
        }
    }
}

float Grid3D::getGridValue(const MeshIndex &mi) const{

    /*
     * bound check to be added!
     */
    size_t ix = mi.getnX();
    size_t iy = mi.getnY();
    size_t iz = mi.getnZ();
    return GDValue[ix][iy][iz];
}

float Grid3D::getGridValue(size_t mi) const{

    /*
     * bound check to be added!
     */
    int plane = this->getBlockGrid().getnX()*this->getBlockGrid().getnY();
    int line = this->getBlockGrid().getnX();
    size_t iz =(size_t) (mi/plane);
    size_t iy =(size_t)((mi-iz*plane)/line);
    size_t ix =(size_t)((mi%plane)%line);
    return GDValue[ix][iy][iz];
}

void Grid3D::setBlockGrid(const GridSize &blockgrid){
    Block::setBlockGrid(blockgrid);
    this->initialization();
}

bool Grid3D::SizeCheck(const Grid3D& gd){
    if(this->getBlockGrid().CheckSameSize(gd.getBlockGrid()))
        return true;
    else
        return false;
}

void Grid3D::conv3(const Grid3D &gd){
    size_t f[3];
    GridSize fg = gd.getBlockGrid();
    f[0] = fg.getnX();
    f[1] = fg.getnY();
    f[2] = fg.getnZ();

    size_t nX = this->getBlockGrid().getnX();
    size_t nY = this->getBlockGrid().getnY();
    size_t nZ = this->getBlockGrid().getnZ();
    Grid3D filtered;
    GridSize fs(nX+f[0],nY+f[1],nZ+f[2]);
    filtered.setBlockGrid(fs);
    //filtered.setBlockSize(this->getBlockSize());

    size_t sX = gd.getBlockGrid().getnX()/2;
    size_t sY = gd.getBlockGrid().getnY()/2;
    size_t sZ = gd.getBlockGrid().getnZ()/2;



    //padding the original image cube.
    for(size_t iX = 0;iX<nX;iX++){
        for(size_t iY = 0;iY<nY;iY++){
            for(size_t iZ = 0; iZ<nZ;iZ++){
                GridSize loc1(iX,iY,iZ);
                GridSize loc2(iX+sX,iY+sY,iZ+sZ);
                float value = this->getGridValue(loc1);
                filtered.setGridValue(loc2,value);
            }
        }
    }

    for(size_t iX = 0;iX<nX;iX++){
        for(size_t iY = 0;iY<nY;iY++){
            for(size_t iZ = 0;iZ<nZ;iZ++){
                //Grid3D batch = gd;
                GridSize loc(iX,iY,iZ);
                float value=0;
                for(size_t ibx = 0;ibx < f[0];ibx++){
                    for(size_t iby = 0;iby < f[1];iby++){
                        for(size_t ibz = 0;ibz < f[2];ibz++){
                            GridSize floc(ibx,iby,ibz);
                            GridSize b(iX+ibx,iY+iby,iZ+ibz);
                            value += filtered.getGridValue(b)*gd.getGridValue(floc);
                        }
                    }
                }
                this->setGridValue(loc,value);
            }
        }
    }
    //return filtered;
}

Grid3D Grid3D::operator+(const Grid3D& gd){
    Grid3D addGD = gd;
    if(this->SizeCheck(gd)){// same size.
        size_t nX = this->getBlockGrid().getnX();
        size_t nY = this->getBlockGrid().getnY();
        size_t nZ = this->getBlockGrid().getnZ();
        for(size_t iX = 0;iX<nX;iX++){
            for(size_t iY = 0;iY<nY;iY++){
                for(size_t iZ = 0; iZ<nZ;iZ++){
                    GridSize loc(iX,iY,iZ);
                    float value = this->getGridValue(loc)+gd.getGridValue(loc);
                    addGD.setGridValue(loc,value);
                }
            }
        }
        return addGD;
    }
    else{
        std::cout<<"Grid3D addtion: two grids have different grid size! "<<std::endl;
        std::exit(-1);
    }
}

Grid3D Grid3D::operator-(const Grid3D& gd){
    Grid3D minusGD = gd;
    if(this->SizeCheck(gd)){// same size.
        size_t nX = this->getBlockGrid().getnX();
        size_t nY = this->getBlockGrid().getnY();
        size_t nZ = this->getBlockGrid().getnZ();
        for(size_t iX = 0;iX<nX;iX++){
            for(size_t iY = 0;iY<nY;iY++){
                for(size_t iZ = 0; iZ<nZ;iZ++){
                    GridSize loc(iX,iY,iZ);
                    float value = this->getGridValue(loc)-gd.getGridValue(loc);
                    minusGD.setGridValue(loc,value);
                }
            }
        }
        return minusGD;
    }
    else{
        std::cout<<"Grid3D substraction: two grids have different grid size! "<<std::endl;
        std::exit(-1);
    }
}

Grid3D Grid3D::operator/(const Grid3D& gd){
    Grid3D divideGD = gd;
    if(this->SizeCheck(gd)){// same size.
        size_t nX = this->getBlockGrid().getnX();
        size_t nY = this->getBlockGrid().getnY();
        size_t nZ = this->getBlockGrid().getnZ();
        for(size_t iX = 0;iX<nX;iX++){
            for(size_t iY = 0;iY<nY;iY++){
                for(size_t iZ = 0; iZ<nZ;iZ++){
                    GridSize loc(iX,iY,iZ);
                    if(gd.getGridValue(loc) == 0.0){ // if the dividor is zero, set the result zero.
                        divideGD.setGridValue(loc,0.0);
                    }
                    else{
                        float value = this->getGridValue(loc)/gd.getGridValue(loc);
                        divideGD.setGridValue(loc,value);
                    }
                }
            }
        }
        return divideGD;
    }
    else{
        std::cout<<"Grid3D division: two grids have different grid size! "<<std::endl;
        std::exit(-1);
    }
}

Grid3D Grid3D::operator*(const Grid3D& gd){
    Grid3D multiGD = gd;
    if(this->SizeCheck(gd)){// same size.
        size_t nX = this->getBlockGrid().getnX();
        size_t nY = this->getBlockGrid().getnY();
        size_t nZ = this->getBlockGrid().getnZ();
        for(size_t iX = 0;iX<nX;iX++){
            for(size_t iY = 0;iY<nY;iY++){
                for(size_t iZ = 0; iZ<nZ;iZ++){
                    GridSize loc(iX,iY,iZ);
                    float value = this->getGridValue(loc)*gd.getGridValue(loc);
                    multiGD.setGridValue(loc,value);
                }
            }
        }
        return multiGD;
    }
    else{
        std::cout<<"Grid3D multiplication: two grids have different grid size! "<<std::endl;
        std::exit(-1);
    }
}

//void Grid3D::operator=(const Grid3D& gd){

//}

void Grid3D::power(float val){
    size_t nX = this->getBlockGrid().getnX();
    size_t nY = this->getBlockGrid().getnY();
    size_t nZ = this->getBlockGrid().getnZ();
    for(size_t iX = 0;iX<nX;iX++){
        for(size_t iY = 0;iY<nY;iY++){
            for(size_t iZ = 0; iZ<nZ;iZ++){
                GridSize loc(iX,iY,iZ);
                float value = std::pow(this->getGridValue(loc),val);
                this->setGridValue(loc,value);
            }
        }
    }
}

void Grid3D::sqrt() {
	size_t nX = this->getBlockGrid().getnX();
	size_t nY = this->getBlockGrid().getnY();
	size_t nZ = this->getBlockGrid().getnZ();
	for (size_t iX = 0; iX<nX; iX++) {
		for (size_t iY = 0; iY<nY; iY++) {
			for (size_t iZ = 0; iZ<nZ; iZ++) {
				GridSize loc(iX, iY, iZ);
				float value = this->getGridValue(loc);
				if (value > 0)
					value = std::sqrt(value);
				else
					value = 0;
				this->setGridValue(loc, value);
			}
		}
	}
}

void Grid3D::abs(){
    size_t nX = this->getBlockGrid().getnX();
    size_t nY = this->getBlockGrid().getnY();
    size_t nZ = this->getBlockGrid().getnZ();
    for(size_t iX = 0;iX<nX;iX++){
        for(size_t iY = 0;iY<nY;iY++){
            for(size_t iZ = 0; iZ<nZ;iZ++){
                GridSize loc(iX,iY,iZ);
                float value = this->getGridValue(loc);
                value = value>=0? value:(-value);
                this->setGridValue(loc,value);
            }
        }
    }
}

float Grid3D::sumGDvalue() const{
    float sumV = 0;
    size_t nX = this->getBlockGrid().getnX();
    size_t nY = this->getBlockGrid().getnY();
    size_t nZ = this->getBlockGrid().getnZ();
    for(size_t iX = 0;iX<nX;iX++){
        for(size_t iY = 0;iY<nY;iY++){
            for(size_t iZ = 0; iZ<nZ;iZ++){
                GridSize loc(iX,iY,iZ);
                sumV += this->getGridValue(loc);
            }
        }
    }
    return sumV;
}

void Grid3D::findMinMax(float &min, float &max) const{
    GridSize p0(0,0,0);
    float fmin,fmax;
    fmin = fmax = this->getGridValue(p0);
    size_t nX = this->getBlockGrid().getnX();
    size_t nY = this->getBlockGrid().getnY();
    size_t nZ = this->getBlockGrid().getnZ();
    for(size_t iX = 0;iX<nX;iX++){
        for(size_t iY = 0;iY<nY;iY++){
            for(size_t iZ = 0; iZ<nZ;iZ++){
                GridSize loc(iX,iY,iZ);
                float val = this->getGridValue(loc);
                fmin = fmin<val? fmin:val;
                fmax = fmax>val? fmax:val;
            }
        }
    }
    min = fmin;
    max = fmax;
}

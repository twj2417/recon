/*
 * @brief:
 * @auther: chengaoyu2013@gmail.com
 * @date:   2017/07/10
*/

#include"ImgFilters.hh"

filter::filter(){

}

filter::~filter(){

}

Grid3D filter::Kaiser(const Block bk, float alpha, int order, float radius){
    Grid3D ft;
    ft.setBlockGrid(bk.getBlockGrid());
    size_t appro = 100; // the larger , the higher accuracy.
    //calculate the modified bessel function.
    double b[appro];
    double interval = radius/appro;
    double Im = boost::math::cyl_bessel_i(order,alpha);
    for(size_t i =0;i<appro;i++){
        double temp = std::sqrt(1-std::pow((interval*(i+0.5))/radius,2));
        double temp1 = temp*alpha;
        double I = boost::math::cyl_bessel_i(order,temp1);
        b[i] =pow(temp,order)*I/Im;
    }

    //size_t len = bk.getBlockGrid().getnX(); // the filter is assumed to be of the same size in all dimensions.
    //generate the 3D kaiser window.
    size_t nX = ft.getBlockGrid().getnX();
    size_t nY = ft.getBlockGrid().getnY();
    size_t nZ = ft.getBlockGrid().getnZ();
    double start = 0;
    if(nX%2 == 0){
        start =0.5;
    }
    if(nX%2 == 1){
        start = 0;
    }
    float gridInterval = 2*radius/nX;
    for(size_t iX = 0;iX<nX;iX++){
        for(size_t iY = 0;iY<nY;iY++){
            for(size_t iZ = 0; iZ<nZ;iZ++){
                GridSize loc(iX,iY,iZ);
                double dis = std::sqrt(std::pow(gridInterval*(std::abs((int)iX-(int)nX/2+start)),2)+std::pow(gridInterval*(std::abs((int)iY-(int)nY/2+start)),2)+std::pow(gridInterval*(std::abs((int)iZ-(int)nZ/2+start)),2));
                int ind = std::round(dis/interval);
                float value = 0;
                if(ind<appro)
                    value = b[ind];
                ft.setGridValue(loc,value);
            }
        }
    }


    return ft;
}

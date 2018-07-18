/*
 * @brief:  this file defines the class TOF
 * @auther: chengaoyu2013@gmail.com
 * @date:   2017/05/18
*/

#include "TOF.hh"
/* class---TOF
 * function definitions
*/

TOF::TOF(){
    FWHM = 30; // set the default full width at half maximum to 30mm (time resolution is 200ps).
    flag = false; // TOF switch is off in default.
    limit  = 9; // set the limit to 3*sigma.
    binsize = 3.66; // the size of a bin.

    Sigma2 = FWHM*FWHM/(8.0*M_LN2);
    std::cout<<"TOF information: "<<std::endl;
    std::cout<<"FWHM: "<<FWHM <<"mm, bin size:"<<binsize<<"mm, limit range:"<<limit<<"mm." <<std::endl;

}

TOF::TOF(float fwhm, float bs,float range){
    flag  = true;
    if( fwhm<0 || binsize<0 || range < 0){
        std::cout<<"TOF initialization: negative value! "<<std::endl;
        std::cout<<"Use default setting! "<<std::endl;
        FWHM =30;
        binsize = 3.66;
        limit = 9;
    }
    else{
        FWHM = fwhm;
        binsize = bs;
        limit = range;
    }
    Sigma2 = FWHM*FWHM/(8.0*M_LN2);
    std::cout<<"TOF information: "<<std::endl;
    std::cout<<"FWHM: "<<FWHM <<"mm, bin size:"<<binsize<<"mm, limit range:"<<limit<<"mm." <<std::endl;
}

TOF::~TOF(){

}

void TOF::setBinSize(float bs){
    binsize = bs;
}

void TOF::setFWHM(float fwhm){
    FWHM = fwhm;
}

void TOF::setLimit(float range){
    limit = range;
}

float TOF::getlimit2() const{
    return limit*limit;
}

float TOF::getSigma2() const{
    return Sigma2;
}


float TOF::getBinSize() const{
    return binsize;
}

float TOF::getFWHM() const{
    return FWHM;
}

void TOF::on(){
    flag = true;
}
void TOF::off(){
    flag = false;
}

bool TOF::isopen() const{
    return flag;
}

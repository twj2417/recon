/*
 * @brief:  the definition of class SphericalPoint3D
 *          declared in "Geometry.hh"
 * @auther: chengaoyu2013@gmail.com
 * @date:   2017/05/24
*/

#include "../GeometryComponents.hh"


/* class---SphericalPoint3D
 * function definitions
*/
SphericalPoint3D::SphericalPoint3D(){
    setSPoint(0.0,0.0,0.0);
}

SphericalPoint3D::SphericalPoint3D(float radius, float theta, float phi){
    setSPoint(radius,theta,phi);
}

SphericalPoint3D::~SphericalPoint3D(){

}

void SphericalPoint3D::setSPoint(float radius, float theta, float phi){

    this->setRadius(radius);
    this->setTheta(theta);
    this->setPhi(phi);
}

void SphericalPoint3D::setRadius(float radius){
    if(radius<0){
        std::cout<<"radius can not be negative!"<<std::endl;
    }
    else{
        this->radius = radius;
    }
}

void SphericalPoint3D::setTheta(float theta){
    //the polar angle ranges from 0~M_PI.
    if(theta>M_PI || theta<0){
        std::cout<<"polar angle out of range!"<<std::endl;
    }
    else{
        this->theta = theta;
    }
}

void SphericalPoint3D::setPhi(float phi){
    // the azimuth angle ranges from 0~2M_PI
    if(phi>2*M_PI || phi<0){
        std::cout<<"azimuth angle out of range!"<<std::endl;
    }
    else {
        this->phi =phi;
    }
}

float SphericalPoint3D::getRadius() const{
    return radius;
}

float SphericalPoint3D::getTheta() const{
    return theta;
}

float SphericalPoint3D::getPhi() const{
    return phi;
}

//calaulate the distance between two points in spherical coordinate system.
float SphericalPoint3D::calculateDist(const SphericalPoint3D &spt1) const{

    Point3D pt = this->SPT2PT();
    Point3D pt1 = spt1.SPT2PT();
    return pt.calculateDist(pt1);
}

//transform a point in spherical coordinate system to a point in rectangular coordinate system.
Point3D SphericalPoint3D::SPT2PT() const{
    Point3D pt(0,0,0);
    pt.setpx(this->getRadius()*sin(this->getTheta())*cos(this->getPhi()));
    pt.setpy(this->getRadius()*sin(this->getTheta())*sin(this->getPhi()));
    pt.setpz(this->getRadius()*cos(this->getTheta()));
    return pt;
}

//transform a point in rectangle coordinate system to a point in spherical coordinate system.
SphericalPoint3D Point3D::PT2SPT() const{
    SphericalPoint3D spt(0,0,0);
    spt.setRadius(sqrt(pow(this->getpx(),2)+pow(this->getpy(),2)+pow(this->getpz(),2)));
    if(spt.getRadius()<0){// the radius shoule not be negative.
        std::cout<<"transfrom error from rectangle coordinate system point to spherical coordinate system !"<<std::endl;
    }
    else{
        spt.setTheta(acos(this->getpz()/spt.getRadius()));
    }
    float phi = atan2(this->getpy(),this->getpx());
    //spt.setPhi(phi);
    //add M_PI to shift the range of phi to 0~2M_PI.
    spt.setPhi(phi= (phi>=0)? phi:phi+2*M_PI);
    return spt;
}
/**********************************************/
/****** class function definition above *******/

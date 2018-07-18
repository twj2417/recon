/*
 * @brief:  definition of class Ray in "Geometry.hh"
 * @auther: chengaoyu2013@gmail.com
 * @date:   2017/05/24
*/

#include"../GeometryComponents.hh"

//class Point;

/* class---Ray
 * function definitions
*/
Ray::Ray():startPt(Point3D(0,0,0)),endPt(Point3D(1,1,1)){

    this->update();
}

Ray::Ray(const Point3D& pt1, const Point3D& pt2){
    this->startPt = pt1;
    this->endPt = pt2;
    this->update();
}

Ray::~Ray(){

}

const direction &Ray::getDir() const{
    return this->dir;
}

const Point3D& Ray::getstartP() const{
    return startPt;
}

const Point3D &Ray::getendP() const{
    return endPt;
}

float Ray::getMint() const{
    return min_t;
}

float Ray::getMaxt() const{
    return max_t;
}

float Ray::getLen() const{
    return len;
}

void Ray::setstartP(const Point3D &pt){
    this->startPt = pt;
    //update();
}

void Ray::setendP(const Point3D &pt){
    this->endPt = pt;
    //update();
}

void Ray::setRay(const Point3D &pt1, const Point3D pt2){
    startPt = pt1;
    endPt = pt2;
    update();
}

void Ray::setTs(float mint, float maxt){
    min_t = mint;
    max_t = maxt;
}

void Ray::update(){
    len = endPt.calculateDist(startPt);
    Point3D temp = endPt-startPt;
    dir = temp/len;
}

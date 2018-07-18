/*
 * @brief:  the definition of class Point3D
 *          declared in "Geometry.hh"
 * @auther: chengaoyu2013@gmail.com
 * @date:   2017/05/24
*/

#include "../GeometryComponents.hh"

/* class---Point3D
 * function definitions
*/
//constructor
Point3D::Point3D(){
    //if no parameter, set the point in origin.
    setPoint(0.0,0.0,0.0);
}
Point3D::~Point3D(){

}

Point3D::Point3D(const Point3D& pt){
    setPoint(pt.getpx(),pt.getpy(),pt.getpz());
}

Point3D::Point3D(float px,float py, float pz){
    setPoint(px, py, pz);
}

Point3D::Point3D(const GridSize& gs){
    setPoint((float)gs.getnX(), (float)gs.getnY(),(float)gs.getnZ());
}

//setter
void Point3D::setPoint(float px, float py, float pz){
    this->setpx(px);
    this->setpy(py);
    this->setpz(pz);
}

void Point3D::setpx(float px){
    this->px = px;
}

void Point3D::setpy(float py){
    this->py = py;
}

void Point3D::setpz(float pz){
    this->pz = pz;
}

//getter
float Point3D::getpx() const{
    return px;
}
float Point3D::getpy() const{
    return py;
}
float Point3D::getpz() const{
    return pz;
}

float Point3D::getp(int index) const{
    switch (index) {
    case 0:
        return px;
    case 1:
        return py;
    case 2:
        return pz;
    default:
        std::cout<<"Point3D::getp(): "<<index<<" is invalid index! "<<std::endl;
        std::exit(-1);
    }
}

void Point3D::setp(float value, int index){
    switch (index) {
    case 0:
        this->setpx(value);
        break;
    case 1:
        this->setpy(value);
        break;
    case 2:
        this->setpz(value);
        break;
    default:
        break;
    }
}

Point3D Point3D::operator+(const Point3D& pt1){
    return Point3D(this->getpx()+pt1.getpx(),this->getpy()+pt1.getpy(),this->getpz()+pt1.getpz());
}

Point3D Point3D::operator-(const Point3D& pt1){
    return Point3D(float(px-pt1.getpx()),float(py-pt1.py),float(pz-pt1.pz));
}

Point3D Point3D::operator*(const Point3D& pt1){
    return Point3D(this->getpx()*pt1.getpx(),this->getpy()*pt1.getpy(),this->getpz()*pt1.getpz());
}

Point3D Point3D::operator/(const Point3D& pt){
    Point3D thePt(0,0,0);
    //when the divisor has 0 component, set the result to be zero.
    if(pt.getpx()==0)
        thePt.setpx(0);
    else
        thePt.setpx(this->getpx()/pt.getpx());
    if(pt.getpy()==0)
        thePt.setpy(0);
    else
        thePt.setpy(this->getpy()/pt.getpy());
    if(pt.getpz()==0)
        thePt.setpz(0);
    else
        thePt.setpz(this->getpz()/pt.getpz());
    //input point is taken as the divisor, and 'this point' is the dividend.
    return thePt;
}

Point3D Point3D::operator/ (float ref){
    if(ref ==0){// divisor can not be 0.
        std::cout<<"invalid value! point3D divided by zero!"<<std::endl;
        Point3D pt(0,0,0);
        return pt;
    }
    Point3D pt(this->getpx()/ref,this->getpy()/ref,this->getpz()/ref);
    return pt;
}

bool Point3D::operator==(const Point3D& pt1)const{
    if(this->getpx()==pt1.getpx()&&this->getpy()==pt1.getpy()&&this->getpz()==pt1.getpz()){
        return true;
    }
    else
        return false;
}


float Point3D::calculateDist(const Point3D &pt1) const{
    Point3D dir(float(px-pt1.getpx()),float(py-pt1.getpy()),float(pz-pt1.getpz()));
    float dist =std::sqrt(std::pow(dir.getpx(),2)+std::pow(dir.getpy(),2)+std::pow(dir.getpz(),2));
    return dist;
}

// calculate the  2 norm
float Point3D::getNorm2() const{
    return std::sqrt(pow(this->getpx(),2)+pow(this->getpy(),2)+pow(this->getpz(),2));
}

/**********************************************/
/****** class function definition above *******/

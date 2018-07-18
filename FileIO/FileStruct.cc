/*
 * @brief:  this file give the definition of methods in "FileStruct.hh"
 * @auther: chengaoyu2013@gmail.com
 * @date:   2017/04/12
*/

//#include "../GeometryComponents/GeometryComponents.hh"
#include "FileStruct.hh"
/* class---LMElement
 * function definitions
*/


Event::Event():time_TOF(-1),ray(Ray()){
}

Event::Event(const Point3D& pt1, const Point3D& pt2):
    time_TOF(-1),ray(pt1,pt2){
}

Event::Event(const Point3D& pt1, const Point3D& pt2, const float weight):
    time_TOF(-1),ray(pt1,pt2),weight(1){

}
Event::~Event(){

}



LMElement::LMElement(){
    this->addr = 0;
    this->dist = 0;
}

LMElement::LMElement(size_t theaddr,float thedist){
    this->addr = theaddr;
    this->dist = thedist;
}

LMElement::~LMElement(){

}

void LMElement::setAddr(size_t theaddr){
    this->addr = theaddr;
}

void LMElement::setDist(float thedist){
    this->dist = thedist;
}

size_t LMElement::getAddr() const{
    return this->addr;
}

float LMElement::getDist() const{
    return this->dist;
}

/* class---LMEvents
 * function definitions
*/

LMEvent::LMEvent(){
    //default is void.
    this->reset();
}

LMEvent::~LMEvent(){

}

void LMEvent::setLMElement(size_t index, const LMElement &lmele){
    if(this->array.size()>index){// valid index
        LMElement ele = lmele;
        array.assign(index,ele);
    }
    else{
        std::cout<<"setLMElement(): invalid index,out of the array! "<<std::endl;
        std::exit(-1);
    }
}

void LMEvent::addLMElement(const LMElement &lmele){
    LMElement ele = lmele;
    array.push_back(ele);
    this->nvoxels += 1;
}

void LMEvent::setNevents(size_t ne){
    nevents = ne;
}

void LMEvent::setNvoxels(size_t nv){
    nvoxels = nv;
}

const LMElement& LMEvent::getLMElement(size_t index) const{
    if(this->array.size()>index){ // valid index
        return this->array.at(index);
    }
    else{
        std::cout<<"getLMElement(): invalid index,out of the array! "<<std::endl;
        std::exit(-1);
    }
}

size_t LMEvent::getNvoxels() const {
    return nvoxels;
}

size_t LMEvent::getNevents() const{
    return nevents;
}

// reset the LMEvent object.
void LMEvent::reset(){
    nvoxels = 0;
    nevents = 1;
    array.clear();
}





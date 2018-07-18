#ifndef FILESTRUCT_HH
#define FILESTRUCT_HH

#include<cstddef>
#include<vector>
#include<iostream>
#include"../GeometryComponents/GeometryComponents.hh"

#define cspeed  3e11; //the speed of light is 3e11 mm/second.
/*
 * Class declaration:
 * LMElement contains a list mode factor (system matrix element)
*/

//typedef Ray Event;
class Event{
public:
    Event();
    Event(const Point3D& pt1, const Point3D& pt2);
    Event(const Point3D& pt1, const Point3D& pt2,const float weight);
    virtual ~Event();
    virtual Ray& getRay(){return ray;}
    virtual void settimeTOF(float dt){time_TOF = dt;}
    virtual void setweight(float wt){weight = wt;}
	virtual float gettimeTOF() const { return time_TOF; }
	virtual float getDistTOF() const {return 0.5*time_TOF*cspeed;}
    virtual float getweight() const {return weight;}
private:
	float     time_TOF;
    Ray       ray;
    float     weight;

};



class LMElement{
public:
    LMElement();
    LMElement(size_t theadds, float thedist);
    virtual ~LMElement();
    virtual void setAddr(size_t theaddr);
    virtual void setDist(float thedist);

    virtual size_t getAddr() const;
    virtual float getDist() const;
private:
    size_t addr;    // location of voxel in line integral
    float dist;         // length of chord through voxel

};


/*
 * Class declaration:
 * LMEVENT contains all of the listmode factors (system matrix elements) for an event.
 * It is row of the system matrix, and each factor is stored as a LMELEMENT.
*/
class LMEvent{
public:
    LMEvent();
    virtual ~LMEvent();

    //setter and getter
    virtual void setNvoxels(size_t nv);
    virtual void setNevents(size_t nv);

    virtual size_t getNvoxels() const;
    virtual size_t getNevents() const;

    virtual const LMElement &getLMElement(size_t index) const;
    virtual void setLMElement(size_t index,const LMElement& lmele);

    // add a LMElement in the end of the array.
    virtual void addLMElement(const LMElement& lmele);

    // reset the LMEvent object.
    virtual void reset();
private:
    size_t nvoxels; // number of in the chord line integral.
    size_t nevents; // number of events for this chord
    std::vector<LMElement> array; // locations of voxels and lengths in line integral
};

#endif //FILESTRUCT_HH

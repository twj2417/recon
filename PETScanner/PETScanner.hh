#ifndef PETSCANNER_HH
#define PETSCANNER_HH

/*
 * @brief:
 * this file declares the different kind of PET scanners.
 * A PET scanner is represented by all of the crystals it owns.
 * These crystals are stored in a crystal List, which is the main data member of a PET scanner.
 * A general PETScanner class is taken as a base class an it can represent any scanner by the crystal List.
 *
 * @auther: chengaoyu2013@gmail.com
 * @date:   2017/03/20
*/

#include <vector>
//#include<cstddef>
#include <fstream>
#include <string>
#include <ctime>

#include "../GeometryComponents/GeometryComponents.hh"
#include "../ImageProcessing/ImgFilters.hh"
#include "../FileIO/FileControl.hh"
#include "TOF.hh"

//define the image volume
typedef Grid3D PETImg;

//use a point to represent a crystal in the PET scanner.
typedef Point3D cryst;
//use a 3D vector to represent the size of a crystal.
typedef Point3D crystalSize;

enum class RayTracingMethod
{
    Huesmann = 1,
    Siddon = 2

    // if other method need to be added, add them below.
};

/*
 * Class declaration:
 * a base class for PET scanner which store the position of the crystals in a list.
 * The center(a 3D point) of a crystal is stored in the cryst list.
*/
class PETScanner
{
  public:
    PETScanner();
    virtual ~PETScanner();
    //
    virtual void initializeScanner();
    virtual void initializeScanner(const std::string &crystListFile) final;

    //add a crystal in the end of the crystal list.
    virtual void addCrystal(const cryst &aCryst);

    //set the a crystal in a given index.
    virtual void setCrystal(size_t inCryst, const cryst &aCryst);

    // set the system name
    virtual void setName(const std::string petname) { PETSystem = petname; }
    // get the system name of this PET Scanner.
    virtual const std::string getName() const { return PETSystem; }

    //set and get the detector attenuation length
    virtual void setDetAtnLen(float dal) { detatn_length = dal; }
    virtual float getDetAtnLen() const { return detatn_length; }

    //set and get the LMEflag which decides if generate a LMEvent in the backprojection process.
    virtual void setLMEFlag(bool flag) { LMEflag = flag; }
    virtual bool getLMEFlag() const { return LMEflag; }

    //set and get the ABFflag which decides if use the after backprojection filtering method.
    virtual void setABFFlag(bool flag) { ABFflag = flag; }
    virtual bool getABFFlag() const { return ABFflag; }

    //set and get the image to be reconstructed.
    virtual void setImage(const PETImg &i) { img = i; }
    virtual PETImg &getImage() { return img; }

    virtual void setLMEvent(const LMEvent &lme) { RTlmevt = lme; }
    virtual const LMEvent &getLMEvent() const { return RTlmevt; }

    // set the total number of crystals in the scanner.
    // virtual void setNumOfCrystals(size_t nCrysts);

    // get the total number of crystals in the scanner.
    virtual size_t getNumOfCrystals() const;

    // clear the Crystal list.
    virtual void clearCrystalList();

    // fill the crystal list.
    virtual void fillCrystalList();

    // calculate the crystal index of a given point.
    virtual size_t locateCrystal(const Point3D &pt) const;

    //given an index, return the corresponding crystal.
    virtual cryst GetCrystal(size_t inCryst) const;

    /// the algorithms to reconstruct an image
    ///
    virtual void setTOFinfo(const TOF &tof) { tofinfo = tof; }
    virtual const TOF &getTOFinfo() const { return tofinfo; }
    // set and get the ray tracing method.
    virtual void setRayTracingMethod(RayTracingMethod rt) { RTMtd = rt; }
    virtual RayTracingMethod getRayTracingMethod() const { return RTMtd; }
    // backproject an event(a LOR line) into the image.

    //calculate the efficiency map of the PET scanner
    virtual void EfficiencyMap(const std::string &vefilename);


    virtual void RayTracing(Event &evt);

    //calculate the voxel length by Huesmann's method.
    virtual void Huesmann(Event &evt);

    //calculate the voxel length by siddon's algorithm.
    virtual void Siddon(Event &evt);

    //intersect a ray with the PET and return the intersect length.
    virtual float IntersectPET(Event &evt);

    //backproject a ray into the image
    virtual void BackProjRay(float t_TOF, float wt, RayCast &raycast);

    //static methods are declared below
    //calculate the intersection points of a ray with a standard ellipsoid(or ellipse)
    static bool IntersectEllipsoid(RegionSize &ShapeSize, Ray &ray, int ndim);

    //calculate the intersection points of a ray with a box.
    static bool IntersectBBox(Block &box, Ray &ray, int boundindex);

    //calculate the min_t and max_t of a ray with two paralell planes
    static bool ClipRay(float toplim, float botlim, int index, Ray &ray);

    virtual void IterRec(const std::string& evtfilename, const std::string &vefilename, const std::string& recfilename,size_t start_index, size_t nIterations);
    virtual void PreMLEM(const std::string& transevtfilename, PETImg& VoxelEff,PETImg& PriorWeight, PETImg& PriorImg, float penalty);
    virtual void TransEvents(const std::string& evtfilename, const std::string& transevtfilename) = 0;
    virtual void MLEM(const std::string &transevtfilename, PETImg &VoxelEff, PETImg& PriorWeight, PETImg& PriorImg,float comfac,size_t start_index, size_t nIterations, const std::string &RecImgFile);
	virtual void MAP(const std::string &transevtfilename, PETImg &VoxelEff, PETImg& PriorWeight, PETImg& PriorImg, float comfac, size_t start_index, size_t nIterations, const std::string &RecImgFile);
	
	
	virtual void ABF(PETImg& img);
    //the EM method embeds the LMF computing peocess.
    virtual float EM(PETImg &ImgNow, PETImg &ImgNext,float comfac, const std::string &evtfilename);
    // LMF method compute the list-mode event of an input evt.
    virtual LMEvent LMF(Event &evt);
    //the EM_Projection forward and back projecte one list mode event.
    virtual float EM_Projection(PETImg& ImgNow, PETImg& ImgNext, LMEvent& lmevt, float weight );


  private:
    //number of crystals in the whole detector
    //size_t nCrystals;
    //name of this PET system.
    std::string PETSystem;
    //crystal list of this PET system.
    std::vector<cryst> crystals;
    PETImg img;
    //the attenuation length of the crystals.
    float detatn_length;
    RayTracingMethod RTMtd;
    //
    TOF tofinfo;

    bool LMEflag = false; // the flag to generate LMEvents in a backprojection process.(generate when the flag is true, default false)
    LMEvent RTlmevt;
    //
    bool ABFflag = false; // the flag to use the after backprojection filter.
};

#endif // PETSCANNER_HH

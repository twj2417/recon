#ifndef SPHEREPET_HH
#define SPHEREPET_HH

#include"PETScanner.hh"
/*
 * @brief:
 * However, this general model will cost much time when given a point and find which crystal it belongs to.
 * Thus, the specific scanner (SpherePET, RingPET) are devoloped to sort the crystals in the list according
 * to their spatial position. It markablely reduce the time to locate a point to the crystals.
 * This file declares the model of spherical PET scanner.
 * @auther: chengaoyu2013@gmail.com
 * @date:   2017/03/22
 *
 *
 * add the RadialRefineFactor for DOI detector at 2017/05/12
*/



/*
 * Class declaration:
 * SpherePET models the discrete scheme of a spherical PET scanner.
 * The sphere is truncated spherical shell which has approximately 3*pi of solid angle (polar angle(0~5/6*pi)).
 * The shell is dicretized into many ribbons along the polar angle first.
 * Then each ribbon is equally divided into pieces with same shape along the azimuth angle.
 *
 * For DOI detector, the RadialRefineFactor is utilized, however the crystal list will not increase, because the
 * the radial value of the crystal is calculated on the fly.
*/
class SpherePET:public PETScanner{
public:
    SpherePET();
    SpherePET(float innerR, float outerR, size_t prf, size_t arf, size_t rrf, float tpv1, float tpv2);
    virtual ~SpherePET();
    virtual void setInnerR(float r1);
    virtual void setOuterR(float r2);
    // refinement factors
    virtual void setPolarRefineFactor(size_t prf);
    virtual void setAzimuthRefineFactor(size_t arf);
    virtual void setRadialRefineFactor(size_t rrf);
    virtual void setNumOfRibbons(size_t nR);
    // set the value of the Ribbon edges
    virtual void setPolarValue(std::vector<float> pvList);
    virtual void setPolarValue(size_t inPv, float pv);



    //getters
    virtual float getInnerR() const{return innerR;}
    virtual float getOuterR() const{return outerR;}
    virtual size_t getPolarRefineFactor() const{ return polarRefineFactor;}
    virtual size_t getAzimuthRefineFactor() const{ return azimuthRefineFactor;}
    virtual size_t getRadialRefineFactor() const{return radialRefineFactor;}
    virtual size_t getNumOfRibbons() const{ return nRibbons;}

    virtual const std::vector<float> &getPolarValue() const;
    virtual float getPolarValue(size_t inPv) const;

    //get the truncated polar angle
    virtual float getTPV1() const {return truncatedPV1;}
    virtual float getTPV2() const {return truncatedPV2;}
    //get the truncated polar angle
    virtual size_t getTL1() const {return truncatedLevel1;}
    virtual size_t getTL2() const {return truncatedLevel2;}

    virtual void initializeScanner() override;
    //virtual void initializeScanner(float innerR, float outerR, size_t polarRefineFactor, size_t azimuthRefineFactor);

    //calculate the number of ribbons of the detector
    virtual void calculateRibbons();
    // set the truncated polar angle
    virtual void setTruncatedPV(float tpv1, float tpv2);
    // set the truncated Level
    virtual void setTruncatedLevel(size_t tl1, size_t tl2);

    // fill the Crystals (std::vector<cryst>)
    virtual void fillCrystalList() override;

    virtual size_t locateCrystal(const Point3D& pt) const override;

    //
    virtual size_t getNumOfCrystals() const override;
    //
    virtual cryst GetCrystal(size_t inCryst) const override;

    virtual void EfficiencyMap(const std::string &vefilename) override;
    //virtual void LMFactors(const std::string& evtfilename, const std::string &lmffilename) override;

    virtual void TransEvents(const std::string& evtfilename, const std::string& transevtfilename);

private:
    float innerR; //circumscribed sphere of inner radius.
    float outerR;
    //the refinement factor in polar and azimuth angle
    size_t polarRefineFactor;
    size_t azimuthRefineFactor;
    //the refinement factor in radial direction.(DOI in sphere)
    size_t radialRefineFactor;
    //number of the ribbons to discretize a simi-sphere.
    size_t nRibbons;
    //the truncated polar angles and detectors fill in this range.(0 =< PV1 <= PI/2, PI/2 <= PV2 <= PI)
    float truncatedPV1;
    float truncatedPV2;
    size_t truncatedLevel1;
    size_t truncatedLevel2;
    //the arccos value bottom edge of the Ribbons.
    std::vector<float> polarValue;
};
#endif // SPHEREPET_HH

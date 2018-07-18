#ifndef RINGPET_HH
#define RINGPET_HH
/*
 * @brief:
 * However, this general model will cost much time when given a point and find which crystal it belongs to.
 * Thus, the specific scanner (SpherePET, RingPET) are devoloped to sort the crystals in the list according
 * to their spatial position. It markablely reduce the time to locate a point to the crystals.
 * This file declares the model of ringy PET scanner.
 * @auther: chengaoyu2013@gmail.com
 * @date:   2017/04/10
*/


#include"PETScanner.hh"


/*
 * Class declaration:
 * RingPET models the discrete scheme of a ringy PET scanner.
 * The whole scanner is made of many rings and each ring is made of amounts of blocks.
 * The cystals are arranged in rings and blocks.
 * If the number of rings are even, the
*/
class RingPET:public PETScanner{
public:
    RingPET();
    virtual ~RingPET();
    RingPET(float inR,float outR, float ht, size_t nring, size_t nbok, const Block& bok);
    //radius setter and getter.
    virtual void setInnerR(float r1){ innerR = r1;}
    virtual void setOuterR(float r2){ outerR = r2;}
    virtual float getInnerR() const{return innerR;}
    virtual float getOuterR() const{return outerR;}
	virtual float getHeight() const { return Height; }

	virtual size_t getNumBlocks() const{ return nBlocks; }
	virtual size_t getNumOfRings() const { return nRings; }
    // the blockorigin setter and getter.
    virtual void setBlockOrigin(const Block& bok){blockorigin = bok;}
    virtual const Block& getBlockOrigin() const{return blockorigin;}

    // ring interval setter and getter
    virtual void setRingInterval(float ri){ ringInterval = ri;}
    virtual float getRingInterval() const{ return ringInterval;}



    // check if the bound of the scanner is valid.
    // give an error message when the blocks are overlapped or out of the cylinder.
    // and then use the default setting
    virtual bool ValidBound();

    //virtual size_t getNumOfCrystals() override;
    // default set of a RingPET.
    virtual void defaultRingPET();
    virtual void initializeScanner() override;

    //fill the crystal list.
    virtual void fillCrystalList() override;
    virtual size_t locateCrystal(const Point3D& pt) const override;
    virtual size_t getNumOfCrystals() const override;
    virtual bool findCrystal( size_t iBlock, size_t iCryst, cryst &pt) const;

	virtual bool findCrystal( cryst &pt)const;

    virtual void EfficiencyMap(const std::string &vefilename) override;
    //virtual void LMFactors(const std::string& evtfilename, const std::string &lmffilename) override;

    virtual void TransEvents(const std::string& evtfilename, const std::string& transevtfilename);
private:
    float innerR;
    float outerR;
    float Height;
    // the number of blocks in a ring.
    size_t nBlocks;

    // the number of rings in the scanner.
    size_t nRings;

    // the block origin(size and grid)
    // all of the blocks can be generated from this origin by rotate and shift.
    Block blockorigin;
    //the interval between two rings
    float ringInterval;
    //crystalSize crystSize;



};
#endif // RINGPET_HH

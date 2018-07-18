#ifndef DODEPET_HH
#define DODEPET_HH
/*
 * @brief:
 * However, this general model will cost much time when given a point and find which crystal it belongs to.
 * Thus, the specific scanner (SpherePET, DodePET) are devoloped to sort the crystals in the list according
 * to their spatial position. It markablely reduce the time to locate a point to the crystals.
 * This file declares the model of ringy PET scanner.
 * @auther: chengaoyu2013@gmail.com
 * @date:   2017/04/10
*/


#include"PETScanner.hh"


/*
 * Class declaration:
 * DodePET models the discrete scheme of a dodecahedral PET scanner.
 * The whole scanner is made of 12 pentagen blocks and the block orgin is the bottom block.
 * The cystals are arranged in blocks.
*/
class DodePET:public PETScanner{
public:
    DodePET();
    virtual ~DodePET();
    DodePET(float inR, float tn, size_t nbok, const Block& bok);

    //there are 9 rotation matrix in the list.
    virtual void setupParas();
    virtual const Rotation& getRotationMatrix(size_t index) const;

    //radius setter and getter.
    virtual void setInnerR(float r){ innerR = r;}
    virtual float getInnerR() const{ return innerR;}

    virtual void setNBlocks(size_t n){nBlocks= n<=12? n:11;}
    virtual size_t getNBlocks() const{ return nBlocks;}

    virtual void setThickness(float t){ thickness = t;}
    virtual float getThickness() const{return thickness;}
    // the blockorigin setter and getter.
    virtual const Block& getBlockOrigin() const {return blockorigin;}


    // check if the bound of the scanner is valid.
    // give an error message when the blocks are overlapped or out of the cylinder.
    // and then use the default setting
    //virtual bool ValidBound();

    // default set of a DodePET.
    virtual void defaultDodePET();
    virtual void initializeScanner() override;
    //fill the crystal list.
    virtual void fillCrystalList() override; //only a block origin of the detector is saved.
    virtual bool findCrystal(Point3D& pt, const size_t &blockIndex) const;
    virtual bool findCrystal(Point3D& pt, const size_t &blockIndex,size_t index) const;

    virtual size_t getNumOfCrystals() const override;

    virtual void EfficiencyMap(const std::string &vefilename) override;
    //virtual void LMFactors(const std::string& evtfilename, const std::string &transevtfilename) override;

    virtual double cross(const Point3D& p0, const Point3D& p1,  const Point3D& p2) const;
    virtual bool  IsInsidePenta(const Point3D& pt) const;


    virtual void TransEvents(const std::string& evtfilename, const std::string& transevtfilename);


private:
    float innerR;    // inradius of DodePET.
    float thickness; // thickness of the detector;
    size_t nBlocks;  // number of blocks in the DodePET, no more than 12.

    // the block origin(size and grid)
    // all of the blocks can be generated from this origin by rotate and shift.
    Block blockorigin;

    std::vector<Rotation> RotationList;
    std::vector<Point3D> PentaBound;


};

#endif // DODEPET_HH

/*
 * @brief:
 * @auther: chengaoyu2013@gmail.com
 * @date:   2017/03/22
*/

#include "SpherePET.hh"

/* class---SpherePET
 * function definitions
*/
SpherePET::SpherePET()
{

    // default initialization
    innerR = 150.0;
    outerR = 170.0;
    //no refinement
    polarRefineFactor = 1;
    azimuthRefineFactor = 1;
    radialRefineFactor = 1;
    //default no truncated
    truncatedPV1 = 0;
    truncatedPV2 = M_PI;
    //SpherePET::fillList();
    this->initializeScanner();
}

SpherePET::SpherePET(float innerR, float outerR, size_t prf, size_t arf, size_t rrf, float tpv1, float tpv2)
{
    this->innerR = innerR;
    this->outerR = outerR;
    this->polarRefineFactor = prf;
    this->azimuthRefineFactor = arf;
    this->radialRefineFactor = rrf;
    this->setTruncatedPV(tpv1, tpv2);
    this->initializeScanner();
}

SpherePET::~SpherePET()
{
}

void SpherePET::initializeScanner()
{
    //PETScanner::initializeScanner();
    this->clearCrystalList();
    //should not be changed!
    this->setName("SpherePET");
    this->clearCrystalList();
    this->calculateRibbons();
    this->fillCrystalList();
}

void SpherePET::calculateRibbons()
{
    //the number of Ribbons is the multiple of innerR and polar refinement factor.
    //the crystals are divided into peices with 2M_PI area.
    //the polarRefineFactor represent the number of peices in polar that a 2M_PI unit will be divided into.
    //calculate a semi-sphere first.
    size_t semi_nRibbons = polarRefineFactor * (size_t)innerR;
    polarValue.clear();
    float pv = 0;
    for (size_t i = 0; i < semi_nRibbons; i++)
    {

        // a common factor of all Ribbons.
        float temp = 1 / innerR / innerR / polarRefineFactor;
        if (i == 0)
        { //for the first ribbon, it is determined by the unit area of crystal
            pv = (1 - temp);
        }
        else
        {
            size_t iRibbon = i / polarRefineFactor + 1;
            pv = polarValue.at(i - 1) - temp * (2 * iRibbon - 1);
        }
        polarValue.push_back(pv);
    }
    // assign the polar value on the other semi-sphere
    //
    /*for (size_t i = simi_nRibbons; i>0; i--){
        // the values are sysmetric with the first semi-sphere, but negative
        pv = -polarValue.at(i-1);
        polarValue.push_back(pv);
    }*/

    for (auto iter = std::begin(polarValue); iter != std::end(polarValue); ++iter)
    {
        //std::cout<<acos(*iter)<<std::endl;

        //transform to angle (0~M_PI/2)
        *iter = acos(*iter);
    }
    this->setNumOfRibbons(polarValue.size());
}

//calculate the positions of the crystals and store them into the List.
void SpherePET::fillCrystalList()
{
    size_t nRibbons = this->getNumOfRibbons();

    // make a list which contains a semi-sphere of crystals .
    // clear the list first

    //record the polar value index.
    size_t inTPV1 = 0;
    size_t inTPV2 = 0;

    //find the ribbon which the truncated edges belong to.
    //this truncated is done in the situation without refinement.
    for (size_t i = 0; i < nRibbons; i++)
    {
        if (polarValue.at(i) >= truncatedPV1)
        { // the detector start below the top edge.
            inTPV1 = i;
            break;
        }
        else
            continue;
    }
    //find the 2nd truncated edge.
    for (size_t i = 0; i < nRibbons; i++)
    {
        if (polarValue.at(i) >= (M_PI - truncatedPV2))
        { // the detector start below the top edge.
            inTPV2 = i;
            break;
        }
        else
            continue;
    }
    // set the truncated Level
    this->setTruncatedLevel(inTPV1, inTPV2);
    // fill in the top part
    for (size_t i = inTPV1; i < nRibbons; i++)
    { //each ribbon
        size_t iRibbon = i / polarRefineFactor + 1;
        for (size_t j = 0; j < (2 * iRibbon - 1) * azimuthRefineFactor; j++)
        { //each refined ribbon
            SphericalPoint3D spt;
            // the radius of the crystal is set to be the middle value of two edges.
            spt.setRadius((innerR + outerR) / 2);
            if (i == 0) //use the average of polar angles, 1st angle should be 0.
                spt.setTheta(0);
            else // use the average of two polar angles.
                spt.setTheta((polarValue.at(i) + polarValue.at(i - 1)) / 2);
            // set the azimuth angle
            spt.setPhi(2 * M_PI / ((2 * iRibbon - 1) * azimuthRefineFactor) * ((float)j + 0.5));

            //transform to a point in rectangular coordinate system.
            Point3D pt = spt.SPT2PT();
            this->addCrystal(pt);
        }
    }
    // fill in the bottom part
    for (size_t i = nRibbons - 1; i >= inTPV2; i--)
    { //each ribbon
        size_t iRibbon = i / polarRefineFactor + 1;
        for (size_t j = 0; j < (2 * iRibbon - 1) * azimuthRefineFactor; j++)
        { //each refined ribbon
            SphericalPoint3D spt;
            // the radius of the crystal is set to be the middle value of two edges.
            spt.setRadius((innerR + outerR) / 2);
            if (i == 0)
            { //use the average of polar angles, 1st angle should be 0.
                spt.setTheta(M_PI);
            }
            else // use the average of two polar angles.
                spt.setTheta(M_PI - (polarValue.at(i - 1) + polarValue.at(i - 1)) / 2);

            // set the azimuth angle
            spt.setPhi(2 * M_PI / ((2 * iRibbon - 1) * azimuthRefineFactor) * ((float)j + 0.5));
            //transform to a point in rectangular coordinate system.
            Point3D pt = spt.SPT2PT();
            this->addCrystal(pt);
        }
        if (i == 0) //avoid i<0, this is a bug!
            break;
    }
    //this->setNumOfCrystals(1);
}

void SpherePET::setInnerR(float r1)
{
    if (r1 >= 0)
        this->innerR = r1;
    else
        std::cout << "invalid value of inner radius for spherical PET!" << std::endl;
}

void SpherePET::setOuterR(float r2)
{
    if (r2 >= 0)
        this->outerR = r2;
    else
        std::cout << "invalid value of outer radius for spherical PET!" << std::endl;
}

void SpherePET::setPolarRefineFactor(size_t prf)
{
    if (prf >= 1)
        this->polarRefineFactor = prf;
    else
        std::cout << "invalid value of polar refinement factor for spherical PET!" << std::endl;
}

void SpherePET::setAzimuthRefineFactor(size_t arf)
{
    if (arf >= 1)
        this->polarRefineFactor = arf;
    else
        std::cout << "invalid value of azimuth refinement factor for spherical PET!" << std::endl;
}

void SpherePET::setRadialRefineFactor(size_t rrf)
{
    if (rrf >= 1)
        this->radialRefineFactor = rrf;
    else
        std::cout << "invalid value of radial refinement factor(DOI) for spherical PET! " << std::endl;
}

void SpherePET::setTruncatedPV(float tpv1, float tpv2)
{ // tpv1 ranges from 0~M_PI/2, tpv2 ranges from M_PI/2~M_PI.
    // tpv1 and tpv2 range from 0~M_PI.
    if (tpv1 < 0 || tpv1 >= M_PI / 2 || tpv2 <= M_PI / 2 || tpv2 > M_PI)
    {
        std::cout << "invalid value! truncated polar angle out of range! " << std::endl;
        std::cout << "Truncated polar angle is set default: PV1 = 0; PV2 = M_PI. " << std::endl;
        //set the range to be 0~M_PI.
        truncatedPV1 = 0;
        truncatedPV2 = M_PI;
    }
    else
    {
        truncatedPV1 = tpv1;
        truncatedPV2 = tpv2;
    }
}

void SpherePET::setTruncatedLevel(size_t tl1, size_t tl2)
{
    // tpv1 and tpv2 range from 0~M_PI.
    if ((int)tl1 < 0 || (int)tl2 < 0)
    {
        std::cout << "invalid value! truncated Level out of range! " << std::endl;
        truncatedLevel1 = 0;
        truncatedLevel2 = 0;
    }
    else
    {
        truncatedLevel1 = tl1;
        truncatedLevel2 = tl2;
    }
}

void SpherePET::setNumOfRibbons(size_t nR)
{
    if (nR >= 1)
        this->nRibbons = nR;
    else
        std::cout << "invalid value of the  for spherical PET!" << std::endl;
}

size_t SpherePET::locateCrystal(const Point3D &pt) const
{
    SphericalPoint3D spt = pt.PT2SPT();
    // set the index of cryst invalid to indicate an error.
    size_t inCryst = this->getNumOfCrystals();
    // indicate a ribbon this point belongs to
    size_t inRibbon = 0;

    int tl1 = this->getTL1();
    //int tl2 = this->getTL2();
    // get the crystal refinement factors.
    size_t prf = this->getPolarRefineFactor();
    size_t arf = this->getAzimuthRefineFactor();
    size_t rrf = this->getRadialRefineFactor();

    float inR = this->getInnerR();
    float outR = this->getOuterR();
    float sptR = spt.getRadius();

    float sptTheta = spt.getTheta();
    float sptPV1 = this->getTPV1();
    float sptPV2 = this->getTPV2();
    //judge if the point is in the detector.
    if (sptR < inR || sptR > outR || sptTheta < sptPV1 || sptTheta > sptPV2)
    {
        //return a num that is larger than the index of crystals to show a mistake!
        //std::cout<<"the point is not in crystal! "<<std::endl;
        return inCryst;
    }
    else
    {
        //calculate the index assume there is no truncation
        /*
         * 0.theta > M_PI/2? or theta < M_PI/2?
         * 1.find the polar angle level
         * 2.find the azimuth angle
         * 3.minus the number of crystals in truncated cap
         * 4.give the index of the crystal where the point is.
         */

        //calcualte the truncated crystal number of a sphere.
        size_t truncatedNum = (size_t)(tl1 / prf) * (size_t)(tl1 / prf) * arf * prf;
        if (sptTheta <= *(polarValue.end() - 1))
        { //point on the top semi-sphere
            for (size_t i = 0; i < polarValue.size(); i++)
            {
                //find which ribbon the point is located
                if (polarValue.at(i) >= sptTheta)
                {
                    inRibbon = i;
                    break;
                }
                else
                {
                    continue;
                }
            }
            size_t inRibbonLarge = (size_t)((float)inRibbon / prf);                         //the ribbon without refinement.
            inCryst = inRibbonLarge * inRibbonLarge * prf * arf - truncatedNum;             // minus the num of truncated num in large ribbon.
            inCryst += (inRibbon % prf) * (2 * inRibbonLarge + 1) * arf;                    //add the small ribbion
            inCryst += (size_t)(spt.getPhi() / (2 * M_PI / (2 * inRibbonLarge + 1) / arf)); //add the
        }
        else if ((M_PI - sptTheta) <= *(polarValue.end() - 1))
        { //point on the bottom semi sphere
            for (size_t i = 0; i < polarValue.size(); i++)
            {
                //std::cout<<"pv"<<i<<":"<<M_PI-sptTheta<<std::endl;
                if (polarValue.at(i) > (M_PI - sptTheta))
                {
                    inRibbon = i;
                    break;
                }
                else
                {
                    continue;
                }
            }
            size_t inRibbonLarge = (size_t)((float)inRibbon / prf);                                                         //the ribbon without refinement.
            inCryst = 2 * inR * inR * prf * arf - truncatedNum;                                                             //
            inCryst -= (inRibbonLarge * inRibbonLarge) * prf * arf;                                                         // minus the num of truncated num in large ribbon.
            inCryst -= (inRibbon % prf) * (2 * inRibbonLarge + 1) * arf;                                                    //add the small ribbion
            inCryst -= (2 * inRibbonLarge + 1) * arf - (size_t)(spt.getPhi() / (2 * M_PI / (2 * inRibbonLarge + 1) / arf)); //add the
        }
        else //ignore the points on zero-Z-plane(Theta = M_PI/2)
            return inCryst;

        // layer indicated the DOI information.
        int layer;
        layer = (sptR - inR) / ((outR - inR) / rrf);
        inCryst = inCryst + (this->getNumOfCrystals() / rrf) * layer;
    }

    return inCryst;
}

void SpherePET::setPolarValue(size_t inPv, float pv)
{
    if (pv >= 0 && pv <= M_PI)
        polarValue[inPv] = pv;
    else
        std::cout << "invalid polar value for the Spherical PET !" << std::endl;
}

void SpherePET::setPolarValue(std::vector<float> pvList)
{
    //clear the old list and copy the input vector.
    polarValue.clear();
    polarValue = pvList;
}

float SpherePET::getPolarValue(size_t inPv) const
{
    if (inPv <= polarValue.size())
        return polarValue[inPv];
    else
        std::cout << "the index of polar value is out of range!" << std::endl;
}

const std::vector<float> &SpherePET::getPolarValue() const
{
    return polarValue;
}

size_t SpherePET::getNumOfCrystals() const
{
    //the crytal list does not store the DOI information, however, calculate the total
    //number of crystals by multiple the radialRefineFactor.
    return PETScanner::getNumOfCrystals() * this->getRadialRefineFactor();
}

cryst SpherePET::GetCrystal(size_t inCryst) const
{
    // get the crystal list size.
    int singlelayer = PETScanner::getNumOfCrystals();
    int crystIndex = inCryst % singlelayer;
    cryst csttemp = PETScanner::GetCrystal(crystIndex);

    SphericalPoint3D scst = csttemp.PT2SPT();

    // judge the DOI layer
    int layer = inCryst / singlelayer;
    //calcualte the radial position center of the crsytal.
    float DOIPosition = (this->getOuterR() - this->getInnerR()) / this->getRadialRefineFactor() * (layer + 0.5) + this->getInnerR();
    scst.setRadius(DOIPosition);
    Point3D cst = scst.SPT2PT();
    return cst;
}
void SpherePET::EfficiencyMap(const std::string &vefilename)
{
    std::clock_t start, end;
    start = clock();
    double diff;
    int nCrysts = this->getNumOfCrystals();
    double totalLOR = (double)nCrysts * (nCrysts - 1) / 2;

    //loop all of the possible LORs to calculate the voxel efficiency map.
    for (int i = 0; i < nCrysts; i++)
    {
        for (int j = i + 1; j < nCrysts; j++)
        {
            //get the two crystals of a LOR
            Point3D pt1(this->GetCrystal(i));
            Point3D pt2(this->GetCrystal(j));
            Event ray(pt1, pt2);
            // backproject the ray into the voxel efficiency map.
            this->RayTracing(ray);
        }
        if (i % 1000 == 1)
        {
            end = clock();
            diff = (double)(end - start) / CLOCKS_PER_SEC;
            std::cout << "cpu time: " << diff / 60 << " minutes,";
            double timetemp = (double)(nCrysts - i) * (nCrysts - i - 1) / 2;
            std::cout << "the time remains:" << diff / 60 / (totalLOR - timetemp) * timetemp << "miuntes." << std::endl;
        }
    }
    FileControl fc1;
    fc1.SaveImg(vefilename, this->getImage());
}

void SpherePET::TransEvents(const std::string &evtfilename, const std::string &transevtfilename)
{
    FileControl fc;
    size_t nChords = fc.ReadEvtHdr(evtfilename);

    //read the event file.
    std::ifstream evtinfile(evtfilename, std::ios_base::in);

    // write the lmfactor file.
    std::ofstream outfile(transevtfilename, std::ios_base::out);

    size_t nEvents = 0;

    for (int i = 0; i < nChords; i++)
    {
        //get the two crystals of a LOR
        Point3D pt1;
        Point3D pt2;
        float time_TOF;
        fc.ReadEvt(evtinfile, pt1, pt2, time_TOF, this->getTOFinfo().isopen());
        size_t ptind1 = (this->locateCrystal(pt1));
        size_t ptind2 = (this->locateCrystal(pt2));
        if (ptind1 >= this->getNumOfCrystals() || ptind2 >= this->getNumOfCrystals())
        {
            //std::cout<<"invalid position!"<<std::endl;
            continue;
        }
        Point3D cryst1(this->GetCrystal(ptind1));
        Point3D cryst2(this->GetCrystal(ptind2));
        Event evt(cryst1, cryst2);
        evt.settimeTOF(time_TOF);
        //backproject the ray into the voxel efficiency map.
        //RayTracing(evt);
        //save the LMevt
        nEvents += fc.SaveEvt(outfile, evt);
    }
    fc.SaveEvtHdr(transevtfilename, nEvents);
}

/**********************************************/
/****** class function definition above *******/

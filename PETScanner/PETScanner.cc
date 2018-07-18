/*
 * @brief: the definitions of the header file "PETScanner.hh"
 * @auther: chengaoyu2013@gmail.com
 * @date:   2017/03/20
*/
#include "PETScanner.hh"

/* class---PETScanner
 * function definitions
*/
PETScanner::PETScanner() : detatn_length(0)
{
    //the name should not be changed!
    this->setName("PETscanner");
    //default settings.
    initializeScanner();
}
PETScanner::~PETScanner()
{
}

//
void PETScanner::initializeScanner()
{
    //if no parameter, set the vector with only one 'zero' element
    cryst acryst(0.0, 0.0, 0.0);
    crystals.clear();
    crystals.assign({acryst});
}

void PETScanner::initializeScanner(const std::string &crystListFile)
{
    std::ifstream inFile(crystListFile, std::ios_base::in);

    //clear the list before initialization.
    crystals.clear();

    /*
     * read the crystal center postion into the crstal list
     * to be added!!!
     *
     */
}

cryst PETScanner::GetCrystal(size_t inCryst) const
{
    //judge if the index is valid.
    if (crystals.size() <= inCryst || int(inCryst) < 0)
    {
        //if the list length is small than the index or the index is negative.
        std::cout << "GetCrystal(): the index " << inCryst << " of crystal is out of range!!!" << std::endl;
        //std::exit(-1);
        return Point3D (0,0,0);
    }
    else
    {
        cryst cst = crystals[inCryst];
        return cst;
    }
}

void PETScanner::addCrystal(const cryst &aCryst)
{
    cryst pt = aCryst;
    this->crystals.push_back(pt);
}

void PETScanner::setCrystal(size_t inCryst, const cryst &aCryst)
{

    //judge if the index is valid.
    if (crystals.size() <= inCryst || int(inCryst) < 0)
    {
        //if the list length is small than the index or the index is negative.
        std::cout << "the index of crystal is out of range!!!" << std::endl;
        std::exit(-1);
    }
    else
    {
        crystals[inCryst] = aCryst;
    }
}

//void PETScanner::setNumOfCrystals(size_t nCrysts){
//    if(nCrysts>=0) //number of crystals can not be negative.
//        this->nCrystals = nCrysts;
//    else
//        std::cout<<"invalid value for number of crystals!"<<std::endl;
//}

size_t PETScanner::getNumOfCrystals() const
{
    return this->crystals.size();
}

size_t PETScanner::locateCrystal(const Point3D &pt) const
{
    // set the minDistance to be maximum value of type float.
    float minDistance = std::numeric_limits<float>::max();
    size_t inCryst = -1;
    for (auto iter = begin(crystals); iter != end(crystals); ++iter)
    {
        //cryst acryst =
        /*
         * find the minimum distance between crystal and the point
         * then break;
         * to be added!!!
         *
         */
    }
    return inCryst;
}

void PETScanner::fillCrystalList()
{
    //to be added.
}

void PETScanner::clearCrystalList()
{
    crystals.clear();
}

void PETScanner::EfficiencyMap(const std::string &vefilename)
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

// void PETScanner::LMFactors(const std::string &evtfilename, const std::string &lmffilename)
// {
//     // enter the list-mode events saving mode.
//     this->setLMEFlag(true);
//     std::clock_t start, end;
//     start = clock();
//     double diff;
//     FileControl fc;
//     size_t nChords = fc.ReadEvtHdr(evtfilename);

//     //read the event file.
//     std::ifstream evtinfile(evtfilename, std::ios_base::in);

//     // write the lmfactor file.
//     std::ofstream outfile(lmffilename, std::ios_base::out);

//     size_t nLMEvents = 0;

//     for (int i = 0; i < nChords; i++)
//     {
//         //get the two crystals of a LOR
//         Point3D pt1;
//         Point3D pt2;
//         float dist_TOF;
//         fc.ReadEvt(evtinfile, pt1, pt2, dist_TOF, this->getTOFinfo().isopen());
//         size_t ptind1 = (this->locateCrystal(pt1));
//         size_t ptind2 = (this->locateCrystal(pt2));
//         if (ptind1 >= this->getNumOfCrystals() || ptind2 >= this->getNumOfCrystals())
//         {
//             //std::cout<<"invalid position!"<<std::endl;
//             continue;
//         }
//         Point3D cryst1(this->GetCrystal(ptind1));
//         Point3D cryst2(this->GetCrystal(ptind2));
//         Event evt(cryst1, cryst2);
//         evt.setDistTOF(dist_TOF);
//         //backproject the ray into the voxel efficiency map.
//         RayTracing(evt);
//         //save the LMevt
//         nLMEvents += fc.SaveLMevt(outfile, this->getLMEvent());
//         if (i % 10000 == 1)
//         {
//             end = clock();
//             diff = (double)(end - start) / CLOCKS_PER_SEC;
//             std::cout << "cpu time: " << diff / 60 << " miuntes,";
//             std::cout << "the time remains:" << diff / i * (nChords - i) / 60 << "miuntes." << std::endl;
//         }
//     }
//     fc.SaveLMevtHdr(lmffilename, nLMEvents);
//     //exit the list mode events saving mode.
//     this->setLMEFlag(false);
// }

void PETScanner::RayTracing(Event &evt)
{
    if (this->LMEflag)
    {
        RTlmevt.reset();
        //this->img.setAllGridValue(0);
    }
    if (this->RTMtd == RayTracingMethod::Huesmann)
    {
        Huesmann(evt);
    }
    else if (this->RTMtd == RayTracingMethod::Siddon)
    {
        Siddon(evt);
    }
    else
    {
        std::cout << "PETScanner::RayTracing(): invalid RayTracing method! " << std::endl;
        std::exit(-1);
    }
}

void PETScanner::Huesmann(Event &evt)
{
    // get the LOR length
    float LORLength = evt.getRay().getLen();

    // the intersect distance in the detectors.
    float IntersectDist = IntersectPET(evt);
    //float natal = IntersectDist/petinfo.getDetAtnLen();
    float natal = 0;

    //upcasting the Grid3D to Block.
    Block &box = this->getImage();
    if (IntersectBBox(box, evt.getRay(), 1))
    {
        Ray &ray = evt.getRay();
        //the ray intersect with the box.(bounded)
        Point3D bbs = box.getBlockSize();   //get the image size
        Point3D bbc = box.getBlockCenter(); //get the image center
        Point3D stp = ray.getstartP();      //get the start point of the ray.

        // get the bottom bound of the image field.
        Point3D p0 = stp - (bbc - (bbs / 2));
        // translate the origin to the corner of the image volume
        ray.setstartP(p0);
        RayCast raycast;
        raycast.SetupRayCast(box, ray);

        //the weight of this LOR line.
        float wt = std::exp(-natal) / std::pow(LORLength, 2);

        float t_TOF = 0.5 * LORLength - evt.getDistTOF() - evt.getRay().getMint();
        BackProjRay(t_TOF, wt, raycast);
    }
    else
    {
        //do nothing!
    }
}
void PETScanner::Siddon(Event &evt)
{
    //to be finished
}

/****************************************************************************
 * IntersectPET()
 * The following function calculates the combined distance of a ray within
 * the PET system.  Endpoints are assumed to be within detector crystals.
 * Return value: length of the ray traveling in the crystals (between endpoints)
****************************************************************************/
float PETScanner::IntersectPET(Event &evt)
{
    // to be added
}

/****************************************************************************
 * IntersectEllipsoid()
 * Calculate the intersection points of a ray with a standard ellipsoid (or ellipse).
 * The ellipsoid is centered on the origin alligned with the coordinate system.
 * The lengths of the axes are given by size.
 * If the ray intersects, min_t and max_t contain the values along the ray,
 * and 1 is returned.  (Otherwise 0 is returned.)
 * the equation of the ellipsoid is:
 *   4r[0]^2/size[0]^2 + 4r[1]^2/size[1]^2 + 4r[2]^2/size[2]^2 = 1
 * and the equation of the line is:
 *   r[i] = p0[i] + t*pdir[i]
 * For ndim=2, we get the correct result for an ellipse.
 * If the ray intersects, min_t and max_t contain the values along the ray,
 * and true is returned.  (Otherwise false is returned.)
 ****************************************************************************/
bool PETScanner::IntersectEllipsoid(RegionSize &ShapeSize, Ray &ray, int ndim)
{
    // a,b,c will represent the coefficients of equation of two degree
    float a, b, c, delta, fac;
    a = b = 0;
    c = -1.0;
    // the axial length of the elliposoid.
    float ElliposoidSize[3];
    ElliposoidSize[0] = ShapeSize.getpx();
    ElliposoidSize[1] = ShapeSize.getpy();
    ElliposoidSize[2] = ShapeSize.getpz();

    // the direction cosine of the ray.
    float dir[3];
    dir[0] = ray.getDir().getpx();
    dir[1] = ray.getDir().getpy();
    dir[2] = ray.getDir().getpz();

    // the start point of the ray.
    float p0[3];
    p0[0] = ray.getstartP().getpx();
    p0[1] = ray.getstartP().getpy();
    p0[2] = ray.getstartP().getpz();

    if (ndim > 1 && ndim < 4)
    { //intersect with an ellipse or an elliposoid.
        for (int i = 0; i < ndim; i++)
        {
            fac = 4.0 / (ElliposoidSize[i] * ElliposoidSize[i]);
            a += dir[i] * dir[i] * fac;
            b += dir[i] * p0[i] * fac;
            c += p0[i] * p0[i] * fac;
        }
    }
    else
    {
        std::cout << "IntersectEllipsoid(): " << ndim << " dimension is invalid!" << std::endl;
        std::exit(-1);
    }

    // calc (B/2)^2 - AC : note b=B/2, a=A, c=C
    delta = b * b - a * c;
    if (delta > 0.)
    {
        delta = std::sqrt(delta);
        ray.setTs((-b - delta) / a, (-b + delta) / a);
        return true;
    }
    else
    {
        return false;
    }
}

/****************************************************************************
 * IntersectBBox()
 * Clip ray->min_t and ray->max_t for the intersection of the ray with the
 * bounding box.  The bounding box is assumed to have all coords for first
 * point less than or equal to coords for second point.
 * If ray intersects bbox return 1, otherwise return 0.
 * "bounded" determines whether we literally clip the ray or whether we
 * locate the real bbox intersections.
 * The bounded version is good for determining the length of a line starting
 * or ending at a point in the box and going to the edge.
****************************************************************************/
bool PETScanner::IntersectBBox(Block &box, Ray &ray, int boundindex)
{
    //check for bounded or unbounded search
    if (boundindex)
    {
        ray.setTs(0, ray.getLen());
    }
    else
    { //unbounded search
        ray.setTs(-HUGE, HUGE);
    }

    //get the bound od the box.
    float bound[2][3];
    bound[0][0] = box.getBottomBound().getpx();
    bound[0][1] = box.getBottomBound().getpy();
    bound[0][2] = box.getBottomBound().getpz();
    bound[1][0] = box.getTopBound().getpx();
    bound[1][1] = box.getTopBound().getpy();
    bound[1][2] = box.getTopBound().getpz();
    /* Keep updating the new min and max t values for the line
     * clipping against each component one at a time */

    for (int i = 0; i < 3; i++)
    {
        if (!ClipRay(bound[1][i], bound[0][i], i, ray))
        {
            return false;
        }
    }

    /* Check that all this clipping has not removed the interval,
        i.e. that the line intersects the bounding box. */
    return ray.getMint() < ray.getMaxt();
}
/****************************************************************************
 * ClipRay()
 * Clip ray->min_t and ray->max_t for the intersection of the ray with the
 * planes given by:
 *   x[index] = toplim
 * and
 *   x[index] = botlim.
 ****************************************************************************/
bool PETScanner::ClipRay(float toplim, float botlim, int index, Ray &ray)
{
    float p0, delt_p, t, tmax, tmin;
    p0 = ray.getstartP().getp(index);
    delt_p = ray.getDir().getp(index);

    tmax = ray.getMaxt();
    tmin = ray.getMint();

    if (delt_p > 0.0)
    {
        t = (toplim - p0) / delt_p;
        if (t < tmax)
            tmax = t;
        t = (botlim - p0) / delt_p;
        if (t > tmin)
            tmin = t;
    }
    else if (delt_p < 0.0)
    {
        t = (botlim - p0) / delt_p;
        if (t < tmax)
            tmax = t;
        t = (toplim - p0) / delt_p;
        if (t > tmin)
            tmin = t;
    }
    else
    {
        if ((p0 < botlim) || (p0 > toplim))
            return 0;
    }
    ray.setTs(tmin, tmax);
    /* Check that all this clipping has not removed the interval. */

    return tmin < tmax;
}
/****************************************************************************
 * BackProjRay()
 * Backproject a ray through a volume, adding "val" times the distance
 * through each voxel to each voxel traversed.

 * If the t_TOF pointer is not NULL, TOF weighting is performed,
 * where *t_TOF is the distance along the ray for the center of the Gaussian
 * weighting (measured from the bounding box intersection)
 * and TOF_slim is the time difference cutoff (in units of sigma) to throw out
 * vanishingly small matrix elements (when non-zero).

 * Called just after SetupRayCast()
****************************************************************************/
void PETScanner::BackProjRay(float t_TOF, float val, RayCast &raycast)
{
    float last_t, l_len;
    int stillin, i;

    //double sigma2_TOF = 0.0, t, sigma2, t2_by_sigma2, t2_by_sigma2_limit;

    /* set up a limit of +- TOF_slim sigmas for the Gaussians for TOF (if non-zero) */
    //t2_by_sigma2_limit = TOF_slim * TOF_slim;

    RayCast &rc = raycast; /* local copy of raycast */
    last_t = 0.0;          /* initialize last distance */

    /* setup TOF sigma-squared if we have a t_TOF pointer */

    //    if (t_TOF) {
    //        sigma2_TOF = pet->dist_TOF_fwhm * pet->dist_TOF_fwhm / (8.0 * M_LN2); /* calculate sigma-squared from FWHM */
    //    }

    /* make a non-zero test value unless we are already out */

    stillin = rc.inBuf[0] * rc.inBuf[1] * rc.inBuf[2];

    while (stillin)
    {

        /* look for next intersection: smallest component of nextT */

        i = (rc.nextT[0] <= rc.nextT[1]) ? ((rc.nextT[0] <= rc.nextT[2]) ? 0 : 2) : ((rc.nextT[1] <= rc.nextT[2]) ? 1 : 2);

        l_len = rc.nextT[i] - last_t;

        /* increment voxel */
        /* no TOF, just the line length */
        //MeshIndex mi;
        //find the mesh in the line raycast pass through the image.
        //        int plane = vol.getBlockGrid().getnZ()*vol.getBlockGrid().getnY();
        //        int line = vol.getBlockGrid().getnZ();
        size_t mi = rc.boffs;
        //if the LMEfalg is true,
        if (this->LMEflag)
        { //compute LMFactors
            LMElement lmele;
            lmele.setAddr(mi);
            TOF tof = this->getTOFinfo();
            if (tof.isopen())
            {                                //TOF switch is on.
                float bs = tof.getBinSize(); // get the binsize of TOF.
                float t = t_TOF - 0.5 * (rc.nextT[i] + last_t);
                float sigma2 = tof.getSigma2() + (l_len * l_len + bs * bs) / 12.0;
                //float weight=sqrt(2.0 * M_PI * tof.getSigma2()) * 0.9973 ;
                float t2_by_sigma2 = t * t / sigma2;
                //float t2_by_sigma2=t*t/tof.getSigma2();
                if (tof.getlimit2() <= 0 || t2_by_sigma2 < tof.getlimit2())
                { // the voxel is valid in the TOF range.
                    //lmele.setDist(val*l_len*l_len*exp(-0.5*t2_by_sigma2)/sqrt(2.0*M_PI*sigma2));
                    lmele.setDist(val * l_len * tof.getBinSize() * exp(-0.5 * t2_by_sigma2) / sqrt(2.0 * M_PI * sigma2));
                    //lmele.setDist(val * l_len  * exp(-0.5 * t2_by_sigma2));
                    RTlmevt.addLMElement(lmele);
                }
            }
            else
            {
                lmele.setDist(val * l_len);
                RTlmevt.addLMElement(lmele);
            }
            //bplme.setNvoxels(bplme.getNvoxels()+1);
        }
        else
        { // voxel efficiency MAP
            float oldvalue = this->img.getGridValue(mi);
            //vol.setGridValue(mi,oldvalue+1);
            this->img.setGridValue(mi, oldvalue + val * l_len);
        }

        /* update */

        last_t = rc.nextT[i];        /* set last distance for next pass */
        rc.nextT[i] += rc.deltaT[i]; /* set next intersection distance */
        rc.boffs += rc.deltaBuf[i];  /* buffer offset for next voxel */
        stillin = --rc.inBuf[i];     /* if we go out this goes to zero */
    }
}

///LMRec_ver1 added

void PETScanner::IterRec(const std::string &evtfilename, const std::string & vefilename, const std::string &recfilename,size_t start_index, size_t nIterations)
{
    PETImg VoxelEff = this->getImage();
    FileControl fc;
    fc.ReadImg(vefilename,VoxelEff);
    PETImg PriorWeight = this->getImage();
    PETImg PriorImg = this->getImage();
    float penalty = 0.002;
    //float comfac = 32.6786;
    float comfac = 1;
    const std::string transevtfilename = "trans_" + evtfilename;
	std::fstream fs;
	fs.open(transevtfilename, std::ios::in);
    if (!fs) {//the event file has not be transformed.
		this->TransEvents(evtfilename, transevtfilename);
    }
    this->PreMLEM(transevtfilename, VoxelEff, PriorWeight, PriorImg, penalty);
    this->MLEM(transevtfilename, VoxelEff, PriorWeight, PriorImg, comfac, start_index, nIterations, recfilename);
}

void PETScanner::PreMLEM(const std::string &evtfilename, PETImg &VoxelEff, PETImg &PriorWeight, PETImg &PriorImg, float penalty)
{
    FileControl fc;
    float nEvents = fc.ReadEvtHdr(evtfilename);
    /* if the initialization volume was not read in
     initialize the reconstruction with the average number
     of events per voxel divided by the average efficiency*/
    float ExpectedVoxelMean = nEvents / VoxelEff.sumGDvalue();
    this->img.setAllGridValue(ExpectedVoxelMean);

    //weighting for the prior is given by 'penalty_strength' divided
    //by the 'expected_voxel_mean' in addition to a power of the efficiency

    PriorWeight.setAllGridValue(ExpectedVoxelMean);
    PriorWeight.power(0.0);
    //penalty function

    float MeanWeight = PriorWeight.sumGDvalue() / PriorWeight.getTotalMeshes();
    PriorWeight.multiAllGridValue(penalty / (ExpectedVoxelMean * MeanWeight));

    //initialize the prior image, same size as the recon, and value to be expected voxel mean.
    PriorImg.setAllGridValue(ExpectedVoxelMean);
    //this->setABFFlag(true);
    if (this->getABFFlag())
    {
        this->ABF(VoxelEff);
    }
}

void PETScanner::MLEM(const std::string &transevtfilename, PETImg &VoxelEff, PETImg &PriorWeight, PETImg &PriorImg, float comfac,size_t start_index, size_t nIterations, const std::string &RecImgFile)
{
	FileControl fc;
    Grid3D ImgNext = this->getImage();
	Grid3D ImgNow = this->getImage();
	
	if(start_index>1) {//intialize the image_now with a input image.

		if (ABFflag)
			fc.ReadImg(RecImgFile + "_ABF" + std::to_string(start_index - 1) + ".rec", ImgNow);
		else
            fc.ReadImg(RecImgFile + "_" + std::to_string(start_index - 1) + ".rec", ImgNow);
		this->setImage(ImgNow);
	}
    
    float loglike = 0;

    float step_test;              /* conv value */
    float vavg, vsig, vmin, vmax; /* for the estimate */
    float gavg, gsig, gmin, gmax; /* for the gradient */
    float savg, ssig, smin, smax; /* for the step */

    
    std::clock_t start, end;
    start = clock();
    double diff;
    std::cout << "MLEM reconstruction start..." << std::endl;

    for (size_t iter = start_index; iter < nIterations+start_index; iter++)
    {
        //'ImgNext' will aaccumulate the backprojection of the reciprocal line integrals of 'this->img'.
        ImgNext.setAllGridValue(0);
        //'loglike' will be the loglikelihood of the data for 'this->img'.
        Grid3D temp = this->img * VoxelEff;
        loglike = -(temp.sumGDvalue());
        ImgNow = this->img;

        loglike += EM(ImgNow, ImgNext, comfac, transevtfilename);

        //what's the matter???
        /*ImgNow = ImgNext - VoxelEff;

        gavg = ImgNow.sumGDvalue() / ImgNow.getTotalMeshes();
        Grid3D temp1 = ImgNow * ImgNow;
        gsig = std::sqrt(temp1.sumGDvalue() / temp1.getTotalMeshes() - gavg * gavg);
        ImgNow.findMinMax(gmin, gmax);

        vavg = this->img.sumGDvalue() / this->img.getTotalMeshes();
        Grid3D temp2 = this->img * this->img;
        vsig = std::sqrt(temp2.sumGDvalue() / temp2.getTotalMeshes() - vavg * vavg);
        this->img.findMinMax(vmin, vmax);*/

        // the update step in MLEM.
        if (ABFflag)
        { //insert the ABF (after backprojection filtering) here!
            //fc.SaveImg(RecimgFile+".rec"+std::to_string(iter+1),(ImgNext*this->img)/VoxelEff);
            this->ABF(ImgNext);
        }
        ImgNext = (ImgNext * this->img) / VoxelEff;
        /*ImgNow = ImgNext - this->img;
        savg = ImgNow.sumGDvalue() / ImgNow.getTotalMeshes();
        Grid3D temp3 = ImgNow * ImgNow;
        ssig = std::sqrt(temp3.sumGDvalue() / temp3.getTotalMeshes() - savg * savg);
        ImgNow.findMinMax(smin, smax);*/

        this->img = ImgNext;

        std::cout << "loglike: " << loglike << std::endl;
        std::cout << "Gavg:" << gavg << ", Gsig:" << gsig << ", Gmin:" << gmin << ", Gmax:" << gmax << std::endl;
        std::cout << "Vavg:" << vavg << ", Vsig:" << vsig << ", Vmin:" << vmin << ", Vmax:" << vmax << std::endl;
        std::cout << "Savg:" << savg << ", Ssig:" << ssig << ", Smin:" << smin << ", Smax:" << smax << std::endl;

        /*ImgNext.abs();
        float crel = 0.00001;
        float cabs = 0.0001;
        ImgNext.multiAllGridValue(crel);
        ImgNext.addAllGridValue(cabs);

        ImgNow.abs();
        ImgNow = ImgNow - ImgNext;
        float NOuse;
        ImgNow.findMinMax(NOuse, step_test);
        if (step_test <= 0.0)
        {
            std::cout << "Convergence achived!" << std::endl;
            break;
        }*/

		end = clock();
		diff = (double)(end - start) / CLOCKS_PER_SEC;
		std::cout << "iteration " << iter << " / " << nIterations + start_index - 1 << std::endl;
		std::cout << "cpu time: " << diff << " s,";
		std::cout << "the time remains:" << diff / (iter - start_index + 1) * (nIterations + start_index - iter - 1) << " s." << std::endl;

        if (iter % 1 == 0)
        {
            if (ABFflag)
                fc.SaveImg(RecImgFile + "_ABF" + "_" + std::to_string(iter)+".rec" , this->getImage());
            else
                fc.SaveImg(RecImgFile + "_" + std::to_string(iter) + ".rec", this->getImage());
        }
    }
}

void PETScanner::MAP(const std::string & transevtfilename, PETImg & VoxelEff, PETImg & PriorWeight, PETImg & PriorImg, float comfac, size_t start_index, size_t nIterations, const std::string & RecImgFile)
{
	FileControl fc;
	Grid3D ImgNext = this->getImage();
	Grid3D ImgNow = this->getImage();
	
	float expect_mean = PriorImg.getGridValue(0);
	float gamma_value = 0.05*expect_mean;

	//intialize the image_now with a input image (for starting reconstruction from an existed image)
	if ((start_index-1) > 0) {
		if (ABFflag)
			fc.ReadImg(RecImgFile + "_ABF" + std::to_string(start_index - 1) + ".rec", ImgNow);
		else
			fc.ReadImg(RecImgFile + std::to_string(start_index - 1) + ".rec", ImgNow);
		this->setImage(ImgNow);
	}

	float loglike = 0;
	//float step_test;              /* conv value */
	//float vavg, vsig, vmin, vmax; /* for the estimate */
	//float gavg, gsig, gmin, gmax; /* for the gradient */
	//float savg, ssig, smin, smax; /* for the step */

	//compute time cost
	std::clock_t start, end;
	start = clock();
	double diff;
	std::cout << "MAP reconstruction start..." << std::endl;

	for (size_t iter = start_index; iter < nIterations + start_index; iter++)
	{
		//'ImgNext' will aaccumulate the backprojection of the reciprocal line integrals of 'this->img'.
		ImgNext.setAllGridValue(0);
		//'loglike' will be the loglikelihood of the data for 'this->img'.
		Grid3D temp = this->img * VoxelEff;
		loglike = -(temp.sumGDvalue());
		ImgNow = this->img;
		loglike += EM(ImgNow, ImgNext, comfac, transevtfilename);

		//what's the matter???
		//ImgNow = ImgNext - VoxelEff;
		//gavg = ImgNow.sumGDvalue() / ImgNow.getTotalMeshes();
		//Grid3D temp1 = ImgNow * ImgNow;
		//gsig = std::sqrt(temp1.sumGDvalue() / temp1.getTotalMeshes() - gavg * gavg);
		//ImgNow.findMinMax(gmin, gmax);
		//vavg = this->img.sumGDvalue() / this->img.getTotalMeshes();
		//Grid3D temp2 = this->img * this->img;
		//vsig = std::sqrt(temp2.sumGDvalue() / temp2.getTotalMeshes() - vavg * vavg);
		//this->img.findMinMax(vmin, vmax);
		// the update step in MLEM.
		if (ABFflag)
		{ //insert the ABF (after backprojection filtering) here!
		  //fc.SaveImg(RecimgFile+".rec"+std::to_string(iter+1),(ImgNext*this->img)/VoxelEff);
			this->ABF(ImgNext);
		}
		
		ImgNext = (ImgNext * this->img);
		ImgNext.multiAllGridValue(4 * gamma_value);

		temp.setAllGridValue(0);
		temp = temp - VoxelEff;
		temp.addAllGridValue(gamma_value*expect_mean);
		temp.power(2);
		ImgNext = ImgNext + temp;
		ImgNext.sqrt();
		ImgNext = ImgNext - VoxelEff;
		ImgNext.addAllGridValue(gamma_value*expect_mean);
		ImgNext.divideAllGridValue(2 * gamma_value);

		//ImgNow = ImgNext - this->img;
		//savg = ImgNow.sumGDvalue() / ImgNow.getTotalMeshes();
		//Grid3D temp3 = ImgNow * ImgNow;
		//ssig = std::sqrt(temp3.sumGDvalue() / temp3.getTotalMeshes() - savg * savg);
		//ImgNow.findMinMax(smin, smax);
		this->img = ImgNext;
		//std::cout << "loglike: " << loglike << std::endl;
		//std::cout << "Gavg:" << gavg << ", Gsig:" << gsig << ", Gmin:" << gmin << ", Gmax:" << gmax << std::endl;
		//std::cout << "Vavg:" << vavg << ", Vsig:" << vsig << ", Vmin:" << vmin << ", Vmax:" << vmax << std::endl;
		//std::cout << "Savg:" << savg << ", Ssig:" << ssig << ", Smin:" << smin << ", Smax:" << smax << std::endl;
		//
		//ImgNext.abs();
		//float crel = 0.00001;
		//float cabs = 0.0001;
		//ImgNext.multiAllGridValue(crel);
		//ImgNext.addAllGridValue(cabs);
		//
		//ImgNow.abs();
		//ImgNow = ImgNow - ImgNext;
		//float NOuse;
		//ImgNow.findMinMax(NOuse, step_test);
		//if (step_test <= 0.0)
		//{
		//	std::cout << "Convergence achived!" << std::endl;
		//	break;
		//}

		end = clock();
		diff = (double)(end - start) / CLOCKS_PER_SEC;
		std::cout << "iteration " << iter << " / " << nIterations + start_index-1 << std::endl;
		std::cout << "cpu time: " << diff << " s,";
		std::cout << "the time remains:" << diff / (iter-start_index+1) * (nIterations + start_index - iter - 1) << " s." << std::endl;

		// reconstructed image output
		if ((iter + 1) % 5 == 0)
		{
			if (ABFflag)
				fc.SaveImg(RecImgFile + "_ABF" + std::to_string(iter + 1) + ".rec", this->getImage());
			else
				fc.SaveImg(RecImgFile + std::to_string(iter + 1) + ".rec", this->getImage());
		}
	}
}

void PETScanner::ABF(PETImg &img)
{
    Grid3D abf;
    filter ft;
    Block bok;
    GridSize gs1(4, 4, 4);
    bok.setBlockGrid(gs1);
    abf = ft.Kaiser(bok, 10.4, 2, 2);
    //abf.setAllGridValue(1);
    //VoxelEff.setAllGridValue(1);
    img.conv3(abf);
}

float PETScanner::EM(PETImg &ImgNow, PETImg &ImgNext, float comfac, const std::string &transevtfilename)
{
    float loglike = 0;
    FileControl fc;
    size_t nEvt = fc.ReadEvtHdr(transevtfilename);
    std::ifstream inEvtStream(transevtfilename, std::ios_base::in);
    this->setLMEFlag(true);
    for (size_t iEvt = 0; iEvt < nEvt; iEvt++)
    {
        Event evt;
        fc.ReadEvt(inEvtStream, evt);
        LMEvent lmevt = LMF(evt);
        if (lmevt.getNevents() > 0 && lmevt.getNvoxels() > 0)
        {
            loglike += lmevt.getNevents() * std::log(comfac * this->EM_Projection(ImgNow, ImgNext, lmevt ,evt.getweight()));
        }
    }
    return loglike;
}

LMEvent PETScanner::LMF(Event &evt)
{
    // enter the list-mode events saving mode.
    //this->setLMEFlag(true);
	
    RayTracing(evt);
	LMEvent lme = this->getLMEvent();
	return lme;
}

float PETScanner::EM_Projection(PETImg &ImgNow, PETImg &ImgNext, LMEvent &lmevt, float weight)
{
    float val = 0;
    float lineItgrl = 0;
    for (size_t i = 0; i < lmevt.getNvoxels(); i++)
    {
        LMElement lmele = lmevt.getLMElement(i);
        size_t gs = lmele.getAddr();
        //forward projection
        lineItgrl += ImgNow.getGridValue(gs) * lmele.getDist();
    }
    if (lineItgrl > 0)
    {
        //back projection
        val = (float)lmevt.getNevents() * weight / lineItgrl;
        //val = (float)lmevt.getNevents() * weight;
        for (size_t i = 0; i < lmevt.getNvoxels(); i++)
        {
            LMElement lmele = lmevt.getLMElement(i);
            size_t gs = lmele.getAddr();
            ImgNext.setGridValue(gs, val * lmele.getDist() + ImgNext.getGridValue(gs));
        }
    }
    return lineItgrl;
}

/**********************************************/
/****** class function definition above *******/

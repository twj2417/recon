/*
 * @brief:  this file gives the definition of the RingPET
 * @auther: chengaoyu2013@gmail.com
 * @date:   2017/04/10
*/

#include "RingPET.hh"

/* class---RingPET
 * function definitions
*/

RingPET::RingPET(){
    this->defaultRingPET();
}

RingPET::RingPET(float inR, float outR, float ht, size_t nring, size_t nbok, const Block& bok){
    this->innerR = inR;
    this->outerR = outR;
    this->Height = ht;
    this->nRings = nring;
    this->nBlocks = nbok;
    this->blockorigin = bok;
    this->setBlockOrigin(bok);
    this->initializeScanner();
}

RingPET::~RingPET(){

}

// check if the bound of the scanner is valid.
// give an error message when the blocks are overlapped or out of the cylinder.
// and then use the default setting.
bool RingPET::ValidBound(){
    //get the ring interval.
   float ri = this->getRingInterval();
   //get the block size  of the block origin.
   BlockSize bs = blockorigin.getBlockSize();

   //get the lengh width and height of the block.
   float bsx = bs.getpx();
   float bsy = bs.getpy();
   float bsz = bs.getpz();

   float heightbound;
   //calcualte the height bound of the cylinder
   heightbound = (bsz+ri)*nRings-ri;

   if (heightbound>this->Height){
       std::cout<<"RingPET initialization: the blocks too high! "<<std::endl;
       this->defaultRingPET();
       return false;
   }
   //the half angle of a block in a ring.
   float phi =M_PI/nBlocks;
   float l1 = ((innerR+outerR-bsx)/2);
   if(std::atan2(bsy/2,l1)> phi){//blocks overlapped
       std::cout<<"RingPET initialization: the blocks overlapped! "<<std::endl;
       this->defaultRingPET();
       return false;
   }
   return true;
}

void RingPET::initializeScanner(){
    //
    this->clearCrystalList();
    //should not be changed;
    this->setName("RingPET");
    //check the bound of the blocks.
    if(!this->ValidBound()){
        std::cout<<"Initialization failed! "<<std::endl;
       // std::exit(-1);
    }
    this->fillCrystalList();
}

// calculate the position of crystal and store them in a vector
// the crystals are stored in a sequnce of rings from bottom(z_min) to top(z_max)
// and in a block sequence from the x axis along the anti-clockwise direction.
void RingPET::fillCrystalList(){

    //get the middle radius of the cylinder.
    float middleR = (this->getInnerR()+this->getOuterR())/2;
    //get the ring interval.
    float ring_interval = this->getRingInterval();
    //get the block size  of the block origin.
    BlockSize bs = blockorigin.getBlockSize();
    // get the num of meshes in a block;
    size_t bnx = blockorigin.getBlockGrid().getnX();
    size_t bny = blockorigin.getBlockGrid().getnY();
    size_t bnz = blockorigin.getBlockGrid().getnZ();

    //get the lengh width and height of the block.
    float bsz = bs.getpz();

    //get the meshsize of the block
    Point3D meshsize = blockorigin.calculateInterval();

	int num_rings = this->getNumOfRings();
	float bottom_ring_offset = -(bsz+ring_interval)*(num_rings - 1)/ 2; //the bottom ring z_value.

	float ring_offset = (bsz + ring_interval);

	for (int iRing = 0; iRing < num_rings; iRing++) {
		float bzoffset = bottom_ring_offset + ring_offset*iRing;
		for (size_t jBlock = 0; jBlock< nBlocks; jBlock++) {
			Rotation rot;
			rot.setMatrix(jBlock * 2 * M_PI / nBlocks, Axis::z);
			// generate a new block.
			Block thebok = blockorigin;
			Point3D center((innerR + outerR) / 2, 0, bzoffset);
			// set the center of this block.
			thebok.setBlockCenter(center);
			//calculate the position of the first mesh.
			Point3D firstmesh = thebok.calculateOffset();
			//loop the meshes(crystals) in a block
			for (size_t bi = 0; bi < bnx; bi++) {
				for (size_t bj = 0; bj < bny; bj++) {
					for (size_t bk = 0; bk < bnz; bk++) {
						Point3D meshoffset(bi*meshsize.getpx(), bj*meshsize.getpy(), bk*meshsize.getpz());
						Point3D pt1 = firstmesh + meshoffset;
						Point3D pt = rot.rotate(pt1);
						this->addCrystal(pt);
					}
				}
			}
		}
	}

}

//find the index of the crystal in the PET scanner.
//only be correct for the single mesh block at present!
size_t RingPET::locateCrystal(const Point3D& pt) const{
    size_t inCryst = this->getNumOfCrystals();
    size_t iRing;
    //get the ring interval;
    float rinter = this->getRingInterval();
    //get the block size in Z direction.
    float bsz = this->getBlockOrigin().getBlockSize().getpz();


    if(nRings%2 == 0){ //even number of rings.
        float ptz = pt.getpz();
        if (ptz>0){
            iRing = std::floor(ptz/(rinter+bsz));
        }
        else{
            iRing = std::floor(std::abs(ptz)/(rinter+bsz))+nRings/2;
        }
    }
    else{ //odd nuber of rings.
        float ptz = pt.getpz();
        if(std::abs(ptz)< bsz/2){
            iRing = 0;
        }
        else if(ptz > 0){
            iRing = std::round(ptz/(rinter+bsz));
        }
        else{
            iRing = std::round(std::abs(ptz)/(rinter+bsz))+(nRings-1)/2;
        }
    }
    //find the block and mesh.
    float phi = std::atan2(pt.getpy(),pt.getpx());
    size_t iBlock = this->nBlocks;
    if(phi>=0){
        iBlock = std::floor(phi/(2*M_PI/this->nBlocks));
    }
    else{
        iBlock = std::floor((2*M_PI+phi)/(2*M_PI/this->nBlocks));
    }
    inCryst = iRing*nBlocks+iBlock;
    return inCryst;
}

size_t RingPET::getNumOfCrystals()const {
	return this->nRings*nBlocks*blockorigin.getTotalMeshes();
}

// get the position of a crystal (given the block number and mesh number)
bool RingPET::findCrystal(size_t iBlock, size_t iCryst, cryst& pt)const {
	size_t inCryst = 0, nBlockY,nBlockZ;
	inCryst += iBlock*this->blockorigin.getTotalMeshes();
	nBlockY = this->blockorigin.getBlockGrid().getnY();
	nBlockZ = this->blockorigin.getBlockGrid().getnZ();
	inCryst += (iCryst%nBlockY)*nBlockZ + (iCryst / nBlockY);
    pt = PETScanner::GetCrystal(inCryst);
	return true;
}

bool RingPET::findCrystal( cryst& pt)const {
	float px = pt.getpx();
	float py = pt.getpy();
	float pz = pt.getpz();
	float height = this->getHeight();
	float ring_interval = this->getRingInterval();
	size_t num_rings = this->getNumOfRings();
	size_t num_blocks = this->getNumBlocks();
	const GridSize gs = blockorigin.getBlockGrid();
	size_t bnx = gs.getnX();
	size_t bny = gs.getnY();
	size_t bnz = gs.getnZ();

	float bsz = blockorigin.getBlockSize().getpz();

	float bottom_ring_offset = -(bsz + ring_interval)*(num_rings - 1) / 2; //the top ring z_value.
	float ring_offset = (bsz + ring_interval);

	int ring_index = (pz - bottom_ring_offset) / ring_offset;

	float phi = std::atan2(py, px);
	phi = (phi >= 0) ? phi : phi + 2 * M_PI;

    float phi_offset = 2*M_PI / nBlocks;
	int phi_index = -1;
	if (phi< (phi_offset / 2) || phi>(2*M_PI - phi_offset / 2)) {
		phi_index = 0;
	}
    else{
        phi_index = (phi+phi_offset/2)/phi_offset;
    }

	Rotation rot;
	rot.setMatrix(-phi_index * 2 * M_PI / nBlocks, Axis::z);
	
	Point3D rpt = rot.rotate(pt);
	Block bok = blockorigin;
	float bzoffset = ring_offset*ring_index + bottom_ring_offset;
	Point3D center((innerR + outerR) / 2, 0, bzoffset);

    bok.setBlockCenter(center);

	GridSize gs1 =  bok.LocateIndex(rpt);
    int inCryst = gs1.getnX()*bny*bnz + gs1.getnY()*bnz + gs1.getnZ();
	inCryst = ring_index*nBlocks*blockorigin.getTotalMeshes()+phi_index * blockorigin.getTotalMeshes()+inCryst;
	pt = PETScanner::GetCrystal(inCryst);
	return true;
}


//default setting is the ECAT.
void RingPET::defaultRingPET(){

    this->innerR = 412;
    this->outerR = 442;
    this->Height = 155;
    this->nRings = 24;
    this->nBlocks = 784;
    //the rings are placed tightly.
    this->ringInterval = 0;
    Block bok;
    //set the block size of 30mm*3.3mm*6.25mm.
    Point3D bs(30,3.3,6.25);
    bok.setBlockSize(bs);
    //set the blockgrid of the block
    GridSize bg(1,1,1);
    bok.setBlockGrid(bg);

    //the center is default at (0,0,0).
    this->setBlockOrigin(bok);
    //initialize the scanner and generate the crystal list.
    this->initializeScanner();
}

void RingPET::EfficiencyMap(const std::string &vefilename){
    std::clock_t start,end;
    start = clock();
    double diff;
    int nCrysts = this->getNumOfCrystals();
    double  totalLOR = (double)nCrysts*(nCrysts-1)/2;

    //loop all of the possible LORs to calculate the voxel efficiency map.
    for (int i = 0; i< nCrysts; i++){
        for(int j = i+1; j<nCrysts; j++){
            //get the two crystals of a LOR
            Point3D pt1(this->GetCrystal(i));
            Point3D pt2(this->GetCrystal(j));
            Event ray(pt1,pt2);
            // backproject the ray into the voxel efficiency map.
            this->RayTracing(ray);
        }
        if(i%1000==1){
            end=clock();
            diff=(double)(end-start)/CLOCKS_PER_SEC;
            std::cout<<"cpu time: "<<diff/60<<" minutes,";
            double timetemp = (double)(nCrysts-i)*(nCrysts-i-1)/2;
            std::cout<<"the time remains:"<<diff/60/(totalLOR-timetemp)*timetemp<< "miuntes."<<std::endl;
        }
    }
    FileControl fc1;
    fc1.SaveImg(vefilename,this->getImage());
}

//     void RingPET::LMFactors(const std::string& evtfilename, const std::string &lmffilename) {
// enter the list-mode events saving mode.
//     this->setLMEFlag(true);
//     std::clock_t start,end;
//     start = clock();
//     double diff;
//     FileControl fc;
//     size_t nChords = fc.ReadEvtHdr(evtfilename);

//     //read the event file.
//     std::ifstream evtinfile(evtfilename,std::ios_base::in);

//     // write the lmfactor file.
//     std::ofstream outfile(lmffilename,std::ios_base::out|std::ios_base::binary);

//     size_t nLMEvents = 0;

//     for (int i = 0; i< nChords; i++){
//         //get the two crystals of a LOR
//         Point3D pt1;
//         Point3D pt2;
//         float dist_TOF;
//         fc.ReadEvt(evtinfile,pt1,pt2,dist_TOF,this->getTOFinfo().isopen());
//         size_t ptind1 = (this->locateCrystal(pt1));
//         size_t ptind2 = (this->locateCrystal(pt2));
//         if(ptind1>=this->getNumOfCrystals()||ptind2>=this->getNumOfCrystals()){
//             //std::cout<<"invalid position!"<<std::endl;
//             continue;
//         }
//         Point3D cryst1(this->GetCrystal(ptind1));
//         Point3D cryst2(this->GetCrystal(ptind2));
//         Event evt(cryst1,cryst2);
//         evt.setDistTOF(dist_TOF);
//         //backproject the ray into the voxel efficiency map.
//         RayTracing(evt);
//         //save the LMevt
//         nLMEvents += fc.SaveLMevt(outfile,this->getLMEvent());
//         if(i%10000==1){
//             end=clock();
//             diff=(double)(end-start)/CLOCKS_PER_SEC;
//             std::cout<<"cpu time: "<<diff/60<<" miuntes,";
//             std::cout<<"the time remains:"<<diff/i*(nChords-i)/60<< "miuntes."<<std::endl;
//         }
//     }
//     fc.SaveLMevtHdr(lmffilename,nLMEvents);
//     //exit the list mode events saving mode.
//     this->setLMEFlag(false);
// }

void RingPET::TransEvents(const std::string& evtfilename, const std::string& transevtfilename){
    // enter the list-mode events saving mode.
    std::clock_t start,end;
    start = clock();
    double diff;
    FileControl fc;
    size_t nChords = fc.ReadEvtHdr(evtfilename);
    //read the event file.
    std::ifstream evtinfile(evtfilename,std::ios_base::in);

    // write the lmfactor file.
    std::ofstream outfile(transevtfilename,std::ios_base::out|std::ios_base::binary);

    size_t nEvents = 0;

    for (int i = 0; i< nChords; i++){
        //get the two crystals of a LOR
        Point3D pt1;
        Point3D pt2;
        double dist_TOF;
		size_t iBlock1, iCryst1, iBlock2, iCryst2;
        double time1, time2, weight;
        //fc.ReadEvt(evtinfile, iBlock1, iCryst1, time1, iBlock2, iCryst2, time2);
        fc.ReadEvt(evtinfile, pt1, pt2,this->getTOFinfo().isopen(), weight, time1, time2);
        //fc.ReadEvt(evtinfile, pt1, pt2,this->getTOFinfo().isopen(), time1, time2);
        if (!(this->findCrystal(pt1))) {
            //find the wrong blocks!
            std::cout << i << ", " << iBlock1 << std::endl;
            continue;
        }
        if (!(this->findCrystal(pt2))) {
            //find the wrong blocks!
            std::cout << i << ", " << iBlock2 << std::endl;
            continue;
        }
		if (time1 > time2) {
			cryst temp = pt1;
			pt1 = pt2;
			pt2 = temp;
        }

        Event evt(pt1,pt2);
		double timediff = std::abs(time2 - time1);
        evt.settimeTOF(timediff);
        evt.setweight(weight);
        //backproject the ray into the voxel efficiency map.
        //RayTracing(evt);
        //save the transformed event
        nEvents += fc.SaveEvt(outfile,evt);
        if(i%10000==1){
            end=clock();
            diff=(double)(end-start)/CLOCKS_PER_SEC;
            std::cout<<"cpu time: "<<diff<<" seconds,";
            std::cout<<"the time remains:"<<diff/i*(nChords-i)<< "seconds."<<std::endl;
        }
    }
    fc.SaveEvtHdr(transevtfilename,nEvents);
    //exit the list mode events saving mode.
}





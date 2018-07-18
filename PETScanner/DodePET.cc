/*
 * @brief:  this file gives the definition of the DodePET
 * @auther: chengaoyu2013@gmail.com
 * @date:   2017/04/10
*/

#include "DodePET.hh"

/* class---DodePET
 * function definitions
*/

DodePET::DodePET(){
    this->defaultDodePET();
}

DodePET::DodePET(float inR, float tn, size_t nbok, const Block& bok){
    this->innerR = inR;
    this->thickness = tn;
    this->setNBlocks(nbok);
    this->blockorigin.setBlockGrid(bok.getBlockGrid());
    this->initializeScanner();
}

DodePET::~DodePET(){

}


void DodePET::initializeScanner(){
    setupParas();
    this->clearCrystalList();
    //should not be changed;
    this->setName("DodePET");
    //check the bound of the blocks.
    this->fillCrystalList();
}


void DodePET::setupParas(){
    float inR = this->innerR;

    // edge length of dodecahedron.
    float edgeLength = inR*20/std::sqrt(250+110*std::sqrt(5));

    // circumradius.
    float circumR = edgeLength*(sqrt(3)+sqrt(15))/4;

    // distance of center to edge.
    float edgeDist = edgeLength*(sqrt(5)+3)/4;

    Rotation r;

    for(size_t i = 0;i <=5; i++){              //rotate i*2*pi/5
       r.setMatrix(i*2*M_PI/5, Axis::z);       //rotate
       RotationList.push_back(r);
    }

    r.setMatrix(2*acos(inR/edgeDist),Axis::x); //bottom to low level
    RotationList.push_back(r);

    r.setMatrix(-2*acos(inR/edgeDist),Axis::x); //low level to bottom
    RotationList.push_back(r);

    r.setMatrix(M_PI,Axis::x);                 // X: pi/2
    RotationList.push_back(r);

    //set the blockorigin
    float blockheight = std::sqrt(edgeDist*edgeDist-inR*inR)+std::sqrt(circumR*circumR-inR*inR);
    float blockwidth = std::sqrt(circumR*circumR-edgeDist*edgeDist+blockheight*blockheight);
    float blockthickness = thickness;

    //block size
    Point3D bs(blockwidth,blockheight,blockthickness);
    this->blockorigin.setBlockSize(bs);

    //block center
    Point3D bc(0,-blockheight/2+std::sqrt(edgeDist*edgeDist-inR*inR),0);
    this->blockorigin.setBlockCenter(bc);

    //set the penta bounds
    this->PentaBound.clear();


    Point3D pentavertex(-edgeLength/2,edgeLength/2/std::tan(M_PI/5),0);
    this->PentaBound.push_back(pentavertex);
    for(size_t i = 1;i<5;i++){
        r.setMatrix(i*M_PI*2/5,Axis::z);
        Point3D pentav = r.rotate(pentavertex);
        this->PentaBound.push_back(pentav);
    }
}

const Rotation &DodePET::getRotationMatrix(size_t index) const{
    return RotationList[index];
}

// crystal are calculated on the fly.
void DodePET::fillCrystalList(){
    //there is no crystal saved!
    // all of the positions are calculated on the fly.
}

//find the position of the crystal in the PET scanner.
// this is for the List-mode Factors.
bool DodePET::findCrystal(Point3D& pt, const size_t& blockIndex) const{
    //size_t inCryst = this->getNumOfCrystals();
    Point3D temp;
    Point3D temp1(0,0,innerR+thickness/2);

    if(blockIndex > nBlocks-1)
    {
        std::cout<<"DodePET::findCrystal(), too large blockIndex! "<<std::endl;
        std::exit(-1);
    }

    if(blockIndex == 0){ // block0
        temp = pt+temp1;

        pt = this->getBlockOrigin().LocatePoint(temp);
        if(!IsInsidePenta(pt)){// point is outside the pentagon
            return false;
        }
        if(pt == this->getBlockOrigin().getBottomBound()){ //point is out of the block
            return false;
        }
        pt = pt-temp1;
        return true;
    }
    else if(blockIndex == 11){
        Point3D ori(0,0,0);
        pt = ori-pt;
        temp = pt+temp1;
        pt = this->getBlockOrigin().LocatePoint(temp);

        if(!IsInsidePenta(pt)){// point is outside the pentagon
            return false;
        }
        if(pt == this->getBlockOrigin().getBottomBound()){ //point is out of the block
            return false;
        }
        pt = pt-temp1;
        pt = ori - pt;
        return true;
    }
    else if(blockIndex <= 5){ //block1 ~ block5

        temp = this->getRotationMatrix(blockIndex-1).rotate(pt); // Z: negative
        pt = this->getRotationMatrix(7).rotate(temp);            // X: negative
        temp = pt+temp1;
        pt = this->getRotationMatrix(8).rotate(temp);            // X: -pi/2
        temp = pt;
        pt = this->getBlockOrigin().LocatePoint(temp);
        if(!IsInsidePenta(pt)){// point is outside the pentagon
            return false;
        }
        if(pt == this->getBlockOrigin().getBottomBound()){
            return false;
        }
        temp = this->getRotationMatrix(8).rotate(pt);            // X: pi/2
        pt = temp-temp1;
        temp = this->getRotationMatrix(6).rotate(pt);            // X: positive
        pt = this->getRotationMatrix(6-blockIndex).rotate(temp); // Z: positive
        return true;
    }
    else if(blockIndex <= 10)
    {
        size_t blockIndexTemp = blockIndex-5;
        Point3D ori(0,0,0);
        pt = ori - pt;                                          // central symmetry
        temp = this->getRotationMatrix(blockIndexTemp-1).rotate(pt); // Z: negative
        pt = this->getRotationMatrix(7).rotate(temp);            // X: negative
        temp = pt+temp1;
        pt = this->getRotationMatrix(8).rotate(temp);            // X: -pi/2
        temp = pt;
        pt = this->getBlockOrigin().LocatePoint(temp);

        if(!IsInsidePenta(pt)){// point is outside the pentagon
            return false;
        }

        if(pt == this->getBlockOrigin().getBottomBound()){
            return false;
        }

        temp = this->getRotationMatrix(8).rotate(pt);            // X: pi/2
        pt = temp-temp1;
        temp = this->getRotationMatrix(6).rotate(pt);            // X: positive
        pt = this->getRotationMatrix(6-blockIndexTemp).rotate(temp); // Z: positive
        pt = ori-pt;
        return true;
    }
    else{
        std::cout<<"DodePET::findCrystal(), invalid blockIndex!"<<std::endl;
        std::exit(-1);
    }
}

//this method is for calculating the normalization map.
bool DodePET::findCrystal(Point3D& pt, const size_t& blockIndex,size_t index) const{
    //size_t inCryst = this->getNumOfCrystals();
    Point3D temp;
    Point3D temp1(0,0,innerR+thickness/2);
    Point3D ori(0,0,0);
    if(blockIndex > nBlocks-1)
    {
        std::cout<<"DodePET::findCrystal(), too large blockIndex! "<<std::endl;
        std::exit(-1);
    }
    if(index > this->getBlockOrigin().getTotalMeshes()-1)
    {
        std::cout<<"DodePET::findCrystal(), too large meshIndex! "<<std::endl;
        std::exit(-1);
    }

    pt = this->getBlockOrigin().LocatePoint(index);
    if(!IsInsidePenta(pt)){// point is outside the pentagon
        return false;
    }
    if(blockIndex == 0){ // block0
        pt = pt-temp1;
        return true;
    }
    else if(blockIndex == 11){
        pt = pt-temp1;
        pt = ori - pt;
        return true;
    }
    else if(blockIndex <= 5){ //block1 ~ block5
        temp = this->getRotationMatrix(8).rotate(pt);            // X: pi/2
        pt = temp-temp1;
        temp = this->getRotationMatrix(6).rotate(pt);            // X: positive
        pt = this->getRotationMatrix(6-blockIndex).rotate(temp); // Z: positive
        return true;
    }
    else if(blockIndex <= 10)
    {
        size_t blockIndexTemp = blockIndex-5;
        temp = this->getRotationMatrix(8).rotate(pt);            // X: pi/2
        pt = temp-temp1;
        temp = this->getRotationMatrix(6).rotate(pt);            // X: positive
        pt = this->getRotationMatrix(6-blockIndexTemp).rotate(temp); // Z: positive
        pt = ori-pt;
        return true;
    }
    else{
        std::cout<<"DodePET::findCrystal(), invalid blockIndex!"<<std::endl;
        std::exit(-1);
    }
}


double DodePET::cross(const Point3D &p0, const Point3D &p1, const Point3D &p2) const{
    return (p1.getpx()-p0.getpx())*(p2.getpy()-p0.getpy()) - (p2.getpx()-p0.getpx())*(p1.getpy()-p0.getpy());
}

bool DodePET::IsInsidePenta(const Point3D &pt) const{
    if(cross(this->PentaBound[0],pt,this->PentaBound[1])>0) return false;
    if(cross(this->PentaBound[0],pt,this->PentaBound[4])<0) return false;
    int i = 2,j = 4;
    int line = -1;
    while(i<=j){
        int mid = (i+j)>>1;
        if(cross(this->PentaBound[0],pt,this->PentaBound[mid])>0){
            line = mid;
            j =mid -1;
        }
        else
            i = mid+1;
    }
    return cross(this->PentaBound[line-1],pt,this->PentaBound[line])<0;
}
// Attention!!! this number is not the real number of the crystals.
// It is set to be the meshes of blockorigin multiply the number of blocks.
// just convenient for locate the crystal.
size_t DodePET::getNumOfCrystals() const{
    return this->blockorigin.getTotalMeshes()*this->getNBlocks();
}

//default setting is the ECAT.
void DodePET::defaultDodePET(){

    this->innerR = 150;
    this->nBlocks = 11;
    Block bok;


    //set the block size of 30mm*3.3mm*6.25mm.
    Point3D bs(2,2,30);
    bok.setBlockSize(bs);
    //set the blockgrid of the block
    GridSize bg(1,1,1);
    bok.setBlockGrid(bg);

    //the center is default at (0,0,0).
    this->blockorigin.setBlockGrid(bok.getBlockGrid());
    //initialize the scanner and generate the crystal list.
    this->initializeScanner();
}

void DodePET::EfficiencyMap(const std::string &vefilename){
    std::clock_t start,end;
    start = clock();
    double diff;

    //loop all of the possible LORs to calculate the voxel efficiency map.
    size_t NBlocks = this->getNBlocks();
    size_t NMeshes = this->getBlockOrigin().getTotalMeshes();
    double  totalBlockPairs = (double)NBlocks*(NBlocks-1)/2;
    for(size_t iBlock = 0; iBlock < NBlocks;iBlock++){
        for(size_t jBlock = iBlock+1;jBlock < NBlocks;jBlock++){
            for(size_t iMesh = 0;iMesh< NMeshes;iMesh++){
                for(size_t jMesh = 0;jMesh<NMeshes;jMesh++){
                    Point3D pt1,pt2;
                    if(!(this->findCrystal(pt1,iBlock,iMesh)&&this->findCrystal(pt2,jBlock,jMesh))){
                        continue;
                    }
                    Event evt(pt1,pt2);
                    // backproject the ray into the voxel efficiency map.
                    this->RayTracing(evt);
                }
            }
            end=clock();
            diff=(double)(end-start)/CLOCKS_PER_SEC;
            std::cout<<"cpu time: "<<diff/60<<" minutes,";
            double timetemp = (double)(NBlocks-iBlock-1)*(NBlocks-iBlock-2)/2+jBlock+1;
            std::cout<<"the time remains:"<<diff/60/(totalBlockPairs-timetemp)*timetemp<< "miuntes."<<std::endl;
        }
    }

    FileControl fc1;
    fc1.SaveImg(vefilename,this->getImage());
}
// void DodePET::LMFactors(const std::string& evtfilename, const std::string &lmffilename){
//     // enter the list-mode events saving mode.
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
//         size_t blockindex1,blockindex2;
//         fc.ReadEvt(evtinfile,pt1,pt2,dist_TOF,this->getTOFinfo().isopen(),blockindex1,blockindex2);

//         if(!(this->findCrystal(pt1,blockindex1))){
//             //find the wrong blocks!
//             std::cout<<i<<", "<<blockindex1<<std::endl;
//             continue;
//         }
//         if(!(this->findCrystal(pt2,blockindex2))){
//             //find the wrong blocks!
//             std::cout<<i<<", "<<blockindex2<<std::endl;
//             continue;
//         }
//         Event evt(pt1,pt2);
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

void DodePET::TransEvents(const std::string& evtfilename, const std::string& transevtfilename ){
   // enter the list-mode events saving mode.
   this->setLMEFlag(true);
   std::clock_t start,end;
   start = clock();
   double diff;
   FileControl fc;
   size_t nChords = fc.ReadEvtHdr(evtfilename);
   //read the event file.
   const std::string evtfilename2=evtfilename;
   const std::string evtfilename1="dode_jaszczak-3_atan.txt";
   std::ifstream evtinfile1(evtfilename1,std::ios_base::in);
   std::ifstream evtinfile2(evtfilename2,std::ios_base::in);
   // write the lmfactor file.
   std::ofstream outfile(transevtfilename,std::ios_base::out|std::ios_base::binary);

   size_t nEvents = 0;

   for (int i = 0; i< nChords; i++){
       //get the two crystals of a LOR
       Point3D pt1;
       Point3D pt2;
       float time_TOF,weight;
       size_t blockindex1,blockindex2;
       fc.ReadEvt(evtinfile1,pt1,pt2,time_TOF,this->getTOFinfo().isopen(),weight);
       fc.ReadEvt(evtinfile2,pt1,pt2,time_TOF,this->getTOFinfo().isopen(),blockindex1,blockindex2);

       if(!(this->findCrystal(pt1,blockindex1))){
           //find the wrong blocks!
           std::cout<<i<<", "<<blockindex1<<std::endl;
           continue;
       }
       if(!(this->findCrystal(pt2,blockindex2))){
           //find the wrong blocks!
           std::cout<<i<<", "<<blockindex2<<std::endl;
           continue;
       }
       Event evt(pt1,pt2,weight);
       evt.settimeTOF(time_TOF);
       evt.setweight(weight);
       //backproject the ray into the voxel efficiency map.
       //RayTracing(evt);
       //save the LMevt
       //nLMEvents += fc.SaveLMevt(outfile,this->getLMEvent());
       nEvents += fc.SaveEvt(outfile,evt);
       if(i%10000==1){
           end=clock();
           diff=(double)(end-start)/CLOCKS_PER_SEC;
           std::cout<<"cpu time: "<<diff/60<<" miuntes,";
           std::cout<<"the time remains:"<<diff/i*(nChords-i)/60<< "miuntes."<<std::endl;
       }
   }
   //fc.SaveEvtHdr(transevtfilename,nLMEvents);
   fc.SaveEvtHdr(transevtfilename,nEvents);
   //exit the list mode events saving mode.
   this->setLMEFlag(false);
}




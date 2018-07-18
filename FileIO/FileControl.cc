/*
 * @brief:  this file define the methods in "FileControl.hh"
 * @auther: chengaoyu2013@gmail.com
 * @date:   2017/04/08
*/

#include "FileControl.hh"
#include <iomanip>

/* class--- FileControl
 * function definitions
*/

FileControl::FileControl(){

}


FileControl::~FileControl(){

}

int FileControl::SaveImg(const std::string filename, const Grid3D &img) const{
    std::ofstream outfile(filename, std::ios_base::out);
    if(! outfile.good()){
        std::cerr<<"SaveImg(): Error while open file!"<<std::endl;
        return -1;
    }
    size_t nX = img.getBlockGrid().getnX();
    size_t nY = img.getBlockGrid().getnY();
    size_t nZ = img.getBlockGrid().getnZ();
    for (size_t k=0; k<nZ;k++){
        for(size_t j=0; j<nY; j++){
            for(size_t i=0; i<nX; i++){
                MeshIndex mi(i,j,k);
                outfile<< img.getGridValue(mi)<<std::endl;
            }
        }
    }
    return 1;
}

//to be added later
int FileControl::ReadImg(const std::string filename, Grid3D &img) const{
    std::ifstream infile(filename, std::ios_base::in);
    if(! infile.good()){
        std::cerr<<"ReadImg(): no file "<<filename <<"!"<<std::endl;
        exit(-1);
    }
    size_t nX = img.getBlockGrid().getnX();
    size_t nY = img.getBlockGrid().getnY();
    size_t nZ = img.getBlockGrid().getnZ();
    for (size_t k=0; k<nZ;k++){
        for(size_t j=0; j<nY; j++){
            for(size_t i=0; i<nX; i++){
                //define a mesh index.
                MeshIndex mi(i,j,k);
                float value;
                //read the value of the image voxels.
                infile>>value;
                img.setGridValue(mi,value);
            }
        }
    }
}

int FileControl::SaveLMevt(std::ofstream& lmevtoutfile, const LMEvent &lmevt) const{
    //std::ofstream outfile(filename,std::ios_base::out);
    if(! lmevtoutfile.good()){
        std::cerr<<"SaveLMevt(): Error while open file!"<<std::endl;
        return 0;
    }
    //get the number of events and voxels this LOR pass through.
    size_t nv = lmevt.getNvoxels();
    size_t ne = lmevt.getNevents();

    lmevtoutfile.write((char*)&ne,sizeof(size_t));
    lmevtoutfile.write((char*)&nv,sizeof(size_t));
    size_t addr =0;
    float dist =0;
    //lmevtoutfile<< ne <<std::endl;
    //lmevtoutfile<< nv <<std::endl;
    for(size_t iv = 0 ;iv<nv;iv++){
        //GridSize3D addr = lmevt.getLMElement(iv).getAddr();
        addr = lmevt.getLMElement(iv).getAddr();
        dist = lmevt.getLMElement(iv).getDist();
        lmevtoutfile.write((char*)&addr,sizeof(size_t));
        lmevtoutfile.write((char*)&dist,sizeof(float));
//        lmevtoutfile<<lmevt.getLMElement(iv).getAddr()<<std::endl;
//        lmevtoutfile<<lmevt.getLMElement(iv).getDist()<<std::endl;
    }
    return 1;
}

int FileControl::SaveLMevtHdr(const std::string filename, size_t nEvents) const{
    std::ofstream outheaderfile(filename+".hdr",std::ios_base::out);
    if(!outheaderfile.good()){
        std::cerr<<"SaveLMevtHdr(): Error while open file!"<<std::endl;
        return 0;
    }
    outheaderfile<<nEvents<<std::endl;
    return 1;
}

size_t FileControl::ReadLMevtHdr(const std::string filename) const{
    std::ifstream infile(filename+".hdr",std::ios_base::in);
    size_t nLMEvents;
    if(!infile.good()){
        std::cerr<<"ReadLMevtHdr(): no file "<<filename<<"!"<<std::endl;
        std::exit(-1);
    }
    infile>>nLMEvents;
    return nLMEvents;
}

int FileControl::ReadLMevt(std::ifstream& lmevtinfile, LMEvent& lmevt) const{

    //reset this LMevent.
    lmevt.reset();

    //std::ifstream infile(filename,std::ios_base::in);

    //check if the file is good
    if(!lmevtinfile.good()){
        std::cerr<<"ReadLMevt():Error while open file!"<<std::endl;
        return -1;
    }
    //
    size_t ne,nv;
    size_t gs;
    float dist;

    lmevtinfile.read((char*)&ne, sizeof(size_t));
    lmevtinfile.read((char*)&nv, sizeof(size_t));

    for(size_t iv = 0; iv < nv;iv++){
        lmevtinfile.read((char*)&gs,sizeof(size_t));
        lmevtinfile.read((char*)&dist,sizeof(float));

        LMElement lmele(gs,dist);
        lmevt.addLMElement(lmele);
    }
    return 1;
}

int FileControl::SaveEvtHdr(const std::string filename, size_t nEvents) const{
    std::ofstream outheaderfile(filename+".hdr",std::ios_base::out);
    if(!outheaderfile.good()){
        std::cerr<<"SaveEvtHdr(): Error while open file!"<<std::endl;
        return 0;
    }
    outheaderfile<<nEvents<<std::endl;
    return 1;
}

size_t FileControl::ReadEvtHdr(const std::string filename) const{
    std::ifstream infile(filename+".hdr",std::ios_base::in);
    size_t nEvents;
    if(!infile.good()){
        std::cerr<<"ReadEvtHdr(): Error while open file!"<<std::endl;
        std::exit(-1);
    }
    infile>>nEvents;
    return nEvents;
}

size_t FileControl::SaveEvt(std::ofstream &evtinfile, Event &evt) const{
    if(!evtinfile.good()){
        std::cerr<<"ReadEvt(): Error while open file!"<<std::endl;
        return -1;
    }
    float pos[6];

	Point3D pt1(evt.getRay().getstartP());
	Point3D pt2(evt.getRay().getendP());
    pos[0] = pt1.getpx();
    pos[1] = pt1.getpy();
    pos[2] = pt1.getpz();
    pos[3] = pt2.getpx();
    pos[4] = pt2.getpy();
    pos[5] = pt2.getpz();

    evtinfile << std::setw(20) << pos[0] << std::setw(20) << pos[1] << std::setw(20) << pos[2] << std::setw(20) << pos[3] << std::setw(20) << pos[4] << std::setw(20) << pos[5] << std::setw(20) << evt.gettimeTOF() << std::setw(20) << evt.getweight()<< std::endl;
    return 1;
}

int FileControl::ReadEvt(std::ifstream& evtinfile, Event &evt)const{
    if(!evtinfile.good()){
        std::cerr<<"ReadEvt(): Error while open file!"<<std::endl;
        return -1;
    }
    float pos[6];
    float time_tof, weight;
    evtinfile >> pos[0] >> pos[1] >> pos[2] >> pos[3] >> pos[4] >> pos[5] >> time_tof >>weight;

    Point3D pt1,pt2;
    pt1.setPoint(pos[0],pos[1],pos[2]);
    pt2.setPoint(pos[3],pos[4],pos[5]);
    evt.getRay().setRay(pt1,pt2);
	evt.settimeTOF(time_tof);
    evt.setweight(weight);
    return 1;
}


int FileControl::ReadEvt(std::ifstream& evtinfile, Point3D& pt1, Point3D& pt2 , float& dist_tof, bool TOFflag) const{
    if(!evtinfile.good()){
        std::cerr<<"ReadEvt(): Error while open file!"<<std::endl;
        return -1;
    }
    if(TOFflag){
        float pos[6];
        //read two point from the input file.
        evtinfile >> pos[0] >> pos[1] >> pos[2] >> pos[3] >> pos[4] >>pos[5]>> dist_tof;
        pt1.setPoint(pos[0],pos[1],pos[2]);
        pt2.setPoint(pos[3],pos[4],pos[5]);
    }
    else{
        float pos[6];
        //read two point from the input file.
        evtinfile >> pos[0] >> pos[1] >> pos[2] >> pos[3] >> pos[4] >>pos[5];
        pt1.setPoint(pos[0],pos[1],pos[2]);
        pt2.setPoint(pos[3],pos[4],pos[5]);
        dist_tof = -1;
    }
    return 1;

}

int FileControl::ReadEvt(std::ifstream& evtinfile, Point3D& pt1, Point3D& pt2 , float& dist_tof, bool TOFflag, size_t &blockindex1, size_t &blockindex2) const{
    if(!evtinfile.good()){
        std::cerr<<"ReadEvt(): Error while open file!"<<std::endl;
        return -1;
    }
    if(TOFflag){
        float pos[6];
        //read two point from the input file.
        evtinfile >> pos[0] >> pos[1] >> pos[2] >> pos[3] >> pos[4] >>pos[5]>> dist_tof >> blockindex1>>blockindex2;
        pt1.setPoint(pos[0],pos[1],pos[2]);
        pt2.setPoint(pos[3],pos[4],pos[5]);
    }
    else{
        float pos[6];
        //read two point from the input file.
        evtinfile >> pos[0] >> pos[1] >> pos[2] >> pos[3] >> pos[4] >>pos[5]>>blockindex1>>blockindex2;
        pt1.setPoint(pos[0],pos[1],pos[2]);
        pt2.setPoint(pos[3],pos[4],pos[5]);
        dist_tof = -1;
    }
    return 1;
}

int FileControl::ReadEvt(std::ifstream& evtinfile, Point3D& pt1, Point3D& pt2 , float& dist_tof, bool TOFflag,  float &weight) const{
    if(!evtinfile.good()){
        std::cerr<<"ReadEvt(): Error while open file!"<<std::endl;
        return -1;
    }
    if(TOFflag){
        float pos[6];
        //read two point from the input file.
        evtinfile >> pos[0] >> pos[1] >> pos[2] >> pos[3] >> pos[4] >>pos[5]>> dist_tof >>weight;
        pt1.setPoint(pos[0],pos[1],pos[2]);
        pt2.setPoint(pos[3],pos[4],pos[5]);
    }
    else{
        float pos[6];
        //read two point from the input file.
        evtinfile >> pos[0] >> pos[1] >> pos[2] >> pos[3] >> pos[4] >>pos[5] >>weight;
        pt1.setPoint(pos[0],pos[1],pos[2]);
        pt2.setPoint(pos[3],pos[4],pos[5]);
        dist_tof = -1;
    }
    return 1;
}

int FileControl::ReadEvt(std::ifstream &evtinfile, size_t &iblock1, size_t &icryst1, double &time1, size_t &iblock2, size_t &icryst2, double &time2) const {
	if (!evtinfile.good()) {
		std::cerr << "ReadEvt(): Error while open file!" << std::endl;
		return -1;
	}
	//read two point from the input file.
	evtinfile >> iblock1 >> icryst1 >> time1 >> iblock2 >> icryst2 >> time2;
	return 1;
}

int FileControl::ReadEvt(std::ifstream &evtinfile, Point3D& pt1, Point3D& pt2, bool TOFflag,double &time1, double &time2)const {
    if (!evtinfile.good()) {
		std::cerr << "ReadEvt(): Error while open file!" << std::endl;
		return -1;
	}
	//read two point from the input file.
	float x1, y1, z1, x2, y2, z2;
    if (TOFflag){
	evtinfile >> x1 >> y1 >> z1 >> x2 >> y2 >> z2>>time1 >> time2;
	pt1.setPoint(x1, y1, z1);
	pt2.setPoint(x2, y2, z2);
    }
    else
    {
        evtinfile >> x1 >> y1 >> z1 >> x2 >> y2 >> z2;
        pt1.setPoint(x1, y1, z1);
        pt2.setPoint(x2, y2, z2);
    }
	return 1;
}

int FileControl::ReadEvt(std::ifstream &evtinfile, Point3D& pt1, Point3D& pt2, bool TOFflag, double &weight, double &time1, double &time2 )const {
    if (!evtinfile.good()) {
        std::cerr << "ReadEvt(): Error while open file!" << std::endl;
        return -1;
    }
    //read two point from the input file.
    float x1, y1, z1, x2, y2, z2;
    //Event evt;
    if (TOFflag){
    evtinfile >> x1 >> y1 >> z1 >> x2 >> y2 >> z2 >>time1 >> time2 >> weight;
    pt1.setPoint(x1, y1, z1);
    pt2.setPoint(x2, y2, z2);
    //evt.setweight(weight);
    }
    else
    {
        evtinfile >> x1 >> y1 >> z1 >> x2 >> y2 >> z2 >> weight;
        pt1.setPoint(x1, y1, z1);
        pt2.setPoint(x2, y2, z2);
        time1 = 0;
        time2 = 0;
        //evt.setweight(weight);
    }
    return 1;
}

int FileControl::savebin(const std::string filename, const Grid3D& img ) const{
    std::ofstream outfile(filename, std::ios_base::out|std::ios_base::binary);
    if(! outfile.good()){
        std::cerr<<"SaveImg(): Error while open file!"<<std::endl;
        return -1;
    }
    size_t nX = img.getBlockGrid().getnX();
    size_t nY = img.getBlockGrid().getnY();
    size_t nZ = img.getBlockGrid().getnZ();
    for (size_t k=0; k<nZ;k++){
        for(size_t j=0; j<nY; j++){
            for(size_t i=0; i<nX; i++){
                MeshIndex mi(i,j,k);
                float value =  img.getGridValue(mi);
                outfile.write((char*)&value,sizeof(float));
            }
        }
    }
    return 1;


    //fout.write((char*)&nNum, sizeof(int));
}

int FileControl::readbin(const std::string filename, Grid3D &img){
    std::ifstream infile(filename, std::ios_base::in|std::ios_base::binary);
    if(! infile.good()){
        std::cerr<<"ReadImg(): no file "<<filename <<"!"<<std::endl;
        exit(-1);
    }
    size_t nX = img.getBlockGrid().getnX();
    size_t nY = img.getBlockGrid().getnY();
    size_t nZ = img.getBlockGrid().getnZ();
    for (size_t k=0; k<nZ;k++){
        for(size_t j=0; j<nY; j++){
            for(size_t i=0; i<nX; i++){
                //define a mesh index.
                MeshIndex mi(i,j,k);
                float value;
                //read the value of the image voxels.
                infile.read((char*)&value, sizeof(float));
                img.setGridValue(mi,value);
            }
        }
    }
}

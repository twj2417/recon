/*
 * @brief:  this file declares the files used in the list-mode reconstruction
 * @auther: chengaoyu2013@gmail.com
 * @date:   2017/04/08
*/

#ifndef FILECONTROL_HH
#define FILECONTROL_HH

#include "../GeometryComponents/GeometryComponents.hh"
#include "FileStruct.hh"

#include<fstream>
/*
 * Class declaration:
 * this class control the data and file used in list-mode reconstruction.
*/
class FileControl{
public:
    FileControl();
    virtual ~FileControl();

    //these two methods is used to save and read images.
    virtual int SaveImg(const std::string filename, const Grid3D& img) const;
    virtual int ReadImg(const std::string filename, Grid3D& img) const;

    // save and read List-Mode events.
    virtual int SaveLMevtHdr(const std::string filename,size_t nEvents) const;
    virtual size_t ReadLMevtHdr(const std::string filename) const;
    virtual int SaveLMevt(std::ofstream &lmevtoutfile, const LMEvent& lmevt) const;
    virtual int ReadLMevt(std::ifstream &lmevtinfile, LMEvent& lmevt) const;


    virtual int SaveEvtHdr(const std::string filename, size_t nEvents) const;
    virtual size_t ReadEvtHdr(const std::string filename) const;

    virtual size_t SaveEvt(std::ofstream &evtinfile, Event& evt) const;
    virtual int ReadEvt(std::ifstream &evtinfile, Event& evt) const;
    virtual int ReadEvt(std::ifstream &evtinfile, Point3D& pt1, Point3D& pt2, float &dist_tof, bool TOFflag) const;
    virtual int ReadEvt(std::ifstream &evtinfile, Point3D& pt1, Point3D& pt2, float &dist_tof, bool TOFflag, size_t& blockindex1,size_t &blockindex2) const;
    virtual int ReadEvt(std::ifstream &evtinfile, Point3D& pt1, Point3D& pt2, float &dist_tof, bool TOFflag, float &weight) const;

	virtual int ReadEvt(std::ifstream &evtinfile, size_t &iBlock1, size_t &iCryst1, double &time1, size_t &iBlock2, size_t &iCryst2, double &time2) const;
    virtual int ReadEvt(std::ifstream &evtinfile, Point3D& pt1, Point3D& pt2, bool TOFflag, double &time1, double &time2) const;
    virtual int ReadEvt(std::ifstream &evtinfile, Point3D& pt1, Point3D& pt2, bool TOFflag, double &weight, double &time1, double &time2) const;

    virtual int savebin(const std::string filename,const Grid3D& img) const;
    virtual int readbin(const std::string filename, Grid3D& img);


//    //these two methods is used to save and read events.
//    virtual int SaveEvents(const std::string filename);
//    virtual int ReadEvents(const std::string filename);

//    //these two methods is used to save and read list-mode factors(SystemMatrix).
//    virtual int SaveLMFactors(const std::string filename);
//    virtual int ReadLMFactors(const std::string filename);

};

#endif // FILECONTROL_HH

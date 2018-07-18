#ifndef TOF_HH
#define TOF_HH
/*
 * @brief:  this file contains the time of flight information of backprojection.
 * @auther: chengaoyu2013@gmail.com
 * @date:   2017/05/20
*/
#include<cmath>
#include<cstddef>
#include<iostream>

class TOF{
public:
    TOF();
    TOF(float fwhm, float bs, float range);
    virtual ~TOF();
    virtual void setFWHM(float fwhm);
    virtual float getFWHM()const;
    virtual void setBinSize(float bs);
    virtual float getBinSize() const;
    virtual void setLimit(float range);
    virtual float getlimit2() const; //get the limit-square.
    virtual float getSigma2() const; //get the sigma-square.

    virtual void on();  //turn on the TOF
    virtual void off(); //turn off the TOF
    virtual bool isopen() const;

private:
    float FWHM;     // measurement error (mm) of annihilation point from TOF
    float Sigma2;   // get the sigma-square.
    float binsize;  // least count for TOF measurement digitization.
    float limit;    // the TOF range limits, the value of points out of this limit is invalid.
    bool  flag;     // the switch to turn on and off the TOF.
};

#endif // TOF_HH

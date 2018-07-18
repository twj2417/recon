/*
 * @brief:  definition of class RayCast in "Geometry.hh"
 * @auther: chengaoyu2013@gmail.com
 * @date:   2017/05/24
*/

#include"../GeometryComponents.hh"


/* class---RayCast
 * function definitions
*/
RayCast::RayCast(){

}
RayCast::~RayCast(){

}


/****************************************************************************
 * SetupRayCast()
 * setup for fast processing of line length calculations
 * IMPORTANT: assumes that the ray has already been intersected with
 * the volume bounding box
****************************************************************************/
void RayCast::SetupRayCast(const Block &imgbox, const Ray &ray){
    int xpos, ypos, zpos;

    xpos = SetupRayCastComponent(imgbox, ray, 0);
    ypos = SetupRayCastComponent(imgbox, ray, 1);
    zpos = SetupRayCastComponent(imgbox, ray, 2);

    int nX = imgbox.getBlockGrid().getnX();
    int nY = imgbox.getBlockGrid().getnY();
    /* return buffer offset for first voxel along ray */
    this->boffs = xpos + nX*(ypos + nY*zpos);
}

/****************************************************************************
 * SetupRayCastComponent()
 * compute distance to next crossing, nextT,
 * distance along ray between successive crossings, deltaT,
 * the increment of the buffer offset for a crossing, deltaBuf,
 * the number of crossings which take us out of the imaging volume, inBuf,
 * and return the index of the voxel component
****************************************************************************/

int RayCast::SetupRayCastComponent(const Block& imgbox, const Ray& ray, int comp)
{
    int pos, i, v_count;
    float pt;
    float v_res;

    float r_pdir = ray.getDir().getp(comp);
    float p0 = ray.getstartP().getp(comp);
    float mint = ray.getMint();
    /* compute component of point at intersection of volume bounding box */

    pt = mint * r_pdir + p0;

    /* get local copies of resolution, and dimension */

    v_res = imgbox.calculateInterval().getp(comp);

    v_count = imgbox.getBlockGrid().getn(comp);

    if (r_pdir > 0.0) {

        /* going to the right, so round down */

        pos = pt / v_res;
//        nextT.setp(((pos+1) * v_res - pt) / r_pdir, comp);
//        deltaT.setp(v_res / r_pdir,comp);
//        deltaBuf.setn(1,comp);
//        inBuf.setn(v_count - pos,comp);
        nextT[comp] = ((pos+1) * v_res - pt) / r_pdir;
        deltaT[comp] = v_res / r_pdir;
        deltaBuf[comp] = 1;
        inBuf[comp] = v_count - pos;

    } else if (r_pdir < 0.0) {

        /* going to the left, so round up and subtract 1 */

        pos = v_count - 1 - (int)(v_count - pt / v_res);
//        nextT.setp((pos * v_res - pt) / r_pdir,comp);
//        deltaT.setp(-v_res / r_pdir,comp);
//        deltaBuf.setn(-1,comp);
//        inBuf.setn(pos+1,comp);
        nextT[comp] = (pos * v_res - pt) / r_pdir;
        deltaT[comp] = -v_res / r_pdir;
        deltaBuf[comp] = -1;
        inBuf[comp] = pos + 1;

    } else {

        pos = pt / v_res;
//        nextT.setp(HUGE,comp);
//        deltaT.setp(HUGE,comp);
//        deltaBuf.setn(0,comp);
//        inBuf.setn(1,comp);
        nextT[comp] = HUGE;
        deltaT[comp] = HUGE;
        deltaBuf[comp] = 0;
        inBuf[comp] = 1;
    }

    /* get the correct spacing for the buffer pointer changes */

    for (i=comp-1; i>=0; i--)
        deltaBuf[comp] *= imgbox.getBlockGrid().getn(i);
//        deltaBuf.setn(deltaBuf.getn(comp)*imgbox.getBlockGrid().getn(i),comp);

    return pos;
}

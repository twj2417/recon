#include <iostream>
//#include <armadillo>

#include "PETScanner/PETScanner.hh"
#include "PETScanner/SpherePET.hh"
#include "PETScanner/RingPET.hh"
#include "PETScanner/DodePET.hh"

//#include "ImageProcessing/ImgFilters.hh"

//using namespace arma;
int paraAnalyse(const std::string parafile);
int main(int argc, char *argv[])
{
    if(argc!=3){
        std::cerr<<"Invalid arguments: "<< argc <<" !"<<std::endl;
        std::cerr<<"there should be two arguments: 1) the ring configure txt file;"<<std::endl;
        std::cerr<<"                               2) the task txt file."<<std::endl;
        std::exit(-1);
    }
    std::string pet_configure_file = argv[1];
    std::string task_file = argv[2];
    std::ifstream configure(pet_configure_file, std::ios_base::in);
    float in_radius, out_radius;
    float height, gap;
    int   nBlocks, nRings;
    int   bnx, bny ,bnz;
    float bgx, bgy, bgz;
    std::string key;
    while(!configure.eof()){
       configure >> key;
       if(key == "PETScanner:")
           configure>>key;
       else if(key =="InnerRadius:")
           configure>>in_radius;
       else if(key =="OuterRadius:")
           configure>>out_radius;
       else if(key =="Height:")
           configure>>height;
       else if(key =="RingGap:")
           configure>>gap;
       else if(key == "NumberOfBlocksPerRing:")
           configure>>nBlocks;
       else if(key == "NumberOfRings:")
           configure>>nRings;
       else if(key =="Grid:")
           configure>>bnx>>bny>>bnz;
       else if(key =="BlockSize:")
           configure>>bgx>>bgy>>bgz;
       else{

       }
    }
    configure.close();
    Block bok;
    Point3D bls(bgx,bgy,bgz);
    bok.setBlockSize(bls);
    //set the blockgrid of the block
    GridSize blg(bnx,bny,bnz);
    bok.setBlockGrid(blg);
    Point3D blc(0,0,0);
    bok.setBlockCenter(blc);
    RingPET rpet(in_radius, out_radius, height, nRings, nBlocks, bok);



    std::ifstream task_stream(task_file,std::ios_base::in);
    std::string task;
    int inx,iny,inz,tof_flag,abf_flag;
    float igx,igy,igz;
    std::string map_name,input_file,output_image_name;
    int iteration_start_index,num_iteration;
    float FWHM,binsize,limit;
    while(!task_stream.eof()){
        task_stream>>key;
        if(key == "Name:"){
            std::string name;
            task_stream>> name;
            if(name == "NormalizationMap"){
                task = "Map";
                std::string key;
                while(key != "End"){
                    task_stream>>key;
                    if(key == "MapGrid:"){
                        task_stream >>inx>>iny>>inz;
                    }
                    else if(key== "MapSize:"){
                        task_stream>>igx>>igy>>igz;
                    }
                    else if(key == "MapName:"){
                        task_stream>>map_name;
                    }
                    else{

                    }
                }
            }
            else if (name == "Reconstruction"){
                task = "Recon";
                while(key != "End"){
                    task_stream>>key;
                    if(key == "Grid:"){
                        task_stream >>inx>>iny>>inz;
                    }
                    else if(key== "Size:"){
                        task_stream>>igx>>igy>>igz;
                    }
                    else if(key == "MapName:"){
                        task_stream>>map_name;
                    }
                    else if(key == "IterationStartIndex:"){
                        task_stream>>iteration_start_index;
                    }
                    else if (key == "NumberOfIteration:"){
                        task_stream>>num_iteration;
                    }
                    else if (key == "InputEventFile:"){
                        task_stream>>input_file;
                    }
                    else if(key == "OutputImage:"){
                        task_stream>>output_image_name;
                    }
                    else if (key == "TOF_Flag:"){
                        task_stream>>tof_flag;
                        if(tof_flag == 1){
                            task_stream>>key>>FWHM; // the full half
                            task_stream>>key>>binsize;
                            task_stream>>key>>limit;
                        }

                    }
                    else if (key == "ABF_Flag:"){
                        task_stream>>abf_flag;
                    }
                    else{

                    }

                }
            }
        }
    }
    task_stream.close();

    Grid3D img;
    BlockSize bs(igx,igy,igz);// image size of the reconstruct region
    GridSize gs(inx,iny,inz); // image grid
    Point3D bc(0,0,0);
    img.setBlockSize(bs);
    img.setBlockGrid(gs);
    img.setBlockCenter(bc);


    if(task=="Map"){  //step1: calculate the voxel efficiency map.
        std::cout<<std::endl;
        std::cout<<"************************************************"<<std::endl;
        std::cout<<"STEP 1 start: calculate the VoxelEfficiency Map!"<<std::endl;
        std::cout<<"************************************************"<<std::endl;
        std::cout<<std::endl;

        const std::string vefilename = map_name;

        rpet.setRayTracingMethod(RayTracingMethod::Huesmann);
        rpet.setImage(img);
        rpet.EfficiencyMap(vefilename);
    }


    if(task=="Recon"){
        if (tof_flag){
            TOF tof(FWHM,binsize,limit); //FHWM = 15mm(100ps),BinSize = 3.66mm, limit = 3*sigma , mCT 527ps
            tof.on();
            rpet.setTOFinfo(tof);
        }
        const std::string vefilename = map_name;
        const std::string GlobalDataName = input_file;
        const std::string GlobalFileName = output_image_name;
        const std::string evtfilename = GlobalDataName+".txt";
        //step2: iterative method to reconstruct the image.
        std::cout<<std::endl;
        std::cout<<"****************************************"<<std::endl;
        std::cout<<"STEP 2 start: reconstruct the PET Image!"<<std::endl;
        std::cout<<"****************************************"<<std::endl;
        std::cout<<std::endl;
        rpet.setABFFlag(abf_flag);
        rpet.setRayTracingMethod(RayTracingMethod::Huesmann);
        rpet.setImage(img);
        rpet.IterRec(evtfilename,vefilename, GlobalFileName,iteration_start_index,num_iteration);
    }
	return 1;
}
int paraAnalyse(const std::string parafile){
    //to be added
    return 1;
}


//Dodecahedron
//Block bok;
//Point3D bls(207.3,217.9,20);
//bok.setBlockSize(bls);
//set the blockgrid of the block
//GridSize blg(103,113,3);
//bok.setBlockGrid(blg);
//Point3D blc(0,0,0);
//bok.setBlockCenter(blc);
//DodePET rpet(150,20,11,bok);
// Grid3D img;
// Point3D bs(256,256,20);
// GridSize gs(256,256,20);    //(240,240,28)
// Point3D bc(0,0,0);
// img.setBlockSize(bs);
// img.setBlockGrid(gs);
// img.setBlockCenter(bc);

//SpherePET rpet(150,170,1,2,2,0,M_PI*5/6);


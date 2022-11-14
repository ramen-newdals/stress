#include <iostream>
#include <vector>
#include <tuple>
#include <string>
#include <stdlib.h>
#include "H5Cpp.h"

class meshReader
{
public:
    std::string file_name {"/home/ramen_newdals/Documents/BSL/stress/mesh.h5"};
    std::string group_name {"Mesh"};
    std::string verticies_dataset {"coordinates"};
    std::string cells_dataset {"topology"};
    int RANK{2}, DIM0{216955}, DIM1{3}, DIM2{1262194}, DIM3{4};
    double *verticies = new double[216955*3];
    double *cells = new double[DIM2*DIM3];

    //verticies = double[216955][3];
    void coolSaying()
    {
        std::cout << "Testing" << std::endl;
    }
    int readVerticies()
    {
        H5::H5File file(file_name, H5F_ACC_RDONLY);
        hsize_t dims[2];
        dims[0] = DIM0;
        dims[1] = DIM1;

        H5::DataSpace dataspace = H5::DataSpace(RANK, dims);
        H5::Group solution;
        H5::DataSet dataset;

        file.openFile(file_name, H5F_ACC_RDONLY);
        solution = file.openGroup(group_name);
        dataset = solution.openDataSet(verticies_dataset);
        dataset.read(verticies, H5T_NATIVE_DOUBLE);
        
        dataset.close();
        solution.close();
        dataspace.close();
        file.close();

        return 0;
    }

    int readCells()
    {
        H5::H5File file(file_name, H5F_ACC_RDONLY);
        hsize_t dims[2];
        dims[0] = DIM2;
        dims[1] = DIM3;

        H5::DataSpace dataspace = H5::DataSpace(RANK, dims);
        H5::Group solution;
        H5::DataSet dataset;

        file.openFile(file_name, H5F_ACC_RDONLY);
        solution = file.openGroup(group_name);
        dataset = solution.openDataSet(cells_dataset);
        dataset.read(cells, H5T_NATIVE_DOUBLE);
        
        dataset.close();
        solution.close();
        dataspace.close();
        file.close();

        return 0;
    }

    void printVerticies()
    {
        for(int i = 0; i<DIM0; i++){
            for(int j = 0; j<DIM1; j++){
                std::cout << verticies[(i*DIM1)+j] << " ";
            }
            std::cout << std::endl;
        }
    }
    
    void printCells()
    {
        for(int i = 0; i<DIM2; i++){
            for(int j = 0; j<DIM3; j++){
                std::cout << cells[(i*DIM1)+j] << " ";
            }
            std::cout << std::endl;
        }
    }

};

class shapeFunctions
{
public:
    shapeFunctions();
    ~shapeFunctions();
    
};

class quadrature
{
public:
    quadrature();
    ~quadrature();
    
};

class gradiant
{
public:
    gradiant();
    ~gradiant();
    
};

class linearSystem
{
public:
    linearSystem();
    ~linearSystem();
    
};

int main(void){
    meshReader mesh1;
    mesh1.coolSaying();
    mesh1.readVerticies();
    mesh1.readCells();
    mesh1.printCells();
    //mesh1.printVerticies();
    std::cout << "HDF5 Api is hell to use" << std::endl;
    return 0; // successfully terminated
}

/*
READ IN MESH VALUES:
====================
1.32717 8.77653 57.4528
0.444988 8.84156 57.832


REAL MESH VALUES:
=================
1.32717, 8.77653, 57.4528,
0.444988, 8.84156, 57.832
*/


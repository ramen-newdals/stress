#include <iostream>
#include <vector>
#include <tuple>
#include <string>
#include <stdlib.h>
#include "H5Cpp.h"

// Deffine type that contains a row of the mesh.h5
typedef struct mytype_t{
    std::vector<double> point;
    hvl_t pointsHandle;
} mytype_t;

class meshReader
{
public:
    std::string file_name {"/home/ramen_newdals/Documents/BSL/stress/mesh.h5"};
    std::string group_name {"Mesh"};
    std::string verticies_dataset {"coordinates"};
    std::string cells_dataset {"topology"};
    int RANK{2}, DIM0{216955}, DIM1{3}, DIM2{126194}, DIM3{4};
    double verticies[216955][3];
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



    void printVerticies()
    {
        for(int i = 0; i<DIM0; i++){
            for(int j = 0; j<DIM1; j++){
                std::cout << verticies[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }

    void creaeCells()
    {
        // Crreates cell mapping from memory
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
    mesh1.printVerticies();
    std::cout << "HDF5 Api is hell to use" << std::endl;
    return 0; // successfully terminated

}
#include <iostream>
#include <vector>
#include <tuple>
#include <string>
#include <stdlib.h>
#include "H5Cpp.h"

class meshReader
{
public:
    std::string file_name {"mesh.h5"};
    std::string group_name {"Mesh"};
    std::string verticies_dataset {"coordinates"};
    std::string cells_dataset {"topology"};
    int RANK{2}, DIM0{216955}, DIM1{3}, DIM2{126194}, DIM3{4};
    double verticies[216955][3];

    void coolSaying()
    {
        std::cout << "Testing" << std::endl;
    }
    int readMesh()
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

    void printRdata()
    {
        for(int i = 0; i<DIM0; i++){
            for(int j = 0; j<DIM1; j++){
                std::cout << verticies[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }
};

int main(void){
    meshReader mesh1;
    mesh1.coolSaying();
    mesh1.readMesh();
    mesh1.printRdata();
    std::cout << "HDF5 Api is hell to use" << std::endl;
    return 0; // successfully terminated

}
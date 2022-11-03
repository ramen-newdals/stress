#include <iostream>
#include <string>
#include "H5Cpp.h"
#include "H5DataType.h"

const H5std_string FILE_NAME("results.h5");
const H5std_string GRORUP_NAME("Solution");
const H5std_string DATASET_NAME("u");

const int RANK     = 2;
const int DIM0     = 216955; // size of dataset
const int DIM1     = 3;

int main(void){
	int i, j;
	double rdata[DIM0][DIM1];;

    int n = 1;
    if(*(char*)&n == 1){
        std::cout << "Little Endian" << std::endl;
    }

	try{
		H5::Exception::dontPrint;

		H5::H5File file(FILE_NAME, H5F_ACC_RDONLY);
		hsize_t dims[2];
		dims[0] = DIM0;
		dims[1] = DIM1;

        H5::DataSpace dataspace = H5::DataSpace(RANK, dims);
        H5::Group solution;
        H5::DataSet dataset;

		file.openFile(FILE_NAME, H5F_ACC_RDONLY);
        solution = file.openGroup(GRORUP_NAME);
		dataset = solution.openDataSet(DATASET_NAME);

        dataset.read(rdata, H5T_NATIVE_DOUBLE);
        std::cout << std::endl << "Data in File:" << std::endl;
        for (i = 0; i < DIM0; i++) {
            for (j = 0; j < DIM1; j++)
                std::cout << " " << rdata[i][j];
            std::cout << std::endl;
        }
        std::cout << std::endl;
    
        dataspace.close();
        solution.close();
        dataset.close();
        file.close();

    } // end of try block

    // catch failure caused by the H5File operations
    catch (H5::FileIException error) {
        error.printErrorStack();
        return -1;
    }

    // catch failure caused by the DataSet operations
    catch (H5::DataSetIException error) {
        error.printErrorStack();
        return -1;
    }

    // catch failure caused by the DataSpace operations
    catch (H5::DataSpaceIException error) {
        error.printErrorStack();
        return -1;
    }

    return 0; // successfully terminated

}
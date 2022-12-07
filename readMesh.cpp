#include <iostream>
#include <vector>
#include <tuple>
#include <string>
#include <stdlib.h>
#include <algorithm>
#include <numbers>
#include <omp.h>
// Timming header
#include "hrtime.h"
// Eigen headders
#include <Eigen/Dense>
// HDF5 Headders
#include "H5Cpp.h"

class meshReader
{
public:
    std::string file_name {"/home/ramen_newdals/Documents/BSL/stress/mesh.h5"};
    std::string group_name {"Mesh"};
    std::string verticies_dataset {"coordinates"};
    std::string cells_dataset {"topology"};

    std::string file_name_velocity {"/home/ramen_newdals/Documents/BSL/stress/solutionData.h5"};
    std::string group_name_velocity {"Solution"};
    std::string velocity_dataset {"u"};

    std::vector<std::vector<int>> cell_vector;
    std::vector<std::vector<double>> verticie_vector;
    std::vector<std::vector<double>> velocity_vector;
    std::vector<std::vector<double>> grad_velocity;

    int RANK, verticies_dim0, verticies_dim1, cells_dim0, cells_dim1;
    double *verticies;
    double *velocity;
    int *cells;

    meshReader()
    {
        //TODO: Change these array dimension names to something more meaningfull
        RANK = 0;
        verticies_dim0 = 216955;
        verticies_dim1 = 3;
        cells_dim0 = 1262194, 
        cells_dim1 = 4;
        verticies = new double[verticies_dim0*verticies_dim1];
        velocity = new double[verticies_dim0*verticies_dim1];
        cells = new int[cells_dim0*cells_dim1];

    };
    ~meshReader()
    {
        delete velocity;
        delete verticies;
        delete cells;
    };


    void coolSaying()
    {
        std::cout << "we are in the beam" << std::endl;
    }
    
    std::vector<double> make_verticie(double x, double y, double z)
    {
        std::vector<double> verticie;
        verticie.push_back(x);
        verticie.push_back(y);
        verticie.push_back(z);
        return verticie;
    }

    std::vector<int> make_cell(int v0, int v1, int v2, int v3)
    {
        std::vector<int> cell;
        cell.push_back(v0);
        cell.push_back(v1);
        cell.push_back(v2);
        cell.push_back(v3);
        return cell;
    }

    int readVerticies()
    {
        H5::H5File file(file_name, H5F_ACC_RDONLY);
        hsize_t dims[2];
        dims[0] = verticies_dim0;
        dims[1] = verticies_dim1;

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

        // Read in all verticies to verticie_vector
        for(int i = 0; i<verticies_dim0; i++)
        {verticie_vector.push_back(make_verticie(verticies[(i*3)], verticies[(i*3)+1],verticies[(i*3)+2]));}
        return 0;
    }

    int readCells()
    {
        H5::H5File file(file_name, H5F_ACC_RDONLY);
        hsize_t dims[2];
        dims[0] = cells_dim0;
        dims[1] = cells_dim1;

        H5::DataSpace dataspace = H5::DataSpace(RANK, dims);
        H5::Group solution;
        H5::DataSet dataset;

        file.openFile(file_name, H5F_ACC_RDONLY);
        solution = file.openGroup(group_name);
        dataset = solution.openDataSet(cells_dataset);
        dataset.read(cells, H5T_NATIVE_INT);
        
        dataset.close();
        solution.close();
        dataspace.close();
        file.close();

        // Read in all cells to the cells_vector
        for(int i = 0; i<cells_dim0; i++)
        {cell_vector.push_back(make_cell(cells[(i*4)], cells[(i*4)+1], cells[(i*4)+2], cells[(i*4)+3]));}
        return 0;
    }

    int readVelocity()
    {   
        H5::H5File file(file_name_velocity, H5F_ACC_RDONLY);
        hsize_t dims[2];
        dims[0] = verticies_dim0;
        dims[1] = verticies_dim1;

        H5::DataSpace dataspace = H5::DataSpace(RANK, dims);
        H5::Group solution;
        H5::DataSet dataset;

        file.openFile(file_name_velocity, H5F_ACC_RDONLY);
        solution = file.openGroup(group_name_velocity);
        dataset = solution.openDataSet(velocity_dataset);
        dataset.read(velocity, H5T_NATIVE_DOUBLE);
        
        dataset.close();
        solution.close();
        dataspace.close();
        file.close();

        // Read in all nodal velocities to velocity_vector
        for(int i = 0; i<verticies_dim0; i++)
        {velocity_vector.push_back(make_verticie(velocity[(i*3)], velocity[(i*3)+1],velocity[(i*3)+2]));}
        return 0;
    }

    std::vector<int> getVetexConnectivity(int vertexNum)
    {
        /* 
        For each vertex go throuhg every cell and find what cells contain this vertex
        if the cell has the vertex it is one of its neightboors, add all of the verticies to
        the vertex connectivity vector such that there are not duplicates
        */
        std::vector<std::vector<int>> vertexConnectivity;
        std::vector<int> connectivityVerticies;
        std::vector<int> connectivityCells;
        connectivityVerticies.push_back(vertexNum);
        // Find all cells that contain a common vertex
        int i, j;
        for(i = 0; i<cells_dim0; i++)
        {
            for(j = 0; j<4; j++)
            {
                if(cell_vector[i][j] == vertexNum)
                {
                    connectivityCells.push_back(i);
                }
            }
        }
        // Using the cells that have a common vertex add all nodes of these cells
        // to the vertex connectivity map
        for(i = 0; i<connectivityCells.size(); i++)
        {
            for(j = 0; j<4; j++)
            {
                if(cell_vector[connectivityCells[i]][j] != vertexNum)
                {
                    connectivityVerticies.push_back(cell_vector[connectivityCells[i]][j]);
                }
            }
        }
        // Remove duplicate verticies from the connectivity list
        std::sort(connectivityVerticies.begin(), connectivityVerticies.end()); 
        auto last = std::unique(connectivityVerticies.begin(), connectivityVerticies.end());
        connectivityVerticies.erase(last, connectivityVerticies.end());
        // for (const auto& i : connectivityVerticies) std::cout << i << " ";
        // std::cout << std::endl;
        return connectivityVerticies;
    }

    void constructLSMatricies(int vertexNum)
    {
        std::vector<int> vertexConnectivity = getVetexConnectivity(vertexNum);
        Eigen::MatrixXd A(vertexConnectivity.size(), 10);
        Eigen::MatrixXd A_T(10, vertexConnectivity.size());
        Eigen::MatrixXd A_F(10, 10);
        Eigen::MatrixXd W(vertexConnectivity.size(), vertexConnectivity.size());
        Eigen::VectorXd F(vertexConnectivity.size());
        Eigen::VectorXd b(10);
        // Eigen::ColPivHouseholderQR<MatrixXd> dec(A_F);
        Eigen::VectorXd x(10);
        double avg_length=0;
        int i, j;
        // Calculate the average length between all verticies
        for(i=0;i<vertexConnectivity.size(); i++)
        {
            avg_length += std::sqrt(std::pow(verticie_vector[vertexNum][0]-verticie_vector[vertexConnectivity[i]][0], 2) + 
                                    std::pow(verticie_vector[vertexNum][1]-verticie_vector[vertexConnectivity[i]][1], 2) + 
                                    std::pow(verticie_vector[vertexNum][2]-verticie_vector[vertexConnectivity[i]][2], 2));
        }
        avg_length = avg_length/vertexConnectivity.size();
        // std::cout<<avg_length<<std::endl;

        double const_1 = 0.398942280401432677, d;
        for(i = 0; i<vertexConnectivity.size(); i++)
        {
            F(i) = velocity_vector[vertexConnectivity[i]][0];
            A(i, 0) = std::pow(verticie_vector[vertexConnectivity[i]][0], 2);
            A(i, 1) = std::pow(verticie_vector[vertexConnectivity[i]][1], 2);
            A(i, 2) = std::pow(verticie_vector[vertexConnectivity[i]][2], 2);
            A(i, 3) = verticie_vector[vertexConnectivity[i]][0]*verticie_vector[vertexConnectivity[i]][1]*verticie_vector[vertexConnectivity[i]][2];
            A(i, 4) = verticie_vector[vertexConnectivity[i]][0]*verticie_vector[vertexConnectivity[i]][1];
            A(i, 5) = verticie_vector[vertexConnectivity[i]][0]*verticie_vector[vertexConnectivity[i]][2];
            A(i, 6) = verticie_vector[vertexConnectivity[i]][1]*verticie_vector[vertexConnectivity[i]][2];
            A(i, 7) = verticie_vector[vertexConnectivity[i]][0];
            A(i, 8) = verticie_vector[vertexConnectivity[i]][1];
            A(i, 9) = verticie_vector[vertexConnectivity[i]][2];
            for(j = 0; j<vertexConnectivity.size(); j++)
            {
                d = std::sqrt(std::pow(verticie_vector[vertexConnectivity[j]][0]-verticie_vector[vertexConnectivity[i]][0], 2) + 
                              std::pow(verticie_vector[vertexConnectivity[j]][1]-verticie_vector[vertexConnectivity[i]][1], 2) + 
                              std::pow(verticie_vector[vertexConnectivity[j]][2]-verticie_vector[vertexConnectivity[i]][2], 2));
                W(i, j) = (1/avg_length)*const_1*std::exp((std::pow(d, 2)/std::pow(avg_length, 2))*-1);
            }
        }
        A_T = A.transpose();
        A_F = A_T*W*A;
        b = A_T*W*F;
        x = A_F.colPivHouseholderQr().solve(b);
        std::cout << "Hello From: " << omp_get_thread_num() << std::endl;
        // std::cout << "A is:" << std::endl << A_F << std::endl;
        // std::cout << "b is:" << std::endl << b << std::endl;
        // std::cout << "The solution is:" << std::endl << x << endl;

        // std::cout << "grad_u[0] ~= " << (2*verticie_vector[vertexNum][0]*x[0]) + (x[3]*verticie_vector[vertexNum][1]*verticie_vector[vertexNum][2]) + (x[4]*verticie_vector[vertexNum][1]) + (x[5]*verticie_vector[vertexNum][2]) + (x[7]) << std::endl;
    }

    int solveLSMatricies()
    {
        return 0;
    }

    int estimateGradient()
    {
        return 0;
    }

    int calculateNormal()
    {
        return 0;
    }

    int checkConnectivity()
    {
        for(int i = 0; i<verticies_dim0; i++)
        {
            getVetexConnectivity(i);
        }
        std::cout << "Checked connectivity done." << std::endl;
        return 0;
    }

    void readMesh()
    {
        readVerticies();
        readCells();
        readVelocity();
    }
};

void printThrread()
{
    std::cout << omp_get_thread_num() << std::endl;
}

void constructLSMatricies()
{
    Eigen::MatrixXd A(5, 5);
    Eigen::VectorXd b(5);
    Eigen::VectorXd x(5);
    double avg_length=0;
    int i, j;
    for(i = 0; i<5; i++)
    {
        b(i) = i;
        for(j = 0; j<5; j++)
        {
            A(i, j) = i*j;
        }
    }
    x = A.colPivHouseholderQr().solve(b);
    //std::cout << "Hello From: " << omp_get_thread_num() << std::endl;
}

double readMesh(int numTrials)
{
    double start_time, end_time, avg_time=0;
    for(int i = 0; i<numTrials; i++)
    {
        meshReader mesh1;
        start_time = getElapsedTime();
        mesh1.readMesh();
        end_time = getElapsedTime();
        avg_time += (end_time - start_time);
    }
    return (avg_time/numTrials);
}


std::vector<float> compareParallelSerial()
{
    double start_time, end_time;
    int i;
    std::vector<float> solutionTimes; 
    start_time = getElapsedTime();
    omp_set_num_threads(4); // Use 4 threads for all consecutive parallel regions
    #pragma omp parallel
    for(i = 0; i<8; i++)
    {constructLSMatricies();}
    end_time = getElapsedTime();
    std::cout << "Parallel Sollution Took: " << end_time - start_time << std::endl;
    solutionTimes.push_back((end_time-start_time));
    start_time = getElapsedTime();
    for(int i=0; i<8; i++)
    {constructLSMatricies();}
    end_time = getElapsedTime();
    std::cout << "Serial Sollution Took: " << end_time - start_time << std::endl;
    solutionTimes.push_back((end_time-start_time));
    return solutionTimes;
}

double constructMatriciesParallel(int matrixSize)
{
    double start_time, end_time;
    Eigen::MatrixXd A(matrixSize, matrixSize);
    start_time = getElapsedTime();
    #pragma omp parallel
    for(int i = 0; i<matrixSize; i++)
    {
        for(int j = 0; j<matrixSize; j++)
        {
            A(i, j) = i*j;
        }
    }
    end_time = getElapsedTime();
    run_time += (end_time - start_time)
    return run_time;
}

int main(void)
{
    // Reading Mesh Performance Analysis
    std::cout << "Reading Mesh Took : " << readMesh(10) << "Seconds" << std::endl;    
    std::cout << "==============================================================" << std::endl;
    // Matrix Construction Performance
    
    return 0; // successfully terminated
}
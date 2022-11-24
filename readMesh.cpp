#include <iostream>
#include <vector>
#include <tuple>
#include <string>
#include <stdlib.h>
#include "H5Cpp.h"
// VTK headders
#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkCellType.h>
#include <vtkDataSetMapper.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkTetra.h>
#include <vtkUnstructuredGrid.h>

class meshReader
{
public:
    std::string file_name {"/home/ramen_newdals/Documents/BSL/stress/mesh.h5"};
    std::string group_name {"Mesh"};
    std::string verticies_dataset {"coordinates"};
    std::string cells_dataset {"topology"};
    std::vector<std::vector<int>> cell_vector;
    std::vector<std::vector<double>> verticie_vector;
    meshReader(){};
    ~meshReader(){};
    //TODO: Change these array dimension names to something more meaningfull
    int RANK{2}, verticies_dim0{216955}, verticies_dim1{3}, cells_dim0{1262194}, cells_dim1{4};
    double *verticies = new double[verticies_dim0*verticies_dim1];
    int *cells = new int[cells_dim0*cells_dim1];

    void coolSaying()
    {std::cout << "we are in the beam" << std::endl;}
    
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
        {
            cell_vector.push_back(make_cell(cells[(i*4)], cells[(i*4)+1], cells[(i*4)+2], cells[(i*4)+3]));
        }
        return 0;
    }

    void printVerticies()
    {
        std::cout << "{ ";
        for (int i = 0; i<verticies_dim0; i++)
        {
            std::cout << "(" << verticie_vector[i][0] <<
                        ", " << verticie_vector[i][1] <<
                        ", " << verticie_vector[i][2] <<
                        ")"  << std::endl;
        }
        std::cout << "}; \n";
    }
    
    void printCells()
    {
        std::cout << "{ ";
        for (int i = 0; i<cells_dim0; i++)
        {
            std::cout << "(" << cell_vector[i][0] <<
                        ", " << cell_vector[i][1] <<
                        ", " << cell_vector[i][2] <<
                        ", " << cell_vector[i][3] <<
                        ")"  << std::endl;
        }
        std::cout << "};" << std::endl;
    }

    std::vector<std::tuple<double, double, double>> extractBoundaryMesh()
    {
        // Extracts all boundary nodes and cells
    }

    int plotMesh()
    {
        vtkNew<vtkNamedColors> colors;

        vtkNew<vtkPoints> points;
        points->SetNumberOfPoints(verticies_dim0);

        vtkNew<vtkUnstructuredGrid> unstructuredGrid;

        // Add verticies to unstructured grid
        double point[3];
        for(int i = 0; i<verticie_vector.size(); i++)
        {
            point[0] = verticie_vector[i][0];
            point[1] = verticie_vector[i][1];
            point[2] = verticie_vector[i][2];
            points->SetPoint(i, point);
        }
        
        unstructuredGrid->SetPoints(points);

        points->Delete();

        std::cout << "Added verticies" << std::endl;
        // Add cells to unstructured grid

        vtkNew<vtkCellArray> cells;
        vtkIdType pointId;

        int *outputCellTypes = new int[cells_dim0];
        int outputCellType = VTK_TETRA;
        int nodesPerTet = 4;

        for(int i =  0; i<cell_vector.size(); i++)
        {
            cells->InsertNextCell(nodesPerTet);
            for(int j = 0; j<nodesPerTet; j++)
            {
                pointId = cell_vector[i][j];
                cells->InsertCellPoint(pointId);
            }
            outputCellTypes[i] = outputCellType;
        }

        unstructuredGrid->SetCells(outputCellTypes, cells);

        cells->Delete();

        std::cout << "Added cells" << std::endl;
        // Create a mapper and actor  

        vtkNew<vtkDataSetMapper> mapper;
        mapper->SetInputData(unstructuredGrid);

        std::cout << "Mapper created" << std::endl;

        vtkNew<vtkActor> actor;
        actor->SetMapper(mapper);
        actor->GetProperty()->SetColor(colors->GetColor3d("Cyan").GetData());

        std::cout << "Actor created" << std::endl;

        // Create a renderer, render window, and interactor
        vtkNew<vtkRenderer> renderer;
        vtkNew<vtkRenderWindow> renderWindow;
        renderWindow->SetWindowName("Tetrahedron");
        renderWindow->AddRenderer(renderer);
        vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
        renderWindowInteractor->SetRenderWindow(renderWindow);

        std::cout << "Render Window Made" << std::endl;

        // Add the actor to the scene
        renderer->AddActor(actor);
        renderer->SetBackground(colors->GetColor3d("DarkGreen").GetData());
        renderer->ResetCamera();
        renderer->GetActiveCamera()->Azimuth(-10);
        renderer->GetActiveCamera()->Elevation(-20);

        std::cout << "Camera and actor set up" << std::endl;

        // Render and interact
        renderWindow->Render();
        renderWindowInteractor->Start();

        std::cout << "Rendered starrted" << std::endl;

        return EXIT_SUCCESS;
    }

};

int main(void){
    meshReader mesh1;
    mesh1.coolSaying();
    mesh1.readVerticies();
    mesh1.readCells();
    //mesh1.printCells();
    //mesh1.printVerticies();
    std::cout << "HDF5 Api is hell to use" << std::endl;
    int testing;
    testing = mesh1.plotMesh();
    
    return 0; // successfully terminated
}
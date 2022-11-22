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
    std::vector<std::tuple<int, int, int, int>> cell_vector;
    std::vector<std::tuple<double, double, double>> verticie_vector;
    meshReader(){};
    ~meshReader(){};
    //TODO: Change these array dimension names to something more meaningfull
    int RANK{2}, verticies_dim0{216955}, verticies_dim1{3}, cells_dim0{1262194}, cells_dim1{4};
    double *verticies = new double[verticies_dim0*verticies_dim1];
    int *cells = new int[cells_dim0*cells_dim1];

    void coolSaying()
    {std::cout << "we are in the beam" << std::endl;}
    
    std::tuple<double, double, double> make_verticie(double x, double y, double z)
    {return std::make_tuple(x, y, z);}

    std::tuple<int, int, int, int> make_cell(int v0, int v1, int v2, int v3)
    {return std::make_tuple(v0, v1, v2, v3);}

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
        for(int i = 0; i<verticies_dim0*verticies_dim1; i+=3)
        {verticie_vector.push_back(make_verticie(verticies[i-2], verticies[i-1],verticies[i]));}
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

        std::cout << "Reading cells into vector" << std::endl;

        // Read in all cells to the cells_vector
        for(int i = 0; i<cells_dim0*cells_dim1; i+=4)
        {cell_vector.push_back(make_cell(cells[i-3], cells[i-2], cells[i-1], cells[i]));}
        return 0;
    }

    void printVerticies()
    {
        std::cout << "{ ";
        for (int i = 0; i<verticies_dim0; i++)
        {
            std::cout << "(" << std::get<0>(verticie_vector[i]) <<
                        ", " << std::get<1>(verticie_vector[i]) <<
                        ", " << std::get<2>(verticie_vector[i]) <<
                        ")"  << std::endl;
        }
        std::cout << "}; \n";
    }
    
    void printCells()
    {
        std::cout << "{ ";
        for (int i = 0; i<cells_dim0; i++)
        {
            std::cout << "(" << std::get<0>(cell_vector[i]) <<
                        ", " << std::get<1>(cell_vector[i]) <<
                        ", " << std::get<2>(cell_vector[i]) <<
                        ", " << std::get<3>(cell_vector[i]) <<
                        ")"  << std::endl;
        }
        std::cout << "};" << std::endl;
    }

    std::vector<std::tuple<double, double, double>> extractBoundaryMesh()
    {
        // Extracts all boundary nodes and cells
    }

    void plotMesh()
    {

    }

};



int plotTet()
{
    vtkNew<vtkNamedColors> colors;

    vtkNew<vtkPoints> points;
    points->InsertNextPoint(0, 0, 0);
    points->InsertNextPoint(1, 0, 0);
    points->InsertNextPoint(1, 1, 0);
    points->InsertNextPoint(0, 1, 1);

    points->InsertNextPoint(2, 2, 2);
    points->InsertNextPoint(3, 2, 2);
    points->InsertNextPoint(3, 3, 2);
    points->InsertNextPoint(2, 3, 3);

    // Method 1
    vtkNew<vtkUnstructuredGrid> unstructuredGrid1;
    unstructuredGrid1->SetPoints(points);

    vtkIdType ptIds[] = {0, 1, 2, 3};
    unstructuredGrid1->InsertNextCell(VTK_TETRA, 4, ptIds);

    // Method 2
    vtkNew<vtkUnstructuredGrid> unstructuredGrid2;
    unstructuredGrid2->SetPoints(points);

    vtkNew<vtkTetra> tetra;

    tetra->GetPointIds()->SetId(0, 4);
    tetra->GetPointIds()->SetId(1, 5);
    tetra->GetPointIds()->SetId(2, 6);
    tetra->GetPointIds()->SetId(3, 7);

    vtkNew<vtkCellArray> cellArray;
    cellArray->InsertNextCell(tetra);
    unstructuredGrid2->SetCells(VTK_TETRA, cellArray);

    // Create a mapper and actor
    vtkNew<vtkDataSetMapper> mapper1;
    mapper1->SetInputData(unstructuredGrid1);

    vtkNew<vtkActor> actor1;
    actor1->SetMapper(mapper1);
    actor1->GetProperty()->SetColor(colors->GetColor3d("Cyan").GetData());

    // Create a mapper and actor
    vtkNew<vtkDataSetMapper> mapper2;
    mapper2->SetInputData(unstructuredGrid2);

    vtkNew<vtkActor> actor2;
    actor2->SetMapper(mapper2);
    actor2->GetProperty()->SetColor(colors->GetColor3d("Yellow").GetData());

    // Create a renderer, render window, and interactor
    vtkNew<vtkRenderer> renderer;
    vtkNew<vtkRenderWindow> renderWindow;
    renderWindow->SetWindowName("Tetrahedron");
    renderWindow->AddRenderer(renderer);
    vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
    renderWindowInteractor->SetRenderWindow(renderWindow);

    // Add the actor to the scene
    renderer->AddActor(actor1);
    renderer->AddActor(actor2);
    renderer->SetBackground(colors->GetColor3d("DarkGreen").GetData());
    renderer->ResetCamera();
    renderer->GetActiveCamera()->Azimuth(-10);
    renderer->GetActiveCamera()->Elevation(-20);

    // Render and interact
    renderWindow->Render();
    renderWindowInteractor->Start();

    return EXIT_SUCCESS;
}


int main(void){
    meshReader mesh1;
    mesh1.coolSaying();
    mesh1.readVerticies();
    mesh1.readCells();
    //mesh1.printCells();
    //mesh1.printVerticies();
    std::cout << "HDF5 Api is hell to use" << std::endl;
    mesh1.plotMesh();
    //int testing;
    //testing = plotTet();
    
    return 0; // successfully terminated
}
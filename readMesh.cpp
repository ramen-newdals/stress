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
#include <vtkPointData.h>
#include <vtkProperty.h>
#include <vtkLookupTable.h>
#include <vtkNamedColors.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkTetra.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDoubleArray.h>
#include <vtkColor.h>
#include <vtkNamedColors.h>
#include <vtkGradientFilter.h>

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

    int RANK, verticies_dim0, verticies_dim1, cells_dim0, cells_dim1;
    double *verticies;
    double *velocity;
    int *cells;

    // vtk objects
    // Data primative objects
    vtkNew<vtkPoints> points;
    vtkNew<vtkUnstructuredGrid> unstructuredGrid;
    vtkNew<vtkDoubleArray> velocityArray;
    vtkNew<vtkDoubleArray> velocityMagnitudeArray;
    vtkNew<vtkCellArray> vtkCells;
    vtkNew<vtkDoubleArray> velocityGradientArray;
    vtkNew<vtkDoubleArray> normalsArray;
    vtkNew<vtkDoubleArray> wallShearRateArray;


    // Rendering objects
    vtkNew<vtkLookupTable> lut1;
    vtkNew<vtkDataSetMapper> mapper;
    vtkNew<vtkActor> actor;
    vtkNew<vtkRenderer> renderer;
    vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
    vtkNew<vtkRenderWindow> renderWindow;

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

        // VTK object setup
        points->SetNumberOfPoints(verticies_dim0);
        velocityArray->SetName("u");
        velocityArray->SetNumberOfComponents(3);
        velocityArray->SetNumberOfTuples(verticies_dim0);

        velocityMagnitudeArray->SetName("u_mag");
        velocityMagnitudeArray->SetNumberOfComponents(1);
        velocityMagnitudeArray->SetNumberOfTuples(verticies_dim0);
    };
    ~meshReader()
    {
        delete velocity;
        delete verticies;
        delete cells;
        points->Delete();
        unstructuredGrid->Delete();
        velocityArray->Delete();
        velocityMagnitudeArray->Delete();
        vtkCells->Delete();
        velocityGradientArray->Delete();
        normalsArray->Delete();
        wallShearRateArray->Delete();

        lut1->Delete();
        mapper->Delete();
        actor->Delete();
        renderer->Delete();
        renderWindowInteractor->Delete();
        renderWindow->Delete();
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

    int deffineMesh()
    {
        // Add verticies to unstructured grid
        double point[3], velocity[3];
        for(int i = 0; i<verticie_vector.size(); i++)
        {
            point[0] = verticie_vector[i][0];
            point[1] = verticie_vector[i][1];
            point[2] = verticie_vector[i][2];
            points->SetPoint(i, point);
            velocity[0] = velocity_vector[i][0];
            velocity[1] = velocity_vector[i][1];
            velocity[2] = velocity_vector[i][2];
            double vmag = std::sqrt(std::pow(velocity_vector[i][0], 2)+std::pow(velocity_vector[i][0], 2)+std::pow(velocity_vector[i][0], 2));
            velocityArray->SetTuple(i, velocity);
            velocityMagnitudeArray->SetValue(i, vmag);
        }

        unstructuredGrid->SetPoints(points);
        unstructuredGrid->GetPointData()->AddArray(velocityArray);
        unstructuredGrid->GetPointData()->AddArray(velocityMagnitudeArray);
        unstructuredGrid->GetPointData()->SetActiveScalars("u_mag");
        
        velocityArray->Delete();
        velocityMagnitudeArray->Delete();
        points->Delete();
        
        
        int *outputCellTypes = new int[cells_dim0];
        int outputCellType = VTK_TETRA;
        int nodesPerTet = 4;
        vtkIdType pointId;
        for(int i =  0; i<cell_vector.size(); i++)
        {
            vtkCells->InsertNextCell(nodesPerTet);
            for(int j = 0; j<nodesPerTet; j++)
            {
                pointId = cell_vector[i][j];
                vtkCells->InsertCellPoint(pointId);
            }
            outputCellTypes[i] = outputCellType;
        }

        unstructuredGrid->SetCells(outputCellTypes, vtkCells);
        vtkCells->Delete();
        delete outputCellTypes;
        
        unstructuredGrid->Print(std::cout);
        return 0;
    }

    int renderMesh()
    {
        // Define a lut
        
        lut1->SetHueRange(.667, 0);

        
        mapper->SetInputData(unstructuredGrid);
        mapper->SetLookupTable(lut1);
        mapper->SetColorModeToMapScalars();

        
        actor->SetMapper(mapper);

        // Create a renderer, render window, and interactor
        
        renderer->AddActor(actor);
        renderer->SetBackground(0.187, 0.808, 0.420);

        
        renderWindow->AddRenderer(renderer);
        renderWindow->SetSize(350, 500);
        renderWindow->SetWindowName("eins_oct_siben_gang");

        // Add the actor to the scene
        
        renderer->ResetCamera();
        vtkCamera *camera = renderer->GetActiveCamera();
        camera->Elevation(-80.0);
        camera->OrthogonalizeViewUp();
        camera->Azimuth(135.0);

        
        renderWindowInteractor->SetRenderWindow(renderWindow);
        renderWindow->Render();
        renderWindowInteractor->Initialize();
        renderWindowInteractor->Start();

        return EXIT_SUCCESS;
    }

    int calculateGradiant()
    {
        vtkIndent indent;
        vtkNew<vtkGradientFilter> gradient;
        //gradient->SetInputConnection(unstructuredGrid->GetPointData());
        gradient->AddInputData(unstructuredGrid);
        gradient->SetInputScalars(3, "u");
        gradient->SetResultArrayName("u_grad");
        //gradient->Print(std::cout);
        gradient->Update();
        std::cout << "Gradiaent Computed" << std::endl;
        gradient->GetUnstructuredGridOutput()->GetPointData()->PrintSelf(std::cout, indent);
        //gradient->Print(std::cout);
        gradient->Delete();
        return 0;
    }

    int calculateStrainRateTensor()
    {
        /**********************************************************************
        Calculate strain rate tensor: E = 0.5 * (\nabla u + (\nabla u)^T)
        Calculate wall shear rate vector: tau = -2 * E*n * (1-n^T*n)
        Reference: Matyka et al., http://dx.doi.org/10.1016/j.compfluid.2012.12.018
        **********************************************************************/
        double velocityGradient[9];
        double normal[3];
        double wallShearRate[3];
        
        double normalShear, shearVector[3], strainRateTensor[9];
        for (int i=0; i<verticies_dim0; i++)
        {
            // compute strain rate tensor
            velocityGradientArray->GetTuple(i,velocityGradient);
            for (int j=0; j<3; j++)
            {
                for (int k=0; k<3; k++)
                {
                strainRateTensor[3*j + k] = 0.5 * (velocityGradient[3*j + k] + velocityGradient[3*k + j]);
                }
            }

            // compute shear rate vector and normal projection
            normalsArray->GetTuple(i,normal);
            normalShear = 0.0;
            for (int j=0; j<3; j++)
            {
                shearVector[j] = 0.0;
                for (int k=0; k<3; k++)
                {
                shearVector[j] += strainRateTensor[3*j + k] * normal[k];
                }
                normalShear += shearVector[j] * normal[j];
            }

            // compute wall shear rate
            for (int j=0; j<3; j++)
            {
                // sign due to normals pointing outwards
                wallShearRate[j] = -2.0 * (shearVector[j] - normalShear*normal[j]);
            }
            wallShearRateArray->SetTuple(i,wallShearRate);
        }
        return 0;
    }
};

int main(void){
    meshReader mesh1;
    mesh1.coolSaying();
    mesh1.readVerticies();
    mesh1.readCells();
    mesh1.readVelocity();
    mesh1.coolSaying();
    std::cout << "HDF5 Api is hell to use" << std::endl;
    mesh1.deffineMesh();
    mesh1.calculateGradiant();
    
    return 0; // successfully terminated
}
/*=========================================================================*/
// .NAME gq - Finite element quadrature rule; Approximation of the definite integral of a function, as a weighted sum of function values at specified points.
// .SECTION Description
// ..

#ifndef __gq_h
#define __gq_h

#include "vtkObject.h"
#include "vtkCell.h"
#include "vtkDoubleArray.h"

class gq
{
public:
  // vtkTypeMacro(gq,vtkObject);
  // static gq* New();

  // vtkGetObjectMacro(QuadraturePoints,vtkDoubleArray);
  // vtkGetObjectMacro(QuadratureWeights,vtkDoubleArray);

  // vtkSetMacro(Order,int);
  // vtkGetMacro(Order,int);
  gq()
  {
    this->QuadraturePoints = NULL;
    this->QuadratureWeights = NULL;

    this->Order = 1;
    this->PreviousOrder = 0;
   
    this->CellType = VTK_EMPTY_CELL;
  }

  ~gq()
  {
    if (this->QuadraturePoints)
    {
      this->QuadraturePoints->Delete();
      this->QuadraturePoints = NULL;
    }

    if (this->QuadratureWeights)
    {
      this->QuadratureWeights->Delete();
      this->QuadratureWeights = NULL;
    }
  }
  
  int GetNumberOfQuadraturePoints()
  {
    return this->QuadraturePoints->GetNumberOfTuples();
  }
 
  double* GetQuadraturePoint(vtkIdType id)
  {
    return this->QuadraturePoints->GetTuple(id);
  }
 
  void GetQuadraturePoint(vtkIdType id, double* quadraturePoint)
  {
    this->QuadraturePoints->GetTuple(id,quadraturePoint);
  }
 
  double GetQuadraturePoint(vtkIdType id, int c)
  {
    return this->QuadraturePoints->GetComponent(id,c);
  }
 
  double GetQuadratureWeight(vtkIdType id)
  {
    return this->QuadratureWeights->GetValue(id);
  }
 
  void Initialize(vtkIdType cellType);

  void Initialize(vtkCell* cell)
  {
    this->Initialize(cell->GetCellType());
  }
 
  void Initialize1DGauss();
  void Initialize1DJacobi(int alpha, int beta);
  void ScaleTo01();
 
protected:

  void TensorProductTetra(gq* gauss1D, gq* jacA1D, gq* jacB1D);
 
  vtkDoubleArray* QuadraturePoints;
  vtkDoubleArray* QuadratureWeights;
 
  int Order;
  int QuadratureType;
  vtkIdType CellType;
  int PreviousOrder;

private:  
  // gq(const gq&);  // Not implemented.
  void operator=(const gq&);  // Not implemented.

};

#endif
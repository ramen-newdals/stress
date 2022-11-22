// .NAME vtkvmtkFEAssembler - Base class on which to build to build specific finite element routines.
// .SECTION Description
// Supports common operations needed during setup of a finite element solver run. 

#ifndef __vtkvmtkFEAssembler_h
#define __vtkvmtkFEAssembler_h

#include "vtkObject.h"
#include "vtkPolyData.h"
#include "vtkvmtkSparseMatrix.h"
#include "vtkvmtkDoubleVector.h"
#include "vtkvmtkWin32Header.h"

class VTK_VMTK_DIFFERENTIAL_GEOMETRY_EXPORT vtkvmtkFEAssembler : public vtkObject
{
public:

  vtkTypeMacro(vtkvmtkFEAssembler,vtkObject);

  vtkSetObjectMacro(DataSet,vtkDataSet);
  vtkGetObjectMacro(DataSet,vtkDataSet);

  vtkSetObjectMacro(Matrix,vtkvmtkSparseMatrix);
  vtkGetObjectMacro(Matrix,vtkvmtkSparseMatrix);

  vtkSetObjectMacro(RHSVector,vtkvmtkDoubleVector);
  vtkGetObjectMacro(RHSVector,vtkvmtkDoubleVector);

  vtkSetObjectMacro(SolutionVector,vtkvmtkDoubleVector);
  vtkGetObjectMacro(SolutionVector,vtkvmtkDoubleVector);

  vtkGetMacro(NumberOfVariables,int);

  vtkSetMacro(QuadratureOrder,int);
  vtkGetMacro(QuadratureOrder,int);

  virtual void Build() = 0;

  void DeepCopy(vtkvmtkFEAssembler *src);
  void ShallowCopy(vtkvmtkFEAssembler *src);

protected:
  vtkvmtkFEAssembler();
  ~vtkvmtkFEAssembler();

  void Initialize(int numberOfVariables);

  vtkDataSet* DataSet;
  vtkvmtkSparseMatrix* Matrix;
  vtkvmtkDoubleVector* RHSVector;
  vtkvmtkDoubleVector* SolutionVector;

  int NumberOfVariables;
  int QuadratureOrder;

private:
  vtkvmtkFEAssembler(const vtkvmtkFEAssembler&);  // Not implemented.
  void operator=(const vtkvmtkFEAssembler&);  // Not implemented.
};

#endif

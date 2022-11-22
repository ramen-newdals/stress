#include "vtkvmtkFEAssembler.h"
#include "vtkvmtkGaussQuadrature.h"
#include "vtkvmtkFEShapeFunctions.h"
#include "vtkObjectFactory.h"

vtkvmtkFEAssembler::vtkvmtkFEAssembler()
{
  this->DataSet = NULL;
  this->Matrix = NULL;
  this->RHSVector = NULL;
  this->SolutionVector = NULL;
  this->NumberOfVariables = 1;
  this->QuadratureOrder = 1;
}

vtkvmtkFEAssembler::~vtkvmtkFEAssembler()
{
  if (this->DataSet)
    {
    this->DataSet->Delete();
    this->DataSet = NULL;
    }
  if (this->Matrix)
    {
    this->Matrix->Delete();
    this->Matrix = NULL;
    }
  if (this->RHSVector)
    {
    this->RHSVector->Delete();
    this->RHSVector = NULL;
    }
  if (this->SolutionVector)
    {
    this->SolutionVector->Delete();
    this->SolutionVector = NULL;
    }
}

void vtkvmtkFEAssembler::Initialize(int numberOfVariables)
{
  this->NumberOfVariables = numberOfVariables;

  this->Matrix->AllocateRowsFromDataSet(this->DataSet,this->NumberOfVariables);

  this->RHSVector->Allocate(this->DataSet->GetNumberOfPoints(),this->NumberOfVariables);
  this->RHSVector->Fill(0.0);

  this->SolutionVector->Allocate(this->DataSet->GetNumberOfPoints(),this->NumberOfVariables);
  this->SolutionVector->Fill(0.0);
}

void vtkvmtkFEAssembler::DeepCopy(vtkvmtkFEAssembler *src)
{
  this->DataSet->DeepCopy(src->DataSet);
  this->Matrix->DeepCopy(src->Matrix);
  this->RHSVector->DeepCopy(src->RHSVector);
  this->SolutionVector->DeepCopy(src->SolutionVector);
  this->NumberOfVariables = src->NumberOfVariables;
  this->QuadratureOrder = src->QuadratureOrder;
}

void vtkvmtkFEAssembler::ShallowCopy(vtkvmtkFEAssembler *src)
{
  this->DataSet = src->DataSet;
  this->DataSet->Register(this);
  this->Matrix = src->Matrix;
  this->Matrix->Register(this);
  this->RHSVector = src->RHSVector;
  this->RHSVector->Register(this);
  this->SolutionVector = src->SolutionVector;
  this->SolutionVector->Register(this);
  this->NumberOfVariables = src->NumberOfVariables;
  this->QuadratureOrder = src->QuadratureOrder;
}


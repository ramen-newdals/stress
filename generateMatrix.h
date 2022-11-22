// .NAME vtkvmtkUnstructuredGridFEGradientAssembler - Construct a gradient based finite element calculation on a mesh.
// .SECTION Description
// ..

#ifndef __vtkvmtkUnstructuredGridFEGradientAssembler_h
#define __vtkvmtkUnstructuredGridFEGradientAssembler_h

#include "vtkvmtkFEAssembler.h"

class VTK_VMTK_DIFFERENTIAL_GEOMETRY_EXPORT vtkvmtkUnstructuredGridFEGradientAssembler : public vtkvmtkFEAssembler
{
public:

  static vtkvmtkUnstructuredGridFEGradientAssembler* New();
  vtkTypeMacro(vtkvmtkUnstructuredGridFEGradientAssembler,vtkvmtkFEAssembler);

  virtual void Build() override;

  vtkSetStringMacro(ScalarsArrayName);
  vtkGetStringMacro(ScalarsArrayName);

  vtkSetMacro(ScalarsComponent,int);
  vtkGetMacro(ScalarsComponent,int);

  vtkSetMacro(Direction,int);
  vtkGetMacro(Direction,int);

  vtkSetMacro(AssemblyMode,int);
  vtkGetMacro(AssemblyMode,int);
  void SetAssemblyModeToGradient()
  { this->SetAssemblyMode(VTKVMTK_GRADIENTASSEMBLY); }
  void SetAssemblyModeToPartialDerivative()
  { this->SetAssemblyMode(VTKVMTK_PARTIALDERIVATIVEASSEMBLY); }

//BTX
  enum {
    VTKVMTK_GRADIENTASSEMBLY,
    VTKVMTK_PARTIALDERIVATIVEASSEMBLY
  };
//ETX

protected:
  vtkvmtkUnstructuredGridFEGradientAssembler();
  ~vtkvmtkUnstructuredGridFEGradientAssembler();

  void BuildGradient();
  void BuildPartialDerivative();

  char* ScalarsArrayName;
  int ScalarsComponent;
  int AssemblyMode;
  int Direction;

private:
  vtkvmtkUnstructuredGridFEGradientAssembler(const vtkvmtkUnstructuredGridFEGradientAssembler&);  // Not implemented.
  void operator=(const vtkvmtkUnstructuredGridFEGradientAssembler&);  // Not implemented.
};

#endif


#ifndef MYDETECTORCONSTRUCTION_HH
#define MYDETECTORCONSTRUCTION_HH

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;

class MyDetectorConstruction : public G4VUserDetectorConstruction {

public:
  
  // Constructor
  MyDetectorConstruction();
  
  // Destructor
  virtual ~MyDetectorConstruction();
  
  // Method
  virtual G4VPhysicalVolume* Construct();
  
private:

  // Helper methods
  void DefineMaterials();
  void SetupGeometry();
  void SetupScoring(G4LogicalVolume*);
  
  // World logical and physical volumes
  G4LogicalVolume*   fpWorldLogical;
  G4VPhysicalVolume* fpWorldPhysical;
  
};

#endif


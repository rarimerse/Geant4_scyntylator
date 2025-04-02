#ifndef MYSCOREPARAMETERISATION_HH
#define MYSCOREPARAMETERISATION_HH

#include "G4VPVParameterisation.hh"


class MyScoreParameterisation : public G4VPVParameterisation {

public:
  
  // Constructor
   MyScoreParameterisation(const G4double, const G4double, const G4int);
  
  // Destructor
  virtual ~MyScoreParameterisation();
  
  // ComputeDimensions
  virtual  void ComputeDimensions(G4Sphere &, 
                                  const G4int, 
       				  const G4VPhysicalVolume *) const; 
  // ComputeTransformation
  void ComputeTransformation(const G4int, G4VPhysicalVolume*)const {}; 
private:
  G4double   finnerRadius;
  G4double   fouterRadius;
  G4int      fnRings;
};
#endif

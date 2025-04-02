#include "MyScoreParameterisation.hh"
#include "G4Sphere.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

MyScoreParameterisation::MyScoreParameterisation(
                         const G4double innerRadius, 
			 const G4double outerRadius,
			 const G4int nRings)
  :finnerRadius(innerRadius)
  ,fouterRadius(outerRadius)
  ,fnRings(nRings)
{}

MyScoreParameterisation::~MyScoreParameterisation() {}

void MyScoreParameterisation::ComputeDimensions(
                         G4Sphere& sphere, 
		       const G4int copyNo, 
		       const G4VPhysicalVolume*) const
{
  G4double deltaTheta = 90.0/fnRings*deg;
  G4double startTheta = copyNo*deltaTheta;
  sphere.SetInnerRadius(finnerRadius);
  sphere.SetOuterRadius(fouterRadius);
  sphere.SetStartPhiAngle(0.*deg);
  sphere.SetDeltaPhiAngle(360.*deg);
  sphere.SetStartThetaAngle(startTheta);
  sphere.SetDeltaThetaAngle(deltaTheta);
}

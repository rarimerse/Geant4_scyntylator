
#include "MyPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"
#include "G4OpticalPhoton.hh"

MyPrimaryGeneratorAction::MyPrimaryGeneratorAction()
{
  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle);
  
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  
  G4ParticleDefinition* particle = particleTable->FindParticle(particleName="proton");
  //G4ParticleDefinition* particle=G4OpticalPhoton::OpticalPhotonDefinition();


  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(-1.,0.,0.));
  particleGun->SetParticleEnergy(100.*GeV);
}

MyPrimaryGeneratorAction::~MyPrimaryGeneratorAction()
{
  delete particleGun;
}

void MyPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  G4double z_min = -5. * cm;
  G4double z_max = 0. * cm;
  G4double y_min = -4. * cm;
  G4double y_max = 4. * cm;

  G4double z = G4UniformRand() * (z_max - z_min) + z_min;
  G4double y = G4UniformRand() * (y_max - y_min) + y_min;
  G4double x = 0. * cm;
  //G4double phi=twopi*G4UniformRand()*radian;
  //G4double cosTheta=-1.0+2.0*G4UniformRand();
  //G4double sinTheta=sqrt(1.0-cosTheta*cosTheta);

particleGun->SetParticlePosition(G4ThreeVector(x,y,z));
  //particleGun->SetParticleMomentumDirection(G4ThreeVector(sinTheta*cos(phi), sinTheta*sin(phi),cosTheta));
  particleGun->GeneratePrimaryVertex(anEvent);
}

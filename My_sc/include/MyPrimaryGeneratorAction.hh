#ifndef MYPRIMARYGENERATORACTION_HH
#define MYPRIMARYGENERATORACTION_HH

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;

class MyPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {

public:

  // Constructor
  MyPrimaryGeneratorAction();    

  // Destructor
  virtual ~MyPrimaryGeneratorAction();
  
  // Method
  void GeneratePrimaries(G4Event*);

private:

  // Data member
  G4ParticleGun* particleGun;	 

};

#endif

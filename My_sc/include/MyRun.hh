#ifndef MYRUN_HH
#define MYRUN_HH

#include "G4Run.hh"
#include "G4THitsMap.hh"
class G4Event;
class MyRun : public G4Run {

public:

  // Constructor 
  MyRun();

  // Destructor
  virtual ~MyRun();

    //Method
  void RecordEvent(const G4Event*);
  //Data
  G4int fMapId;
  G4THitsMap<G4double> frunHitsMap;
  G4int fCollectionId;
  G4double photonZ;
  G4double photonY;
};

#endif
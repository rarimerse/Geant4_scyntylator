#ifndef MYRUNACTION_HH
#define MYRUNACTION_HH

#include "G4UserRunAction.hh"
class G4Run;
class MyRunAction : public G4UserRunAction {

public:

  // Constructor
  MyRunAction();

  // Destructor
  virtual ~MyRunAction();

// Method
   void BeginOfRunAction(const G4Run* );
   void EndOfRunAction(const G4Run* );
   G4Run* GenerateRun();
};

#endif
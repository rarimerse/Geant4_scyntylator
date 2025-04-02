#include "MyRunAction.hh"
#include "G4Run.hh"
#include "globals.hh"
#include "MyRun.hh"
#include "MyAnalysis.hh"
MyRunAction::MyRunAction() 
{}

MyRunAction::~MyRunAction() 
{}

void MyRunAction::BeginOfRunAction(const G4Run* aRun)
{
 G4cout << "MyLog:   Begin of run" << G4endl;
 // Create analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetFirstHistoId(1);
 
// Open an output file
  //analysisManager->OpenFile("hist.root");

// Creating histograms
 //G4int idPhotonFlux = analysisManager->CreateH1("PhotonFlux","photon flow",15,0.,90.);
   analysisManager->SetFirstHistoId(6);
  analysisManager->OpenFile("Histograms.root");
  analysisManager->CreateH2("Photon Z and time", "Photon Z and time", 10, -5., 0., 25, 0., 50.);
  analysisManager->CreateH2("Photon Y and time", "Photon Y and time", 10, -4., 4., 25, 0., 50.);
  analysisManager->SetFirstNtupleId(1);
  analysisManager->OpenFile("tree.root");
  analysisManager->CreateNtuple("PhotonCollection","Kolekcja fotonow");
  analysisManager->CreateNtupleDColumn("E");
  analysisManager->CreateNtupleDColumn("X");
  analysisManager->CreateNtupleDColumn("Y");
  analysisManager->CreateNtupleDColumn("Z");
  analysisManager->CreateNtupleDColumn("t");
  analysisManager->CreateNtupleDColumn("lambda");
  analysisManager->FinishNtuple();

}

void MyRunAction::EndOfRunAction(const G4Run* aRun)
{
    // Get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();
 G4cout <<"MyLog:   Number of processed events:" << aRun->GetNumberOfEvent()
               << " events. " <<G4endl;
 G4cout << "MyLog:   End of run" << G4endl;
}

G4Run* MyRunAction::GenerateRun()
{
 return new MyRun();
}
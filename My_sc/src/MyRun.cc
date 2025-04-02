#include "MyRun.hh"
#include "G4Event.hh"
#include "globals.hh"
#include "G4TrajectoryContainer.hh"
#include "G4VTrajectory.hh"
#include <iomanip>
#include "G4SDManager.hh"
#include "G4THitsMap.hh"
#include "G4THitsCollection.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "MyPhotonHit.hh"
#include "MyAnalysis.hh"

MyRun::MyRun()
{
    G4SDManager* manager = G4SDManager::GetSDMpointer();
fMapId = manager->GetCollectionID("MyDetector/MyScorer");
G4cout << "MyLog:  MyRun constructor: index of photon scorer map: " << fMapId << G4endl;

fCollectionId =  manager->GetCollectionID("Photon/PhotonCollection");
 G4cout << "MyLog:  MyRun constructor: index of photon collection : " << fCollectionId << G4endl;
}

MyRun::~MyRun()
{}

void MyRun::RecordEvent(const G4Event* evt)
{
 numberOfEvent++;
G4HCofThisEvent* hce = evt->GetHCofThisEvent();

G4int eventID = evt->GetEventID(); 
    G4cout << "MyLog: Processing event number: " << eventID << G4endl;
 
  G4THitsMap<double>* hitsMap =  (G4THitsMap<double>*) (hce->GetHC(fMapId));
  G4cout << "MyLog: number of entries " << hitsMap->entries() << G4endl;
     frunHitsMap += *hitsMap;

    G4THitsCollection<MyPhotonHit>* photonCollection =
      dynamic_cast<G4THitsCollection<MyPhotonHit>*> (hce->GetHC(fCollectionId));

   G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
   G4double num=0;
   G4double t_sum=0;
   G4PrimaryVertex* primaryVertex = evt->GetPrimaryVertex(0);
  photonY = primaryVertex->GetY0()/cm;
  photonZ = primaryVertex->GetZ0()/cm;
  for( G4int i = 0; i< photonCollection->entries(); i++) {
     float e = (*photonCollection)[i]->GetEnergy()/MeV;
     float dlugosc_fali=1239.8419/e;
     float x = (*photonCollection)[i]->GetPosition().x()/cm;
     float y = (*photonCollection)[i]->GetPosition().y()/cm;
     float z = (*photonCollection)[i]->GetPosition().z()/cm;
     float t = (*photonCollection)[i]->GetTime();



     analysisManager->FillNtupleDColumn(0,e);
     analysisManager->FillNtupleDColumn(1,x);
     analysisManager->FillNtupleDColumn(2,y);
     analysisManager->FillNtupleDColumn(3,z);
     analysisManager->FillNtupleDColumn(4,t);
     analysisManager->FillNtupleDColumn(5,dlugosc_fali);
     analysisManager->AddNtupleRow();
     num+=1;
     t_sum+=t;
  }
   G4double t_mean=t_sum/num;
   analysisManager->FillH2(6, photonZ,t_mean); 
   analysisManager->FillH2(7, photonY,t_mean); 
   /*
  for( G4int i = 0; i< photonCollection->entries(); i++)
   G4cout     << "MyLog: energy and position of photon: "
              << std::setprecision(3) << (*photonCollection)[i]->GetEnergy()/MeV << "     "
              << std::setprecision(3) << (*photonCollection)[i]->GetPosition().x()/cm << " "
              << std::setprecision(3) << (*photonCollection)[i]->GetPosition().y()/cm << " "
              << std::setprecision(3) << (*photonCollection)[i]->GetPosition().z()/cm << "      "
              << std::setprecision(3) << photonZ <<"  " << photonY << "      "
              << std::setprecision(3) << (*photonCollection)[i]->GetTime() << G4endl;
  std::map<int,double*>::iterator iter;
  iter = hitsMap->GetMap()->begin();
  while( iter != hitsMap->GetMap()->end() )
   {
   G4cout<< "MyLog:      value of HitsMap: "<< *iter->second << " index of value: "<< iter->first << G4endl;
   iter++;
   }
   */
}
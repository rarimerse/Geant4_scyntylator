#include "MyDetectorConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4PVParameterised.hh"
#include "G4PVPlacement.hh"
#include "G4RunManager.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4Trap.hh"
#include "G4VisAttributes.hh"
#include "G4NistManager.hh"
#include "G4Sphere.hh"
#include "MyScoreParameterisation.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4SDManager.hh"
#include "G4PSFlatSurfaceCurrent.hh"
#include "MyPhotonSD.hh"
#include "G4Para.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include <cmath>

MyDetectorConstruction::MyDetectorConstruction()
  :fpWorldLogical(0)
  ,fpWorldPhysical(0)
{}

MyDetectorConstruction::~MyDetectorConstruction() {}

G4VPhysicalVolume* MyDetectorConstruction::Construct()
{
  // Material Definition
  DefineMaterials();  

  // Geometry Definition
  SetupGeometry();   

  // Return world volume
  return fpWorldPhysical;  
}

void MyDetectorConstruction::DefineMaterials()
{
  G4String symbol;             
  G4double a, z, density;     
  G4int ncomponents;
  G4double fractionmass;

  G4NistManager* man = G4NistManager::Instance();
  man->SetVerbose(1);
  
  // Define simple materials
  // Define beryllium, silicon and iron 

  new G4Material("Titanium",  z=22., a=47.90*g/mole,    density=4.540*g/cm3);
  new G4Material ("Beryllium", z=4., a=9.012182*g/mole, density=1.8480*g/cm3);
  new G4Material("Iron", z=26., a=55.845*g/mole, density=7.87*g/cm3);
  new G4Material("Silicon", z=14., a=28.0855*g/mole, density=2.33*g/cm3);

  // Define elements
  G4Element* N = new G4Element("Nitrogen", symbol="N", z=7., a=14.01*g/mole);
  G4Element* O = new G4Element("Oxygen",   symbol="O", z=8., a=16.00*g/mole);
  G4Element* C = man->FindOrBuildElement("C");
  G4Element* H = man->FindOrBuildElement("H");

  // Define air
  G4Material* air = new G4Material("Air", density= 1.290*mg/cm3, ncomponents=2);

  // Define polyvinyltoluene
  G4Material* my_polyvinyltoluene = new G4Material("Polyvinyltoluene", density=1.032*g/cm3, ncomponents=2);
  my_polyvinyltoluene->AddElement(H, fractionmass=0.5247148289);
  my_polyvinyltoluene->AddElement(C, fractionmass=1.-0.5247148289);
  air->AddElement(N, fractionmass=0.7);
  air->AddElement(O, fractionmass=0.3);
  G4cout<< my_polyvinyltoluene << G4endl;

  // Define vacuum
  G4Material* vacuum = new G4Material("Vacuum", density= 1.e-5*g/cm3, 
				      ncomponents=1, kStateGas, STP_Temperature, 
				      2.e-2*bar);
  
  vacuum->AddMaterial(air, fractionmass=1.);
  
  
  // Dump material information
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
  
  // optical properties of scintilator
  const G4int NUMENTRIES = 9;
  const G4int NUMENTRIES2 = 7;
  const G4double Decay_Time=2.1*ns; //fast time constant
  // const G4double Light_Output=64; // %

  //G4double wavelengths_sc[NUMENTRIES2] = {400*nm,420*nm,430*nm,440*nm,460*nm};
  //G4double Relative_Light_Output[NUMENTRIES2]={42, 84, 100, 84, 20};
  G4double energies[NUMENTRIES2] = {0.0000000001*eV,2.70*eV, 2.82*eV, 2.89*eV, 2.92*eV, 3.10*eV,100000000000000000.10*eV};
  G4double Light_Output_FAST[NUMENTRIES2]={0.002,0.125, 0.255, 0.303, 0.255, 0.058,0.003};
  G4double Scnt_SLOW_poly[NUMENTRIES2] = { 0.000010, 0.000020, 0.000030, 0.004000, 0.008000, 0.005000, 0.020000};
  G4double refractiveIndexBC[NUMENTRIES2]={1.58,1.58, 1.58, 1.58, 1.58, 1.58,1.58};
  G4double absorption[NUMENTRIES2]= {15. * m,15. * m,  15. * m,  15. * m,  15. * m,  15. * m,15. * m};

 // mirror
 G4Material * mirror = man->FindOrBuildMaterial("G4_Al");
 G4MaterialPropertiesTable* mirrorProperties = new G4MaterialPropertiesTable();
  G4double REFLECTIVITY[NUMENTRIES2] = { 
     0.99999999,0.99999999,0.99999999,0.99999999,0.99999999,0.99999999,0.99999999};
  G4double EFFICIENCY[NUMENTRIES2] = { 
     1.,1.,1.,1.,1.,1.,1.};
     G4double RINDEX_mirror[NUMENTRIES2]={0.00000001,0.000000001,0.000000001,0.000000001,0.0000000001,0.000000001,0.000000001};
    mirrorProperties->AddProperty("REFLECTIVITY", energies, EFFICIENCY, NUMENTRIES2);
    mirrorProperties->AddProperty("EFFICIENCY", energies, EFFICIENCY, NUMENTRIES2);
    mirrorProperties->AddProperty("RINDEX", energies, RINDEX_mirror, NUMENTRIES2);
    mirror->SetMaterialPropertiesTable(mirrorProperties);

  std::vector<G4double> refractiveIndexPS0 = { 1.59, 1.62, 1.63 };
std::vector<G4double> energyPS0    = { 2.00 *eV, 3.09 *eV, 3.47 *eV };
std::vector<G4double> AbsLengthPS0 = { 4.67 * m, 4.67 * m, 0.09 * mm };
 
std::vector<G4double> energy = {
    2.00*eV, 2.03*eV, 2.06*eV, 2.09*eV, 2.12*eV, 2.15*eV, 2.18*eV, 2.21*eV, 2.24*eV, 2.27*eV,
    2.30*eV, 2.33*eV, 2.36*eV, 2.39*eV, 2.42*eV, 2.45*eV, 2.48*eV, 2.51*eV, 2.54*eV, 2.57*eV,
    2.60*eV, 2.63*eV, 2.66*eV, 2.69*eV, 2.72*eV, 2.75*eV, 2.78*eV, 2.81*eV, 2.84*eV, 2.87*eV,
    2.90*eV, 2.93*eV, 2.96*eV, 2.99*eV, 3.02*eV, 3.05*eV, 3.08*eV, 3.11*eV, 3.14*eV, 3.17*eV,
    3.20*eV, 3.23*eV, 3.26*eV, 3.29*eV, 3.32*eV, 3.35*eV, 3.38*eV, 3.41*eV, 3.44*eV, 3.47*eV
  };
 
 std::vector<G4double> emissionFast = {
     0.05, 0.10, 0.30, 0.50,  0.75,  1.00,  1.26,  1.50,  1.81,  2.09,
     2.39, 2.63, 2.78, 2.98,  3.17,  3.49,  4.11,  4.70,  5.25,  5.79,
     6.42, 7.78, 8.69, 9.78, 10.61, 11.88, 13.68, 15.48, 16.56, 16.42,
     15., 12.54, 5.41, 2.13,  0.88,  0.40, 0.028, 0.13,  0.05,   0.02,
     0.01,  0.01,   0.01,    0,     0,     0,     0,     0,     0,  0
  };
 
// Add entries into properties table
  auto mptPolystyrene = new G4MaterialPropertiesTable();
  mptPolystyrene->AddProperty("RINDEX", energyPS0, refractiveIndexPS0);
  mptPolystyrene->AddProperty("ABSLENGTH", energyPS0, AbsLengthPS0);
  mptPolystyrene->AddProperty("SCINTILLATIONCOMPONENT1", energy, emissionFast,false,true);
  mptPolystyrene->AddConstProperty("SCINTILLATIONYIELD", 10. / keV); //10/keV
  mptPolystyrene->AddConstProperty("RESOLUTIONSCALE", 1.);
  mptPolystyrene->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 10. * ns);
  my_polyvinyltoluene->SetMaterialPropertiesTable(mptPolystyrene);
   // Set the Birk's Constant for the Polystyrene scintillator
  my_polyvinyltoluene->GetIonisation()->SetBirksConstant(0.126 * mm / MeV);
 
 
std::vector<G4double> energySmall = { 1.9 * eV, 6.7 * eV };
  std::vector<G4double> clad1RIndex = { 1.49, 1.49 };

std::vector<G4double> clad1AbsLen = {
    5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m,
    5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m,
    5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m,
    5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m, 5.40 * m,
    5.40 * m, 1.10 * m, 1.10 * m, 1.10 * m, 1.10 * m, 1.10 * m, 1.10 * m,
    1.10 * m, 1. * mm,  1. * mm,  1. * mm,  1. * mm,  1. * mm,  1. * mm,
    1. * mm,  1. * mm,  1. * mm,  1. * mm,  1. * mm,  1. * mm,  1. * mm,
    1. * mm
  };
 G4Material* pmma = man->FindOrBuildMaterial("G4_PLEXIGLASS");
// Add entries into properties table
  auto mptClad1 = new G4MaterialPropertiesTable();
  mptClad1->AddProperty("RINDEX", energySmall, clad1RIndex);
  mptClad1->AddProperty("ABSLENGTH", energy, clad1AbsLen);
 
  pmma->SetMaterialPropertiesTable(mptClad1);

 
std::vector<G4double> refractiveIndex = { 1.0, 1.0 };
 G4Material* fAir = man->FindOrBuildMaterial("Air");
  auto mpt = new G4MaterialPropertiesTable();
  mpt->AddProperty("RINDEX", energySmall, refractiveIndex);
 
  fAir->SetMaterialPropertiesTable(mpt);
}
void MyDetectorConstruction::SetupScoring(G4LogicalVolume* scoringVolume)
{
 // Create a new sensitive detector named "MyDetector"
   G4MultiFunctionalDetector* detector =
    new G4MultiFunctionalDetector("MyDetector");

  // Get pointer to detector manager
  G4SDManager* manager = G4SDManager::GetSDMpointer();  

  // Register detector with manager
  manager->AddNewDetector(detector);

  // Attach detector to scoring volume
  //scoringVolume->SetSensitiveDetector(detector);

  // Create a primitive Scorer named MyScorer
  G4PSFlatSurfaceCurrent* scorer =
    new G4PSFlatSurfaceCurrent("MyScorer",fCurrent_In);
// Create a new sensitive detector named "Photon"
  MyPhotonSD* photonSD =  new MyPhotonSD("Photon");
  manager->AddNewDetector(photonSD);
  scoringVolume->SetSensitiveDetector(photonSD);
  // Register scorer with detector
  detector->RegisterPrimitive(scorer);
}
void MyDetectorConstruction::SetupGeometry()
{
  // NIST definition of air
  G4NistManager* man = G4NistManager::Instance();
  man->SetVerbose(1);

  // definitions


  G4Material* air = G4Material::GetMaterial("Air");
  G4cout<< "air material properties table:" << G4endl;
  //air->GetMaterialPropertiesTable()->DumpTable();

  G4Material* PbWO4 = man->FindOrBuildMaterial("G4_PbWO4");
  G4Material* mirrorMaterial = man->FindOrBuildMaterial("G4_Al");
  G4Material* stykMaterial = air;
  G4Material* pmma = man->FindOrBuildMaterial("G4_PLEXIGLASS");
  G4Material* plastic_poly = G4Material::GetMaterial("Polyvinyltoluene");
G4cout<< "polyvinyltoluene materials table:" << G4endl;
  plastic_poly->GetMaterialPropertiesTable()->DumpTable();


  // World volume
  G4Box* worldSolid = new G4Box("World_Solid",           // Name
				0.5*m, 0.5*m, 0.5*m);    // Half lengths
  
  fpWorldLogical = new G4LogicalVolume(worldSolid,	 // Solid
				       air,	         // Material
				       "World_Logical"); // Name
  
  fpWorldPhysical = new G4PVPlacement(0,	         // Rotation matrix pointer
				      G4ThreeVector(),   // Translation vector
				      fpWorldLogical,	 // Logical volume
				      "World_Physical",	 // Name
				      0,		 // Mother volume
				      false,		 // Unused boolean parameter
				      0);		 // Copy number


  ////////////////////////////////////////////////////////////////////////
  // Scintillator
  G4Box* ScintillatorSolid = new G4Box("Scintillator_Solid",           // Name
				8.0*mm/2.0, 8.0*cm/2.0, 5.0*cm/2.0);    // Half lengths
  
  G4LogicalVolume* fpScincillatorLogical = new G4LogicalVolume(ScintillatorSolid,	 // Solid
				       plastic_poly,	         // Material
				       "Scintillator_Logical"); // Name
  
  G4VPhysicalVolume* fpScintillatorPhysical=new G4PVPlacement(0,	         // Rotation matrix pointer
				      G4ThreeVector(0.,0.,-2.5*cm),   // Translation vector
				      fpScincillatorLogical,	 // Logical volume
				      "Scintillator_Physical",	 // Name
				      fpWorldLogical,		 // Mother volume
				      false,		 // Unused boolean parameter
				      0);		 // Copy number
  
  ////////////////////////////////////////////////////////////////////////
  // OpticalFiber
G4double theta=std::atan(8.0/125.0);
G4double theta_deg = theta*180.0/M_PI;   // Kst theta (pochylenie osi trapezoidu)
G4Trap* OpticalFiberSolid = new G4Trap("Optical_Fiber_Solid",
  12.5 * cm / 2.0,        // Polowa wysokosci w z
   theta_deg*deg,     // Theta
   0*deg,       // Phi
   4 * cm / 2.0,        // Wysokosc dolnej podstawy
   8.0 * mm / 2.0,       // Szerokosc dolnej podstawy (x1)
   8.0 * mm / 2.0,
   0*deg,
   6.0 * mm / 2.0,       // Szerokosc dolnej podstawy (x2)
   8.0 * mm / 2.0,    // Kat dolnej podstawy
   8.0 * mm / 2.0,        // Wysokosc gornej podstawy       // Szerokosc gornej podstawy (x1)       // Szerokosc grnej podstawy (x2)
   0*deg);


  G4LogicalVolume* OpticalFiberLogical = new G4LogicalVolume(
        OpticalFiberSolid,          // Ksztalt geometryczny
        pmma,       // Material
        "Optical_Fiber_Logical"
        ); // Name
  
  G4VPhysicalVolume* OpticalFiberPhysical1=new G4PVPlacement(0,	         // Rotation matrix pointer
				      G4ThreeVector(8.0*mm/2.0,4*cm/2.0,12.5*cm/2.),   // Translation vector
				      OpticalFiberLogical,	 // Logical volume
				      "Optical_Fiber_Physical",	 // Name
				      fpWorldLogical,		 // Mother volume
				      false,		 // Unused boolean parameter
				      0);		 // Copy number
  ////////////////////////////////////////////////////////////////////////
  // Second optic fiber
  G4VPhysicalVolume* OpticalFiberPhysical2=new G4PVPlacement(0,	         // Rotation matrix pointer
				      G4ThreeVector(8.0*mm/2.0,-4*cm/2.0,12.5*cm/2.),   // Translation vector
				      OpticalFiberLogical,	 // Logical volume
				      "Optical_Fiber_Physical",	 // Name
				      fpWorldLogical,		 // Mother volume
				      false,		 // Unused boolean parameter
				      0);		 // Copy number
 ////////////////////////////////////////////////////////////////

  //// Mirrors - optical fiber

  G4Trap* OpticalFiberSolid_mirror = new G4Trap("Optical_Fiber_Solid",
  12.5 * cm / 2.0,      
   theta_deg*deg,    
   0*deg,      
   4 * cm / 2.0,       
   8.0 * nm / 2.0,     
   8.0 * nm / 2.0,
   0*deg,
   6.0 * mm / 2.0,     
   8.0 * nm / 2.0,   
   8.0 * nm / 2.0,           
   0*deg);


  G4LogicalVolume* OpticalFiberLogical_mirror = new G4LogicalVolume(
        OpticalFiberSolid_mirror,         
        mirrorMaterial,    
        "Optical_Fiber_Logical_mirror"
        ); // Name   

  G4VPhysicalVolume* OpticalFiberPhysical_mirror1=new G4PVPlacement(0,	         // Rotation matrix pointer
				      G4ThreeVector(8.0*mm/2.0+8.0/2*mm+8.0/2.0*nm,4*cm/2.0,12.5*cm/2.),   // Translation vector
				      OpticalFiberLogical_mirror,	 // Logical volume
				      "Optical_Fiber_Physical_mirror",	 // Name
				      fpWorldLogical,		 // Mother volume
				      false,		 // Unused boolean parameter
				      0);		 // Copy number
  G4VPhysicalVolume* OpticalFiberPhysical_mirror2=new G4PVPlacement(0,	         // Rotation matrix pointer
				      G4ThreeVector(8.0*mm/2.0-8.0/2*mm-8.0/2.0*nm,4*cm/2.0,12.5*cm/2.),   // Translation vector
				      OpticalFiberLogical_mirror,	 // Logical volume
				      "Optical_Fiber_Physical_mirror",	 // Name
				      fpWorldLogical,		 // Mother volume
				      false,		 // Unused boolean parameter
				      1);		 // Copy number
  G4VPhysicalVolume* OpticalFiberPhysical_mirror3=new G4PVPlacement(0,	         // Rotation matrix pointer
				      G4ThreeVector(8.0*mm/2.0+8.0/2*mm+8.0/2.0*nm,-4*cm/2.0,12.5*cm/2.),   // Translation vector
				      OpticalFiberLogical_mirror,	 // Logical volume
				      "Optical_Fiber_Physical_mirror",	 // Name
				      fpWorldLogical,		 // Mother volume
				      false,		 // Unused boolean parameter
				      2);		 // Copy number
              
  G4VPhysicalVolume* OpticalFiberPhysical_mirror4=new G4PVPlacement(0,	         // Rotation matrix pointer
				      G4ThreeVector(8.0*mm/2.0-8.0/2*mm-8.0/2.0*nm,-4*cm/2.0,12.5*cm/2.),   // Translation vector
				      OpticalFiberLogical_mirror,	 // Logical volume
				      "Optical_Fiber_Physical_mirror",	 // Name
				      fpWorldLogical,		 // Mother volume
				      false,		 // Unused boolean parameter
				      3);		 // Copy number

  G4LogicalVolume* OpticalFiberLogical_mirror2 = new G4LogicalVolume(
        OpticalFiberSolid,   
        mirrorMaterial,      
        "Optical_Fiber_Logical_mirror2"
        ); // Name   

G4RotationMatrix* rotation = new G4RotationMatrix();
rotation->rotateY(180.0 * deg); 
  G4double theta_mirr_rad=std::atan(1.7/12.50);
  G4double theta_mirr_deg = theta_mirr_rad * (180.0 / M_PI);
           
  G4VPhysicalVolume* OpticalFiberPhysical_mirror5=new G4PVPlacement(rotation,	         // Rotation matrix pointer
				      G4ThreeVector(8.0*mm/2.0,8*cm/2+6.0*mm/2,12.5*cm/2.),   // Translation vector
				      OpticalFiberLogical_mirror2,	 // Logical volume
				      "Optical_Fiber_Physical_mirror5",	 // Name
				      fpWorldLogical,		 // Mother volume
				      false,		 // Unused boolean parameter
				      0);		 // Copy number
  G4VPhysicalVolume* OpticalFiberPhysical_mirror6=new G4PVPlacement(rotation,	         // Rotation matrix pointer
				      G4ThreeVector(8.0*mm/2.0,-8*cm/2-6.0*mm/2,12.5*cm/2.),   // Translation vector
				      OpticalFiberLogical_mirror2,	 // Logical volume
				      "Optical_Fiber_Physical_mirror6",	 // Name
				      fpWorldLogical,		 // Mother volume
				      false,		 // Unused boolean parameter
				      1);		 // Copy number


  G4Trap* OpticalFiber_mirror_between = new G4Trap("Optical_Fiber_mirror_between",
  12.5 * cm / 2.0,      
   theta_deg*deg,     
   0*deg,     
   3.4 * cm / 2.0,     
   8.0 * mm / 2.0,    
   8.0 * mm / 2.0,
   0*deg,
   0.0000001 * nm / 2.0,      
   8.0 * mm / 2.0,   
   8.0 * mm / 2.0,         
   0*deg);


  G4LogicalVolume* OpticalFiberLogical_between = new G4LogicalVolume(
        OpticalFiber_mirror_between,       
        mirrorMaterial,     
        "Optical_Fiber_Logical_mirror"
        ); // Name  

  G4VPhysicalVolume* OpticalFiberPhysical_mirror7=new G4PVPlacement(rotation,	         // Rotation matrix pointer
				      G4ThreeVector(8.0*mm/2.0,0,12.5*cm/2+1*nm),   // Translation vector
				      OpticalFiberLogical_between,	 // Logical volume
				      "Optical_Fiber_Physical_mirror7",	 // Name
				      fpWorldLogical,		 // Mother volume
				      false,		 // Unused boolean parameter
				      0);		 // Copy number
  
  G4OpticalSurface* opticalMirror = new G4OpticalSurface("MirrorSurface");
  opticalMirror->SetType(dielectric_metal);
  opticalMirror->SetFinish(polished); 
  opticalMirror->SetModel(unified); 
  opticalMirror->SetMaterialPropertiesTable(mirrorMaterial->GetMaterialPropertiesTable());

  G4OpticalSurface* opticalstyk = new G4OpticalSurface("stykSurface");
  opticalMirror->SetType(dielectric_dielectric);
  opticalMirror->SetFinish(ground); 
  opticalMirror->SetModel(unified); 
  opticalstyk->SetMaterialPropertiesTable(stykMaterial->GetMaterialPropertiesTable());

  
  G4Box* mirror = new G4Box("Mirror_Solid",           // Name
				2.0*nm/2.0, 8.0*cm/2.0, 5.0*cm/2.0);

  G4LogicalVolume* mirrorLogical = new G4LogicalVolume(
        mirror,         
        mirrorMaterial,     
        "mirror_Logical"
        );
  G4VPhysicalVolume* fpMirror1Physical = new G4PVPlacement(0,	         // Rotation matrix pointer
				      G4ThreeVector(-4.0*mm-2.0*nm/2.0,0.*cm,-2.5*cm),   // Translation vector
				      mirrorLogical,	 // Logical volume
				      "mirror_Physical",	 // Name
				      fpWorldLogical,		 // Mother volume
				      false,		 // Unused boolean parameter
				      0);		 // Copy number
              
  G4VPhysicalVolume* fpMirror2Physical =new G4PVPlacement(0,	         // Rotation matrix pointer
				      G4ThreeVector(4.0*mm+2.0*nm/2.0,0.*cm,-2.5*cm),   // Translation vector
				      mirrorLogical,	 // Logical volume
				      "mirror_Physical",	 // Name
				      fpWorldLogical,		 // Mother volume
				      false,		 // Unused boolean parameter
				      1);		 // Copy number
              
  G4Box* mirror_shorter = new G4Box("Mirror_shorter_Solid",           // Name
				8.0*mm/2.0, 8.0*cm/2.0, 4.0*nm/2.0);

  G4LogicalVolume* mirrorShorterLogical = new G4LogicalVolume(
        mirror_shorter,        
        mirrorMaterial,   
        "mirror_shorter_Logical"
        );
  G4VPhysicalVolume* fpMirrorShorterPhysical = new G4PVPlacement(0,	         // Rotation matrix pointer
				      G4ThreeVector(0.*cm,0.*cm,-5*cm - 4.0/2.0*nm),   // Translation vector
				      mirrorShorterLogical,	 // Logical volume
				      "mirror_Shorter_Physical",	 // Name
				      fpWorldLogical,		 // Mother volume
				      false,		 // Unused boolean parameter
				      0);		 // Copy number


  G4Box* mirror_top = new G4Box("Mirror_top_Solid",           // Name
				8.0*mm/2.0, 1.0*nm/2.0, 5.0*cm/2.0);

  G4LogicalVolume* mirrorTopLogical = new G4LogicalVolume(
        mirror_top,        
        mirrorMaterial, 
        "mirror_shorter_Logical"
        );
  G4VPhysicalVolume* fpMirrorTopPhysical = new G4PVPlacement(0,	         // Rotation matrix pointer
				      G4ThreeVector(0.*cm,4.*cm+1.0/2.0*nm,-2.5*cm),   // Translation vector
				      mirrorTopLogical,	 // Logical volume
				      "mirror_Top_Physical",	 // Name
				      fpWorldLogical,		 // Mother volume
				      false,		 // Unused boolean parameter
				      0);		 // Copy number
  G4VPhysicalVolume* fpMirrorTopPhysical2 = new G4PVPlacement(0,	         // Rotation matrix pointer
				      G4ThreeVector(0.*cm,-4*cm-1.0/2.0*nm,-2.5*cm),   // Translation vector
				      mirrorTopLogical,	 // Logical volume
				      "mirror_Top_Physical",	 // Name
				      fpWorldLogical,		 // Mother volume
				      false,		 // Unused boolean parameter
				      1);		 // Copy number

  G4Box* SD = new G4Box("SD_Solid",           // Name
				8.0*mm/2.0, 6.0*mm/2.0, 0.000001*cm/2.0);

  G4LogicalVolume* SDLogical = new G4LogicalVolume(
        SD,       
        pmma,     
        "SD_Logical"
        );
  SetupScoring(SDLogical); 
  
  G4VPhysicalVolume* SDPhysical = new G4PVPlacement(0,	         // Rotation matrix pointer
				      G4ThreeVector(0.8*cm,2.*cm,12.500001*cm),   // Translation vector
				      SDLogical,	 // Logical volume
				      "SD_Physical",	 // Name
				      fpWorldLogical,		 // Mother volume
				      false,		 // Unused boolean parameter
				      0);		 // Copy number  
             
  G4VPhysicalVolume* SDPhysical2 = new G4PVPlacement(0,	         // Rotation matrix pointer
				      G4ThreeVector(0.8*cm,-2.*cm,12.500001*cm),   // Translation vector
				      SDLogical,	 // Logical volume
				      "SD_Physical",	 // Name
				      fpWorldLogical,		 // Mother volume
				      false,		 // Unused boolean parameter
				      1);	
                       
// sphere
  G4double innerRadius = 40.0*cm;
  G4double outerRadius = 40.01*cm;

  G4VSolid* sphereSolid =
    new G4Sphere("Sphere_Solid",     //its name
                 innerRadius,                     //its inner radius
                 outerRadius,                     //its outer radius
                 0.0*deg,                           //its starting phi angle
                 360.0*deg,                       //its total phi angle
                 0.0*deg,                           //its starting theta angle
                 90.0*deg);                       //its total theta angle

  G4LogicalVolume* sphereLogical =
    new G4LogicalVolume(sphereSolid,              //its solid
                        air,                                              //its material
                        "Sphere_logical");                      //its name

  new G4PVPlacement(0,                         //no rotation
                    G4ThreeVector(),                 //no transition
                    sphereLogical,                     //its logical volume
                    "Sphere_physical",               //its name
                    fpWorldLogical,                   //its mother volume
                    false,                                    //true not implemented
                    0);                                        // copy number                            
////////////////////////////////////////////////////////////////////////
// score
G4int nRings=15;
G4double deltaTheta=90.0*deg/nRings;
                
  G4VSolid * scoreSolid =
     new G4Sphere("Score_Solid",          //its name
                 innerRadius,            //its inner radius
                 outerRadius,            //its outer radius
                 0.0*deg,                   //its starting phi angle
                 360.0*deg,               //its total phi angle
                 0.0*deg,                   //its starting theta angle
                 deltaTheta);             //its total theta angle

  G4LogicalVolume* scoreLogical =
    new G4LogicalVolume(scoreSolid,      //its solid
                        air,                                    //its material 
                        "Score_Logical");             //its name
                   
  MyScoreParameterisation* param =
    new MyScoreParameterisation(innerRadius,outerRadius,nRings);

  new G4PVParameterised("Score_Physical",   //its name
                        scoreLogical,                             //its logical volume
               sphereLogical,                                    //its mother logical
                       kZAxis,                                        //axis
                       nRings,                                        //number of copies
                       param);                                       //parameterisation
////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  // Visualisation attributes
  
  // Invisible world volume.
  fpWorldLogical->SetVisAttributes(G4VisAttributes::GetInvisible());

  // Sphere - white
  G4VisAttributes* sphereAttributes = new G4VisAttributes(G4Colour(1.0,1.0,1.0,0.5));
  sphereAttributes->SetVisibility(false);
  sphereLogical->SetVisAttributes(sphereAttributes);

  // Score - green
  G4VisAttributes* scoreAttributes = new G4VisAttributes(G4Colour(0.0,1.0,0.0,0.5));
  scoreAttributes->SetVisibility(false);
  scoreLogical->SetVisAttributes(scoreAttributes);
  // scintillator - purple
  G4VisAttributes* scintillatorAttributes = new G4VisAttributes(G4Colour(G4Colour(0.745, 0.431, 0.961, 0.5)));
  scintillatorAttributes->SetForceSolid(true);
  fpScincillatorLogical->SetVisAttributes(scintillatorAttributes);
  // optical fibers - pink
  G4VisAttributes* opticalAttributes = new G4VisAttributes(G4Colour(0.886, 0.659, 0.945, 0.5));
  opticalAttributes->SetForceSolid(true);
  OpticalFiberLogical->SetVisAttributes(opticalAttributes);
  // mirrors of optical fiber - invisible
  G4VisAttributes* opticalAttributes_fibers = new G4VisAttributes();
  opticalAttributes_fibers->SetForceSolid(false);
  opticalAttributes_fibers->SetVisibility(false);
  OpticalFiberLogical_mirror2->SetVisAttributes(opticalAttributes_fibers);
  OpticalFiberLogical_between->SetVisAttributes(opticalAttributes_fibers);
  
}



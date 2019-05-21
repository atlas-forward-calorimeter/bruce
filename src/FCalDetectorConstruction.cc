//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
//This is the source code of the detctor constuction                                  |
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
//Author: Zhaoyuan Cui                                                                |
//        2016 Summer                                                                 |
//        Physics department, The University of Arizona                               |
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
//This code generates the geometry of FCal1, including liguid Argon gaps and rods.    |
//Arrangement of the rods is accomplished by following the actual region arrangement. |
//Sensetive detector should be provided by the user for his or her own purpose.       |
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|

//particle distance on z: 4750mm

#include "FCalDetectorConstruction.hh"
#include "FCalEmCalorimeterSD.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Polyhedra.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4PVReplica.hh"

#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
#include "G4AssemblyVolume.hh"
#include "G4Material.hh"
#include "G4SubtractionSolid.hh"
#include "G4UserLimits.hh"



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

FCalDetectorConstruction::FCalDetectorConstruction()
  : //G4VUserDetectorConstruction(),
    //gapLogical(0),
    fpWorldLogical(0),
    fpWorldPhysical(0)
{}

FCalDetectorConstruction::~FCalDetectorConstruction()
{

}
/*
void FCalDetectorConstruction::ConstructSDandField()
{
  G4SDManager *SDman=G4SDManager::GetSDMpointer();
  G4String SDname;

  G4VSensitiveDetector* emCal
    =new FCalEmCalorimeterSD(SDname="/EMCalorimeter",
  			     "CalorimeterHitsCollection");
  SDman->AddNewDetector(emCal);
  gapLogical->SetSensitiveDetector(emCal);
}
*/
//new SD method start
void FCalDetectorConstruction::ConstructSDandField()
{
  // Sensitive detectors

  G4String trackerChamberSDname = "FCalSD";
  FCalEmCalorimeterSD* aTrackerSD = new FCalEmCalorimeterSD(trackerChamberSDname,
                                            "FCalHitsCollection");
  G4SDManager::GetSDMpointer()->AddNewDetector(aTrackerSD);
  // Setting aTrackerSD to all logical volumes with the same name 
  // of "Chamber_LV".
  SetSensitiveDetector("Gap_Logical", aTrackerSD, true);

  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  /*
  G4ThreeVector fieldValue = G4ThreeVector();
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);
  
  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
  */
}

//new SD method end

void FCalDetectorConstruction::ConstructMaterials()
{
  G4String symbol;             
  G4double a, z, n, density;     
  G4int ncomponents,nisotopes, natoms;
  G4double fractionmass;

  new G4Material("Copper",z=29., a=63.546*CLHEP::g/CLHEP::mole, density=8.96*CLHEP::g/CLHEP::cm3);
  new G4Material("LAr",   z=18., a=39.948*CLHEP::g/CLHEP::mole, density=1.396*CLHEP::g/CLHEP::cm3);
  new G4Material("Tungsten",z=74.,a=183.84*CLHEP::g/CLHEP::mole,density=19.25*CLHEP::g/CLHEP::cm3);
  new G4Material("Aluminum",z=13.,a=26.981*CLHEP::g/CLHEP::mole,density=2.7*CLHEP::g/CLHEP::cm3);
  new G4Material("Lead",   z=82.,a=207.2*CLHEP::g/CLHEP::mole, density=11.34*CLHEP::g/CLHEP::cm3);  //density original 11.34
    
  // Define elements
  G4Element* N  = new G4Element("Nitrogen", symbol="N",  z=7.  , a=14.01*CLHEP::g/CLHEP::mole);
  G4Element* O  = new G4Element("Oxygen",   symbol="O",  z=8.  , a=16.00*CLHEP::g/CLHEP::mole);
  G4Element* Ar = new G4Element("Argon",    symbol="Ar", z=18. , a=39.948*CLHEP::g/CLHEP::mole);
  G4Element* Cu = new G4Element("Copper",   symbol="Cu", z=29. , a=63.546*CLHEP::g/CLHEP::mole);
  G4Element* W = new G4Element("Tungsten",  symbol="W",  z=74. , a=183.84*CLHEP::g/CLHEP::mole);//shouldn't need this
  G4Element* H = new G4Element("Hydrogen",  symbol="H",  z=1.  , a=1.0078*CLHEP::g/CLHEP::mole);
  G4Element* C = new G4Element("Carbon",  symbol="C",  z=6.  , a=12*CLHEP::g/CLHEP::mole);
  

  // Define air
  G4Material* air = new G4Material("Air", density= 1.290*CLHEP::mg/CLHEP::cm3, ncomponents=2);
  air->AddElement(N, fractionmass=0.7);
  air->AddElement(O, fractionmass=0.3);

  // Define vacuum
  G4Material* vacuum = new G4Material("Vacuum", density= 1.e-5*CLHEP::g/CLHEP::cm3, 
				      ncomponents=1, kStateGas, CLHEP::STP_Temperature, 
				      2.e-2*CLHEP::bar);
  
  vacuum->AddMaterial(air, fractionmass=1.);

  // Define PEEK
  G4Material* PEEK = new G4Material("PEEK", density= 1.32*CLHEP::g/CLHEP::cm3, ncomponents=3);
  PEEK->AddElement(C, natoms=19);
  PEEK->AddElement(H, natoms=12);
  PEEK->AddElement(O, natoms=3);


  // Dump material information
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;


  
}

G4VPhysicalVolume *FCalDetectorConstruction::Construct()
{
  //All materials
  ConstructMaterials();
  //Geometry Definition
  
  SetupGeometry();   




    
  return fpWorldPhysical;
}

void FCalDetectorConstruction::SetupGeometry()
{
  G4Material* air    = G4Material::GetMaterial("Air");
  G4Material* vacuum = G4Material::GetMaterial("Vacuum");
  G4Material* copper = G4Material::GetMaterial("Aluminum");  // copper is aluminum
  G4Material* Rcopper = G4Material::GetMaterial("Copper");  //Rcopper is copper

  G4Material* lar    = G4Material::GetMaterial("LAr");
  G4Material* tungsten=G4Material::GetMaterial("Tungsten");
  G4Material* PEEK   = G4Material::GetMaterial("PEEK");
  G4Material* lead   = G4Material::GetMaterial("Lead");
  
  // World volume
  G4Box* worldSolid = new G4Box("World_Solid",           // Name
				20*CLHEP::cm, 20*CLHEP::cm, 20*CLHEP::cm);    // Half lengths
  
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
  // HandsOn3: Scoring volumes

  // HandsOn3: Construct scoring components


  ////////////////////////////////////////////////////////////////////////
  // Shaft (the inner copper part of the source)
  
  G4Tubs* shaftTubs = new G4Tubs("Shaft", // Name
				  0.*CLHEP::cm,          // Inner radius
				  0.15*CLHEP::cm,        // Outer radius
				  0.55*CLHEP::cm,        // Half length in z
				  0.*CLHEP::deg,         // Starting phi angle         
				  360.*CLHEP::deg);      // Segment angle
  
  G4LogicalVolume* shaftLogical =
    new G4LogicalVolume(shaftTubs,        // Solid
			copper,           // Material
			"Shaft_Logical"); // Name
  //left Shaft
  G4Tubs* shaftTubsL = new G4Tubs("ShaftL", // Name
				  0.*CLHEP::cm,          // Inner radius
				  0.41/2.*CLHEP::cm,        // Outer radius
				  0.6*CLHEP::cm,        // Half length in z
				  0.*CLHEP::deg,         // Starting phi angle         
				  360.*CLHEP::deg);      // Segment angle
  
  G4LogicalVolume* shaftLLogical =
    new G4LogicalVolume(shaftTubsL,        // Solid
			copper,           // Material
			"ShaftL_Logical"); // Name
  //right Shaft
  G4Tubs* shaftTubsR = new G4Tubs("ShaftR", // Name
				  0.*CLHEP::cm,          // Inner radius
				  0.41/2.*CLHEP::cm,        // Outer radius
				  0.6*CLHEP::cm,        // Half length in z
				  0.*CLHEP::deg,         // Starting phi angle         
				  360.*CLHEP::deg);      // Segment angle
  
  G4LogicalVolume* shaftRLogical =
    new G4LogicalVolume(shaftTubsR,        // Solid
			copper,           // Material
			"ShaftR_Logical"); // Name
  /*
  new G4PVPlacement(0,                                // Rotation matrix pointer
		    G4ThreeVector(0.,0., 0.),         // Translation vector
		    shaftLogical,                     // Logical volume
		    "Shaft_Physical",                 // Name
		    fpWorldLogical,                   // Mother volume
		    false,                            // Unused boolean 
		    0);                               // Copy number
  */
  ////////////////////////////////////////////////////////////////////////
  // Cavity (the space between shaft and copper foil)
  
  G4Tubs* cavityTubsOut = new G4Tubs("CavityOut", // Name
				  0.15*CLHEP::cm,        // Inner radius
				  0.205*CLHEP::cm,      // Outer radius  //changed from 0.2026 to 0.205
				  0.55*CLHEP::cm,        // Half length in z
				  0.*CLHEP::deg,         // Starting phi angle         
				  360.*CLHEP::deg);      // Segment angle
  G4Tubs* foilTubs = new G4Tubs("Foil", // Name
				  0.199*CLHEP::cm,      // Inner radius //changed from 0.2026 to 0.199
				  0.205*CLHEP::cm,      // Outer radius  //changed from 0.2086 to 0.205
				  0.5*CLHEP::cm,        // Half length in z
				  0.*CLHEP::deg,         // Starting phi angle         
				  360.*CLHEP::deg);      // Segment angle
  G4SubtractionSolid *cavityTubs = new G4SubtractionSolid("Cavity",
							cavityTubsOut,
							  foilTubs);

  
  G4LogicalVolume* cavityLogical =
    new G4LogicalVolume(cavityTubs,        // Solid
			vacuum,            // Material
			"Cavity_Logical"); // Name
  /*
  new G4PVPlacement(0,                                // Rotation matrix pointer
		    G4ThreeVector(0.,0., 0.),         // Translation vector
		    cavityLogical,                    // Logical volume
		    "Cavity_Physical",                // Name
		    fpWorldLogical,                   // Mother volume
		    false,                            // Unused boolean 
		    0);                               // Copy number
  */
  ////////////////////////////////////////////////////////////////////////
  // Foil (the copper foil on which the Sr90 is placed)
  

  
  G4LogicalVolume* foilLogical =
    new G4LogicalVolume(foilTubs,        // Solid
			copper,          // Material
			"Foil_Logical"); // Name
  /*
  new G4PVPlacement(0,                                // Rotation matrix pointer
		    G4ThreeVector(0.,0., 0.),         // Translation vector
		    foilLogical,                    // Logical volume
		    "Foil_Physical",                // Name
		    fpWorldLogical,                   // Mother volume
		    false,                            // Unused boolean 
		    0);                               // Copy number
  */
  ////////////////////////////////////////////////////////////////////////
  // Space (the space between copper foil and wall)
  
  G4Tubs* spaceTubs = new G4Tubs("Space", // Name
				  0.2086*CLHEP::cm,      // Inner radius
				  0.2106*CLHEP::cm,      // Outer radius
				  0.70*CLHEP::cm,        // Half length in z
				  0.*CLHEP::deg,         // Starting phi angle         
				  360.*CLHEP::deg);      // Segment angle
  
  G4LogicalVolume* spaceLogical =
    new G4LogicalVolume(spaceTubs,        // Solid
			vacuum,           // Material
			"Space_Logical"); // Name
  /*
  new G4PVPlacement(0,                                // Rotation matrix pointer
		    G4ThreeVector(0.,0., 0.),         // Translation vector
		    spaceLogical,                     // Logical volume
		    "Space_Physical",                 // Name
		    fpWorldLogical,                   // Mother volume
		    false,                            // Unused boolean 
		    0);                               // Copy number
  */
  ////////////////////////////////////////////////////////////////////////
  // Wall (the cavity wall)
  
  G4Tubs* wallTubs = new G4Tubs("Wall", // Name
				  0.205*CLHEP::cm,      // Inner radius changed from 0.2106 to 0.205
				  0.2356*CLHEP::cm,      // Outer radius
				  1.75*CLHEP::cm,        // Half length in z
				  0.*CLHEP::deg,         // Starting phi angle         
				  360.*CLHEP::deg);      // Segment angle
  
  G4LogicalVolume* wallLogical =
    new G4LogicalVolume(wallTubs,        // Solid
			copper,          // Material
			"Wall_Logical"); // Name
  /*
  new G4PVPlacement(0,                                // Rotation matrix pointer
		    G4ThreeVector(0.,0., 0.),         // Translation vector
		    wallLogical,                      // Logical volume
		    "Wall_Physical",                  // Name
		    fpWorldLogical,                   // Mother volume
		    false,                            // Unused boolean 
		    0);                               // Copy number
  */
  ////////////////////////////////////////////////////////////////////////
  // Gap (the LAr gap)
  
  G4Tubs* gapTubs = new G4Tubs("Gap", // Name
			       0.2356*CLHEP::cm,      // Inner radius
			       0.2625*CLHEP::cm,      // Outer radius
			       3.5/2.*CLHEP::cm,         // Half length in z
			       0.*CLHEP::deg,         // Starting phi angle         
			       360.*CLHEP::deg);      // Segment angle
  
  G4LogicalVolume* gapLogical =
    new G4LogicalVolume(gapTubs,        // Solid
			lar,            // Material
			"Gap_Logical"); // Name
  /*
  new G4PVPlacement(0,                                // Rotation matrix pointer
		    G4ThreeVector(0.,0., 0.),         // Translation vector
		    gapLogical,                       // Logical volume
		    "Gap_Physical",                   // Name
		    fpWorldLogical,                   // Mother volume
		    false,                            // Unused boolean 
		    0);                               // Copy number
  */
  ////////////////////////////////////////////////////////////////////////
  // Tube (the FCal tube)
  
  G4Tubs* tubeTubs = new G4Tubs("Tube", // Name
			       0.2625*CLHEP::cm,      // Inner radius
			       0.2875*CLHEP::cm,      // Outer radius  changed from 0.2878 to 0.2875
			       0.4*CLHEP::cm,         // Half length in z changed from 1.0 to 0.4
			       0.*CLHEP::deg,         // Starting phi angle         
			       360.*CLHEP::deg);      // Segment angle
  
  G4LogicalVolume* tubeLogical =
    new G4LogicalVolume(tubeTubs,        // Solid
			copper,          // Material
			"Tube_Logical"); // Name
  //Tube left
  G4Tubs* tubeTubsL = new G4Tubs("TubeL", // Name
			       0.2625*CLHEP::cm,      // Inner radius
			       0.2875*CLHEP::cm,      // Outer radius  changed from 0.2878 to 0.2875
			       1.35/2.*CLHEP::cm,         // Half length in z changed from 1.0 to 4.0
			       0.*CLHEP::deg,         // Starting phi angle         
			       360.*CLHEP::deg);      // Segment angle
  
  G4LogicalVolume* tubeLLogical =
    new G4LogicalVolume(tubeTubsL,        // Solid
			copper,          // Material
			"TubeL_Logical"); // Name
  //Tube right
  G4Tubs* tubeTubsR = new G4Tubs("TubeR", // Name
			       0.2625*CLHEP::cm,      // Inner radius
			       0.2875*CLHEP::cm,      // Outer radius  changed from 0.2878 to 0.2875
			       1.35/2.*CLHEP::cm,         // Half length in z changed from 1.0 to 4.0
			       0.*CLHEP::deg,         // Starting phi angle         
			       360.*CLHEP::deg);      // Segment angle
  
  G4LogicalVolume* tubeRLogical =
    new G4LogicalVolume(tubeTubsR,        // Solid
			copper,          // Material
			"TubeR_Logical"); // Name



  /*
  new G4PVPlacement(0,                                // Rotation matrix pointer
		    G4ThreeVector(0.,0., 0.),         // Translation vector
		    tubeLogical,                      // Logical volume
		    "Tube_Physical",                  // Name
		    fpWorldLogical,                   // Mother volume
		    false,                            // Unused boolean 
		    0);                               // Copy number
  */
  ////////////////////////////////////////////////////////////////////////
  // PlugL (the left end)
  
  G4Tubs* plugLTubs = new G4Tubs("PlugL", // Name
			       0.0*CLHEP::cm,         // Inner radius
			       0.2356*CLHEP::cm,      // Outer radius
			       0.3*CLHEP::cm/2,       // Half length in z
			       0.*CLHEP::deg,         // Starting phi angle         
			       360.*CLHEP::deg);      // Segment angle
  
  G4LogicalVolume* plugLLogical =
    new G4LogicalVolume(plugLTubs,        // Solid
			copper,           // Material
			"PlugL_Logical"); // Name
  /*
  new G4PVPlacement(0,                                     // Rotation matrix pointer
		    G4ThreeVector(0.,0., -0.85*CLHEP::cm), // Translation vector
		    plugLLogical,                          // Logical volume
		    "PlugL_Physical",                      // Name
		    fpWorldLogical,                        // Mother volume
		    false,                                 // Unused boolean 
		    0);                               // Copy number
  */
  ////////////////////////////////////////////////////////////////////////
  // PlugR (the right end)
  
  G4Tubs* plugRTubs = new G4Tubs("PlugR", // Name
			       0.0*CLHEP::cm,         // Inner radius
			       0.2356*CLHEP::cm,      // Outer radius
			       0.3*CLHEP::cm/2,       // Half length in z
			       0.*CLHEP::deg,         // Starting phi angle         
			       360.*CLHEP::deg);      // Segment angle
  
  G4LogicalVolume* plugRLogical =
    new G4LogicalVolume(plugRTubs,        // Solid
			copper,           // Material
			"PlugR_Logical"); // Name
  /*
  new G4PVPlacement(0,                                     // Rotation matrix pointer
		    G4ThreeVector(0.,0., 0.85*CLHEP::cm),  // Translation vector
		    plugRLogical,                          // Logical volume
		    "PlugR_Physical",                      // Name
		    fpWorldLogical,                        // Mother volume
		    false,                                 // Unused boolean 
		    0);                                    // Copy number
  */
  //PEEK box dimensions
  G4double boxhx=1.5;
  G4double boxhy=1.5;
  G4double boxhz=4.7/2.;
  
  ////////////////////////////////////////////////////////////////////////
  //Tungsten plate series
  G4double tunghzTot=1.675;  //changed from 1.0 to 1.34  for new structure
  G4double larGThz=0.1;
  int tungPN=10;  //changed from 6 to 8 for new structure
  int tungBPN=6;  //regulate the back plate number 4 for new structure
  G4Box* tungPlate = new G4Box("tungPlate", //Name
			       boxhx*CLHEP::cm, //half x thick
			       boxhy*CLHEP::cm, //half y thick
			       tunghzTot/(double)tungPN*CLHEP::cm); //half z thick
  G4Box* larGapTp = new G4Box("larGapTp", //Name
			       boxhx*CLHEP::cm, //half x thick
			       boxhy*CLHEP::cm, //half y thick
			       larGThz*CLHEP::cm); //half z thick
  G4LogicalVolume* tungPlateLogical =
    new G4LogicalVolume(tungPlate,  //Solid
			tungsten,   //Material
			"tungPlate_Logical");  //changed for edep investivation should be tungPlate_Logical it is Gap_Logical if we want to include tungsten for edep investigation
  G4LogicalVolume* larGapTpLogical =
    new G4LogicalVolume(larGapTp,  //Solid
			lar,   //Material
			"Gap_Logical");
  G4double LarGapcenter;
  for(int itunp=0; itunp<tungPN-1;itunp++)
    {
      //front
      new G4PVPlacement(0,                                        //Rotation matrix
			G4ThreeVector(0.,0.,(boxhz+tunghzTot/tungPN+(tunghzTot/tungPN+larGThz)*2*itunp)*CLHEP::cm),//Translation vector
			tungPlateLogical,                         //Logical volume
			"tungPlate_Physical",
			fpWorldLogical,
			false,
			0);
      LarGapcenter=(boxhz+tunghzTot/tungPN*2+larGThz+(tunghzTot/tungPN+larGThz)*2*itunp);
      new G4PVPlacement(0,                                        //Rotation matrix
			G4ThreeVector(0.,0.,LarGapcenter*CLHEP::cm),//Translation vector
			larGapTpLogical,                         //Logical volume
			"larGapTp_Physical",
			fpWorldLogical,
			false,
			0);
      //back
      if(itunp<tungBPN-1)
	{
      new G4PVPlacement(0,                                        //Rotation matrix
			G4ThreeVector(0.,0.,-(boxhz+tunghzTot/tungPN+(tunghzTot/tungPN+larGThz)*2*itunp)*CLHEP::cm),//Translation vector
			tungPlateLogical,                         //Logical volume
			"tungPlate_Physical",
			fpWorldLogical,
			false,
			0);
      new G4PVPlacement(0,                                        //Rotation matrix
			G4ThreeVector(0.,0.,-LarGapcenter*CLHEP::cm),//Translation vector
			larGapTpLogical,                         //Logical volume
			"larGapTp_Physical",
			fpWorldLogical,
			false,
			0);
	}
  G4cout<<"Tungplate series: center for Lar gap is"<<LarGapcenter<<"cm"<<G4endl;
  G4cout<<"    Lar gap is from "<<LarGapcenter-0.1<<" cm to "<<LarGapcenter+0.1<<" cm"<<G4endl;

    }
  // last piece front
      new G4PVPlacement(0,                                        //Rotation matrix
			G4ThreeVector(0.,0.,(boxhz+tunghzTot/tungPN+(tunghzTot/tungPN+larGThz)*2*(tungPN-1))*CLHEP::cm),//Translation vector
			tungPlateLogical,                         //Logical volume
			"tungPlate_Physical",
			fpWorldLogical,
			false,
			0);
  // last piece back

      new G4PVPlacement(0,                                        //Rotation matrix
			G4ThreeVector(0.,0.,-(boxhz+tunghzTot/tungPN+(tunghzTot/tungPN+larGThz)*2*(tungBPN-1))*CLHEP::cm),//Translation vector
			tungPlateLogical,                         //Logical volume
			"tungPlate_Physical",
			fpWorldLogical,
			false,
			0);  

  ////////////////////////////////////////////////////////////////////////
  //PEEK box
  const int tubcount=4;
  double tubtriside=0.75;
  double tshifty[tubcount];
  double tshiftx[tubcount];
  double tshiftz[tubcount];
  G4RotationMatrix Rth;
  G4Transform3D Trh[tubcount];
  
  for(int it=0; it<tubcount;it++)
    {
      tshifty[it]=-tubtriside/4.0*sqrt(3.)*pow(-1,it);
      tshiftx[it]=tubtriside/4.0*(it*2-3);
      tshiftz[it]=0;
      Trh[it] = G4Transform3D(Rth, G4ThreeVector(tshiftx[it]*CLHEP::cm,
						 tshifty[it]*CLHEP::cm,
						 tshiftz[it]*CLHEP::cm));

    }

  G4Box* outBox = new G4Box("outBox", //Name
			    boxhx*CLHEP::cm, //half x thick
			    boxhy*CLHEP::cm, //half y thick
			    boxhz*CLHEP::cm);//half z thick
  G4Tubs* tubHole = new G4Tubs("tubHole", //Name
			       0.0*CLHEP::cm,         // Inner radius
			       0.58/2.*CLHEP::cm,      // Outer radius
			       boxhz*CLHEP::cm,       // Half length in z
			       0.*CLHEP::deg,         // Starting phi angle         
			       360.*CLHEP::deg);      // Segment angle	
  G4SubtractionSolid *peekBox0 = new G4SubtractionSolid("peekBox0",
							outBox,
							tubHole,
							Trh[0]);
  G4SubtractionSolid *peekBox1 = new G4SubtractionSolid("peekBox1",
							peekBox0,
							tubHole,
							Trh[1]);
  G4SubtractionSolid *peekBox2 = new G4SubtractionSolid("peekBox2",
							peekBox1,
							tubHole,
							Trh[2]);
  G4SubtractionSolid *peekBox3 = new G4SubtractionSolid("peekBox3",
							peekBox2,
							tubHole,
							Trh[3]);
  G4LogicalVolume* peekBoxLogical =
    new G4LogicalVolume(peekBox3,  //Solid
			PEEK,   //Material
			"peekBox_Logical");
  
  new G4PVPlacement(0,                                        //Rotation matrix
		    G4ThreeVector(0.,0.,0.),                  //Translation vector
		    peekBoxLogical,                           //Logical volume
		    "peekBox_Physical",
		    fpWorldLogical,
		    false,
		    0);  
  					

  ////////////////////////////////////////////////////////////////////////
  //assembly start for entire tube
  G4AssemblyVolume* assemblyTube = new G4AssemblyVolume();
  G4Transform3D Tr0;
  G4RotationMatrix Ro;
  Tr0 = G4Transform3D(Ro, G4ThreeVector(0.,0.,0.));

  assemblyTube->AddPlacedVolume(
				shaftLogical,Tr0);
  Tr0 = G4Transform3D(Ro, G4ThreeVector(0.,0.,-(0.6+0.55)*CLHEP::cm));
 
  assemblyTube->AddPlacedVolume(
				shaftLLogical,Tr0);
  Tr0 = G4Transform3D(Ro, G4ThreeVector(0.,0.,(0.6+0.55)*CLHEP::cm));
 
  assemblyTube->AddPlacedVolume(
				shaftRLogical,Tr0);

  Tr0 = G4Transform3D(Ro, G4ThreeVector(0.,0.,0.));

  assemblyTube->AddPlacedVolume(
				cavityLogical,Tr0);
  assemblyTube->AddPlacedVolume(
				foilLogical, Tr0);
  assemblyTube->AddPlacedVolume(
				wallLogical,Tr0);
  assemblyTube->AddPlacedVolume(
				gapLogical,Tr0);
  assemblyTube->AddPlacedVolume(
				tubeLogical,Tr0);
  Tr0 = G4Transform3D(Ro, G4ThreeVector(0.,0.,-((0.8+1.35)/2+0.0025)*CLHEP::cm));
  assemblyTube->AddPlacedVolume(
				tubeLLogical,Tr0);
  Tr0 = G4Transform3D(Ro, G4ThreeVector(0.,0.,((0.8+1.35)/2+0.0025)*CLHEP::cm));
  assemblyTube->AddPlacedVolume(
				tubeRLogical,Tr0);

  


  for(int it=0; it<tubcount;it++)
    {
      //four tube
      Tr0 = G4Transform3D(Ro, G4ThreeVector(tshiftx[it]*CLHEP::cm,tshifty[it]*CLHEP::cm,tshiftz[it]*CLHEP::cm));
      assemblyTube->MakeImprint(fpWorldLogical, Tr0);
    }

  //single tube
  //Tr0 = G4Transform3D(Ro, G4ThreeVector(0.,0.,0.));
  //assemblyTube->MakeImprint(fpWorldLogical, Tr0);


  
  //ConstructSDandField(gapLogical);

  ////////////////////////////////////////////////////////////////////////
  // Visualisation attributes
  
  // Invisible world volume.
  fpWorldLogical->SetVisAttributes(G4VisAttributes::Invisible);

  // LAr gap - light blue
  G4VisAttributes* gapAttributes = new G4VisAttributes(G4Colour(0.0,0.5,0.5,1.0));
  gapAttributes->SetForceSolid(true);
  gapLogical->SetVisAttributes(gapAttributes);

  // copper parts - orange
  G4VisAttributes* copperAttributes = new G4VisAttributes(G4Colour(0.7,0.7,0.0,1.0));
  copperAttributes->SetForceSolid(true);

  shaftLogical->SetVisAttributes(copperAttributes);
  wallLogical->SetVisAttributes(copperAttributes);
  tubeLogical->SetVisAttributes(copperAttributes);

  // Sr90 coated copper foil - red
  G4VisAttributes* foilAttributes = new G4VisAttributes(G4Colour(1.0,0.0,0.0,1.0));
  foilAttributes->SetForceSolid(true);
  foilLogical->SetVisAttributes(foilAttributes);

  // tracking lenght in each volume

  G4double maxStep = 0.0001*CLHEP::cm; // 1 mu tracking
  G4UserLimits * limits = new G4UserLimits(maxStep);

  fpWorldLogical->SetUserLimits(limits);

  shaftLogical->SetUserLimits(limits);
  cavityLogical->SetUserLimits(limits);
  foilLogical->SetUserLimits(limits);
  spaceLogical->SetUserLimits(limits);
  wallLogical->SetUserLimits(limits);
  gapLogical->SetUserLimits(limits);
  tubeLogical->SetUserLimits(limits);  
  plugLLogical->SetUserLimits(limits);  
  plugRLogical->SetUserLimits(limits);  

}




//Bot

//Geant4 headers
#include "lampslowDetectorConstruction.hh"
#include "SiCsISD.hh"
#include "NeutronSD.hh"

//Geant4 headers
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4Trap.hh"
#include "G4TwistTrapAlphaSide.hh"
#include "G4TwistedTrap.hh"
#include "G4Tubs.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4FieldManager.hh"
#include <TMath.h>
#include <sstream>
#include <iostream>

//ReadOut
//#include "SiCsIROGeometry.hh"

using std::stringstream;
using namespace std;
using namespace CLHEP;

lampslowDetectorConstruction::lampslowDetectorConstruction()
//:sicsi_PVName("SiCsIPV"), block_PVName("BlockPV")
{
  ConstructMaterials();
  DefineDimensions();
}

lampslowDetectorConstruction::~lampslowDetectorConstruction()
{
  DestructMaterials();
}

void lampslowDetectorConstruction::DefineDimensions()
{
  //Lab
  labX = 1500*CLHEP::mm;
  labY = 1500*CLHEP::mm;
  labZ = 2500*CLHEP::mm;

  // Size of Detector
  Sisize = 9.00*CLHEP::cm;
  CsIsize = 9.00*CLHEP::cm;
  Vetosize = 9.00*CLHEP::cm;

  // Thickness of Detector
  SiThick = 100*CLHEP::um;
  CsIThick = 5*CLHEP::cm;
  VetoThick = 300*CLHEP::um;

  testblock_X = 150*CLHEP::mm;
  testblock_Y = 150*CLHEP::mm;
  testblock_Z = 200*CLHEP::mm;
}

G4VPhysicalVolume* lampslowDetectorConstruction::Construct()
{
  G4VSolid* labSolid = new G4Box("labSolid", labX/2, labY/2, labZ/2);
  G4LogicalVolume* labLV = new G4LogicalVolume(labSolid, Vacuum, "labLV");
  G4VPhysicalVolume* labPV = new G4PVPlacement(0, G4ThreeVector(), "labPV", labLV, 0, false, 0);

  // ######### Test block detector #########  - jaebeom
 
  //******Tracking Sphere******

  //Tracking Vacuum
/*  G4Sphere *Track_Sphere = new G4Sphere("Track_Sphere", 110*CLHEP::cm, 110.2*CLHEP::cm, 0, CLHEP::twopi, 0, pi);
  G4LogicalVolume *Track_LVSphere = new G4LogicalVolume(Track_Sphere, Vacuum, "Track_LVSphere");
  G4VPhysicalVolume* Track_PV = new G4PVPlacement(0,G4ThreeVector(0,0,400*CLHEP::mm+110*CLHEP::um+1.25*CLHEP::cm),"Track_PV",Track_LVSphere,labPV,false,0);

  G4VisAttributes* Track_VisAttrib_LVSphere = new G4VisAttributes(G4Colour(1., 0., 0.));
  Track_LVSphere -> SetVisAttributes(Track_VisAttrib_LVSphere);
    G4VSolid* Track_Block = new G4Box("labSolid", 100*CLHEP::cm, 100*CLHEP::cm, 100*CLHEP::cm);
    G4LogicalVolume *Track_LVSphere = new G4LogicalVolume(Track_Block, Vacuum, "Track_LVSphere");
    G4VPhysicalVolume* Track_PV = new G4PVPlacement(0,G4ThreeVector(),"Track_PV",Track_LVSphere,labPV,false,0);

	  G4VisAttributes* Track_VisAttrib_LVSphere = new G4VisAttributes(G4Colour(1., 0., 0.));
  	Track_LVSphere -> SetVisAttributes(Track_VisAttrib_LVSphere);

  //Neutron Block test
  G4Box* BlockSolid = new G4Box("Block", testblock_X/2, testblock_Y/2, testblock_Z/2);
  G4LogicalVolume* BlockLV = new G4LogicalVolume(BlockSolid, Scint, "BlockLV");
  new G4PVPlacement(0, G4ThreeVector(0,0,400*CLHEP::mm+testblock_Z/2+110*um+50.1*CLHEP::mm), "testblock_PV", BlockLV, labPV, false, 0);

  G4VisAttributes* BlockTest = new G4VisAttributes(G4Colour(0.,1.,0.));
  BlockLV -> SetVisAttributes(BlockTest);
*/

  //veto test
  G4Box* TestVetoBlock = new G4Box("TestVeto", 9.0*CLHEP::cm/2, 9.0*CLHEP::cm/2, 300*CLHEP::um/2);
  G4LogicalVolume* TestVetoBlockLV = new G4LogicalVolume(TestVetoBlock, Scint, "TestVetoBlockLV");
  new G4PVPlacement(0, G4ThreeVector(0,0,40*cm+CsIThick+SiThick+VetoThick+6*cm), "TestVetoBlock_PV", TestVetoBlockLV, labPV, false,0);

  G4VisAttributes* TestVetoBlockSD = new G4VisAttributes(G4Colour(1.,1.,0.));
  TestVetoBlockLV -> SetVisAttributes(TestVetoBlockSD);


  //Silicon test
  G4Box* TestSiBlock_1 = new G4Box("TestSi_1", 9.0*CLHEP::cm/2, 9.0*CLHEP::cm/2, 100*CLHEP::um/2);
  G4LogicalVolume* TestSiBlockLV_1 = new G4LogicalVolume(TestSiBlock_1, Silicon, "TestSiBlockLV_1");
  new G4PVPlacement(0, G4ThreeVector(0,0,400*mm), "TestSiBlock_PV_1", TestSiBlockLV_1, labPV, false,0);

  G4VisAttributes* TestSiBlockSD = new G4VisAttributes(G4Colour(1.,0.,0.));
  TestSiBlockLV_1 -> SetVisAttributes(TestSiBlockSD);
/*  
  G4Box* TestSiBlock_2 = new G4Box("TestSi_2", 9.0*CLHEP::cm/2, 9.0*CLHEP::cm/2, 400*CLHEP::um/2);
  G4LogicalVolume* TestSiBlockLV_2 = new G4LogicalVolume(TestSiBlock_2, Silicon, "TestSiBlockLV_2");
  new G4PVPlacement(0, G4ThreeVector(0,0,400*mm+1*cm), "TestSiBlock_PV_2", TestSiBlockLV_2, labPV, false,0);

  TestSiBlockLV_2 -> SetVisAttributes(TestSiBlockSD);

  G4Box* TestSiBlock_3 = new G4Box("TestSi_3", 9.0*CLHEP::cm/2, 9.0*CLHEP::cm/2, 400*CLHEP::um/2);
  G4LogicalVolume* TestSiBlockLV_3 = new G4LogicalVolume(TestSiBlock_3, Silicon, "TestSiBlockLV_3");
  new G4PVPlacement(0, G4ThreeVector(0,0,400*mm+2*cm), "TestSiBlock_PV_3", TestSiBlockLV_3, labPV, false,0);

  TestSiBlockLV_3 -> SetVisAttributes(TestSiBlockSD);
*/  
//CsI test
  G4Box* TestCsIBlock = new G4Box("TestCsI", 9.0*CLHEP::cm/2, 9.0*CLHEP::cm/2, 5.0*CLHEP::cm/2);
  G4LogicalVolume* TestCsIBlockLV = new G4LogicalVolume(TestCsIBlock, CsI, "TestCsIBlockLV");
  new G4PVPlacement(0, G4ThreeVector(0,0,400*mm+CsIThick+SiThick+1*cm), "TestCsIBlock_PV", TestCsIBlockLV, labPV, false,0);

  G4VisAttributes* TestCsIBlockSD = new G4VisAttributes(G4Colour(0.,0.,1.));
  TestCsIBlockLV -> SetVisAttributes(TestCsIBlockSD);

  //Det Register
  G4SDManager* sdManager = G4SDManager::GetSDMpointer();
/*
  NeutronSD* NSD = new NeutronSD("/Neutron/NeutronSD","Neutron");
  sdManager -> AddNewDetector(NSD);
  BlockLV -> SetSensitiveDetector(NSD);
*/
  SiCsISD* VetoSD = new SiCsISD("/SiCsI/VetoSD","Veto");
  sdManager -> AddNewDetector(VetoSD);
  TestVetoBlockLV -> SetSensitiveDetector(VetoSD);
  
  SiCsISD* SiSD = new SiCsISD("/SiCsI/SiSD","Si");
  sdManager -> AddNewDetector(SiSD);
  TestSiBlockLV_1 -> SetSensitiveDetector(SiSD);
 // TestSiBlockLV_2 -> SetSensitiveDetector(SiSD);
  //TestSiBlockLV_3 -> SetSensitiveDetector(SiSD);
  
  SiCsISD* CsISD = new SiCsISD("/SiCsI/CsISD","CsI");
  sdManager -> AddNewDetector(CsISD);
  TestCsIBlockLV -> SetSensitiveDetector(CsISD);
 // Track_LVSphere -> SetSensitiveDetector(CsISD);

  

  //Define detectors
//  ConstructSiCsI(labPV, sdManager);
//  ConstructBlock(labPV, sdManager);

  return labPV;
}

//SiCsI
void lampslowDetectorConstruction::ConstructSiCsI(G4VPhysicalVolume* labPV, G4SDManager* sdManager)
{
	


  //############ Define Parameters ################# - jaebeom
   
	// Detector Size -  polar angle 17.5 ~ 77.5 degree
	const G4double SiSize_1 = 9.00*CLHEP::cm;
	const G4double CsISize_1 = 9.00*CLHEP::cm;

	// Detector Size -  polar angle 77.5 ~ 150 degree
	const G4double SiSize_2 = 15.00*CLHEP::cm;
	const G4double CsISize_2 = 15.00*CLHEP::cm;



	// Si, CsI Thickness
	const G4double SiThick_Z = 100*CLHEP::um;
	const G4double CsIThick_Z = 5*CLHEP::cm;
	
	// Angle
	G4double Angle1 = 25*CLHEP::degree;
	G4double Angle2 = 40*CLHEP::degree;
	G4double Angle3 = 55*CLHEP::degree;
	G4double Angle4 = 70*CLHEP::degree;
	G4double Angle5 = 90*CLHEP::degree;
	G4double Angle6 = 114*CLHEP::degree;
	G4double Angle7 = 138*CLHEP::degree;

	// Number of Detectors per unit polar angle
	const G4int CsIPhiN_1 = 8;
	const G4int CsIPhiN_2 = 12;
	const G4int CsIPhiN_3 = 18;
	const G4int CsIPhiN_4 = 20;
	const G4int CsIPhiN_5 = 15;
	const G4int CsIPhiN_6 = 12;
	const G4int CsIPhiN_7 = 8;
  

	// Voxel Size
  const G4double X_VOXEL_SIZE_SiCsI = 10*CLHEP::mm;
  const G4double Y_VOXEL_SIZE_SiCsI = 10*CLHEP::mm;
  const G4double Z_VOXEL_SIZE_SiCsI = 10*CLHEP::mm;
	
  
  
  //############ Define Detectors ################# - jaebeom


  //============================1st ring==============================================
	//Si1
    G4Box *Si_Solid_1 = new G4Box("Si_Solid_1", SiSize_1/2, SiSize_1/2, SiThick_Z/2);
	  G4LogicalVolume *Si_LV_1 = new G4LogicalVolume(Si_Solid_1, Silicon, "Si_LV_1");
	
	for(G4int copyN=0; copyN<CsIPhiN_1; copyN++)
	{
		new G4PVPlacement(new G4RotationMatrix((360/CsIPhiN_1)*copyN*CLHEP::degree, Angle1, 0*CLHEP::degree), G4ThreeVector(0*CLHEP::mm, 0*CLHEP::mm, 400*CLHEP::mm).rotateX(Angle1).rotateZ((360/CsIPhiN_1)*copyN*CLHEP::degree),"Si_PV", Si_LV_1, labPV, false, 111+copyN);
	}

	G4VisAttributes* Si_VisAttrib = new G4VisAttributes(G4Colour(1., 0., 0.));
	Si_LV_1 -> SetVisAttributes(Si_VisAttrib);


	//CsI1
	
	G4Box *CsI_Solid_1 = new G4Box("CsI_Solid_1", CsISize_1/2, CsISize_1/2, CsIThick_Z/2);
	G4LogicalVolume *CsI_LV_1 = new G4LogicalVolume(CsI_Solid_1, CsI, "CsI_LV_1");
	
	for(G4int copyN=0; copyN<CsIPhiN_1; copyN++)
	{
		new G4PVPlacement(new G4RotationMatrix((360/CsIPhiN_1)*copyN*CLHEP::degree, Angle1, 0*CLHEP::degree), G4ThreeVector(0*CLHEP::mm, 0*CLHEP::mm, 410*CLHEP::mm+SiThick_Z/2+CsIThick_Z/2).rotateX(Angle1).rotateZ((360/CsIPhiN_1)*copyN*CLHEP::degree), "CsI_PV", CsI_LV_1, labPV, false, 121+copyN);
	}

	G4VisAttributes* CsI_VisAttrib = new G4VisAttributes(G4Colour(0., 0., 1.));
	CsI_LV_1 -> SetVisAttributes(CsI_VisAttrib);
	

  
  //*******Trapezoid*********
/*
  //Si
    G4Trap *Si_Trd_1 = new G4Trap("Si_Trd_1", SiThick_Z/2, 0, 0, 5*CLHEP::cm, 4.5*CLHEP::cm, 8*CLHEP::cm, 0, 5*CLHEP::cm, 4.5*CLHEP::cm, 8*CLHEP::cm, 0);
    G4LogicalVolume *Si_LVTrd_1 = new G4LogicalVolume(Si_Trd_1, Silicon, "Si_LVTrd_1");

	
  for(G4int copyN=0; copyN<CsIPhiN_1; copyN++)
	{
		new G4PVPlacement(new G4RotationMatrix((360/CsIPhiN_1)*copyN*degree, -Angle1, 0*degree), G4ThreeVector(0*CLHEP::mm, 0*CLHEP::mm, 400*CLHEP::mm).rotateX(-Angle1).rotateZ((360/CsIPhiN_1)*copyN*degree),"Si_PV", Si_LVTrd_1, labPV, false, 111+copyN);
	}

	G4VisAttributes* Si_VisAttrib_Trd = new G4VisAttributes(G4Colour(1., 0., 0.));
	Si_LVTrd_1 -> SetVisAttributes(Si_VisAttrib_Trd);


  //CsI1

	G4Trap *CsI_Trd_1 = new G4Trap("CsI_Trd_1", CsIThick_Z/2, 0, 0, 5*CLHEP::cm, 4.5*CLHEP::cm, 8*CLHEP::cm, 0, 5*CLHEP::cm, 4.5*CLHEP::cm, 8*CLHEP::cm, 0);
	G4LogicalVolume *CsI_LVTrd_1 = new G4LogicalVolume(CsI_Trd_1, CsI, "CsI_LVTrd_1");
	
	for(G4int copyN=0; copyN<CsIPhiN_1; copyN++)
	{
		new G4PVPlacement(new G4RotationMatrix((360/CsIPhiN_1)*copyN*degree, -Angle1, 0*degree), G4ThreeVector(0*CLHEP::mm, 0*CLHEP::mm, 425*CLHEP::mm+SiThick_Z/2+CsIThick_Z/2).rotateX(-Angle1).rotateZ((360/CsIPhiN_1)*copyN*degree), "CsI_PV", CsI_LVTrd_1, labPV, false, 121+copyN);
	}

	G4VisAttributes* CsI_VisAttribTrd = new G4VisAttributes(G4Colour(0., 0., 1.));
	CsI_LVTrd_1 -> SetVisAttributes(CsI_VisAttribTrd);

*/
  //******Sphere******
/* 
  //Si
    G4Sphere *Si_Sphere_1 = new G4Sphere("Si_Sphere_1", 40.*CLHEP::cm, 40.01*CLHEP::cm, 0, twopi, 0, pi);
    G4LogicalVolume *Si_LVSphere_1 = new G4LogicalVolume(Si_Sphere_1, Silicon, "Si_LVSphere_1");
    new G4PVPlacement(0,G4ThreeVector(),"Si_PV",Si_LVSphere_1,labPV,false,0);

	  G4VisAttributes* Si_VisAttrib_LVSphere = new G4VisAttributes(G4Colour(1., 0., 0.));
  	Si_LVSphere_1 -> SetVisAttributes(Si_VisAttrib_LVSphere);


  //Csi
    G4Sphere *CsI_Sphere_1 = new G4Sphere("CsI_Sphere_1", 41.01*CLHEP::cm, 46.01*CLHEP::cm, 0, twopi, 0, pi);
    G4LogicalVolume *CsI_LVSphere_1 = new G4LogicalVolume(CsI_Sphere_1, CsI, "CsI_LVSphere_1");
    new G4PVPlacement(0,G4ThreeVector(),"CsI_PV",CsI_LVSphere_1,labPV,false,0);

	  G4VisAttributes* CsI_VisAttrib_LVSphere = new G4VisAttributes(G4Colour(1., 0., 0.));
  	CsI_LVSphere_1 -> SetVisAttributes(CsI_VisAttrib_LVSphere);
*/

 
  //=========================2nd ring================================================
	//Si2
    G4Box *Si_Solid_2 = new G4Box("Si_Solid_2", SiSize_1/2, SiSize_1/2, SiThick_Z/2);
	  G4LogicalVolume *Si_LV_2 = new G4LogicalVolume(Si_Solid_2, Silicon, "Si_LV_2");
	
	for(G4int copyN=0; copyN<CsIPhiN_2; copyN++)
	{
		new G4PVPlacement(new G4RotationMatrix((360/CsIPhiN_2)*copyN*CLHEP::degree, Angle2, 0*CLHEP::degree), G4ThreeVector(0*CLHEP::mm, 0*CLHEP::mm, 400*CLHEP::mm).rotateX(Angle2).rotateZ((360/CsIPhiN_2)*copyN*CLHEP::degree),"Si_PV", Si_LV_2, labPV, false, 211+copyN);
	}

	Si_LV_2 -> SetVisAttributes(Si_VisAttrib);


	//CsI2
	
	G4Box *CsI_Solid_2 = new G4Box("CsI_Solid_2", CsISize_1/2, CsISize_1/2, CsIThick_Z/2);
	G4LogicalVolume *CsI_LV_2 = new G4LogicalVolume(CsI_Solid_2, CsI, "CsI_LV_2");
	
	for(G4int copyN=0; copyN<CsIPhiN_2; copyN++)
	{
		new G4PVPlacement(new G4RotationMatrix((360/CsIPhiN_2)*copyN*CLHEP::degree, Angle2, 0*CLHEP::degree), G4ThreeVector(0*CLHEP::mm, 0*CLHEP::mm, 410*CLHEP::mm+SiThick_Z/2+CsIThick_Z/2+1*CLHEP::mm).rotateX(Angle2).rotateZ((360/CsIPhiN_2)*copyN*CLHEP::degree), "CsI_PV", CsI_LV_2, labPV, false, 221+copyN);
	}

	CsI_LV_2 -> SetVisAttributes(CsI_VisAttrib);
  

  //=========================3rd ring================================================
	//Si3
    G4Box *Si_Solid_3 = new G4Box("Si_Solid_3", SiSize_1/2, SiSize_1/2, SiThick_Z/2);
	  G4LogicalVolume *Si_LV_3 = new G4LogicalVolume(Si_Solid_3, Silicon, "Si_LV_3");
	
	for(G4int copyN=0; copyN<CsIPhiN_3; copyN++)
	{
		new G4PVPlacement(new G4RotationMatrix((360/CsIPhiN_3)*copyN*CLHEP::degree, Angle3, 0*CLHEP::degree), G4ThreeVector(0*CLHEP::mm, 0*CLHEP::mm, 400*CLHEP::mm).rotateX(Angle3).rotateZ((360/CsIPhiN_3)*copyN*CLHEP::degree),"Si_PV", Si_LV_3, labPV, false, 311+copyN);
	}

	Si_LV_3 -> SetVisAttributes(Si_VisAttrib);


	//CsI3
	
	G4Box *CsI_Solid_3 = new G4Box("CsI_Solid_3", CsISize_1/2, CsISize_1/2, CsIThick_Z/2);
	G4LogicalVolume *CsI_LV_3 = new G4LogicalVolume(CsI_Solid_3, CsI, "CsI_LV_3");
	
	for(G4int copyN=0; copyN<CsIPhiN_3; copyN++)
	{
		new G4PVPlacement(new G4RotationMatrix((360/CsIPhiN_3)*copyN*CLHEP::degree, Angle3, 0*CLHEP::degree), G4ThreeVector(0*CLHEP::mm, 0*CLHEP::mm, 410*CLHEP::mm+SiThick_Z/2+CsIThick_Z/2+1*CLHEP::mm).rotateX(Angle3).rotateZ((360/CsIPhiN_3)*copyN*CLHEP::degree), "CsI_PV", CsI_LV_3, labPV, false, 321+copyN);
	}

	CsI_LV_3 -> SetVisAttributes(CsI_VisAttrib);
  
  
  //=========================4th ring================================================
	//Si4
    G4Box *Si_Solid_4 = new G4Box("Si_Solid_4", SiSize_1/2, SiSize_1/2, SiThick_Z/2);
	  G4LogicalVolume *Si_LV_4 = new G4LogicalVolume(Si_Solid_4, Silicon, "Si_LV_4");
	
	for(G4int copyN=0; copyN<CsIPhiN_4; copyN++)
	{
		new G4PVPlacement(new G4RotationMatrix((360/CsIPhiN_4)*copyN*CLHEP::degree, Angle4, 0*CLHEP::degree), G4ThreeVector(0*CLHEP::mm, 0*CLHEP::mm, 400*CLHEP::mm).rotateX(Angle4).rotateZ((360/CsIPhiN_4)*copyN*CLHEP::degree),"Si_PV", Si_LV_4, labPV, false, 411+copyN);
	}

	Si_LV_4 -> SetVisAttributes(Si_VisAttrib);


	//CsI4
	
	G4Box *CsI_Solid_4 = new G4Box("CsI_Solid_4", CsISize_1/2, CsISize_1/2, CsIThick_Z/2);
	G4LogicalVolume *CsI_LV_4 = new G4LogicalVolume(CsI_Solid_4, CsI, "CsI_LV_4");
	
	for(G4int copyN=0; copyN<CsIPhiN_4; copyN++)
	{
		new G4PVPlacement(new G4RotationMatrix((360/CsIPhiN_4)*copyN*CLHEP::degree, Angle4, 0*CLHEP::degree), G4ThreeVector(0*CLHEP::mm, 0*CLHEP::mm, 410*CLHEP::mm+SiThick_Z/2+CsIThick_Z/2+1*CLHEP::mm).rotateX(Angle4).rotateZ((360/CsIPhiN_4)*copyN*CLHEP::degree), "CsI_PV", CsI_LV_4, labPV, false, 421+copyN);
	}

	CsI_LV_4 -> SetVisAttributes(CsI_VisAttrib);
 

  //=========================5th ring================================================
	//Si5
    G4Box *Si_Solid_5 = new G4Box("Si_Solid_5", SiSize_2/2, SiSize_2/2, SiThick_Z/2);
	  G4LogicalVolume *Si_LV_5 = new G4LogicalVolume(Si_Solid_5, Silicon, "Si_LV_5");
	
	for(G4int copyN=0; copyN<CsIPhiN_5; copyN++)
	{
		new G4PVPlacement(new G4RotationMatrix((360/CsIPhiN_5)*copyN*CLHEP::degree, Angle5, 0*CLHEP::degree), G4ThreeVector(0*CLHEP::mm, 0*CLHEP::mm, 400*CLHEP::mm).rotateX(Angle5).rotateZ((360/CsIPhiN_5)*copyN*CLHEP::degree),"Si_PV", Si_LV_5, labPV, false, 511+copyN);
	}

	Si_LV_5 -> SetVisAttributes(Si_VisAttrib);


	//CsI5
	
	G4Box *CsI_Solid_5 = new G4Box("CsI_Solid_5", CsISize_2/2, CsISize_2/2, CsIThick_Z/2);
	G4LogicalVolume *CsI_LV_5 = new G4LogicalVolume(CsI_Solid_5, CsI, "CsI_LV_5");
	
	for(G4int copyN=0; copyN<CsIPhiN_5; copyN++)
	{
		new G4PVPlacement(new G4RotationMatrix((360/CsIPhiN_5)*copyN*CLHEP::degree, Angle5, 0*CLHEP::degree), G4ThreeVector(0*CLHEP::mm, 0*CLHEP::mm, 410*CLHEP::mm+SiThick_Z/2+CsIThick_Z/2+1*CLHEP::mm).rotateX(Angle5).rotateZ((360/CsIPhiN_5)*copyN*CLHEP::degree), "CsI_PV", CsI_LV_5, labPV, false, 521+copyN);
	}

	CsI_LV_5 -> SetVisAttributes(CsI_VisAttrib);
 

  //=========================6th ring================================================
	//Si6
    G4Box *Si_Solid_6 = new G4Box("Si_Solid_6", SiSize_2/2, SiSize_2/2, SiThick_Z/2);
	  G4LogicalVolume *Si_LV_6 = new G4LogicalVolume(Si_Solid_6, Silicon, "Si_LV_6");
	
	for(G4int copyN=0; copyN<CsIPhiN_6; copyN++)
	{
		new G4PVPlacement(new G4RotationMatrix((360/CsIPhiN_6)*copyN*CLHEP::degree, Angle6, 0*CLHEP::degree), G4ThreeVector(0*CLHEP::mm, 0*CLHEP::mm, 400*CLHEP::mm).rotateX(Angle6).rotateZ((360/CsIPhiN_6)*copyN*CLHEP::degree),"Si_PV", Si_LV_6, labPV, false, 611+copyN);
	}

	Si_LV_6 -> SetVisAttributes(Si_VisAttrib);


	//CsI6
	
	G4Box *CsI_Solid_6 = new G4Box("CsI_Solid_6", CsISize_2/2, CsISize_2/2, CsIThick_Z/2);
	G4LogicalVolume *CsI_LV_6 = new G4LogicalVolume(CsI_Solid_6, CsI, "CsI_LV_6");
	
	for(G4int copyN=0; copyN<CsIPhiN_6; copyN++)
	{
		new G4PVPlacement(new G4RotationMatrix((360/CsIPhiN_6)*copyN*CLHEP::degree, Angle6, 0*CLHEP::degree), G4ThreeVector(0*CLHEP::mm, 0*CLHEP::mm, 410*CLHEP::mm+SiThick_Z/2+CsIThick_Z/2+1*CLHEP::mm).rotateX(Angle6).rotateZ((360/CsIPhiN_6)*copyN*CLHEP::degree), "CsI_PV", CsI_LV_6, labPV, false, 621+copyN);
	}

	CsI_LV_6 -> SetVisAttributes(CsI_VisAttrib);
  
  //=========================7th ring================================================
	//Si7
    G4Box *Si_Solid_7 = new G4Box("Si_Solid_7", SiSize_2/2, SiSize_2/2, SiThick_Z/2);
	  G4LogicalVolume *Si_LV_7 = new G4LogicalVolume(Si_Solid_7, Silicon, "Si_LV_7");
	
	for(G4int copyN=0; copyN<CsIPhiN_7; copyN++)
	{
		new G4PVPlacement(new G4RotationMatrix((360/CsIPhiN_7)*copyN*CLHEP::degree, Angle7, 0*CLHEP::degree), G4ThreeVector(0*CLHEP::mm, 0*CLHEP::mm, 400*CLHEP::mm).rotateX(Angle7).rotateZ((360/CsIPhiN_7)*copyN*CLHEP::degree),"Si_PV", Si_LV_7, labPV, false, 711+copyN);
	}

	Si_LV_7 -> SetVisAttributes(Si_VisAttrib);


	//CsI7
	
	G4Box *CsI_Solid_7 = new G4Box("CsI_Solid_7", CsISize_2/2, CsISize_2/2, CsIThick_Z/2);
	G4LogicalVolume *CsI_LV_7 = new G4LogicalVolume(CsI_Solid_7, CsI, "CsI_LV_7");
	
	for(G4int copyN=0; copyN<CsIPhiN_7; copyN++)
	{
		new G4PVPlacement(new G4RotationMatrix((360/CsIPhiN_7)*copyN*CLHEP::degree, Angle7, 0*CLHEP::degree), G4ThreeVector(0*CLHEP::mm, 0*CLHEP::mm, 410*CLHEP::mm+SiThick_Z/2+CsIThick_Z/2+1*CLHEP::mm).rotateX(Angle7).rotateZ((360/CsIPhiN_7)*copyN*CLHEP::degree), "CsI_PV", CsI_LV_7, labPV, false, 721+copyN);
	}

	CsI_LV_7 -> SetVisAttributes(CsI_VisAttrib);

  /*
  SiCsIROGeometry* SICSI_RO_GEOMETRY = new SiCsIROGeometry("SICSI_RO_GEOMETRY",
	                                            							labX,
							                                            	labY,
							                                            	labZ,
							                                            	SiThick_Z,
							                                            	CsIThick_Z,
								                                            X_VOXEL_SIZE_SiCsI,
			                                            					Y_VOXEL_SIZE_SiCsI,
			                                            					Z_VOXEL_SIZE_SiCsI,
		                                            						SiSize_1,
					                                            			CsISize_1,
						                                            		Angle1,
								                                            CsIPhiN_1);

	SICSI_RO_GEOMETRY -> BuildROGeometry();
*/



  // ############# Detector Registration ################# - jaebeom

  SiCsISD* SiSD = new SiCsISD("/SiCsI/SiSD","Si");
  sdManager -> AddNewDetector(SiSD);

  SiCsISD* CsISD = new SiCsISD("/SiCsI/CsISD","CsI");
  sdManager -> AddNewDetector(CsISD);

  // Silicon   
  Si_LV_1->SetSensitiveDetector(SiSD);
  Si_LV_2->SetSensitiveDetector(SiSD);
  Si_LV_3->SetSensitiveDetector(SiSD);
  Si_LV_4->SetSensitiveDetector(SiSD);
  Si_LV_5->SetSensitiveDetector(SiSD);
  Si_LV_6->SetSensitiveDetector(SiSD);
  Si_LV_7->SetSensitiveDetector(SiSD);

 
  // CsI
  CsI_LV_1->SetSensitiveDetector(CsISD);
  CsI_LV_2->SetSensitiveDetector(CsISD);
  CsI_LV_3->SetSensitiveDetector(CsISD);
  CsI_LV_4->SetSensitiveDetector(CsISD);
  CsI_LV_5->SetSensitiveDetector(CsISD);
  CsI_LV_6->SetSensitiveDetector(CsISD);
  CsI_LV_7->SetSensitiveDetector(CsISD);


  //Sphere
 // Si_LVSphere_1->SetSensitiveDetector(SiSD);

  //CsI_LVSphere_1->SetSensitiveDetector(CsISD);

  //trapezoid
  //Si_LVTrd_1->SetSensitiveDetector(SiSD);

  //trapezoid
  //CsI_LVTrd_1->SetSensitiveDetector(CsISD);


}


void lampslowDetectorConstruction::ConstructBlock(G4VPhysicalVolume* labPV, G4SDManager* sdManager)
{

  //########### Define Parameters ################ - jaebeom
  
  const G4double TargetRadiusDistance = 1500*CLHEP::mm;

  G4double NeutronSize    = 150*CLHEP::mm;
  G4double NeutronThick_Z = 200*CLHEP::mm;
  
  G4double NeutronPhiNum[10]={20,30,36,43,49,53,57,60,61,61};
  G4double Angle[10]={24*CLHEP::degree,32*CLHEP::degree,40*CLHEP::degree,48*CLHEP::degree,56*CLHEP::degree,64*CLHEP::degree,72*CLHEP::degree,80*CLHEP::degree,88*CLHEP::degree,96*CLHEP::degree};

  G4int n;

  //check
  for(int i=0; i<10; i++)
    {
     cout << "Neutron Phi Number " << i << " : " << NeutronPhiNum[i] << endl;
     cout << "Center Angle " << i << " : " << Angle[i] << endl;
    }


  
/*

  //########### Define Detectors ################  - jaebeom
  
  
  //============================ 1st ring : theta = 20 ~ 28 degree ==============================================

  n=0;

	G4Box *Neutron_Solid_1 = new G4Box("Neutron_Solid_1", NeutronSize/2, NeutronSize/2, NeutronThick_Z/2);
	G4LogicalVolume *Neutron_LV_1 = new G4LogicalVolume(Neutron_Solid_1, Scint, "Neutron_LV_1");
	
	for(G4int copyN=0; copyN<NeutronPhiNum[n]; copyN++)
	{
		new G4PVPlacement(new G4RotationMatrix((360/NeutronPhiNum[n])*copyN*CLHEP::degree, Angle[n], 0*CLHEP::degree), G4ThreeVector(0*CLHEP::mm, 0*CLHEP::mm, TargetRadiusDistance+NeutronThick_Z/2).rotateX(Angle[n]).rotateZ((360/NeutronPhiNum[n])*copyN*CLHEP::degree), "Neutron_PV", Neutron_LV_1, labPV, false, 11+copyN);
	}

	G4VisAttributes* Neutron_VisAttrib = new G4VisAttributes(G4Colour(0., 1., 0.));
	Neutron_LV_1 -> SetVisAttributes(Neutron_VisAttrib);
    
  
  
  
  
  //============================ 2rd ring : theta = 28 ~ 36 degree ==============================================

  n=1;

	G4Box *Neutron_Solid_2 = new G4Box("Neutron_Solid_2", NeutronSize/2, NeutronSize/2, NeutronThick_Z/2);
	G4LogicalVolume *Neutron_LV_2 = new G4LogicalVolume(Neutron_Solid_2, Scint, "Neutron_LV_2");
	
	for(G4int copyN=0; copyN<NeutronPhiNum[n]; copyN++)
	{
		new G4PVPlacement(new G4RotationMatrix((360/NeutronPhiNum[n])*copyN*CLHEP::degree, Angle[n], 0*CLHEP::degree), G4ThreeVector(0*CLHEP::mm, 0*CLHEP::mm, TargetRadiusDistance+NeutronThick_Z/2).rotateX(Angle[n]).rotateZ((360/NeutronPhiNum[n])*copyN*CLHEP::degree), "Neutron_PV", Neutron_LV_2, labPV, false, 21+copyN);
	}

	Neutron_LV_2 -> SetVisAttributes(Neutron_VisAttrib);
  



  //============================ 3rd ring : theta = 36 ~ 44 degree ==============================================

  n=2;

	G4Box *Neutron_Solid_3 = new G4Box("Neutron_Solid_3", NeutronSize/2, NeutronSize/2, NeutronThick_Z/2);
	G4LogicalVolume *Neutron_LV_3 = new G4LogicalVolume(Neutron_Solid_3, Scint, "Neutron_LV_3");
	
	for(G4int copyN=0; copyN<NeutronPhiNum[n]; copyN++)
	{
		new G4PVPlacement(new G4RotationMatrix((360/NeutronPhiNum[n])*copyN*CLHEP::degree, Angle[n], 0*CLHEP::degree), G4ThreeVector(0*CLHEP::mm, 0*CLHEP::mm, TargetRadiusDistance+NeutronThick_Z/2).rotateX(Angle[n]).rotateZ((360/NeutronPhiNum[n])*copyN*CLHEP::degree), "Neutron_PV", Neutron_LV_3, labPV, false, 31+copyN);
	}

	Neutron_LV_3 -> SetVisAttributes(Neutron_VisAttrib);
    
  


  //============================ 4th ring : theta = 44 ~ 52 degree ==============================================

  n=3;

	G4Box *Neutron_Solid_4 = new G4Box("Neutron_Solid_4", NeutronSize/2, NeutronSize/2, NeutronThick_Z/2);
	G4LogicalVolume *Neutron_LV_4 = new G4LogicalVolume(Neutron_Solid_4, Scint, "Neutron_LV_4");
	
	for(G4int copyN=0; copyN<NeutronPhiNum[n]; copyN++)
	{
		new G4PVPlacement(new G4RotationMatrix((360/NeutronPhiNum[n])*copyN*CLHEP::degree, Angle[n], 0*CLHEP::degree), G4ThreeVector(0*CLHEP::mm, 0*CLHEP::mm, TargetRadiusDistance+NeutronThick_Z/2).rotateX(Angle[n]).rotateZ((360/NeutronPhiNum[n])*copyN*CLHEP::degree), "Neutron_PV", Neutron_LV_4, labPV, false, 41+copyN);
	}

	Neutron_LV_4 -> SetVisAttributes(Neutron_VisAttrib);
    
  


  //============================ 5th ring : theta = 52 ~ 60 degree ==============================================

  n=4;

	G4Box *Neutron_Solid_5 = new G4Box("Neutron_Solid_5", NeutronSize/2, NeutronSize/2, NeutronThick_Z/2);
	G4LogicalVolume *Neutron_LV_5 = new G4LogicalVolume(Neutron_Solid_5, Scint, "Neutron_LV_5");
	
	for(G4int copyN=0; copyN<NeutronPhiNum[n]; copyN++)
	{
		new G4PVPlacement(new G4RotationMatrix((360/NeutronPhiNum[n])*copyN*CLHEP::degree, Angle[n], 0*CLHEP::degree), G4ThreeVector(0*CLHEP::mm, 0*CLHEP::mm, TargetRadiusDistance+NeutronThick_Z/2).rotateX(Angle[n]).rotateZ((360/NeutronPhiNum[n])*copyN*CLHEP::degree), "Neutron_PV", Neutron_LV_5, labPV, false, 51+copyN);
	}

	Neutron_LV_5 -> SetVisAttributes(Neutron_VisAttrib);
    
  



  //============================ 6th ring : theta = 60 ~ 68 degree ==============================================

  n=5;

	G4Box *Neutron_Solid_6 = new G4Box("Neutron_Solid_6", NeutronSize/2, NeutronSize/2, NeutronThick_Z/2);
	G4LogicalVolume *Neutron_LV_6 = new G4LogicalVolume(Neutron_Solid_6, Scint, "Neutron_LV_6");
	
	for(G4int copyN=0; copyN<NeutronPhiNum[n]; copyN++)
	{
		new G4PVPlacement(new G4RotationMatrix((360/NeutronPhiNum[n])*copyN*CLHEP::degree, Angle[n], 0*CLHEP::degree), G4ThreeVector(0*CLHEP::mm, 0*CLHEP::mm, TargetRadiusDistance+NeutronThick_Z/2).rotateX(Angle[n]).rotateZ((360/NeutronPhiNum[n])*copyN*CLHEP::degree), "Neutron_PV", Neutron_LV_6, labPV, false, 61+copyN);
	}

	Neutron_LV_6 -> SetVisAttributes(Neutron_VisAttrib);
    
  



  //============================ 7th ring : theta = 68 ~ 76 degree ==============================================

  n=6;

	G4Box *Neutron_Solid_7 = new G4Box("Neutron_Solid_7", NeutronSize/2, NeutronSize/2, NeutronThick_Z/2);
	G4LogicalVolume *Neutron_LV_7 = new G4LogicalVolume(Neutron_Solid_7, Scint, "Neutron_LV_7");
	
	for(G4int copyN=0; copyN<NeutronPhiNum[n]; copyN++)
	{
		new G4PVPlacement(new G4RotationMatrix((360/NeutronPhiNum[n])*copyN*CLHEP::degree, Angle[n], 0*CLHEP::degree), G4ThreeVector(0*CLHEP::mm, 0*CLHEP::mm, TargetRadiusDistance+NeutronThick_Z/2).rotateX(Angle[n]).rotateZ((360/NeutronPhiNum[n])*copyN*CLHEP::degree), "Neutron_PV", Neutron_LV_7, labPV, false, 71+copyN);
	}

	Neutron_LV_7 -> SetVisAttributes(Neutron_VisAttrib);
    
  



  //============================ 8th ring : theta = 76 ~ 84 CLHEP::degree ==============================================

  n=7;

	G4Box *Neutron_Solid_8 = new G4Box("Neutron_Solid_8", NeutronSize/2, NeutronSize/2, NeutronThick_Z/2);
	G4LogicalVolume *Neutron_LV_8 = new G4LogicalVolume(Neutron_Solid_8, Scint, "Neutron_LV_8");
	
	for(G4int copyN=0; copyN<NeutronPhiNum[n]; copyN++)
	{
		new G4PVPlacement(new G4RotationMatrix((360/NeutronPhiNum[n])*copyN*CLHEP::degree, Angle[n], 0*CLHEP::degree), G4ThreeVector(0*CLHEP::mm, 0*CLHEP::mm, TargetRadiusDistance+NeutronThick_Z/2).rotateX(Angle[n]).rotateZ((360/NeutronPhiNum[n])*copyN*CLHEP::degree), "Neutron_PV", Neutron_LV_8, labPV, false, 81+copyN);
	}

	Neutron_LV_8 -> SetVisAttributes(Neutron_VisAttrib);
    


  

  //============================ 9th ring : theta = 84 ~ 92 CLHEP::degree ==============================================

  n=8;

	G4Box *Neutron_Solid_9 = new G4Box("Neutron_Solid_9", NeutronSize/2, NeutronSize/2, NeutronThick_Z/2);
	G4LogicalVolume *Neutron_LV_9 = new G4LogicalVolume(Neutron_Solid_9, Scint, "Neutron_LV_9");
	
	for(G4int copyN=0; copyN<NeutronPhiNum[n]; copyN++)
	{
		new G4PVPlacement(new G4RotationMatrix((360/NeutronPhiNum[n])*copyN*degree, Angle[n], 0*CLHEP::degree), G4ThreeVector(0*CLHEP::mm, 0*CLHEP::mm, TargetRadiusDistance+NeutronThick_Z/2).rotateX(Angle[n]).rotateZ((360/NeutronPhiNum[n])*copyN*CLHEP::degree), "Neutron_PV", Neutron_LV_9, labPV, false, 91+copyN);
	}

	Neutron_LV_9 -> SetVisAttributes(Neutron_VisAttrib);
    


  

  //============================ 10th ring : theta = 92 ~ 100 degree ==============================================

  n=9;

	G4Box *Neutron_Solid_10 = new G4Box("Neutron_Solid_10", NeutronSize/2, NeutronSize/2, NeutronThick_Z/2);
	G4LogicalVolume *Neutron_LV_10 = new G4LogicalVolume(Neutron_Solid_10, Scint, "Neutron_LV_10");
	
	for(G4int copyN=0; copyN<NeutronPhiNum[n]; copyN++)
	{
		new G4PVPlacement(new G4RotationMatrix((360/NeutronPhiNum[n])*copyN*CLHEP::degree, Angle[n], 0*CLHEP::degree), G4ThreeVector(0*CLHEP::mm, 0*CLHEP::CLHEP::mm, TargetRadiusDistance+NeutronThick_Z/2).rotateX(Angle[n]).rotateZ((360/NeutronPhiNum[n])*copyN*CLHEP::degree), "Neutron_PV", Neutron_LV_10, labPV, false, 101+copyN);
	}

	Neutron_LV_10 -> SetVisAttributes(Neutron_VisAttrib);
    
  
    
  
  NeutronSD* NSD = new NeutronSD("/Neutron/NeutronSD","Neutron");
  sdManager -> AddNewDetector(NSD);
  Neutron_LV_1->SetSensitiveDetector(NSD);
  Neutron_LV_2->SetSensitiveDetector(NSD);
  Neutron_LV_3->SetSensitiveDetector(NSD);
  Neutron_LV_4->SetSensitiveDetector(NSD);
  Neutron_LV_5->SetSensitiveDetector(NSD);
  Neutron_LV_6->SetSensitiveDetector(NSD);
  Neutron_LV_7->SetSensitiveDetector(NSD);
  Neutron_LV_8->SetSensitiveDetector(NSD);
  Neutron_LV_9->SetSensitiveDetector(NSD);
  Neutron_LV_10->SetSensitiveDetector(NSD);
*/
}

//Neutron Block Test Detector
/*
void lampslowDetectorConstruction::ConstructBlock(G4VPhysicalVolume* labPV, G4SDManager* sdManager)
{
  G4VSolid* BlockSolid = new G4Box("Block", block_dxy/2, block_dxy/2, block_dz/2);
  G4LogicalVolume* BlockLV = new G4LogicalVolume(BlockSolid, Scint, "BlockLV");
  new G4PVPlacement(0, G4ThreeVector(0,0,block_z+block_dz/2), block_PVName, BlockLV, labPV, false, 0);

  DetSD* blockSD = new DetSD("/block", block_PVName, 1, 0, 0);
  sdManager -> AddNewDetector(blockSD);
  BlockLV -> SetSensitiveDetector(blockSD);
}

*/
//Materials
void lampslowDetectorConstruction::ConstructMaterials()
{
  const G4double labTemp = STP_Temperature + 20.*kelvin;

  G4double atomicNumber = 1.;
  G4double massOfMole = 1.008*g/mole;
  G4double densityUniverse = 1.e-25*g/CLHEP::cm3;
  G4double tempUniverse= 2.73 *kelvin;
  G4double pressureUniverse = 3.e-18*pascal;

  G4double density_PET = 1.4*g/CLHEP::cm3;

  // Elements
  elN  = new G4Element("Nitrogen", "N",  7,  14.000674*g/mole);
  elO  = new G4Element("Oxygen",   "O",  8,  15.9994*g/mole);
  elC  = new G4Element("Carbon",   "C",  6,  12.011*g/mole);
  elH  = new G4Element("Hydrogen", "H",  1,  1.00794*g/mole);
  elAr = new G4Element("Argon",    "Ar", 18, 39.938*g/mole);
  elSi = new G4Element("Silicon",    "Si", 14, 28.0855*g/mole);
  elCs = new G4Element("Cesium",    "Cs", 55, 132.90543*g/mole);
  elI = new G4Element("Iodine",    "I", 53, 126.90447*g/mole);

  // Materials
  Air = new G4Material("Air", 1.2929e-03*g/CLHEP::cm3, 3, kStateGas, labTemp);
  Air -> AddElement(elN, 75.47/99.95);
  Air -> AddElement(elO, 23.20/99.95);
  Air -> AddElement(elAr, 1.28/99.95);

  //Vacuum
  G4double temperature = tempUniverse;
  G4double pressure = pressureUniverse;
  G4double density = densityUniverse;
   
  Vacuum =  new G4Material("Vacuum", density, 2, kStateGas, temperature);
  Vacuum -> AddElement(elO,.3);
  Vacuum -> AddElement(elN,.7);


  PET = new G4Material("PET", density_PET, 3, kStateSolid);
  PET -> AddElement(elC,10);
  PET -> AddElement(elH,8);
  PET -> AddElement(elO,4);

  Scint = new G4Material("Scintillator", 1.05*g/cm3, 2, kStateSolid, labTemp);
  Scint -> AddElement(elC, 10);
  Scint -> AddElement(elH, 11);

  //SI
  Silicon = new G4Material("Silicon", 2.33*g/cm3, 1, kStateSolid, labTemp);
  Silicon -> AddElement(elSi, 1);

  //CsI
  CsI = new G4Material("CsI", 4.51*g/cm3, 2, kStateSolid, labTemp);
  CsI -> AddElement(elCs, 1);
  CsI -> AddElement(elI, 1);
 // CsI -> GetIonization() -> SetMeanExcitationEnergy(533.1*eV);


}

void lampslowDetectorConstruction::DestructMaterials()
{
  delete Air;
  delete Vacuum;

  delete elH;
  delete elC;
  delete elAr;
  delete elO;
  delete elN;
  delete Silicon;
  delete CsI;
}

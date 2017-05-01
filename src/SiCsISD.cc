//Geant4 headers
#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"

//lampslow headers
#include "SiCsISD.hh"
#include "SiCsIHit.hh"
#include "lampslowData.hh"

//C++ headers
#include <iostream>
#include <sstream>

//Root Headers
#include "TMath.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TStyle.h"

using namespace std;


TH1D *hist = new TH1D("hist_scattering","#theta_{scatt}; #theta (degree);Counts",100,0,180);
TH2D *hist2D = new TH2D("EnergyScatt","Energy vs #theta_{scatt}; #theta (degree);Energy",300,0,180,150,0,13);
TH1D *hist1 = new TH1D("hist_scattering1","#theta_{scatt}; #theta (degree);Counts",100,0,18000);
//TCanvas *cvs = new TCanvas("scattering","",700,700);
//TCanvas *cvs2D = new TCanvas("EnergyScatt","",700,700);


SiCsISD::SiCsISD(const G4String &name, const G4String pv_name)
:G4VSensitiveDetector(name), PVName(pv_name)
{
  //for unique name ...
  /*
  G4int fullID = detTypeID*10000+sectID*100+detID;
  std::ostringstream ss_fullID;
  ss_fullID << fullID;
  G4String str_collectionName = "DetHitsColletion" + G4String(ss_fullID.str());
  collectionName.insert(str_collectionName);
  */
  G4String hitCollectionName = "SiCsIHitCollection_" + PVName; 
  collectionName.insert(hitCollectionName);
}

SiCsISD::~SiCsISD() {}

void SiCsISD::Initialize(G4HCofThisEvent *HCTE)
{
  hitsCollection = new G4THitsCollection<SiCsIHit>(SensitiveDetectorName, collectionName[0]);
  G4int hcid = G4SDManager::GetSDMpointer() -> GetCollectionID(collectionName[0]);
  HCTE -> AddHitsCollection(hcid, hitsCollection);
}
  
G4bool SiCsISD::ProcessHits(G4Step *aStep, G4TouchableHistory *ROHist)
{
  
  gStyle -> SetOptStat(0);
  gStyle -> SetTitleYOffset(1.85);
  gStyle -> SetTitleXOffset(1.45);
  gStyle -> SetTitleSize(0.04,"xy");
  gStyle -> SetPadLeftMargin(0.16);
  gStyle -> SetPadRightMargin(0.05);
  gStyle -> SetPadTopMargin(0.12);
  gStyle -> SetPadBottomMargin(0.15);
  
//  if(!ROHist) return false;

  energyDeposit = aStep -> GetTotalEnergyDeposit();
//  if(energyDeposit == 0.) return false;
 /* 
  if((aStep -> GetPreStepPoint() -> GetPhysicalVolume() -> GetName() == "Si_PV") &&(energyDeposit!=0.))
  {
    parentID    = aStep -> GetTrack() -> GetParentID();
    pID         = aStep -> GetTrack() -> GetDefinition() -> GetPDGEncoding();
    prePosition = aStep -> GetPreStepPoint() -> GetPosition();
    preMomentum = aStep -> GetPreStepPoint() -> GetMomentum();
    time        = aStep -> GetPreStepPoint() -> GetGlobalTime();
    copyNum     = aStep -> GetPreStepPoint() -> GetTouchableHandle() -> GetCopyNumber();
    detID       = 1;   

    SiCsIHit *aHit = new SiCsIHit(parentID, pID, prePosition, preMomentum, time, energyDeposit,
                              detTypeID, sectID, copyNum, detID);
    hitsCollection -> insert(aHit);
  }
  
  if((aStep -> GetPreStepPoint() -> GetPhysicalVolume() -> GetName() == "CsI_PV") &&(energyDeposit!=0.))
  {
    parentID    = aStep -> GetTrack() -> GetParentID();
    pID         = aStep -> GetTrack() -> GetDefinition() -> GetPDGEncoding();
    prePosition = aStep -> GetPreStepPoint() -> GetPosition();
    preMomentum = aStep -> GetPreStepPoint() -> GetMomentum();
    time        = aStep -> GetPreStepPoint() -> GetGlobalTime();
    copyNum     = aStep -> GetPreStepPoint() -> GetTouchableHandle() -> GetCopyNumber();
    detID       = 2;

    SiCsIHit *aHit = new SiCsIHit(parentID, pID, prePosition, preMomentum, time, energyDeposit,
                              detTypeID, sectID, copyNum, detID);
    hitsCollection -> insert(aHit);
    }
*/  
  if((aStep -> GetPreStepPoint() -> GetPhysicalVolume() -> GetName() == "TestVetoBlock_PV") &&(energyDeposit!=0.))
  {
    parentID    = aStep -> GetTrack() -> GetParentID();
    pID         = aStep -> GetTrack() -> GetDefinition() -> GetPDGEncoding();
    prePosition = aStep -> GetPreStepPoint() -> GetPosition();
    preMomentum = aStep -> GetPreStepPoint() -> GetMomentum();
    time        = aStep -> GetPreStepPoint() -> GetGlobalTime();
    copyNum     = aStep -> GetPreStepPoint() -> GetTouchableHandle() -> GetCopyNumber();
    detID       = 5;

    SiCsIHit *aHit = new SiCsIHit(parentID, pID, prePosition, preMomentum, time, energyDeposit,
                              detTypeID, sectID, copyNum, detID);
    hitsCollection -> insert(aHit);
    }

  if((aStep -> GetPreStepPoint() -> GetPhysicalVolume() -> GetName() == "TestSiBlock_PV_1") &&(energyDeposit!=0.))
  {
    parentID    = aStep -> GetTrack() -> GetParentID();
    pID         = aStep -> GetTrack() -> GetDefinition() -> GetPDGEncoding();
    prePosition = aStep -> GetPreStepPoint() -> GetPosition();
    preMomentum = aStep -> GetPreStepPoint() -> GetMomentum();
    time        = aStep -> GetPreStepPoint() -> GetGlobalTime();
    copyNum     = aStep -> GetPreStepPoint() -> GetTouchableHandle() -> GetCopyNumber();
    detID       = 7;

    SiCsIHit *aHit = new SiCsIHit(parentID, pID, prePosition, preMomentum, time, energyDeposit,
                              detTypeID, sectID, copyNum, detID);
    hitsCollection -> insert(aHit);
    }
/*
  if((aStep -> GetPreStepPoint() -> GetPhysicalVolume() -> GetName() == "TestSiBlock_PV_2") &&(energyDeposit!=0.))
  {
    parentID    = aStep -> GetTrack() -> GetParentID();
    pID         = aStep -> GetTrack() -> GetDefinition() -> GetPDGEncoding();
    prePosition = aStep -> GetPreStepPoint() -> GetPosition();
    preMomentum = aStep -> GetPreStepPoint() -> GetMomentum();
    time        = aStep -> GetPreStepPoint() -> GetGlobalTime();
    copyNum     = aStep -> GetPreStepPoint() -> GetTouchableHandle() -> GetCopyNumber();
    detID       = 8;

    SiCsIHit *aHit = new SiCsIHit(parentID, pID, prePosition, preMomentum, time, energyDeposit,
                              detTypeID, sectID, copyNum, detID);
    hitsCollection -> insert(aHit);
    }
 
  if((aStep -> GetPreStepPoint() -> GetPhysicalVolume() -> GetName() == "TestSiBlock_PV_3") &&(energyDeposit!=0.))
  {
    parentID    = aStep -> GetTrack() -> GetParentID();
    pID         = aStep -> GetTrack() -> GetDefinition() -> GetPDGEncoding();
    prePosition = aStep -> GetPreStepPoint() -> GetPosition();
    preMomentum = aStep -> GetPreStepPoint() -> GetMomentum();
    time        = aStep -> GetPreStepPoint() -> GetGlobalTime();
    copyNum     = aStep -> GetPreStepPoint() -> GetTouchableHandle() -> GetCopyNumber();
    detID       = 9;

    SiCsIHit *aHit = new SiCsIHit(parentID, pID, prePosition, preMomentum, time, energyDeposit,
                              detTypeID, sectID, copyNum, detID);
    hitsCollection -> insert(aHit);
    }
  */
  if((aStep -> GetPreStepPoint() -> GetPhysicalVolume() -> GetName() == "TestCsIBlock_PV") &&(energyDeposit!=0.))
  {
    parentID    = aStep -> GetTrack() -> GetParentID();
    pID         = aStep -> GetTrack() -> GetDefinition() -> GetPDGEncoding();
    prePosition = aStep -> GetPreStepPoint() -> GetPosition();
    postPosition = aStep -> GetPostStepPoint() -> GetPosition();
    preMomentum = aStep -> GetPreStepPoint() -> GetMomentum();
    time        = aStep -> GetPreStepPoint() -> GetGlobalTime();
    copyNum     = aStep -> GetPreStepPoint() -> GetTouchableHandle() -> GetCopyNumber();
    detID       = 11;

    SiCsIHit *aHit = new SiCsIHit(parentID, pID, prePosition, preMomentum, time, energyDeposit,
                              detTypeID, sectID, copyNum, detID);
    hitsCollection -> insert(aHit);
/*
    cout << "*** HitCsI ***" << endl;
    cout << "pID : " << pID << endl;
    cout << "energy dep : " << energyDeposit << endl;
/    cout << endl;
    cout << endl;
    cout << endl;
  
    if(energyDeposit==0)
    {
    hist1->Fill(pID);
    }
    cvs->cd();
    hist1->Draw();
    cvs->Write();*/
   }
  
 /* 
  if((aStep -> GetPreStepPoint() -> GetPhysicalVolume() -> GetName() == "TestCsIBlock_PV"))
  {
    parentID    = aStep -> GetTrack() -> GetParentID();
    pID         = aStep -> GetTrack() -> GetDefinition() -> GetPDGEncoding();
    prePosition = aStep -> GetPreStepPoint() -> GetPosition();
    postPosition = aStep -> GetPostStepPoint() -> GetPosition();
    preMomentum = aStep -> GetPreStepPoint() -> GetMomentum();
    time        = aStep -> GetPreStepPoint() -> GetGlobalTime();
    copyNum     = aStep -> GetPreStepPoint() -> GetTouchableHandle() -> GetCopyNumber();
    detID       = 8;

    SiCsIHit *aHit = new SiCsIHit(parentID, pID, prePosition, postPosition, preMomentum, time, energyDeposit,
                              detTypeID, sectID, copyNum, detID);
    hitsCollection -> insert(aHit);

    cout << "*** HitCsI ***" << endl;
    cout << "pID : " << pID << endl;
    cout << "energy dep : " << energyDeposit << endl;
    cout << endl;
    cout << endl;
    cout << endl;
  
   }*/
  /*
  if((aStep -> GetPreStepPoint() -> GetPhysicalVolume() -> GetName() == "TestCsIBlock_PV") &&(aStep -> GetPostStepPoint() -> GetPhysicalVolume() -> GetName() == "Track_PV"))
  {
    parentID    = aStep -> GetTrack() -> GetParentID();
    pID         = aStep -> GetTrack() -> GetDefinition() -> GetPDGEncoding();
    prePosition = aStep -> GetPreStepPoint() -> GetPosition();
    postPosition = aStep -> GetPostStepPoint() -> GetPosition();
   // postPosition = prePosition;
    preMomentum = aStep -> GetPreStepPoint() -> GetMomentum();
    time        = aStep -> GetPreStepPoint() -> GetGlobalTime();
    copyNum     = aStep -> GetPreStepPoint() -> GetTouchableHandle() -> GetCopyNumber();
    detID       = 0;
    G4int sectionID = 100;

    SiCsIHit *aHit = new SiCsIHit(parentID, pID, prePosition, postPosition, preMomentum, time, energyDeposit,
                              detTypeID, sectID, copyNum, detID);
    hitsCollection -> insert(aHit);

 
    G4double Length_PreT = TMath::Sqrt(prePosition(0)*prePosition(0)+prePosition(1)*prePosition(1));
    G4double Length_PostT = TMath::Sqrt(postPosition(0)*postPosition(0)+postPosition(1)*postPosition(1));
    G4double Z_diff = postPosition(2)-prePosition(2);
    G4double Trans_diff = Length_PostT-Length_PreT;
    G4double Tangent_Theta=0;
    if(Z_diff>0) Tangent_Theta = TMath::ATan(TMath::Abs(Trans_diff/Z_diff));
    else if(Z_diff<0) Tangent_Theta = TMath::Pi()-TMath::ATan(TMath::Abs(Trans_diff/Z_diff));
    else if(Z_diff==0) Tangent_Theta = TMath::Pi()/2;

    G4double ScatteringAngle = Tangent_Theta*180/TMath::Pi();

    cout << "*** Hit ***" << endl;
    cout << "pID : " << pID << endl;
    cout << "Pre Positiion : " << prePosition << endl;
    cout << "Pos Positiion : " << postPosition << endl;
    cout << "Z Difference : " << Z_diff  << endl;
    cout << "T Difference : " << Trans_diff  << endl;
    cout << "Angle : " << ScatteringAngle << endl;
    cout << "Time : " << time << endl;
    cout << "energy dep : " << energyDeposit << endl;
    cout << endl;
    cout << endl;
    cout << endl;

    G4double pTot = TMath::Sqrt(preMomentum(0)*preMomentum(0)+preMomentum(1)*preMomentum(1)+preMomentum(2)*preMomentum(2));
    G4double pZ = preMomentum(2);
    G4double pX = preMomentum(0);
    G4double pY = preMomentum(1);
    G4double pT = TMath::Sqrt(pX*pX+pY*pY);
    G4double Angle = TMath::ASin(pT/pTot);
    G4double ScatteringAngle=0;
    if(pZ>=0) ScatteringAngle = Angle*180/TMath::Pi();
    else if(pZ<0) ScatteringAngle = 180 - Angle*180/TMath::Pi();

    G4double mass = 939.5;
    G4double Etot = TMath::Sqrt(pTot*pTot+mass*mass);
    G4double Ek = Etot-mass;



    G4double Pos = prePosition(2);
    cout << "????" << endl;
    cout << "ps : " << Pos << endl;
    cout << endl;
    cout << endl;
    cout << endl;

    if(energyDeposit!=0) 
    {
      cout << "### Energy Deposit exist!! ###" << endl;
      cout << endl;
      cout << "eDep " << ": " << energyDeposit <<endl;
      cout << "Scatt Angle : " << ScatteringAngle << endl;
      cout << "pID : " << pID << endl;
      cout << endl;
      cout << endl;
      cout << endl;
    }
  hist -> SetMinimum(0.1001);
  hist -> GetXaxis() -> CenterTitle();
  hist -> GetYaxis() -> CenterTitle();
  if(ScatteringAngle>10)
  {
   hist -> Fill(ScatteringAngle);
  }
  hist2D -> GetXaxis() -> CenterTitle();
  hist2D -> GetYaxis() -> CenterTitle();
  hist2D -> Fill(ScatteringAngle,Ek);
// if(pID!=2112)  hist1->Fill(Z_diff);
  cvs->cd();
  cvs->SetLogy();
  hist->Draw(); 
//  hist1->Draw();
  cvs->Write();
  cvs2D->cd();
  hist2D->Draw("colz");
  cvs2D->Write();
  }*/

  if((aStep -> GetPreStepPoint() -> GetPhysicalVolume() -> GetName() == "Track_PV"))
  {
    parentID    = aStep -> GetTrack() -> GetParentID();
    pID         = aStep -> GetTrack() -> GetDefinition() -> GetPDGEncoding();
    prePosition = aStep -> GetPreStepPoint() -> GetPosition();
    postPosition = aStep -> GetPostStepPoint() -> GetPosition();
    preMomentum = aStep -> GetPreStepPoint() -> GetMomentum();
    time        = aStep -> GetPreStepPoint() -> GetGlobalTime();
    copyNum     = aStep -> GetPreStepPoint() -> GetTouchableHandle() -> GetCopyNumber();
    detID       = 0;
    G4int sectionID = 100;

    SiCsIHit *aHit = new SiCsIHit(parentID, pID, prePosition, preMomentum, time, energyDeposit,
                              detTypeID, sectID, copyNum, detID);
    hitsCollection -> insert(aHit);
    } 
  return true;
}

void SiCsISD::EndOfEvent(G4HCofThisEvent *HCTE)
{
}

///____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Afroditi Papadopoulou <apapadop \at mit.edu>
         Massachusetts Institute of Technology - October 04, 2019

 @ October 4, 2019 - Afroditi Papadopoulou (AP)
   Created this new module that controls the addition of the recoil nucleon in the event record 
   and extracts its kinematics 
*/
//____________________________________________________________________________

#include <cstdlib>

#include <TLorentzVector.h>
#include <TVector3.h>
#include <TParticlePDG.h>
#include <TMath.h>

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Physics/NuclearState/SRCNuclearRecoil.h"

#include "Physics/NuclearState/NuclearModel.h"
#include "Physics/NuclearState/NuclearModelI.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/NuclearState/FermiMomentumTablePool.h"
#include "Physics/NuclearState/FermiMomentumTable.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Utils/KineUtils.h"
#include "Physics/NuclearState/NuclearUtils.h"

using namespace genie;
using namespace genie::constants;

//___________________________________________________________________________
SRCNuclearRecoil::SRCNuclearRecoil() :
SecondNucleonEmissionI("genie::SRCNuclearRecoil")
{

}
//___________________________________________________________________________
SRCNuclearRecoil::SRCNuclearRecoil(string config) :
  SecondNucleonEmissionI("genie::SRCNuclearRecoil", config )
{

}

//___________________________________________________________________________

SRCNuclearRecoil::~SRCNuclearRecoil()
{

}

//___________________________________________________________________________

void SRCNuclearRecoil::ProcessEventRecord(GHepRecord * evrec) const
{

  Interaction *  interaction = evrec       -> Summary();
  InitialState * init_state  = interaction -> InitStatePtr();
  Target *       tgt         = init_state  -> TgtPtr();

  // do nothing for non-nuclear targets
  if(!tgt->IsNucleus()) return;

  // access the hit nucleon and target nucleus at the GHEP record
  GHepParticle * nucleon = evrec->HitNucleon();
  GHepParticle * nucleus = evrec->TargetNucleus();
  assert(nucleon);
  assert(nucleus);

  TVector3 p3 = fNuclModel->Momentum3();
  double pF2 = p3.Mag2(); // (fermi momentum)^2

  // Set this to either a proton or neutron to eject a secondary particle
  int eject_nucleon_pdg = this->SRCRecoilPDG(nucleon, nucleus, tgt, pF2);

  // Ejection of secondary particle
  if (eject_nucleon_pdg != 0) { 
	if (fGaussianEmission) { EmitSecondNucleon(evrec,eject_nucleon_pdg); } 
        else { SecondNucleonEmissionI::EmitSecondNucleon(evrec,eject_nucleon_pdg); }
  }

}

//___________________________________________________________________________

bool SRCNuclearRecoil::EmitSecondNucleon(GHepRecord * evrec, const int eject_nucleon_pdg) const {

  LOG("SRCNuclearRecoil", pINFO) << "Adding a recoil nucleon with PDG " << eject_nucleon_pdg ;

  GHepParticle * nucleon = evrec->HitNucleon();

  GHepStatus_t status = kIStHadronInTheNucleus;
  int imom = evrec->TargetNucleusPosition();

  //-- Has opposite momentum from the struck nucleon
  double vx = nucleon->Vx();
  double vy = nucleon->Vy();
  double vz = nucleon->Vz();
  // recoil nucleon not exactly in the opposite direction of the hit nucleon 
  // Use gaussian distribution for center-of-mass motion
  double px = gRandom->Gaus(0,fGaussianSigma) - nucleon->Px();
  double py = gRandom->Gaus(0,fGaussianSigma) - nucleon->Py();
  double pz = gRandom->Gaus(0,fGaussianSigma) - nucleon->Pz();

  double M  = PDGLibrary::Instance()->Find(eject_nucleon_pdg)->Mass();
  double E  = TMath::Sqrt(px*px+py*py+pz*pz+M*M);

  evrec->AddParticle( eject_nucleon_pdg, status, imom, -1, -1, -1, px, py, pz, E, vx, vy, vz, 0 );

  return true ;

}
//___________________________________________________________________________

int SRCNuclearRecoil::SRCRecoilPDG(GHepParticle * nucleon, GHepParticle * nucleus, Target* tgt, double pF2) const {

      int eject_nucleon_pdg = 0;

      const int nucleus_pdgc = nucleus->Pdg();
      const int nucleon_pdgc = nucleon->Pdg();

      // Calculate the Fermi momentum, using a local Fermi gas if the
      // nuclear model is LocalFGM, and RFG otherwise
      double kF;
      if(fNuclModel->ModelType(*tgt) == kNucmLocalFermiGas){
	assert(pdg::IsProton(nucleon_pdgc) || pdg::IsNeutron(nucleon_pdgc));
	int A = tgt->A();
	bool is_p = pdg::IsProton(nucleon_pdgc);
	double numNuc = (is_p) ? (double)tgt->Z():(double)tgt->N();
	double radius = nucleon->X4()->Vect().Mag();
	double hbarc = kLightSpeed*kPlankConstant/genie::units::fermi;
	kF= TMath::Power(3*kPi2*numNuc*
		  genie::utils::nuclear::Density(radius,A),1.0/3.0) *hbarc;
      }else{
        kF = fKFTable->FindClosestKF(nucleus_pdgc, nucleon_pdgc);
      }
      if (TMath::Sqrt(pF2) > kF) {

        double Pp = (nucleon->Pdg() == kPdgProton) ? fPPPairPercentage : fPNPairPercentage;
        RandomGen * rnd = RandomGen::Instance();
        double prob = rnd->RndGen().Rndm();
        eject_nucleon_pdg = (prob > Pp) ? kPdgNeutron : kPdgProton;

      }

      return eject_nucleon_pdg;
}
//___________________________________________________________________________
void SRCNuclearRecoil::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void SRCNuclearRecoil::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void SRCNuclearRecoil::LoadConfig(void)
{

  SecondNucleonEmissionI::LoadConfig() ;

  this->GetParamDef("PNPairPercentage",       fPNPairPercentage,    0.95);

  if (fPNPairPercentage < 0. || fPNPairPercentage > 1.) { 

	LOG("SRCNuclearRecoil", pFATAL)
	<< "PNPairPercentage either less than 0 or greater than 1: Exiting" ;

	exit(78); 
  }

  fPPPairPercentage = 1. - fPNPairPercentage;

  // Get the Fermi momentum table for relativistic Fermi gas
  GetParam( "FermiMomentumTable", fKFTableName ) ;
  fKFTable = 0;
  FermiMomentumTablePool * kftp = FermiMomentumTablePool::Instance();
  fKFTable = kftp->GetTable(fKFTableName);
  assert(fKFTable);

  this->GetParam("SRC-GaussianEmission",fGaussianEmission);
  this->GetParamDef("SRC-GaussianSigma",fGaussianSigma,0.);

}
//____________________________________________________________________________

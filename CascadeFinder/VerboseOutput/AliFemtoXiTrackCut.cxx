#include "AliFemtoXiTrackCut.h"
#include "AliESDtrack.h"
#include <cstdio>
#include <bitset>

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoXiTrackCut);
  /// \endcond
#endif



AliFemtoXiTrackCut::AliFemtoXiTrackCut() : AliFemtoV0TrackCut(), fMaxEtaXi(100), fMinPtXi(0), fMaxPtXi(100), fChargeXi(1.0), fMaxEtaBac(100), fMinPtBac(0), fMaxPtBac(100), fTPCNclsBac(0), fNdofBac(100), fStatusBac(0), fMaxDcaXi(0), fMinDcaXiBac(0), fMaxDcaXiDaughters(0), fMinCosPointingAngleXi(0), fMaxDecayLengthXi(100.0), fInvMassXiMin(0), fInvMassXiMax(1000), fFailedXi1(0), fFailedXi2(0), fFailedXi4(0), fFailedXi5(0)
{
  // Default constructor
  fTotal = 0;
  fTotalFailedInvMass = 0;
  fTotalPassedInvMass = 0;
  for(int i=0; i<20; i++)
  {
    fFailedTracker[i] = 0;
    fTrueFailedTracker[i] = 0;
  }
}
 //------------------------------
AliFemtoXiTrackCut::~AliFemtoXiTrackCut(){
  /* noop */
}
//------------------------------
bool AliFemtoXiTrackCut::Pass(const AliFemtoXi* aXi)
{
  //cout << endl << "---------- ENTERING AliFemtoXiTrackCut::Pass!!!!!!" << endl;
  // test the particle and return 
  // true if it meets all the criteria
  // false if it doesn't meet at least one of the criteria

  fTotal++;

   //invariant mass Xi
  if(aXi->MassXi()<fInvMassXiMin || aXi->MassXi()>fInvMassXiMax)
     {
       fTotalFailedInvMass++;
//       cout << "fTotalFailedInvMass = " << fTotalFailedInvMass << endl;
       return false;
     }
   
  fTotalPassedInvMass++;
//  cout << "fTotalPassedInvMass = " << fTotalPassedInvMass << endl;
//  cout << endl << "PASSED invariant mass Xi" << endl;

  cout << "fTotal = " << fTotal << endl;
  cout << "\tfTotalFailedInvMass = " << fTotalFailedInvMass << endl;
  cout << "\tfTotalPassedInvMass = " << fTotalPassedInvMass << endl;


  for(int i=0; i<20; i++)
  {
    cout << "fFailedTracker[" << i << "] = " << fFailedTracker[i] << endl;
  }

  for(int i=0; i<20; i++)
  {
    cout << "fTrueFailedTracker[" << i << "] = " << fTrueFailedTracker[i] << endl;
  }


  Float_t pt = aXi->PtXi();
  Float_t eta = aXi->EtaXi();


//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  if(aXi->ChargeXi()==0)
  {
    fTrueFailedTracker[0]++;
  }

  //ParticleType selection
  //If fParticleTypeXi=kAll, any charged candidate will pass 
  if(fParticleTypeXi == kXiPlus && aXi->ChargeXi() == -1) 
  {
    return false;
  }

  if(fParticleTypeXi == kXiMinus && aXi->ChargeXi() == 1) 
  {
    return false;
  }


  //kinematic cuts
  if(TMath::Abs(eta) > fMaxEtaXi) 
  {
    fTrueFailedTracker[1]++;
  }

  if(pt < fMinPtXi)
  {
    fTrueFailedTracker[2]++;
  }

  if(pt > fMaxPtXi)
  {
    fTrueFailedTracker[3]++;
  }

  if(TMath::Abs(aXi->EtaBac()) > fMaxEtaBac)
  {
    fTrueFailedTracker[4]++;
  }

  if(aXi->PtBac()< fMinPtBac)
  {
    fTrueFailedTracker[5]++;
  }

  if(aXi->PtBac()> fMaxPtBac)
  {
    fTrueFailedTracker[6]++;
  }

    //quality cuts
    if(aXi->StatusBac() == 999)
    {
      fTrueFailedTracker[8]++;
    }

    if(aXi->TPCNclsBac()<fTPCNclsBac)
    {
      fTrueFailedTracker[9]++;
    }

    if(aXi->NdofBac()>fNdofBac)
    {
      fTrueFailedTracker[10]++;
    }


//WHAT IS THE POINT OF THIS CUT?!!?!?!?!?!
    if(!(aXi->StatusBac()&fStatusBac))
    {
      fTrueFailedTracker[11]++;
    }

    //DCA Xi to prim vertex
    if(TMath::Abs(aXi->DcaXiToPrimVertex())>fMaxDcaXi)
    {
      fTrueFailedTracker[12]++;
    }


    //DCA Xi bachelor to prim vertex
    if(TMath::Abs(aXi->DcaBacToPrimVertex())<fMinDcaXiBac)
    {
      fTrueFailedTracker[13]++;
    }



    //DCA Xi daughters
    if(TMath::Abs(aXi->DcaXiDaughters())>fMaxDcaXiDaughters)
    {
      fTrueFailedTracker[14]++;
    }    

    //cos pointing angle
    if(aXi->CosPointingAngleXi()<fMinCosPointingAngleXi)
    {
      fTrueFailedTracker[15]++;
    }
    
    //decay length
    if(aXi->DecayLengthXi()>fMaxDecayLengthXi)
    {
      fTrueFailedTracker[16]++;
    }
    

  bool pid_check1=false;
  // Looking for Xi
  if (fParticleTypeXi == kXiMinus || fParticleTypeXi == kXiPlus) {
    if (IsPionNSigmaBac(aXi->PtBac(), aXi->BacNSigmaTPCPi(), aXi->BacNSigmaTOFPi())) //pion
	{
	  pid_check1=true;
	}

  }

  if (!pid_check1) 
  {
    fTrueFailedTracker[17]++;
  }  

  if(!AliFemtoV0TrackCut::Pass(aXi))
  {
    fTrueFailedTracker[18]++;
  }
 
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


  
  if(aXi->ChargeXi()==0)
  {
    fFailedTracker[0]++;
    return false;
  }

  //ParticleType selection
  //If fParticleTypeXi=kAll, any charged candidate will pass 
  if(fParticleTypeXi == kXiPlus && aXi->ChargeXi() == -1) 
  {
    return false;
  }

  if(fParticleTypeXi == kXiMinus && aXi->ChargeXi() == 1) 
  {
    return false;
  }

/*
  //kinematic cuts
  if(TMath::Abs(eta) > fMaxEtaXi) return false;  //put in kinematic cuts by hand
  if(pt < fMinPtXi) return false;
  if(pt > fMaxPtXi) return false;
  if(TMath::Abs(aXi->EtaBac()) > fMaxEtaBac) return false;
//  if(aXi->PtPos()< fMinPtBac) return false;  //WRONG!
  if(aXi->PtBac()< fMinPtBac) return false; 
//  if(aXi->PtPos()> fMaxPtBac) return false; //WRONG!
  if(aXi->PtBac()> fMaxPtBac) return false; 
*/

  //kinematic cuts
  if(TMath::Abs(eta) > fMaxEtaXi) 
  {
    fFailedTracker[1]++;

    fFailedXi1.push_back(TMath::Abs(eta));
    cout << "--------------------- fFailedXi1: TMath::Abs(eta) > fMaxEtaXi --------------------" << endl;
    for(unsigned int i=0; i<fFailedXi1.size(); i++)
    {
      cout << fFailedXi1[i] << endl;
    }

    return false;
  }

  if(pt < fMinPtXi)
  {
    fFailedTracker[2]++;

    fFailedXi2.push_back(pt);
    cout << "--------------------- fFailedXi2: pt<fMinPtXi --------------------" << endl;
    for(unsigned int i=0; i<fFailedXi2.size(); i++)
    {
      cout << fFailedXi2[i] << endl;
    }

    return false;
  }

  if(pt > fMaxPtXi)
  {
    fFailedTracker[3]++;
    return false;
  }

  if(TMath::Abs(aXi->EtaBac()) > fMaxEtaBac)
  {
    fFailedTracker[4]++;

    fFailedXi4.push_back(TMath::Abs(aXi->EtaBac()));
    cout << "--------------------- fFailedXi4: TMath::Abs(aXi->EtaBac()) > fMaxEtaBac --------------------" << endl;
    for(unsigned int i=0; i<fFailedXi4.size(); i++)
    {
      cout << fFailedXi4[i] << endl;
    }


    return false;
  }

  if(aXi->PtBac()< fMinPtBac)
  {
    fFailedTracker[5]++;

    fFailedXi5.push_back(aXi->PtBac());
    cout << "--------------------- fFailedXi5: aXi->PtBac()< fMinPtBac --------------------" << endl;
    for(unsigned int i=0; i<fFailedXi5.size(); i++)
    {
      cout << fFailedXi5[i] << endl;
    }

    return false;
  }

  if(aXi->PtBac()> fMaxPtBac)
  {
    fFailedTracker[6]++;
    return false;
  }


    //Xi from kinematics information
    //WHAT IS THE POINT OF THIS BLOCK OF CODE?!
    //The InvMass was already tested above, and I'm not sure why we would want to require BacNSigmaTPCPi to be zero
    if (fParticleTypeXi == kXiMinusMC || fParticleTypeXi == kXiPlusMC) {
      if(!(aXi->MassXi()>fInvMassXiMin && aXi->MassXi()<fInvMassXiMax) || !(aXi->BacNSigmaTPCPi()==0))  //WHY DO WE REQUIRE BacNSigmaTPCPi = 0?!?!?!?!?
        {
          fFailedTracker[7]++;
	  return false; 
        }

      else
	{
	  return true;  //WHY ARE WE RETURNING TRUE HERE?!?!?!?!
	} 

    }

/*
    //quality cuts
    if(aXi->StatusBac() == 999) return false;
    if(aXi->TPCNclsBac()<fTPCNclsBac) return false;
    if(aXi->NdofBac()>fNdofBac) return false;
    if(!(aXi->StatusBac()&fStatusBac)) return false;
*/
    //quality cuts
    if(aXi->StatusBac() == 999)
    {
      fFailedTracker[8]++;
      return false;
    }

    if(aXi->TPCNclsBac()<fTPCNclsBac)
    {
      fFailedTracker[9]++;
      return false;
    }

    if(aXi->NdofBac()>fNdofBac)
    {
      fFailedTracker[10]++;
      return false;
    }


//WHAT IS THE POINT OF THIS CUT?!!?!?!?!?!
    if(!(aXi->StatusBac()&fStatusBac))
    {
      fFailedTracker[11]++;
      return false;
    }

    //DCA Xi to prim vertex
    if(TMath::Abs(aXi->DcaXiToPrimVertex())>fMaxDcaXi)
    {
      fFailedTracker[12]++;
      return false;
    }


    //DCA Xi bachelor to prim vertex
    if(TMath::Abs(aXi->DcaBacToPrimVertex())<fMinDcaXiBac)
    {
      fFailedTracker[13]++;
      return false;
    }



    //DCA Xi daughters
    if(TMath::Abs(aXi->DcaXiDaughters())>fMaxDcaXiDaughters)
    {
      fFailedTracker[14]++;
      return false;
    }    

    //cos pointing angle
    if(aXi->CosPointingAngleXi()<fMinCosPointingAngleXi)
    {
      fFailedTracker[15]++;
      return false; 
    }
    
    //decay length
    if(aXi->DecayLengthXi()>fMaxDecayLengthXi)
    {
      fFailedTracker[16]++;
      return false;
    }
    
 
  if(fParticleTypeXi == kAll)
    return true;



  bool pid_check=false;
  // Looking for Xi
  if (fParticleTypeXi == kXiMinus || fParticleTypeXi == kXiPlus) {
    if (IsPionNSigmaBac(aXi->PtBac(), aXi->BacNSigmaTPCPi(), aXi->BacNSigmaTOFPi())) //pion
	{
	  pid_check=true;
	}

  }

  if (!pid_check) 
  {
    fFailedTracker[17]++;
    return false;
  }  

  if(!AliFemtoV0TrackCut::Pass(aXi))
  {
    fFailedTracker[18]++;
    return false;
  }

cout << endl << "\t\t\t\t\t pt = " << pt << endl;
cout << "\t\t\t\t\t aXi->StatusBac() = " << std::bitset<16>(aXi->StatusBac()) << endl;
cout << "\t\t\t\t\t fStatusBac       = " << std::bitset<16>(fStatusBac) << endl << endl;
  
  fFailedTracker[19]++;
  return true;
    
    
}
//------------------------------
AliFemtoString AliFemtoXiTrackCut::Report()
{
  // Prepare report from the execution
  string tStemp;
  char tCtemp[100];


  AliFemtoString returnThis = tStemp;
  return returnThis;
}
TList *AliFemtoXiTrackCut::ListSettings()
{
  // return a list of settings in a writable form
  TList *tListSetttings = new TList();
 
  return tListSetttings;
}


bool AliFemtoXiTrackCut::IsPionNSigmaBac(float mom, float nsigmaTPCPi, float nsigmaTOFPi)
{
  if(TMath::Abs(nsigmaTPCPi)<3.0) return true;


  return false;
}



void AliFemtoXiTrackCut::SetEtaXi(double x){
  fMaxEtaXi = x;
}

void AliFemtoXiTrackCut::SetPtXi(double min, double max){
  fMinPtXi = min;
  fMaxPtXi = max;
}

void AliFemtoXiTrackCut::SetChargeXi(int x){
  fChargeXi = x;
}

void AliFemtoXiTrackCut::SetMaxDecayLengthXi(double x){
  fMaxDecayLengthXi = x;
}

void AliFemtoXiTrackCut::SetEtaBac(double x){
  fMaxEtaBac = x;
}

void AliFemtoXiTrackCut::SetPtBac(double min, double max){
  fMinPtBac = min;
  fMaxPtBac = max;
}


void AliFemtoXiTrackCut::SetTPCnclsBac(int x){
  fTPCNclsBac = x; 
}

void AliFemtoXiTrackCut::SetNdofBac(double x){
  fNdofBac = x;
}


void AliFemtoXiTrackCut::SetStatusBac(unsigned long x) {
  fStatusBac = x;
}

void AliFemtoXiTrackCut::SetInvariantMassXi(double min, double max)
{
  fInvMassXiMin = min;
  fInvMassXiMax = max;

}

void AliFemtoXiTrackCut::SetMaxDcaXi(double x)
{
  fMaxDcaXi = x;
}

void AliFemtoXiTrackCut::SetMinDcaXiBac(double x)
{
  fMinDcaXiBac = x;
}

void AliFemtoXiTrackCut::SetMaxDcaXiDaughters(double x)
{
  fMaxDcaXiDaughters = x;
}

void AliFemtoXiTrackCut::SetMinCosPointingAngleXi(double max)
{
  fMinCosPointingAngleXi = max;
}

void AliFemtoXiTrackCut::SetParticleTypeXi(short x)
{
  fParticleTypeXi = x;
}

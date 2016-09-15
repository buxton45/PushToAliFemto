#include "AliFemtoV0TrackCut.h"
#include "AliESDtrack.h"
#include <cstdio>

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoV0TrackCut);
  /// \endcond
#endif


AliFemtoV0TrackCut::AliFemtoV0TrackCut():
  fInvMassLambdaMin(0.0),
  fInvMassLambdaMax(99.0),
  fInvMassK0sMin(0.0),
  fInvMassK0sMax(99.0),
  fMinDcaDaughterPosToVert(0.0),
  fMinDcaDaughterNegToVert(0.0),
  fMaxDcaV0Daughters(99.0),
  fMaxDcaV0(99.0),
  fMinDcaV0(0.0),
  fMaxCosPointingAngle(0.0),
  fMinCosPointingAngle(0.0),
  fParticleType(99.0),
  fEta(0.8),
  fPtMin(0.0),
  fPtMax(100.0),
  fOnFlyStatus(kFALSE),
  fMaxEtaDaughters(100.0),
  fTPCNclsDaughters(0),
  fNdofDaughters(10),
  fStatusDaughters(0),
  fPtMinPosDaughter(0.0),
  fPtMaxPosDaughter(99.0),
  fPtMinNegDaughter(0.0),
  fPtMaxNegDaughter(99.0),
  fMinAvgSepDaughters(0)
{
  // Default constructor
  // Default constructor
  fTotalV0 = 0;
  fPassedV0 = 0;
  for(int i=0; i<29; i++)
  {
    fFailedV0Tracker[i] = 0;
    fTrueFailedV0Tracker[i] = 0;
  }
}
//------------------------------
AliFemtoV0TrackCut::~AliFemtoV0TrackCut()
{
  /* noop */
}
//------------------------------
bool AliFemtoV0TrackCut::Pass(const AliFemtoV0* aV0)
{
  // test the particle and return
  // true if it meets all the criteria
  // false if it doesn't meet at least one of the criteria

cout << "&&&&&&&&&&&&&&&&&&&&&& V0 &&&&&&&&&&&&&&&&&&&" << endl;
cout << "fTotalV0 = " << fTotalV0 << endl;
cout << "fPassedV0 - " << fPassedV0 << endl;

  for(int i=0; i<29; i++)
  {
    cout << "fFailedV0Tracker[" << i << "] = " << fFailedV0Tracker[i] << endl;
  }
  for(int i=0; i<29; i++)
  {
    cout << "fTrueFailedV0Tracker[" << i << "] = " << fTrueFailedV0Tracker[i] << endl;
  }

  fTotalV0++;

  Float_t pt = aV0->PtV0();
  Float_t eta = aV0->EtaV0();


//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //kinematic cuts
  if (TMath::Abs(eta) > fEta) 
  {
    fTrueFailedV0Tracker[0]++;
  }
  if (pt < fPtMin || fPtMax < pt) 
  {
    fTrueFailedV0Tracker[1]++;
  }

  if (TMath::Abs(aV0->EtaPos()) > fMaxEtaDaughters)
  {
    fTrueFailedV0Tracker[2]++;
  }
  if (TMath::Abs(aV0->EtaNeg()) > fMaxEtaDaughters)
  {
    fTrueFailedV0Tracker[3]++;
  }

  if (aV0->PtPos() < fPtMinPosDaughter)
  {
    fTrueFailedV0Tracker[4]++;
  }
  if (aV0->PtNeg() < fPtMinNegDaughter)
  {
    fTrueFailedV0Tracker[5]++;
  }
  if (aV0->PtPos() > fPtMaxPosDaughter)
  {
    fTrueFailedV0Tracker[6]++;
  }
  if (aV0->PtNeg() > fPtMaxNegDaughter)
  {
    fTrueFailedV0Tracker[7]++;
  }

  //V0 from kinematics information
  if (fParticleType == kLambdaMC ) {
    if (!(aV0->MassLambda() > fInvMassLambdaMin && aV0->MassLambda() < fInvMassLambdaMax) || !(aV0->PosNSigmaTPCP() == 0)) {
      fTrueFailedV0Tracker[8]++;
    }
  } else if (fParticleType == kAntiLambdaMC) {
    if (!(aV0->MassLambda() > fInvMassLambdaMin && aV0->MassLambda() < fInvMassLambdaMax) || !(aV0->NegNSigmaTPCP() == 0)) {
      fTrueFailedV0Tracker[9]++;
    }
  } else if (fParticleType == kAll) {
    if (!(aV0->MassK0Short() > fInvMassK0sMin && aV0->MassK0Short() < fInvMassK0sMax) || !(aV0->NegNSigmaTPCP() == 0)) {
      fTrueFailedV0Tracker[10]++;
    }
  }

  //quality cuts
  if (aV0->OnFlyStatusV0() != fOnFlyStatus)
  {
    fTrueFailedV0Tracker[11]++;
  }
  if (aV0->StatusNeg() == 999 || aV0->StatusPos() == 999)
  {
    fTrueFailedV0Tracker[12]++;
  }
  if (aV0->TPCNclsPos() < fTPCNclsDaughters)
  {
    fTrueFailedV0Tracker[13]++;
  }
  if (aV0->TPCNclsNeg() < fTPCNclsDaughters)
  {
    fTrueFailedV0Tracker[14]++;
  }
  if (aV0->NdofPos() > fNdofDaughters)
  {
    fTrueFailedV0Tracker[15]++;
  }
  if (aV0->NdofNeg() > fNdofDaughters)
  {
    fTrueFailedV0Tracker[16]++;
  }
  if (!(aV0->StatusNeg() & fStatusDaughters))
  {
    fTrueFailedV0Tracker[17]++;
  }
  if (!(aV0->StatusPos() & fStatusDaughters))
  {
    fTrueFailedV0Tracker[18]++;
  }


  //DCA between daughter particles
  if (TMath::Abs(aV0->DcaV0Daughters()) > fMaxDcaV0Daughters)
  {
    fTrueFailedV0Tracker[19]++;
  }

  //DCA of daughters to primary vertex
  if (TMath::Abs(aV0->DcaPosToPrimVertex()) < fMinDcaDaughterPosToVert || TMath::Abs(aV0->DcaNegToPrimVertex()) < fMinDcaDaughterNegToVert)
  {
    fTrueFailedV0Tracker[20]++;
  }

  //DCA V0 to prim vertex
  if (TMath::Abs(aV0->DcaV0ToPrimVertex()) > fMaxDcaV0 || TMath::Abs(aV0->DcaV0ToPrimVertex()) < fMinDcaV0)
  {
    fTrueFailedV0Tracker[21]++;
  }

  //becomes obsolete - wrong name of the data memeber and the corresnponding methods (by default is set to fMaxCosPointingAngle = 0.0)
  //cos pointing angle
  if (aV0->CosPointingAngle() < fMaxCosPointingAngle)
  {
    fTrueFailedV0Tracker[22]++;
  }

  //this is the correct name of the data member and the corresponding methods (we are accepting cos(pointing angle bigger than certain minimum)
  //cos pointing angle
  if (aV0->CosPointingAngle() < fMinCosPointingAngle)
  {
    fTrueFailedV0Tracker[23]++;
  }

  //decay length
  if (aV0->DecayLengthV0() > fMaxDecayLength)
  {
    fTrueFailedV0Tracker[24]++;
  }


  bool pid_check1 = false;
  // Looking for lambdas = proton + pim
  if (fParticleType == kLambda) {
    if (IsProtonNSigma(aV0->PtPos(), aV0->PosNSigmaTPCP(), aV0->PosNSigmaTOFP())) //proton
      if (IsPionNSigma(aV0->PtNeg(), aV0->NegNSigmaTPCPi(), aV0->NegNSigmaTOFPi())) { //pion
        pid_check1 = true;
        //invariant mass lambda
        if (aV0->MassLambda() < fInvMassLambdaMin || aV0->MassLambda() > fInvMassLambdaMax)
        {
          fTrueFailedV0Tracker[25]++;
        }
      }

  }//Looking for antilambdas =  antiproton + pip
  else if (fParticleType == kAntiLambda) {
    if (IsProtonNSigma(aV0->PtNeg(), aV0->NegNSigmaTPCP(), aV0->NegNSigmaTOFP())) //proton
      if (IsPionNSigma(aV0->PtPos(), aV0->PosNSigmaTPCPi(), aV0->PosNSigmaTOFPi())) { //pion
        pid_check1 = true;
        //invariant mass antilambda
        if (aV0->MassAntiLambda() < fInvMassLambdaMin || aV0->MassAntiLambda() > fInvMassLambdaMax)
        {
          fTrueFailedV0Tracker[26]++;
        }
      }
  }//Looking for K0s = pip + pim
  else if (fParticleType == kK0s) {
    if (IsPionNSigma(aV0->PtNeg(), aV0->NegNSigmaTPCPi(), aV0->NegNSigmaTOFPi())) //pion-
      if (IsPionNSigma(aV0->PtPos(), aV0->PosNSigmaTPCPi(), aV0->PosNSigmaTOFPi())) { //pion+
        pid_check1 = true;
        //invariant mass K0s
        if (aV0->MassK0Short() < fInvMassK0sMin || aV0->MassK0Short() > fInvMassK0sMax)
        {
          fTrueFailedV0Tracker[27]++;
        }
      }
  }

//cout << "aV0->PtNeg(), aV0->NegNSigmaTPCP(), aV0->NegNSigmaTOFP() = " << aV0->PtNeg() << ", " << aV0->NegNSigmaTPCP() << ", " << aV0->NegNSigmaTOFP() << endl;
//cout << "aV0->PtPos(), aV0->PosNSigmaTPCPi(), aV0->PosNSigmaTOFPi() = " << aV0->PtPos() << ", " << aV0->PosNSigmaTPCPi() << ", " << aV0->PosNSigmaTOFPi() << endl;

  if (!pid_check1) 
  {
    fTrueFailedV0Tracker[28]++;
  }


//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&








  //kinematic cuts
  if (TMath::Abs(eta) > fEta) 
  {
    fFailedV0Tracker[0]++;
    return false;  //put in kinematic cuts by hand
  }
  if (pt < fPtMin || fPtMax < pt) 
  {
    fFailedV0Tracker[1]++;
    return false;
  }

  if (TMath::Abs(aV0->EtaPos()) > fMaxEtaDaughters)
  {
    fFailedV0Tracker[2]++;
    return false;
  }
  if (TMath::Abs(aV0->EtaNeg()) > fMaxEtaDaughters)
  {
    fFailedV0Tracker[3]++;
    return false;
  }

  if (aV0->PtPos() < fPtMinPosDaughter)
  {
//    cout << "aV0->PtPos() = " << aV0->PtPos() << endl;
//    cout << "fPtMinPosDaughter = " << fPtMinPosDaughter << endl;
    fFailedV0Tracker[4]++;
    return false;
  }
  if (aV0->PtNeg() < fPtMinNegDaughter)
  {
//    cout << "aV0->PtNeg() = " << aV0->PtNeg() << endl;
//    cout << "fPtMinNegDaughter = " << fPtMinNegDaughter << endl;
    fFailedV0Tracker[5]++;
    return false;
  }
  if (aV0->PtPos() > fPtMaxPosDaughter)
  {
    fFailedV0Tracker[6]++;
    return false;
  }
  if (aV0->PtNeg() > fPtMaxNegDaughter)
  {
    fFailedV0Tracker[7]++;
    return false;
  }

  //V0 from kinematics information
  if (fParticleType == kLambdaMC ) {
    if (!(aV0->MassLambda() > fInvMassLambdaMin && aV0->MassLambda() < fInvMassLambdaMax) || !(aV0->PosNSigmaTPCP() == 0)) {
      fFailedV0Tracker[8]++;
      return false;
    } else {
      fPassedV0++;
      return true;
    }
  } else if (fParticleType == kAntiLambdaMC) {
    if (!(aV0->MassLambda() > fInvMassLambdaMin && aV0->MassLambda() < fInvMassLambdaMax) || !(aV0->NegNSigmaTPCP() == 0)) {
      fFailedV0Tracker[9]++;
      return false;
    } else {
      fPassedV0++;
      return true;
    }
  } else if (fParticleType == kAll) {
    if (!(aV0->MassK0Short() > fInvMassK0sMin && aV0->MassK0Short() < fInvMassK0sMax) || !(aV0->NegNSigmaTPCP() == 0)) {
      fFailedV0Tracker[10]++;
      return false;
    } else {
      fPassedV0++;
      return true;
    }
  }

  //quality cuts
  if (aV0->OnFlyStatusV0() != fOnFlyStatus)
  {
    fFailedV0Tracker[11]++;
    return false;
  }
  if (aV0->StatusNeg() == 999 || aV0->StatusPos() == 999)
  {
    fFailedV0Tracker[12]++;
    return false;
  }
  if (aV0->TPCNclsPos() < fTPCNclsDaughters)
  {
    fFailedV0Tracker[13]++;
    return false;
  }
  if (aV0->TPCNclsNeg() < fTPCNclsDaughters)
  {
    fFailedV0Tracker[14]++;
    return false;
  }
  if (aV0->NdofPos() > fNdofDaughters)
  {
    fFailedV0Tracker[15]++;
    return false;
  }
  if (aV0->NdofNeg() > fNdofDaughters)
  {
    fFailedV0Tracker[16]++;
    return false;
  }
  if (!(aV0->StatusNeg() & fStatusDaughters))
  {
    fFailedV0Tracker[17]++;
    return false;
  }
  if (!(aV0->StatusPos() & fStatusDaughters))
  {
    fFailedV0Tracker[18]++;
    return false;
  }


  //DCA between daughter particles
  if (TMath::Abs(aV0->DcaV0Daughters()) > fMaxDcaV0Daughters)
  {
    fFailedV0Tracker[19]++;
    return false;
  }

  //DCA of daughters to primary vertex
  if (TMath::Abs(aV0->DcaPosToPrimVertex()) < fMinDcaDaughterPosToVert || TMath::Abs(aV0->DcaNegToPrimVertex()) < fMinDcaDaughterNegToVert)
  {
    fFailedV0Tracker[20]++;
    return false;
  }

  //DCA V0 to prim vertex
  if (TMath::Abs(aV0->DcaV0ToPrimVertex()) > fMaxDcaV0 || TMath::Abs(aV0->DcaV0ToPrimVertex()) < fMinDcaV0)
  {
    fFailedV0Tracker[21]++;
    return false;
  }

  //becomes obsolete - wrong name of the data memeber and the corresnponding methods (by default is set to fMaxCosPointingAngle = 0.0)
  //cos pointing angle
  if (aV0->CosPointingAngle() < fMaxCosPointingAngle)
  {
    fFailedV0Tracker[22]++;
    return false;
  }

  //this is the correct name of the data member and the corresponding methods (we are accepting cos(pointing angle bigger than certain minimum)
  //cos pointing angle
  if (aV0->CosPointingAngle() < fMinCosPointingAngle)
  {
    fFailedV0Tracker[23]++;
    return false;
  }

  //decay length
  if (aV0->DecayLengthV0() > fMaxDecayLength)
  {
    fFailedV0Tracker[24]++;
    return false;
  }


  if (fParticleType == kAll)
  {
    fPassedV0++;
    return true;
  }

  bool pid_check = false;
  // Looking for lambdas = proton + pim
  if (fParticleType == kLambda) {
    if (IsProtonNSigma(aV0->PtPos(), aV0->PosNSigmaTPCP(), aV0->PosNSigmaTOFP())) //proton
      if (IsPionNSigma(aV0->PtNeg(), aV0->NegNSigmaTPCPi(), aV0->NegNSigmaTOFPi())) { //pion
        pid_check = true;
        //invariant mass lambda
        if (aV0->MassLambda() < fInvMassLambdaMin || aV0->MassLambda() > fInvMassLambdaMax)
        {
          fFailedV0Tracker[25]++;
          return false;
        }
      }

  }//Looking for antilambdas =  antiproton + pip
  else if (fParticleType == kAntiLambda) {
    if (IsProtonNSigma(aV0->PtNeg(), aV0->NegNSigmaTPCP(), aV0->NegNSigmaTOFP())) //proton
      if (IsPionNSigma(aV0->PtPos(), aV0->PosNSigmaTPCPi(), aV0->PosNSigmaTOFPi())) { //pion
        pid_check = true;
        //invariant mass antilambda
        if (aV0->MassAntiLambda() < fInvMassLambdaMin || aV0->MassAntiLambda() > fInvMassLambdaMax)
        {
          fFailedV0Tracker[26]++;
          return false;
        }
      }
  }//Looking for K0s = pip + pim
  else if (fParticleType == kK0s) {
    if (IsPionNSigma(aV0->PtNeg(), aV0->NegNSigmaTPCPi(), aV0->NegNSigmaTOFPi())) //pion-
      if (IsPionNSigma(aV0->PtPos(), aV0->PosNSigmaTPCPi(), aV0->PosNSigmaTOFPi())) { //pion+
        pid_check = true;
        //invariant mass K0s
        if (aV0->MassK0Short() < fInvMassK0sMin || aV0->MassK0Short() > fInvMassK0sMax)
        {
          fFailedV0Tracker[27]++;
          return false;
        }
      }
  }

//cout << "aV0->PtNeg(), aV0->NegNSigmaTPCP(), aV0->NegNSigmaTOFP() = " << aV0->PtNeg() << ", " << aV0->NegNSigmaTPCP() << ", " << aV0->NegNSigmaTOFP() << endl;
//cout << "aV0->PtPos(), aV0->PosNSigmaTPCPi(), aV0->PosNSigmaTOFPi() = " << aV0->PtPos() << ", " << aV0->PosNSigmaTPCPi() << ", " << aV0->PosNSigmaTOFPi() << endl;

  if (!pid_check) 
  {
    fFailedV0Tracker[28]++;
    return false;
  }


  fPassedV0++;
  return true;
}
//------------------------------
AliFemtoString AliFemtoV0TrackCut::Report()
{
  // Prepare report from the execution
  string tStemp;
  char tCtemp[100];
  snprintf(tCtemp, 100, "Minimum of Invariant Mass assuming Lambda:\t%lf\n", fInvMassLambdaMin);
  tStemp += tCtemp;
  snprintf(tCtemp, 100, "Maximum of Invariant Mass assuming Lambda:\t%lf\n", fInvMassLambdaMax);
  tStemp += tCtemp;
  snprintf(tCtemp, 100, "Minimum DCA of positive daughter to primary vertex:\t%lf\n", fMinDcaDaughterPosToVert);
  tStemp += tCtemp;
  snprintf(tCtemp, 100, "Minimum DCA of negative daughter to primary vertex:\t%lf\n", fMinDcaDaughterNegToVert);
  tStemp += tCtemp;
  snprintf(tCtemp, 100, "Max DCA of daughters:\t%lf\n", fMaxDcaV0Daughters);
  tStemp += tCtemp;

  AliFemtoString returnThis = tStemp;
  return returnThis;
}
TList *AliFemtoV0TrackCut::ListSettings()
{
  // return a list of settings in a writable form
  TList *tListSetttings = new TList();
  char buf[200];
  snprintf(buf, 200, "AliFemtoV0TrackCut.InvMassLambdaMin=%lf", fInvMassLambdaMin);
  tListSetttings->AddLast(new TObjString(buf));
  return tListSetttings;
}

void AliFemtoV0TrackCut::SetMinDaughtersToPrimVertex(double minPos, double minNeg)
{
  fMinDcaDaughterPosToVert = minPos;
  fMinDcaDaughterNegToVert = minNeg;
}

void AliFemtoV0TrackCut::SetMaxDcaV0Daughters(double max)
{
  fMaxDcaV0Daughters = max;
}

void AliFemtoV0TrackCut::SetMaxV0DecayLength(double max)
{
  fMaxDecayLength = max;
}

void AliFemtoV0TrackCut::SetMaxDcaV0(double max)
{
  fMaxDcaV0 = max;
}

void AliFemtoV0TrackCut::SetMinDcaV0(double min)
{
  fMinDcaV0 = min;
}

void AliFemtoV0TrackCut::SetMaxCosPointingAngle(double max) //obsolete
{
  fMaxCosPointingAngle = max;
}

void AliFemtoV0TrackCut::SetMinCosPointingAngle(double min) //correct
{
  fMinCosPointingAngle = min;
}

void AliFemtoV0TrackCut::SetParticleType(short x)
{
  fParticleType = x;
}

void AliFemtoV0TrackCut::SetEta(double x)
{
  fEta = x;
}

void AliFemtoV0TrackCut::SetPt(double min, double max)
{
  fPtMin = min;
  fPtMax = max;
}

void AliFemtoV0TrackCut::SetEtaDaughters(float x)
{
  fMaxEtaDaughters = x;
}
void AliFemtoV0TrackCut::SetTPCnclsDaughters(int x)
{
  fTPCNclsDaughters = x;
}
void AliFemtoV0TrackCut::SetNdofDaughters(int x)
{
  fNdofDaughters = x;
}
void AliFemtoV0TrackCut::SetStatusDaughters(unsigned long x)
{
  fStatusDaughters = x;
}
void AliFemtoV0TrackCut::SetPtPosDaughter(float min, float max)
{
  fPtMinPosDaughter = min;
  fPtMaxPosDaughter = max;
}

void AliFemtoV0TrackCut::SetPtNegDaughter(float min, float max)
{
  fPtMinNegDaughter = min;
  fPtMaxNegDaughter = max;
}

void AliFemtoV0TrackCut::SetOnFlyStatus(bool x)
{
  fOnFlyStatus = x;
}

void AliFemtoV0TrackCut::SetInvariantMassLambda(double min, double max)
{
  fInvMassLambdaMin = min;
  fInvMassLambdaMax = max;
}

void AliFemtoV0TrackCut::SetInvariantMassK0s(double min, double max)
{
  fInvMassK0sMin = min;
  fInvMassK0sMax = max;
}

void AliFemtoV0TrackCut::SetMinAvgSeparation(double minSep)
{
  fMinAvgSepDaughters = minSep;
}

//---------------------PID n Sigma ---------------------------------//
bool AliFemtoV0TrackCut::IsKaonTPCdEdxNSigma(float mom, float nsigmaK)
{
  if (mom < 0.35 && TMath::Abs(nsigmaK) < 5.0) return true;
  if (mom >= 0.35 && mom < 0.5 && TMath::Abs(nsigmaK) < 3.0) return true;
  if (mom >= 0.5 && mom < 0.7 && TMath::Abs(nsigmaK) < 2.0) return true;

  return false;
}


bool AliFemtoV0TrackCut::IsKaonTOFNSigma(float mom, float nsigmaK)
{
  if (mom >= 1.5 && TMath::Abs(nsigmaK) < 2.0) return true;
  return false;
}

bool AliFemtoV0TrackCut::IsKaonNSigma(float mom, float nsigmaTPCK, float nsigmaTOFK)
{
  if (mom < 0.4) {
    if (nsigmaTOFK < -999.) {
      if (TMath::Abs(nsigmaTPCK) < 1.0) return true;
    } else if (TMath::Abs(nsigmaTOFK) < 3.0 && TMath::Abs(nsigmaTPCK) < 3.0) return true;
  } else if (mom >= 0.4 && mom <= 0.6) {
    if (nsigmaTOFK < -999.) {
      if (TMath::Abs(nsigmaTPCK) < 2.0) return true;
    } else if (TMath::Abs(nsigmaTOFK) < 3.0 && TMath::Abs(nsigmaTPCK) < 3.0) return true;
  } else if (nsigmaTOFK < -999.) {
    return false;
  } else if (TMath::Abs(nsigmaTOFK) < 3.0 && TMath::Abs(nsigmaTPCK) < 3.0) return true;

  return false;
}



bool AliFemtoV0TrackCut::IsPionNSigma(float mom, float nsigmaTPCPi, float nsigmaTOFPi)
{
  if (TMath::Abs(nsigmaTPCPi) < 3.0) return true;

  /*if (nsigmaTOFPi<-999.)
    {
      if (TMath::Abs(nsigmaTPCPi)<3.0) return true;
    }
  else
    {
      if (TMath::Abs(nsigmaTPCPi)<3.0 && TMath::Abs(nsigmaTOFPi)<4.0) return true;
      }*/

  /*if (mom<0.65)
    {
      if (nsigmaTOFPi<-999.)
  {
    if (mom<0.35 && TMath::Abs(nsigmaTPCPi)<3.0) return true;
    else if (mom<0.5 && mom>=0.35 && TMath::Abs(nsigmaTPCPi)<3.0) return true;
    else if (mom>=0.5 && TMath::Abs(nsigmaTPCPi)<2.0) return true;
    else return false;
  }
  else if (TMath::Abs(nsigmaTOFPi)<3.0 && TMath::Abs(nsigmaTPCPi)<3.0) return true;
  if (TMath::Abs(nsigmaTPCPi)<3.0) return true;
      else return false;
    }
  else if (nsigmaTOFPi<-999.)
    {
      return false;
    }
  else if (mom<1.5 && TMath::Abs(nsigmaTOFPi)<3.0 && TMath::Abs(nsigmaTPCPi)<5.0) return true;
  else if (mom>=1.5 && TMath::Abs(nsigmaTOFPi)<2.0 && TMath::Abs(nsigmaTPCPi)<5.0) return true;*/

  return false;
}

bool AliFemtoV0TrackCut::IsProtonNSigma(float mom, float nsigmaTPCP, float nsigmaTOFP)
{
  if (mom < 0.8) {
    if (TMath::Abs(nsigmaTPCP) < 3.0) return true;
  } else {
    if (nsigmaTOFP < -999.) {
      if (TMath::Abs(nsigmaTPCP) < 3.0) return true;
    } else {
      if (TMath::Abs(nsigmaTPCP) < 3.0 && TMath::Abs(nsigmaTOFP) < 3.0) return true;
    }
  }


  /*if (nsigmaTOFP<-999.)
    {
      if (TMath::Abs(nsigmaTPCP)<3.0) return true;
    }
  else
    {
      if (TMath::Abs(nsigmaTPCP)<3.0 && TMath::Abs(nsigmaTOFP)<3.0) return true;
      }*/


  /*if (mom<0.8)
    {
      if (nsigmaTOFP<-999.)
  {
    if (TMath::Abs(nsigmaTPCP)<3.0) return true;
    else return false;
  }
  else if (TMath::Abs(nsigmaTOFP)<3.0 && TMath::Abs(nsigmaTPCP)<3.0) return true;
  if (TMath::Abs(nsigmaTPCP)<3.0) return true;
      else return false;
    }
  else if (nsigmaTOFP<-999.)
    {
      return false;
    }
    else if (TMath::Abs(nsigmaTOFP)<3.0 && TMath::Abs(nsigmaTPCP)<3.0) return true;*/

  return false;
}

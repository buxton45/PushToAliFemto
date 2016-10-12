////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  myAliFemtoModelCorrFctnKStar                                              //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(myAliFemtoModelCorrFctnKStar);
  /// \endcond
#endif

#include "myAliFemtoModelCorrFctnKStar.h"
#include "AliFemtoModelHiddenInfo.h"
#include <TH1D.h>
#include <TH2D.h>

//_________________________________________
myAliFemtoModelCorrFctnKStar::myAliFemtoModelCorrFctnKStar() :
  AliFemtoModelCorrFctn()
{
  //default constructor
}

myAliFemtoModelCorrFctnKStar::myAliFemtoModelCorrFctnKStar(const char *title, int aNbins, double aKStarLo, double aKStarHi) :
  AliFemtoModelCorrFctn(title,aNbins,aKStarLo,aKStarHi)
{
  //normal constructor
}

//_________________________________________
myAliFemtoModelCorrFctnKStar::~myAliFemtoModelCorrFctnKStar()
{
  //destructor
}

//_________________________________________
myAliFemtoModelCorrFctnKStar::myAliFemtoModelCorrFctnKStar(const myAliFemtoModelCorrFctnKStar& aCorrFctn) :
  AliFemtoModelCorrFctn(aCorrFctn)
{
  //copy constructor
}

//_________________________________________
myAliFemtoModelCorrFctnKStar& myAliFemtoModelCorrFctnKStar::operator=(const myAliFemtoModelCorrFctnKStar& aCorrFctn)
{
  //assignment operator
  if (this == &aCorrFctn) return *this;

  AliFemtoModelCorrFctn::operator=(aCorrFctn);

  return *this;
}

//_________________________________________
double myAliFemtoModelCorrFctnKStar::GetKStarTrue(AliFemtoPair* aPair)
{
  // Generate a simple femtoscopic weight coming simple Lednicky equation
  // The weight is generated using the TrueMomentum in the HiddenInfo

  AliFemtoParticle *tPart1 = (AliFemtoParticle*)aPair->Track1();
  AliFemtoParticle *tPart2 = (AliFemtoParticle*)aPair->Track2();

  if(tPart1 != NULL && tPart2 != NULL)
  {
    AliFemtoModelHiddenInfo *tPart1HiddenInfo = (AliFemtoModelHiddenInfo*)tPart1->GetHiddenInfo();
    AliFemtoModelHiddenInfo *tPart2HiddenInfo = (AliFemtoModelHiddenInfo*)tPart2->GetHiddenInfo();

    if(tPart1HiddenInfo != NULL && tPart2HiddenInfo != NULL)
    {
      double px1 = tPart1HiddenInfo->GetTrueMomentum()->x();
      double py1 = tPart1HiddenInfo->GetTrueMomentum()->y();
      double pz1 = tPart1HiddenInfo->GetTrueMomentum()->z();
      double mass1 = tPart1HiddenInfo->GetMass();
      double E1 = std::sqrt(mass1*mass1 + px1*px1 + py1*py1 + pz1*pz1);

      double px2 = tPart2HiddenInfo->GetTrueMomentum()->x();
      double py2 = tPart2HiddenInfo->GetTrueMomentum()->y();
      double pz2 = tPart2HiddenInfo->GetTrueMomentum()->z();
      double mass2 = tPart2HiddenInfo->GetMass();
      double E2 = std::sqrt(mass2*mass2 + px2*px2 + py2*py2 + pz2*pz2);
      //------------------------------------------------------------

      double tMinvSq = (E1+E2)*(E1+E2) - (px1+px2)*(px1+px2) - (py1+py2)*(py1+py2) - (pz1+pz2)*(pz1+pz2);

      double tQinvSq = ((mass1*mass1 - mass2*mass2)*(mass1*mass1 - mass2*mass2))/tMinvSq + tMinvSq - 2.0*(mass1*mass1 + mass2*mass2);

      double tKStar = 0.5*sqrt(tQinvSq);
      return tKStar;
    }

  }

  return -999;
}

//_________________________________________
void myAliFemtoModelCorrFctnKStar::AddRealPair(AliFemtoPair* aPair)
{
  double tWeight = fManager->GetWeight(aPair);
  fNumeratorTrue->Fill(aPair->KStar(),tWeight);

  double tKStarTrue = GetKStarTrue(aPair);

  fNumeratorTrueIdeal->Fill(tKStarTrue,tWeight);
}

//_________________________________________
void myAliFemtoModelCorrFctnKStar::AddMixedPair(AliFemtoPair* aPair)
{
  double tWeight = fManager->GetWeight(aPair);
  fNumeratorFake->Fill(aPair->KStar(),tWeight);
  fDenominator->Fill(aPair->KStar(), 1.0);

  double tKStarTrue = GetKStarTrue(aPair);

  fNumeratorFakeIdeal->Fill(tKStarTrue,tWeight);
  fDenominatorIdeal->Fill(tKStarTrue, 1.0);

  fQgenQrec->Fill(tKStarTrue,aPair->KStar());
}




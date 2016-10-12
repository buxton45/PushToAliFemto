/// \class AliFemtoV0TrackCutNSigmaFilter

#ifndef ALIFEMTOV0TRACKCUTNSIGMAFILTER_H
#define ALIFEMTOV0TRACKCUTNSIGMAFILTER_H

#include "AliFemtoTrackCut.h"

#include "AliFemtoNSigmaFilter.h"
class AliFemtoNSigmaFilter;

#include "TH1D.h"

/**
 * \class AliFemtoV0TrackCutNSigmaFilter

 */
class AliFemtoV0TrackCutNSigmaFilter : public AliFemtoParticleCut {
public:
  /// Enumerated type to easily select correct algorithm for each particle type
  enum V0Type {
    kLambda = 0,
    kAntiLambda = 1,
    kK0s = 2,
    kAll = 99,
    kLambdaMC = 101,
    kAntiLambdaMC = 102,
    kK0sMC = 3
  };
  typedef enum V0Type AliFemtoV0Type;

  AliFemtoV0TrackCutNSigmaFilter();
  virtual ~AliFemtoV0TrackCutNSigmaFilter();

  virtual bool Pass(const AliFemtoV0* aV0);

  virtual AliFemtoString Report();
  virtual TList *ListSettings();
  virtual AliFemtoParticleType Type(){return hbtV0;}

  void SetInvariantMassLambda(double,double);
  void SetInvariantMassK0s(double,double);
  void SetMinDaughtersToPrimVertex(double,double);
  void SetMaxDcaV0Daughters(double);
  void SetMaxDcaV0(double);
  void SetMinDcaV0(double);
  void SetMaxCosPointingAngle(double);
  void SetMinCosPointingAngle(double);
  void SetMaxV0DecayLength(double);
  void SetParticleType(short);
  void SetEta(double);
  void SetPt(double,double);
  void SetEtaDaughters(float);
  void SetTPCnclsDaughters(int);
  void SetNdofDaughters(int);
  void SetStatusDaughters(unsigned long);
  void SetPtPosDaughter(float,float);
  void SetPtNegDaughter(float,float);
  void SetOnFlyStatus(bool);
  void SetMinAvgSeparation(double);

  //----n sigma----
  bool IsKaonTPCdEdxNSigma(float mom, float nsigmaK);
  bool IsKaonTOFNSigma(float mom, float nsigmaK);
  bool IsKaonNSigma(float mom, float nsigmaTPCK, float nsigmaTOFK);  //tweaked 14/12/2015
  bool IsPionNSigma(float mom, float nsigmaTPCPi, float nsigmaTOFPi);  //tweaked 14/12/2015
  bool IsProtonNSigma(float mom, float nsigmaTPCP, float nsigmaTOFP);  //tweaked 14/12/2015


// !!!!!----- 14/12/2015 ----------------------------------------------------------
  enum DaughterParticleType {kPion=0, kKaon=1, kProton=2};

  void CreateCustomNSigmaFilter(DaughterParticleType aDaughterType);
  void CreateCustomPionNSigmaFilter();
  void CreateCustomKaonNSigmaFilter();
  void CreateCustomProtonNSigmaFilter();

  void AddTPCAndTOFNSigmaCut(DaughterParticleType aDaughterType, double aMomMin, double aMomMax, double aNSigmaValueTPC, double aNSigmaValueTOF);
  void AddPionTPCAndTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC, double aNSigmaValueTOF);
  void AddKaonTPCAndTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC, double aNSigmaValueTOF);
  void AddProtonTPCAndTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC, double aNSigmaValueTOF);

  void AddTPCNSigmaCut(DaughterParticleType aDaughterType, double aMomMin, double aMomMax, double aNSigmaValueTPC);
  void AddPionTPCNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC);
  void AddKaonTPCNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC);
  void AddProtonTPCNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC);

  void AddTOFNSigmaCut(DaughterParticleType aDaughterType, double aMomMin, double aMomMax, double aNSigmaValueTOF);
  void AddPionTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTOF);
  void AddKaonTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTOF);
  void AddProtonTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTOF);

  //-----The MinvHisto is built immediately before the (final) Minv cut, and thus may be used to calculate the purity of the V0 collection
  // However, the use MUST manually add this histogram to the output list of the analysis.
  //   i.e. add something similar to:  tOutputList->Add(p1cut->GetMinvHisto());
  //   where p1cut is a AliFemtoV0TrackCut object
  void SetMinvHisto(const char* title, const int& nbins, const float& aInvMassMin, const float& aInvMassMax);   //set the Minv histogram attributes and
														//set flag fBuildMinvHisto=true
  TH1D* GetMinvHisto();

 private:   // here are the quantities I want to cut on...

  double fInvMassLambdaMin;        ///< invariant mass lambda min
  double fInvMassLambdaMax;        ///< invariant mass lambda max
  double fInvMassK0sMin;           ///< invariant mass lambda min
  double fInvMassK0sMax;           ///< invariant mass lambda max
  double fMinDcaDaughterPosToVert; ///< DCA of positive daughter to primary vertex
  double fMinDcaDaughterNegToVert; ///< DCA of negative daughter to primary vertex
  double fMaxDcaV0Daughters;       ///< Max DCA of v0 daughters at Decay vertex
  double fMaxDcaV0;
  double fMinDcaV0;
  double fMaxDecayLength;

  double fMaxCosPointingAngle;    //obsolete
  double fMinCosPointingAngle;    //correct
  short fParticleType;             ///< 0-lambda
  double fEta;
  double fPtMin;
  double fPtMax;
  bool fOnFlyStatus;

  float fMaxEtaDaughters;         ///< Eta of positive daughter
  int fTPCNclsDaughters;          ///< No. of cls of pos daughter
  int fNdofDaughters;             ///< No. of degrees of freedom of the pos. daughter track
  unsigned long fStatusDaughters; ///< Status (tpc refit, its refit...)
  float fPtMinPosDaughter;
  float fPtMaxPosDaughter;
  float fPtMinNegDaughter;
  float fPtMaxNegDaughter;
  double fMinAvgSepDaughters;

// !!!!!----- 14/12/2015 ----------------------------------------------------------
  bool fUseCustomPionNSigmaFilter;
  bool fUseCustomKaonNSigmaFilter;
  bool fUseCustomProtonNSigmaFilter;

  AliFemtoNSigmaFilter *fPionNSigmaFilter;
  AliFemtoNSigmaFilter *fKaonNSigmaFilter;
  AliFemtoNSigmaFilter *fProtonNSigmaFilter;

  bool fBuildMinvHisto;
  TH1D* fMinvHisto;


#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoV0TrackCutNSigmaFilter, 1);
  /// \endcond
#endif

};


#endif

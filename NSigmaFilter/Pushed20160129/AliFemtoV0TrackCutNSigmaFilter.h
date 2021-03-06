///
/// \file AliFemtoV0TrackCutNSigmaFilter.h
///

#ifndef ALIFEMTOV0TRACKCUTNSIGMAFILTER_H
#define ALIFEMTOV0TRACKCUTNSIGMAFILTER_H

#pragma once

#include "AliFemtoTrackCut.h"
#include "AliFemtoV0TrackCut.h"
#include "AliFemtoNSigmaFilter.h"

class TH1D;


/// \class AliFemtoV0TrackCutNSigmaFilter
/// \brief A class designed to help track cuts test particle calculate NSigma values
///
/// \author Jesse Buxton <jesse.thomas.buxton@cern.ch>
///
class AliFemtoV0TrackCutNSigmaFilter : public AliFemtoV0TrackCut {
public:

  /// Default Constructor - wide ranges and
  AliFemtoV0TrackCutNSigmaFilter();

  /// Destructor - deletes filters and histograms
  virtual ~AliFemtoV0TrackCutNSigmaFilter();

  /// Main method of the filter. Returns true if V0 passes NSigma
  virtual bool Pass(const AliFemtoV0* aV0);

  virtual AliFemtoString Report();
  virtual TList* ListSettings();
  virtual TList* AppendSettings(TList *settings, const TString &prefix = "") const;

  //----n sigma----
  bool IsKaonNSigma(float mom, float nsigmaTPCK, float nsigmaTOFK);
  bool IsPionNSigma(float mom, float nsigmaTPCPi, float nsigmaTOFPi);
  bool IsProtonNSigma(float mom, float nsigmaTPCP, float nsigmaTOFP);

  enum DaughterParticleType {
    kPion = 0,
    kKaon = 1,
    kProton = 2
  };

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

  //DO NOT set these to true unless you completely understand the consequences.  See AliFemtoNSigmaFilter.h for more information
  void SetOverrideImproperPionNSigmaFilter(bool aOverride);
  void SetOverrideImproperKaonNSigmaFilter(bool aOverride);
  void SetOverrideImproperProtonNSigmaFilter(bool aOverride);

protected:

  bool fUseCustomPionNSigmaFilter;
  bool fUseCustomKaonNSigmaFilter;
  bool fUseCustomProtonNSigmaFilter;

  AliFemtoNSigmaFilter *fPionNSigmaFilter;
  AliFemtoNSigmaFilter *fKaonNSigmaFilter;
  AliFemtoNSigmaFilter *fProtonNSigmaFilter;


#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoV0TrackCutNSigmaFilter, 1);
  /// \endcond
#endif

};



inline void AliFemtoV0TrackCutNSigmaFilter::SetOverrideImproperPionNSigmaFilter(bool aOverride)
{
  fPionNSigmaFilter->SetOverrideImproperConfig(aOverride);
}

inline void AliFemtoV0TrackCutNSigmaFilter::SetOverrideImproperKaonNSigmaFilter(bool aOverride)
{
  fKaonNSigmaFilter->SetOverrideImproperConfig(aOverride);
}

inline void AliFemtoV0TrackCutNSigmaFilter::SetOverrideImproperProtonNSigmaFilter(bool aOverride)
{
  fProtonNSigmaFilter->SetOverrideImproperConfig(aOverride);
}

inline TList* AliFemtoV0TrackCutNSigmaFilter::ListSettings()
{
  return AppendSettings(new TList());
}

inline void AliFemtoV0TrackCutNSigmaFilter::CreateCustomPionNSigmaFilter() {
  CreateCustomNSigmaFilter(kPion);
}

inline void AliFemtoV0TrackCutNSigmaFilter::CreateCustomKaonNSigmaFilter()
{
  CreateCustomNSigmaFilter(kKaon);
}

inline void AliFemtoV0TrackCutNSigmaFilter::CreateCustomProtonNSigmaFilter()
{
  CreateCustomNSigmaFilter(kProton);
}

inline void AliFemtoV0TrackCutNSigmaFilter::AddPionTPCAndTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC, double aNSigmaValueTOF)
{
  AddTPCAndTOFNSigmaCut(kPion, aMomMin, aMomMax, aNSigmaValueTPC, aNSigmaValueTOF);
}

inline void AliFemtoV0TrackCutNSigmaFilter::AddKaonTPCAndTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC, double aNSigmaValueTOF)
{
  AddTPCAndTOFNSigmaCut(kKaon, aMomMin, aMomMax, aNSigmaValueTPC, aNSigmaValueTOF);
}

inline void AliFemtoV0TrackCutNSigmaFilter::AddProtonTPCAndTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC, double aNSigmaValueTOF)
{
  AddTPCAndTOFNSigmaCut(kProton, aMomMin, aMomMax, aNSigmaValueTPC, aNSigmaValueTOF);
}


inline void AliFemtoV0TrackCutNSigmaFilter::AddPionTPCNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC)
{
  AddTPCNSigmaCut(kPion, aMomMin, aMomMax, aNSigmaValueTPC);
}

inline void AliFemtoV0TrackCutNSigmaFilter::AddKaonTPCNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC)
{
  AddTPCNSigmaCut(kKaon, aMomMin, aMomMax, aNSigmaValueTPC);
}

inline void AliFemtoV0TrackCutNSigmaFilter::AddProtonTPCNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC)
{
  AddTPCNSigmaCut(kProton, aMomMin, aMomMax, aNSigmaValueTPC);
}



inline void AliFemtoV0TrackCutNSigmaFilter::AddPionTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTOF)
{
  AddTOFNSigmaCut(kPion, aMomMin, aMomMax, aNSigmaValueTOF);
}

inline void AliFemtoV0TrackCutNSigmaFilter::AddKaonTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTOF)
{
  AddTOFNSigmaCut(kKaon, aMomMin, aMomMax, aNSigmaValueTOF);
}

inline void AliFemtoV0TrackCutNSigmaFilter::AddProtonTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTOF)
{
  AddTOFNSigmaCut(kProton, aMomMin, aMomMax, aNSigmaValueTOF);
}




#endif

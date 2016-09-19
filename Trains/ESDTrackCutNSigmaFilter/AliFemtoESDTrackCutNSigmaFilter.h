////////////////////////////////////////////////////////////////////////////
//                                                                        //
// AliFemtoESDTrackCutNSigmaFilter:                                       //
//  This will give me more control over the Kch track cuts                //
//  For now, the only this I want to change is the method IsKaonNSigma    //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOESDTRACKCUTNSIGMAFILTER_H
#define ALIFEMTOESDTRACKCUTNSIGMAFILTER_H

#include "AliFemtoESDTrackCut.h"

class AliFemtoESDTrackCutNSigmaFilter : public AliFemtoESDTrackCut {

public:

  AliFemtoESDTrackCutNSigmaFilter();
  AliFemtoESDTrackCutNSigmaFilter(const AliFemtoESDTrackCutNSigmaFilter& aCut);
  AliFemtoESDTrackCutNSigmaFilter& operator=(const AliFemtoESDTrackCutNSigmaFilter& aCut);
  virtual ~AliFemtoESDTrackCutNSigmaFilter();

  virtual bool Pass(const AliFemtoTrack* aTrack);

  void SetPionRejection(bool aReject);

private:
  bool fPionRejection;  // 26/02/2016

  bool IsKaonNSigma(float mom, float nsigmaTPC, float nsigmaTOF);

  bool IsElectron(float nsigmaTPCE);  // 26/02/2016


#ifdef __ROOT__
  ClassDef(AliFemtoESDTrackCutNSigmaFilter, 1)
#endif

};

inline void AliFemtoESDTrackCutNSigmaFilter::SetPionRejection(bool aReject) {fPionRejection = aReject;}

#endif

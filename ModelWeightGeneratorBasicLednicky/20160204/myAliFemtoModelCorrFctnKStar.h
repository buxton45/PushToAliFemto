////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  myAliFemtoModelCorrFctnKStar                                              //
//                                                                            //
//  Author: Jesse Buxton jesse.thomas.buxton@cern.ch                          //
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOMODELCORRFCTNKSTAR_H
#define ALIFEMTOMODELCORRFCTNKSTAR_H

#include "AliFemtoCorrFctn.h"
#include "AliFemtoModelCorrFctn.h"
#include "AliFemtoModelHiddenInfo.h"
#include "AliFemtoPair.h"
#include "AliFemtoModelManager.h"


class myAliFemtoModelCorrFctnKStar : public AliFemtoModelCorrFctn 
{
public:
  myAliFemtoModelCorrFctnKStar();
  myAliFemtoModelCorrFctnKStar(const char *title, int aNbins, double aKStarLo, double aKStarHi);
  virtual ~myAliFemtoModelCorrFctnKStar();

  myAliFemtoModelCorrFctnKStar(const myAliFemtoModelCorrFctnKStar& aCorrFctn);
  myAliFemtoModelCorrFctnKStar& operator=(const myAliFemtoModelCorrFctnKStar& aCorrFctn);

  double GetKStarTrue(AliFemtoPair* aPair);

  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPair);




protected:








#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(myAliFemtoModelCorrFctnKStar, 1);
  /// \endcond
#endif

};

#endif

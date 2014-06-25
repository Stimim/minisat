#ifndef STIMIM_GB_ANALYZER_H
#define STIMIM_GB_ANALYZER_H

#include "polybori.h"
#include "polybori/grobner/groebner_alg.h"

namespace GroebnerBasis {
  struct Analyzer {
    enum ReturnCode {
      SKIPPED,
      CONST_1,
      CONST_0,
      OTHER,
    };

    polybori::BoolePolyRing ring;
    polybori::groebner::CacheManager cacheMgr;

    Analyzer(size_t nVar) : ring(nVar, polybori::COrderEnums::dlex), cacheMgr() {}

    bool checkClock(unsigned conflicts, unsigned starts) const {
      return false;
    }
  };
};

#endif

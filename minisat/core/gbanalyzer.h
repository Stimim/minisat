#ifndef STIMIM_GB_ANALYZER_H
#define STIMIM_GB_ANALYZER_H

#include "polybori.h"
#include "polybori/groebner/groebner_alg.h"

namespace GroebnerBasis {

typedef polybori::BoolePolynomial Poly;
typedef polybori::BooleVariable   Var;

struct PolyComp {
  bool operator() (const Poly& x, const Poly& y) const {
    return x.compare(y) < 0;
  }
};

struct VarComp {
  const Minisat::VMap<double>& solver_act;
  const std::vector<double>& activity;

  VarComp(const Minisat::VMap<double>& solver_act,
          const std::vector<double>& activity) :
      solver_act(solver_act), activity(activity) { }

  bool operator() (const Minisat::Var& u, const Minisat::Var& v) {
    return solver_act[v] * activity[v] < solver_act[u] * activity[u];
  }
};

struct Analyzer {
  enum ReturnCode {
    SKIPPED,
    CONST_1,
    CONST_0,
    OTHER,
  };

  polybori::BoolePolyRing ring;
  polybori::groebner::CacheManager cacheMgr;

  std::vector<double> activity;

  Analyzer(size_t nVar) : ring(nVar, polybori::COrderEnums::dlex), cacheMgr(), activity(nVar, 1) {}

  void reorder_vars(std::vector<Minisat::Var>& vars,
                    const Minisat::VMap<double>& solver_act) {
    std::sort(vars.begin(), vars.end(), VarComp(solver_act, activity));
  }

  bool check_clock(unsigned conflicts, unsigned starts) const {
    unsigned long long n = 50000;
    return (1+conflicts) % n == 0;
  }

  void collect_vars(const std::vector<std::vector<Minisat::Lit> >& reasons,
                    std::vector<Minisat::Var>& vars) {
    std::set<Minisat::Var> _vars;
    for (size_t i = 0; i < reasons.size(); ++ i) {
      for (size_t j = 0; j < reasons[i].size(); ++ j) {
        _vars.insert(var(reasons[i][j]));
      }
    }
    vars.clear();
    for (std::set<Minisat::Var>::const_iterator it = _vars.begin();
         it != _vars.end(); ++ it) {
      vars.push_back(*it);
    }
  }

  ReturnCode analyze(const std::vector<std::vector<Minisat::Lit> >& reasons,
                     std::vector<std::vector<Minisat::Lit> >& out,
                     const Minisat::VMap<double>& solver_act) {
    std::vector<Minisat::Var> vars;
    collect_vars(reasons, vars);

    if (skip_this_time(reasons, vars)) return SKIPPED;

    reorder_vars(vars, solver_act);

    std::map<Minisat::Var, unsigned> rev_map;
    for (size_t i = 0; i < vars.size(); ++ i) {
      rev_map[vars[i]] = i;
    }

    polybori::groebner::GroebnerStrategy strat(ring);
    initStrategy(strat);

    std::vector<Poly> polys;
    std::set<Poly, PolyComp> origs;

    convert_to_polys(reasons, rev_map, polys, origs);

    std::sort(polys.begin(), polys.end(), PolyComp());

    for (std::vector<Poly>::const_iterator i = polys.begin();
         i != polys.end(); ++ i) {
      std::vector<Poly> v = full_implication_gb(*i, cacheMgr, strat);
      for (std::vector<polybori::BoolePolynomial>::const_iterator j = v.begin();
           j != v.end(); ++ j) {
        const polybori::BoolePolynomial& p = strat.nf(*j);
        if (p.isZero()) {
          continue; // skip
        } else if (p.isOne()) {
          return CONST_1;
        }
        strat.addGenerator(p);
      }
    }

    strat.symmGB_F2();
    const std::vector<polybori::BoolePolynomial>& basis =
        strat.minimalizeAndTailReduce();
    // std::cout << "----- basis -----" << std::endl;
    for (std::vector<polybori::BoolePolynomial>::const_iterator it = basis.begin();
         it != basis.end(); ++ it) {
      // it->print(std::cout) << std::endl;
      if (it->isOne()) {
        return CONST_1;
      } else if (it->isZero()) {
        continue;
      } else if (origs.count(*it) > 0) {
        continue;
      }

      const polybori::BooleMonomial& m = it->usedVariables();
      std::vector<Minisat::Lit> cla;
      bool good = true;
      for (polybori::BooleMonomial::variable_iterator v = m.variableBegin();
           v != m.variableEnd(); ++ v) {
        if ((*it % *v).isZero()) {
          cla.push_back(Minisat::mkLit(vars[v->index()], false));
        } else if ((*it % (*v+1)).isZero()) {
          cla.push_back(Minisat::mkLit(vars[v->index()], true));
        } else {
          std::cout << "Can't convert ";
          it->print(std::cout) << std::endl;
          good = false;
          break;
        }
      }

      if (good) {
        out.push_back(cla);
      }

      if (cla.size() <= 2) {
        for (unsigned i = 0; i < cla.size(); ++ i) {
          activity[var(cla[i])] += 1;
        }
      }
    }

    return OTHER;
  }

  void initStrategy(polybori::groebner::GroebnerStrategy& strat) const {
    strat.optLazy = true;
    strat.optDrawMatrices = false;
    strat.optModifiedLinearAlgebra = true;
    strat.enabledLog = false;
    strat.optExchange = true;
    strat.optAllowRecursion = false;
    strat.optLinearAlgebraInLastBlock = true;
  }




  static const double THRES = 2.0;

  bool skip_this_time(const std::vector<std::vector<Minisat::Lit> >& reasons,
                      const std::vector<Minisat::Var>& vars) const {
    if ((double) vars.size() / reasons.size() > THRES) {
      return true;
    }
    return false;
  }

  void convert_to_polys(const std::vector<std::vector<Minisat::Lit> >& reasons,
                        std::map<Minisat::Var, unsigned>& rev_map,
                        std::vector<Poly>& polys,
                        std::set<Poly, PolyComp>& origs) {
    for (size_t i = 0; i < reasons.size(); ++ i) {
      Poly poly(true, ring);
      const std::vector<Minisat::Lit>& c = reasons[i];
      for (int j = 0, sz = c.size(); j < sz; ++ j) {
        Minisat::Lit p = c[j];
        if (sign(p)) {
          poly *= Var(rev_map[var(p)], ring) + 1;
        } else {
          poly *= Var(rev_map[var(p)], ring);
        }
      }
      if (poly.isZero()) continue;

      polys.push_back(poly);
      origs.insert(poly);
    }
  }
};
};

#endif

#include <cassert>
/* assert macro */
#include <fstream>
using std::ifstream;
#include <sstream>
using std::ostringstream;
#include <string>
using std::string;
#include <vector>
using std::vector;

#include <iostream>
using std::cout;

#include "selection.h"
#include "small_powers.h"
#include "settings.h"
#include "utilities.h"
#include "stablepoint.h"

#define require(condition) if (!(condition)) return false;
#define veto(condition) if (condition) return false;

Selection::Selection()
{}

// helper: ------------------------------------------------------------

namespace {
const string getCodePartFromThisFile(const string& marker_tag)
{
  const string path = __FILE__;     // preprocessor macro
  if (!path_exists(path)) {      // utilities.cpp
    return "The selection cuts cannot be displayed because the source file\nselection.cpp used during compilation cannot be opened!\n";
  }
  ifstream ifstr(path.c_str());
  if (!ifstr) {
    fatal_error("Cannot open input file ", path);
  }
  ostringstream oss;
  char ch;
  while (ifstr.get(ch)) {
    oss.put(ch);
  }
  assert(ifstr.eof());
  const string s = oss.str();

  const string marker_begin = "  // MARKER BEGIN " + marker_tag + "\n";
  const string::size_type begin_pos = s.find(marker_begin);
  assert(begin_pos != string::npos);     // begin marker found
  const string marker_end = "  // MARKER END " + marker_tag + "\n";
  const string::size_type end_pos = s.rfind(marker_end);
  assert(end_pos != string::npos);     // end marker found
  const string::size_type begin_length = marker_begin.size();
  const string::size_type code_part_begin = begin_pos + begin_length;
  const int code_part_length = end_pos - code_part_begin;
  assert(code_part_length >= 0);
  const string code_part = s.substr(code_part_begin, code_part_length);
  return code_part;
}
}     // unnamed namespace

// --------------------------------------------------------------------

const bool MinimalCutsSelection::includes(const PhaseSpace& measuredPS, const PhaseSpace& truePS) const
{
  // MARKER BEGIN MinimalCutsSelection
  // amplitude/formcalc/stablepoint.cpp:
#if defined QUADPREC || defined MIXEDPREC
  require(quadruple_precision_sufficient(truePS));
#else
  require(double_precision_sufficient(truePS));
#endif
  require(Mapping::minInvMassPhotonPropCutsPassed(truePS));
  // MARKER END MinimalCutsSelection
  return true;
}

MinimalCutsSelection::MinimalCutsSelection()
  : Selection()
{}

void MinimalCutsSelection::writeSelectionCode(ostream& os) const
{
  os << getCodePartFromThisFile("MinimalCutsSelection");     // helper
}

// --------------------------------------------------------------------

const bool GeneratorCutsSelection::includes(const PhaseSpace& measuredPS, const PhaseSpace& truePS) const
{
  // MARKER BEGIN GeneratorCutsSelection
  // amplitude/formcalc/stablepoint.cpp:
#if defined QUADPREC || defined MIXEDPREC
  require(quadruple_precision_sufficient(truePS));
#else
  require(double_precision_sufficient(truePS));
#endif
  require(Mapping::minInvMassPhotonPropCutsPassed(truePS));
  assert(3.9999 < minInvMassPhotonProp && minInvMassPhotonProp < 4.0001);
#if defined PROCMACRO_WW2l2v || defined PROCMACRO_WWZAZ_2l2v
  const FourMomentum& l1 = measuredPS.aa();
  const FourMomentum& l2 = measuredPS.pp();
  require(min(l1.pt(), l2.pt()) > 3.);
  require(max(fabs(l1.eta()), fabs(l2.eta())) < 2.6);
#elif defined PROCMACRO_ZAZA2l2l || defined PROCMACRO_ZAZA4l
  const FourMomentum& l1 = truePS.aa();
  const FourMomentum& l2 = truePS.ap();
  const FourMomentum& l3 = truePS.pa();
  const FourMomentum& l4 = truePS.pp();
    require(max(l1.pt(), l2.pt(), l3.pt(), l4.pt()) > 20.);
    require(min(l1.pt(), l2.pt(), l3.pt(), l4.pt()) > 5.);
    require(NextToLeading(l1.pt(), l2.pt(), l3.pt(), l4.pt()) > 10. );
  require(max(fabs(l1.eta()), fabs(l2.eta())) < 2.5);
  require(max(fabs(l3.eta()), fabs(l4.eta())) < 2.4);

  ConstValue<FourMomentum> ll1(l1 + l2);
  require(ll1().m() > 4.);
  ConstValue<FourMomentum> ll2(l3 + l4);
  require(ll2().m() > 4.);

  const FourMomentum& ml1 = measuredPS.aa();
  const FourMomentum& ml2 = measuredPS.ap();
  const FourMomentum& ml3 = measuredPS.pa();
  const FourMomentum& ml4 = measuredPS.pp();

  ConstValue<FourMomentum> llll(ml1 + ml2 + ml3 + ml4);
  require(llll().m() > 100.);

#elif defined PROCMACRO_ZAZ_2l2v
  const ConstValue<FourMomentum>& l1 = truePS.aa;
  const ConstValue<FourMomentum>& l2 = truePS.ap;
  require(min(l1().pt(), l2().pt()) > 20.);

  const ConstValue<FourMomentum>& n1 = truePS.pa();
  const ConstValue<FourMomentum>& n2 = truePS.pp();
  //require(min(n1().pt(), n2().pt()) > 5.);

  require(max(fabs(l1().eta()), fabs(l2().eta())) < 2.5); //, fabs(n1().eta()), fabs(n2().eta())) < 2.5);

  ConstValue<FourMomentum> ll(l1() + l2());
  require(ll().m() > 40.);

  //ConstValue<FourMomentum> nn(n1() + n2());
  //require(ll().m() > 5.);

  //require(truePS.p().m() > 60.);

    //ConstValue<FourMomentum> llll(l1() + l2() + n1() + n2());
    //require(llll().m() > 170.); 

#endif    
  // MARKER END GeneratorCutsSelection
  return true;
}

GeneratorCutsSelection::GeneratorCutsSelection()
  : Selection()
{}

const float GeneratorCutsSelection::NextToLeading(const float pt1, const float pt2, const float pt3, const float pt4) const
{
    if (pt1 == max(pt1,pt2,pt3,pt4))
        return max(pt2,pt3,pt4);
    else if (pt2 == max(pt1,pt2,pt3,pt4))
        return max(pt1,pt3,pt4);
    else if (pt3 == max(pt1,pt2,pt3,pt4))
        return max(pt1,pt2,pt4);
    else return max(pt1,pt2,pt3);
}

void GeneratorCutsSelection::writeSelectionCode(ostream& os) const
{
  os << getCodePartFromThisFile("GeneratorCutsSelection");     // helper
}

// --------------------------------------------------------------------

const bool StandardCutsSelection::includes(const PhaseSpace& measuredPS, const PhaseSpace& truePS) const
{
  // MARKER BEGIN StandardCutsSelection
  // amplitude/formcalc/stablepoint.cpp:
#if defined QUADPREC || defined MIXEDPREC
  require(quadruple_precision_sufficient(truePS));
#else
  require(double_precision_sufficient(truePS));
#endif
  require(Mapping::minInvMassPhotonPropCutsPassed(truePS));

  const ConstValue<FourMomentum>& l1 = measuredPS.aa;
  const ConstValue<FourMomentum>& l2 = measuredPS.pp;
  ConstValue<FourMomentum> pmiss(measuredPS.ap() + measuredPS.pa());

  require(min(l1().pt(), l2().pt()) > 20.);
  require(max(fabs(l1().eta()), fabs(l2().eta())) < 2.5);
  require(pmiss().pt() > 30.);

  ConstValue<FourMomentum> ll(l1() + l2());
  require(ll().m() > 12.);
  // MARKER END StandardCutsSelection
  return true;
}

StandardCutsSelection::StandardCutsSelection()
  : Selection()
{}

void StandardCutsSelection::writeSelectionCode(ostream& os) const
{
  os << getCodePartFromThisFile("StandardCutsSelection");     // helper
}

// --------------------------------------------------------------------

const bool HiggsSearchSelection::includes(const PhaseSpace& measuredPS, const PhaseSpace& truePS) const
{
  // MARKER BEGIN HiggsSearchSelection
  // amplitude/formcalc/stablepoint.cpp:
#if defined QUADPREC || defined MIXEDPREC
  require(quadruple_precision_sufficient(truePS));
#else
  require(double_precision_sufficient(truePS));
#endif
  require(Mapping::minInvMassPhotonPropCutsPassed(truePS));

  const ConstValue<FourMomentum>& l1 = truePS.aa;
  const ConstValue<FourMomentum>& l2 = truePS.pp;
  require(min(l1().pt(), l2().pt()) > 20.);
  require(max(fabs(l1().eta()), fabs(l2().eta())) < 2.5);

  ConstValue<FourMomentum> pmiss(truePS.ap()+truePS.pa());
  require(pmiss().pt() > 30.);

  ConstValue<FourMomentum> ll(l1() + l2());
  require(12. < ll().m() && ll().m() < 50.);

  const FourMomentumPair llpair(l1, l2);
  require(llpair.Dphi() < 1.8);
  // MARKER END HiggsSearchSelection
  return true;
}

HiggsSearchSelection::HiggsSearchSelection()
  : Selection()
{}

void HiggsSearchSelection::writeSelectionCode(ostream& os) const
{
  os << getCodePartFromThisFile("HiggsSearchSelection");     // helper
}

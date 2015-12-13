#include <cmath>
using std::log;
using std::sqrt;
#include <cstdlib>
using std::getenv;
#include <fstream>
using std::ifstream;

#include "settings.h"
#include "phasespace.h"
#include "amp_formcalc.h"
#include "pdf.h"
#include "pi.h"
#include "small_powers.h"
#include "rng.h"
#include "utilities.h"

// Description =============================================================

const string settingsDescription = "LHC Higgs Cross Section WG, arXiv:1101.0593 [hep-ph], App. A, G_\\mu scheme";

// MC INTEGRATION MODE =====================================================

int maxNumberOfGridAdaptationIterations = 5;     // zero: skip grid adaptation
const Int64 gridAdaptationNumberOfShots = 10000; // before 50000
const Int64 mainRunNumberOfShots = 100000;       // zero: skip main run (before was 100000long int seedForRandomNumberGenerator = 1;
long int seedForRandomNumberGenerator = 1;
// EVENT GENERATION MODE ===================================================

const bool generateEvents = false;
const bool unweightEvents = true;
const int numberOfEvents = 10000;

#if defined PROCMACRO_WWZAZ_2l2v
const int wwzaz_2l2v_intermediate_state_in_events = W_PAIR;   // W_PAIR or Z_PAIR
#endif

// STANDARD MODEL PARAMETERS ===============================================

const double readHiggsMassFromFile()
{
  ifstream fin("higgsmass");
  double higgsmass = 0.;
  if (fin) {
    higgsmass = getObj<double>(fin);
    fin.close();
  }
  else {
    fatal_error("error: could not read file 'higgsmass'");
  }
  return higgsmass;
}

// input parameters
const double SMP::gfermi = 1.16637e-5;
const double SMP::mW = 80.398;
const double SMP::mZ = 91.1876;
// derived parameters
const double SMP::sin2w = 1. - pow2(SMP::mW/SMP::mZ);
const double SMP::alpha = sqrt(2.)/Pi*SMP::gfermi*pow2(SMP::mW)*SMP::sin2w;
// W and Z widths
const double SMP::wW = 2.08872;   // NLO
const double SMP::wZ = 2.49595;   // NLO
// unstable lepton masses (for event generation only)
const double SMP::mmuon = 0.1057;
const double SMP::mtau = 1.777;
// massive quark loop amplitude parameters
const double SMP::mt = 172.5;
const double SMP::mb = 4.75;
// Higgs mass and width
const double SMP::mH = 125.5;
//const double SMP::mH = readHiggsMassFromFile();
const double SMP::wH = SMP::gamma_H_HTO();
// Higgs rescaling relative to SM couplings 
const double BSM::HggCouplingRescalingFactor = 1.;
const double BSM::HWWCouplingRescalingFactor = 1.;
const double BSM::HZZCouplingRescalingFactor = 1.;

// process dependent vector boson mass and width
#if defined PROCMACRO_WW2l2v
const double SMP::mV = SMP::mW;
const double SMP::wV = SMP::wW;
#endif
#if defined PROCMACRO_ZAZA2l2l
const double SMP::mV = SMP::mZ;
const double SMP::wV = SMP::wZ;
#endif
#if defined PROCMACRO_ZAZ_2l2v
const double SMP::mV = SMP::mZ;
const double SMP::wV = SMP::wZ;
#endif
#if defined PROCMACRO_ZAZA4l
const double SMP::mV = SMP::mZ;
const double SMP::wV = SMP::wZ;
#endif
#if defined PROCMACRO_WWZAZ_2l2v
// not defined because of V=W,Z ambiguity
#endif

// HIGGS PROPAGATOR ========================================================

const bool applyingHiggsNWA = false;   // narrow(zero)-width approximation

//const int bwshape = 1;   // running-width Breit-Wigner
//const int bwshape = 2;   // fixed-width Breit-Wigner
const int bwshape = 3;   // Passarino OFFP prescription (complex pole scheme)

// AMPLITUDE ===============================================================

const unsigned int amplitude_selection = (
+AmpSqr::CONTINUUMBACKGROUND
+AmpSqr::HIGGSSIGNAL

);

// SCALES ==================================================================

vector<Scales*> scales_ptr_vector;

namespace {
RootshatScales central_scale_1(0.5);
RootshatScales central_scale_2(1.);
RootshatScales central_scale_3(0.25);
FixedScales central_scale_4(SMP::mH/2.);
FixedScales central_scale_5(SMP::mV);
}     // unnamed namespace

void initialize_scales()
{
  scales_ptr_vector.push_back(&central_scale_1);
  scales_ptr_vector.push_back(&central_scale_2);
  scales_ptr_vector.push_back(&central_scale_3);
  scales_ptr_vector.push_back(&central_scale_4);
  scales_ptr_vector.push_back(&central_scale_5);
}

const int scale_variation_factor = 2;

// SELECTION CUTS ==========================================================

const double minInvMassPhotonProp = 4.;     // GeV

//namespace { const MinimalCutsSelection _selection; }
namespace { const GeneratorCutsSelection _selection; }
//namespace { const StandardCutsSelection _selection; }
//namespace { const HiggsSearchSelection _selection; }

const Selection& selection = _selection;

// COLLIDER ================================================================

//ConstValue<double> PhaseSpace::E_CMS = 8.e3;     // GeV
ConstValue<double> PhaseSpace::E_CMS = 13.e3;     // GeV
//ConstValue<double> PhaseSpace::E_CMS = 14.e3;     // GeV

// PARTON DISTRIBUTION FUNCTIONS and ALPHAS ================================
//const double alpha_s = 05;
// select PDF set:
//const string lhapdf_file = "MSTW2008lo68cl.LHgrid";
//const string lhapdf_file = "CT14llo.LHgrid";
//const string lhapdf_file = "CT14lo.LHgrid";
const string lhapdf_file = "PDF4LHC15_nlo_mc.LHgrid";
//const string lhapdf_file = "PDF4LHC15_nnlo_mc.LHgrid";
//const string lhapdf_file = "cteq6ll.LHpdf";
//const string lhapdf_file = "MSTW2008nnlo68cl.LHgrid";
//const string lhapdf_file = "CT10nnlo.LHgrid";

// calculate PDF error (if supported by PDF set):
const bool requestPDFError = false;

// using set's default alpha_s (LHAPDF)
// =========================================================================

// Some settings should only be modified by advanced users:
#include "advanced_settings.inc"

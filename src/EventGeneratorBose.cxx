#include "RhoPrime/EventGeneratorBose.h"

#include <TAxis.h>
#include <TFile.h>
#include <TParameter.h>
#include <TTree.h>
#include <TVector3.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <functional>
#include <memory>
#include <queue>
#include <stdexcept>
#include <utility>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace rhoprime {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kTiny = 1.0e-300;

void Require(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error(message);
}

TH2D* CloneRequiredHist(TFile& file, const char* name) {
  auto* h = dynamic_cast<TH2D*>(file.Get(name));
  if (!h) throw std::runtime_error(std::string("missing required TH2D: ") + name);
  auto* clone = dynamic_cast<TH2D*>(h->Clone((std::string(name) + "_bose_clone").c_str()));
  if (!clone) throw std::runtime_error(std::string("failed to clone TH2D: ") + name);
  clone->SetDirectory(nullptr);
  return clone;
}

std::string HistNameForYBin(const char* base, int iy) {
  return std::string(base) + "_y" + std::to_string(iy);
}

double SafeInterpolate(const TH2D* h, double x, double y) {
  if (!h) return 0.0;
  const TAxis* ax = h->GetXaxis();
  const TAxis* ay = h->GetYaxis();
  if (x < ax->GetXmin() || x > ax->GetXmax() || y < ay->GetXmin() || y > ay->GetXmax()) {
    return 0.0;
  }
  const double v = h->Interpolate(x, y);
  return std::isfinite(v) ? v : 0.0;
}

double TwoBodyMomentum(double parent_m, double m1, double m2) {
  if (!(parent_m > m1 + m2) || !(parent_m > 0.0)) return 0.0;
  const double m2p = parent_m * parent_m;
  const double plus = m1 + m2;
  const double minus = m1 - m2;
  return std::sqrt(std::max(0.0, (m2p - plus * plus) * (m2p - minus * minus))) /
         (2.0 * parent_m);
}

void DecayTwoBodyIsotropic(const TLorentzVector& parent, double m1, double m2, TLorentzVector& p1,
                           TLorentzVector& p2, TRandom3& rng) {
  const double p = TwoBodyMomentum(parent.M(), m1, m2);
  const double cos_theta = rng.Uniform(-1.0, 1.0);
  const double sin_theta = std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));
  const double phi = rng.Uniform(-kPi, kPi);
  const double px = p * sin_theta * std::cos(phi);
  const double py = p * sin_theta * std::sin(phi);
  const double pz = p * cos_theta;
  p1.SetPxPyPzE(px, py, pz, std::sqrt(p * p + m1 * m1));
  p2.SetPxPyPzE(-px, -py, -pz, std::sqrt(p * p + m2 * m2));
  p1.Boost(parent.BoostVector());
  p2.Boost(parent.BoostVector());
}

double CauchyCdf(double x, double loc, double scale) {
  return std::atan((x - loc) / scale) / kPi + 0.5;
}

double CauchyPdf(double x, double loc, double scale) {
  const double u = (x - loc) / scale;
  return 1.0 / (kPi * scale * (1.0 + u * u));
}

double SampleTruncatedCauchy(TRandom3& rng, double loc, double width, double lo, double hi) {
  Require(hi > lo, "invalid truncated Cauchy interval");
  const double scale = std::max(width * 0.5, 1.0e-12);
  const double c_lo = CauchyCdf(lo, loc, scale);
  const double c_hi = CauchyCdf(hi, loc, scale);
  const double u = rng.Uniform(c_lo, c_hi);
  return loc + scale * std::tan(kPi * (u - 0.5));
}

double TruncatedCauchyPdf(double x, double loc, double width, double lo, double hi) {
  if (x < lo || x > hi || !(hi > lo)) return 0.0;
  const double scale = std::max(width * 0.5, 1.0e-12);
  const double norm = CauchyCdf(hi, loc, scale) - CauchyCdf(lo, loc, scale);
  if (!(norm > 0.0)) return 0.0;
  return CauchyPdf(x, loc, scale) / norm;
}

double DecayMomentum(double m, double daughter_m) {
  if (m <= 2.0 * daughter_m) return 0.0;
  return std::sqrt(std::max(0.0, 0.25 * m * m - daughter_m * daughter_m));
}

double CascadePhaseSpaceWeight(double parent_m, double m_rho, double m_sigma) {
  if (!(parent_m > m_rho + m_sigma)) return 0.0;
  const double p_rho_sigma = TwoBodyMomentum(parent_m, m_rho, m_sigma);
  const double q_rho = DecayMomentum(m_rho, kMPi);
  const double q_sigma = DecayMomentum(m_sigma, kMPi);
  if (!(p_rho_sigma > 0.0) || !(q_rho > 0.0) || !(q_sigma > 0.0)) return 0.0;

  // Constants common to all events are omitted. The retained factors are the
  // variable part of dPhi4 in a sequential two-body parametrization, including
  // dm_rho^2 dm_sigma^2 = 4 m_rho m_sigma dm_rho dm_sigma.
  return (2.0 * m_rho) * (2.0 * m_sigma) * (p_rho_sigma / parent_m) *
         (q_rho / m_rho) * (q_sigma / m_sigma);
}

double MassDependentWidth(double m, double pole, double width0, int L, double daughter_m) {
  if (!(m > 2.0 * daughter_m) || !(pole > 2.0 * daughter_m)) return 0.0;
  const double q = DecayMomentum(m, daughter_m);
  const double q0 = DecayMomentum(pole, daughter_m);
  if (!(q > 0.0) || !(q0 > 0.0)) return 0.0;
  const int power = 2 * L + 1;
  return width0 * std::pow(q / q0, power) * (pole / m);
}

std::complex<double> BreitWigner(double m, double pole, double width0, int L) {
  if (!(m > 2.0 * kMPi)) return {0.0, 0.0};
  const double gamma = MassDependentWidth(m, pole, width0, L, kMPi);
  if (!(gamma > 0.0)) return {0.0, 0.0};
  return 1.0 / std::complex<double>(m * m - pole * pole, pole * gamma);
}

double InvariantMass(const TLorentzVector& a, const TLorentzVector& b) {
  return std::max(0.0, (a + b).M());
}

ComplexVector3 AmpPair(const TLorentzVector& rho_plus, const TLorentzVector& rho_minus,
                       const TLorentzVector& sig_plus, const TLorentzVector& sig_minus) {
  const std::complex<double> bw_rho =
      BreitWigner(InvariantMass(rho_plus, rho_minus), kMRho, kGRho, 1);
  const std::complex<double> bw_sigma =
      BreitWigner(InvariantMass(sig_plus, sig_minus), kMSigma, kGSigma, 0);
  const std::complex<double> coeff = bw_rho * bw_sigma;
  const TVector3 q = rho_plus.Vect() - rho_minus.Vect();
  return {coeff * q.X(), coeff * q.Y(), coeff * q.Z()};
}

double Norm2(const std::complex<double>& z) {
  return std::norm(z);
}

}  // namespace

ComplexVector3 ComplexVector3::operator+(const ComplexVector3& o) const {
  return {x + o.x, y + o.y, z + o.z};
}

ComplexVector3& ComplexVector3::operator+=(const ComplexVector3& o) {
  x += o.x;
  y += o.y;
  z += o.z;
  return *this;
}

void HistSampler2D::Build(TH2D* h) {
  Require(h != nullptr, "cannot build sampler from null histogram");
  hist.reset(dynamic_cast<TH2D*>(h->Clone((std::string(h->GetName()) + "_sampler").c_str())));
  Require(hist != nullptr, "failed to clone histogram sampler");
  hist->SetDirectory(nullptr);

  cdf.clear();
  total_sum = 0.0;
  cdf.reserve(hist->GetNbinsX() * hist->GetNbinsY());
  for (int ix = 1; ix <= hist->GetNbinsX(); ++ix) {
    const double dx = hist->GetXaxis()->GetBinWidth(ix);
    for (int iy = 1; iy <= hist->GetNbinsY(); ++iy) {
      const double dy = hist->GetYaxis()->GetBinWidth(iy);
      const double content = hist->GetBinContent(ix, iy);
      const double w = std::isfinite(content) ? std::max(0.0, content) * dx * dy : 0.0;
      total_sum += w;
      cdf.push_back(total_sum);
    }
  }
  Require(total_sum > 0.0, std::string("histogram has no positive entries: ") + h->GetName());
  for (double& v : cdf) v /= total_sum;
  cdf.back() = 1.0;
}

bool HistSampler2D::Sample(TRandom3& rng, double& x, double& y) const {
  if (!hist || cdf.empty()) return false;
  const double u = rng.Uniform(0.0, 1.0);
  auto it = std::lower_bound(cdf.begin(), cdf.end(), u);
  int flat = static_cast<int>(it - cdf.begin());
  flat = std::clamp(flat, 0, static_cast<int>(cdf.size()) - 1);
  const int ny = hist->GetNbinsY();
  const int ix = flat / ny + 1;
  const int iy = flat % ny + 1;
  x = rng.Uniform(hist->GetXaxis()->GetBinLowEdge(ix), hist->GetXaxis()->GetBinUpEdge(ix));
  y = rng.Uniform(hist->GetYaxis()->GetBinLowEdge(iy), hist->GetYaxis()->GetBinUpEdge(iy));
  return Interpolate(x, y) > 0.0;
}

double HistSampler2D::Interpolate(double x, double y) const {
  return SafeInterpolate(hist.get(), x, y);
}

struct EventGeneratorBose::Candidate {
  double M = 0.0;
  double y = 0.0;
  double px = 0.0;
  double py = 0.0;
  double m_rho = 0.0;
  double m_sigma = 0.0;
  int iy_bin = 0;
  double proposal_density = 0.0;
  double phase_space_weight = 0.0;
  double pol_xx = 0.0;
  double pol_yy = 0.0;
  double pol_rexy = 0.0;
  double pol_imxy = 0.0;
  TLorentzVector rhoprime_lab;
  TLorentzVector rho_rest;
  TLorentzVector sigma_rest;
  TLorentzVector rho_lab;
  TLorentzVector sigma_lab;
  TLorentzVector p1_rest;
  TLorentzVector p2_rest;
  TLorentzVector p3_rest;
  TLorentzVector p4_rest;
  TLorentzVector p1_lab;
  TLorentzVector p2_lab;
  TLorentzVector p3_lab;
  TLorentzVector p4_lab;
};

EventGeneratorBose::EventGeneratorBose(BoseGeneratorConfig cfg) : cfg_(std::move(cfg)) {}

void EventGeneratorBose::LoadHistograms() {
  std::unique_ptr<TFile> file(TFile::Open(cfg_.input_root.c_str(), "READ"));
  Require(file && !file->IsZombie(), "cannot open input ROOT file: " + cfg_.input_root);

  std::unique_ptr<TH2D> h_my(CloneRequiredHist(*file, "h2_dM_dy"));
  cross_section_mb_ = h_my->Integral("width");
  Require(cross_section_mb_ > 0.0, "h2_dM_dy has no positive width integral");
  sampler_M_y_.Build(h_my.get());

  auto* n_y_param = dynamic_cast<TParameter<int>*>(file->Get("n_rapidity_bins"));
  n_y_bins_ = n_y_param ? n_y_param->GetVal() : 1;
  Require(n_y_bins_ > 0, "n_rapidity_bins must be positive");

  sampler_px_py_.clear();
  h_pol_xx_.clear();
  h_pol_yy_.clear();
  h_pol_rexy_.clear();
  h_pol_imxy_.clear();
  sampler_px_py_.resize(n_y_bins_);
  h_pol_xx_.reserve(n_y_bins_);
  h_pol_yy_.reserve(n_y_bins_);
  h_pol_rexy_.reserve(n_y_bins_);
  h_pol_imxy_.reserve(n_y_bins_);

  for (int iy = 0; iy < n_y_bins_; ++iy) {
    const bool use_legacy_names = (n_y_param == nullptr);
    const std::string sigma_name =
        use_legacy_names ? "h2_px_py_sigma" : HistNameForYBin("h2_px_py_sigma", iy);
    const std::string ax2_name =
        use_legacy_names ? "h2_px_py_Ax2" : HistNameForYBin("h2_px_py_Ax2", iy);
    const std::string ay2_name =
        use_legacy_names ? "h2_px_py_Ay2" : HistNameForYBin("h2_px_py_Ay2", iy);
    const std::string re_name =
        use_legacy_names ? "h2_px_py_ReAxAy" : HistNameForYBin("h2_px_py_ReAxAy", iy);
    const std::string im_name =
        use_legacy_names ? "h2_px_py_ImAxAy" : HistNameForYBin("h2_px_py_ImAxAy", iy);

    std::unique_ptr<TH2D> h_pxpy(CloneRequiredHist(*file, sigma_name.c_str()));
    sampler_px_py_[iy].Build(h_pxpy.get());
    h_pol_xx_.emplace_back(CloneRequiredHist(*file, ax2_name.c_str()));
    h_pol_yy_.emplace_back(CloneRequiredHist(*file, ay2_name.c_str()));
    h_pol_rexy_.emplace_back(CloneRequiredHist(*file, re_name.c_str()));
    h_pol_imxy_.emplace_back(CloneRequiredHist(*file, im_name.c_str()));
  }
}

void EventGeneratorBose::Initialize() {
  Require(!cfg_.input_root.empty(), "input ROOT path is required");
  Require(!cfg_.output_root.empty(), "output ROOT path is required");
  Require(cfg_.n_events > 0, "nEvents must be positive");
  Require(cfg_.decay_norm_trials_per_mass > 0, "decay_norm_trials_per_mass must be positive");
  Require(cfg_.unweighting_mode >= 1 && cfg_.unweighting_mode <= 3,
          "unweighting_mode must be 1 (accept-reject), 2 (resample), or 3 (reservoir)");
  Require(cfg_.output_compression_level >= 0 && cfg_.output_compression_level <= 9,
          "output_compression_level must be in [0,9]");
  Require(cfg_.mode1_threads > 0, "mode1_threads must be positive");
  LoadHistograms();
  BuildDecayWeightNormalization();
  initialized_ = true;
}

EventGeneratorBose::Candidate EventGeneratorBose::MakeCandidate(TRandom3& rng) const {
  const double min_pair_m = 2.0 * kMPi;

  for (int tries = 0; tries < 10000; ++tries) {
    Candidate c;
    if (!sampler_M_y_.Sample(rng, c.M, c.y)) continue;
    if (c.M <= 2.0 * min_pair_m) continue;
    c.iy_bin = std::clamp(sampler_M_y_.hist->GetYaxis()->FindBin(c.y) - 1, 0, n_y_bins_ - 1);
    if (!sampler_px_py_[c.iy_bin].Sample(rng, c.px, c.py)) continue;
    const double px_lookup = c.px;
    const double py_lookup = c.py;

    const double mr_lo = min_pair_m;
    const double mr_hi = c.M - min_pair_m;
    c.m_rho = SampleTruncatedCauchy(rng, kMRho, kGRho, mr_lo, mr_hi);
    const double ms_lo = min_pair_m;
    const double ms_hi = c.M - c.m_rho;
    if (!(ms_hi > ms_lo)) continue;
    c.m_sigma = SampleTruncatedCauchy(rng, kMSigma, kGSigma, ms_lo, ms_hi);
    if (c.m_rho + c.m_sigma >= c.M) continue;

    const double mt = std::sqrt(c.M * c.M + c.px * c.px + c.py * c.py);
    c.rhoprime_lab.SetPxPyPzE(c.px, c.py, mt * std::sinh(c.y), mt * std::cosh(c.y));

    const TLorentzVector rhoprime_rest(0.0, 0.0, 0.0, c.M);
    DecayTwoBodyIsotropic(rhoprime_rest, c.m_rho, c.m_sigma, c.rho_rest, c.sigma_rest, rng);
    DecayTwoBodyIsotropic(c.rho_rest, kMPi, kMPi, c.p1_rest, c.p2_rest, rng);
    DecayTwoBodyIsotropic(c.sigma_rest, kMPi, kMPi, c.p3_rest, c.p4_rest, rng);

    const double pdf_mr = TruncatedCauchyPdf(c.m_rho, kMRho, kGRho, mr_lo, mr_hi);
    const double pdf_ms = TruncatedCauchyPdf(c.m_sigma, kMSigma, kGSigma, ms_lo, ms_hi);
    c.proposal_density = pdf_mr * pdf_ms;
    if (!(c.proposal_density > kTiny) || !std::isfinite(c.proposal_density)) continue;
    c.phase_space_weight = CascadePhaseSpaceWeight(c.M, c.m_rho, c.m_sigma);
    if (!(c.phase_space_weight > 0.0) || !std::isfinite(c.phase_space_weight)) continue;

    c.pol_xx = SafeInterpolate(h_pol_xx_[c.iy_bin].get(), px_lookup, py_lookup);
    c.pol_yy = SafeInterpolate(h_pol_yy_[c.iy_bin].get(), px_lookup, py_lookup);
    c.pol_rexy = SafeInterpolate(h_pol_rexy_[c.iy_bin].get(), px_lookup, py_lookup);
    c.pol_imxy = SafeInterpolate(h_pol_imxy_[c.iy_bin].get(), px_lookup, py_lookup);

    c.rho_lab = c.rho_rest;
    c.sigma_lab = c.sigma_rest;
    c.p1_lab = c.p1_rest;
    c.p2_lab = c.p2_rest;
    c.p3_lab = c.p3_rest;
    c.p4_lab = c.p4_rest;
    const TVector3 beta = c.rhoprime_lab.BoostVector();
    c.rho_lab.Boost(beta);
    c.sigma_lab.Boost(beta);
    c.p1_lab.Boost(beta);
    c.p2_lab.Boost(beta);
    c.p3_lab.Boost(beta);
    c.p4_lab.Boost(beta);
    return c;
  }
  throw std::runtime_error("failed to create a physical rho-prime Bose candidate");
}

double EventGeneratorBose::RawCandidateWeight(const Candidate& c) const {
  ComplexVector3 mtot;
  // A: (1,2) -> rho, (3,4) -> sigma
  mtot += AmpPair(c.p1_rest, c.p2_rest, c.p3_rest, c.p4_rest);
  if (cfg_.bose_symmetrize) {
    // Both pi+ and both pi- labels are interchangeable. The rho and sigma are
    // distinct intermediate states, so these assignments enter with the same sign.
    // B: (3,2) -> rho, (1,4) -> sigma
    mtot += AmpPair(c.p3_rest, c.p2_rest, c.p1_rest, c.p4_rest);
    // C: (1,4) -> rho, (3,2) -> sigma
    mtot += AmpPair(c.p1_rest, c.p4_rest, c.p3_rest, c.p2_rest);
    // D: (3,4) -> rho, (1,2) -> sigma
    mtot += AmpPair(c.p3_rest, c.p4_rest, c.p1_rest, c.p2_rest);
  }

  const std::complex<double> rho_xy(c.pol_rexy, c.pol_imxy);
  // Normalize SDM by total intensity to get dimensionless spin density matrix
  const double I_tot = sampler_px_py_[c.iy_bin].Interpolate(c.px, c.py);
  if (!(I_tot > 0.0)) return 0.0;
  const double rho_xx = c.pol_xx / I_tot;
  const double rho_yy = c.pol_yy / I_tot;
  const std::complex<double> rho_xy_norm = rho_xy / I_tot;
  const double intensity =
      rho_xx * Norm2(mtot.x) + rho_yy * Norm2(mtot.y) +
      2.0 * std::real(rho_xy_norm * std::conj(mtot.y) * mtot.x);

  if (!std::isfinite(intensity) || intensity <= 0.0) return 0.0;
  const double w = intensity * c.phase_space_weight / c.proposal_density;
  return std::isfinite(w) ? std::max(0.0, w) : 0.0;
}

double EventGeneratorBose::CandidateWeight(const Candidate& c) const {
  const double raw = RawCandidateWeight(c);
  if (!(raw > 0.0)) return 0.0;
  if (decay_weight_norm_by_mass_bin_.empty() || !sampler_M_y_.hist) return raw;
  const TAxis* ax = sampler_M_y_.hist->GetXaxis();
  int bin = std::clamp(ax->FindBin(c.M), 1, ax->GetNbins());
  double norm = decay_weight_norm_by_mass_bin_[bin - 1];
  const double center = ax->GetBinCenter(bin);
  if (c.M > center && bin < ax->GetNbins()) {
    const double next_center = ax->GetBinCenter(bin + 1);
    const double t = (c.M - center) / (next_center - center);
    norm = norm * (1.0 - t) + decay_weight_norm_by_mass_bin_[bin] * t;
  } else if (c.M < center && bin > 1) {
    const double prev_center = ax->GetBinCenter(bin - 1);
    const double t = (c.M - prev_center) / (center - prev_center);
    norm = decay_weight_norm_by_mass_bin_[bin - 2] * (1.0 - t) + norm * t;
  }
  if (!(norm > 0.0) || !std::isfinite(norm)) return 0.0;
  return raw / norm;
}

double EventGeneratorBose::EstimateDecayWeightNorm(double M, TRandom3& rng,
                                                   long long trials) const {
  const double min_pair_m = 2.0 * kMPi;
  if (!(M > 2.0 * min_pair_m)) return 1.0;

  double sum = 0.0;
  long long accepted = 0;
  for (long long i = 0; i < trials; ++i) {
    Candidate c;
    c.M = M;

    const double mr_lo = min_pair_m;
    const double mr_hi = c.M - min_pair_m;
    c.m_rho = SampleTruncatedCauchy(rng, kMRho, kGRho, mr_lo, mr_hi);
    const double ms_lo = min_pair_m;
    const double ms_hi = c.M - c.m_rho;
    if (!(ms_hi > ms_lo)) continue;
    c.m_sigma = SampleTruncatedCauchy(rng, kMSigma, kGSigma, ms_lo, ms_hi);
    if (c.m_rho + c.m_sigma >= c.M) continue;

    const TLorentzVector rhoprime_rest(0.0, 0.0, 0.0, c.M);
    DecayTwoBodyIsotropic(rhoprime_rest, c.m_rho, c.m_sigma, c.rho_rest, c.sigma_rest, rng);
    DecayTwoBodyIsotropic(c.rho_rest, kMPi, kMPi, c.p1_rest, c.p2_rest, rng);
    DecayTwoBodyIsotropic(c.sigma_rest, kMPi, kMPi, c.p3_rest, c.p4_rest, rng);

    const double pdf_mr = TruncatedCauchyPdf(c.m_rho, kMRho, kGRho, mr_lo, mr_hi);
    const double pdf_ms = TruncatedCauchyPdf(c.m_sigma, kMSigma, kGSigma, ms_lo, ms_hi);
    c.proposal_density = pdf_mr * pdf_ms;
    c.phase_space_weight = CascadePhaseSpaceWeight(c.M, c.m_rho, c.m_sigma);
    if (!(c.proposal_density > kTiny) || !(c.phase_space_weight > 0.0)) continue;

    ComplexVector3 mtot;
    mtot += AmpPair(c.p1_rest, c.p2_rest, c.p3_rest, c.p4_rest);
    if (cfg_.bose_symmetrize) {
      mtot += AmpPair(c.p3_rest, c.p2_rest, c.p1_rest, c.p4_rest);
      mtot += AmpPair(c.p1_rest, c.p4_rest, c.p3_rest, c.p2_rest);
      mtot += AmpPair(c.p3_rest, c.p4_rest, c.p1_rest, c.p2_rest);
    }

    const double unpolarized_transverse_intensity = 0.5 * (Norm2(mtot.x) + Norm2(mtot.y));
    const double w =
        unpolarized_transverse_intensity * c.phase_space_weight / c.proposal_density;
    if (std::isfinite(w) && w > 0.0) {
      sum += w;
      ++accepted;
    }
  }

  Require(accepted > 0, "failed to estimate Bose decay normalization");
  const double norm = sum / static_cast<double>(accepted);
  Require(norm > 0.0 && std::isfinite(norm), "invalid Bose decay normalization");
  return norm;
}

void EventGeneratorBose::BuildDecayWeightNormalization() {
  Require(sampler_M_y_.hist != nullptr, "mass sampler must be initialized");
  const int n_mass_bins = sampler_M_y_.hist->GetNbinsX();
  decay_weight_norm_by_mass_bin_.assign(n_mass_bins, 1.0);

  TRandom3 rng(cfg_.random_seed + 0x9e3779b9U);
  for (int i = 1; i <= n_mass_bins; ++i) {
    const double M = sampler_M_y_.hist->GetXaxis()->GetBinCenter(i);
    decay_weight_norm_by_mass_bin_[i - 1] =
        EstimateDecayWeightNorm(M, rng, cfg_.decay_norm_trials_per_mass);
  }
}

double EventGeneratorBose::EstimateUnweightingEnvelope(TRandom3& rng, long long trials) {
  Require(trials > 0, "unweighting_trials must be positive");
  double max_weight = 0.0;
  long long positive = 0;
  for (long long i = 0; i < trials; ++i) {
    Candidate c;
    try {
      c = MakeCandidate(rng);
    } catch (const std::exception&) {
      continue;
    }
    const double w = CandidateWeight(c);
    if (std::isfinite(w) && w > 0.0) {
      max_weight = std::max(max_weight, w);
      ++positive;
    }
  }
  Require(positive > 0 && max_weight > 0.0, "failed to estimate unweighting envelope");
  const double safety = std::max(1.0, cfg_.unweighting_safety_factor);
  return max_weight * safety;
}

BoseGenerationSummary EventGeneratorBose::Generate() {
  if (!initialized_) Initialize();
  TRandom3 rng(cfg_.random_seed);

  long long rejected = 0;
  long long overweight_candidates = 0;
  double unweighting_envelope = 0.0;
  int mode1_threads_used = 1;
  if (cfg_.unweighted_events && cfg_.unweighting_mode == 1) {
    TRandom3 envelope_rng(cfg_.random_seed + 0x85ebca6bU);
    unweighting_envelope = EstimateUnweightingEnvelope(envelope_rng, cfg_.unweighting_trials);
    std::cout << "Unweighted accept-reject mode: envelope=" << unweighting_envelope
              << " from trials=" << cfg_.unweighting_trials
              << " safety=" << cfg_.unweighting_safety_factor << "\n";
  } else if (cfg_.unweighted_events && cfg_.unweighting_mode == 2) {
    std::cout << "Fast resampled unit-weight mode: pool_trials=" << cfg_.unweighting_trials
              << "\n";
  } else if (cfg_.unweighted_events && cfg_.unweighting_mode == 3) {
    std::cout << "Weighted-reservoir unit-weight mode: pool_trials=" << cfg_.unweighting_trials
              << " target_events=" << cfg_.n_events << "\n";
  }

  std::unique_ptr<TFile> out(TFile::Open(cfg_.output_root.c_str(), "RECREATE"));
  Require(out && !out->IsZombie(), "cannot create output ROOT file: " + cfg_.output_root);
  out->SetCompressionLevel(cfg_.output_compression_level);

  TTree tree("Events", "Bose-symmetrized rho-prime -> rho0 sigma -> pi+ pi- pi+ pi-");

  double event_weight = 0.0;
  double M_rhoprime = 0.0;
  double y_rhoprime = 0.0;
  double pt_rhoprime = 0.0;
  double m_rho = 0.0;
  double m_sigma = 0.0;
  double pol_xx = 0.0;
  double pol_yy = 0.0;
  double pol_rexy = 0.0;
  double pol_imxy = 0.0;
  TLorentzVector pi1;
  TLorentzVector pi2;
  TLorentzVector pi3;
  TLorentzVector pi4;
  TLorentzVector rho_momentum;
  TLorentzVector sigma_momentum;
  int pi_parent[4] = {1, 1, 2, 2};
  int pi_charge[4] = {1, -1, 1, -1};

  tree.Branch("event_weight", &event_weight, "event_weight/D");
  tree.Branch("M_rhoprime", &M_rhoprime, "M_rhoprime/D");
  tree.Branch("y_rhoprime", &y_rhoprime, "y_rhoprime/D");
  tree.Branch("pt_rhoprime", &pt_rhoprime, "pt_rhoprime/D");
  tree.Branch("m_rho", &m_rho, "m_rho/D");
  tree.Branch("m_sigma", &m_sigma, "m_sigma/D");
  tree.Branch("pol_xx", &pol_xx, "pol_xx/D");
  tree.Branch("pol_yy", &pol_yy, "pol_yy/D");
  tree.Branch("pol_rexy", &pol_rexy, "pol_rexy/D");
  tree.Branch("pol_imxy", &pol_imxy, "pol_imxy/D");
  tree.Branch("pi1", &pi1);
  tree.Branch("pi2", &pi2);
  tree.Branch("pi3", &pi3);
  tree.Branch("pi4", &pi4);
  tree.Branch("rho_momentum", &rho_momentum);
  tree.Branch("sigma_momentum", &sigma_momentum);
  tree.Branch("pi_parent", pi_parent, "pi_parent[4]/I");
  tree.Branch("pi_charge", pi_charge, "pi_charge[4]/I");

  auto fill_candidate = [&](const Candidate& c, double weight_to_store) {
    event_weight = weight_to_store;
    M_rhoprime = c.M;
    y_rhoprime = c.y;
    pt_rhoprime = std::sqrt(c.px * c.px + c.py * c.py);
    m_rho = c.m_rho;
    m_sigma = c.m_sigma;
    pol_xx = c.pol_xx;
    pol_yy = c.pol_yy;
    pol_rexy = c.pol_rexy;
    pol_imxy = c.pol_imxy;
    pi1 = c.p1_lab;
    pi2 = c.p2_lab;
    pi3 = c.p3_lab;
    pi4 = c.p4_lab;
    rho_momentum = c.rho_lab;
    sigma_momentum = c.sigma_lab;

    const TLorentzVector residual = c.rhoprime_lab - pi1 - pi2 - pi3 - pi4;
    if (std::abs(residual.E()) > 1.0e-8 || residual.Vect().Mag() > 1.0e-8) {
      throw std::runtime_error("internal four-momentum conservation failure");
    }
    tree.Fill();
  };

  long long filled = 0;
  if (cfg_.unweighted_events && cfg_.unweighting_mode == 2) {
    std::vector<Candidate> pool;
    std::vector<double> cdf;
    pool.reserve(static_cast<size_t>(std::min<long long>(cfg_.unweighting_trials, 10000000LL)));
    cdf.reserve(pool.capacity());
    double sum_w = 0.0;
    for (long long trials = 0; trials < cfg_.unweighting_trials; ++trials) {
      Candidate c;
      try {
        c = MakeCandidate(rng);
      } catch (const std::exception&) {
        ++rejected;
        continue;
      }
      const double w = CandidateWeight(c);
      if (!(w > 0.0) || !std::isfinite(w)) {
        ++rejected;
        continue;
      }
      pool.push_back(c);
      sum_w += w;
      cdf.push_back(sum_w);
    }
    Require(!pool.empty() && sum_w > 0.0, "resampling pool has no positive-weight candidates");
    for (; filled < cfg_.n_events; ++filled) {
      const double u = rng.Uniform(0.0, sum_w);
      auto it = std::lower_bound(cdf.begin(), cdf.end(), u);
      const size_t idx = static_cast<size_t>(std::clamp<long long>(it - cdf.begin(), 0,
                                            static_cast<long long>(pool.size()) - 1));
      fill_candidate(pool[idx], 1.0);
    }
  } else if (cfg_.unweighted_events && cfg_.unweighting_mode == 3) {
    struct ReservoirItem {
      double key = -std::numeric_limits<double>::infinity();
      Candidate candidate;
      bool operator>(const ReservoirItem& other) const { return key > other.key; }
    };
    std::priority_queue<ReservoirItem, std::vector<ReservoirItem>, std::greater<ReservoirItem>> reservoir;
    long long positive = 0;
    for (long long trials = 0; trials < cfg_.unweighting_trials; ++trials) {
      Candidate c;
      try {
        c = MakeCandidate(rng);
      } catch (const std::exception&) {
        ++rejected;
        continue;
      }
      const double w = CandidateWeight(c);
      if (!(w > 0.0) || !std::isfinite(w)) {
        ++rejected;
        continue;
      }
      ++positive;
      const double u = std::max(rng.Uniform(0.0, 1.0), std::numeric_limits<double>::min());
      const double key = std::log(u) / w;  // Efraimidis-Spirakis weighted reservoir key.
      ReservoirItem item{key, c};
      if (static_cast<long long>(reservoir.size()) < cfg_.n_events) {
        reservoir.push(std::move(item));
      } else if (key > reservoir.top().key) {
        reservoir.pop();
        reservoir.push(std::move(item));
      }
    }
    Require(positive >= cfg_.n_events && static_cast<long long>(reservoir.size()) == cfg_.n_events,
            "reservoir mode needs unweighting_trials >= requested events with positive weights");
    std::vector<Candidate> selected;
    selected.reserve(static_cast<size_t>(cfg_.n_events));
    while (!reservoir.empty()) {
      selected.push_back(reservoir.top().candidate);
      reservoir.pop();
    }
    std::reverse(selected.begin(), selected.end());
    for (const auto& c : selected) {
      fill_candidate(c, 1.0);
      ++filled;
    }
  } else if (cfg_.unweighted_events && cfg_.unweighting_mode == 1 && cfg_.mode1_threads > 1) {
#ifdef _OPENMP
    const int n_threads = std::max(1, std::min(cfg_.mode1_threads, omp_get_max_threads()));
    mode1_threads_used = n_threads;
    const long long trials_per_batch = std::max<long long>(262144LL, 16384LL * n_threads);
    std::vector<TRandom3> thread_rngs;
    thread_rngs.reserve(n_threads);
    for (int tid = 0; tid < n_threads; ++tid) {
      thread_rngs.emplace_back(cfg_.random_seed + 0xc2b2ae35U * static_cast<unsigned int>(tid + 1));
    }
    std::cout << "Parallel exact accept-reject: threads=" << n_threads
              << " trials_per_batch=" << trials_per_batch << "\n";

    while (filled < cfg_.n_events) {
      std::vector<std::vector<Candidate>> accepted_by_thread(n_threads);
      std::vector<long long> rejected_by_thread(n_threads, 0);
      std::vector<long long> overweight_by_thread(n_threads, 0);
      std::vector<double> max_overweight_by_thread(n_threads, 0.0);

#pragma omp parallel num_threads(n_threads)
      {
        const int tid = omp_get_thread_num();
        TRandom3& thread_rng = thread_rngs[tid];
        auto& accepted = accepted_by_thread[tid];
        accepted.reserve(static_cast<size_t>(
            std::max<long long>(16, trials_per_batch / n_threads / 1000)));

#pragma omp for schedule(static)
        for (long long trial = 0; trial < trials_per_batch; ++trial) {
          Candidate c;
          try {
            c = MakeCandidate(thread_rng);
          } catch (const std::exception&) {
            ++rejected_by_thread[tid];
            continue;
          }
          const double w = CandidateWeight(c);
          if (!(w > 0.0) || !std::isfinite(w)) {
            ++rejected_by_thread[tid];
            continue;
          }
          if (w > unweighting_envelope) {
            ++overweight_by_thread[tid];
            max_overweight_by_thread[tid] = std::max(max_overweight_by_thread[tid], w);
            continue;
          }
          if (thread_rng.Uniform(0.0, 1.0) > w / unweighting_envelope) {
            ++rejected_by_thread[tid];
            continue;
          }
          accepted.push_back(std::move(c));
        }
      }

      double max_overweight = 0.0;
      for (int tid = 0; tid < n_threads; ++tid) {
        rejected += rejected_by_thread[tid];
        overweight_candidates += overweight_by_thread[tid];
        max_overweight = std::max(max_overweight, max_overweight_by_thread[tid]);
      }
      if (overweight_candidates > 0) {
        throw std::runtime_error(
            "unweighting envelope was exceeded in parallel mode; max observed weight=" +
            std::to_string(max_overweight) +
            "; increase unweighting_trials/safety or run weighted mode");
      }
      for (const auto& accepted : accepted_by_thread) {
        for (const auto& c : accepted) {
          if (filled >= cfg_.n_events) break;
          fill_candidate(c, 1.0);
          ++filled;
        }
        if (filled >= cfg_.n_events) break;
      }
    }
#else
    throw std::runtime_error("mode1_threads > 1 requires an OpenMP-enabled build");
#endif
  }

  const long long max_trials_per_event = 2000000LL;
  const long long max_trials =
      cfg_.n_events > std::numeric_limits<long long>::max() / max_trials_per_event
          ? std::numeric_limits<long long>::max()
          : std::max(max_trials_per_event, cfg_.n_events * max_trials_per_event);
  for (long long trials = 0; filled < cfg_.n_events && trials < max_trials; ++trials) {
    Candidate c;
    try {
      c = MakeCandidate(rng);
    } catch (const std::exception&) {
      ++rejected;
      continue;
    }

    const double w = CandidateWeight(c);
    if (!(w > 0.0)) {
      ++rejected;
      continue;
    }
    if (cfg_.unweighted_events) {
      const double envelope = unweighting_envelope;
      if (w > envelope) {
        ++overweight_candidates;
        throw std::runtime_error(
            "unweighting envelope was exceeded; increase unweighting_trials/safety or run weighted mode");
      }
      const double accept_prob = w / envelope;
      if (rng.Uniform(0.0, 1.0) > accept_prob) {
        ++rejected;
        continue;
      }
      event_weight = 1.0;
    } else {
      // Weighted event mode: accept all valid candidates, store weight in TTree.
      // No accept-reject — user normalizes with weights in analysis.
      event_weight = w;
    }
    M_rhoprime = c.M;
    y_rhoprime = c.y;
    pt_rhoprime = std::sqrt(c.px * c.px + c.py * c.py);
    m_rho = c.m_rho;
    m_sigma = c.m_sigma;
    pol_xx = c.pol_xx;
    pol_yy = c.pol_yy;
    pol_rexy = c.pol_rexy;
    pol_imxy = c.pol_imxy;
    pi1 = c.p1_lab;
    pi2 = c.p2_lab;
    pi3 = c.p3_lab;
    pi4 = c.p4_lab;
    rho_momentum = c.rho_lab;
    sigma_momentum = c.sigma_lab;

    const TLorentzVector residual = c.rhoprime_lab - pi1 - pi2 - pi3 - pi4;
    if (std::abs(residual.E()) > 1.0e-8 || residual.Vect().Mag() > 1.0e-8) {
      throw std::runtime_error("internal four-momentum conservation failure");
    }

    tree.Fill();
    ++filled;
  }

  if (filled < cfg_.n_events) {
    std::cerr << "Warning: generated " << filled << "/" << cfg_.n_events
              << " events before max trial limit.\n";
  }

  TParameter<long long>("requested_events", cfg_.n_events).Write();
  TParameter<long long>("filled_events", filled).Write();
  TParameter<long long>("rejected_candidates", rejected).Write();
  TParameter<long long>("overweight_candidates", overweight_candidates).Write();
  TParameter<int>("weighted_event_mode", cfg_.unweighted_events ? 0 : 1).Write();
  TParameter<int>("unweighted_event_mode", cfg_.unweighted_events ? 1 : 0).Write();
  TParameter<int>("unweighting_mode", cfg_.unweighting_mode).Write();
  TParameter<double>("unweighting_envelope", unweighting_envelope).Write();
  TParameter<long long>("unweighting_trials", cfg_.unweighting_trials).Write();
  TParameter<double>("unweighting_safety_factor", cfg_.unweighting_safety_factor).Write();
  TParameter<int>("output_compression_level", cfg_.output_compression_level).Write();
  TParameter<int>("mode1_threads_requested", cfg_.mode1_threads).Write();
  TParameter<int>("mode1_threads_used", mode1_threads_used).Write();
  TParameter<int>("bose_symmetrize", cfg_.bose_symmetrize ? 1 : 0).Write();
  TParameter<double>("cross_section_mb", cross_section_mb_).Write();
  TParameter<long long>("decay_norm_trials_per_mass",
                        cfg_.decay_norm_trials_per_mass).Write();
  TParameter<double>("m_pi_GeV", kMPi).Write();
  TParameter<double>("m_rho_GeV", kMRho).Write();
  TParameter<double>("gamma_rho_GeV", kGRho).Write();
  TParameter<double>("m_sigma_GeV", kMSigma).Write();
  TParameter<double>("gamma_sigma_GeV", kGSigma).Write();
  tree.Write();
  out->Close();

  return BoseGenerationSummary{cfg_.n_events, filled, rejected, overweight_candidates,
                               unweighting_envelope, cfg_.output_root};
}

}  // namespace rhoprime

#include "RhoPrime/UpcRhoPrimeModelOpt.h"

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TString.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace rhoprime {
namespace {

struct SystemConfig {
  std::string name;
  int Z = 0;
  int A = 0;
  double R_WS = 0.0;
  double a_WS = 0.0;
  double rho_0 = 0.0;
  double R_A_hard = 0.0;
};

struct GridPoint {
  double x = 0.0;
  double y = 0.0;
  double r = 0.0;
  int idx_rt = 0;
  double TA = 0.0;
};

constexpr int kLutBins = 3000;
constexpr double kDrLutFm = 0.1;
constexpr int kNLoc = 80;
constexpr double kLLocFm = 40.0;
constexpr double kDxLocFm = kLLocFm / kNLoc;

constexpr double kAlpha = 1.0 / 137.036;
constexpr double kHbarC = 0.197327;
constexpr double kHbarC2Fm2 = 0.0389379;
constexpr double kMP = 0.93827;
constexpr double kBV = 9.4;
constexpr double kFV2_4PI = 8.4;
constexpr double kBR4PI = 0.2;
constexpr double kXP = 0.00026;
constexpr double kEps = 0.28;
constexpr double kYP = 0.02078;
constexpr double kEta = 1.21;
constexpr double kMPi = 0.13957;
constexpr double kMThreshold = 4.0 * kMPi;
constexpr double kRapidityMirrorEps = 1.0e-12;

constexpr double kFitARes = 1.11569;
constexpr double kFitMRes = 1.57773;
constexpr double kFitGamma0 = 0.611252;
constexpr double kFitANr = 0.151507;
constexpr double kFitBNr = 0.0836185;
constexpr double kFitPhi = 5.17036e-5;

struct Workspace {
  SystemConfig sys;
  UpcGridConfig cfg;
  TH1D* h_prob_surv = nullptr;
  TH1D* h_prob_tag = nullptr;
  std::vector<double> ta_table;
  std::vector<GridPoint> grid;
  double dp_global = 0.0;
  int n_bins_total = 0;
  std::vector<std::vector<std::complex<double>>> phase_x;
  std::vector<std::vector<std::complex<double>>> phase_y;
  double mass_norm = 1.0;
  bool probability_file_loaded = false;
  bool probability_fallback_used = false;
  std::string probability_file_path;
};

struct ThreadLocalData {
  ThreadLocalData(int n_bins, int n_y_bins)
      : n_y_bins(n_y_bins),
        n_bins(n_bins),
        hist_sigma(n_y_bins, std::vector<double>(n_bins * n_bins, 0.0)),
        hist_Ax2(n_y_bins, std::vector<double>(n_bins * n_bins, 0.0)),
        hist_Ay2(n_y_bins, std::vector<double>(n_bins * n_bins, 0.0)),
        hist_ReAxAy(n_y_bins, std::vector<double>(n_bins * n_bins, 0.0)),
        hist_ImAxAy(n_y_bins, std::vector<double>(n_bins * n_bins, 0.0)),
        sum_y_A1x(n_bins, std::vector<std::complex<double>>(kNLoc)),
        sum_y_A1y(n_bins, std::vector<std::complex<double>>(kNLoc)),
        sum_y_A2x(n_bins, std::vector<std::complex<double>>(kNLoc)),
        sum_y_A2y(n_bins, std::vector<std::complex<double>>(kNLoc)) {}

  int n_y_bins;
  int n_bins;
  std::vector<std::vector<double>> hist_sigma;
  std::vector<std::vector<double>> hist_Ax2;
  std::vector<std::vector<double>> hist_Ay2;
  std::vector<std::vector<double>> hist_ReAxAy;
  std::vector<std::vector<double>> hist_ImAxAy;
  double sum_spatial = 0.0;
  double sum_momentum = 0.0;
  std::vector<std::vector<std::complex<double>>> sum_y_A1x;
  std::vector<std::vector<std::complex<double>>> sum_y_A1y;
  std::vector<std::vector<std::complex<double>>> sum_y_A2x;
  std::vector<std::vector<std::complex<double>>> sum_y_A2y;
};

struct OmegaLUT {
  std::vector<double> omega_values;
  std::vector<double> sig_vp;
  std::vector<std::vector<double>> photon_field;
  int n_bins = 0;
};

int MaxThreads() {
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
}

int ThreadNum() {
#ifdef _OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
}

void Require(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error(message);
}

void SetupSystem(Workspace& ws) {
  if (ws.cfg.system == "AuAu" || ws.cfg.system == "Au" || ws.cfg.system == "Au197") {
    ws.sys = SystemConfig{"AuAu", 79, 197, 6.38, 0.535, 0.169, 0.0};
  } else if (ws.cfg.system == "PbPb" || ws.cfg.system == "Pb" || ws.cfg.system == "Pb208") {
    ws.sys = SystemConfig{"PbPb", 82, 208, 6.62, 0.546, 0.161, 0.0};
  } else {
    throw std::invalid_argument("system must be AuAu or PbPb");
  }
  ws.sys.R_A_hard = 1.2 * std::pow(static_cast<double>(ws.sys.A), 1.0 / 3.0);
}

double ThicknessFunctionWS(const Workspace& ws, double r) {
  const double z_max = 15.0;
  const int steps = 100;
  const double dz = z_max / steps;
  double sum = 0.0;
  for (int i = 0; i < steps; ++i) {
    const double z = (i + 0.5) * dz;
    sum += ws.sys.rho_0 /
           (1.0 + std::exp((std::sqrt(r * r + z * z) - ws.sys.R_WS) / ws.sys.a_WS));
  }
  return 2.0 * sum * dz;
}

double CoherentLengthFactor(const Workspace& ws, double r, double q_z) {
  double sum_cos = 0.0;
  double sum_norm = 0.0;
  const double z_max = 15.0;
  const int steps = 100;
  const double dz = z_max / steps;
  for (int i = 0; i < steps; ++i) {
    const double z = (i + 0.5) * dz;
    const double rho = ws.sys.rho_0 /
                       (1.0 + std::exp((std::sqrt(r * r + z * z) - ws.sys.R_WS) /
                                        ws.sys.a_WS));
    sum_cos += rho * std::cos(q_z * z / kHbarC) * dz;
    sum_norm += rho * dz;
  }
  return sum_norm > 0.0 ? sum_cos / sum_norm : 0.0;
}

void PrecomputeGlobals(Workspace& ws) {
  ws.ta_table.assign(kLutBins, 0.0);
  for (int i = 0; i < kLutBins; ++i) ws.ta_table[i] = ThicknessFunctionWS(ws, i * kDrLutFm);

  ws.grid.clear();
  ws.grid.reserve(kNLoc * kNLoc);
  for (int ix = 0; ix < kNLoc; ++ix) {
    const double x = -kLLocFm / 2.0 + (ix + 0.5) * kDxLocFm;
    for (int iy = 0; iy < kNLoc; ++iy) {
      const double y = -kLLocFm / 2.0 + (iy + 0.5) * kDxLocFm;
      const double r = std::sqrt(x * x + y * y);
      int idx = static_cast<int>(r / kDrLutFm);
      idx = std::clamp(idx, 0, kLutBins - 1);
      ws.grid.push_back(GridPoint{x, y, r, idx, ws.ta_table[idx]});
    }
  }

  ws.dp_global = (2.0 * TMath::Pi() / ws.cfg.global_box_fm) * kHbarC;
  ws.n_bins_total = 2 * ws.cfg.n_p_bins_each_side + 1;
  ws.phase_x.assign(ws.n_bins_total, std::vector<std::complex<double>>(kNLoc));
  ws.phase_y.assign(ws.n_bins_total, std::vector<std::complex<double>>(kNLoc));
  for (int ip = -ws.cfg.n_p_bins_each_side; ip <= ws.cfg.n_p_bins_each_side; ++ip) {
    const double p = ip * ws.dp_global;
    const int ip_idx = ip + ws.cfg.n_p_bins_each_side;
    for (int k = 0; k < kNLoc; ++k) {
      const double coord = -kLLocFm / 2.0 + (k + 0.5) * kDxLocFm;
      const std::complex<double> val = std::exp(std::complex<double>(0.0, -p * coord / kHbarC));
      ws.phase_x[ip_idx][k] = val;
      ws.phase_y[ip_idx][k] = val;
    }
  }
}

void LoadProbabilityHistograms(Workspace& ws) {
  const auto exists = [](const std::string& path) {
    if (std::FILE* fp = std::fopen(path.c_str(), "r")) {
      std::fclose(fp);
      return true;
    }
    return false;
  };
  const auto canonical_name = [&]() {
    return "UPC_Probabilities_" + ws.sys.name + "_" +
           std::to_string(static_cast<int>(std::round(ws.cfg.sqrt_s_NN_GeV))) +
           "GeV_Glauber_Final.root";
  };
  std::vector<std::string> names;
  names.push_back(canonical_name());
  if (ws.sys.name == "PbPb" && std::abs(ws.cfg.sqrt_s_NN_GeV - 5360.0) < 1.0) {
    names.push_back("UPC_Probabilities_PbPb_5.36TeV_Glauber_Final.root");
  }
  if (ws.sys.name == "PbPb" && std::abs(ws.cfg.sqrt_s_NN_GeV - 5020.0) < 50.0) {
    names.push_back("UPC_Probabilities_PbPb_5.36TeV_Glauber_Final.root");
  }
  if (ws.sys.name == "AuAu" && std::abs(ws.cfg.sqrt_s_NN_GeV - 200.0) < 1.0) {
    names.push_back("UPC_Probabilities_AuAu_200GeV_Glauber_Final.root");
  }

  std::vector<std::string> prefixes = {"", "rhoprime/", "doublerho/", "../rhoprime/",
                                       "../doublerho/", "../../rhoprime/", "../../doublerho/"};
  std::string selected;
  for (const auto& prefix : prefixes) {
    for (const auto& name : names) {
      const std::string path = prefix + name;
      if (exists(path)) {
        selected = path;
        break;
      }
    }
    if (!selected.empty()) break;
  }

  std::unique_ptr<TFile> f(selected.empty() ? nullptr : TFile::Open(selected.c_str(), "READ"));
  if (!f || f->IsZombie()) {
    ws.probability_fallback_used = true;
    if (ws.cfg.require_probability_file) {
      throw std::runtime_error("cannot load required UPC probability ROOT for " + ws.sys.name +
                               " @ " + std::to_string(ws.cfg.sqrt_s_NN_GeV));
    }
    std::cerr << "Warning: cannot load UPC probability ROOT for " << ws.sys.name << " @ "
              << ws.cfg.sqrt_s_NN_GeV
              << "; using survival=1 and tagged probability=0 fallback.\n";
    return;
  }
  ws.probability_file_loaded = true;
  ws.probability_file_path = selected;
  std::cout << "Loaded UPC probability ROOT: " << selected << "\n";
  if (auto* h = dynamic_cast<TH1D*>(f->Get("hProb_Surv"))) {
    ws.h_prob_surv = dynamic_cast<TH1D*>(h->Clone("hProb_Surv_rhoprime"));
    ws.h_prob_surv->SetDirectory(nullptr);
  }
  if (ws.cfg.trigger != "NoTag") {
    std::string tag = ws.cfg.trigger;
    if (tag.rfind("EMD_", 0) == 0) tag = tag.substr(4);
    const std::vector<std::string> tag_candidates = {"hProb_EMD_" + tag, "hProb_" + tag};
    for (const auto& tag_name : tag_candidates) {
      if (auto* h = dynamic_cast<TH1D*>(f->Get(tag_name.c_str()))) {
        ws.h_prob_tag = dynamic_cast<TH1D*>(h->Clone("hProb_Tag_rhoprime"));
        ws.h_prob_tag->SetDirectory(nullptr);
        std::cout << "Loaded UPC tag probability histogram: " << tag_name << "\n";
        break;
      }
    }
    if (!ws.h_prob_tag) {
      throw std::runtime_error("missing pure EMD tag histogram for trigger " + ws.cfg.trigger +
                               " in " + selected +
                               "; expected hProb_EMD_1n, hProb_EMD_Xn, "
                               "hProb_EMD_1n1n, or hProb_EMD_XnXn");
    }
  }
}

double GammaM(double M, double M0, double Gamma0) {
  if (M <= kMThreshold) return 0.0;
  const double q_M = std::sqrt(M * M - kMThreshold * kMThreshold);
  const double q_M0 = std::sqrt(M0 * M0 - kMThreshold * kMThreshold);
  return Gamma0 * (M0 / M) * std::pow(q_M / q_M0, 3);
}

double MassDistUnnormalized(double M) {
  if (M <= kMThreshold) return 0.0;
  const double s = M * M;
  const double gamma = GammaM(M, kFitMRes, kFitGamma0);
  const double denom = std::pow(kFitMRes * kFitMRes - s, 2) + std::pow(kFitMRes * gamma, 2);
  const double num_res = kFitARes * std::sqrt(kFitMRes * gamma);
  const double re_res = num_res * (kFitMRes * kFitMRes - s) / denom;
  const double im_res = num_res * (kFitMRes * gamma) / denom;
  const double phase_space = std::sqrt(1.0 - std::pow(kMThreshold / M, 2));
  const double amp_nr = kFitANr * phase_space * std::exp(-kFitBNr * (M - kMThreshold));
  const double re_nr = amp_nr * std::cos(kFitPhi);
  const double im_nr = amp_nr * std::sin(kFitPhi);
  return std::pow(re_res + re_nr, 2) + std::pow(im_res + im_nr, 2);
}

void InitMassNormalization(Workspace& ws) {
  double sum = 0.0;
  const int steps = 5000;
  const double dM = (3.5 - kMThreshold) / steps;
  for (int i = 0; i < steps; ++i) {
    const double M = kMThreshold + (i + 0.5) * dM;
    sum += MassDistUnnormalized(M) * dM;
  }
  Require(sum > 0.0, "rho-prime mass normalization failed");
  ws.mass_norm = sum;
}

double NormalizedMassDist(const Workspace& ws, double M) {
  return MassDistUnnormalized(M) / ws.mass_norm;
}

double SigmaGammaP(double W) {
  if (W < 4.0 * kMPi + kMP) return 0.0;
  const double sig_mb = (kXP * std::pow(W, kEps) + kYP * std::pow(W, -kEta)) / kBR4PI;
  return sig_mb * 0.1;
}

double SigmaTotVp(double W) {
  const double dsig_vp_dt0 = kBV * SigmaGammaP(W);
  return std::sqrt(16.0 * TMath::Pi() * kHbarC2Fm2 * dsig_vp_dt0 * (kFV2_4PI / kAlpha));
}

double CalcPhotonField(double omega, double r, double gamma_L, double R_A) {
  if (r < 1.0e-3) r = 1.0e-3;
  const double xi = omega * r / (gamma_L * kHbarC);
  if (xi >= 20.0) return 0.0;
  if (r >= R_A) return xi * TMath::BesselK1(xi) / r;
  const double xi_R = omega * R_A / (gamma_L * kHbarC);
  return (xi_R * TMath::BesselK1(xi_R) / R_A) * (r / R_A);
}

double CalcAbsorptionProfile(double TA_r, double sig_vp) {
  return 2.0 * (1.0 - std::exp(-sig_vp * TA_r / 2.0));
}

double InterpProbability(TH1D* h, double b, double fallback) {
  if (!h) return fallback;
  if (b <= h->GetXaxis()->GetXmax()) return std::max(0.0, h->Interpolate(b));
  return fallback;
}

OmegaLUT PrecomputeOmegaLUT(const Workspace& ws, int n_omega) {
  Require(n_omega >= 2, "n_omega must be at least 2");
  const double y_abs = std::max(std::abs(ws.cfg.rapidity_min), std::abs(ws.cfg.rapidity_max));
  const double omega_min = 0.5 * ws.cfg.mass_min_GeV * std::exp(-y_abs);
  const double omega_max = 0.5 * ws.cfg.mass_max_GeV * std::exp(y_abs);
  Require(omega_min > 0.0 && omega_max > omega_min, "invalid omega lookup range");

  OmegaLUT lut;
  lut.omega_values.resize(n_omega);
  lut.sig_vp.assign(n_omega, 0.0);
  const double photon_r_max =
      std::max(ws.cfg.b_tail_max_notag_fm, ws.cfg.b_tail_max_tagged_fm) + 16.0;
  const int photon_bins = static_cast<int>(std::ceil(photon_r_max / kDrLutFm)) + 1;
  lut.n_bins = photon_bins;
  lut.photon_field.assign(n_omega, std::vector<double>(photon_bins, 0.0));

  const double gamma_L = ws.cfg.sqrt_s_NN_GeV / (2.0 * kMP);
  const double log_min = std::log(omega_min);
  const double dlog = (std::log(omega_max) - log_min) / (n_omega - 1);

  for (int iw = 0; iw < n_omega; ++iw) {
    const double omega = std::exp(log_min + iw * dlog);
    lut.omega_values[iw] = omega;
    const double sig_vp =
        SigmaTotVp(std::sqrt(kMP * kMP + 4.0 * omega * gamma_L * kMP));
    lut.sig_vp[iw] = sig_vp;

    for (int i = 0; i < photon_bins; ++i) {
      lut.photon_field[iw][i] =
          CalcPhotonField(omega, i * kDrLutFm, gamma_L, ws.sys.R_WS);
    }
  }
  return lut;
}

void OmegaBracket(const OmegaLUT& lut, double omega, int& lo, int& hi, double& t) {
  Require(!lut.omega_values.empty(), "empty omega lookup table");
  if (omega <= lut.omega_values.front()) {
    lo = hi = 0;
    t = 0.0;
    return;
  }
  if (omega >= lut.omega_values.back()) {
    lo = hi = static_cast<int>(lut.omega_values.size()) - 1;
    t = 0.0;
    return;
  }
  const double log_omega = std::log(omega);
  auto upper = std::upper_bound(lut.omega_values.begin(), lut.omega_values.end(), omega);
  hi = static_cast<int>(upper - lut.omega_values.begin());
  lo = hi - 1;
  t = (log_omega - std::log(lut.omega_values[lo])) /
      (std::log(lut.omega_values[hi]) - std::log(lut.omega_values[lo]));
}

double InterpolateOmegaScalar(const OmegaLUT& lut, double omega, const std::vector<double>& src) {
  int lo = 0, hi = 0;
  double t = 0.0;
  OmegaBracket(lut, omega, lo, hi, t);
  return src[lo] * (1.0 - t) + src[hi] * t;
}

void InterpolateOmegaVector(const OmegaLUT& lut, double omega,
                            const std::vector<std::vector<double>>& src,
                            std::vector<double>& out) {
  int lo = 0, hi = 0;
  double t = 0.0;
  OmegaBracket(lut, omega, lo, hi, t);
  out.assign(src[lo].size(), 0.0);
  for (size_t i = 0; i < out.size(); ++i) out[i] = src[lo][i] * (1.0 - t) + src[hi][i] * t;
}

double LookupRadial(const std::vector<double>& values, double r) {
  if (r <= 0.0) return values.front();
  const double pos = r / kDrLutFm;
  const int lo = static_cast<int>(pos);
  if (lo >= static_cast<int>(values.size()) - 1) return values.back();
  const double t = pos - lo;
  return values[lo] * (1.0 - t) + values[lo + 1] * t;
}

double D2SigmaDMdyOpt(const Workspace& ws, const OmegaLUT& lut, double M, double y, double dM,
                      double dy, int iy_bin, ThreadLocalData& tl) {
  const double gamma_L = ws.cfg.sqrt_s_NN_GeV / (2.0 * kMP);
  const double P_M = NormalizedMassDist(ws, M);
  if (P_M < 1.0e-6) return 0.0;

  const double omega1 = 0.5 * M * std::exp(y);
  const double omega2 = 0.5 * M * std::exp(-y);
  const double sig_vp1 = InterpolateOmegaScalar(lut, omega1, lut.sig_vp);
  const double sig_vp2 = InterpolateOmegaScalar(lut, omega2, lut.sig_vp);
  if (sig_vp1 <= 0.0 && sig_vp2 <= 0.0) return 0.0;

  const double qz1 = (M * M) / (4.0 * gamma_L * omega1);
  const double qz2 = (M * M) / (4.0 * gamma_L * omega2);
  std::vector<double> clf1(kLutBins), clf2(kLutBins);
  for (int i = 0; i < kLutBins; ++i) {
    clf1[i] = CoherentLengthFactor(ws, i * kDrLutFm, qz1);
    clf2[i] = CoherentLengthFactor(ws, i * kDrLutFm, qz2);
  }

  std::vector<double> photon_field1, photon_field2;
  InterpolateOmegaVector(lut, omega1, lut.photon_field, photon_field1);
  InterpolateOmegaVector(lut, omega2, lut.photon_field, photon_field2);

  const double fft_norm = (ws.sys.Z * ws.sys.Z * kAlpha * kAlpha * kBR4PI * 10.0) /
                          (16.0 * std::pow(TMath::Pi(), 4) * kHbarC2Fm2 * kFV2_4PI);
  double d2sigma = 0.0;

  auto eval_b = [&](double b, double db_step) {
    const double hard_survival = (b < 2.0 * ws.sys.R_WS) ? 0.0 : 1.0;
    const double p_surv = InterpProbability(ws.h_prob_surv, b, hard_survival);
    if (p_surv <= 0.0) return 0.0;
    double p_tag = 1.0;
    if (ws.cfg.trigger != "NoTag") p_tag = InterpProbability(ws.h_prob_tag, b, 0.0);
    const double total_prob = p_surv * p_tag;
    if (total_prob <= 0.0) return 0.0;

    double sum_spatial_local = 0.0;
    for (int ix = 0; ix < kNLoc; ++ix) {
      std::vector<double> A1x(kNLoc, 0.0), A1y(kNLoc, 0.0), A2x(kNLoc, 0.0), A2y(kNLoc, 0.0);
      for (int iy = 0; iy < kNLoc; ++iy) {
        const auto& pt = ws.grid[ix * kNLoc + iy];
        if (pt.r < 0.2 || pt.r > 16.0) continue;
        const double r1 = std::sqrt((pt.x - b) * (pt.x - b) + pt.y * pt.y);
        const double r2 = std::sqrt((pt.x + b) * (pt.x + b) + pt.y * pt.y);
        if (r1 > 0.2) {
          const double mag = LookupRadial(photon_field1, r1) *
                             CalcAbsorptionProfile(pt.TA, sig_vp1) * clf1[pt.idx_rt];
          A1x[iy] = mag * (pt.x - b) / r1;
          A1y[iy] = mag * pt.y / r1;
          sum_spatial_local += A1x[iy] * A1x[iy] + A1y[iy] * A1y[iy];
        }
        if (r2 > 0.2) {
          const double mag = LookupRadial(photon_field2, r2) *
                             CalcAbsorptionProfile(pt.TA, sig_vp2) * clf2[pt.idx_rt];
          A2x[iy] = mag * (pt.x + b) / r2;
          A2y[iy] = mag * pt.y / r2;
          sum_spatial_local += A2x[iy] * A2x[iy] + A2y[iy] * A2y[iy];
        }
      }
      for (int ipy = 0; ipy < ws.n_bins_total; ++ipy) {
        std::complex<double> sy1x(0, 0), sy1y(0, 0), sy2x(0, 0), sy2y(0, 0);
        for (int iy = 0; iy < kNLoc; ++iy) {
          const auto pf = ws.phase_y[ipy][iy];
          sy1x += A1x[iy] * pf;
          sy1y += A1y[iy] * pf;
          sy2x += A2x[iy] * pf;
          sy2y += A2y[iy] * pf;
        }
        tl.sum_y_A1x[ipy][ix] = sy1x;
        tl.sum_y_A1y[ipy][ix] = sy1y;
        tl.sum_y_A2x[ipy][ix] = sy2x;
        tl.sum_y_A2y[ipy][ix] = sy2y;
      }
    }

    const double phase_space_vol = (2.0 * TMath::Pi() * b * db_step) * dM * dy;
    const double weight_mapping = total_prob * P_M * phase_space_vol;
    const double spatial_mb = sum_spatial_local * (kDxLocFm * kDxLocFm) * fft_norm *
                              (4.0 * TMath::Pi() * TMath::Pi() * kHbarC * kHbarC);
    tl.sum_spatial += spatial_mb * weight_mapping;

    const double d2sig_factor = std::pow(kDxLocFm, 4) * fft_norm * weight_mapping;
    double sum_momentum_raw = 0.0;
    const int zero_idx = ws.cfg.n_p_bins_each_side;
    auto add_momentum_bin = [&](int ix, int iy, double Ax2, double Ay2, double Re, double Im,
                                double sign_re, double sign_im) {
      const int idx = ix * ws.n_bins_total + iy;
      tl.hist_sigma[iy_bin][idx] += (Ax2 + Ay2) * d2sig_factor;
      tl.hist_Ax2[iy_bin][idx] += Ax2 * d2sig_factor;
      tl.hist_Ay2[iy_bin][idx] += Ay2 * d2sig_factor;
      tl.hist_ReAxAy[iy_bin][idx] += sign_re * Re * d2sig_factor;
      tl.hist_ImAxAy[iy_bin][idx] += sign_im * Im * d2sig_factor;
      sum_momentum_raw += Ax2 + Ay2;
    };

    for (int ipx = zero_idx; ipx < ws.n_bins_total; ++ipx) {
      const double px = (ipx - ws.cfg.n_p_bins_each_side) * ws.dp_global;
      const auto shift1 = std::exp(std::complex<double>(0.0, px * b / (2.0 * kHbarC)));
      const auto shift2 = std::exp(std::complex<double>(0.0, -px * b / (2.0 * kHbarC)));
      for (int ipy = zero_idx; ipy < ws.n_bins_total; ++ipy) {
        std::complex<double> A1x(0, 0), A1y(0, 0), A2x(0, 0), A2y(0, 0);
        for (int ix = 0; ix < kNLoc; ++ix) {
          const auto pf = ws.phase_x[ipx][ix];
          A1x += tl.sum_y_A1x[ipy][ix] * pf;
          A1y += tl.sum_y_A1y[ipy][ix] * pf;
          A2x += tl.sum_y_A2x[ipy][ix] * pf;
          A2y += tl.sum_y_A2y[ipy][ix] * pf;
        }
        const auto Ax = A1x * shift1 + A2x * shift2;
        const auto Ay = A1y * shift1 + A2y * shift2;
        const double Ax2 = std::norm(Ax);
        const double Ay2 = std::norm(Ay);
        const auto AxAy = Ax * std::conj(Ay);
        const double re = AxAy.real();
        const double im = AxAy.imag();

        add_momentum_bin(ipx, ipy, Ax2, Ay2, re, im, 1.0, 1.0);
        if (ipx != zero_idx) {
          add_momentum_bin(ws.n_bins_total - ipx - 1, ipy, Ax2, Ay2, re, im, -1.0, 1.0);
        }
        if (ipy != zero_idx) {
          add_momentum_bin(ipx, ws.n_bins_total - ipy - 1, Ax2, Ay2, re, im, -1.0, -1.0);
        }
        if (ipx != zero_idx && ipy != zero_idx) {
          add_momentum_bin(ws.n_bins_total - ipx - 1, ws.n_bins_total - ipy - 1, Ax2, Ay2, re,
                           im, 1.0, -1.0);
        }
      }
    }

    const double momentum_mb =
        sum_momentum_raw * std::pow(kDxLocFm, 4) * fft_norm * ws.dp_global * ws.dp_global;
    tl.sum_momentum += momentum_mb * weight_mapping;
    return momentum_mb * total_prob;
  };

  const double db_core = (ws.cfg.b_core_max_fm - ws.cfg.b_core_min_fm) / ws.cfg.b_core_steps;
  for (int j = 0; j < ws.cfg.b_core_steps; ++j) {
    const double b = ws.cfg.b_core_min_fm + (j + 0.5) * db_core;
    d2sigma += eval_b(b, db_core) * (2.0 * TMath::Pi() * b * db_core);
  }
  const double b_tail_max =
      (ws.cfg.trigger == "NoTag") ? ws.cfg.b_tail_max_notag_fm : ws.cfg.b_tail_max_tagged_fm;
  const int tail_steps =
      std::max(0, static_cast<int>((b_tail_max - ws.cfg.b_tail_min_fm) / ws.cfg.b_tail_step_fm));
  if (tail_steps > 0) {
    const double db_tail = (b_tail_max - ws.cfg.b_tail_min_fm) / tail_steps;
    for (int j = 0; j < tail_steps; ++j) {
      const double b = ws.cfg.b_tail_min_fm + (j + 0.5) * db_tail;
      d2sigma += eval_b(b, db_tail) * (2.0 * TMath::Pi() * b * db_tail);
    }
  }
  return d2sigma * P_M;
}

double D2SigmaDMdy(Workspace& ws, double M, double y, double dM, double dy, int iy_bin,
                   ThreadLocalData& tl) {
  const double gamma_L = ws.cfg.sqrt_s_NN_GeV / (2.0 * kMP);
  const double P_M = NormalizedMassDist(ws, M);
  if (P_M < 1.0e-6) return 0.0;

  const double omega1 = 0.5 * M * std::exp(y);
  const double omega2 = 0.5 * M * std::exp(-y);
  const double sig_vp1 = SigmaTotVp(std::sqrt(kMP * kMP + 4.0 * omega1 * gamma_L * kMP));
  const double sig_vp2 = SigmaTotVp(std::sqrt(kMP * kMP + 4.0 * omega2 * gamma_L * kMP));
  if (sig_vp1 <= 0.0 && sig_vp2 <= 0.0) return 0.0;

  const double qz1 = (M * M) / (4.0 * gamma_L * omega1);
  const double qz2 = (M * M) / (4.0 * gamma_L * omega2);
  std::vector<double> clf1(kLutBins), clf2(kLutBins);
  for (int i = 0; i < kLutBins; ++i) {
    clf1[i] = CoherentLengthFactor(ws, i * kDrLutFm, qz1);
    clf2[i] = CoherentLengthFactor(ws, i * kDrLutFm, qz2);
  }

  const double fft_norm = (ws.sys.Z * ws.sys.Z * kAlpha * kAlpha * kBR4PI * 10.0) /
                          (16.0 * std::pow(TMath::Pi(), 4) * kHbarC2Fm2 * kFV2_4PI);
  double d2sigma = 0.0;

  auto eval_b = [&](double b, double db_step) {
    const double hard_survival = (b < 2.0 * ws.sys.R_WS) ? 0.0 : 1.0;
    const double p_surv = InterpProbability(ws.h_prob_surv, b, hard_survival);
    if (p_surv <= 0.0) return 0.0;
    double p_tag = 1.0;
    if (ws.cfg.trigger != "NoTag") p_tag = InterpProbability(ws.h_prob_tag, b, 0.0);
    const double total_prob = p_surv * p_tag;
    if (total_prob <= 0.0) return 0.0;

    double sum_spatial_local = 0.0;
    for (int ix = 0; ix < kNLoc; ++ix) {
      std::vector<double> A1x(kNLoc, 0.0), A1y(kNLoc, 0.0), A2x(kNLoc, 0.0), A2y(kNLoc, 0.0);
      for (int iy = 0; iy < kNLoc; ++iy) {
        const auto& pt = ws.grid[ix * kNLoc + iy];
        if (pt.r < 0.2 || pt.r > 16.0) continue;
        const double r1 = std::sqrt((pt.x - b) * (pt.x - b) + pt.y * pt.y);
        if (r1 > 0.2) {
          const double mag = CalcPhotonField(omega1, r1, gamma_L, ws.sys.R_WS) *
                             CalcAbsorptionProfile(pt.TA, sig_vp1) * clf1[pt.idx_rt];
          A1x[iy] = mag * (pt.x - b) / r1;
          A1y[iy] = mag * pt.y / r1;
          sum_spatial_local += A1x[iy] * A1x[iy] + A1y[iy] * A1y[iy];
        }
        const double r2 = std::sqrt((pt.x + b) * (pt.x + b) + pt.y * pt.y);
        if (r2 > 0.2) {
          const double mag = CalcPhotonField(omega2, r2, gamma_L, ws.sys.R_WS) *
                             CalcAbsorptionProfile(pt.TA, sig_vp2) * clf2[pt.idx_rt];
          A2x[iy] = mag * (pt.x + b) / r2;
          A2y[iy] = mag * pt.y / r2;
          sum_spatial_local += A2x[iy] * A2x[iy] + A2y[iy] * A2y[iy];
        }
      }
      for (int ipy = 0; ipy < ws.n_bins_total; ++ipy) {
        std::complex<double> sy1x(0, 0), sy1y(0, 0), sy2x(0, 0), sy2y(0, 0);
        for (int iy = 0; iy < kNLoc; ++iy) {
          const auto pf = ws.phase_y[ipy][iy];
          sy1x += A1x[iy] * pf;
          sy1y += A1y[iy] * pf;
          sy2x += A2x[iy] * pf;
          sy2y += A2y[iy] * pf;
        }
        tl.sum_y_A1x[ipy][ix] = sy1x;
        tl.sum_y_A1y[ipy][ix] = sy1y;
        tl.sum_y_A2x[ipy][ix] = sy2x;
        tl.sum_y_A2y[ipy][ix] = sy2y;
      }
    }

    const double phase_space_vol = (2.0 * TMath::Pi() * b * db_step) * dM * dy;
    const double weight_mapping = total_prob * P_M * phase_space_vol;
    const double spatial_mb = sum_spatial_local * (kDxLocFm * kDxLocFm) * fft_norm *
                              (4.0 * TMath::Pi() * TMath::Pi() * kHbarC * kHbarC);
    tl.sum_spatial += spatial_mb * weight_mapping;

    const double d2sig_factor = std::pow(kDxLocFm, 4) * fft_norm * weight_mapping;
    double sum_momentum_raw = 0.0;
    const int zero_idx = ws.cfg.n_p_bins_each_side;
    auto add_momentum_bin = [&](int ix, int iy, double Ax2, double Ay2, double Re, double Im,
                                double sign_re, double sign_im) {
      const int idx = ix * ws.n_bins_total + iy;
      tl.hist_sigma[iy_bin][idx] += (Ax2 + Ay2) * d2sig_factor;
      tl.hist_Ax2[iy_bin][idx] += Ax2 * d2sig_factor;
      tl.hist_Ay2[iy_bin][idx] += Ay2 * d2sig_factor;
      tl.hist_ReAxAy[iy_bin][idx] += sign_re * Re * d2sig_factor;
      tl.hist_ImAxAy[iy_bin][idx] += sign_im * Im * d2sig_factor;
      sum_momentum_raw += Ax2 + Ay2;
    };

    for (int ipx = zero_idx; ipx < ws.n_bins_total; ++ipx) {
      const double px = (ipx - ws.cfg.n_p_bins_each_side) * ws.dp_global;
      const auto shift1 = std::exp(std::complex<double>(0.0, px * b / (2.0 * kHbarC)));
      const auto shift2 = std::exp(std::complex<double>(0.0, -px * b / (2.0 * kHbarC)));
      for (int ipy = zero_idx; ipy < ws.n_bins_total; ++ipy) {
        std::complex<double> A1x(0, 0), A1y(0, 0), A2x(0, 0), A2y(0, 0);
        for (int ix = 0; ix < kNLoc; ++ix) {
          const auto pf = ws.phase_x[ipx][ix];
          A1x += tl.sum_y_A1x[ipy][ix] * pf;
          A1y += tl.sum_y_A1y[ipy][ix] * pf;
          A2x += tl.sum_y_A2x[ipy][ix] * pf;
          A2y += tl.sum_y_A2y[ipy][ix] * pf;
        }
        const auto Ax = A1x * shift1 + A2x * shift2;
        const auto Ay = A1y * shift1 + A2y * shift2;
        const double Ax2 = std::norm(Ax);
        const double Ay2 = std::norm(Ay);
        const auto AxAy = Ax * std::conj(Ay);
        const double re = AxAy.real();
        const double im = AxAy.imag();

        add_momentum_bin(ipx, ipy, Ax2, Ay2, re, im, 1.0, 1.0);
        if (ipx != zero_idx) {
          add_momentum_bin(ws.n_bins_total - ipx - 1, ipy, Ax2, Ay2, re, im, -1.0, 1.0);
        }
        if (ipy != zero_idx) {
          add_momentum_bin(ipx, ws.n_bins_total - ipy - 1, Ax2, Ay2, re, im, -1.0, -1.0);
        }
        if (ipx != zero_idx && ipy != zero_idx) {
          add_momentum_bin(ws.n_bins_total - ipx - 1, ws.n_bins_total - ipy - 1, Ax2, Ay2, re,
                           im, 1.0, -1.0);
        }
      }
    }

    const double momentum_mb =
        sum_momentum_raw * std::pow(kDxLocFm, 4) * fft_norm * ws.dp_global * ws.dp_global;
    tl.sum_momentum += momentum_mb * weight_mapping;
    return momentum_mb * total_prob;
  };

  const double db_core = (ws.cfg.b_core_max_fm - ws.cfg.b_core_min_fm) / ws.cfg.b_core_steps;
  for (int j = 0; j < ws.cfg.b_core_steps; ++j) {
    const double b = ws.cfg.b_core_min_fm + (j + 0.5) * db_core;
    d2sigma += eval_b(b, db_core) * (2.0 * TMath::Pi() * b * db_core);
  }
  const double b_tail_max =
      (ws.cfg.trigger == "NoTag") ? ws.cfg.b_tail_max_notag_fm : ws.cfg.b_tail_max_tagged_fm;
  const int tail_steps =
      std::max(0, static_cast<int>((b_tail_max - ws.cfg.b_tail_min_fm) / ws.cfg.b_tail_step_fm));
  if (tail_steps > 0) {
    const double db_tail = (b_tail_max - ws.cfg.b_tail_min_fm) / tail_steps;
    for (int j = 0; j < tail_steps; ++j) {
      const double b = ws.cfg.b_tail_min_fm + (j + 0.5) * db_tail;
      d2sigma += eval_b(b, db_tail) * (2.0 * TMath::Pi() * b * db_tail);
    }
  }
  return d2sigma * P_M;
}

std::string DefaultOutputName(const UpcGridConfig& cfg) {
  return "UPC_CrossSection_" + cfg.system + "_" +
         std::to_string(static_cast<int>(std::round(cfg.sqrt_s_NN_GeV))) + "GeV_" +
         cfg.trigger + "_Optimized.root";
}

}  // namespace

UpcGridSummary GenerateUpcRhoPrimeGridOpt(const UpcGridConfig& cfg) {
  Workspace ws;
  ws.cfg = cfg;
  if (ws.cfg.output_root.empty()) ws.cfg.output_root = DefaultOutputName(ws.cfg);
  Require(ws.cfg.n_mass_bins > 0, "n_mass_bins must be positive");
  Require(ws.cfg.n_rapidity_bins > 0, "n_rapidity_bins must be positive");
  Require(ws.cfg.n_p_bins_each_side > 0, "n_p_bins_each_side must be positive");

  std::cout << "Generating rho-prime UPC grid: system=" << ws.cfg.system
            << " sqrt_s_NN=" << ws.cfg.sqrt_s_NN_GeV << " GeV trigger=" << ws.cfg.trigger
            << "\n";
  SetupSystem(ws);
  LoadProbabilityHistograms(ws);
  InitMassNormalization(ws);
  PrecomputeGlobals(ws);
  std::cout << "Precomputing optimized omega LUT (50 log-spaced points)\n";
  const OmegaLUT omega_lut = PrecomputeOmegaLUT(ws, 50);

  const double dM = (ws.cfg.mass_max_GeV - ws.cfg.mass_min_GeV) / ws.cfg.n_mass_bins;
  const double dy = (ws.cfg.rapidity_max - ws.cfg.rapidity_min) / ws.cfg.n_rapidity_bins;
  const double p_edge = (ws.cfg.n_p_bins_each_side + 0.5) * ws.dp_global;

  TH2D h2_dM_dy("h2_dM_dy", "d^{2}#sigma/(dM dy) (mb/GeV);M;y", ws.cfg.n_mass_bins,
                ws.cfg.mass_min_GeV, ws.cfg.mass_max_GeV, ws.cfg.n_rapidity_bins,
                ws.cfg.rapidity_min, ws.cfg.rapidity_max);
  TH1D h_dM("h_dM", "d#sigma/dM (mb/GeV);M", ws.cfg.n_mass_bins, ws.cfg.mass_min_GeV,
            ws.cfg.mass_max_GeV);
  TH1D h_dy("h_dy", "d#sigma/dy (mb);y", ws.cfg.n_rapidity_bins, ws.cfg.rapidity_min,
            ws.cfg.rapidity_max);
  TH2D h_sigma("h2_px_py_sigma", "d^{2}#sigma/dpxdpy;px;py", ws.n_bins_total, -p_edge,
               p_edge, ws.n_bins_total, -p_edge, p_edge);
  TH2D h_Ax2("h2_px_py_Ax2", "Ax2;px;py", ws.n_bins_total, -p_edge, p_edge, ws.n_bins_total,
             -p_edge, p_edge);
  TH2D h_Ay2("h2_px_py_Ay2", "Ay2;px;py", ws.n_bins_total, -p_edge, p_edge, ws.n_bins_total,
             -p_edge, p_edge);
  TH2D h_Re("h2_px_py_ReAxAy", "ReAxAy;px;py", ws.n_bins_total, -p_edge, p_edge,
            ws.n_bins_total, -p_edge, p_edge);
  TH2D h_Im("h2_px_py_ImAxAy", "ImAxAy;px;py", ws.n_bins_total, -p_edge, p_edge,
            ws.n_bins_total, -p_edge, p_edge);
  std::vector<std::unique_ptr<TH2D>> h_sigma_y;
  std::vector<std::unique_ptr<TH2D>> h_Ax2_y;
  std::vector<std::unique_ptr<TH2D>> h_Ay2_y;
  std::vector<std::unique_ptr<TH2D>> h_Re_y;
  std::vector<std::unique_ptr<TH2D>> h_Im_y;
  h_sigma_y.reserve(ws.cfg.n_rapidity_bins);
  h_Ax2_y.reserve(ws.cfg.n_rapidity_bins);
  h_Ay2_y.reserve(ws.cfg.n_rapidity_bins);
  h_Re_y.reserve(ws.cfg.n_rapidity_bins);
  h_Im_y.reserve(ws.cfg.n_rapidity_bins);
  for (int iy = 0; iy < ws.cfg.n_rapidity_bins; ++iy) {
    h_sigma_y.push_back(std::make_unique<TH2D>(
        Form("h2_px_py_sigma_y%d", iy), "d^{2}#sigma/dpxdpy;px;py", ws.n_bins_total, -p_edge,
        p_edge, ws.n_bins_total, -p_edge, p_edge));
    h_Ax2_y.push_back(std::make_unique<TH2D>(Form("h2_px_py_Ax2_y%d", iy), "Ax2;px;py",
                                             ws.n_bins_total, -p_edge, p_edge, ws.n_bins_total,
                                             -p_edge, p_edge));
    h_Ay2_y.push_back(std::make_unique<TH2D>(Form("h2_px_py_Ay2_y%d", iy), "Ay2;px;py",
                                             ws.n_bins_total, -p_edge, p_edge, ws.n_bins_total,
                                             -p_edge, p_edge));
    h_Re_y.push_back(std::make_unique<TH2D>(Form("h2_px_py_ReAxAy_y%d", iy), "ReAxAy;px;py",
                                            ws.n_bins_total, -p_edge, p_edge, ws.n_bins_total,
                                            -p_edge, p_edge));
    h_Im_y.push_back(std::make_unique<TH2D>(Form("h2_px_py_ImAxAy_y%d", iy), "ImAxAy;px;py",
                                            ws.n_bins_total, -p_edge, p_edge, ws.n_bins_total,
                                            -p_edge, p_edge));
  }

  const int n_threads = MaxThreads();
  std::vector<std::unique_ptr<ThreadLocalData>> thread_data;
  for (int i = 0; i < n_threads; ++i) {
    thread_data.push_back(
        std::make_unique<ThreadLocalData>(ws.n_bins_total, ws.cfg.n_rapidity_bins));
  }

#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
  for (int i = 1; i <= ws.cfg.n_mass_bins; ++i) {
    for (int j = 1; j <= ws.cfg.n_rapidity_bins; ++j) {
      const int tid = ThreadNum();
      const double M = h2_dM_dy.GetXaxis()->GetBinCenter(i);
      const double y = h2_dM_dy.GetYaxis()->GetBinCenter(j);
      if (y < -kRapidityMirrorEps) continue;
      const double val = D2SigmaDMdyOpt(ws, omega_lut, M, y, dM, dy, j - 1, *thread_data[tid]);
#ifdef _OPENMP
#pragma omp critical
#endif
      {
        h2_dM_dy.SetBinContent(i, j, val);
        const int j_mirror = h2_dM_dy.GetYaxis()->FindBin(-y);
        if (j_mirror != j) h2_dM_dy.SetBinContent(i, j_mirror, val);
        const double y_weight = (j_mirror != j) ? 2.0 : 1.0;
        h_dM.SetBinContent(i, h_dM.GetBinContent(i) + y_weight * val * dy);
        h_dy.SetBinContent(j, h_dy.GetBinContent(j) + val * dM);
        if (j_mirror != j) h_dy.SetBinContent(j_mirror, h_dy.GetBinContent(j_mirror) + val * dM);
      }
    }
  }

  UpcGridSummary summary;
  summary.output_root = ws.cfg.output_root;
  for (const auto& td : thread_data) {
    summary.spatial_integral_mb += td->sum_spatial;
    summary.momentum_integral_mb += td->sum_momentum;
    for (int iy_bin = 0; iy_bin < ws.cfg.n_rapidity_bins; ++iy_bin) {
      for (int k = 0; k < ws.n_bins_total * ws.n_bins_total; ++k) {
        const int ix = k / ws.n_bins_total + 1;
        const int iy = k % ws.n_bins_total + 1;
        const double x = h_sigma.GetXaxis()->GetBinCenter(ix);
        const double y = h_sigma.GetYaxis()->GetBinCenter(iy);
        h_sigma_y[iy_bin]->Fill(x, y, td->hist_sigma[iy_bin][k]);
        h_Ax2_y[iy_bin]->Fill(x, y, td->hist_Ax2[iy_bin][k]);
        h_Ay2_y[iy_bin]->Fill(x, y, td->hist_Ay2[iy_bin][k]);
        h_Re_y[iy_bin]->Fill(x, y, td->hist_ReAxAy[iy_bin][k]);
        h_Im_y[iy_bin]->Fill(x, y, td->hist_ImAxAy[iy_bin][k]);
        h_sigma.Fill(x, y, td->hist_sigma[iy_bin][k]);
        h_Ax2.Fill(x, y, td->hist_Ax2[iy_bin][k]);
        h_Ay2.Fill(x, y, td->hist_Ay2[iy_bin][k]);
        h_Re.Fill(x, y, td->hist_ReAxAy[iy_bin][k]);
        h_Im.Fill(x, y, td->hist_ImAxAy[iy_bin][k]);
      }
    }
  }

  for (int iy_bin = 0; iy_bin < ws.cfg.n_rapidity_bins; ++iy_bin) {
    const double y = h2_dM_dy.GetYaxis()->GetBinCenter(iy_bin + 1);
    if (y <= kRapidityMirrorEps) continue;
    const int iy_mirror = h2_dM_dy.GetYaxis()->FindBin(-y) - 1;
    if (iy_mirror < 0 || iy_mirror >= ws.cfg.n_rapidity_bins || iy_mirror == iy_bin) continue;

    for (int ix = 1; ix <= ws.n_bins_total; ++ix) {
      const int ix_mirror = ws.n_bins_total - ix + 1;
      const double px = h_sigma.GetXaxis()->GetBinCenter(ix);
      for (int iy = 1; iy <= ws.n_bins_total; ++iy) {
        const int iy_momentum_mirror = ws.n_bins_total - iy + 1;
        const double py = h_sigma.GetYaxis()->GetBinCenter(iy);

        const double sigma = h_sigma_y[iy_bin]->GetBinContent(ix_mirror, iy_momentum_mirror);
        const double Ax2 = h_Ax2_y[iy_bin]->GetBinContent(ix_mirror, iy_momentum_mirror);
        const double Ay2 = h_Ay2_y[iy_bin]->GetBinContent(ix_mirror, iy_momentum_mirror);
        const double Re = h_Re_y[iy_bin]->GetBinContent(ix_mirror, iy_momentum_mirror);
        const double Im = h_Im_y[iy_bin]->GetBinContent(ix_mirror, iy_momentum_mirror);

        h_sigma_y[iy_mirror]->SetBinContent(ix, iy, sigma);
        h_Ax2_y[iy_mirror]->SetBinContent(ix, iy, Ax2);
        h_Ay2_y[iy_mirror]->SetBinContent(ix, iy, Ay2);
        h_Re_y[iy_mirror]->SetBinContent(ix, iy, Re);
        h_Im_y[iy_mirror]->SetBinContent(ix, iy, Im);

        h_sigma.Fill(px, py, sigma);
        h_Ax2.Fill(px, py, Ax2);
        h_Ay2.Fill(px, py, Ay2);
        h_Re.Fill(px, py, Re);
        h_Im.Fill(px, py, Im);
      }
    }
  }

  const double half_momentum_integral = summary.momentum_integral_mb;
  const double full_momentum_integral =
      h_sigma.Integral() * ws.dp_global * ws.dp_global;
  if (half_momentum_integral > 0.0) {
    const double mirror_scale = full_momentum_integral / half_momentum_integral;
    summary.spatial_integral_mb *= mirror_scale;
    summary.momentum_integral_mb = full_momentum_integral;
  }

  std::unique_ptr<TFile> out(TFile::Open(ws.cfg.output_root.c_str(), "RECREATE"));
  Require(out && !out->IsZombie(), "cannot create UPC grid ROOT file");
  h2_dM_dy.Write();
  h_dM.Write();
  h_dy.Write();
  h_sigma.Write();
  h_Ax2.Write();
  h_Ay2.Write();
  h_Re.Write();
  h_Im.Write();
  for (int iy = 0; iy < ws.cfg.n_rapidity_bins; ++iy) {
    h_sigma_y[iy]->Write();
    h_Ax2_y[iy]->Write();
    h_Ay2_y[iy]->Write();
    h_Re_y[iy]->Write();
    h_Im_y[iy]->Write();
  }
  TParameter<int>("n_rapidity_bins", ws.cfg.n_rapidity_bins).Write();
  TParameter<int>("probability_file_loaded", ws.probability_file_loaded ? 1 : 0).Write();
  TParameter<int>("probability_fallback_used", ws.probability_fallback_used ? 1 : 0).Write();
  TParameter<int>("require_probability_file", ws.cfg.require_probability_file ? 1 : 0).Write();
  TNamed("probability_file_path", ws.probability_file_path.c_str()).Write();
  TParameter<int>("local_grid_n", kNLoc).Write();
  TParameter<double>("local_grid_box_fm", kLLocFm).Write();
  TParameter<double>("local_grid_dx_fm", kDxLocFm).Write();
  TParameter<int>("omega_lut_points", 50).Write();
  out->Close();

  std::cout << "UPC grid written: " << summary.output_root
            << "\nspatial integral=" << summary.spatial_integral_mb
            << " mb momentum integral=" << summary.momentum_integral_mb << " mb\n";
  return summary;
}

}  // namespace rhoprime

#pragma once

#include <TH2D.h>
#include <TLorentzVector.h>
#include <TRandom3.h>

#include <complex>
#include <memory>
#include <string>
#include <vector>

namespace rhoprime {

constexpr double kMPi = 0.13957;
constexpr double kMRho = 0.77526;
constexpr double kGRho = 0.1491;
constexpr double kMSigma = 0.500;
constexpr double kGSigma = 0.400;

struct ComplexVector3 {
  std::complex<double> x{0.0, 0.0};
  std::complex<double> y{0.0, 0.0};
  std::complex<double> z{0.0, 0.0};

  ComplexVector3 operator+(const ComplexVector3& o) const;
  ComplexVector3& operator+=(const ComplexVector3& o);
};

struct HistSampler2D {
  std::unique_ptr<TH2D> hist;
  std::vector<double> cdf;
  double total_sum = 0.0;

  void Build(TH2D* h);
  bool Sample(TRandom3& rng, double& x, double& y) const;
  double Interpolate(double x, double y) const;
};

struct BoseGeneratorConfig {
  std::string input_root;
  std::string output_root;
  long long n_events = 100000;
  unsigned int random_seed = 42;
  long long decay_norm_trials_per_mass = 20000;
  bool bose_symmetrize = true;
  bool unweighted_events = false;
  int unweighting_mode = 1;  // 1: exact accept-reject, 2: resample with replacement, 3: weighted reservoir
  long long unweighting_trials = 200000;
  double unweighting_safety_factor = 1.25;
  int output_compression_level = 1;
  int mode1_threads = 1;
};

struct BoseGenerationSummary {
  long long requested = 0;
  long long filled = 0;
  long long rejected = 0;
  long long overweight_candidates = 0;
  double unweighting_envelope = 0.0;
  std::string output_root;
};

class EventGeneratorBose {
 public:
  explicit EventGeneratorBose(BoseGeneratorConfig cfg);

  void Initialize();
  BoseGenerationSummary Generate();

 private:
  struct Candidate;

  Candidate MakeCandidate(TRandom3& rng) const;
  double EstimateDecayWeightNorm(double M, TRandom3& rng, long long trials) const;
  double EstimateUnweightingEnvelope(TRandom3& rng, long long trials);
  void BuildDecayWeightNormalization();
  double CandidateWeight(const Candidate& c) const;
  double RawCandidateWeight(const Candidate& c) const;
  void LoadHistograms();

  BoseGeneratorConfig cfg_;
  bool initialized_ = false;

  HistSampler2D sampler_M_y_;
  std::vector<HistSampler2D> sampler_px_py_;
  std::vector<std::unique_ptr<TH2D>> h_pol_xx_;
  std::vector<std::unique_ptr<TH2D>> h_pol_yy_;
  std::vector<std::unique_ptr<TH2D>> h_pol_rexy_;
  std::vector<std::unique_ptr<TH2D>> h_pol_imxy_;
  int n_y_bins_ = 0;
  double cross_section_mb_ = 0.0;
  std::vector<double> decay_weight_norm_by_mass_bin_;
};

}  // namespace rhoprime

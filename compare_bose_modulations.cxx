#include <TFile.h>
#include <TLorentzVector.h>
#include <TParameter.h>
#include <TTree.h>
#include <TVector2.h>
#include <TVector3.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double kMRho = 0.77526;
constexpr double kMSigma = 0.500;
constexpr double kPi = 3.14159265358979323846;
constexpr std::array<double, 7> kMassEdges{{0.55828, 0.9, 1.2, 1.5, 1.8, 2.2, 3.0}};
constexpr std::array<double, 7> kPtEdges{{0.0, 0.02, 0.04, 0.06, 0.08, 0.12, 0.25}};

struct Harmonics {
  long long count = 0;
  std::array<double, 5> cos_sum{};
  std::array<double, 5> sin_sum{};
  std::array<double, 5> cos2_sum{};
  std::array<double, 5> sin2_sum{};

  void Fill(double phi) {
    ++count;
    for (int n = 1; n <= 4; ++n) {
      const double c = std::cos(n * phi);
      const double s = std::sin(n * phi);
      cos_sum[n] += c;
      sin_sum[n] += s;
      cos2_sum[n] += c * c;
      sin2_sum[n] += s * s;
    }
  }

  double Cos(int n) const { return count ? 2.0 * cos_sum[n] / count : 0.0; }
  double Sin(int n) const { return count ? 2.0 * sin_sum[n] / count : 0.0; }
  double CosError(int n) const {
    if (count < 2) return 0.0;
    const double mean = cos_sum[n] / count;
    const double variance = std::max(0.0, cos2_sum[n] / count - mean * mean);
    return 2.0 * std::sqrt(variance / count);
  }
  double SinError(int n) const {
    if (count < 2) return 0.0;
    const double mean = sin_sum[n] / count;
    const double variance = std::max(0.0, sin2_sum[n] / count - mean * mean);
    return 2.0 * std::sqrt(variance / count);
  }
};

struct AnglePair {
  float a = 0.0F;
  float b = 0.0F;
};

struct MethodResults {
  Harmonics phi_a;
  Harmonics phi_b;
  Harmonics dphi;
  std::vector<AnglePair> pairs;

  void Fill(double a, double b) {
    phi_a.Fill(a);
    phi_b.Fill(b);
    dphi.Fill(TVector2::Phi_mpi_pi(a - b));
    pairs.push_back({static_cast<float>(a), static_cast<float>(b)});
  }
};

struct SampleResults {
  std::string label;
  bool bose = false;
  long long events = 0;
  long long valid = 0;
  long long best_matches_tag = 0;
  MethodResults tagged;
  MethodResults best;
  Harmonics all_pairs_phi;
  Harmonics all_matchings_dphi;
  std::array<Harmonics, 6> tagged_dphi_mass;
  std::array<Harmonics, 6> best_dphi_mass;
  std::array<Harmonics, 6> tagged_dphi_pt;
  std::array<Harmonics, 6> best_dphi_pt;
};

template <size_t N>
int FindBin(double value, const std::array<double, N>& edges) {
  for (size_t i = 0; i + 1 < edges.size(); ++i) {
    if (value >= edges[i] && value < edges[i + 1]) return static_cast<int>(i);
  }
  return value == edges.back() ? static_cast<int>(edges.size() - 2) : -1;
}

double PairPhi(const TLorentzVector& plus, const TLorentzVector& minus,
               const TVector3& x_axis, const TVector3& y_axis) {
  const TLorentzVector pair = plus + minus;
  TLorentzVector plus_rest = plus;
  plus_rest.Boost(-pair.BoostVector());
  return std::atan2(plus_rest.Vect().Dot(y_axis), plus_rest.Vect().Dot(x_axis));
}

int SignCos(int n, double phi, double setting) {
  return std::cos(n * (phi - setting)) >= 0.0 ? 1 : -1;
}

double Correlation(const std::vector<AnglePair>& pairs, int n, double a, double b) {
  double sum = 0.0;
  for (const auto& pair : pairs) {
    sum += SignCos(n, pair.a, a) * SignCos(n, pair.b, b);
  }
  return pairs.empty() ? 0.0 : sum / pairs.size();
}

std::pair<double, double> SymmetricChshScan(const std::vector<AnglePair>& pairs, int n) {
  double best_s = 0.0;
  double best_theta = 0.0;
  for (int i = 0; i <= 180; ++i) {
    const double theta = (kPi / n) * i / 180.0;
    const double e_ab = Correlation(pairs, n, 0.0, theta);
    const double e_abp = Correlation(pairs, n, 0.0, -theta);
    const double e_apb = Correlation(pairs, n, 2.0 * theta, theta);
    const double e_apbp = Correlation(pairs, n, 2.0 * theta, -theta);
    const double s = std::abs(e_ab + e_abp + e_apb - e_apbp);
    if (s > best_s) {
      best_s = s;
      best_theta = theta;
    }
  }
  return {best_s, best_theta};
}

SampleResults Analyze(const char* path, const char* label) {
  TFile file(path, "READ");
  if (file.IsZombie()) throw std::runtime_error(std::string("cannot open ") + path);
  auto* tree = dynamic_cast<TTree*>(file.Get("Events"));
  if (!tree) throw std::runtime_error(std::string("missing Events in ") + path);

  SampleResults result;
  result.label = label;
  result.events = tree->GetEntries();
  if (auto* p = dynamic_cast<TParameter<int>*>(file.Get("bose_symmetrize"))) {
    result.bose = p->GetVal() != 0;
  }
  result.tagged.pairs.reserve(result.events);
  result.best.pairs.reserve(result.events);

  TLorentzVector *pi1 = nullptr, *pi2 = nullptr, *pi3 = nullptr, *pi4 = nullptr;
  int charge[4]{};
  int parent[4]{};
  tree->SetBranchAddress("pi1", &pi1);
  tree->SetBranchAddress("pi2", &pi2);
  tree->SetBranchAddress("pi3", &pi3);
  tree->SetBranchAddress("pi4", &pi4);
  tree->SetBranchAddress("pi_charge", charge);
  tree->SetBranchAddress("pi_parent", parent);

  for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry) {
    tree->GetEntry(entry);
    TLorentzVector* pion[4] = {pi1, pi2, pi3, pi4};
    int plus[2]{}, minus[2]{};
    int nplus = 0, nminus = 0;
    int tagged_rho_plus = -1, tagged_rho_minus = -1;
    int tagged_sigma_plus = -1, tagged_sigma_minus = -1;
    for (int i = 0; i < 4; ++i) {
      if (charge[i] > 0) plus[nplus++] = i;
      if (charge[i] < 0) minus[nminus++] = i;
      if (parent[i] == 1 && charge[i] > 0) tagged_rho_plus = i;
      if (parent[i] == 1 && charge[i] < 0) tagged_rho_minus = i;
      if (parent[i] == 2 && charge[i] > 0) tagged_sigma_plus = i;
      if (parent[i] == 2 && charge[i] < 0) tagged_sigma_minus = i;
    }
    if (nplus != 2 || nminus != 2 || tagged_rho_plus < 0 || tagged_rho_minus < 0 ||
        tagged_sigma_plus < 0 || tagged_sigma_minus < 0) {
      continue;
    }

    const TLorentzVector p4pi = *pi1 + *pi2 + *pi3 + *pi4;
    TVector3 x_axis = p4pi.Vect();
    x_axis.SetZ(0.0);
    if (x_axis.Mag() < 1.0e-12) continue;
    x_axis = x_axis.Unit();
    const TVector3 y_axis = TVector3(0.0, 0.0, 1.0).Cross(x_axis).Unit();

    double phi[2][2]{};
    double mass[2][2]{};
    for (int ip = 0; ip < 2; ++ip) {
      for (int im = 0; im < 2; ++im) {
        phi[ip][im] = PairPhi(*pion[plus[ip]], *pion[minus[im]], x_axis, y_axis);
        mass[ip][im] = (*pion[plus[ip]] + *pion[minus[im]]).M();
        result.all_pairs_phi.Fill(phi[ip][im]);
      }
    }
    result.all_matchings_dphi.Fill(TVector2::Phi_mpi_pi(phi[0][0] - phi[1][1]));
    result.all_matchings_dphi.Fill(TVector2::Phi_mpi_pi(phi[0][1] - phi[1][0]));

    const double tagged_rho_phi =
        PairPhi(*pion[tagged_rho_plus], *pion[tagged_rho_minus], x_axis, y_axis);
    const double tagged_sigma_phi =
        PairPhi(*pion[tagged_sigma_plus], *pion[tagged_sigma_minus], x_axis, y_axis);
    result.tagged.Fill(tagged_rho_phi, tagged_sigma_phi);
    const double tagged_dphi = TVector2::Phi_mpi_pi(tagged_rho_phi - tagged_sigma_phi);

    struct Assignment {
      int rp;
      int rm;
      int sp;
      int sm;
      double score;
    };
    std::array<Assignment, 4> assignments{{
        {0, 0, 1, 1, std::abs(mass[0][0] - kMRho) + std::abs(mass[1][1] - kMSigma)},
        {1, 1, 0, 0, std::abs(mass[1][1] - kMRho) + std::abs(mass[0][0] - kMSigma)},
        {0, 1, 1, 0, std::abs(mass[0][1] - kMRho) + std::abs(mass[1][0] - kMSigma)},
        {1, 0, 0, 1, std::abs(mass[1][0] - kMRho) + std::abs(mass[0][1] - kMSigma)},
    }};
    const auto best = *std::min_element(
        assignments.begin(), assignments.end(),
        [](const Assignment& a, const Assignment& b) { return a.score < b.score; });
    result.best.Fill(phi[best.rp][best.rm], phi[best.sp][best.sm]);
    const double best_dphi = TVector2::Phi_mpi_pi(phi[best.rp][best.rm] - phi[best.sp][best.sm]);
    const int mass_bin = FindBin(p4pi.M(), kMassEdges);
    const int pt_bin = FindBin(p4pi.Pt(), kPtEdges);
    if (mass_bin >= 0) {
      result.tagged_dphi_mass[mass_bin].Fill(tagged_dphi);
      result.best_dphi_mass[mass_bin].Fill(best_dphi);
    }
    if (pt_bin >= 0) {
      result.tagged_dphi_pt[pt_bin].Fill(tagged_dphi);
      result.best_dphi_pt[pt_bin].Fill(best_dphi);
    }
    if (plus[best.rp] == tagged_rho_plus && minus[best.rm] == tagged_rho_minus &&
        plus[best.sp] == tagged_sigma_plus && minus[best.sm] == tagged_sigma_minus) {
      ++result.best_matches_tag;
    }
    ++result.valid;
  }
  pi1 = pi2 = pi3 = pi4 = nullptr;
  return result;
}

void PrintHarmonics(const std::string& sample, const std::string& method,
                    const std::string& variable, const Harmonics& h) {
  for (int n = 1; n <= 4; ++n) {
    std::cout << "MOMENT," << sample << ',' << method << ',' << variable << ',' << n << ','
              << h.Cos(n) << ',' << h.CosError(n) << ',' << h.Sin(n) << ',' << h.SinError(n)
              << ',' << h.count << '\n';
  }
}

void PrintMethod(const SampleResults& sample, const char* name, const MethodResults& method) {
  PrintHarmonics(sample.label, name, "phi_a", method.phi_a);
  PrintHarmonics(sample.label, name, "phi_b", method.phi_b);
  PrintHarmonics(sample.label, name, "dphi", method.dphi);
  for (int n = 1; n <= 2; ++n) {
    const auto [s, theta] = SymmetricChshScan(method.pairs, n);
    const double harmonic_proxy = 2.0 * std::sqrt(2.0) * std::abs(method.dphi.Cos(n) / 2.0);
    std::cout << "BELL_LIKE," << sample.label << ',' << name << ',' << n << ',' << s << ','
              << theta << ',' << harmonic_proxy << ',' << method.pairs.size() << '\n';
  }
}

void Print(const SampleResults& sample) {
  std::cout << "SAMPLE," << sample.label << ",bose=" << sample.bose << ",events=" << sample.events
            << ",valid=" << sample.valid << ",best_matches_generator_tag_fraction="
            << (sample.valid ? static_cast<double>(sample.best_matches_tag) / sample.valid : 0.0)
            << '\n';
  PrintMethod(sample, "generator_tag", sample.tagged);
  PrintMethod(sample, "best_mass_no_truth", sample.best);
  PrintHarmonics(sample.label, "all_pairs_no_truth", "phi", sample.all_pairs_phi);
  PrintHarmonics(sample.label, "all_matchings_no_truth", "dphi", sample.all_matchings_dphi);
  for (size_t i = 0; i + 1 < kMassEdges.size(); ++i) {
    for (int n = 1; n <= 2; ++n) {
      std::cout << "BINNED," << sample.label << ",mass," << kMassEdges[i] << ','
                << kMassEdges[i + 1] << ',' << n << ',' << sample.tagged_dphi_mass[i].Cos(n)
                << ',' << sample.tagged_dphi_mass[i].CosError(n) << ','
                << sample.best_dphi_mass[i].Cos(n) << ',' << sample.best_dphi_mass[i].CosError(n)
                << ',' << sample.tagged_dphi_mass[i].count << '\n';
    }
  }
  for (size_t i = 0; i + 1 < kPtEdges.size(); ++i) {
    for (int n = 1; n <= 2; ++n) {
      std::cout << "BINNED," << sample.label << ",pt," << kPtEdges[i] << ',' << kPtEdges[i + 1]
                << ',' << n << ',' << sample.tagged_dphi_pt[i].Cos(n) << ','
                << sample.tagged_dphi_pt[i].CosError(n) << ',' << sample.best_dphi_pt[i].Cos(n)
                << ',' << sample.best_dphi_pt[i].CosError(n) << ','
                << sample.tagged_dphi_pt[i].count << '\n';
    }
  }
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 3) {
    std::cerr << "Usage: compare_bose_modulations bose.root nonbose.root\n";
    return 2;
  }
  try {
    std::cout << std::setprecision(10);
    std::cout << "MOMENT_HEADER,sample,method,variable,n,2cos,error,2sin,error,count\n";
    std::cout << "BELL_LIKE_HEADER,sample,method,n,S_symmetric_max,theta_at_max,"
                 "harmonic_visibility_proxy,count\n";
    std::cout << "BINNED_HEADER,sample,axis,low,high,n,tagged_2cos,error,best_2cos,error,count\n";
    Print(Analyze(argv[1], "bose"));
    Print(Analyze(argv[2], "nonbose"));
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "compare_bose_modulations: " << e.what() << '\n';
    return 1;
  }
}

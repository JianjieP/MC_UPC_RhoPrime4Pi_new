// Compare weighted reference events to unit-weight/resampled events.
// Usage:
//   root -l -b -q 'review/scripts/compare_event_samples.C("weighted.root","unit.root")'

#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>
#include <TMath.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

TTree* GetTree(TFile& file) {
  auto* tree = dynamic_cast<TTree*>(file.Get("Events"));
  if (!tree) throw std::runtime_error("missing TTree Events");
  return tree;
}

struct Moments {
  double sum_w = 0.0;
  double mean = 0.0;
  double rms = 0.0;
};

Moments FillHist(TTree* tree, const char* branch, TH1D& hist, bool weighted) {
  double x = 0.0;
  double w = 1.0;
  tree->SetBranchAddress(branch, &x);
  tree->SetBranchAddress("event_weight", &w);
  double sum_w = 0.0;
  double sum_x = 0.0;
  double sum_x2 = 0.0;
  for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);
    const double ww = weighted ? w : 1.0;
    if (!std::isfinite(x) || !std::isfinite(ww) || ww <= 0.0) continue;
    hist.Fill(x, ww);
    sum_w += ww;
    sum_x += ww * x;
    sum_x2 += ww * x * x;
  }
  Moments m;
  m.sum_w = sum_w;
  if (sum_w > 0.0) {
    m.mean = sum_x / sum_w;
    m.rms = std::sqrt(std::max(0.0, sum_x2 / sum_w - m.mean * m.mean));
  }
  return m;
}

double ShapeChi2(const TH1D& ref, const TH1D& test, int& ndof) {
  TH1D r(ref), t(test);
  if (r.Integral() > 0.0) r.Scale(1.0 / r.Integral());
  if (t.Integral() > 0.0) t.Scale(1.0 / t.Integral());
  double chi2 = 0.0;
  ndof = 0;
  for (int i = 1; i <= r.GetNbinsX(); ++i) {
    const double a = r.GetBinContent(i);
    const double b = t.GetBinContent(i);
    const double ea = ref.GetBinError(i) / std::max(1.0, ref.Integral());
    const double eb = test.GetBinError(i) / std::max(1.0, test.Integral());
    const double err2 = ea * ea + eb * eb;
    if (err2 <= 0.0 && (a == 0.0 && b == 0.0)) continue;
    chi2 += (a - b) * (a - b) / std::max(err2, 1.0e-12);
    ++ndof;
  }
  if (ndof > 1) --ndof;
  return chi2;
}

void CompareOne(TTree* ref_tree, TTree* test_tree, const char* branch, int nbins, double lo,
                double hi) {
  TH1D h_ref(Form("h_ref_%s", branch), branch, nbins, lo, hi);
  TH1D h_test(Form("h_test_%s", branch), branch, nbins, lo, hi);
  h_ref.Sumw2();
  h_test.Sumw2();
  const auto ref_m = FillHist(ref_tree, branch, h_ref, true);
  const auto test_m = FillHist(test_tree, branch, h_test, false);
  int ndof = 0;
  const double chi2 = ShapeChi2(h_ref, h_test, ndof);
  const double rel_mean = ref_m.mean != 0.0 ? (test_m.mean - ref_m.mean) / ref_m.mean : 0.0;
  std::cout << branch << " ref_mean=" << ref_m.mean << " test_mean=" << test_m.mean
            << " rel_mean=" << rel_mean << " ref_rms=" << ref_m.rms
            << " test_rms=" << test_m.rms << " shape_chi2_ndof="
            << (ndof > 0 ? chi2 / ndof : 0.0) << " ndof=" << ndof << "\n";
}

}  // namespace

void compare_event_samples(const char* weighted_ref, const char* unit_weight_test) {
  TFile f_ref(weighted_ref, "READ");
  TFile f_test(unit_weight_test, "READ");
  if (f_ref.IsZombie()) throw std::runtime_error(std::string("cannot open ") + weighted_ref);
  if (f_test.IsZombie()) throw std::runtime_error(std::string("cannot open ") + unit_weight_test);
  TTree* ref_tree = GetTree(f_ref);
  TTree* test_tree = GetTree(f_test);
  std::cout << "COMPARE weighted_ref=" << weighted_ref << " unit_test=" << unit_weight_test
            << " ref_entries=" << ref_tree->GetEntries()
            << " test_entries=" << test_tree->GetEntries() << "\n";
  CompareOne(ref_tree, test_tree, "M_rhoprime", 40, 0.55, 3.0);
  CompareOne(ref_tree, test_tree, "y_rhoprime", 40, -1.0, 1.0);
  CompareOne(ref_tree, test_tree, "pt_rhoprime", 40, 0.0, 0.25);
}

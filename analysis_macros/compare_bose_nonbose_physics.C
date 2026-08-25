// Compare Bose and non-Bose rho' -> pi+ pi- pi+ pi- samples.
//
// Run with all events:
//   root -l -b -q 'analysis_macros/compare_bose_nonbose_physics.C()'
//
// Quick test with 100000 events per sample:
//   root -l -b -q 'analysis_macros/compare_bose_nonbose_physics.C(-1,-1,-1,100000)'
//
// The angular observable is defined identically for both samples and does not
// use pi_parent truth information:
//   * x axis: rho' transverse-momentum direction in the lab
//   * y axis: z_lab x x
//   * boost pi+ to the pi+pi- pair rest frame
//   * test all disjoint rho/sigma assignments and minimize
//       |m(rho candidate)-m_rho| + |m(sigma candidate)-m_sigma|
//   * use the pi+ azimuth from that best rho candidate
//
// Thus every accepted event has one entry in the moment profiles.

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TParameter.h>
#include <TProfile.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TVector2.h>
#include <TVector3.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace BoseNonBoseComparison {

constexpr double kRhoMass = 0.77526;
constexpr double kSigmaMass = 0.500;
constexpr double kPairMassMin = 0.6;
constexpr double kPairMassMax = 0.9;

struct EventData {
  double event_weight = 1.0;
  TLorentzVector* pi1 = nullptr;
  TLorentzVector* pi2 = nullptr;
  TLorentzVector* pi3 = nullptr;
  TLorentzVector* pi4 = nullptr;
  int pi_charge[4] = {1, -1, 1, -1};
};

struct SampleResult {
  std::string tag;
  std::string label;
  Long64_t entries = 0;
  Long64_t processed = 0;
  Long64_t valid_weight_events = 0;
  Long64_t angular_events = 0;
  double sum_weight = 0.0;
  double sum_weight2 = 0.0;

  std::unique_ptr<TH1D> mass;
  std::unique_ptr<TH1D> rapidity;
  std::unique_ptr<TH1D> pt;
  std::unique_ptr<TProfile> two_cos_phi;
  std::unique_ptr<TProfile> two_cos_phi_vs_pt;

  // Additional useful physics observables.
  std::unique_ptr<TProfile> two_cos_2phi;
  std::unique_ptr<TProfile> two_cos_2phi_vs_pt;
  std::unique_ptr<TH1D> all_os_pair_mass;
  std::unique_ptr<TH1D> best_rho_mass;
  std::unique_ptr<TH1D> best_sigma_mass;
  std::unique_ptr<TH1D> best_rho_phi;
  std::unique_ptr<TH1D> best_rho_sigma_dphi;
  std::unique_ptr<TH1D> selected_pair_multiplicity;
};

void RequireBranch(TTree& tree, const char* name) {
  if (!tree.GetBranch(name)) {
    throw std::runtime_error(std::string("missing branch '") + name +
                             "' in tree " + tree.GetName());
  }
}

void BindBranches(TTree& tree, EventData& ev) {
  tree.SetBranchStatus("*", 0);
  for (const char* name : {"event_weight", "pi1", "pi2", "pi3", "pi4"}) {
    RequireBranch(tree, name);
    tree.SetBranchStatus(name, 1);
  }
  tree.SetBranchAddress("event_weight", &ev.event_weight);
  tree.SetBranchAddress("pi1", &ev.pi1);
  tree.SetBranchAddress("pi2", &ev.pi2);
  tree.SetBranchAddress("pi3", &ev.pi3);
  tree.SetBranchAddress("pi4", &ev.pi4);

  if (tree.GetBranch("pi_charge")) {
    tree.SetBranchStatus("pi_charge", 1);
    tree.SetBranchAddress("pi_charge", ev.pi_charge);
  }
}

double PairPhi(const TLorentzVector& pi_plus, const TLorentzVector& pi_minus,
               const TVector3& x_axis, const TVector3& y_axis) {
  const TLorentzVector pair = pi_plus + pi_minus;
  TLorentzVector pi_plus_rest = pi_plus;
  pi_plus_rest.Boost(-pair.BoostVector());
  const TVector3 p = pi_plus_rest.Vect();
  return std::atan2(p.Dot(y_axis), p.Dot(x_axis));
}

template <class HistType, class... Args>
std::unique_ptr<HistType> MakeHist(const std::string& name, Args&&... args) {
  auto hist =
      std::make_unique<HistType>(name.c_str(), std::forward<Args>(args)...);
  hist->SetDirectory(nullptr);
  hist->Sumw2();
  return hist;
}

SampleResult MakeResult(const std::string& tag, const std::string& label) {
  SampleResult r;
  r.tag = tag;
  r.label = label;
  r.mass = MakeHist<TH1D>(
      tag + "_rho_prime_mass", ";M_{4#pi} [GeV];Weighted events", 250, 0.5, 3.0);
  r.rapidity = MakeHist<TH1D>(
      tag + "_rho_prime_y", ";y_{4#pi};Weighted events", 96, -2.4, 2.4);
  r.pt = MakeHist<TH1D>(
      tag + "_rho_prime_pt", ";p_{T,4#pi} [GeV];Weighted events", 40, 0.0, 0.2);

  r.two_cos_phi = MakeHist<TProfile>(
      tag + "_two_cos_phi",
      ";angular moment;2#LTcos#phi#GT", 1, 0.0, 1.0, -2.0, 2.0);
  r.two_cos_phi->GetXaxis()->SetBinLabel(1, "best #rho/#sigma pairing");
  r.two_cos_phi_vs_pt = MakeHist<TProfile>(
      tag + "_two_cos_phi_vs_pt",
      ";p_{T,4#pi} [GeV];2#LTcos#phi#GT", 20, 0.0, 0.2, -2.0, 2.0);

  r.two_cos_2phi = MakeHist<TProfile>(
      tag + "_two_cos_2phi",
      ";angular moment;2#LTcos(2#phi)#GT", 1, 0.0, 1.0, -2.0, 2.0);
  r.two_cos_2phi->GetXaxis()->SetBinLabel(1, "best #rho/#sigma pairing");
  r.two_cos_2phi_vs_pt = MakeHist<TProfile>(
      tag + "_two_cos_2phi_vs_pt",
      ";p_{T,4#pi} [GeV];2#LTcos(2#phi)#GT", 20, 0.0, 0.2, -2.0, 2.0);

  r.all_os_pair_mass = MakeHist<TH1D>(
      tag + "_all_os_pair_mass",
      ";m_{#pi^{+}#pi^{-}} [GeV];Weighted pairs", 100, 0.2, 1.5);
  r.best_rho_mass = MakeHist<TH1D>(
      tag + "_best_rho_mass",
      ";m_{#pi^{+}#pi^{-}}^{best #rho} [GeV];Weighted events", 100, 0.2, 1.5);
  r.best_sigma_mass = MakeHist<TH1D>(
      tag + "_best_sigma_mass",
      ";m_{#pi^{+}#pi^{-}}^{best #sigma} [GeV];Weighted events", 100, 0.2, 1.5);
  r.best_rho_phi = MakeHist<TH1D>(
      tag + "_best_rho_phi",
      ";#phi_{#pi^{+}}^{best #rho} [rad];Weighted events", 64, -M_PI, M_PI);
  r.best_rho_sigma_dphi = MakeHist<TH1D>(
      tag + "_best_rho_sigma_dphi",
      ";#Delta#phi_{#rho,#sigma} [rad];Weighted events", 64, -M_PI, M_PI);
  r.selected_pair_multiplicity = MakeHist<TH1D>(
      tag + "_selected_pair_multiplicity",
      ";N_{OS pairs}(0.6 < m_{#pi#pi} < 0.9 GeV);Weighted events",
      5, -0.5, 4.5);
  return r;
}

SampleResult ProcessSample(const char* input_file, const std::string& tag,
                           const std::string& label, Long64_t max_events) {
  std::unique_ptr<TFile> file(TFile::Open(input_file, "READ"));
  if (!file || file->IsZombie()) {
    throw std::runtime_error(std::string("cannot open input file: ") + input_file);
  }
  auto* tree = dynamic_cast<TTree*>(file->Get("Events"));
  if (!tree) {
    throw std::runtime_error(std::string("missing Events tree in: ") + input_file);
  }

  SampleResult r = MakeResult(tag, label);
  r.entries = tree->GetEntries();
  r.processed =
      (max_events > 0) ? std::min(max_events, r.entries) : r.entries;
  EventData ev;
  BindBranches(*tree, ev);

  std::cout << "\n[" << label << "] processing " << r.processed << " / "
            << r.entries << " events from\n  " << input_file << std::endl;
  const Long64_t progress_step = std::max<Long64_t>(1, r.processed / 20);

  for (Long64_t i = 0; i < r.processed; ++i) {
    tree->GetEntry(i);
    if (i % progress_step == 0) {
      std::cout << "  " << (100 * i / std::max<Long64_t>(1, r.processed))
                << "%\r" << std::flush;
    }
    const double w = ev.event_weight;
    if (!std::isfinite(w) || !(w > 0.0) || !ev.pi1 || !ev.pi2 ||
        !ev.pi3 || !ev.pi4) {
      continue;
    }
    ++r.valid_weight_events;
    r.sum_weight += w;
    r.sum_weight2 += w * w;

    TLorentzVector* pions[4] = {ev.pi1, ev.pi2, ev.pi3, ev.pi4};
    const TLorentzVector rho_prime = *ev.pi1 + *ev.pi2 + *ev.pi3 + *ev.pi4;
    r.mass->Fill(rho_prime.M(), w);
    r.rapidity->Fill(rho_prime.Rapidity(), w);
    r.pt->Fill(rho_prime.Pt(), w);

    std::vector<int> plus;
    std::vector<int> minus;
    for (int j = 0; j < 4; ++j) {
      if (ev.pi_charge[j] > 0) plus.push_back(j);
      if (ev.pi_charge[j] < 0) minus.push_back(j);
    }
    if (plus.size() != 2 || minus.size() != 2) continue;

    TVector3 x_axis(rho_prime.Px(), rho_prime.Py(), 0.0);
    if (x_axis.Mag2() < 1.0e-24) continue;
    x_axis = x_axis.Unit();
    const TVector3 y_axis = TVector3(0.0, 0.0, 1.0).Cross(x_axis).Unit();

    double pair_mass[2][2] = {};
    double pair_phi[2][2] = {};
    int selected_pairs = 0;

    for (int ip = 0; ip < 2; ++ip) {
      for (int im = 0; im < 2; ++im) {
        const TLorentzVector& pplus = *pions[plus[ip]];
        const TLorentzVector& pminus = *pions[minus[im]];
        pair_mass[ip][im] = (pplus + pminus).M();
        pair_phi[ip][im] = PairPhi(pplus, pminus, x_axis, y_axis);
        r.all_os_pair_mass->Fill(pair_mass[ip][im], w);
        if (pair_mass[ip][im] > kPairMassMin &&
            pair_mass[ip][im] < kPairMassMax) {
          ++selected_pairs;
        }
      }
    }
    r.selected_pair_multiplicity->Fill(selected_pairs, w);

    // Choose the disjoint rho/sigma assignment closest to their pole masses.
    struct Assignment {
      int rp;
      int rm;
      int sp;
      int sm;
    };
    const Assignment assignments[4] = {
        {0, 0, 1, 1}, {1, 1, 0, 0}, {0, 1, 1, 0}, {1, 0, 0, 1}};
    int best = 0;
    double best_score = 1.0e100;
    for (int ia = 0; ia < 4; ++ia) {
      const Assignment& a = assignments[ia];
      const double score =
          std::abs(pair_mass[a.rp][a.rm] - kRhoMass) +
          std::abs(pair_mass[a.sp][a.sm] - kSigmaMass);
      if (score < best_score) {
        best_score = score;
        best = ia;
      }
    }
    const Assignment& a = assignments[best];
    const double rho_phi = pair_phi[a.rp][a.rm];
    const double event_two_cos_phi = 2.0 * std::cos(rho_phi);
    const double event_two_cos_2phi = 2.0 * std::cos(2.0 * rho_phi);
    r.two_cos_phi->Fill(0.5, event_two_cos_phi, w);
    r.two_cos_phi_vs_pt->Fill(rho_prime.Pt(), event_two_cos_phi, w);
    r.two_cos_2phi->Fill(0.5, event_two_cos_2phi, w);
    r.two_cos_2phi_vs_pt->Fill(rho_prime.Pt(), event_two_cos_2phi, w);
    ++r.angular_events;

    r.best_rho_mass->Fill(pair_mass[a.rp][a.rm], w);
    r.best_sigma_mass->Fill(pair_mass[a.sp][a.sm], w);
    r.best_rho_phi->Fill(rho_phi, w);
    r.best_rho_sigma_dphi->Fill(
        TVector2::Phi_mpi_pi(rho_phi - pair_phi[a.sp][a.sm]), w);
  }
  std::cout << "  100%\n  valid-weight events: " << r.valid_weight_events
            << ", angular events: " << r.angular_events
            << ", sum weights: " << r.sum_weight << std::endl;

  // ROOT may otherwise try to manage branch objects while closing the file.
  ev.pi1 = ev.pi2 = ev.pi3 = ev.pi4 = nullptr;
  tree->ResetBranchAddresses();
  return r;
}

double EffectiveEntries(const SampleResult& r) {
  return (r.sum_weight2 > 0.0)
             ? r.sum_weight * r.sum_weight / r.sum_weight2
             : 0.0;
}

void SetHistStyle(TH1* hist, int color, int marker) {
  hist->SetLineColor(color);
  hist->SetMarkerColor(color);
  hist->SetMarkerStyle(marker);
  hist->SetMarkerSize(1.0);
  hist->SetLineWidth(2);
  hist->GetXaxis()->SetTitleSize(0.048);
  hist->GetYaxis()->SetTitleSize(0.048);
  hist->GetXaxis()->SetLabelSize(0.042);
  hist->GetYaxis()->SetLabelSize(0.042);
  hist->GetXaxis()->SetTitleOffset(1.05);
  hist->GetYaxis()->SetTitleOffset(1.35);
}

std::unique_ptr<TH1D> DensityClone(const TH1D& source,
                                   const std::string& name) {
  auto result = std::unique_ptr<TH1D>(
      static_cast<TH1D*>(source.Clone(name.c_str())));
  result->SetDirectory(nullptr);
  const double integral = result->Integral("width");
  if (integral > 0.0) result->Scale(1.0 / integral);
  return result;
}

std::string PhysicsTitle(const std::string& basename) {
  if (basename == "rho_prime_mass") return "M_{#rho'}";
  if (basename == "rho_prime_y") return "y_{#rho'}";
  if (basename == "rho_prime_pt") return "p_{T,#rho'}";
  if (basename == "two_cos_phi") return "2#LTcos#phi#GT";
  if (basename == "two_cos_phi_vs_pt") return "2#LTcos#phi#GT vs. p_{T,#rho'}";
  if (basename == "two_cos_2phi") return "2#LTcos(2#phi)#GT";
  if (basename == "two_cos_2phi_vs_pt") return "2#LTcos(2#phi)#GT vs. p_{T,#rho'}";
  if (basename == "all_os_pair_mass") return "All opposite-sign pair masses";
  if (basename == "best_rho_mass") return "Best #rho-candidate mass";
  if (basename == "best_sigma_mass") return "Best #sigma-candidate mass";
  if (basename == "best_rho_phi") return "Best #rho-candidate #phi";
  if (basename == "best_rho_sigma_dphi") return "Best #rho-#sigma #Delta#phi";
  if (basename == "selected_pair_multiplicity") return "Selected-pair multiplicity";
  return basename;
}

std::string NormalizedYAxis(const std::string& basename) {
  if (basename == "rho_prime_mass") return "(1/N) dN/dM_{#rho'} [GeV^{-1}]";
  if (basename == "rho_prime_pt") return "(1/N) dN/dp_{T,#rho'} [GeV^{-1}]";
  if (basename == "rho_prime_y") return "(1/N) dN/dy_{#rho'}";
  if (basename == "best_rho_mass" || basename == "best_sigma_mass" ||
      basename == "all_os_pair_mass") {
    return "(1/N) dN/dm [GeV^{-1}]";
  }
  return "Normalized density";
}

void DrawHeader(const std::string& title, const char* extra = nullptr) {
  TLatex text;
  text.SetNDC();
  text.SetTextFont(42);
  // Standard centered plot title.
  text.SetTextAlign(22);
  text.SetTextSize(0.045);
  text.DrawLatex(0.55, 0.95, title.c_str());
  if (extra) {
    text.SetTextSize(0.028);
    text.DrawLatex(0.55, 0.895, extra);
  }
  // Place the CMS label immediately inside the upper-left plot-frame corner.
  // Keeping it inside avoids ROOT clipping it for TProfile drawing modes.
  text.SetTextAlign(13);
  text.SetTextSize(0.040);
  text.DrawLatex(0.15, 0.805, "#bf{CMS} #it{Work in progress}");
}

void SaveCanvas(TCanvas& canvas, const std::string& output_dir,
                const std::string& basename) {
  canvas.SaveAs((output_dir + "/" + basename + ".png").c_str());
  canvas.SaveAs((output_dir + "/" + basename + ".pdf").c_str());
}

void DrawShapeComparison(const TH1D& bose_source, const TH1D& nonbose_source,
                         const std::string& output_dir,
                         const std::string& basename, const char* extra = nullptr) {
  auto bose = DensityClone(bose_source, basename + "_bose_draw");
  auto nonbose = DensityClone(nonbose_source, basename + "_nonbose_draw");
  SetHistStyle(bose.get(), kRed + 1, 20);
  SetHistStyle(nonbose.get(), kBlue + 1, 24);
  bose->GetYaxis()->SetTitle(NormalizedYAxis(basename).c_str());
  nonbose->GetYaxis()->SetTitle(NormalizedYAxis(basename).c_str());
  const double ymax = 1.25 * std::max(bose->GetMaximum(), nonbose->GetMaximum());
  bose->SetMaximum(ymax > 0.0 ? ymax : 1.0);
  bose->SetMinimum(0.0);

  TCanvas canvas((basename + "_canvas").c_str(), "", 900, 760);
  canvas.SetLeftMargin(0.14);
  canvas.SetRightMargin(0.04);
  canvas.SetBottomMargin(0.12);
  canvas.SetTopMargin(0.18);
  canvas.SetTicks(1, 1);
  bose->Draw("HIST");
  nonbose->Draw("HIST SAME");
  TLegend legend(0.66, 0.68, 0.92, 0.80);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.AddEntry(bose.get(), "Bose", "l");
  legend.AddEntry(nonbose.get(), "Non-Bose", "l");
  legend.Draw();
  DrawHeader(PhysicsTitle(basename), extra);
  SaveCanvas(canvas, output_dir, basename);
}

void SetMomentRange(TProfile& bose, TProfile& nonbose) {
  double ymin = 0.0;
  double ymax = 0.0;
  bool found = false;
  for (TProfile* p : {&bose, &nonbose}) {
    for (int bin = 1; bin <= p->GetNbinsX(); ++bin) {
      if (p->GetBinEntries(bin) <= 0.0) continue;
      const double value = p->GetBinContent(bin);
      const double error = p->GetBinError(bin);
      if (!found) {
        ymin = value - error;
        ymax = value + error;
        found = true;
      } else {
        ymin = std::min(ymin, value - error);
        ymax = std::max(ymax, value + error);
      }
    }
  }
  ymin = std::min(ymin, 0.0);
  ymax = std::max(ymax, 0.0);
  double margin = 0.18 * std::max(0.02, ymax - ymin);
  bose.SetMinimum(std::max(-2.0, ymin - margin));
  bose.SetMaximum(std::min(2.0, ymax + margin));
}

void DrawMomentComparison(TProfile& bose, TProfile& nonbose,
                          const std::string& output_dir,
                          const std::string& basename, const char* extra) {
  SetHistStyle(&bose, kRed + 1, 20);
  SetHistStyle(&nonbose, kBlue + 1, 24);
  SetMomentRange(bose, nonbose);

  TCanvas canvas((basename + "_canvas").c_str(), "", 900, 760);
  canvas.SetLeftMargin(0.14);
  canvas.SetRightMargin(0.04);
  canvas.SetBottomMargin(0.12);
  canvas.SetTopMargin(0.18);
  canvas.SetTicks(1, 1);
  canvas.SetGridy();
  bose.Draw("E1");
  nonbose.Draw("E1 SAME");
  TLegend legend(0.66, 0.68, 0.92, 0.80);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.AddEntry(&bose, "Bose", "lep");
  legend.AddEntry(&nonbose, "Non-Bose", "lep");
  legend.Draw();
  DrawHeader(PhysicsTitle(basename), extra);
  SaveCanvas(canvas, output_dir, basename);
}

void WriteSample(TFile& output, SampleResult& r) {
  TDirectory* dir = output.mkdir(r.tag.c_str());
  dir->cd();
  for (TH1* hist :
       {static_cast<TH1*>(r.mass.get()),
        static_cast<TH1*>(r.rapidity.get()),
        static_cast<TH1*>(r.pt.get()),
        static_cast<TH1*>(r.two_cos_phi.get()),
        static_cast<TH1*>(r.two_cos_phi_vs_pt.get()),
        static_cast<TH1*>(r.two_cos_2phi.get()),
        static_cast<TH1*>(r.two_cos_2phi_vs_pt.get()),
        static_cast<TH1*>(r.all_os_pair_mass.get()),
        static_cast<TH1*>(r.best_rho_mass.get()),
        static_cast<TH1*>(r.best_sigma_mass.get()),
        static_cast<TH1*>(r.best_rho_phi.get()),
        static_cast<TH1*>(r.best_rho_sigma_dphi.get()),
        static_cast<TH1*>(r.selected_pair_multiplicity.get())}) {
    hist->Write();
  }
  TParameter<Long64_t>("input_entries", r.entries).Write();
  TParameter<Long64_t>("processed_entries", r.processed).Write();
  TParameter<Long64_t>("valid_weight_events", r.valid_weight_events).Write();
  TParameter<Long64_t>("angular_events", r.angular_events).Write();
  TParameter<double>("sum_event_weight", r.sum_weight).Write();
  TParameter<double>("effective_entries", EffectiveEntries(r)).Write();
  output.cd();
}

void PrintSummary(const SampleResult& r) {
  std::cout << "\n=== " << r.label << " ===\n"
            << "processed             = " << r.processed << "\n"
            << "valid weighted events = " << r.valid_weight_events << "\n"
            << "angular events        = " << r.angular_events << "\n"
            << "sum weights           = " << r.sum_weight << "\n"
            << "effective entries     = " << EffectiveEntries(r) << "\n"
            << "M(4pi) mean [GeV]     = " << r.mass->GetMean() << "\n"
            << "y(4pi) mean           = " << r.rapidity->GetMean() << "\n"
            << "pT(4pi) mean [GeV]    = " << r.pt->GetMean() << "\n"
            << "2<cos(phi)>           = " << r.two_cos_phi->GetBinContent(1)
            << " +/- " << r.two_cos_phi->GetBinError(1) << "\n"
            << "2<cos(2phi)>          = " << r.two_cos_2phi->GetBinContent(1)
            << " +/- " << r.two_cos_2phi->GetBinError(1) << std::endl;
}

}  // namespace BoseNonBoseComparison

void compare_bose_nonbose_physics(
    const char* bose_file =
        "/eos/cms/store/group/phys_heavyions/jianjie/MC_UPC_RhoPrime4Pi_new/result/events_bose_5000000.root",
    const char* nonbose_file =
        "/eos/cms/store/group/phys_heavyions/jianjie/MC_UPC_RhoPrime4Pi_new/result/events_nonbose_5000000.root",
    const char* output_dir =
        "/eos/cms/store/group/phys_heavyions/jianjie/MC_UPC_RhoPrime4Pi_new/result/bose_nonbose_physics_comparison",
    Long64_t max_events = -1) {
  using namespace BoseNonBoseComparison;

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetLegendFont(42);
  if (gSystem->mkdir(output_dir, kTRUE) != 0 &&
      gSystem->AccessPathName(output_dir)) {
    std::cerr << "ERROR: cannot create output directory " << output_dir
              << std::endl;
    return;
  }

  try {
    SampleResult bose =
        ProcessSample(bose_file, "bose", "Bose", max_events);
    SampleResult nonbose =
        ProcessSample(nonbose_file, "nonbose", "Non-Bose", max_events);

    // Five required plots.
    DrawShapeComparison(*bose.mass, *nonbose.mass, output_dir,
                        "rho_prime_mass");
    DrawShapeComparison(*bose.rapidity, *nonbose.rapidity, output_dir,
                        "rho_prime_y");
    DrawShapeComparison(*bose.pt, *nonbose.pt, output_dir,
                        "rho_prime_pt");
    DrawMomentComparison(
        *bose.two_cos_phi, *nonbose.two_cos_phi, output_dir,
        "two_cos_phi",
        "best disjoint #rho/#sigma mass assignment");
    DrawMomentComparison(
        *bose.two_cos_phi_vs_pt, *nonbose.two_cos_phi_vs_pt, output_dir,
        "two_cos_phi_vs_pt",
        "best disjoint #rho/#sigma mass assignment");

    // Additional physics comparisons.
    DrawMomentComparison(
        *bose.two_cos_2phi, *nonbose.two_cos_2phi, output_dir,
        "two_cos_2phi",
        "best disjoint #rho/#sigma mass assignment");
    DrawMomentComparison(
        *bose.two_cos_2phi_vs_pt, *nonbose.two_cos_2phi_vs_pt, output_dir,
        "two_cos_2phi_vs_pt",
        "best disjoint #rho/#sigma mass assignment");
    DrawShapeComparison(*bose.all_os_pair_mass, *nonbose.all_os_pair_mass,
                        output_dir, "all_os_pair_mass");
    DrawShapeComparison(*bose.best_rho_mass, *nonbose.best_rho_mass,
                        output_dir, "best_rho_mass",
                        "best disjoint #rho/#sigma mass assignment");
    DrawShapeComparison(*bose.best_sigma_mass, *nonbose.best_sigma_mass,
                        output_dir, "best_sigma_mass",
                        "best disjoint #rho/#sigma mass assignment");
    DrawShapeComparison(*bose.best_rho_phi, *nonbose.best_rho_phi,
                        output_dir, "best_rho_phi",
                        "best disjoint #rho/#sigma mass assignment");
    DrawShapeComparison(
        *bose.best_rho_sigma_dphi, *nonbose.best_rho_sigma_dphi,
        output_dir, "best_rho_sigma_dphi",
        "best disjoint #rho/#sigma mass assignment");
    DrawShapeComparison(
        *bose.selected_pair_multiplicity,
        *nonbose.selected_pair_multiplicity, output_dir,
        "selected_pair_multiplicity");

    std::unique_ptr<TFile> output(TFile::Open(
        (std::string(output_dir) + "/bose_nonbose_comparison.root").c_str(),
        "RECREATE"));
    if (!output || output->IsZombie()) {
      throw std::runtime_error("cannot create output ROOT file");
    }
    WriteSample(*output, bose);
    WriteSample(*output, nonbose);
    output->Write();
    output->Close();

    PrintSummary(bose);
    PrintSummary(nonbose);
    std::cout << "\nOutput directory: " << output_dir
              << "\nROOT histograms : " << output_dir
              << "/bose_nonbose_comparison.root"
              << "\nPNG and PDF plots are in the same directory." << std::endl;
  } catch (const std::exception& error) {
    std::cerr << "ERROR: " << error.what() << std::endl;
  }
}

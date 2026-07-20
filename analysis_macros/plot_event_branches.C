#include <algorithm>
#include <cmath>
#include <cctype>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

namespace {

struct EventData {
    Double_t event_weight = 1.0;
    Double_t M_rhoprime = 0.0;
    Double_t y_rhoprime = 0.0;
    Double_t pt_rhoprime = 0.0;
    Double_t m_rho = 0.0;
    Double_t m_sigma = 0.0;
    Double_t pol_xx = 0.0;
    Double_t pol_yy = 0.0;
    Double_t pol_rexy = 0.0;
    Double_t pol_imxy = 0.0;
    TLorentzVector* pi1 = nullptr;
    TLorentzVector* pi2 = nullptr;
    TLorentzVector* pi3 = nullptr;
    TLorentzVector* pi4 = nullptr;
    TLorentzVector* rho_momentum = nullptr;
    TLorentzVector* sigma_momentum = nullptr;
    Int_t pi_parent[4] = {0, 0, 0, 0};
    Int_t pi_charge[4] = {0, 0, 0, 0};
};

struct VarSpec {
    std::string name;
    std::string title;
    int bins = 100;
    bool fixed_range = false;
    double fixed_min = 0.0;
    double fixed_max = 1.0;
    std::function<void(const EventData&, std::vector<double>&)> values;
};

struct Range {
    double min = std::numeric_limits<double>::infinity();
    double max = -std::numeric_limits<double>::infinity();
    Long64_t count = 0;
};

bool EnsureDir(const std::string& dir) {
    if (gSystem->AccessPathName(dir.c_str())) {
        return gSystem->mkdir(dir.c_str(), kTRUE) == 0;
    }
    return true;
}

std::string Sanitize(const std::string& in) {
    std::string out;
    out.reserve(in.size());
    for (unsigned char c : in) {
        if (std::isalnum(c) || c == '_' || c == '-' || c == '.') out.push_back(c);
        else out.push_back('_');
    }
    return out.empty() ? "hist" : out;
}

void UpdateRange(Range& r, double x) {
    if (!std::isfinite(x)) return;
    r.min = std::min(r.min, x);
    r.max = std::max(r.max, x);
    ++r.count;
}

void ExpandedRange(const Range& r, double& lo, double& hi) {
    lo = r.min;
    hi = r.max;
    if (r.count == 0 || !std::isfinite(lo) || !std::isfinite(hi)) {
        lo = 0.0;
        hi = 1.0;
        return;
    }
    if (lo == hi) {
        const double pad = std::max(1.0e-3, std::abs(lo) * 0.1);
        lo -= pad;
        hi += pad;
        return;
    }
    const double pad = 0.08 * (hi - lo);
    lo -= pad;
    hi += pad;
}

bool HasBranch(TTree* tree, const char* name) {
    return tree && tree->GetBranch(name);
}

void SetBranchIfPresent(TTree* tree, const char* name, void* address) {
    if (HasBranch(tree, name)) tree->SetBranchAddress(name, address);
    else std::cerr << "Warning: missing branch " << name << std::endl;
}

void BindBranches(TTree* tree, EventData& ev) {
    SetBranchIfPresent(tree, "event_weight", &ev.event_weight);
    SetBranchIfPresent(tree, "M_rhoprime", &ev.M_rhoprime);
    SetBranchIfPresent(tree, "y_rhoprime", &ev.y_rhoprime);
    SetBranchIfPresent(tree, "pt_rhoprime", &ev.pt_rhoprime);
    SetBranchIfPresent(tree, "m_rho", &ev.m_rho);
    SetBranchIfPresent(tree, "m_sigma", &ev.m_sigma);
    SetBranchIfPresent(tree, "pol_xx", &ev.pol_xx);
    SetBranchIfPresent(tree, "pol_yy", &ev.pol_yy);
    SetBranchIfPresent(tree, "pol_rexy", &ev.pol_rexy);
    SetBranchIfPresent(tree, "pol_imxy", &ev.pol_imxy);
    SetBranchIfPresent(tree, "pi1", &ev.pi1);
    SetBranchIfPresent(tree, "pi2", &ev.pi2);
    SetBranchIfPresent(tree, "pi3", &ev.pi3);
    SetBranchIfPresent(tree, "pi4", &ev.pi4);
    SetBranchIfPresent(tree, "rho_momentum", &ev.rho_momentum);
    SetBranchIfPresent(tree, "sigma_momentum", &ev.sigma_momentum);
    SetBranchIfPresent(tree, "pi_parent", ev.pi_parent);
    SetBranchIfPresent(tree, "pi_charge", ev.pi_charge);
}

void AddScalar(std::vector<VarSpec>& vars, const std::string& name,
               const std::string& title, double EventData::* member) {
    VarSpec v;
    v.name = name;
    v.title = title + ";value;Events";
    v.bins = 120;
    v.values = [member](const EventData& ev, std::vector<double>& out) {
        out.push_back(ev.*member);
    };
    vars.push_back(v);
}

void AddArray(std::vector<VarSpec>& vars, const std::string& name,
              const std::string& title, Int_t (EventData::* member)[4], int index,
              bool all_values) {
    VarSpec v;
    v.name = name;
    v.title = title + ";value;Entries";
    v.bins = 7;
    v.fixed_range = true;
    v.fixed_min = (name.find("charge") != std::string::npos) ? -3.5 : -0.5;
    v.fixed_max = (name.find("charge") != std::string::npos) ? 3.5 : 6.5;
    v.values = [member, index, all_values](const EventData& ev, std::vector<double>& out) {
        const Int_t* arr = (ev.*member);
        if (all_values) {
            for (int i = 0; i < 4; ++i) out.push_back(arr[i]);
        } else {
            out.push_back(arr[index]);
        }
    };
    vars.push_back(v);
}

void AddLV(std::vector<VarSpec>& vars, const std::string& name,
           const std::string& label, TLorentzVector* EventData::* member) {
    auto add = [&](const std::string& suffix, const std::string& axis,
                   std::function<double(const TLorentzVector&)> f,
                   bool fixed = false, double lo = 0.0, double hi = 1.0) {
        VarSpec v;
        v.name = name + "_" + suffix;
        v.title = label + " " + suffix + ";" + axis + ";Events";
        v.bins = 120;
        v.fixed_range = fixed;
        v.fixed_min = lo;
        v.fixed_max = hi;
        v.values = [member, f](const EventData& ev, std::vector<double>& out) {
            const TLorentzVector* p = ev.*member;
            if (p) out.push_back(f(*p));
        };
        vars.push_back(v);
    };

    add("E", "E [GeV]", [](const TLorentzVector& p) { return p.E(); });
    add("Px", "p_{x} [GeV]", [](const TLorentzVector& p) { return p.Px(); });
    add("Py", "p_{y} [GeV]", [](const TLorentzVector& p) { return p.Py(); });
    add("Pz", "p_{z} [GeV]", [](const TLorentzVector& p) { return p.Pz(); });
    add("Pt", "p_{T} [GeV]", [](const TLorentzVector& p) { return p.Pt(); });
    add("P", "|p| [GeV]", [](const TLorentzVector& p) { return p.P(); });
    add("Eta", "#eta", [](const TLorentzVector& p) { return p.Eta(); });
    add("Phi", "#phi [rad]", [](const TLorentzVector& p) { return p.Phi(); },
        true, -M_PI, M_PI);
    add("Rapidity", "y", [](const TLorentzVector& p) { return p.Rapidity(); });
    add("M", "mass [GeV]", [](const TLorentzVector& p) { return p.M(); });
}

std::vector<VarSpec> MakeVariables() {
    std::vector<VarSpec> vars;
    AddScalar(vars, "event_weight", "event_weight", &EventData::event_weight);
    AddScalar(vars, "M_rhoprime", "M_{#rho'}", &EventData::M_rhoprime);
    AddScalar(vars, "y_rhoprime", "y_{#rho'}", &EventData::y_rhoprime);
    AddScalar(vars, "pt_rhoprime", "p_{T,#rho'} [GeV]", &EventData::pt_rhoprime);
    AddScalar(vars, "m_rho", "m_{#rho} [GeV]", &EventData::m_rho);
    AddScalar(vars, "m_sigma", "m_{#sigma} [GeV]", &EventData::m_sigma);
    AddScalar(vars, "pol_xx", "pol_xx", &EventData::pol_xx);
    AddScalar(vars, "pol_yy", "pol_yy", &EventData::pol_yy);
    AddScalar(vars, "pol_rexy", "pol_rexy", &EventData::pol_rexy);
    AddScalar(vars, "pol_imxy", "pol_imxy", &EventData::pol_imxy);

    AddArray(vars, "pi_parent_all", "pi_parent all pions", &EventData::pi_parent, 0, true);
    AddArray(vars, "pi_charge_all", "pi_charge all pions", &EventData::pi_charge, 0, true);
    for (int i = 0; i < 4; ++i) {
        AddArray(vars, "pi_parent_" + std::to_string(i + 1),
                 "pi_parent[" + std::to_string(i) + "]", &EventData::pi_parent, i, false);
        AddArray(vars, "pi_charge_" + std::to_string(i + 1),
                 "pi_charge[" + std::to_string(i) + "]", &EventData::pi_charge, i, false);
    }

    AddLV(vars, "pi1", "pi1", &EventData::pi1);
    AddLV(vars, "pi2", "pi2", &EventData::pi2);
    AddLV(vars, "pi3", "pi3", &EventData::pi3);
    AddLV(vars, "pi4", "pi4", &EventData::pi4);
    AddLV(vars, "rho_momentum", "#rho momentum", &EventData::rho_momentum);
    AddLV(vars, "sigma_momentum", "#sigma momentum", &EventData::sigma_momentum);
    return vars;
}

void StyleHist(TH1D* h, int color) {
    h->SetLineColor(color);
    h->SetMarkerColor(color);
    h->SetLineWidth(2);
    h->GetYaxis()->SetTitleOffset(1.25);
}

void DrawOne(TH1D* h, const std::string& label, const std::string& path) {
    TCanvas c("c_one_branch", "c_one_branch", 900, 750);
    c.SetLeftMargin(0.13);
    c.SetBottomMargin(0.12);
    c.SetGrid();
    StyleHist(h, kBlue + 1);
    h->Draw("HIST");
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.16, 0.92, label.c_str());
    c.SaveAs(path.c_str());
}

void DrawCompare(TH1D* hb, TH1D* hn, const std::string& path) {
    std::unique_ptr<TH1D> b(static_cast<TH1D*>(hb->Clone("b_norm")));
    std::unique_ptr<TH1D> n(static_cast<TH1D*>(hn->Clone("n_norm")));
    b->SetDirectory(nullptr);
    n->SetDirectory(nullptr);
    if (b->Integral() > 0.0) b->Scale(1.0 / b->Integral());
    if (n->Integral() > 0.0) n->Scale(1.0 / n->Integral());
    StyleHist(b.get(), kRed + 1);
    StyleHist(n.get(), kBlue + 1);
    b->SetTitle((std::string(hb->GetTitle()) + ";normalized entries").c_str());
    n->SetTitle(b->GetTitle());
    const double ymax = 1.25 * std::max(b->GetMaximum(), n->GetMaximum());
    b->SetMaximum(ymax > 0.0 ? ymax : 1.0);

    TCanvas c("c_compare_branch", "c_compare_branch", 900, 750);
    c.SetLeftMargin(0.13);
    c.SetBottomMargin(0.12);
    c.SetGrid();
    b->Draw("HIST");
    n->Draw("HIST SAME");
    TLegend leg(0.62, 0.76, 0.88, 0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry(b.get(), "Bose", "l");
    leg.AddEntry(n.get(), "Non-Bose", "l");
    leg.Draw();
    c.SaveAs(path.c_str());
}

std::map<std::string, std::unique_ptr<TH1D>> ProcessSample(
    const char* inputFile, const std::string& sampleName, const std::string& outputDir,
    const std::vector<VarSpec>& vars, Long64_t maxEvents, bool useEventWeights) {

    std::map<std::string, std::unique_ptr<TH1D>> hists;
    std::unique_ptr<TFile> file(TFile::Open(inputFile, "READ"));
    if (!file || file->IsZombie()) {
        std::cerr << "Error: cannot open " << inputFile << std::endl;
        return hists;
    }
    TTree* tree = dynamic_cast<TTree*>(file->Get("Events"));
    if (!tree) {
        std::cerr << "Error: no Events tree in " << inputFile << std::endl;
        return hists;
    }

    const Long64_t nEntries = tree->GetEntries();
    const Long64_t nProcess = (maxEvents > 0) ? std::min(maxEvents, nEntries) : nEntries;
    std::cout << "Scanning " << sampleName << ": " << nProcess << " / " << nEntries << " events" << std::endl;

    EventData ev;
    BindBranches(tree, ev);
    std::vector<Range> ranges(vars.size());
    std::vector<double> values;
    for (Long64_t i = 0; i < nProcess; ++i) {
        tree->GetEntry(i);
        for (size_t iv = 0; iv < vars.size(); ++iv) {
            if (vars[iv].fixed_range) continue;
            values.clear();
            vars[iv].values(ev, values);
            for (double x : values) UpdateRange(ranges[iv], x);
        }
    }

    for (size_t iv = 0; iv < vars.size(); ++iv) {
        double lo = vars[iv].fixed_min;
        double hi = vars[iv].fixed_max;
        if (!vars[iv].fixed_range) ExpandedRange(ranges[iv], lo, hi);
        std::string hname = sampleName + "_" + vars[iv].name;
        auto h = std::make_unique<TH1D>(hname.c_str(), vars[iv].title.c_str(),
                                        vars[iv].bins, lo, hi);
        h->Sumw2();
        h->SetDirectory(nullptr);
        hists[vars[iv].name] = std::move(h);
    }

    std::cout << "Filling " << sampleName << " histograms" << std::endl;
    for (Long64_t i = 0; i < nProcess; ++i) {
        tree->GetEntry(i);
        const double w = (useEventWeights && std::isfinite(ev.event_weight) && ev.event_weight > 0.0)
                             ? ev.event_weight
                             : 1.0;
        for (const auto& v : vars) {
            values.clear();
            v.values(ev, values);
            TH1D* h = hists[v.name].get();
            const double fillWeight = (v.name == "event_weight") ? 1.0 : w;
            for (double x : values) {
                if (std::isfinite(x)) h->Fill(x, fillWeight);
            }
        }
    }

    const std::string sampleDir = outputDir + "/" + sampleName;
    EnsureDir(sampleDir);
    for (const auto& item : hists) {
        DrawOne(item.second.get(), sampleName.c_str(),
                sampleDir + "/" + Sanitize(item.first) + ".png");
    }
    return hists;
}

}  // namespace

void plot_event_branches(
    const char* boseFile =
        "/eos/cms/store/group/phys_heavyions/jianjie/MC_UPC_RhoPrime4Pi_new/result/events_bose_5000000.root",
    const char* nonboseFile =
        "/eos/cms/store/group/phys_heavyions/jianjie/MC_UPC_RhoPrime4Pi_new/result/events_nonbose_5000000.root",
    const char* outputDir =
        "/eos/cms/store/group/phys_heavyions/jianjie/MC_UPC_RhoPrime4Pi_new/result/event_branch_hists",
    Long64_t maxEvents = -1,
    bool useEventWeights = true) {

    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(1110);
    gStyle->SetTitleSize(0.045, "xyz");
    gStyle->SetLabelSize(0.04, "xyz");

    if (!EnsureDir(outputDir)) {
        std::cerr << "Error: cannot create output directory " << outputDir << std::endl;
        return;
    }
    EnsureDir(std::string(outputDir) + "/compare");

    const auto vars = MakeVariables();
    auto bose = ProcessSample(boseFile, "bose", outputDir, vars, maxEvents, useEventWeights);
    auto nonbose = ProcessSample(nonboseFile, "nonbose", outputDir, vars, maxEvents, useEventWeights);

    std::unique_ptr<TFile> fout(TFile::Open((std::string(outputDir) + "/event_branch_hists.root").c_str(),
                                            "RECREATE"));
    if (fout && !fout->IsZombie()) {
        fout->mkdir("bose");
        fout->mkdir("nonbose");
        fout->cd("bose");
        for (const auto& item : bose) item.second->Write(item.first.c_str());
        fout->cd("nonbose");
        for (const auto& item : nonbose) item.second->Write(item.first.c_str());
        fout->Close();
    }

    int nCompare = 0;
    for (const auto& v : vars) {
        auto ib = bose.find(v.name);
        auto in = nonbose.find(v.name);
        if (ib == bose.end() || in == nonbose.end()) continue;
        DrawCompare(ib->second.get(), in->second.get(),
                    std::string(outputDir) + "/compare/" + Sanitize(v.name) + ".png");
        ++nCompare;
    }

    std::cout << "Done." << std::endl;
    std::cout << "  Individual PNGs: " << outputDir << "/bose and " << outputDir << "/nonbose" << std::endl;
    std::cout << "  Comparison PNGs: " << outputDir << "/compare (" << nCompare << " plots)" << std::endl;
    std::cout << "  ROOT hist file : " << outputDir << "/event_branch_hists.root" << std::endl;
}

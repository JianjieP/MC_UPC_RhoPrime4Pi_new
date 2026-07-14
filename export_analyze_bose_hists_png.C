#include <algorithm>
#include <cctype>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <string>

#include <TCanvas.h>
#include <TClass.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TKey.h>
#include <TList.h>
#include <TObject.h>
#include <TProfile.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TH1.h>
#include <TH2.h>

namespace {

std::string SanitizeFileName(const std::string& name) {
    std::string out;
    out.reserve(name.size());
    for (unsigned char ch : name) {
        if (std::isalnum(ch) || ch == '_' || ch == '-' || ch == '.') {
            out.push_back(static_cast<char>(ch));
        } else {
            out.push_back('_');
        }
    }
    while (out.find("__") != std::string::npos) {
        out.replace(out.find("__"), 2, "_");
    }
    if (out.empty()) out = "hist";
    return out;
}

bool EnsureDir(const std::string& dir) {
    if (gSystem->AccessPathName(dir.c_str())) {
        return gSystem->mkdir(dir.c_str(), kTRUE) == 0;
    }
    return true;
}

void SetTH1YRange(TH1* h) {
    double minY = std::numeric_limits<double>::infinity();
    double maxY = -std::numeric_limits<double>::infinity();
    bool hasNegative = false;
    bool hasBin = false;

    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        const double y = h->GetBinContent(i);
        if (!std::isfinite(y)) continue;
        minY = std::min(minY, y);
        maxY = std::max(maxY, y);
        hasNegative = hasNegative || (y < 0.0);
        hasBin = true;
    }

    if (!hasBin || !std::isfinite(minY) || !std::isfinite(maxY)) return;
    if (!hasNegative && maxY > 0.0) {
        h->SetMinimum(0.0);
        h->SetMaximum(maxY * 1.25);
        return;
    }

    const double span = maxY - minY;
    const double pad = span > 0.0 ? 0.15 * span : std::max(1.0, std::abs(maxY)) * 0.15;
    h->SetMinimum(minY - pad);
    h->SetMaximum(maxY + pad);
}

void SetProfileYRange(TProfile* p) {
    double minY = std::numeric_limits<double>::infinity();
    double maxY = -std::numeric_limits<double>::infinity();
    bool hasBin = false;

    for (int i = 1; i <= p->GetNbinsX(); ++i) {
        if (p->GetBinEntries(i) <= 0.0) continue;
        const double y = p->GetBinContent(i);
        const double e = p->GetBinError(i);
        if (!std::isfinite(y) || !std::isfinite(e)) continue;
        minY = std::min(minY, y - e);
        maxY = std::max(maxY, y + e);
        hasBin = true;
    }

    if (!hasBin || !std::isfinite(minY) || !std::isfinite(maxY)) return;
    const double span = maxY - minY;
    const double pad = span > 0.0 ? 0.15 * span : std::max(0.05, std::abs(maxY) * 0.2);
    p->SetMinimum(minY - pad);
    p->SetMaximum(maxY + pad);
}

void DrawAndSave(TH1* h, const std::string& outputPath) {
    const bool is2D = h->InheritsFrom(TH2::Class()) && h->GetDimension() >= 2;
    const bool isProfile = h->InheritsFrom(TProfile::Class());

    TCanvas canvas("c_export_hist", h->GetName(), is2D ? 1000 : 900, 800);
    canvas.SetLeftMargin(0.12);
    canvas.SetBottomMargin(0.12);
    if (is2D) canvas.SetRightMargin(0.16);
    if (!is2D) canvas.SetGrid();

    h->SetLineWidth(2);
    h->SetMarkerStyle(20);
    h->SetMarkerSize(0.7);
    h->GetYaxis()->SetTitleOffset(1.25);

    if (is2D) {
        h->Draw("COLZ");
    } else if (isProfile) {
        SetProfileYRange(static_cast<TProfile*>(h));
        h->Draw("PE");
    } else {
        SetTH1YRange(h);
        h->Draw("HIST E");
    }

    canvas.SaveAs(outputPath.c_str());
}

void ExportDirectory(TDirectory* dir, const std::string& outputDir,
                     const std::string& relDir, int& saved, int& skipped) {
    TIter next(dir->GetListOfKeys());
    while (TKey* key = static_cast<TKey*>(next())) {
        std::unique_ptr<TObject> obj(key->ReadObj());
        if (!obj) {
            ++skipped;
            continue;
        }

        if (obj->InheritsFrom(TDirectory::Class())) {
            const std::string nextRel = relDir.empty() ? obj->GetName() : relDir + "/" + obj->GetName();
            ExportDirectory(static_cast<TDirectory*>(obj.get()), outputDir, nextRel, saved, skipped);
            continue;
        }

        if (!obj->InheritsFrom(TH1::Class())) {
            ++skipped;
            continue;
        }

        TH1* hist = static_cast<TH1*>(obj.get());
        hist->SetDirectory(nullptr);

        std::string baseName = relDir.empty() ? hist->GetName() : relDir + "_" + hist->GetName();
        if (key->GetCycle() > 1) baseName += "_cycle" + std::to_string(key->GetCycle());
        const std::string outputPath = outputDir + "/" + SanitizeFileName(baseName) + ".png";
        DrawAndSave(hist, outputPath);
        ++saved;
    }
}

}  // namespace

void export_analyze_bose_hists_png(
    const char* inputFile = "/eos/user/j/jianjie/STARlight/CMSSW_15_1_0_patch3/src/MC_UPC_RhoPrime4Pi_new/build/analyze_bose_5000000.root",
    const char* outputDir = "/eos/user/j/jianjie/STARlight/CMSSW_15_1_0_patch3/src/MC_UPC_RhoPrime4Pi_new/analyze_bose_5000000_hist_png") {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    std::unique_ptr<TFile> file(TFile::Open(inputFile, "READ"));
    if (!file || file->IsZombie()) {
        std::cerr << "Error: cannot open input ROOT file: " << inputFile << std::endl;
        return;
    }

    if (!EnsureDir(outputDir)) {
        std::cerr << "Error: cannot create output directory: " << outputDir << std::endl;
        return;
    }

    int saved = 0;
    int skipped = 0;
    ExportDirectory(file.get(), outputDir, "", saved, skipped);

    std::cout << "Saved " << saved << " histogram PNG files to " << outputDir << std::endl;
    std::cout << "Skipped " << skipped << " non-hist objects." << std::endl;
}

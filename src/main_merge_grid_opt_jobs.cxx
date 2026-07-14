#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TAxis.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void Require(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error(message);
}

bool FileExists(const std::string& path) {
  if (std::FILE* fp = std::fopen(path.c_str(), "r")) {
    std::fclose(fp);
    return true;
  }
  return false;
}

std::string RootPathWithSuffix(const std::string& base, const std::string& suffix) {
  Require(!base.empty(), "base output ROOT path must not be empty");
  const std::string root_ext = ".root";
  const auto slash_pos = base.find_last_of("/\\");
  const auto root_pos = base.rfind(root_ext);
  if (root_pos != std::string::npos && root_pos + root_ext.size() == base.size() &&
      (slash_pos == std::string::npos || root_pos > slash_pos)) {
    return base.substr(0, root_pos) + suffix + root_ext;
  }
  return base + suffix + root_ext;
}

int JobFieldWidth(int n_jobs) {
  int width = 2;
  for (int x = std::max(1, n_jobs - 1); x >= 100; x /= 10) ++width;
  return width;
}

std::string JobSuffix(int job_id, int n_jobs) {
  std::ostringstream out;
  out << "_split" << n_jobs << "_job" << std::setw(JobFieldWidth(n_jobs))
      << std::setfill('0') << job_id << "of" << n_jobs;
  return out.str();
}

template <typename T>
T* GetOrThrow(TFile& file, const char* name) {
  auto* obj = dynamic_cast<T*>(file.Get(name));
  if (!obj) throw std::runtime_error(std::string("missing object ") + name + " in " + file.GetName());
  return obj;
}

int ReadIntParam(TFile& file, const char* name) {
  auto* param = dynamic_cast<TParameter<int>*>(file.Get(name));
  if (!param) throw std::runtime_error(std::string("missing int parameter ") + name + " in " + file.GetName());
  return param->GetVal();
}

double ReadDoubleParam(TFile& file, const char* name) {
  auto* param = dynamic_cast<TParameter<double>*>(file.Get(name));
  if (!param) throw std::runtime_error(std::string("missing double parameter ") + name + " in " + file.GetName());
  return param->GetVal();
}

std::string ReadNamedTitle(TFile& file, const char* name) {
  auto* named = dynamic_cast<TNamed*>(file.Get(name));
  if (!named) throw std::runtime_error(std::string("missing named object ") + name + " in " + file.GetName());
  return named->GetTitle();
}

void CheckAxis(const TAxis& axis, int n_bins, double x_min, double x_max, const std::string& label) {
  const double scale = std::max({1.0, std::abs(x_min), std::abs(x_max)});
  Require(axis.GetNbins() == n_bins, label + " has inconsistent bin count");
  Require(std::abs(axis.GetXmin() - x_min) <= 1.0e-8 * scale, label + " has inconsistent minimum");
  Require(std::abs(axis.GetXmax() - x_max) <= 1.0e-8 * scale, label + " has inconsistent maximum");
}

std::unique_ptr<TH1D> MakeEmptyTH1(const TH1D& src, const char* name) {
  auto hist = std::make_unique<TH1D>(name, src.GetTitle(), src.GetXaxis()->GetNbins(),
                                    src.GetXaxis()->GetXmin(), src.GetXaxis()->GetXmax());
  hist->SetDirectory(nullptr);
  return hist;
}

std::unique_ptr<TH2D> MakeEmptyTH2(const TH2D& src, const char* name) {
  auto hist = std::make_unique<TH2D>(name, src.GetTitle(), src.GetXaxis()->GetNbins(),
                                    src.GetXaxis()->GetXmin(), src.GetXaxis()->GetXmax(),
                                    src.GetYaxis()->GetNbins(), src.GetYaxis()->GetXmin(),
                                    src.GetYaxis()->GetXmax());
  hist->SetDirectory(nullptr);
  return hist;
}

void AddChecked(TH1D& target, TH1D& src, const std::string& label) {
  CheckAxis(*src.GetXaxis(), target.GetXaxis()->GetNbins(), target.GetXaxis()->GetXmin(),
            target.GetXaxis()->GetXmax(), label + " x-axis");
  target.Add(&src);
}

void AddChecked(TH2D& target, TH2D& src, const std::string& label) {
  CheckAxis(*src.GetXaxis(), target.GetXaxis()->GetNbins(), target.GetXaxis()->GetXmin(),
            target.GetXaxis()->GetXmax(), label + " x-axis");
  CheckAxis(*src.GetYaxis(), target.GetYaxis()->GetNbins(), target.GetYaxis()->GetXmin(),
            target.GetYaxis()->GetXmax(), label + " y-axis");
  target.Add(&src);
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc < 2) {
      std::cerr << "Usage: merge_grid_opt_jobs base_output.root [n_jobs=20] [merged_output.root]\n";
      return 2;
    }

    const std::string base_output = argv[1];
    const int n_jobs = argc > 2 ? std::stoi(argv[2]) : 20;
    const std::string merged_output =
        argc > 3 ? argv[3] : RootPathWithSuffix(base_output, "_split" + std::to_string(n_jobs) + "_merged");
    Require(n_jobs > 0, "n_jobs must be positive");
    Require(!FileExists(merged_output), "refusing to overwrite existing merged output: " + merged_output);

    const std::string first_path = RootPathWithSuffix(base_output, JobSuffix(0, n_jobs));
    std::unique_ptr<TFile> first(TFile::Open(first_path.c_str(), "READ"));
    Require(first && !first->IsZombie(), "cannot open first job file: " + first_path);

    const int global_n_mass_bins = ReadIntParam(*first, "split_global_n_mass_bins");
    const double global_mass_min = ReadDoubleParam(*first, "split_global_mass_min_GeV");
    const double global_mass_max = ReadDoubleParam(*first, "split_global_mass_max_GeV");
    const int n_rapidity_bins = ReadIntParam(*first, "n_rapidity_bins");
    const int probability_file_loaded = ReadIntParam(*first, "probability_file_loaded");
    const int probability_fallback_used = ReadIntParam(*first, "probability_fallback_used");
    const int require_probability_file = ReadIntParam(*first, "require_probability_file");
    const int local_grid_n = ReadIntParam(*first, "local_grid_n");
    const double local_grid_box_fm = ReadDoubleParam(*first, "local_grid_box_fm");
    const double local_grid_dx_fm = ReadDoubleParam(*first, "local_grid_dx_fm");
    const int omega_lut_points = ReadIntParam(*first, "omega_lut_points");
    const std::string probability_file_path = ReadNamedTitle(*first, "probability_file_path");

    auto* first_h2_dMdy = GetOrThrow<TH2D>(*first, "h2_dM_dy");
    auto* first_h_dM = GetOrThrow<TH1D>(*first, "h_dM");
    auto* first_h_dy = GetOrThrow<TH1D>(*first, "h_dy");
    Require(n_rapidity_bins == first_h2_dMdy->GetYaxis()->GetNbins(),
            "n_rapidity_bins metadata does not match h2_dM_dy");
    const double rapidity_min = first_h2_dMdy->GetYaxis()->GetXmin();
    const double rapidity_max = first_h2_dMdy->GetYaxis()->GetXmax();

    TH2D h2_dMdy("h2_dM_dy", first_h2_dMdy->GetTitle(), global_n_mass_bins, global_mass_min,
                 global_mass_max, n_rapidity_bins, rapidity_min, rapidity_max);
    TH1D h_dM("h_dM", first_h_dM->GetTitle(), global_n_mass_bins, global_mass_min,
              global_mass_max);
    auto h_dy = MakeEmptyTH1(*first_h_dy, "h_dy");
    auto h_sigma = MakeEmptyTH2(*GetOrThrow<TH2D>(*first, "h2_px_py_sigma"), "h2_px_py_sigma");
    auto h_Ax2 = MakeEmptyTH2(*GetOrThrow<TH2D>(*first, "h2_px_py_Ax2"), "h2_px_py_Ax2");
    auto h_Ay2 = MakeEmptyTH2(*GetOrThrow<TH2D>(*first, "h2_px_py_Ay2"), "h2_px_py_Ay2");
    auto h_Re = MakeEmptyTH2(*GetOrThrow<TH2D>(*first, "h2_px_py_ReAxAy"), "h2_px_py_ReAxAy");
    auto h_Im = MakeEmptyTH2(*GetOrThrow<TH2D>(*first, "h2_px_py_ImAxAy"), "h2_px_py_ImAxAy");

    std::vector<std::unique_ptr<TH2D>> h_sigma_y;
    std::vector<std::unique_ptr<TH2D>> h_Ax2_y;
    std::vector<std::unique_ptr<TH2D>> h_Ay2_y;
    std::vector<std::unique_ptr<TH2D>> h_Re_y;
    std::vector<std::unique_ptr<TH2D>> h_Im_y;
    h_sigma_y.reserve(n_rapidity_bins);
    h_Ax2_y.reserve(n_rapidity_bins);
    h_Ay2_y.reserve(n_rapidity_bins);
    h_Re_y.reserve(n_rapidity_bins);
    h_Im_y.reserve(n_rapidity_bins);
    for (int iy = 0; iy < n_rapidity_bins; ++iy) {
      const std::string suffix = std::to_string(iy);
      h_sigma_y.push_back(MakeEmptyTH2(*GetOrThrow<TH2D>(*first, ("h2_px_py_sigma_y" + suffix).c_str()),
                                       ("h2_px_py_sigma_y" + suffix).c_str()));
      h_Ax2_y.push_back(MakeEmptyTH2(*GetOrThrow<TH2D>(*first, ("h2_px_py_Ax2_y" + suffix).c_str()),
                                     ("h2_px_py_Ax2_y" + suffix).c_str()));
      h_Ay2_y.push_back(MakeEmptyTH2(*GetOrThrow<TH2D>(*first, ("h2_px_py_Ay2_y" + suffix).c_str()),
                                     ("h2_px_py_Ay2_y" + suffix).c_str()));
      h_Re_y.push_back(MakeEmptyTH2(*GetOrThrow<TH2D>(*first, ("h2_px_py_ReAxAy_y" + suffix).c_str()),
                                    ("h2_px_py_ReAxAy_y" + suffix).c_str()));
      h_Im_y.push_back(MakeEmptyTH2(*GetOrThrow<TH2D>(*first, ("h2_px_py_ImAxAy_y" + suffix).c_str()),
                                    ("h2_px_py_ImAxAy_y" + suffix).c_str()));
    }
    first->Close();
    first.reset();

    std::vector<int> mass_bin_seen(global_n_mass_bins, 0);
    double spatial_integral_mb = 0.0;
    double momentum_integral_mb = 0.0;

    for (int job = 0; job < n_jobs; ++job) {
      const std::string job_path = RootPathWithSuffix(base_output, JobSuffix(job, n_jobs));
      std::unique_ptr<TFile> file(TFile::Open(job_path.c_str(), "READ"));
      Require(file && !file->IsZombie(), "cannot open job file: " + job_path);
      Require(ReadIntParam(*file, "split_job_id") == job, "job_id metadata mismatch in " + job_path);
      Require(ReadIntParam(*file, "split_n_jobs") == n_jobs, "n_jobs metadata mismatch in " + job_path);
      Require(ReadIntParam(*file, "split_global_n_mass_bins") == global_n_mass_bins,
              "global nM mismatch in " + job_path);
      Require(std::abs(ReadDoubleParam(*file, "split_global_mass_min_GeV") - global_mass_min) < 1.0e-10,
              "global mass min mismatch in " + job_path);
      Require(std::abs(ReadDoubleParam(*file, "split_global_mass_max_GeV") - global_mass_max) < 1.0e-10,
              "global mass max mismatch in " + job_path);

      const int first_mass_bin0 = ReadIntParam(*file, "split_first_mass_bin0");
      const int n_mass_bins_this = ReadIntParam(*file, "split_n_mass_bins");
      auto* h2_sub = GetOrThrow<TH2D>(*file, "h2_dM_dy");
      auto* h_dM_sub = GetOrThrow<TH1D>(*file, "h_dM");
      auto* h_dy_sub = GetOrThrow<TH1D>(*file, "h_dy");
      Require(h2_sub->GetXaxis()->GetNbins() == n_mass_bins_this,
              "split nM metadata mismatch in " + job_path);
      Require(h_dM_sub->GetXaxis()->GetNbins() == n_mass_bins_this,
              "split h_dM nM mismatch in " + job_path);
      CheckAxis(*h2_sub->GetYaxis(), n_rapidity_bins, rapidity_min, rapidity_max,
                "h2_dM_dy y-axis in " + job_path);

      for (int ix = 1; ix <= n_mass_bins_this; ++ix) {
        const int global_bin = first_mass_bin0 + ix;
        Require(global_bin >= 1 && global_bin <= global_n_mass_bins,
                "mass bin range outside global grid in " + job_path);
        const double expected_center = h2_dMdy.GetXaxis()->GetBinCenter(global_bin);
        const double actual_center = h2_sub->GetXaxis()->GetBinCenter(ix);
        Require(std::abs(expected_center - actual_center) < 1.0e-8,
                "mass bin center mismatch in " + job_path);
        ++mass_bin_seen[global_bin - 1];
        h_dM.SetBinContent(global_bin, h_dM_sub->GetBinContent(ix));
        for (int iy = 1; iy <= n_rapidity_bins; ++iy) {
          h2_dMdy.SetBinContent(global_bin, iy, h2_sub->GetBinContent(ix, iy));
        }
      }

      AddChecked(*h_dy, *h_dy_sub, "h_dy in " + job_path);
      AddChecked(*h_sigma, *GetOrThrow<TH2D>(*file, "h2_px_py_sigma"), "h2_px_py_sigma in " + job_path);
      AddChecked(*h_Ax2, *GetOrThrow<TH2D>(*file, "h2_px_py_Ax2"), "h2_px_py_Ax2 in " + job_path);
      AddChecked(*h_Ay2, *GetOrThrow<TH2D>(*file, "h2_px_py_Ay2"), "h2_px_py_Ay2 in " + job_path);
      AddChecked(*h_Re, *GetOrThrow<TH2D>(*file, "h2_px_py_ReAxAy"), "h2_px_py_ReAxAy in " + job_path);
      AddChecked(*h_Im, *GetOrThrow<TH2D>(*file, "h2_px_py_ImAxAy"), "h2_px_py_ImAxAy in " + job_path);
      for (int iy = 0; iy < n_rapidity_bins; ++iy) {
        const std::string suffix = std::to_string(iy);
        AddChecked(*h_sigma_y[iy], *GetOrThrow<TH2D>(*file, ("h2_px_py_sigma_y" + suffix).c_str()),
                   "h2_px_py_sigma_y" + suffix + " in " + job_path);
        AddChecked(*h_Ax2_y[iy], *GetOrThrow<TH2D>(*file, ("h2_px_py_Ax2_y" + suffix).c_str()),
                   "h2_px_py_Ax2_y" + suffix + " in " + job_path);
        AddChecked(*h_Ay2_y[iy], *GetOrThrow<TH2D>(*file, ("h2_px_py_Ay2_y" + suffix).c_str()),
                   "h2_px_py_Ay2_y" + suffix + " in " + job_path);
        AddChecked(*h_Re_y[iy], *GetOrThrow<TH2D>(*file, ("h2_px_py_ReAxAy_y" + suffix).c_str()),
                   "h2_px_py_ReAxAy_y" + suffix + " in " + job_path);
        AddChecked(*h_Im_y[iy], *GetOrThrow<TH2D>(*file, ("h2_px_py_ImAxAy_y" + suffix).c_str()),
                   "h2_px_py_ImAxAy_y" + suffix + " in " + job_path);
      }

      spatial_integral_mb += ReadDoubleParam(*file, "split_spatial_integral_mb");
      momentum_integral_mb += ReadDoubleParam(*file, "split_momentum_integral_mb");
      std::cout << "Merged split job " << job << "/" << n_jobs << ": " << job_path << "\n";
    }

    for (int i = 0; i < global_n_mass_bins; ++i) {
      Require(mass_bin_seen[i] == 1, "global mass bin " + std::to_string(i) +
                                      " was filled " + std::to_string(mass_bin_seen[i]) + " times");
    }

    std::unique_ptr<TFile> out(TFile::Open(merged_output.c_str(), "CREATE"));
    Require(out && !out->IsZombie(), "cannot create merged output ROOT: " + merged_output);
    h2_dMdy.Write();
    h_dM.Write();
    h_dy->Write();
    h_sigma->Write();
    h_Ax2->Write();
    h_Ay2->Write();
    h_Re->Write();
    h_Im->Write();
    for (int iy = 0; iy < n_rapidity_bins; ++iy) {
      h_sigma_y[iy]->Write();
      h_Ax2_y[iy]->Write();
      h_Ay2_y[iy]->Write();
      h_Re_y[iy]->Write();
      h_Im_y[iy]->Write();
    }
    TParameter<int>("n_rapidity_bins", n_rapidity_bins).Write();
    TParameter<int>("probability_file_loaded", probability_file_loaded).Write();
    TParameter<int>("probability_fallback_used", probability_fallback_used).Write();
    TParameter<int>("require_probability_file", require_probability_file).Write();
    TNamed("probability_file_path", probability_file_path.c_str()).Write();
    TParameter<int>("local_grid_n", local_grid_n).Write();
    TParameter<double>("local_grid_box_fm", local_grid_box_fm).Write();
    TParameter<double>("local_grid_dx_fm", local_grid_dx_fm).Write();
    TParameter<int>("omega_lut_points", omega_lut_points).Write();
    TNamed("split_base_output_root", base_output.c_str()).Write();
    TParameter<int>("split_n_jobs", n_jobs).Write();
    TParameter<int>("split_global_n_mass_bins", global_n_mass_bins).Write();
    TParameter<double>("split_global_mass_min_GeV", global_mass_min).Write();
    TParameter<double>("split_global_mass_max_GeV", global_mass_max).Write();
    TParameter<double>("merged_spatial_integral_mb", spatial_integral_mb).Write();
    TParameter<double>("merged_momentum_integral_mb", momentum_integral_mb).Write();
    out->Close();

    const int y0_bin = h_dy->GetXaxis()->FindBin(0.0);
    const double dsigma_dy = h_dy->GetBinContent(y0_bin);
    const double y_center = h_dy->GetXaxis()->GetBinCenter(y0_bin);
    const double sigma_y_window = h_dy->Integral("width");
    const double sigma_mass_integral = h_dM.Integral("width");
    const double sigma_2d_integral = h2_dMdy.Integral("width");
    std::cout << "MERGE_RESULT output=" << merged_output << " y_center=" << y_center
              << " dsigma_dy_y0_mb=" << dsigma_dy
              << " sigma_y_window_mb=" << sigma_y_window
              << " h_dM_integral_mb=" << sigma_mass_integral
              << " h2_dM_dy_integral_mb=" << sigma_2d_integral
              << " spatial_integral_mb=" << spatial_integral_mb
              << " momentum_integral_mb=" << momentum_integral_mb << "\n";
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "merge_grid_opt_jobs: " << e.what() << "\n";
    return 1;
  }
}

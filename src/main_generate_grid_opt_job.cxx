#include "RhoPrime/UpcRhoPrimeModelOpt.h"

#include <TFile.h>
#include <TNamed.h>
#include <TParameter.h>

#include <algorithm>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

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

void WriteSplitMetadata(const rhoprime::UpcGridSummary& summary, const std::string& base_output,
                        int job_id, int n_jobs, int global_n_mass_bins,
                        double global_mass_min, double global_mass_max, int first_mass_bin0,
                        int n_mass_bins_this) {
  std::unique_ptr<TFile> out(TFile::Open(summary.output_root.c_str(), "UPDATE"));
  Require(out && !out->IsZombie(), "cannot reopen job output ROOT for metadata");

  TNamed("split_base_output_root", base_output.c_str()).Write();
  TParameter<int>("split_job_id", job_id).Write();
  TParameter<int>("split_n_jobs", n_jobs).Write();
  TParameter<int>("split_global_n_mass_bins", global_n_mass_bins).Write();
  TParameter<double>("split_global_mass_min_GeV", global_mass_min).Write();
  TParameter<double>("split_global_mass_max_GeV", global_mass_max).Write();
  TParameter<int>("split_first_mass_bin0", first_mass_bin0).Write();
  TParameter<int>("split_n_mass_bins", n_mass_bins_this).Write();
  TParameter<double>("split_spatial_integral_mb", summary.spatial_integral_mb).Write();
  TParameter<double>("split_momentum_integral_mb", summary.momentum_integral_mb).Write();
  out->Close();
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc < 14) {
      std::cerr << "Usage: generate_grid_opt_job AuAu|PbPb sqrt_s_NN_GeV tag "
                   "base_output.root nM nP bCoreSteps bTailStep yHalfWidth nY "
                   "global_box_fm require_probability_file job_id [n_jobs=20]\n";
      return 2;
    }

    rhoprime::UpcGridConfig cfg;
    cfg.system = argv[1];
    cfg.sqrt_s_NN_GeV = std::stod(argv[2]);
    cfg.trigger = argv[3];
    const std::string base_output = argv[4];
    cfg.n_mass_bins = std::stoi(argv[5]);
    cfg.n_p_bins_each_side = std::stoi(argv[6]);
    cfg.b_core_steps = std::stoi(argv[7]);
    cfg.b_tail_step_fm = std::stod(argv[8]);
    const double y_half_width = std::stod(argv[9]);
    cfg.n_rapidity_bins = std::stoi(argv[10]);
    cfg.global_box_fm = std::stod(argv[11]);
    cfg.require_probability_file = (std::stoi(argv[12]) != 0);
    const int job_id = std::stoi(argv[13]);
    const int n_jobs = argc > 14 ? std::stoi(argv[14]) : 20;
    cfg.rapidity_min = -y_half_width;
    cfg.rapidity_max = y_half_width;

    Require(cfg.n_mass_bins > 0, "global nM must be positive");
    Require(n_jobs > 0, "n_jobs must be positive");
    Require(job_id >= 0 && job_id < n_jobs, "job_id must be in [0, n_jobs)");

    const int global_n_mass_bins = cfg.n_mass_bins;
    const double global_mass_min = cfg.mass_min_GeV;
    const double global_mass_max = cfg.mass_max_GeV;
    const int base_bins = global_n_mass_bins / n_jobs;
    const int extra_bins = global_n_mass_bins % n_jobs;
    const int first_mass_bin0 = job_id * base_bins + std::min(job_id, extra_bins);
    const int n_mass_bins_this = base_bins + (job_id < extra_bins ? 1 : 0);
    Require(n_mass_bins_this > 0, "this job receives no mass bins; reduce n_jobs or increase nM");

    const double dM_global = (global_mass_max - global_mass_min) / global_n_mass_bins;
    cfg.mass_min_GeV = global_mass_min + first_mass_bin0 * dM_global;
    cfg.mass_max_GeV = cfg.mass_min_GeV + n_mass_bins_this * dM_global;
    cfg.n_mass_bins = n_mass_bins_this;
    cfg.output_root = RootPathWithSuffix(base_output, JobSuffix(job_id, n_jobs));
    Require(!FileExists(cfg.output_root), "refusing to overwrite existing job output: " + cfg.output_root);

    std::cout << "Split grid job: job=" << job_id << "/" << n_jobs
              << " global_mass_bins=" << first_mass_bin0 << ".."
              << (first_mass_bin0 + n_mass_bins_this - 1)
              << " mass_range_GeV=[" << cfg.mass_min_GeV << "," << cfg.mass_max_GeV
              << "] output=" << cfg.output_root << "\n";

    const auto summary = rhoprime::GenerateUpcRhoPrimeGridOpt(cfg);
    WriteSplitMetadata(summary, base_output, job_id, n_jobs, global_n_mass_bins,
                       global_mass_min, global_mass_max, first_mass_bin0, n_mass_bins_this);

    std::cout << "JOB_RESULT job=" << job_id << " nJobs=" << n_jobs
              << " output=" << summary.output_root
              << " spatial_integral_mb=" << summary.spatial_integral_mb
              << " momentum_integral_mb=" << summary.momentum_integral_mb << "\n";
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "generate_grid_opt_job: " << e.what() << "\n";
    return 1;
  }
}

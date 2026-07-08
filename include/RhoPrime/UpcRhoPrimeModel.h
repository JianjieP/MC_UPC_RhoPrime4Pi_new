#pragma once

#include <string>

namespace rhoprime {

struct UpcGridConfig {
  std::string system = "PbPb";
  double sqrt_s_NN_GeV = 5360.0;
  std::string trigger = "NoTag";
  std::string output_root = "";

  int n_mass_bins = 30;
  double mass_min_GeV = 4.0 * 0.13957;
  double mass_max_GeV = 3.0;
  int n_rapidity_bins = 10;
  double rapidity_min = -1.0;
  double rapidity_max = 1.0;

  int n_p_bins_each_side = 40;
  double global_box_fm = 40.0;

  int b_core_steps = 38;
  double b_core_min_fm = 12.0;
  double b_core_max_fm = 50.0;
  double b_tail_min_fm = 50.0;
  double b_tail_max_notag_fm = 1000.0;
  double b_tail_max_tagged_fm = 100.0;
  double b_tail_step_fm = 4.0;
  bool require_probability_file = false;
};

struct UpcGridSummary {
  std::string output_root;
  double spatial_integral_mb = 0.0;
  double momentum_integral_mb = 0.0;
};

UpcGridSummary GenerateUpcRhoPrimeGrid(const UpcGridConfig& cfg);

}  // namespace rhoprime

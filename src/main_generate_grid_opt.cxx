#include "RhoPrime/UpcRhoPrimeModelOpt.h"

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

namespace {

template <typename T>
T* GetOrThrow(TFile& file, const char* name) {
  auto* obj = dynamic_cast<T*>(file.Get(name));
  if (!obj) throw std::runtime_error(std::string("missing histogram ") + name);
  return obj;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc < 4) {
      std::cerr << "Usage: generate_grid_opt AuAu|PbPb sqrt_s_NN_GeV tag "
                   "[output_grid.root] [nM=30] [nP=100] [bCoreSteps=38] "
                   "[bTailStep=4] [yHalfWidth=0.1] [nY=10] [global_box_fm=40] "
                   "[require_probability_file=0]\n";
      return 2;
    }

    rhoprime::UpcGridConfig cfg;
    cfg.system = argv[1];
    cfg.sqrt_s_NN_GeV = std::stod(argv[2]);
    cfg.trigger = argv[3];
    cfg.output_root = argc > 4 ? argv[4] : "";
    cfg.n_mass_bins = argc > 5 ? std::stoi(argv[5]) : cfg.n_mass_bins;
    cfg.n_p_bins_each_side = argc > 6 ? std::stoi(argv[6]) : cfg.n_p_bins_each_side;
    cfg.b_core_steps = argc > 7 ? std::stoi(argv[7]) : cfg.b_core_steps;
    cfg.b_tail_step_fm = argc > 8 ? std::stod(argv[8]) : cfg.b_tail_step_fm;
    const double y_half_width = argc > 9 ? std::stod(argv[9]) : 0.1;
    cfg.n_rapidity_bins = argc > 10 ? std::stoi(argv[10]) : cfg.n_rapidity_bins;
    cfg.global_box_fm = argc > 11 ? std::stod(argv[11]) : cfg.global_box_fm;
    cfg.require_probability_file = argc > 12 ? (std::stoi(argv[12]) != 0) : cfg.require_probability_file;
    cfg.rapidity_min = -y_half_width;
    cfg.rapidity_max = y_half_width;

    const auto summary = rhoprime::GenerateUpcRhoPrimeGridOpt(cfg);

    std::unique_ptr<TFile> file(TFile::Open(summary.output_root.c_str(), "READ"));
    if (!file || file->IsZombie()) throw std::runtime_error("cannot open output grid");
    const auto* h_dy = GetOrThrow<TH1D>(*file, "h_dy");
    const auto* h_dM = GetOrThrow<TH1D>(*file, "h_dM");
    const auto* h_dMdy = GetOrThrow<TH2D>(*file, "h2_dM_dy");
    const int y0_bin = h_dy->GetXaxis()->FindBin(0.0);
    const double dsigma_dy = h_dy->GetBinContent(y0_bin);
    const double y_center = h_dy->GetXaxis()->GetBinCenter(y0_bin);
    const double sigma_y_window = h_dy->Integral("width");
    const double sigma_mass_integral = h_dM->Integral("width");
    const double sigma_2d_integral = h_dMdy->Integral("width");

    std::cout << "RESULT system=" << cfg.system << " sqrt_s_NN_GeV=" << cfg.sqrt_s_NN_GeV
              << " trigger=" << cfg.trigger << " y_center=" << y_center
              << " dsigma_dy_y0_mb=" << dsigma_dy
              << " sigma_y_window_mb=" << sigma_y_window
              << " h_dM_integral_mb=" << sigma_mass_integral
              << " h2_dM_dy_integral_mb=" << sigma_2d_integral
              << " spatial_integral_mb=" << summary.spatial_integral_mb
              << " momentum_integral_mb=" << summary.momentum_integral_mb
              << " nM=" << cfg.n_mass_bins << " nP=" << cfg.n_p_bins_each_side
              << " bCoreSteps=" << cfg.b_core_steps << " bTailStep=" << cfg.b_tail_step_fm
              << " nY=" << cfg.n_rapidity_bins << " global_box_fm=" << cfg.global_box_fm
              << " output=" << summary.output_root << "\n";
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "generate_grid_opt: " << e.what() << "\n";
    return 1;
  }
}

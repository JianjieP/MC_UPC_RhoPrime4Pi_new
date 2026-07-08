#include "RhoPrime/EventGeneratorBose.h"

#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
  try {
    if (argc < 5) {
      std::cerr << "Usage: generate_bose input.root output.root nEvents seed "
                   "[decay_norm_trials_per_mass] [bose_symmetrize=1] "
                   "[unweighted=0] [unweighting_trials=200000] [safety=1.25] "
                   "[unweighting_mode=1: accept-reject, 2: resample, 3: reservoir] "
                   "[compression_level=1] [mode1_threads=1]\n";
      return 2;
    }

    rhoprime::BoseGeneratorConfig cfg;
    cfg.input_root = argv[1];
    cfg.output_root = argv[2];
    cfg.n_events = std::stoll(argv[3]);
    cfg.random_seed = static_cast<unsigned int>(std::stoul(argv[4]));
    if (argc > 5) cfg.decay_norm_trials_per_mass = std::stoll(argv[5]);
    if (argc > 6) cfg.bose_symmetrize = std::stoi(argv[6]) != 0;
    if (argc > 7) cfg.unweighted_events = std::stoi(argv[7]) != 0;
    if (argc > 8) cfg.unweighting_trials = std::stoll(argv[8]);
    if (argc > 9) cfg.unweighting_safety_factor = std::stod(argv[9]);
    if (argc > 10) cfg.unweighting_mode = std::stoi(argv[10]);
    if (argc > 11) cfg.output_compression_level = std::stoi(argv[11]);
    if (argc > 12) cfg.mode1_threads = std::stoi(argv[12]);

    rhoprime::EventGeneratorBose generator(cfg);
    const auto summary = generator.Generate();

    std::cout << "Generated "
              << (cfg.bose_symmetrize ? "Bose-symmetrized" : "non-Bose")
              << " rho-prime events: " << summary.filled << "/"
              << summary.requested << "\nRejected candidates: " << summary.rejected
              << "\nOverweight candidates: " << summary.overweight_candidates
              << "\nUnweighting envelope: " << summary.unweighting_envelope
              << "\nOutput ROOT: " << summary.output_root << "\n";
    return summary.filled == summary.requested ? 0 : 1;
  } catch (const std::exception& e) {
    std::cerr << "generate_bose: " << e.what() << "\n";
    return 1;
  }
}

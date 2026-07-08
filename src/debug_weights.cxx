#include "RhoPrime/EventGeneratorBose.h"

#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
  try {
    if (argc < 2) {
      std::cerr << "Usage: debug_weights input.root [output.root] [nEvents=100] "
                   "[decay_norm_trials_per_mass=1000]\n";
      return 2;
    }

    rhoprime::BoseGeneratorConfig cfg;
    cfg.input_root = argv[1];
    cfg.output_root = argc > 2 ? argv[2] : "debug_out.root";
    cfg.n_events = argc > 3 ? std::stoll(argv[3]) : 100;
    cfg.decay_norm_trials_per_mass = argc > 4 ? std::stoll(argv[4]) : 1000;

    rhoprime::EventGeneratorBose gen(cfg);
    const auto summary = gen.Generate();
    std::cout << "Generated " << summary.filled << "/" << summary.requested
              << " events, rejected " << summary.rejected
              << ", output=" << summary.output_root << "\n";
    return summary.filled == summary.requested ? 0 : 1;
  } catch (const std::exception& e) {
    std::cerr << "debug_weights: " << e.what() << "\n";
    return 1;
  }
}

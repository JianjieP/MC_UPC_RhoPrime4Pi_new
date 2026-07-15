#include "RhoPrime/EventGeneratorBose.h"

#include <TFile.h>
#include <TLorentzVector.h>
#include <TTree.h>

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unistd.h>

namespace {

struct Options {
  std::string grid;
  std::string config_card;
  std::string output_lhe;
  long long events = 0;
  unsigned int seed = 0;
  int threads = 1;
  double beam_energy_gev = 2680.0;
  int process_id = 81;
  rhoprime::BoseGeneratorConfig generator;
};

void PrintUsage(std::ostream& out) {
  out << "Usage: rhoprime_lhe --grid data/grid.root --config cards/rhoprime.card "
         "--events N --seed SEED --threads N --output cmsgrid_final.lhe\n"
      << "\n"
      << "Optional overrides:\n"
      << "  --beam-energy GEV              LHE beam energy per beam (default: 2680)\n"
      << "  --weighted                     write weighted events instead of unit-weight events\n"
      << "  --unweighting-mode MODE         1=accept-reject, 2=resample, 3=reservoir\n"
      << "  --unweighting-trials N          trials/pool size for unit-weight generation\n"
      << "  --unweighting-safety FACTOR     envelope safety factor for mode 1\n"
      << "  --decay-norm-trials N           decay normalization trials per mass bin\n";
}

[[noreturn]] void ThrowUsage(const std::string& message = "") {
  if (!message.empty()) std::cerr << "rhoprime_lhe: " << message << "\n";
  PrintUsage(std::cerr);
  throw std::invalid_argument("invalid command line");
}

std::string Trim(const std::string& input) {
  size_t first = 0;
  while (first < input.size() &&
         std::isspace(static_cast<unsigned char>(input[first]))) {
    ++first;
  }
  size_t last = input.size();
  while (last > first &&
         std::isspace(static_cast<unsigned char>(input[last - 1]))) {
    --last;
  }
  return input.substr(first, last - first);
}

std::string CanonicalKey(std::string key) {
  key = Trim(key);
  std::transform(key.begin(), key.end(), key.begin(), [](unsigned char c) {
    if (c == '-') return '_';
    return static_cast<char>(std::tolower(c));
  });
  return key;
}

std::map<std::string, std::string> ReadConfigCard(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open config card: " + path);

  std::map<std::string, std::string> values;
  std::string line;
  int line_number = 0;
  while (std::getline(in, line)) {
    ++line_number;
    const size_t comment = line.find_first_of("#;");
    if (comment != std::string::npos) line.erase(comment);
    line = Trim(line);
    if (line.empty()) continue;

    std::string key;
    std::string value;
    size_t sep = line.find('=');
    if (sep == std::string::npos) sep = line.find(':');
    if (sep != std::string::npos) {
      key = line.substr(0, sep);
      value = line.substr(sep + 1);
    } else {
      std::istringstream iss(line);
      iss >> key >> value;
      if (key.empty() || value.empty()) {
        throw std::runtime_error("cannot parse config card line " +
                                 std::to_string(line_number) + ": " + line);
      }
    }

    key = CanonicalKey(key);
    value = Trim(value);
    if (!key.empty()) values[key] = value;
  }
  return values;
}

bool ParseBool(const std::string& value) {
  std::string v = CanonicalKey(value);
  if (v == "1" || v == "true" || v == "yes" || v == "on") return true;
  if (v == "0" || v == "false" || v == "no" || v == "off") return false;
  throw std::runtime_error("invalid boolean value: " + value);
}

long long ParseLongLong(const std::string& value, const std::string& name) {
  size_t pos = 0;
  const long long parsed = std::stoll(value, &pos);
  if (pos != value.size()) throw std::runtime_error("invalid integer for " + name);
  return parsed;
}

int ParseInt(const std::string& value, const std::string& name) {
  size_t pos = 0;
  const int parsed = std::stoi(value, &pos);
  if (pos != value.size()) throw std::runtime_error("invalid integer for " + name);
  return parsed;
}

double ParseDouble(const std::string& value, const std::string& name) {
  size_t pos = 0;
  const double parsed = std::stod(value, &pos);
  if (pos != value.size()) throw std::runtime_error("invalid number for " + name);
  return parsed;
}

std::string RequiredValue(int& i, int argc, char** argv) {
  if (i + 1 >= argc) ThrowUsage(std::string("missing value for ") + argv[i]);
  return argv[++i];
}

void ApplyCardValue(Options& opts, const std::string& key,
                    const std::string& value) {
  auto& cfg = opts.generator;
  if (key == "bose_symmetrize" || key == "bose_symmetrization") {
    cfg.bose_symmetrize = ParseBool(value);
  } else if (key == "event_mode") {
    const std::string mode = CanonicalKey(value);
    if (mode == "weighted") cfg.unweighted_events = false;
    else if (mode == "unweighted" || mode == "unit_weight" || mode == "unit_weighted") cfg.unweighted_events = true;
    else throw std::runtime_error("event_mode must be weighted or unweighted");
  } else if (key == "grid_file") {
    if (opts.grid.empty()) opts.grid = value;
  } else if (key == "unweighted_events" || key == "unweighted") {
    cfg.unweighted_events = ParseBool(value);
  } else if (key == "unweighting_mode") {
    cfg.unweighting_mode = ParseInt(value, key);
  } else if (key == "unweighting_trials") {
    cfg.unweighting_trials = ParseLongLong(value, key);
  } else if (key == "unweighting_safety_factor" ||
             key == "unweighting_safety") {
    cfg.unweighting_safety_factor = ParseDouble(value, key);
  } else if (key == "decay_norm_trials_per_mass" ||
             key == "decay_norm_trials") {
    cfg.decay_norm_trials_per_mass = ParseLongLong(value, key);
  } else if (key == "output_compression_level" ||
             key == "compression_level") {
    cfg.output_compression_level = ParseInt(value, key);
  } else if (key == "beam_energy_gev") {
    opts.beam_energy_gev = ParseDouble(value, key);
  } else if (key == "process_id") {
    opts.process_id = ParseInt(value, key);
  } else if (key == "sqrt_snn_gev") {
    opts.beam_energy_gev = 0.5 * ParseDouble(value, key);
  }
}

Options ParseOptions(int argc, char** argv) {
  Options opts;
  opts.generator.bose_symmetrize = true;
  opts.generator.unweighted_events = true;
  opts.generator.unweighting_mode = 1;

  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--help" || arg == "-h") {
      PrintUsage(std::cout);
      std::exit(0);
    }
    if (arg == "--grid") {
      opts.grid = RequiredValue(i, argc, argv);
    } else if (arg == "--config") {
      opts.config_card = RequiredValue(i, argc, argv);
    } else if (arg == "--events") {
      opts.events = ParseLongLong(RequiredValue(i, argc, argv), "events");
    } else if (arg == "--seed") {
      opts.seed =
          static_cast<unsigned int>(ParseLongLong(RequiredValue(i, argc, argv), "seed"));
    } else if (arg == "--threads") {
      opts.threads = ParseInt(RequiredValue(i, argc, argv), "threads");
    } else if (arg == "--output") {
      opts.output_lhe = RequiredValue(i, argc, argv);
    } else if (arg == "--beam-energy") {
      opts.beam_energy_gev = ParseDouble(RequiredValue(i, argc, argv), "beam-energy");
    } else if (arg == "--weighted") {
      opts.generator.unweighted_events = false;
    } else if (arg == "--unweighted") {
      opts.generator.unweighted_events = true;
    } else if (arg == "--unweighting-mode") {
      opts.generator.unweighting_mode =
          ParseInt(RequiredValue(i, argc, argv), "unweighting-mode");
    } else if (arg == "--unweighting-trials") {
      opts.generator.unweighting_trials =
          ParseLongLong(RequiredValue(i, argc, argv), "unweighting-trials");
    } else if (arg == "--unweighting-safety") {
      opts.generator.unweighting_safety_factor =
          ParseDouble(RequiredValue(i, argc, argv), "unweighting-safety");
    } else if (arg == "--decay-norm-trials") {
      opts.generator.decay_norm_trials_per_mass =
          ParseLongLong(RequiredValue(i, argc, argv), "decay-norm-trials");
    } else {
      ThrowUsage("unknown option: " + arg);
    }
  }

  if (opts.config_card.empty()) ThrowUsage("--config is required");
  if (opts.output_lhe.empty()) ThrowUsage("--output is required");
  if (opts.events <= 0) ThrowUsage("--events must be positive");
  if (opts.threads <= 0) ThrowUsage("--threads must be positive");

  for (const auto& [key, value] : ReadConfigCard(opts.config_card)) {
    ApplyCardValue(opts, key, value);
  }

  if (opts.grid.empty()) ThrowUsage("--grid is required unless grid_file is set in the card");
  if (!(opts.beam_energy_gev > 0.0)) ThrowUsage("--beam-energy must be positive");

  opts.generator.input_root = opts.grid;
  opts.generator.n_events = opts.events;
  opts.generator.random_seed = opts.seed;
  opts.generator.mode1_threads = opts.threads;
  return opts;
}

std::string MakeTemporaryRootPath() {
  const char* tmpdir_env = std::getenv("TMPDIR");
  const std::string tmpdir = (tmpdir_env && tmpdir_env[0] != '\0') ? tmpdir_env : "/tmp";
  std::string pattern = tmpdir + "/rhoprime_lhe_XXXXXX.root";
  const int fd = mkstemps(pattern.data(), 5);
  if (fd < 0) throw std::runtime_error("failed to create temporary ROOT path");
  close(fd);
  return pattern;
}


std::string XmlEscape(const std::string& input) {
  std::string out;
  out.reserve(input.size());
  for (char c : input) {
    switch (c) {
      case '&':
        out += "&amp;";
        break;
      case '<':
        out += "&lt;";
        break;
      case '>':
        out += "&gt;";
        break;
      case '"':
        out += "&quot;";
        break;
      case '\'':
        out += "&apos;";
        break;
      default:
        out += c;
        break;
    }
  }
  return out;
}

void RequireBranch(TTree& tree, const char* name) {
  if (!tree.GetBranch(name)) {
    throw std::runtime_error(std::string("missing Events branch: ") + name);
  }
}

void WriteParticle(std::ostream& out, int pdg_id, const TLorentzVector& p) {
  out << pdg_id << " 1 0 0 0 0 "
      << p.Px() << ' ' << p.Py() << ' ' << p.Pz() << ' ' << p.E() << ' '
      << rhoprime::kMPi << " 0.0 9.0\n";
}

long long WriteLheFromRoot(const std::string& root_path, const std::string& lhe_path,
                           const Options& opts) {
  std::unique_ptr<TFile> file(TFile::Open(root_path.c_str(), "READ"));
  if (!file || file->IsZombie()) {
    throw std::runtime_error("cannot open intermediate ROOT file: " + root_path);
  }

  auto* tree = dynamic_cast<TTree*>(file->Get("Events"));
  if (!tree) throw std::runtime_error("missing Events tree in " + root_path);

  RequireBranch(*tree, "event_weight");
  RequireBranch(*tree, "pi1");
  RequireBranch(*tree, "pi2");
  RequireBranch(*tree, "pi3");
  RequireBranch(*tree, "pi4");

  double event_weight = 1.0;
  TLorentzVector* pi1 = nullptr;
  TLorentzVector* pi2 = nullptr;
  TLorentzVector* pi3 = nullptr;
  TLorentzVector* pi4 = nullptr;
  tree->SetBranchAddress("event_weight", &event_weight);
  tree->SetBranchAddress("pi1", &pi1);
  tree->SetBranchAddress("pi2", &pi2);
  tree->SetBranchAddress("pi3", &pi3);
  tree->SetBranchAddress("pi4", &pi4);

  std::ofstream out(lhe_path);
  if (!out) throw std::runtime_error("cannot create LHE file: " + lhe_path);
  out << std::scientific << std::setprecision(10);

  out << "<LesHouchesEvents version=\"1.0\">\n";
  out << "<header>\n";
  out << "<rhoprime generator=\"rhoprime_lhe\" grid=\"" << XmlEscape(opts.grid)
      << "\" config=\"" << XmlEscape(opts.config_card) << "\" events=\""
      << opts.events << "\" seed=\"" << opts.seed << "\" unit_weight=\""
      << (opts.generator.unweighted_events ? 1 : 0) << "\" />\n";
  out << "</header>\n";
  out << "<init>\n";
  out << "22 22 " << opts.beam_energy_gev << ' ' << opts.beam_energy_gev
      << " 0 0 0 0 3 1\n";
  out << "1.0 0.0 3.0 " << opts.process_id << "\n";
  out << "</init>\n";

  const long long entries = tree->GetEntries();
  for (long long i = 0; i < entries; ++i) {
    tree->GetEntry(i);
    if (!pi1 || !pi2 || !pi3 || !pi4) {
      throw std::runtime_error("failed to read pion four-vectors from Events tree");
    }

    const double lhe_weight = opts.generator.unweighted_events ? 1.0 : event_weight;
    out << "<event>\n";
    out << "4 " << opts.process_id << ' ' << lhe_weight << " -1.0 -1.0 -1.0\n";
    WriteParticle(out, 211, *pi1);
    WriteParticle(out, -211, *pi2);
    WriteParticle(out, 211, *pi3);
    WriteParticle(out, -211, *pi4);
    out << "</event>\n";
  }

  out << "</LesHouchesEvents>\n";
  if (!out) throw std::runtime_error("failed while writing LHE file: " + lhe_path);
  return entries;
}

}  // namespace

int main(int argc, char** argv) {
  std::string intermediate_root;
  try {
    Options opts = ParseOptions(argc, argv);
    intermediate_root = MakeTemporaryRootPath();

    opts.generator.output_root = intermediate_root;
    rhoprime::EventGeneratorBose generator(opts.generator);
    const auto summary = generator.Generate();
    if (summary.filled != summary.requested) {
      throw std::runtime_error("generated " + std::to_string(summary.filled) +
                               "/" + std::to_string(summary.requested) +
                               " events in the intermediate ROOT file");
    }

    const long long written = WriteLheFromRoot(intermediate_root, opts.output_lhe, opts);
    if (written != opts.events) {
      throw std::runtime_error("wrote " + std::to_string(written) + "/" +
                               std::to_string(opts.events) + " LHE events");
    }

    std::remove(intermediate_root.c_str());
    std::cout << "Generated rho-prime LHE events: " << written << "\n"
              << "Output LHE: " << opts.output_lhe << "\n";
    return 0;
  } catch (const std::exception& e) {
    if (!intermediate_root.empty()) {
      std::remove(intermediate_root.c_str());
    }
    std::cerr << "rhoprime_lhe: " << e.what() << "\n";
    return 1;
  }
}

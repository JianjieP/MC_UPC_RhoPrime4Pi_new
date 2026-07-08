// analyze_nonbose.cxx — spin alignment of rho meson in non-Bose events
// Uses pi_parent to identify pi+pi- from rho decay directly
// Axes: z = lab-z, x = rho' pT direction (no boost needed)

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TVector2.h>
#include <TParameter.h>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <string>

namespace {

constexpr double kMRhoNom = 0.77526;
constexpr double kMSigmaNom = 0.500;

struct FourierMoment {
    double sw = 0.0;
    double c[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    double s[5] = {0.0, 0.0, 0.0, 0.0, 0.0};

    void Fill(double x, double w) {
        sw += w;
        for (int n = 1; n <= 4; ++n) {
            c[n] += w * std::cos(n * x);
            s[n] += w * std::sin(n * x);
        }
    }

    double CosCoeff(int n) const { return (sw > 0.0) ? 2.0 * c[n] / sw : 0.0; }
    double SinCoeff(int n) const { return (sw > 0.0) ? 2.0 * s[n] / sw : 0.0; }
};

double PairPhi(TLorentzVector* pi_plus, TLorentzVector* pi_minus,
               const TVector3& x_axis, const TVector3& y_axis) {
    TLorentzVector pair = *pi_plus + *pi_minus;
    TLorentzVector pi_plus_rest = *pi_plus;
    pi_plus_rest.Boost(-pair.BoostVector());
    TVector3 p3 = pi_plus_rest.Vect();
    return std::atan2(p3.Dot(y_axis), p3.Dot(x_axis));
}

void PrintMoment(const char* label, const FourierMoment& m) {
    std::cout << "\n--- " << label << "  sumw=" << m.sw << std::endl;
    for (int n = 1; n <= 4; ++n) {
        std::cout << "    2<cos" << n << "x>=" << m.CosCoeff(n)
                  << "    2<sin" << n << "x>=" << m.SinCoeff(n) << std::endl;
    }
}

void WriteMomentParams(const char* prefix, const FourierMoment& m) {
    TParameter<double>((std::string(prefix) + "_sumw").c_str(), m.sw).Write();
    for (int n = 1; n <= 4; ++n) {
        TParameter<double>((std::string(prefix) + "_2cos" + std::to_string(n)).c_str(),
                           m.CosCoeff(n)).Write();
        TParameter<double>((std::string(prefix) + "_2sin" + std::to_string(n)).c_str(),
                           m.SinCoeff(n)).Write();
    }
}

}  // namespace

int main(int argc, char** argv) {
    const char* input_file  = (argc > 1) ? argv[1] : "BoseEvents_AuAu200_XnXn_y1_dp001_nonbose_100M.root";
    const char* output_file = (argc > 2) ? argv[2] : "analyze_nonbose_output.root";
    Long64_t max_events     = (argc > 3) ? std::atoll(argv[3]) : 1000000;

    TFile* fin = TFile::Open(input_file, "READ");
    if (!fin || fin->IsZombie()) { std::cerr << "ERROR open input" << std::endl; return 1; }
    TTree* t = dynamic_cast<TTree*>(fin->Get("Events"));
    if (!t) { std::cerr << "ERROR no Events" << std::endl; return 1; }

    Long64_t nentries = t->GetEntries();
    Long64_t nprocess = std::min(nentries, max_events);
    std::cout << "Processing " << nprocess << " / " << nentries << " events" << std::endl;
    std::cout << "Input: " << input_file << std::endl;

    // Branch addresses
    Double_t event_weight, M_stored, y_stored, pt_stored, m_rho, m_sigma;
    TLorentzVector *pi1 = nullptr, *pi2 = nullptr, *pi3 = nullptr, *pi4 = nullptr;
    TLorentzVector *rho_momentum = nullptr, *sigma_momentum = nullptr;
    Int_t pi_parent[4], pi_charge[4];

    t->SetBranchAddress("event_weight", &event_weight);
    t->SetBranchAddress("M_rhoprime",   &M_stored);
    t->SetBranchAddress("y_rhoprime",   &y_stored);
    t->SetBranchAddress("pt_rhoprime",  &pt_stored);
    t->SetBranchAddress("m_rho",        &m_rho);
    t->SetBranchAddress("m_sigma",      &m_sigma);
    t->SetBranchAddress("pi1",          &pi1);
    t->SetBranchAddress("pi2",          &pi2);
    t->SetBranchAddress("pi3",          &pi3);
    t->SetBranchAddress("pi4",          &pi4);
    t->SetBranchAddress("rho_momentum", &rho_momentum);
    t->SetBranchAddress("sigma_momentum", &sigma_momentum);
    t->SetBranchAddress("pi_parent",    pi_parent);
    t->SetBranchAddress("pi_charge",    pi_charge);

    // Histograms for rho pair characteristics
    TH1D* h_M_rho_pair = new TH1D("h_M_rho_pair", ";m_{#pi+#pi- from #rho} [GeV];Events", 100, 0.2, 1.5);
    TH1D* h_M_sig_pair = new TH1D("h_M_sig_pair", ";m_{#pi+#pi- from #sigma} [GeV];Events", 100, 0.2, 1.5);
    TH1D* h_pt_rho     = new TH1D("h_pt_rho",     ";p_{T}(#rho pair) [GeV];Events",        20, 0.0, 0.2);

    // rho' reconstruction from 4 pions
    TH1D* h_M_reco   = new TH1D("h_M_reco",   ";M_{4#pi} [GeV];Events",           200, 0.5, 5.0);
    TH1D* h_y_reco   = new TH1D("h_y_reco",   ";y_{4#pi};Events",                100, -5.0, 5.0);
    TH1D* h_pt_reco  = new TH1D("h_pt_reco",  ";p_{T,4#pi} [GeV];Events",        20, 0.0, 0.2);
    TH1D* h_M_store  = new TH1D("h_M_store",  ";M_{stored} [GeV];Events",        200, 0.5, 5.0);
    TH1D* h_y_store  = new TH1D("h_y_store",  ";y_{stored};Events",             100, -5.0, 5.0);
    TH1D* h_pt_store = new TH1D("h_pt_store", ";p_{T,stored} [GeV];Events",      20, 0.0, 0.2);
    TH1D* h_M_diff   = new TH1D("h_M_diff",   ";M_{stored} - M_{4#pi} [GeV];Events", 200, -0.01, 0.01);

    TH1D* h_cos_theta_rho = new TH1D("h_cos_theta_rho",
        ";cos#theta_{#pi+ in #rho rest};Events", 40, -1.0, 1.0);
    TH1D* h_cos_theta_rho_win = new TH1D("h_cos_theta_rho_win",
        ";cos#theta_{#pi+ in #rho rest} (0.6<M<0.9);Events", 40, -1.0, 1.0);
    TH1D* h_cos_theta_sigma = new TH1D("h_cos_theta_sigma",
        ";cos#theta_{#pi+ in #sigma rest};Events", 40, -1.0, 1.0);

    TH1D* h_phi_rho = new TH1D("h_phi_rho",
        ";#phi_{#pi+ in #rho rest} [rad];Events", 40, -M_PI, M_PI);
    TH1D* h_phi_sigma = new TH1D("h_phi_sigma",
        ";#phi_{#pi+ in #sigma rest} [rad];Events", 40, -M_PI, M_PI);
    TH1D* h_dphi_rho_sigma = new TH1D("h_dphi_rho_sigma",
        ";#Delta#phi_{#rho-#sigma} [rad];Events", 40, -M_PI, M_PI);
    TH2D* h_phi_rho_sigma = new TH2D("h_phi_rho_sigma",
        ";#phi_{#rho} [rad];#phi_{#sigma} [rad]", 40, -M_PI, M_PI, 40, -M_PI, M_PI);

    TH1D* h_phi_truth_rho = new TH1D("h_phi_truth_rho",
        ";#phi_{#pi+}^{truth #rho} [rad];Events", 80, -M_PI, M_PI);
    TH1D* h_phi_truth_sigma = new TH1D("h_phi_truth_sigma",
        ";#phi_{#pi+}^{truth #sigma} [rad];Events", 80, -M_PI, M_PI);
    TH1D* h_phi_all_pairs = new TH1D("h_phi_all_pairs",
        ";#phi_{#pi+}^{all #pi+#pi- pairs} [rad];Events", 80, -M_PI, M_PI);
    TH1D* h_phi_best_rho = new TH1D("h_phi_best_rho",
        ";#phi_{#pi+}^{best #rho pair} [rad];Events", 80, -M_PI, M_PI);
    TH1D* h_dphi_truth_rho_sigma = new TH1D("h_dphi_truth_rho_sigma",
        ";#Delta#phi^{truth}_{#rho-#sigma} [rad];Events", 80, -M_PI, M_PI);
    TH1D* h_dphi_all_matchings = new TH1D("h_dphi_all_matchings",
        ";#Delta#phi^{all OS matchings} [rad];Events", 80, -M_PI, M_PI);
    TH1D* h_dphi_best_rho_sigma = new TH1D("h_dphi_best_rho_sigma",
        ";#Delta#phi^{best #rho/#sigma} [rad];Events", 80, -M_PI, M_PI);

    // Spin alignment: 2<cos(phi)> and 2<cos(2phi)> vs rho' pT
    // Only for rho pairs with mass in [0.6, 0.9]
    TProfile* h_cos2phi_vs_pt = new TProfile("h_cos2phi_vs_pt",
        ";p_{T}(#rho') [GeV];2#upoint<cos2#phi>_{#rho}", 20, 0.0, 0.2);
    TProfile* h_cos1phi_vs_pt = new TProfile("h_cos1phi_vs_pt",
        ";p_{T}(#rho') [GeV];2#upoint<cos#phi>_{#rho}",  20, 0.0, 0.2);

    // Also vs m_rho (rho pair mass)
    TProfile* h_cos2phi_vs_m = new TProfile("h_cos2phi_vs_m",
        ";m_{#pi+#pi-} [GeV];2#upoint<cos2#phi>_{#rho}", 40, 0.5, 1.2);
    TProfile* h_cos1phi_vs_m = new TProfile("h_cos1phi_vs_m",
        ";m_{#pi+#pi-} [GeV];2#upoint<cos#phi>_{#rho}",  40, 0.5, 1.2);

    TProfile* h_sigma_cos2phi_vs_pt = new TProfile("h_sigma_cos2phi_vs_pt",
        ";p_{T}(#rho') [GeV];2#upoint<cos2#phi>_{#sigma}", 20, 0.0, 0.2);
    TProfile* h_sigma_cos1phi_vs_pt = new TProfile("h_sigma_cos1phi_vs_pt",
        ";p_{T}(#rho') [GeV];2#upoint<cos#phi>_{#sigma}",  20, 0.0, 0.2);

    FourierMoment m_phi_truth_rho, m_phi_truth_sigma, m_phi_all_pairs, m_phi_best_rho;
    FourierMoment m_dphi_truth_rho_sigma, m_dphi_all_matchings, m_dphi_best_rho_sigma;

    // Weight range scan
    double wmax = 0.0;
    Long64_t n_scan = std::min(nprocess, 100000LL);
    for (Long64_t i = 0; i < n_scan; ++i) {
        t->GetEntry(i);
        if (std::isfinite(event_weight) && event_weight > wmax) wmax = event_weight;
    }
    TH1D* h_weight = new TH1D("h_weight", ";event weight;Events", 200, 0.0, wmax * 1.1);

    // Main loop
    double sum_w = 0.0, sum_w2 = 0.0;
    Long64_t step = nprocess / 10 + 1;

    for (Long64_t i = 0; i < nprocess; ++i) {
        t->GetEntry(i);
        if (i % step == 0)
            std::cout << "  " << i << "/" << nprocess << " (" << 100.0*i/nprocess << "%)" << std::endl;

        double w = event_weight;
        if (!std::isfinite(w) || w <= 0.0) continue;
        sum_w += w; sum_w2 += w * w;
        h_weight->Fill(w);

        // Find rho pair: pions with pi_parent == 1
        TLorentzVector* pions[4] = {pi1, pi2, pi3, pi4};
        TLorentzVector *pi_plus_rho = nullptr, *pi_minus_rho = nullptr;
        int idx_plus[2] = {-1, -1};
        int idx_minus[2] = {-1, -1};
        int nplus = 0, nminus = 0;

        for (int j = 0; j < 4; ++j) {
            if (pi_charge[j] > 0 && nplus < 2) idx_plus[nplus++] = j;
            if (pi_charge[j] < 0 && nminus < 2) idx_minus[nminus++] = j;
            if (pi_parent[j] == 1) {
                if (pi_charge[j] > 0)  pi_plus_rho  = pions[j];
                else                   pi_minus_rho = pions[j];
            }
        }
        if (!pi_plus_rho || !pi_minus_rho) continue;

        TLorentzVector P_rho = *pi_plus_rho + *pi_minus_rho;
        double M_rho_pair = P_rho.M();
        h_M_rho_pair->Fill(M_rho_pair, w);
        h_pt_rho->Fill(P_rho.Pt(), w);

        // Sigma pair: pi_parent == 2
        TLorentzVector *pi_plus_sig = nullptr, *pi_minus_sig = nullptr;
        for (int j = 0; j < 4; ++j) {
            if (pi_parent[j] == 2) {
                if (pi_charge[j] > 0)  pi_plus_sig  = pions[j];
                else                   pi_minus_sig = pions[j];
            }
        }
        if (pi_plus_sig && pi_minus_sig) {
            TLorentzVector P_sig = *pi_plus_sig + *pi_minus_sig;
            h_M_sig_pair->Fill(P_sig.M(), w);
        }

        // cos(theta) of pi+ in rho rest frame (z = lab-z)
        {
            TVector3 beta_rho = P_rho.BoostVector();
            TLorentzVector pi_plus_rest = *pi_plus_rho;
            pi_plus_rest.Boost(-beta_rho);
            double cos_theta = pi_plus_rest.Vect().Unit().Z();
            h_cos_theta_rho->Fill(cos_theta, w);
            if (M_rho_pair >= 0.6 && M_rho_pair <= 0.9)
                h_cos_theta_rho_win->Fill(cos_theta, w);
        }

        // rho' 4-momentum
        TLorentzVector p_rho_prime = *pi1 + *pi2 + *pi3 + *pi4;
        double M_reco  = p_rho_prime.M();
        double y_reco  = p_rho_prime.Rapidity();
        double pt_prime = p_rho_prime.Pt();

        h_M_reco->Fill(M_reco, w);
        h_y_reco->Fill(y_reco, w);
        h_pt_reco->Fill(pt_prime, w);
        h_M_store->Fill(M_stored, w);
        h_y_store->Fill(y_stored, w);
        h_pt_store->Fill(pt_stored, w);
        h_M_diff->Fill(M_stored - M_reco, w);

        // rho' pT direction in lab
        TVector3 pt_dir = p_rho_prime.Vect();
        pt_dir.SetZ(0);
        if (pt_dir.Mag() < 1e-12) continue;
        pt_dir = pt_dir.Unit();

        TVector3 z_axis(0, 0, 1);
        TVector3 x_axis = pt_dir;
        TVector3 y_axis = z_axis.Cross(x_axis).Unit();

        const double phi_truth_rho = PairPhi(pi_plus_rho, pi_minus_rho, x_axis, y_axis);
        h_phi_truth_rho->Fill(phi_truth_rho, w);
        m_phi_truth_rho.Fill(phi_truth_rho, w);

        if (pi_plus_sig && pi_minus_sig) {
            const double phi_truth_sigma = PairPhi(pi_plus_sig, pi_minus_sig, x_axis, y_axis);
            const double dphi_truth = TVector2::Phi_mpi_pi(phi_truth_rho - phi_truth_sigma);
            h_phi_truth_sigma->Fill(phi_truth_sigma, w);
            h_dphi_truth_rho_sigma->Fill(dphi_truth, w);
            m_phi_truth_sigma.Fill(phi_truth_sigma, w);
            m_dphi_truth_rho_sigma.Fill(dphi_truth, w);
        }

        if (nplus == 2 && nminus == 2) {
            double phi_pair[2][2];
            double mass_pair[2][2];
            for (int ip = 0; ip < 2; ++ip) {
                for (int im = 0; im < 2; ++im) {
                    phi_pair[ip][im] = PairPhi(pions[idx_plus[ip]], pions[idx_minus[im]], x_axis, y_axis);
                    mass_pair[ip][im] = (*pions[idx_plus[ip]] + *pions[idx_minus[im]]).M();
                    h_phi_all_pairs->Fill(phi_pair[ip][im], w);
                    m_phi_all_pairs.Fill(phi_pair[ip][im], w);
                }
            }

            const double dphi_match_a = TVector2::Phi_mpi_pi(phi_pair[0][0] - phi_pair[1][1]);
            const double dphi_match_b = TVector2::Phi_mpi_pi(phi_pair[0][1] - phi_pair[1][0]);
            h_dphi_all_matchings->Fill(dphi_match_a, w);
            h_dphi_all_matchings->Fill(dphi_match_b, w);
            m_dphi_all_matchings.Fill(dphi_match_a, w);
            m_dphi_all_matchings.Fill(dphi_match_b, w);

            const double score0 = std::abs(mass_pair[0][0] - kMRhoNom) + std::abs(mass_pair[1][1] - kMSigmaNom);
            const double score1 = std::abs(mass_pair[1][1] - kMRhoNom) + std::abs(mass_pair[0][0] - kMSigmaNom);
            const double score2 = std::abs(mass_pair[0][1] - kMRhoNom) + std::abs(mass_pair[1][0] - kMSigmaNom);
            const double score3 = std::abs(mass_pair[1][0] - kMRhoNom) + std::abs(mass_pair[0][1] - kMSigmaNom);

            int best = 0;
            double best_score = score0;
            if (score1 < best_score) { best = 1; best_score = score1; }
            if (score2 < best_score) { best = 2; best_score = score2; }
            if (score3 < best_score) { best = 3; best_score = score3; }

            double phi_best_rho = phi_pair[0][0];
            double phi_best_sigma = phi_pair[1][1];
            if (best == 1) { phi_best_rho = phi_pair[1][1]; phi_best_sigma = phi_pair[0][0]; }
            if (best == 2) { phi_best_rho = phi_pair[0][1]; phi_best_sigma = phi_pair[1][0]; }
            if (best == 3) { phi_best_rho = phi_pair[1][0]; phi_best_sigma = phi_pair[0][1]; }

            const double dphi_best = TVector2::Phi_mpi_pi(phi_best_rho - phi_best_sigma);
            h_phi_best_rho->Fill(phi_best_rho, w);
            h_dphi_best_rho_sigma->Fill(dphi_best, w);
            m_phi_best_rho.Fill(phi_best_rho, w);
            m_dphi_best_rho_sigma.Fill(dphi_best, w);
        }

        // Check rho mass window for spin analysis
        if (M_rho_pair < 0.6 || M_rho_pair > 0.9) continue;

        // Boost pi+ from rho to rho pair rest frame
        TVector3 beta_rho = P_rho.BoostVector();
        TLorentzVector pi_plus_rest = *pi_plus_rho;
        pi_plus_rest.Boost(-beta_rho);

        // Axes: pure spatial, no boost needed
        TVector3 p3 = pi_plus_rest.Vect();
        double phi = std::atan2(p3.Dot(y_axis), p3.Dot(x_axis));
        double cos2phi = std::cos(2.0 * phi);
        double cos1phi = std::cos(phi);

        h_cos2phi_vs_pt->Fill(pt_prime, 2.0 * cos2phi, w);
        h_cos1phi_vs_pt->Fill(pt_prime, 2.0 * cos1phi, w);
        h_cos2phi_vs_m->Fill(M_rho_pair, 2.0 * cos2phi, w);
        h_cos1phi_vs_m->Fill(M_rho_pair, 2.0 * cos1phi, w);

        h_phi_rho->Fill(phi, w);

        if (pi_plus_sig && pi_minus_sig) {
            TLorentzVector P_sig = *pi_plus_sig + *pi_minus_sig;
            TVector3 beta_sig = P_sig.BoostVector();
            TLorentzVector pi_plus_sig_rest = *pi_plus_sig;
            pi_plus_sig_rest.Boost(-beta_sig);

            TVector3 p3_sig = pi_plus_sig_rest.Vect();
            if (p3_sig.Mag() > 0.0) {
                const double cos_theta_sig = p3_sig.Unit().Z();
                const double phi_sig = std::atan2(p3_sig.Dot(y_axis), p3_sig.Dot(x_axis));
                const double dphi = TVector2::Phi_mpi_pi(phi - phi_sig);

                h_cos_theta_sigma->Fill(cos_theta_sig, w);
                h_phi_sigma->Fill(phi_sig, w);
                h_dphi_rho_sigma->Fill(dphi, w);
                h_phi_rho_sigma->Fill(phi, phi_sig, w);
                h_sigma_cos2phi_vs_pt->Fill(pt_prime, 2.0 * std::cos(2.0 * phi_sig), w);
                h_sigma_cos1phi_vs_pt->Fill(pt_prime, 2.0 * std::cos(phi_sig), w);
            }
        }
    }
    std::cout << "Loop done." << std::endl;

    double Neff = (sum_w2 > 0) ? sum_w*sum_w/sum_w2 : 0;
    std::cout << "Sum weights: " << sum_w << "  Neff: " << Neff << std::endl;

    // Null pointers
    pi1 = pi2 = pi3 = pi4 = rho_momentum = sigma_momentum = nullptr;

    // Write output
    std::cout << "Writing " << output_file << "..." << std::endl;
    TFile* fout = TFile::Open(output_file, "RECREATE");
    if (!fout || fout->IsZombie()) { std::cerr << "ERROR creating output" << std::endl; return 1; }

    h_M_rho_pair->Write();
    h_M_sig_pair->Write();
    h_cos_theta_rho->Write();  h_cos_theta_rho_win->Write();
    h_pt_rho->Write();
    h_M_reco->Write();  h_y_reco->Write();  h_pt_reco->Write();
    h_M_store->Write(); h_y_store->Write(); h_pt_store->Write(); h_M_diff->Write();
    h_cos_theta_sigma->Write();
    h_phi_rho->Write(); h_phi_sigma->Write(); h_dphi_rho_sigma->Write(); h_phi_rho_sigma->Write();
    h_phi_truth_rho->Write(); h_phi_truth_sigma->Write();
    h_phi_all_pairs->Write(); h_phi_best_rho->Write();
    h_dphi_truth_rho_sigma->Write(); h_dphi_all_matchings->Write(); h_dphi_best_rho_sigma->Write();
    h_cos2phi_vs_pt->Write();  h_cos1phi_vs_pt->Write();
    h_cos2phi_vs_m->Write();   h_cos1phi_vs_m->Write();
    h_sigma_cos2phi_vs_pt->Write(); h_sigma_cos1phi_vs_pt->Write();
    h_weight->Write();
    TParameter<double>("sum_event_weight", sum_w).Write();
    TParameter<double>("Neff", Neff).Write();
    WriteMomentParams("mod_phi_truth_rho", m_phi_truth_rho);
    WriteMomentParams("mod_phi_truth_sigma", m_phi_truth_sigma);
    WriteMomentParams("mod_phi_all_pairs", m_phi_all_pairs);
    WriteMomentParams("mod_phi_best_rho", m_phi_best_rho);
    WriteMomentParams("mod_dphi_truth_rho_sigma", m_dphi_truth_rho_sigma);
    WriteMomentParams("mod_dphi_all_matchings", m_dphi_all_matchings);
    WriteMomentParams("mod_dphi_best_rho_sigma", m_dphi_best_rho_sigma);
    fout->Close();
    std::cout << "Output saved: " << output_file << std::endl;

    // Print summary
    std::cout << "\n--- Rho' reco M [GeV]  mean=" << h_M_reco->GetMean()
              << "  rms=" << h_M_reco->GetRMS() << std::endl;
    std::cout << "--- Rho' reco y          mean=" << h_y_reco->GetMean()
              << "  rms=" << h_y_reco->GetRMS() << std::endl;
    std::cout << "--- Rho' reco pT [GeV]   mean=" << h_pt_reco->GetMean()
              << "  rms=" << h_pt_reco->GetRMS() << std::endl;
    std::cout << "--- M_diff [GeV]         mean=" << h_M_diff->GetMean()
              << "  rms=" << h_M_diff->GetRMS() << std::endl;
    std::cout << "\n--- Rho pair mass   mean=" << h_M_rho_pair->GetMean()
              << "  rms=" << h_M_rho_pair->GetRMS() << std::endl;
    std::cout << "--- Sigma pair mass  mean=" << h_M_sig_pair->GetMean()
              << "  rms=" << h_M_sig_pair->GetRMS() << std::endl;
    std::cout << "--- cos_theta_rho    mean=" << h_cos_theta_rho->GetMean()
              << "  rms=" << h_cos_theta_rho->GetRMS() << std::endl;
    std::cout << "--- cos_theta_rho(win) mean=" << h_cos_theta_rho_win->GetMean()
              << "  rms=" << h_cos_theta_rho_win->GetRMS() << std::endl;
    std::cout << "--- cos_theta_sigma mean=" << h_cos_theta_sigma->GetMean()
              << "  rms=" << h_cos_theta_sigma->GetRMS() << std::endl;
    std::cout << "--- dphi(rho-sigma) mean=" << h_dphi_rho_sigma->GetMean()
              << "  rms=" << h_dphi_rho_sigma->GetRMS() << std::endl;

    std::cout << "\n=== Fourier modulation moments: f(x) ~ 1 + a_n cos(nx) + b_n sin(nx) ===" << std::endl;
    PrintMoment("phi truth rho", m_phi_truth_rho);
    PrintMoment("phi truth sigma", m_phi_truth_sigma);
    PrintMoment("phi all pi+pi- pairs, no truth", m_phi_all_pairs);
    PrintMoment("phi best rho-like pair, no truth", m_phi_best_rho);
    PrintMoment("Delta phi truth rho-sigma", m_dphi_truth_rho_sigma);
    PrintMoment("Delta phi all opposite-sign matchings, no truth", m_dphi_all_matchings);
    PrintMoment("Delta phi best rho/sigma assignment, no truth", m_dphi_best_rho_sigma);

    std::cout << "--- rho_cos2phi/cos1phi vs pT(ρ'):" << std::endl;
    for (int b = 1; b <= h_cos2phi_vs_pt->GetNbinsX(); ++b)
        std::cout << "    pT[" << h_cos2phi_vs_pt->GetBinLowEdge(b) << "-"
                  << h_cos2phi_vs_pt->GetBinLowEdge(b)+h_cos2phi_vs_pt->GetBinWidth(b)
                  << "]  2<cos2phi>=" << h_cos2phi_vs_pt->GetBinContent(b)
                  << " ± " << h_cos2phi_vs_pt->GetBinError(b)
                  << "    2<cosphi>=" << h_cos1phi_vs_pt->GetBinContent(b)
                  << " ± " << h_cos1phi_vs_pt->GetBinError(b) << std::endl;

    std::cout << "\nDone." << std::endl;
    return 0;
}

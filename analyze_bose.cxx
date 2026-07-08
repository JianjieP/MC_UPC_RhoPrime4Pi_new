// analyze_bose.cxx — reconstruct rho' from 4-pion, TLorentzVector pointers fixed
// Build: cmake + make
// Added: cos(2phi) spin alignment analysis vs pT

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
#include <vector>

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
    const char* input_file  = (argc > 1) ? argv[1] : "BoseEvents_AuAu200_XnXn_y1_dp001_bose_100M.root";
    const char* output_file = (argc > 2) ? argv[2] : "analyze_bose_output.root";
    Long64_t max_events     = (argc > 3) ? std::atoll(argv[3]) : 1000000;

    TFile* fin = TFile::Open(input_file, "READ");
    if (!fin || fin->IsZombie()) { std::cerr << "ERROR open input" << std::endl; return 1; }
    TTree* t = dynamic_cast<TTree*>(fin->Get("Events"));
    if (!t) { std::cerr << "ERROR no Events" << std::endl; return 1; }

    Long64_t nentries = t->GetEntries();
    Long64_t nprocess = std::min(nentries, max_events);
    std::cout << "Processing " << nprocess << " / " << nentries << " events" << std::endl;

    // Branch addresses — TLorentzVector objects must use pointer-to-pointer
    Double_t event_weight, M_stored, y_stored, pt_stored, m_rho, m_sigma;
    TLorentzVector *pi1 = nullptr, *pi2 = nullptr, *pi3 = nullptr, *pi4 = nullptr;

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

    Int_t pi_charge[4];
    t->SetBranchAddress("pi_charge",  pi_charge);

    // Histograms
    TH1D* h_M_reco   = new TH1D("h_M_reco",   ";M_{4#pi} [GeV];Events",           200, 0.5, 5.0);
    TH1D* h_y_reco   = new TH1D("h_y_reco",   ";y_{4#pi};Events",                100, -5.0, 5.0);
    TH1D* h_pt_reco  = new TH1D("h_pt_reco",  ";p_{T,4#pi} [GeV];Events",        20, 0.0, 0.2);
    TH1D* h_M_store  = new TH1D("h_M_store",  ";M_{stored} [GeV];Events",        200, 0.5, 5.0);
    TH1D* h_y_store  = new TH1D("h_y_store",  ";y_{stored};Events",             100, -5.0, 5.0);
    TH1D* h_pt_store = new TH1D("h_pt_store", ";p_{T,stored} [GeV];Events",     20, 0.0, 0.2);
    TH1D* h_M_diff   = new TH1D("h_M_diff",   ";M_{stored} - M_{4#pi} [GeV];Events", 200, -0.01, 0.01);
    TH2D* h_M_y_reco = new TH2D("h_M_y_reco", ";M_{4#pi} [GeV];y_{4#pi}", 100, 0.5, 5.0, 60, -3.0, 3.0);
    TH2D* h_M_pt_reco= new TH2D("h_M_pt_reco",";M_{4#pi} [GeV];p_{T} [GeV]", 100, 0.5, 5.0, 20, 0.0, 0.2);
    TH1D* h_m_rho    = new TH1D("h_m_rho",    ";m_{#rho^{0}} [GeV];Events",      100, 0.2, 3.5);
    TH1D* h_m_sigma  = new TH1D("h_m_sigma",  ";m_{#sigma} [GeV];Events",       100, 0.2, 3.5);
    TH2D* h_mrho_msig= new TH2D("h_mrho_msig",";m_{#rho^{0}} [GeV];m_{#sigma} [GeV]", 60, 0.2, 2.5, 60, 0.2, 2.5);
    TH1D* h_M_os_pairs = new TH1D("h_M_os_pairs",
        ";m_{#pi+#pi-}^{all OS pairs} [GeV];Pairs", 100, 0.2, 1.5);
    TH1D* h_M_best_rho_pair = new TH1D("h_M_best_rho_pair",
        ";m_{#pi+#pi-}^{best #rho pair} [GeV];Events", 100, 0.2, 1.5);
    TH1D* h_M_best_sigma_pair = new TH1D("h_M_best_sigma_pair",
        ";m_{#pi+#pi-}^{best #sigma pair} [GeV];Events", 100, 0.2, 1.5);

    // cos(2phi) and cos(1phi) spin alignment: TProfile
    TProfile* h_cos2phi_vs_pt = new TProfile("h_cos2phi_vs_pt",
        ";p_{T}(#rho') [GeV];2#upoint<cos2#phi>", 20, 0.0, 0.2);
    TProfile* h_cos1phi_vs_pt = new TProfile("h_cos1phi_vs_pt",
        ";p_{T}(#rho') [GeV];2#upoint<cos#phi>", 20, 0.0, 0.2);
    TProfile* h_npair_vs_pt = new TProfile("h_npair_vs_pt",
        ";p_{T}(#rho') [GeV];<npairs> per event", 20, 0.0, 0.2);

    // cos(theta) for accepted pi+pi- pairs (0.6<M<0.9)
    TH1D* h_cos_theta_bose = new TH1D("h_cos_theta_bose",
        ";cos#theta_{#pi+ in pair rest} (0.6<M<0.9);Events", 40, -1.0, 1.0);
    TH1D* h_phi_bose_pairs = new TH1D("h_phi_bose_pairs",
        ";#phi_{#pi+ in pair rest} (0.6<M<0.9) [rad];Pairs", 80, -M_PI, M_PI);
    TH1D* h_phi_all_pairs = new TH1D("h_phi_all_pairs",
        ";#phi_{#pi+}^{all #pi+#pi- pairs} [rad];Pairs", 80, -M_PI, M_PI);
    TH1D* h_phi_best_rho = new TH1D("h_phi_best_rho",
        ";#phi_{#pi+}^{best #rho pair} [rad];Events", 80, -M_PI, M_PI);
    TH1D* h_dphi_all_matchings = new TH1D("h_dphi_all_matchings",
        ";#Delta#phi^{all OS matchings} [rad];Matchings", 80, -M_PI, M_PI);
    TH1D* h_dphi_best_rho_sigma = new TH1D("h_dphi_best_rho_sigma",
        ";#Delta#phi^{best #rho/#sigma} [rad];Events", 80, -M_PI, M_PI);

    FourierMoment m_phi_all_pairs, m_phi_best_rho;
    FourierMoment m_dphi_all_matchings, m_dphi_best_rho_sigma;

    // Weight range
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

        TLorentzVector p4pi = *pi1 + *pi2 + *pi3 + *pi4;
        double M_reco  = p4pi.M();
        double y_reco  = p4pi.Rapidity();
        double pt_reco = p4pi.Pt();

        h_weight->Fill(w);
        h_M_reco->Fill(M_reco, w);   h_y_reco->Fill(y_reco, w);   h_pt_reco->Fill(pt_reco, w);
        h_M_store->Fill(M_stored, w); h_y_store->Fill(y_stored, w); h_pt_store->Fill(pt_stored, w);
        h_M_diff->Fill(M_stored - M_reco, w);
        h_M_y_reco->Fill(M_reco, y_reco, w);
        h_M_pt_reco->Fill(M_reco, pt_reco, w);
        h_m_rho->Fill(m_rho, w);
        h_m_sigma->Fill(m_sigma, w);
        h_mrho_msig->Fill(m_rho, m_sigma, w);

        // — cos(2phi) spin alignment analysis —
        // 4 pairings of (pi+, pi-) — pi_charge = {1,-1,1,-1} for pi1..pi4
        TLorentzVector* pions[4] = {pi1, pi2, pi3, pi4};
        std::vector<int> idx_plus, idx_minus;
        for (int j = 0; j < 4; ++j) {
            if (pi_charge[j] > 0) idx_plus.push_back(j);
            else                 idx_minus.push_back(j);
        }
        if (idx_plus.size() != 2 || idx_minus.size() != 2) continue;

        TLorentzVector p_rho_prime = *pi1 + *pi2 + *pi3 + *pi4;
        double pt_prime = p_rho_prime.Pt();
        // rho' transverse direction in lab
        TVector3 pt_dir_lab = p_rho_prime.Vect();
        pt_dir_lab.SetZ(0);
        if (pt_dir_lab.Mag() < 1e-12) continue;
        pt_dir_lab = pt_dir_lab.Unit();

        TVector3 z_axis(0, 0, 1);
        TVector3 x_axis = pt_dir_lab;
        TVector3 y_axis = z_axis.Cross(x_axis).Unit();

        double phi_pair[2][2];
        double mass_pair[2][2];
        for (int ip = 0; ip < 2; ++ip) {
            for (int im = 0; im < 2; ++im) {
                phi_pair[ip][im] = PairPhi(pions[idx_plus[ip]], pions[idx_minus[im]], x_axis, y_axis);
                mass_pair[ip][im] = (*pions[idx_plus[ip]] + *pions[idx_minus[im]]).M();
                h_M_os_pairs->Fill(mass_pair[ip][im], w);
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
        double mass_best_rho = mass_pair[0][0];
        double mass_best_sigma = mass_pair[1][1];
        if (best == 1) {
            phi_best_rho = phi_pair[1][1]; phi_best_sigma = phi_pair[0][0];
            mass_best_rho = mass_pair[1][1]; mass_best_sigma = mass_pair[0][0];
        }
        if (best == 2) {
            phi_best_rho = phi_pair[0][1]; phi_best_sigma = phi_pair[1][0];
            mass_best_rho = mass_pair[0][1]; mass_best_sigma = mass_pair[1][0];
        }
        if (best == 3) {
            phi_best_rho = phi_pair[1][0]; phi_best_sigma = phi_pair[0][1];
            mass_best_rho = mass_pair[1][0]; mass_best_sigma = mass_pair[0][1];
        }

        const double dphi_best = TVector2::Phi_mpi_pi(phi_best_rho - phi_best_sigma);
        h_M_best_rho_pair->Fill(mass_best_rho, w);
        h_M_best_sigma_pair->Fill(mass_best_sigma, w);
        h_phi_best_rho->Fill(phi_best_rho, w);
        h_dphi_best_rho_sigma->Fill(dphi_best, w);
        m_phi_best_rho.Fill(phi_best_rho, w);
        m_dphi_best_rho_sigma.Fill(dphi_best, w);

        double sum_cos2phi = 0.0;
        double sum_cos1phi = 0.0;
        int npair_acc = 0;

        for (int ip : idx_plus) {
            for (int im : idx_minus) {
                TLorentzVector P_pair = *pions[ip] + *pions[im];
                double M_pair = P_pair.M();
                if (M_pair < 0.6 || M_pair > 0.9) continue;

                // Boost to pair rest frame
                TVector3 beta_pair = P_pair.BoostVector();
                TLorentzVector pi_plus_rest = *pions[ip];
                pi_plus_rest.Boost(-beta_pair);

                // Axes defined in lab frame — pure spatial directions, no boost needed
                // phi of pi+ in this frame
                TVector3 p3 = pi_plus_rest.Vect();
                double phi = std::atan2(p3.Dot(y_axis), p3.Dot(x_axis));
                sum_cos2phi += std::cos(2.0 * phi);
                sum_cos1phi += std::cos(phi);
                // cos(theta): z-axis = lab-z
                double cos_theta = p3.Unit().Z();
                h_cos_theta_bose->Fill(cos_theta, w);
                h_phi_bose_pairs->Fill(phi, w);
                npair_acc++;
            }
        }

        if (npair_acc > 0) {
            double cos2phi_avg = sum_cos2phi / npair_acc;
            double cos1phi_avg = sum_cos1phi / npair_acc;
            h_cos2phi_vs_pt->Fill(pt_prime, 2.0 * cos2phi_avg, w);
            h_cos1phi_vs_pt->Fill(pt_prime, 2.0 * cos1phi_avg, w);
            h_npair_vs_pt->Fill(pt_prime, npair_acc, w);
        }
    }
    std::cout << "Loop done." << std::endl;

    double Neff = (sum_w2 > 0) ? sum_w*sum_w/sum_w2 : 0;
    std::cout << "Sum weights: " << sum_w << "  Neff: " << Neff << std::endl;

    // Null TLorentzVector pointers so fin->Close() doesn't try to delete them
    pi1 = pi2 = pi3 = pi4 = nullptr;

    // Write output (don't call fin->Close, let it happen on exit)
    std::cout << "Writing " << output_file << "..." << std::endl;
    TFile* fout = TFile::Open(output_file, "RECREATE");
    if (!fout || fout->IsZombie()) { std::cerr << "ERROR creating output" << std::endl; return 1; }

    h_M_reco->Write();   h_y_reco->Write();   h_pt_reco->Write();
    h_M_store->Write();  h_y_store->Write();  h_pt_store->Write();
    h_M_diff->Write();   h_M_y_reco->Write(); h_M_pt_reco->Write();
    h_m_rho->Write();    h_m_sigma->Write();  h_mrho_msig->Write();
    h_M_os_pairs->Write(); h_M_best_rho_pair->Write(); h_M_best_sigma_pair->Write();
    h_weight->Write();
    h_cos2phi_vs_pt->Write();  h_cos1phi_vs_pt->Write();  h_npair_vs_pt->Write();
    h_cos_theta_bose->Write(); h_phi_bose_pairs->Write();
    h_phi_all_pairs->Write(); h_phi_best_rho->Write();
    h_dphi_all_matchings->Write(); h_dphi_best_rho_sigma->Write();
    TParameter<double>("sum_event_weight", sum_w).Write();
    TParameter<double>("Neff", Neff).Write();
    WriteMomentParams("mod_phi_all_pairs", m_phi_all_pairs);
    WriteMomentParams("mod_phi_best_rho", m_phi_best_rho);
    WriteMomentParams("mod_dphi_all_matchings", m_dphi_all_matchings);
    WriteMomentParams("mod_dphi_best_rho_sigma", m_dphi_best_rho_sigma);
    fout->Close();
    std::cout << "Output saved: " << output_file << std::endl;

    std::cout << "\n--- Reco M [GeV]  mean=" << h_M_reco->GetMean() << "  rms=" << h_M_reco->GetRMS() << std::endl;
    std::cout << "--- Reco y          mean=" << h_y_reco->GetMean() << "  rms=" << h_y_reco->GetRMS() << std::endl;
    std::cout << "--- Reco pT [GeV]   mean=" << h_pt_reco->GetMean() << "  rms=" << h_pt_reco->GetRMS() << std::endl;
    std::cout << "--- M_diff [GeV]    mean=" << h_M_diff->GetMean() << "  rms=" << h_M_diff->GetRMS() << std::endl;
    std::cout << "--- all OS pair mass mean=" << h_M_os_pairs->GetMean() << "  rms=" << h_M_os_pairs->GetRMS() << std::endl;
    std::cout << "--- best rho mass    mean=" << h_M_best_rho_pair->GetMean() << "  rms=" << h_M_best_rho_pair->GetRMS() << std::endl;
    std::cout << "--- best sigma mass  mean=" << h_M_best_sigma_pair->GetMean() << "  rms=" << h_M_best_sigma_pair->GetRMS() << std::endl;
    std::cout << "--- cos_theta_bose   mean=" << h_cos_theta_bose->GetMean() << "  rms=" << h_cos_theta_bose->GetRMS() << std::endl;

    std::cout << "\n=== Fourier modulation moments: f(x) ~ 1 + a_n cos(nx) + b_n sin(nx) ===" << std::endl;
    PrintMoment("phi all pi+pi- pairs, no truth", m_phi_all_pairs);
    PrintMoment("phi best rho-like pair, no truth", m_phi_best_rho);
    PrintMoment("Delta phi all opposite-sign matchings, no truth", m_dphi_all_matchings);
    PrintMoment("Delta phi best rho/sigma assignment, no truth", m_dphi_best_rho_sigma);

    std::cout << "--- cos2phi/cos1phi vs pT:" << std::endl;
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

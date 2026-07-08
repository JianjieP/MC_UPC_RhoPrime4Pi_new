// analyze_nonbose_dphi.C
// Full analysis of non-Bose Δφ distribution
// Step-by-step documentation of how Δφ is computed

void analyze_nonbose_dphi(
    const char* input = "build/BoseEvents_AuAu200_XnXn_nonbose_benchmark_100k.root"
) {
    // ── 0. Open file ──
    TFile* f = TFile::Open(input);
    TTree* t = (TTree*)f->Get("Events");

    Double_t event_weight, M_stored, pt_stored;
    TLorentzVector *pi1=nullptr, *pi2=nullptr, *pi3=nullptr, *pi4=nullptr;
    TLorentzVector *rho_mom=nullptr, *sigma_mom=nullptr;
    Int_t pi_parent[4], pi_charge[4];

    t->SetBranchAddress("event_weight", &event_weight);
    t->SetBranchAddress("M_rhoprime", &M_stored);
    t->SetBranchAddress("pt_rhoprime", &pt_stored);
    t->SetBranchAddress("pi1",&pi1); t->SetBranchAddress("pi2",&pi2);
    t->SetBranchAddress("pi3",&pi3); t->SetBranchAddress("pi4",&pi4);
    t->SetBranchAddress("rho_momentum", &rho_mom);
    t->SetBranchAddress("sigma_momentum", &sigma_mom);
    t->SetBranchAddress("pi_parent", pi_parent);
    t->SetBranchAddress("pi_charge", pi_charge);

    // ── Histograms ──
    TH1D* h_phi_rho  = new TH1D("h_phi_rho", ";#phi_{#pi^{+} in #rho RF} (rad);weighted", 40, -TMath::Pi(), TMath::Pi());
    TH1D* h_phi_sig  = new TH1D("h_phi_sig", ";#phi_{#pi^{+} in #sigma RF} (rad);weighted", 40, -TMath::Pi(), TMath::Pi());
    TH1D* h_dphi     = new TH1D("h_dphi",    ";#Delta#phi = #phi_{#rho} - #phi_{#sigma} (rad);weighted", 40, -TMath::Pi(), TMath::Pi());
    TH2D* h2d        = new TH2D("h2d","#phi_{#rho} vs #phi_{#sigma};#phi_{#rho};#phi_{#sigma}",
                                40, -TMath::Pi(), TMath::Pi(), 40, -TMath::Pi(), TMath::Pi());

    Long64_t n = t->GetEntries();
    cout << "=== Computing Δφ for " << n << " events ===\n";
    cout << "\nΔφ computation steps:\n";
    cout << "  1. Identify π⁺ from ρ (pi_parent==1, charge>0)\n";
    cout << "  2. Identify π⁺ from σ (pi_parent==2, charge>0)\n";
    cout << "  3. Boost π⁺_ρ to ρ RF:   pip_rf = boost(-rho.BoostVector())\n";
    cout << "  4. Boost π⁺_σ to σ RF:   same for sigma\n";
    cout << "  5. Define axes in LAB (no boost):\n";
    cout << "       Z = (0,0,1)  beam direction\n";
    cout << "       X = ρ' pT direction (unit vector in transverse plane)\n";
    cout << "       Y = Z × X\n";
    cout << "  6. φ_ρ = atan2(p_ρ·Y, p_ρ·X)  where p_ρ = 3-momentum of π⁺ in ρ RF\n";
    cout << "  7. φ_σ = atan2(p_σ·Y, p_σ·X)  where p_σ = 3-momentum of π⁺ in σ RF\n";
    cout << "  8. Δφ = φ_ρ - φ_σ  wrapped to [-π,π]\n\n";

    // Print one event in detail
    bool printed = false;

    for (Long64_t i = 0; i < n; ++i) {
        t->GetEntry(i);
        if (!isfinite(event_weight) || event_weight <= 0) continue;

        // Step 1-2: identify pions
        TLorentzVector* pions[4] = {pi1, pi2, pi3, pi4};
        TLorentzVector *pi_rho_plus=nullptr, *pi_sig_plus=nullptr;
        TLorentzVector *pi_rho_minus=nullptr, *pi_sig_minus=nullptr;

        for (int j = 0; j < 4; ++j) {
            if (pi_parent[j] == 1) {
                if (pi_charge[j] > 0) pi_rho_plus  = pions[j];
                else                  pi_rho_minus = pions[j];
            } else {
                if (pi_charge[j] > 0) pi_sig_plus  = pions[j];
                else                  pi_sig_minus = pions[j];
            }
        }
        if (!pi_rho_plus || !pi_sig_plus) continue;

        // Step 5: axes
        TLorentzVector p_rhoprime = *pi1 + *pi2 + *pi3 + *pi4;
        TVector3 vX = p_rhoprime.Vect(); vX.SetZ(0);
        if (vX.Mag() < 1e-12) continue;
        vX = vX.Unit();
        TVector3 vZ(0,0,1);
        TVector3 vY = vZ.Cross(vX).Unit();

        // Step 3-4, 6: φ_rho
        TVector3 beta_rho = rho_mom->BoostVector();
        TLorentzVector pip_rho_rf = *pi_rho_plus;
        pip_rho_rf.Boost(-beta_rho);
        TVector3 p3_rho = pip_rho_rf.Vect();
        double phi_rho = TMath::ATan2(p3_rho.Dot(vY), p3_rho.Dot(vX));

        // Step 7: φ_sigma
        TVector3 beta_sig = sigma_mom->BoostVector();
        TLorentzVector pip_sig_rf = *pi_sig_plus;
        pip_sig_rf.Boost(-beta_sig);
        TVector3 p3_sig = pip_sig_rf.Vect();
        double phi_sig = TMath::ATan2(p3_sig.Dot(vY), p3_sig.Dot(vX));

        // Step 8: Δφ
        double dphi = TVector2::Phi_mpi_pi(phi_rho - phi_sig);

        // Print first event
        if (!printed) {
            cout << "=== Example event 0 ===\n";
            cout << "  ρ' M=" << M_stored << " pT=" << pt_stored << "\n";
            cout << "  rho 4-mom: (" << rho_mom->Px() << "," << rho_mom->Py() << ","
                 << rho_mom->Pz() << "," << rho_mom->E() << ")\n";
            cout << "  sigma 4-mom: (" << sigma_mom->Px() << "," << sigma_mom->Py() << ","
                 << sigma_mom->Pz() << "," << sigma_mom->E() << ")\n";
            cout << "  X axis (ρ' pT dir): (" << vX.X() << "," << vX.Y() << "," << vX.Z() << ")\n";
            cout << "  Y axis: (" << vY.X() << "," << vY.Y() << "," << vY.Z() << ")\n";
            cout << "  π⁺_ρ in ρ RF: p3=(" << p3_rho.X() << "," << p3_rho.Y() << "," << p3_rho.Z() << ")\n";
            cout << "  π⁺_σ in σ RF: p3=(" << p3_sig.X() << "," << p3_sig.Y() << "," << p3_sig.Z() << ")\n";
            cout << "  φ_ρ = " << phi_rho << "  φ_σ = " << phi_sig << "  Δφ = " << dphi << "\n";
            cout << "  event_weight = " << event_weight << "\n\n";
            printed = true;
        }

        // Fill
        h_phi_rho->Fill(phi_rho, event_weight);
        h_phi_sig->Fill(phi_sig, event_weight);
        h_dphi->Fill(dphi, event_weight);
        h2d->Fill(phi_rho, phi_sig, event_weight);
    }
    pi1=pi2=pi3=pi4=rho_mom=sigma_mom=nullptr;

    // ── Fit ──
    auto doFit = [](TH1D* h, const char* label) {
        TF1* ff = new TF1("ff","[0]*(1+[1]*cos(x)+[2]*cos(2*x)+[3]*cos(3*x)+[4]*cos(4*x))",
                          -TMath::Pi(), TMath::Pi());
        double n0 = fabs(h->Integral())/(2*TMath::Pi())*h->GetBinWidth(1);
        ff->SetParameters(TMath::Max(n0,1.0), 0,0,0,0);
        h->Fit(ff,"RQ");
        cout << "\n=== " << label << " ===" << endl;
        for (int p=0;p<5;++p) cout << Form("  B%d = %+.4f ± %.4f",p,ff->GetParameter(p),ff->GetParError(p))<<endl;
    };

    doFit(h_dphi,    "Δφ (non-Bose)");
    doFit(h_phi_rho, "φ_rho (non-Bose)");
    doFit(h_phi_sig, "φ_sigma (non-Bose)");

    // ── Draw ──
    TCanvas* c = new TCanvas("c","non-bose dphi",1400,1000);
    c->Divide(2,2);
    c->cd(1); h_phi_rho->Draw("HIST E");
    c->cd(2); h_phi_sig->Draw("HIST E");
    c->cd(3); h_dphi->Draw("HIST E");
    c->cd(4); h2d->Draw("COLZ");
    c->SaveAs("plots_old/nonbose_dphi_full.png");

    f->Close();
    cout << "\nDone." << endl;
}

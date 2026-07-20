// analyze_nontruth_dphi.C
// Δφ WITHOUT using pi_parent — emulate realistic analysis
// Pairing: assign each π⁺ to a π⁻ to minimize |m(π⁺π⁻) - m_ρ|

void analyze_nontruth_dphi() {
    const double m_rho_nom = 0.77526;  // nominal rho mass

    TFile* f = TFile::Open("build/BoseEvents_AuAu200_XnXn_nonbose_benchmark_100k.root");
    TTree* t = (TTree*)f->Get("Events");

    Double_t event_weight;
    TLorentzVector *pi1=nullptr, *pi2=nullptr, *pi3=nullptr, *pi4=nullptr;
    Int_t pi_charge[4];

    t->SetBranchAddress("event_weight", &event_weight);
    t->SetBranchAddress("pi1",&pi1); t->SetBranchAddress("pi2",&pi2);
    t->SetBranchAddress("pi3",&pi3); t->SetBranchAddress("pi4",&pi4);
    t->SetBranchAddress("pi_charge", pi_charge);

    // Method 1: pair each π⁺ with its "best" π⁻ (closest to m_ρ)
    TH1D* hDphi_mass = new TH1D("hDphi_mass","#Delta#phi (best-mass pairing);#Delta#phi",
                                40, -TMath::Pi(), TMath::Pi());
    // Method 2: always use pairing (pi1,pi2)-(pi3,pi4) i.e. index-order
    TH1D* hDphi_idx = new TH1D("hDphi_idx","#Delta#phi (index-order);#Delta#phi",
                               40, -TMath::Pi(), TMath::Pi());

    // Also: just use lab φ of each π⁺ (no boost, no pairing)
    TH1D* hDphi_lab = new TH1D("hDphi_lab","#Delta#phi (lab frame, no pair RF);#Delta#phi",
                               40, -TMath::Pi(), TMath::Pi());

    auto getPhiInPairRF = [](TLorentzVector* pi_p, TLorentzVector* pi_m,
                              const TVector3& vX, const TVector3& vY) -> double {
        TLorentzVector pair = *pi_p + *pi_m;
        TVector3 beta = pair.BoostVector();
        TLorentzVector pi_rf = *pi_p;
        pi_rf.Boost(-beta);
        TVector3 p3 = pi_rf.Vect();
        return TMath::ATan2(p3.Dot(vY), p3.Dot(vX));
    };

    Long64_t n = t->GetEntries();
    for (Long64_t i = 0; i < n; ++i) {
        t->GetEntry(i);
        if (!isfinite(event_weight) || event_weight <= 0) continue;

        TLorentzVector* pions[4] = {pi1, pi2, pi3, pi4};

        // Find π⁺ and π⁻
        int ip_idx[2], im_idx[2];
        int nip=0, nim=0;
        for (int j=0; j<4; ++j) {
            if (pi_charge[j] > 0) ip_idx[nip++] = j;
            else im_idx[nim++] = j;
        }
        if (nip != 2 || nim != 2) continue;

        // Axes
        TLorentzVector p4pi = *pi1 + *pi2 + *pi3 + *pi4;
        TVector3 vX = p4pi.Vect(); vX.SetZ(0);
        if (vX.Mag() < 1e-12) continue;
        vX = vX.Unit();
        TVector3 vY = TVector3(0,0,1).Cross(vX).Unit();

        // ── Method 1: best-mass pairing ──
        // Try both assignments of π⁻ to π⁺, pick the one with best mass match
        // Assignment A: π⁺_a ↔ π⁻_c, π⁺_b ↔ π⁻_d
        // Assignment B: π⁺_a ↔ π⁻_d, π⁺_b ↔ π⁻_c
        int a=ip_idx[0], b=ip_idx[1], c=im_idx[0], d=im_idx[1];

        double m_ac = (*pions[a] + *pions[c]).M();
        double m_bd = (*pions[b] + *pions[d]).M();
        double m_ad = (*pions[a] + *pions[d]).M();
        double m_bc = (*pions[b] + *pions[c]).M();

        // Score: sum of |m - m_ρ| for both pairs
        double score_A = fabs(m_ac - m_rho_nom) + fabs(m_bd - m_rho_nom);
        double score_B = fabs(m_ad - m_rho_nom) + fabs(m_bc - m_rho_nom);

        int pA, mA, pB, mB;  // (pi+_A, pi-_A) and (pi+_B, pi-_B)
        if (score_A <= score_B) {
            pA=a; mA=c; pB=b; mB=d;
        } else {
            pA=a; mA=d; pB=b; mB=c;
        }

        double phi_A = getPhiInPairRF(pions[pA], pions[mA], vX, vY);
        double phi_B = getPhiInPairRF(pions[pB], pions[mB], vX, vY);
        // Order by pT
        double dphi_mass;
        if (pions[pA]->Pt() > pions[pB]->Pt())
            dphi_mass = TVector2::Phi_mpi_pi(phi_A - phi_B);
        else
            dphi_mass = TVector2::Phi_mpi_pi(phi_B - phi_A);
        hDphi_mass->Fill(dphi_mass, event_weight);

        // ── Method 2: index-order pairing (pi1,pi2)-(pi3,pi4) ──
        // pi1 and pi2 are π⁺ from ρ, pi3 and pi4 are π⁻ from σ in non-Bose
        // But we don't know that — we just use charge
        TLorentzVector *pip1=pions[a], *pip2=pions[b], *pim1=pions[c], *pim2=pions[d];
        double phi_i1 = getPhiInPairRF(pip1, pim1, vX, vY);
        double phi_i2 = getPhiInPairRF(pip2, pim2, vX, vY);
        double dphi_idx;
        if (pip1->Pt() > pip2->Pt())
            dphi_idx = TVector2::Phi_mpi_pi(phi_i1 - phi_i2);
        else
            dphi_idx = TVector2::Phi_mpi_pi(phi_i2 - phi_i1);
        hDphi_idx->Fill(dphi_idx, event_weight);

        // ── Method 3: lab-frame φ, no pair RF ──
        double phi_lab1 = pip1->Phi();
        double phi_lab2 = pip2->Phi();
        double dphi_lab = TVector2::Phi_mpi_pi(phi_lab1 - phi_lab2);
        hDphi_lab->Fill(dphi_lab, event_weight);
    }
    pi1=pi2=pi3=pi4=nullptr;

    auto doFit = [](TH1D* h, const char* label) {
        TF1* ff = new TF1("ff","[0]*(1+[1]*cos(x)+[2]*cos(2*x)+[3]*cos(3*x)+[4]*cos(4*x))",
                          -TMath::Pi(), TMath::Pi());
        double n0 = fabs(h->Integral())/(2*TMath::Pi())*h->GetBinWidth(1);
        ff->SetParameters(TMath::Max(n0,1.0),0,0,0,0);
        h->Fit(ff,"RQ");
        cout << "\n=== " << label << " ===" << endl;
        for (int p=0;p<5;++p) cout << Form("  B%d = %+.4f ± %.4f",p,ff->GetParameter(p),ff->GetParError(p))<<endl;
    };

    doFit(hDphi_mass, "Δφ (best-mass pairing, no MC truth)");
    doFit(hDphi_idx,  "Δφ (index-order pairing)");
    doFit(hDphi_lab,  "Δφ (lab frame, no boost)");

    TCanvas* c = new TCanvas("c","no truth",1200,400);
    c->Divide(3,1);
    c->cd(1); hDphi_mass->Draw("HIST E");
    c->cd(2); hDphi_idx->Draw("HIST E");
    c->cd(3); hDphi_lab->Draw("HIST E");
    c->SaveAs("plots_old/nontruth_dphi.png");

    f->Close();
    cout << "\nDone." << endl;
}

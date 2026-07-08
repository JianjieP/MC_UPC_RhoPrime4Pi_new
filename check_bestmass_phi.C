// check_bestmass_phi.C
// For each event: pick the charge pairing that minimizes |m(π⁺π⁻) - 0.770|
// Look at φ and Δφ in the chosen pairing

void check_bestmass_phi(const char* input, const char* label) {
    const double m_rho_nom = 0.770;  // target rho mass

    TFile* f = TFile::Open(input);
    TTree* t = (TTree*)f->Get("Events");

    Double_t event_weight;
    TLorentzVector *pi1=nullptr, *pi2=nullptr, *pi3=nullptr, *pi4=nullptr;
    Int_t pi_charge[4];

    t->SetBranchAddress("event_weight", &event_weight);
    t->SetBranchAddress("pi1",&pi1); t->SetBranchAddress("pi2",&pi2);
    t->SetBranchAddress("pi3",&pi3); t->SetBranchAddress("pi4",&pi4);
    t->SetBranchAddress("pi_charge", pi_charge);

    TH1D* h_phi_a  = new TH1D("h_phi_a", ";#phi of #pi^{+}_{A} in pair-A RF;", 40, -TMath::Pi(), TMath::Pi());
    TH1D* h_phi_b  = new TH1D("h_phi_b", ";#phi of #pi^{+}_{B} in pair-B RF;", 40, -TMath::Pi(), TMath::Pi());
    TH1D* h_dphi   = new TH1D("h_dphi",  ";#Delta#phi = #phi_A - #phi_B;", 40, -TMath::Pi(), TMath::Pi());

    // Also show which pairing was chosen: correct or wrong (using pi_parent if available)
    bool has_parent = (t->GetBranch("pi_parent") != nullptr);
    Int_t pi_parent[4];
    if (has_parent) t->SetBranchAddress("pi_parent", pi_parent);

    TH1D* h_dphi_correct = new TH1D("h_dphi_correct","#Delta#phi (correct pairing chosen)",40,-TMath::Pi(),TMath::Pi());
    TH1D* h_dphi_wrong   = new TH1D("h_dphi_wrong","#Delta#phi (wrong pairing chosen)",40,-TMath::Pi(),TMath::Pi());
    TH1D* h_choice = new TH1D("h_choice",";chosen correct? (1=yes 0=no)",2,0,2);

    auto getPhi = [](TLorentzVector* pip, TLorentzVector* pim,
                      const TVector3& vX, const TVector3& vY) -> double {
        TLorentzVector pair = *pip + *pim;
        TVector3 beta = pair.BoostVector();
        TLorentzVector pi_rf = *pip; pi_rf.Boost(-beta);
        TVector3 p3 = pi_rf.Vect();
        return TMath::ATan2(p3.Dot(vY), p3.Dot(vX));
    };

    Long64_t n = t->GetEntries();
    for (Long64_t i = 0; i < n; ++i) {
        t->GetEntry(i);
        if (!isfinite(event_weight)) continue;

        TLorentzVector* pions[4] = {pi1, pi2, pi3, pi4};
        int ip[2], im[2], nip=0, nim=0;
        for (int j=0; j<4; ++j) {
            if (pi_charge[j]>0) ip[nip++]=j; else im[nim++]=j;
        }
        if (nip!=2||nim!=2) continue;

        TLorentzVector p4pi = *pi1+*pi2+*pi3+*pi4;
        TVector3 vX = p4pi.Vect(); vX.SetZ(0);
        if (vX.Mag() < 1e-12) continue;
        vX = vX.Unit();
        TVector3 vY = TVector3(0,0,1).Cross(vX).Unit();

        int a=ip[0], b=ip[1], c=im[0], d=im[1];

        double m_ac = (*pions[a] + *pions[c]).M();
        double m_bd = (*pions[b] + *pions[d]).M();
        double m_ad = (*pions[a] + *pions[d]).M();
        double m_bc = (*pions[b] + *pions[c]).M();

        double score_A = fabs(m_ac - m_rho_nom) + fabs(m_bd - m_rho_nom);
        double score_B = fabs(m_ad - m_rho_nom) + fabs(m_bc - m_rho_nom);

        bool use_A = (score_A <= score_B);

        int pA, mA, pB, mB;
        if (use_A) { pA=a; mA=c; pB=b; mB=d; }
        else       { pA=a; mA=d; pB=b; mB=c; }

        double phi_A = getPhi(pions[pA], pions[mA], vX, vY);
        double phi_B = getPhi(pions[pB], pions[mB], vX, vY);
        double dphi = TVector2::Phi_mpi_pi(phi_A - phi_B);

        h_phi_a->Fill(phi_A, event_weight);
        h_phi_b->Fill(phi_B, event_weight);
        h_dphi->Fill(dphi, event_weight);

        // Check if pairing is "correct" (using MC truth)
        if (has_parent) {
            bool is_correct = 
                (pi_parent[pA]==1 && pi_parent[pB]==2) ||
                (pi_parent[pA]==2 && pi_parent[pB]==1);
            h_choice->Fill(is_correct ? 1.5 : 0.5, event_weight);
            if (is_correct) h_dphi_correct->Fill(dphi, event_weight);
            else h_dphi_wrong->Fill(dphi, event_weight);
        }
    }
    pi1=pi2=pi3=pi4=nullptr;

    auto doFit = [&](TH1D* h, const char* tag) {
        TF1* ff = new TF1("ff","[0]*(1+[1]*cos(x)+[2]*cos(2*x)+[3]*cos(3*x)+[4]*cos(4*x))",
                          -TMath::Pi(), TMath::Pi());
        double n0 = fabs(h->Integral())/(2*TMath::Pi())*h->GetBinWidth(1);
        ff->SetParameters(TMath::Max(n0,1.0),0,0,0,0);
        h->Fit(ff,"RQ");
        cout << "\n=== " << label << " " << tag << " ===" << endl;
        for (int p=0;p<5;++p) cout << Form("  B%d = %+.4f ± %.4f",p,ff->GetParameter(p),ff->GetParError(p))<<endl;
    };

    doFit(h_phi_a, "φ_A");
    doFit(h_phi_b, "φ_B");
    doFit(h_dphi,  "Δφ");

    if (has_parent) {
        cout << "\n--- Pairing choice correctness ---" << endl;
        cout << "  Fraction correct = " << h_choice->GetBinContent(2)/(h_choice->GetBinContent(1)+h_choice->GetBinContent(2)) << endl;
        doFit(h_dphi_correct, "Δφ (correct pairing)");
        doFit(h_dphi_wrong,   "Δφ (wrong pairing)");
    }

    TCanvas* c = new TCanvas("c",label, 1200, 400);
    c->Divide(3,1);
    c->cd(1); h_phi_a->Draw("HIST E");
    c->cd(2); h_phi_b->Draw("HIST E");
    c->cd(3); h_dphi->Draw("HIST E");
    c->SaveAs(Form("plots_old/bestmass_%s.png", label));

    f->Close();
}

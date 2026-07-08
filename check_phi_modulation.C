// check_phi_modulation.C
// Single-pion φ in pair RF: check if φ_rho and φ_sigma have modulation
// Z = lab Z, X = ρ' pT direction

void check_phi_modulation(const char* input, const char* label) {
    TFile* f = TFile::Open(input);
    TTree* t = (TTree*)f->Get("Events");

    Double_t event_weight;
    TLorentzVector *pi1=nullptr, *pi2=nullptr, *pi3=nullptr, *pi4=nullptr;
    Int_t pi_parent[4], pi_charge[4];

    bool has_parent = (t->GetBranch("pi_parent") != nullptr);
    if (has_parent) {
        t->SetBranchAddress("pi_parent", pi_parent);
        t->SetBranchAddress("pi_charge", pi_charge);
    } else {
        t->SetBranchAddress("pi_charge", pi_charge);
    }
    t->SetBranchAddress("event_weight", &event_weight);
    t->SetBranchAddress("pi1",&pi1); t->SetBranchAddress("pi2",&pi2);
    t->SetBranchAddress("pi3",&pi3); t->SetBranchAddress("pi4",&pi4);

    TH1D* h_phi_rho = new TH1D("h_phi_rho","#phi of #pi^{+}_{#rho} in #rho RF", 40, -TMath::Pi(), TMath::Pi());
    TH1D* h_phi_sig = new TH1D("h_phi_sig","#phi of #pi^{+}_{#sigma} in #sigma RF", 40, -TMath::Pi(), TMath::Pi());
    TH1D* h_phi_all = new TH1D("h_phi_all","#phi of all #pi^{+} in pair RF (both pairings)", 40, -TMath::Pi(), TMath::Pi());

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

        TLorentzVector p4pi = *pi1+*pi2+*pi3+*pi4;
        TVector3 vX = p4pi.Vect(); vX.SetZ(0);
        if (vX.Mag() < 1e-12) continue;
        vX = vX.Unit();
        TVector3 vY = TVector3(0,0,1).Cross(vX).Unit();

        if (has_parent) {
            // Use MC truth
            TLorentzVector *pi_rho_p=nullptr, *pi_rho_m=nullptr;
            TLorentzVector *pi_sig_p=nullptr, *pi_sig_m=nullptr;
            for (int j=0; j<4; ++j) {
                if (pi_parent[j] == 1) {
                    if (pi_charge[j]>0) pi_rho_p=pions[j]; else pi_rho_m=pions[j];
                } else {
                    if (pi_charge[j]>0) pi_sig_p=pions[j]; else pi_sig_m=pions[j];
                }
            }
            if (pi_rho_p && pi_rho_m)
                h_phi_rho->Fill(getPhi(pi_rho_p, pi_rho_m, vX, vY), event_weight);
            if (pi_sig_p && pi_sig_m)
                h_phi_sig->Fill(getPhi(pi_sig_p, pi_sig_m, vX, vY), event_weight);
        }

        // All pairings (model-independent)
        int ip[2], im[2], nip=0, nim=0;
        for (int j=0; j<4; ++j) {
            if (pi_charge[j]>0) ip[nip++]=j; else im[nim++]=j;
        }
        if (nip!=2||nim!=2) continue;

        // Both pairings
        h_phi_all->Fill(getPhi(pions[ip[0]], pions[im[0]], vX, vY), event_weight);
        h_phi_all->Fill(getPhi(pions[ip[0]], pions[im[1]], vX, vY), event_weight);
        h_phi_all->Fill(getPhi(pions[ip[1]], pions[im[0]], vX, vY), event_weight);
        h_phi_all->Fill(getPhi(pions[ip[1]], pions[im[1]], vX, vY), event_weight);
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

    if (has_parent) {
        doFit(h_phi_rho, "φ_rho (MC truth)");
        doFit(h_phi_sig, "φ_sigma (MC truth)");
    }
    doFit(h_phi_all, "φ all (4 pairs, no truth)");

    TCanvas* c = new TCanvas("c",label, has_parent?1200:600, 500);
    c->Divide(has_parent?3:1,1);
    int pad=1;
    if (has_parent) { c->cd(pad++); h_phi_rho->Draw("HIST E"); }
    if (has_parent) { c->cd(pad++); h_phi_sig->Draw("HIST E"); }
    c->cd(pad); h_phi_all->Draw("HIST E");
    c->SaveAs(Form("plots_old/phi_%s.png", label));

    f->Close();
}

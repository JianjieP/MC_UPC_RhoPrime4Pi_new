// analyze_allpairings.C
// Fill BOTH charge pairings for each event (no selection)
// Compare non-Bose vs Bose

void analyze_allpairings(const char* input, const char* label) {
    TFile* f = TFile::Open(input);
    TTree* t = (TTree*)f->Get("Events");

    Double_t event_weight;
    TLorentzVector *pi1=nullptr, *pi2=nullptr, *pi3=nullptr, *pi4=nullptr;
    Int_t pi_charge[4];

    t->SetBranchAddress("event_weight", &event_weight);
    t->SetBranchAddress("pi1",&pi1); t->SetBranchAddress("pi2",&pi2);
    t->SetBranchAddress("pi3",&pi3); t->SetBranchAddress("pi4",&pi4);
    t->SetBranchAddress("pi_charge", pi_charge);

    TH1D* hDphi = new TH1D("hDphi","#Delta#phi (both pairings);#Delta#phi", 40, -TMath::Pi(), TMath::Pi());

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
            if (pi_charge[j] > 0) ip[nip++] = j;
            else im[nim++] = j;
        }
        if (nip!=2 || nim!=2) continue;

        TLorentzVector p4pi = *pi1+*pi2+*pi3+*pi4;
        TVector3 vX = p4pi.Vect(); vX.SetZ(0);
        if (vX.Mag() < 1e-12) continue;
        vX = vX.Unit();
        TVector3 vY = TVector3(0,0,1).Cross(vX).Unit();

        // Pairing 1: (ip0,im0)-(ip1,im1)
        double phi_a1 = getPhi(pions[ip[0]], pions[im[0]], vX, vY);
        double phi_b1 = getPhi(pions[ip[1]], pions[im[1]], vX, vY);
        hDphi->Fill(TVector2::Phi_mpi_pi(phi_a1 - phi_b1), event_weight);

        // Pairing 2: (ip0,im1)-(ip1,im0)
        double phi_a2 = getPhi(pions[ip[0]], pions[im[1]], vX, vY);
        double phi_b2 = getPhi(pions[ip[1]], pions[im[0]], vX, vY);
        hDphi->Fill(TVector2::Phi_mpi_pi(phi_a2 - phi_b2), event_weight);
    }
    pi1=pi2=pi3=pi4=nullptr;

    TF1* ff = new TF1("ff","[0]*(1+[1]*cos(x)+[2]*cos(2*x)+[3]*cos(3*x)+[4]*cos(4*x))",
                      -TMath::Pi(), TMath::Pi());
    double n0 = fabs(hDphi->Integral())/(2*TMath::Pi())*hDphi->GetBinWidth(1);
    ff->SetParameters(TMath::Max(n0,1.0),0,0,0,0);
    hDphi->Fit(ff,"RQ");
    cout << "\n=== " << label << ": both pairings ===" << endl;
    for (int p=0;p<5;++p) cout << Form("  B%d = %+.4f ± %.4f",p,ff->GetParameter(p),ff->GetParError(p))<<endl;

    TCanvas* c = new TCanvas("c",label,600,500);
    hDphi->Draw("HIST E");
    c->SaveAs(Form("plots_old/dphi_allpairs_%s.png", label));

    f->Close();
}

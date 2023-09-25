void finder() {
    fstream file;
	file.open("finder_check_BG10_cov15.txt", ios::in);

    TH1F *h_mom_all = new TH1F("mom_focus", "Any #chi^{2}; p_{fit}-p_{real} [MeV]; Tracks", 100, -40, 40);
    TH1F *h_mom_focus = new TH1F("mom_all", "#chi^{2} < 50; p_{fit}-p_{real} [MeV]; Tracks", 100, -40, 40);
    TH1F *h_ultra = new TH1F("mom", "#chi^{2} < 7; p_{fit}-p_{real} [MeV]; Tracks", 100, -40, 40);
    auto hcol1 = new TH2F("Fit","Fitted momentum; p_{fit}-p_{real} [MeV/c]; Fitting #chi^{2}", 100, -10, 10, 200, 0, 300);
    
    auto hcol2 = new TH2F("Recon","Reconstructed measurements; p_{fit}-p_{real} [MeV/c]; Reconstructed nHits", 100, -20, 20, 17, 15, 32);
    auto hcol3 = new TH2F("Recon2","Reconstructed measurements; Reconstructed nHits; Fitting #chi^{2}", 6, 15, 32, 200, 0, 300);
    
    auto hcol4 = new TH2F("Fit2","Fitted momentum; p_{fit} [MeV/c]; p_{estim} [MeV/c]", 100, 30, 160, 200, 30, 160);
    
    TH1F *h_estim = new TH1F("estim", "Any #chi^{2}; p_{estim}-p_{real} [MeV]; Tracks", 100, -80, 80);
    
    TH1F *h_q = new TH1F("Recon3", "Amount of reconstructed hits; nHits; Entries", 32, 0, 32);

    auto hcol5 = new TH2F("Purity","Purity; p_{fit}-p_{estim} [MeV/c]; Purity", 100, -40, 40, 25, 0, 1);
    
    int i, ok;
    double p, p_og, ndf, chi2, p_init, purity;
    
    int nOk    = 0;
    int nHitsOK = 0;
    int nFocus = 0;
    int nTotal = 0;
    int nPure = 0;
    int n7 = 0;

    while(1){
		file >> i >> ok >> ndf >> chi2 >> p >> p_init >> purity;
    p_og = p_init;
		
		double quality = chi2 / ndf;
		nTotal++;
		if (ok) {
			nOk++;
      if (i != 800000) {
			if (ndf >= 30) {
				h_mom_all->Fill(p - p_og);
				nHitsOK++;
				if (chi2 / ndf < 50) {
					nFocus++;
					h_mom_focus->Fill(p - p_og);
					if (chi2 / ndf < 7) {
						h_ultra->Fill(p-p_og);
            n7++;
					}
				}
				hcol1->Fill(p-p_og, quality);
			}
			h_q->Fill(ndf / 2);
			hcol2->Fill(p - p_og, ndf/2);
			hcol3->Fill(ndf/2, quality);
			hcol4->Fill(p, p_init*1000);
			h_estim->Fill(p_init*1000 - p_og);
      if (i != 800000) {
        hcol5->Fill(p-p_og, purity);
        if (purity > 0.5) {nPure++;}
      }
      else {std::cout << "BG has been seeded" << std::endl;}
		}
    }

		if ( file.eof()) break;
	}
	
	std::cout << "Total number of events: " << nTotal << "\n" <<
		"Total number of successfully fitted events: " << nOk << "\n" <<
		"Total number of events with nHits > 15: " << nHitsOK << "\n" <<
		"Total number of events with chi2 < 50 and nHits > 15: " << nFocus << std::endl;
    std::cout << nPure << std::endl;

    TCanvas *c1 = new TCanvas("c1", "c1", 1500, 1000);
    h_mom_all->Draw();
    h_mom_all->SetStats(0);
     TF1 *fit = new TF1("fit", "gaus", -30, 30);
    //h_mom_all->Fit("fit", "R");
    h_mom_focus->Draw("SAME");
    h_mom_focus->SetLineColor(kRed);
    h_ultra->Draw("SAME");
    h_ultra->SetLineColor(6);
    gPad->BuildLegend();
    gPad->SetLogy();
    //h_mom_focus->SetLineWidth(3);
    
    TCanvas *c2 = new TCanvas("c2", "c2", 1500, 1000);
    hcol2->Draw("COLZ");
    hcol2->SetStats(0);
    
    TCanvas *c3 = new TCanvas("c3", "c3", 1500, 1000);
    h_q->Draw();
    h_q->SetStats(0);
    
    TCanvas *c4 = new TCanvas("c4", "c4", 1500, 1000);
    hcol1->Draw("COLZ");
    hcol1->SetStats(0);
    
    TCanvas *c5 = new TCanvas("c5", "c5", 1500, 1000);
    h_ultra->Draw();
    h_ultra->SetStats(0);
    
    TCanvas *c6 = new TCanvas("c6", "c6", 1500, 1000);
    hcol3->Draw("COLZ");
    hcol3->SetStats(0);
    
    TCanvas *c7 = new TCanvas("c7", "c7", 1500, 1000);
    hcol4->Draw("COLZ");
    hcol4->SetStats(0);
    
    TCanvas *c8 = new TCanvas("c8", "c8", 1500, 1000);
    h_estim->Draw();
    h_estim->SetStats(0);

    TCanvas *c9 = new TCanvas("c9", "c9", 1500, 1000);
    hcol5->Draw("COLZ");
    hcol5->SetStats(0);
}

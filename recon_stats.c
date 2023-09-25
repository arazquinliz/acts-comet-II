void recon_stats() {
    fstream file, file2;
	file.open("recon_mom.txt", ios::in);

    TH1F *h_mom = new TH1F("Momentum", "Reconstructed momentum; p [MeV]; Tracks", 100, 40, 160);
    TH1F *h_init = new TH1F("Initial momentum", "Estimated initial momentum; |p_{estim} - 105| [MeV]; , Entries", 100, 0, 100);
	auto hcol1 = new TH2F("Recon","Reconstructed measurements; Reconstructed momentum [MeV/c]; Reconstructed nHits", 100, 40, 160, 17, 15, 32);
	auto hcol2 = new TH2F("Recon","Reconstructed measurements; Reconstructed momentum [MeV/c]; Initial momentum", 100, 40, 160, 100, 40, 160);
	auto hcol3 = new TH2F("Recon","Reconstructed measurements; Reconstructed momentum [MeV/c]; Initial momentum", 100, 40, 160, 10, 0, 1.1);
	// Find another way to represent the purity

    double p, chi2;
    double ndf, init_p, purity;
	int nHits, event;

	double puri = 0.;
	int n = 0;

    while(1){
		file >> p >> ndf >> chi2 >> init_p >> nHits >> purity >> event;
		//if (ndf > 40) {
            if ((init_p > 40) && (init_p < 160)) {
            		h_mom->Fill(p);
            		h_init->Fill(init_p - 105);
					hcol1->Fill(p, nHits);
					hcol2->Fill(p, init_p);
					puri += purity;
					n += 1;
            }
        //}
		if ( file.eof()) break;
	}

	// Fitting functions
	// Gauss
    TF1 *fit = new TF1("fit", "gaus", 40, 160);
    
    // Crystal ball
    TF1 *fitFunction = new TF1("fitFunction", "crystalball", 40, 160);
	fitFunction->SetParameters(105, 100, 20, 0.1, 50); // x, mean, sigma, alpha, N,
	//fitFunction->Print();
	
	// Exponentially convoluted gaussian
	TF1Convolution *f_conv = new TF1Convolution("expo", "gaus", 40., 160., true); // f1, f2, xmin, xmax, use fast fourier transform
	f_conv->SetNofPointsFFT(1000);
	TF1 *f = new TF1("f", *f_conv, 40, 160, f_conv->GetNpar());
	f->SetParameters(0.5, -0.3, 2., 10.);
	
	// Double gaussian 
	double par[6];
	TF1 *G1 = new TF1("G1", "gaus", 90, 110);
	TF1 *G2 = new TF1("G2", "gaus", 40, 160);
	G1->SetLineColor(kBlue);
	G2->SetLineColor(kGreen);
	
	//h_mom->Fit(G1, "R");
	//h_mom->Fit(G2);
	
	G1->GetParameters(&par[0]);
	G2->GetParameters(&par[3]);
	
	TF1 *fitDouble = new TF1("fitDouble", "gaus(0)+gaus(3)", 40, 160);
	fitDouble->SetLineColor(kRed);
	fitDouble->SetParameters(par);

	// Crystal ball + gaussian
	TF1 *fitFunc = new TF1("fitFunc", "gaus(0) + crystalball(3)", 40, 160);
	fitFunc->SetLineColor(kBlue);
	fitFunc->SetParameters(80, 93, 30, 105, 100, 20, 0.1, 50); // x, mean, sigma, alpha, N,
	
	// Plot and fitting
    TCanvas *c1 = new TCanvas("c1", "c1", 1500, 1000);
    h_mom->Draw();
    //h_mom->Fit("fitFunction", "R");
    //h_mom->Fit("f", "R");
	//h_mom->Fit("fit", "R");
	//h_mom->Fit(fitDouble);
	h_mom->Fit("fitFunc", "R");
	
	
    //TCanvas *c2 = new TCanvas("c2", "c2", 1500, 1000);
    //h_init->Draw();

	/*TCanvas *c2 = new TCanvas("c2", "c2", 1500, 1000);
    hcol2->Draw("COLZ");
	hcol2->Draw("COLZ");
	TLine *lX = new TLine(105, 40, 105, 160);
	lX->Draw();
	lX->SetLineColor(kRed);
	TLine *lY = new TLine(40, 105, 160, 105);
	lY->Draw();
	lY->SetLineColor(kRed);*/

    //file.close();
}

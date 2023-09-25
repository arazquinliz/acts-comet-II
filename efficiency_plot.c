void efficiency_plot() {
    // 9900 total number of events with enough hits
    // files were created only with the first station
    fstream noCut, tCut, pCut, allCut;
	noCut.open("seedEfficiency_allCut_0.txt", ios::in);
    tCut.open("seedEfficiency_noBG_095.txt", ios::in);
    pCut.open("seedEfficiency_pCut_all.txt", ios::in);
    allCut.open("seedEfficiency_allCut_095.txt", ios::in);

    TEfficiency* noEff  = new TEfficiency("Efficiency base cuts","ACTS seeder w/ backgrounds; pT [MeV/c]; Efficiency", 50, 0, 105);
    TEfficiency* tEff   = new TEfficiency("Efficiency t cut","Seeder w/o backgrounds; pT [MeV/c]; Efficiency", 50, 0.1, 105.1);
    TEfficiency* pEff   = new TEfficiency("Efficiency p cut","Base + p cuts; pT [MeV/c]; Efficiency", 50, 0, 105);
    TEfficiency* allEff = new TEfficiency("Efficiency all cuts","All cuts w/ backgrounds; pT [MeV/c]; Efficiency", 50, 0.1, 105.1);

    TH1F *allFake       = new TH1F("allFake", "All cuts; pT [MeV/c] signal e; Fake seeds per run", 50, 0, 105);
    TH1F *noFake        = new TH1F("noFake", "Base cuts; pT [MeV/c] signal e; Fake seeds per run", 50, 0, 105);
    TH1F *tFake         = new TH1F("tFake", "Base + t cuts; pT [MeV/c] signal e; Fake seeds per run", 50, 0, 105);
    TH1F *pFake         = new TH1F("pFake", "Base + p cuts; pT [MeV/c] signal e; Fake seeds per run", 50, 0, 105);

    TH1F *allBG       = new TH1F("allBG", "All cuts; pT [MeV/c] signal electron; True BG seeds per run", 50, 0, 105);
    TH1F *noBG        = new TH1F("noBG", "Base cuts; pT [MeV/c] signal electron; True BG seeds per run", 50, 0.1, 105.1);
    TH1F *tBG         = new TH1F("tBG", "Base + t cuts; pT [MeV/c] signal electron; True BG seeds per run", 50, 0, 105);
    TH1F *pBG         = new TH1F("pBG", "Base + p cuts; pT [MeV/c] signal electron; True BG seeds per run", 50, 0.1, 105.1);
    
    int i, type, eff;
    double pT, mu, p;

    int nFake = 0, nBG = 0, nSignal = 0, nEmpty = 0;
    int prev = -1;
    while(1) {
		noCut >> i >> type >> eff >> p >> pT >> mu;
        if (type == 0) {
            if (eff == 1) {
            		if (i != prev) {
				allEff->Fill(true, pT);
				prev = i;
			}
			nSignal++;
            }
            else if (eff == 0) {
                allEff->Fill(false, pT);
                nEmpty++;
            }
        }
        if (type == 2) {
            allFake->Fill(pT);
            nFake++;
        }
        if (type == 1) {
            allBG->Fill(pT);
            nBG++;
        }
		if ( noCut.eof()) {break;}
	}
    std::cout << nSignal << " " << nEmpty << " " << nBG << " " << nFake << std::endl;
    
    prev = -1;
    while(1) {
		noCut >> i >> type >> eff >> p >> pT >> mu;
        if (type == 0) {
            if (eff == 1) {
		        if (i != prev) {
				noEff->Fill(true, pT);
				prev = i;
			}
			nSignal++;
            }
            else if (eff == 0) {
                noEff->Fill(false, pT);
            }
        }
        if (type == 2) {
            noFake->Fill(pT);
        }
        if (type == 1) {
            noBG->Fill(pT);
        }
		if ( noCut.eof()) {break;}
	}
	
	prev = -1;
    while(1) {
		tCut >> i >> type >> eff >> p >> pT >> mu;
        if (type == 0) {
            if (eff == 1) {
		        if (i != prev) {
				tEff->Fill(true, pT);
				prev = i;
			}
			nSignal++;
            }
            else if (eff == 0) {
                tEff->Fill(false, pT);
            }
        }
        if (type == 2) {
            tFake->Fill(pT);
        }
        if (type == 1) {
            tBG->Fill(pT);
        }
		if ( tCut.eof()) {break;}
	}
    /*while(1) {
		pCut >> i >> type >> eff >> p >> pT >> mu;
        if (type == 0) {
            if (eff == 1) {
		        pEff->Fill(true, pT);
            }
            else if (eff == 0) {
                pEff->Fill(false, pT);
            }
        }
        if (type == 2) {
            pFake->Fill(pT);
        }
        if (type == 1) {
            pBG->Fill(pT);
        }
		if ( pCut.eof()) {break;}
	}*/
    
    TCanvas* c1 = new TCanvas("c1","",600,400);
    c1->SetFillStyle(1001);
    c1->SetFillColor(kWhite);
    tEff->Draw("AP");
    noEff->Draw("SAME");
    noEff->SetLineColor(2);
    //pEff->Draw("SAME");
    //pEff->SetLineColor(8);
    allEff->Draw("SAME");
    allEff->SetLineColor(4);
    gPad->BuildLegend();

    /*TCanvas* c2 = new TCanvas("c2", "", 600, 400);
    noFake->Scale(1./9900);
    tFake->Scale(1./9900);
    pFake->Scale(1./9900);
    allFake->Scale(1./9900);
    noFake->Draw("HIST");
    noFake->SetStats(0);
    noFake->SetLineColor(1);
    tFake->Draw("SAME HIST");
    tFake->SetLineColor(2);
    pFake->Draw("SAME HIST");
    pFake->SetLineColor(8);
    allFake->Draw("SAME HIST");
    allFake->SetLineColor(4);
    gPad->BuildLegend();
    c2->SetLogy();

    TCanvas* c3 = new TCanvas("c3", "", 600, 400);
    noBG->Scale(1./9900);
    tBG->Scale(1./9900);
    pBG->Scale(1./9900);
    allBG->Scale(1./9900);
    noBG->Draw("HIST");
    noBG->SetStats(0);
    noBG->SetLineColor(1);
    tBG->Draw("SAME HIST");
    tBG->SetLineColor(2);
    pBG->Draw("SAME HIST");
    pBG->SetLineColor(8);
    allBG->Draw("SAME HIST");
    allBG->SetLineColor(4);
    gPad->BuildLegend();
    c3->SetLogy();*/
}

// Plot Kalman Fitting with measurements
void plot_fitting() {
        std::fstream file;
        std::fstream fitfile;
        file.open("linked_meas.txt", std::ios::in);
        fitfile.open("predicted_meas.txt", std::ios::in);

        TCanvas *c0 = new TCanvas();
        c0->Divide(2,1);
        
        TGraph *g1 = new TGraph();
        g1->SetTitle("Hits in StrawTrk; x; y");
        TGraph *g2 = new TGraph();
        g2->SetTitle("Hits in StrawTrk; z; y");
        TGraph *g3 = new TGraph();
        g3->SetTitle("Hits in StrawTrk; z; y");
        TGraph *g4 = new TGraph();
        g4->SetTitle("Hits in StrawTrk; x; y");

        //TCanvas *c1 = new TCanvas();
        TGraph2D *gr = new TGraph2D();
        gr->SetName("hits_ST");
        gr->SetTitle("Hits in StrawTrk; x; z; y");
        TGraph2D *gr2 = new TGraph2D();
        gr2->SetName("hits_ST");
        gr2->SetTitle("Hits in StrawTrk; x; z; y");

        int j = 0;
        double x, y, z, t, r;
        double fit_x, fit_y, fit_z, dir_x, dir_y, dir_z, mom;
        while (1) {
                /*file >> x >> y >> z >> t;
                g1->SetPoint(j, x-10800, y); 
                g2->SetPoint(j, z-3334.5, y); */
                file >> x >> y >> z;
                g1->SetPoint(j, x, y); 
                g2->SetPoint(j, z, x);
                gr->SetPoint(j, x, y, z);
                ++j;
                if ( file.eof()) {break;}
        }

        g1->GetXaxis()->SetRangeUser(-695, 695); // in mm
        //g1->GetYaxis()->SetRangeUser(-695, 695);

        //g1->SetMarkerStyle(5);
        g1->SetMarkerColor(kBlue);

        g2->GetXaxis()->SetRangeUser(-1400, 1400);
        g2->GetYaxis()->SetRangeUser(-695, 695);

        //g2->SetMarkerStyle(5);
        g2->SetMarkerColor(kBlue);
        //g2->SetMarkerSize(2);

        gr->GetXaxis()->SetRangeUser(-695, 695);
        gr->GetZaxis()->SetRangeUser(-1400, 1400);
        gr->GetYaxis()->SetRangeUser(-695, 695);
        //gr->Draw();

        //TLine *l = new TLine(0,-695,0,695);

        //c0->cd(1);
        //g1->Draw("AP");
        //l->Draw();

        //gr->Draw("p0");
        
        int i = 0;
        while (1) {
                fitfile >> fit_x >> fit_y >> fit_z;
                //std::cout << fit_x << std::endl;
                g3->SetPoint(i, fit_z, fit_x); 
                g4->SetPoint(i, fit_x, fit_y);
                gr2->SetPoint(i, fit_x, fit_y, fit_z);
                i++;
                if ( fitfile.eof()) {break;}
        }
        
        g3->GetXaxis()->SetRangeUser(-1400, 1400);
        g3->GetYaxis()->SetRangeUser(-695, 695);

        g3->SetMarkerStyle(4);
        g3->SetMarkerColor(kRed);
        
        //g4->GetXaxis()->SetRangeUser(-695, 695); // in mm
        //g4->GetYaxis()->SetRangeUser(-695, 695);

        g4->SetMarkerStyle(4);
        g4->SetMarkerColor(kRed);
        
        gr2->GetXaxis()->SetRangeUser(-695, 695);
        gr2->GetZaxis()->SetRangeUser(-1400, 1400);
        gr2->GetYaxis()->SetRangeUser(-695, 695);
        
        //gr2->Draw("p0");
        c0->cd(1);
        g2->Draw("A*");
        g3->Draw("p");
        
        int mani, lay;
        int NofStations  = 8;
        double StationD  = 390;
        double manifoldD = 5.9;
        double layerD    = 9.35;
        for (int i = 0; i < NofStations; i++) {
        	int manifold = 0;
        	while (manifold < 2) {
        		if (manifold == 0) {mani = -1; lay = -1;}
            		else if (manifold == 1) {mani = 1; lay = 1;}
		    	int layer = 0;
		    	while (layer < 2) {
        			TLine *l = new TLine(-1365 + i * StationD + mani * manifoldD + lay*layerD*layer, -695, -1365 + i * StationD + mani * manifoldD + lay*layerD*layer, 695);
        		l->Draw();
        			//std::cout << -1365 + i * StationD + mani * manifoldD + lay*layerD*0.5 << std::endl;
        			layer++;
        		} // layer in manifold
        		manifold++;
        	} // manifold in station
        } // station
        
        c0->cd(2); 
        g1->Draw("A*");
        g4->Draw("p");
}

// Plot measurements on surface

void plot_sp() {
        std::fstream file;
        file.open("spacepoints.txt", std::ios::in);
        //file.open("event_meas.txt", std::ios::in);

        TCanvas *c0 = new TCanvas();
        c0->Divide(2,1);
        TGraph *g1 = new TGraph();
        g1->SetTitle("Hits in StrawTrk; x; y");
        TGraph *g2 = new TGraph();
        g2->SetTitle("Hits in StrawTrk; z; x");

        //TCanvas *c1 = new TCanvas();
        TGraph2D *gr = new TGraph2D();
        gr->SetName("hits_ST");
        gr->SetTitle("Hits in StrawTrk; x; z; y");

        int j = 0;
        double x, y, z, t, r;
        int vol, layer, surf;
        while (1) {
                /*file >> x >> y >> z >> t;
                g1->SetPoint(j, x-10800, y); 
                g2->SetPoint(j, z-3334.5, y); */
                file >> x >> y >> z;
                g1->SetPoint(j, x, y); 
                g2->SetPoint(j, z, x);
                gr->SetPoint(j, x, z, y);
                ++j;
                //if ((x == 0.147148) & (y == -0.0889063) & (z == -1380)) {break;}
                if ( file.eof()) {break;}
        }

        g1->GetXaxis()->SetRangeUser(-695, 695); // in mm
        g1->GetYaxis()->SetRangeUser(-695, 695);

        g1->SetMarkerStyle(5);
        g1->SetMarkerColor(kBlue);

        g2->GetXaxis()->SetRangeUser(-1400, 1400);
        g2->GetYaxis()->SetRangeUser(-695, 695);

        //g2->SetMarkerStyle(5);
        g2->SetMarkerColor(kBlue);
        //g2->SetMarkerSize(2);

        gr->GetXaxis()->SetRangeUser(0, 2760);
        gr->GetZaxis()->SetRangeUser(-695, 695);
        gr->GetYaxis()->SetRangeUser(-695, 695);
        //gr->Draw();

        //TLine *l = new TLine(0,-695,0,695);

        c0->cd(1);
        g2->Draw("A*");
        
        
        int mani, lay;
        int NofStations  = 8;
        double StationD  = 390;
        double manifoldD = 10.575;
        double layerD    = 9.35;
        for (int i = 0; i < NofStations; i++) {
        	int manifold = 0;
        	while (manifold < 2) {
        		if (manifold == 0) {mani = -1;}
            		else if (manifold == 1) {mani = 1;}
		    	int layer = 0;
		    	while (layer < 2) {
                		if (layer == 0) {lay = -1;}
                		else if (layer == 1) {lay = 1;}
        			TLine *l = new TLine(-1365 + i * StationD + mani * manifoldD + lay*layerD*0.5, -695, -1365 + i * StationD + mani * manifoldD + lay*layerD*0.5, 695);
        		l->Draw();
        			//std::cout << -1365 + i * StationD + mani * manifoldD + lay*layerD*0.5 << std::endl;
        			layer++;
        		} // layer in manifold
        		manifold++;
        	} // manifold in station
        } // station
        
        c0->cd(2);
        g1->Draw("A*");
}
//sed -i '/nan/d' propagation_steps.txt

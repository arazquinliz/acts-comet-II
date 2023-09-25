void mapSeeds() {
    auto f = TFile::Open("seedMap.root");
    auto tree = (TTree*)f->Get("seeds");

    std::vector<ROOT::Math::XYZVector> *originalPos = 0;
    std::vector<ROOT::Math::XYZVector> *newPos = 0;
    std::vector<ROOT::Math::XYZVector> *comformalPos = 0;
    std::vector<int> *layer = 0;
    int event;

    tree->SetBranchAddress("originalPos", &originalPos);
    tree->SetBranchAddress("newPos", &newPos);
    tree->SetBranchAddress("comformalPos", &comformalPos);
    tree->SetBranchAddress("layer", &layer);
    tree->SetBranchAddress("event", &event);

    // ORIGINAL POSITION
    TCanvas *c1 = new TCanvas();

    TGraph *g[tree->GetEntries()];
    TMultiGraph *mg = new TMultiGraph();
    mg->SetTitle("Seeded SpacePoints; x; y");
    
    int i = 1;
    int prev_ev = -1;
    for (int ev = 0; ev < tree->GetEntries(); ev++) {
        tree->GetEntry(ev);
        if (event != prev_ev) {i = 1; prev_ev = event;}
        g[ev] = new TGraph();
        Double_t x[3] = {(*originalPos)[0].X(), (*originalPos)[1].X(), (*originalPos)[2].X()};
        Double_t y[3] = {(*originalPos)[0].Y(), (*originalPos)[1].Y(), (*originalPos)[2].Y()};
        
        /*for (auto point : *originalPos) {
            std::cout << point.X() << " " << point.Y() << " " << point.Z() << std::endl;
        }
        for (auto lay : *layer) {
            std::cout << lay << std::endl;
        }*/
        g[ev] = new TGraph(3, x, y);
        g[ev]->SetLineColor(i);
        //g[ev]->SetMarkerColor(event);
        g[ev]->SetMarkerSize(4);
        g[ev]->SetLineWidth(3);
        mg->Add(g[ev]);
        i++;
    }
    prev_ev = -1;
    
    mg->Draw("AC*");
    // Change the axis limits
    /*gPad->Modified();
    mg->GetXaxis()->SetLimits(-20, 30);
    mg->SetMinimum(75);
    mg->SetMaximum(125);*/
    c1->SetGrid();

    // NEW FRAME
    TCanvas *c2 = new TCanvas();

    TGraph *g2[tree->GetEntries()];
    TMultiGraph *mg2 = new TMultiGraph();
    mg2->SetTitle("Seeded SpacePoints in the new reference frame; x'; y'");
    
    i = 1;
    for (int ev = 0; ev < tree->GetEntries(); ev++) {
        tree->GetEntry(ev);
        if (event != prev_ev) {i = 1; prev_ev = event;}
        g2[ev] = new TGraph();
        Double_t x[3] = {(*newPos)[0].X(), (*newPos)[1].X(), (*newPos)[2].X()};
        Double_t y[3] = {(*newPos)[0].Y(), (*newPos)[1].Y(), (*newPos)[2].Y()};
        
        g2[ev] = new TGraph(3, x, y);
        g2[ev]->SetLineColor(i);
        //g2[ev]->SetMarkerColor(event);
        g2[ev]->SetLineWidth(3);
        mg2->Add(g2[ev]);
        i++;
    }
    prev_ev = -1;
    
    mg2->Draw("AC*");
    // Change the axis limits
    /*gPad->Modified();
    mg2->GetXaxis()->SetLimits(-1, 55);
    mg2->SetMinimum(-4);
    mg2->SetMaximum(1);*/
    c2->SetGrid();

    // COMFORMAL POSITIONS
    TCanvas *c3 = new TCanvas();

    TGraph *g3[tree->GetEntries()];
    TMultiGraph *mg3 = new TMultiGraph();
    mg3->SetTitle("Seeded SpacePoints in comformal space; u; v");
    
    i = 1;
    for (int ev = 0; ev < tree->GetEntries(); ev++) {
        tree->GetEntry(ev);
        if (event != prev_ev) {i = 1; prev_ev = event;}
        g3[ev] = new TGraph();
        Double_t x[2] = {(*comformalPos)[1].X(), (*comformalPos)[2].X()};
        Double_t y[2] = {(*comformalPos)[1].Y(), (*comformalPos)[2].Y()};
        
        g3[ev] = new TGraph(2, x, y);
        g3[ev]->SetLineColor(i);
        //g3[ev]->SetMarkerColor(event);
        g3[ev]->SetLineWidth(3);
        mg3->Add(g3[ev]);
        i++;
    }
    
    mg3->Draw("AL*");
    // Change the axis limits
    /*gPad->Modified();
    mg3->GetXaxis()->SetLimits(0, 0.08);
    mg3->SetMinimum(-0.0015);
    mg3->SetMaximum(0.0002);*/
    c3->SetGrid();
}
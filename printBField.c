// Visualization of magnetic field from printBField.c
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TProfile3D.h"
#include "TROOT.h"
#include "TTree.h"

void printBField(){
	// Parameters
	std::string outFile = "magnetic_plots.root"; // This file should be made following BFieldExample.cpp->RootBFieldWriter.cpp
	
	float xmin = -0.88, xmax = 0.9;
	float ymin = -0.9, ymax = 0.9;
	float zmin = -3.3345, zmax = 1.6655;
	int nBinsx = 90; // 91
	int nBinsy = 91; // 251
	int nBinsz = 227; // 251 90

	// Visualization settings
	const Int_t NRGBs        = 5;
	const Int_t NCont        = 225;
	Double_t    stops[NRGBs] = {0.00, 0.34, 0.61, 0.84, 1.00};
	Double_t    red[NRGBs]   = {0.00, 0.00, 0.87, 1.00, 0.51};
	Double_t    green[NRGBs] = {0.00, 0.81, 1.00, 0.20, 0.00};
	Double_t    blue[NRGBs]  = {0.51, 1.00, 0.12, 0.00, 0.00};
	TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
	gStyle->SetNumberContours(NCont);
	gStyle->SetOptStat(0);
	
	// Input information
	TFile* inputFile = TFile::Open("magnetic_field.root");
	if (inputFile == nullptr) {
		throw std::runtime_error("file does not exist");
	}
	
	TTree* tree = inputFile->Get<TTree>("bField"); // Name of the TTree in the ROOT file
	if (tree == nullptr) {
		throw std::runtime_error("tree object not found in file");
	}
	Int_t entries = tree->GetEntries();
	
	double x = 0, y = 0, z = 0;
	double Bx = 0, By = 0, Bz = 0;
	
	tree->SetBranchAddress("x", &x);
	tree->SetBranchAddress("y", &y);
	tree->SetBranchAddress("z", &z);

	tree->SetBranchAddress("Bx", &Bx);
	tree->SetBranchAddress("By", &By);
	tree->SetBranchAddress("Bz", &Bz);
	
	// Output information 
	std::cout << "Creating new output file: " << outFile << " and writing out histograms. " << std::endl;
	TFile outputFile(outFile.c_str(), "recreate");

	TProfile2D* bField_xy = new TProfile2D("BField_xy", "Magnetic Field", nBinsx, xmin, xmax, nBinsy, ymin, ymax);
	bField_xy->GetXaxis()->SetTitle("x [m]");
	bField_xy->GetYaxis()->SetTitle("y [m]");
	
	TProfile2D* bField_yz = new TProfile2D("BField_yz", "Magnetic Field", nBinsz, zmin, zmax, nBinsy, ymin, ymax);
	bField_yz->GetXaxis()->SetTitle("z [m]");
	bField_yz->GetYaxis()->SetTitle("y [m]");
	
	TProfile2D* bField_xz = new TProfile2D("BField_xz", "Magnetic Field", nBinsz, zmin, zmax, nBinsx, xmin, xmax);
	bField_xz->GetXaxis()->SetTitle("z [m]");
	bField_xz->GetYaxis()->SetTitle("x [m]");
	
	// This one doesn't work for some reason
	TProfile3D* bField_xyz = new TProfile3D("BField_xyz", "Magnetic Field", nBinsz, zmin, zmax, nBinsx, xmin, xmax, nBinsy, ymin, ymax);
	bField_xyz->GetXaxis()->SetTitle("z [m]");
	bField_xyz->GetYaxis()->SetTitle("x [m]");
	bField_xyz->GetZaxis()->SetTitle("y [m]");
 
	for (int i = 0; i < entries; i++) {
		tree->GetEvent(i);
		float bFieldValue = sqrt(Bx * Bx + By * By + Bz * Bz);

		bField_xy->Fill(x / 1000., y / 1000., bFieldValue);
		bField_yz->Fill(z / 1000., y / 1000., bFieldValue);
		bField_xz->Fill(z / 1000., x / 1000., bFieldValue);
		bField_xyz->Fill(z / 1000., x / 1000., y / 1000., bFieldValue);
	}
	bField_xy->Write();
	bField_yz->Write();
	bField_xz->Write();
	bField_xyz->Write();
	
	delete bField_xy;
	delete bField_yz;
	delete bField_xz;
	delete bField_xyz;
 
 	outputFile.Close();
 }

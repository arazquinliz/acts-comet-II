// ---------------------------------------------------------------- //
//                Creator: Amaia Razquin Lizarraga                  //
//                Contact: amaiafisika@gmail.com                    //
//                Date:    7 September 2023                         //
// ---------------------------------------------------------------- //
/* This module gets the magnetic field vector in each position from a 
   fieldmap root file (in get_mf()). The resulting map is printed out 
   via writeout_mf() as magnetic_field.root, then it can be visualized
   with printBField.c. The magnetic field is introduced in the provider
   in aGetMagneticField, the main function. */

#include <iostream>

// Definitions & utilities
#include <Acts/Definitions/Units.hpp>

// Root
#include <TGeoManager.h>
#include <TFile.h>
#include <TTree.h>

#include "MagField.hpp"

using namespace Acts::UnitLiterals;
using InterpolatedMagneticField3 =
    Acts::InterpolatedBFieldMap<Acts::detail::Grid<
        Acts::Vector3, Acts::detail::EquidistantAxis,
        Acts::detail::EquidistantAxis, Acts::detail::EquidistantAxis>>;     
        
// Get interpolated fieldmap from file
InterpolatedMagneticField3 get_mf(
	const std::function<size_t(std::array<size_t, 3> binsXYZ,
    std::array<size_t, 3> nBinsXYZ)>& localToGlobalBin) {
      	
    Acts::ActsScalar lengthUnit = 1_cm; 
    Acts::ActsScalar BFieldUnit = 0.1_T; 
    bool firstOctant = false; // Is field only described for the first octant?
    const std::string fieldMapFile = "fieldmap.root";
    const std::string treeName     = "field_description";
      	
    /// [1] Read in field map file
	// Grid position points in x, y and z
	std::vector<double> xPos;
	std::vector<double> yPos;
	std::vector<double> zPos;
	// components of magnetic field on grid points
	std::vector<Acts::Vector3> bField;
	// [1] Read in file and fill values
	TFile* inputFile = TFile::Open(fieldMapFile.c_str());
	if (inputFile == nullptr) {
		throw std::runtime_error("file does not exist");
	}
	TTree* tree = inputFile->Get<TTree>(treeName.c_str());
	if (tree == nullptr) {
		throw std::runtime_error("object not found in file");
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

	// reserve size
	xPos.reserve(entries);
	yPos.reserve(entries);
	zPos.reserve(entries);
	bField.reserve(entries);

	for (int i = 0; i < entries; i++) {
		tree->GetEvent(i);
		// Transport to Detector Solenoid frame
	    xPos.push_back(x - 1080);
	    yPos.push_back(y);
	   	zPos.push_back(z - 333.45);
	    bField.push_back(Acts::Vector3(Bx, By, Bz));
	}
	delete inputFile;
	
	// Setup the magnetic field
	auto BFieldmap = Acts::fieldMapXYZ(localToGlobalBin, xPos, yPos, zPos, bField, lengthUnit, 
		BFieldUnit, firstOctant); // map is InterpolatedBFieldMap
  	
	return BFieldmap;
}

// Write the interpolated magnetic field in a root file to plot
void writeout_mf(InterpolatedMagneticField3 BFieldmap) {
  	// Write out magnetic field (code from RootBFieldWriter.cpp)
 	TFile* outputFile = TFile::Open("magnetic_field.root", "recreate");
	if (outputFile == nullptr) {
		throw std::runtime_error("file does not exist");
	}
	
	TTree* outputTree = new TTree("bField", "bField"); // Name of the TTree in the ROOT file
	
	// Access the minima and maxima of all axes
  	auto minima = BFieldmap.getMin();
  	auto maxima = BFieldmap.getMax();
  	auto nBins  = BFieldmap.getNBins();
  	
  	// Write out the interpolated magnetic field map
    double minX = 0., minY = 0., minZ = 0.;
    double maxX = 0., maxY = 0., maxZ = 0.;
    size_t nBinsX = 0, nBinsY = 0, nBinsZ = 0;

    // The position values in xyz
    double x = 0;
    outputTree->Branch("x", &x);
    double y = 0;
    outputTree->Branch("y", &y);
    double z = 0;
    outputTree->Branch("z", &z);

    // The BField values in xyz
    double Bx = 0;
    outputTree->Branch("Bx", &Bx);
    double By = 0;
    outputTree->Branch("By", &By);
    double Bz = 0;
    outputTree->Branch("Bz", &Bz);
 	
 	// check dimension of Bfieldmap
    if (minima.size() == 3 && maxima.size() == 3) {
        minX = minima.at(0);
        minY = minima.at(1);
        minZ = minima.at(2);

        maxX = maxima.at(0);
        maxY = maxima.at(1);
        maxZ = maxima.at(2);

        nBinsX = nBins.at(0);
        nBinsY = nBins.at(1);
        nBinsZ = nBins.at(2);

    } 
    else if (minima.size() == 2 && maxima.size() == 2) {
        minX = -maxima.at(0);
        minY = -maxima.at(0);
        minZ = minima.at(1);

        maxX = maxima.at(0);
        maxY = maxima.at(0);
        maxZ = maxima.at(1);

        nBinsX = nBins.at(0);
        nBinsY = nBins.at(0);
        nBinsZ = nBins.at(1);
    } 
    else {
        throw std::invalid_argument(
			"BField has wrong dimension. The dimension needs to be "
            "either 2 (r,z,Br,Bz) or 3(x,y,z,Bx,By,Bz) in order to be "
            "written out by this writer.");
    }
 	
	assert(maxX > minX);
    assert(maxY > minY);
    assert(maxZ > minZ);

    double stepX = (maxX - minX) / (nBinsX - 1);
    double stepY = (maxY - minY) / (nBinsY - 1);
    double stepZ = (maxZ - minZ) / (nBinsZ - 1);

    for (size_t i = 0; i < nBinsX; i++) {
    	double raw_x = minX + i * stepX;
  		for (size_t j = 0; j < nBinsY; j++) {
        	double raw_y = minY + j * stepY;
        	for (size_t k = 0; k < nBinsZ; k++) {
      			double raw_z = minZ + k * stepZ;
      			Acts::Vector3 position(raw_x, raw_y, raw_z);
          		Acts::Vector3 mbField = BFieldmap.getFieldUnchecked(position);

          		x = raw_x / Acts::UnitConstants::mm;
          		y = raw_y / Acts::UnitConstants::mm;
          		z = raw_z / Acts::UnitConstants::mm;
          		Bx = mbField.x() / Acts::UnitConstants::T;
          		By = mbField.y() / Acts::UnitConstants::T;
          		Bz = mbField.z() / Acts::UnitConstants::T;
          			
          		outputTree->Fill();
        	}  // for z
      	}    // for y
    }      // for x
    outputTree->Write();
	delete outputFile; 
}

// Main function: get field and write it out
std::shared_ptr<Acts::MagneticFieldProvider> aGetMagneticField(bool debug,
	Acts::MagneticFieldContext mctx,
	const std::function<size_t(std::array<size_t, 3> binsXYZ,
    std::array<size_t, 3> nBinsXYZ)>& localToGlobalBin) {

	auto BFieldmap = get_mf(localToGlobalBin);
	
	auto cache = BFieldmap.makeCache(mctx);
	
	// Access the minima and maxima of all axes
  	auto minima = BFieldmap.getMin();
  	auto maxima = BFieldmap.getMax();
  	auto nBins  = BFieldmap.getNBins();
  	
  	if (debug) {
		std::cout << "Minimum:";
    		for (auto m : minima) {
      			std::cout << " " << m;
    		}
    		std::cout << " Maximum:";
    		for (auto m : maxima) {
      			std::cout << " " << m;
    		}
    		std::cout << " nBins:";
    		for (auto m : nBins) {
      			std::cout << " " << m;
    		}
    		std::cout << std::endl;
  	}
  	
  	writeout_mf(BFieldmap);
  	
  	return std::make_shared<InterpolatedMagneticField3>(std::move(BFieldmap));
}
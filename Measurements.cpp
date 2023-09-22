// ---------------------------------------------------------------- //
//                Creator: Amaia Razquin Lizarraga                  //
//                Contact: amaiafisika@gmail.com                    //
//                Date:    7 September 2023                         //
// ---------------------------------------------------------------- //
/* This module gets the measurements from premade measurement files
   the geometry layout has to be created and introduce from the main
   program. The measurements are extracted from the file for straw and
   disk geometry in different functions. In the straw case, the measurements 
   can be 2D or 1D (New). A standard covariance generator function is also included */

#include <iostream>
#include <ostream>
#include <vector>

// ACTS
#include <Acts/Definitions/Units.hpp>

// ROOT
#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>

// AMAIA
#include "Measurements.hpp"
#include "base/Calibrator.hpp"

using namespace Acts::UnitLiterals;

Acts::ActsSquareMatrix<2> aGetDiskCovariance(const Acts::GeometryContext& gctx, 
	std::shared_ptr<const Acts::Surface> surface, 
	const Acts::Vector2 locCov, double posZ, int k) {

	double radius    = 0.49_cm;

	Acts::ActsSquareMatrix<2> cova = Acts::ActsSquareMatrix<2>::Zero();
	cova(Acts::eBoundLoc0, Acts::eBoundLoc1) = 0._mm;
	cova(Acts::eBoundLoc1, Acts::eBoundLoc0) = 0._mm;

	Acts::Vector2 uniCov = surface->globalToLocal(gctx, 
	Acts::Vector3{surface->center(gctx).x() - 0.15, //abs((1.0/12) * std::pow(2*radius, 2)
		surface->center(gctx).y() - 0.15, posZ},
		Acts::Vector3::Zero()).value(); 
	cova(Acts::eBoundLoc0, Acts::eBoundLoc0) = abs(uniCov.x());
	cova(Acts::eBoundLoc1, Acts::eBoundLoc1) = abs(uniCov.y());

	return cova;
}

// Makes a layout of the hits by classifying them per station
std::map<std::string, std::vector<Acts::Vector3>> aReadDiskLayout(
	std::vector<ROOT::Math::XYZVector> spOfEvent,
	std::shared_ptr<const Acts::TrackingGeometry> tGeometry, 
	const Acts::GeometryContext& gctx, bool debug) {
      	
	double hitX, hitY, hitZ, hitT;
	double StationD    = 39.0_cm;    // Distance between stations
	double NofStations = 8;          // Number of stations
	double thickness   = 20.15_mm;   // Half thickness of station (position+radius+clearance)     

    std::map<std::string, std::vector<Acts::Vector3>> hit_layout; // link hits to layers
	double prev = 0; // Helper to make sure that the measurements are in order
    for (size_t it = 0; it != spOfEvent.size(); ++it) {
		if (prev > spOfEvent[it].Z()) {continue;}
		else {prev = spOfEvent[it].Z();}
		hitX = spOfEvent[it].X(); // in mm
		hitY = spOfEvent[it].Y();
		hitZ = spOfEvent[it].Z();
      		Acts::Vector3 globalPos(hitX - 1080._cm, hitY, 
      			hitZ - 333.45_cm); // translated
      		int sttn; // Organize hits per station	
      		const auto laynum = 
      			tGeometry->associatedLayer(gctx, globalPos)->geometryId().layer();
	      	if (laynum % 2 == 0) {
	      		sttn = laynum / 2 - 1;
	      	}
	      	else if (laynum % 2 != 0) {
	      		sttn = (laynum + 1) / 2 - 1;
	      	}

		int stat;
		for (int i = 0; i < NofStations; i++) {
			double tVol_z = -136.5_cm + i * StationD;
			if ((globalPos.z() <= tVol_z + thickness) && (globalPos.z() >= tVol_z - thickness)) {
				stat = i;
			}
		}
      		
	    auto station = "Station_" + std::to_string(stat);
	    hit_layout[station].push_back(globalPos);
		
	    if (debug){
	      	std::cout << station << " x = " << 
 	      		globalPos.x() << " y = " << 
	      		globalPos.y() << " z = " << 
	      		globalPos.z() << std::endl;
	    }
    }
	return hit_layout;
}

// Makes a layout of the hits by classifying them per station
// For straw and 2D measurements
std::map<std::string, std::vector<std::pair<Acts::Vector3, Acts::Vector3>>> aReadStrawLayout(
	std::vector<ROOT::Math::XYZVector> spOfEvent,
	std::vector<ROOT::Math::XYZVector> momOfEvent,
	std::shared_ptr<const Acts::TrackingGeometry> tGeometry, 
	const Acts::GeometryContext& gctx, bool debug) {
    
	std::cout << "Reading straw hit layout.\n";
	double hitX, hitY, hitZ, hitT;
	double StationD    = 39.0_cm;    // Distance between stations
	double NofStations = 8;          // Number of stations
	double thickness   = 20.15_mm;   // Half thickness of station (position+radius+clearance)     

    std::map<std::string, std::vector<std::pair<Acts::Vector3, Acts::Vector3>>> hit_layout; // link hits to layers
	double prev = 0; // Helper to make sure that the measurements are in order
    for (size_t it = 0; it != spOfEvent.size(); ++it) {
		if (prev > spOfEvent[it].Z()) {continue;}
		else {prev = spOfEvent[it].Z();}
		hitX = spOfEvent[it].X(); // in mm
		hitY = spOfEvent[it].Y();
		hitZ = spOfEvent[it].Z();
      	
		Acts::Vector3 globalPos(hitX - 1080._cm, hitY, 
      			hitZ - 333.45_cm); // translated
		Acts::Vector3 globalMom(momOfEvent[it].X() * 0.001, momOfEvent[it].Y() * 0.001,
			momOfEvent[it].Z() * 0.001);
			
		// Get station
		int stat;
		for (int i = 0; i < NofStations; i++) {
			double tVol_z = -136.5_cm + i * StationD;
			if ((globalPos.z() <= tVol_z + thickness) && (globalPos.z() >= tVol_z - thickness)) {
				stat = i;
			}
		}
      	
		// Get layer 
		int laynum;
		double centerPos = -136.5_cm + stat * StationD;
		if (globalPos.z() < centerPos) {
			if (globalPos.z() < centerPos - 10.075) {laynum = 2;}
			else if (globalPos.z() > centerPos - 10.075) {laynum = 4;}
		}
		else if (globalPos.z() >= centerPos) {
			if (globalPos.z() < centerPos + 10.075) {laynum = 6;}
			else if (globalPos.z() >= centerPos + 10.075) {laynum = 8;}
		}

	    auto station_layer = "Station_" + std::to_string(stat) + "_" + std::to_string(laynum);
	    hit_layout[station_layer].push_back(std::make_pair(globalPos, globalMom));
		bool yes = false;
	    if (yes){
	      	std::cout << station_layer << " x = " << 
 	      		globalPos.x() << " y = " << 
	      		globalPos.y() << " z = " << 
	      		globalPos.z() << std::endl;
	    }
    }
	return hit_layout;
}

// Makes a layout of the hits by classifying them per station
// For straw and 1D measurements
std::map<std::string, std::vector<std::tuple<Acts::Vector3, Acts::Vector3, double>>> aReadStrawLayoutNew(
	std::vector<ROOT::Math::XYZVector> spOfEvent,
	std::vector<ROOT::Math::XYZVector> momOfEvent,
	std::vector<double> driftD,
	std::shared_ptr<const Acts::TrackingGeometry> tGeometry, 
	const Acts::GeometryContext& gctx, bool debug) {
    
	std::cout << "Reading straw hit layout.\n";
	double hitX, hitY, hitZ, hitT;
	double StationD    = 39.0_cm;    // Distance between stations
	double NofStations = 8;          // Number of stations
	double thickness   = 20.15_mm;   // Half thickness of station (position+radius+clearance)     
    std::map<std::string, std::vector<std::tuple<Acts::Vector3, Acts::Vector3, double>>> hit_layout; // link hits to layers
	double prev = 0; // Helper to make sure that the measurements are in order
    for (size_t it = 0; it != spOfEvent.size(); ++it) {
		if (prev > spOfEvent[it].Z()) {continue;}
		else {prev = spOfEvent[it].Z();}
		hitX = spOfEvent[it].X(); // in mm
		hitY = spOfEvent[it].Y();
		hitZ = spOfEvent[it].Z();
      	
		Acts::Vector3 globalPos(hitX - 1080._cm, hitY, 
      			hitZ - 333.45_cm); // translated
		Acts::Vector3 globalMom(momOfEvent[it].X() * 0.001, momOfEvent[it].Y() * 0.001,
			momOfEvent[it].Z() * 0.001);
			
		// Get station
		int stat;
		for (int i = 0; i < NofStations; i++) {
			double tVol_z = -136.5_cm + i * StationD;
			if ((globalPos.z() <= tVol_z + thickness) && (globalPos.z() >= tVol_z - thickness)) {
				stat = i;
			}
		}
      	
		// Get layer 
		int laynum;
		double centerPos = -136.5_cm + stat * StationD;
		if (globalPos.z() < centerPos) {
			if (globalPos.z() < centerPos - 10.075) {laynum = 2;}
			else if (globalPos.z() > centerPos - 10.075) {laynum = 4;}
		}
		else if (globalPos.z() >= centerPos) {
			if (globalPos.z() < centerPos + 10.075) {laynum = 6;}
			else if (globalPos.z() >= centerPos + 10.075) {laynum = 8;}
		}

	    auto station_layer = "Station_" + std::to_string(stat) + "_" + std::to_string(laynum);
	    hit_layout[station_layer].push_back(std::make_tuple(globalPos, globalMom, driftD[it]));
		bool yes = false;
	    if (yes){
	      	std::cout << station_layer << " x = " << 
 	      		globalPos.x() << " y = " << 
	      		globalPos.y() << " z = " << 
	      		globalPos.z() << std::endl;
	    }
    }
	return hit_layout;
}

// Reads through the events to create sourcelinks and measurements
// Straw measurements in 2D
std::vector<result> aCreateStrawLinks(
	std::map<std::string, std::vector<std::shared_ptr<const Acts::Surface>>> tracker_layout,
	std::shared_ptr<const Acts::TrackingGeometry> tGeometry,
	const Acts::GeometryContext& gctx, bool debug) {

	int NofStations  = 8;            // Number of stations
	double radius    = 0.49_cm;      // StrawTubeOuterRadius

	// get the tree and iterate over the branches
	auto *f = TFile::Open("detector_oa_all_t.root");
	if (!f) {
		std::cout << "Couldn't open measurement file\n";
	}

	auto *tree = (TTree*)f->Get("detections");
	TVector3 *initialPos = 0; 
	TVector3 *initialMom = 0;
	std::vector<ROOT::Math::XYZVector> *spOfEvent = 0;
	std::vector<ROOT::Math::XYZVector> *momOfEvent = 0;
	tree->SetBranchAddress("spOfEvent", &spOfEvent);
	tree->SetBranchAddress("momOfEvent", &momOfEvent);
	tree->SetBranchAddress("initialPos", &initialPos);	
	tree->SetBranchAddress("initialMom", &initialMom);

	std::vector<result> result_vector;

	std::fstream linkedfile;
	linkedfile.open("linked_meas.txt", std::ios::out);
	for (size_t ev = 0; ev < 20000; ++ev) {
	//for (size_t ev = 0; ev < tree->GetEntries(); ++ev) {
		
		std::vector<Acts::SourceLink> SourceLinkVec;
		std::vector<Acts::BoundVariantMeasurement> new_measurements;
		double ix, iy, iz;

		tree->GetEntry(ev);
		std::map<std::string, std::vector<std::pair<Acts::Vector3, Acts::Vector3>>> hit_layout = 
			aReadStrawLayout(*spOfEvent, *momOfEvent, tGeometry, gctx, debug);

	    // Average the hits around each surface to get a single point in the surface
		std::string station; 
	    for (int i = 0; i < NofStations; ++i) {
			for (int j = 1; j < 5; ++j) {
				bool found = false;
	      		auto station = "Station_" + std::to_string(i) + "_" + std::to_string(2 * j);
				for (auto& surface : tracker_layout[station]) {
					if (found) {continue;}
					double posZ = 0, posY = 0, posX = 0;
					double momX = 0, momY = 0, momZ = 0;
					int k = 0;
					// Value vector for covariance 
					std::vector<Acts::Vector3> chi;
					chi.clear();
					for (auto hit : hit_layout[station]) {
						Acts::Vector3 surf_center = surface->center(gctx);
						if ((j == 1 || j == 2) && (abs(hit.first.x() - surf_center.x()) <= radius)) {
							posX += hit.first.x();
							posY += hit.first.y();
							posZ  = surf_center.z();
							k    += 1; 
							momX += hit.second.x();
							momY += hit.second.y();
							momZ += hit.second.z();
							chi.push_back(Acts::Vector3{hit.first.x(), hit.first.y(), hit.first.z()});
							found = true;
						}
						else if ((j == 3 || j == 4) && (abs(hit.first.y() - surf_center.y()) <= radius)) {
							posX += hit.first.x();
							posY += hit.first.y();
							posZ  = surf_center.z();
							k    += 1; 
							chi.push_back(Acts::Vector3{hit.first.x(), hit.first.y(), hit.first.z()});
							found = true;
						}
					}
						
					if (posY == 0 & posX == 0) {continue;} // do not add sp without hits
					
					int kx, ky;
					if (j == 1 || j == 2) {kx = k; ky = k;}
					else if (j == 3 || j == 4) {kx = k; ky = k;}
					Acts::Vector2 loc = surface->globalToLocal(gctx, 
						Acts::Vector3( posX/kx, posY/ky, posZ ),
						Acts::Vector3{momX/k, momY/k, momY/k}).value();
						//Acts::Vector3::Zero()).value();
						
					auto geoID        = surface->geometryId();
					
					Acts::Vector3 glob = surface->localToGlobal(gctx, loc, Acts::Vector3::Zero());
					linkedfile << posX/kx << " " << posY/ky << " " << posZ << "\n";
					// Caclulate covariace explicitly
					Acts::Vector3 cov = Acts::Vector3::Zero();  
					for (int j = 0; j < k; ++j) {
						cov[0] += std::pow(chi[j].x() - posX/kx, 2);  
						cov[1] += std::pow(chi[j].y() - posY/ky, 2);
					}
					cov[0] = surface->center(gctx).x() - cov[0];
					cov[1] = surface->center(gctx).y() - cov[1];
					cov[2] = posZ;

					Acts::Vector2 locCov = surface->globalToLocal(gctx, cov,
						Acts::Vector3::Zero()).value(); 
						
					bool yes = true;
					if (yes){
						//bool isit = surface->Acts::Surface::isOnSurface(*gctx,
							//Acts::Vector3( posX, posY/k, posZ/k ),
							//Acts::Vector3::Zero()); 
						std::cout << geoID << " x = "   << 
							posX/k  << " y = "    << 
							posY/k  << " z = "    << 
							posZ    << " loc0 = " <<
							loc.x() << " loc1 = " << 
							loc.y() << std::endl;
					}
					
					// NEW
					std::array<Acts::BoundIndices, 2> indices;
					indices[0] = Acts::BoundIndices::eBoundLoc0;
					indices[1] = Acts::BoundIndices::eBoundLoc1;

					Acts::ActsVector<2> loca;
					loca[Acts::eBoundLoc0] = loc.x();
					loca[Acts::eBoundLoc1] = loc.y();

					Acts::ActsSquareMatrix<2> cova = Acts::ActsSquareMatrix<2>::Zero();
					cova = aGetDiskCovariance(gctx, surface, locCov, posZ, k);
					
					ActsSourceLink::Index index = new_measurements.size();
					ActsSourceLink sl(surface->geometryId(), index);
					Acts::Measurement<Acts::BoundIndices, 2> meas(Acts::SourceLink{sl}, indices, loca, cova);
					SourceLinkVec.push_back(Acts::SourceLink{sl});
					new_measurements.push_back(meas);
				} // surface
			} // layer
		} // station
		std::cout << "Size of measurements is: " << new_measurements.size() << std::endl;

		auto iPos = Acts::Vector3{initialPos->X() - 1080._cm, initialPos->Y(), 
			initialPos->Z() - 333.45_cm};
		auto iMom = Acts::Vector3{initialMom->X(), initialMom->Y(), initialMom->Z()};
		result_vector.push_back(result{new_measurements, SourceLinkVec, iPos, iMom});
	} // event loop

	linkedfile.close();

	return result_vector;
}

// Reads through the events to create sourcelinks and measurements
// Straw measurements in 1D
std::vector<result> aCreateStrawLinksNew(
	std::map<std::string, std::vector<std::shared_ptr<const Acts::Surface>>> tracker_layout,
	std::shared_ptr<const Acts::TrackingGeometry> tGeometry,
	const Acts::GeometryContext& gctx, bool debug) {

	int NofStations  = 8;            // Number of stations
	double radius    = 0.49_cm;      // StrawTubeOuterRadius

	// get the tree and iterate over the branches
	auto *f = TFile::Open("detector_oa_all_td.root"); // detector_oa_all.root
	if (!f) {
		std::cout << "Couldn't open measurement file\n";
	}

	auto *tree = (TTree*)f->Get("detections");
	TVector3 *initialPos = 0; 
	TVector3 *initialMom = 0;
	std::vector<ROOT::Math::XYZVector> *spOfEvent = 0;
	std::vector<ROOT::Math::XYZVector> *momOfEvent = 0;
	std::vector<double> *driftD = 0;
	tree->SetBranchAddress("spOfEvent", &spOfEvent);
	tree->SetBranchAddress("momOfEvent", &momOfEvent);
	tree->SetBranchAddress("initialPos", &initialPos);	
	tree->SetBranchAddress("initialMom", &initialMom);
	tree->SetBranchAddress("driftD", &driftD);

	std::vector<result> result_vector;

	std::fstream linkedfile;
	linkedfile.open("linked_meas.txt", std::ios::out);
	for (size_t ev = 0; ev < 8000; ++ev) {
	//for (size_t ev = 0; ev < tree->GetEntries(); ++ev) {
		
		std::vector<Acts::SourceLink> SourceLinkVec;
		std::vector<Acts::BoundVariantMeasurement> new_measurements;
		double ix, iy, iz;

		tree->GetEntry(ev);
		std::map<std::string, std::vector<std::tuple<Acts::Vector3, Acts::Vector3, double>>> hit_layout = 
			aReadStrawLayoutNew(*spOfEvent, *momOfEvent, *driftD, tGeometry, gctx, debug);

	    // Average the hits around each surface to get a single point in the surface
		std::string station; 
	    for (int i = 0; i < NofStations; ++i) {
			for (int j = 1; j < 5; ++j) {
				bool found = false;
	      		auto station = "Station_" + std::to_string(i) + "_" + std::to_string(2 * j);
				for (auto& surface : tracker_layout[station]) {
					if (found) {continue;}
					double posZ = 0, posY = 0, posX = 0;
					double momX = 0, momY = 0, momZ = 0;
					int k = 0;
					double drift = 5._mm;
					// Value vector for covariance 
					std::vector<Acts::Vector3> chi;
					chi.clear();
					for (auto hit : hit_layout[station]) {
						Acts::Vector3 surf_center = surface->center(gctx);
						if ((j == 1 || j == 2) && (abs(std::get<0>(hit).x() - surf_center.x()) <= radius)) {
							posX += std::get<0>(hit).x();
							posY += std::get<0>(hit).y();
							posZ  = surf_center.z();
							k    += 1; 
							momX += std::get<1>(hit).x();
							momY += std::get<1>(hit).y();
							momZ += std::get<1>(hit).z();
							chi.push_back(Acts::Vector3{std::get<0>(hit).x(), std::get<0>(hit).y(), std::get<0>(hit).z()});
							found = true;
						}
						else if ((j == 3 || j == 4) && (abs(std::get<0>(hit).y() - surf_center.y()) <= radius)) {
							posX += std::get<0>(hit).x();
							posY += std::get<0>(hit).y();
							posZ  = std::get<0>(hit).z();
							k    += 1; 
							chi.push_back(Acts::Vector3{std::get<0>(hit).x(), std::get<0>(hit).y(), std::get<0>(hit).z()});
							found = true;
						}
						if (std::get<2>(hit) < drift) {drift = std::get<2>(hit);}
					}
						
					if (posY == 0 & posX == 0) {continue;} // do not add sp without hits
					
					int kx, ky;
					if (j == 1 || j == 2) {kx = k; ky = k;}
					else if (j == 3 || j == 4) {kx = k; ky = k;}
					Acts::Vector2 loc = surface->globalToLocal(gctx, 
						Acts::Vector3( posX/kx, posY/ky, posZ ),
						Acts::Vector3{momX/k, momY/k, momY/k}).value();
						//Acts::Vector3::Zero()).value();
						
					auto geoID        = surface->geometryId();
					
					Acts::Vector3 glob = surface->localToGlobal(gctx, loc, Acts::Vector3::Zero());
					linkedfile << posX/kx << " " << posY/ky << " " << posZ << "\n";
					// Caclulate covariace explicitly
					Acts::Vector3 cov = Acts::Vector3::Zero();  
					for (int j = 0; j < k; ++j) {
						cov[0] += std::pow(chi[j].x() - posX/kx, 2);  
						cov[1] += std::pow(chi[j].y() - posY/ky, 2);
					}
					cov[0] = surface->center(gctx).x() - cov[0];
					cov[1] = surface->center(gctx).y() - cov[1];
					cov[2] = posZ;
					//std::cout << k << kx << ky << " " << cov[0] << " " << cov[1] << " " << cov[2] << std::endl;

					Acts::Vector2 locCov = surface->globalToLocal(gctx, cov,
						Acts::Vector3::Zero()).value(); 
					//std::cout << locCov.x() << " " << locCov.y() << std::endl;
					bool yes = true;
					if (yes){
						//bool isit = surface->Acts::Surface::isOnSurface(*gctx,
						//Acts::Vector3( posX, posY/k, posZ/k ),
						//Acts::Vector3::Zero()); 
						std::cout << geoID << " x = "   << 
							posX/k  << " y = "    << 
							posY/k  << " z = "    << 
							posZ    << " loc0 = " <<
							loc.x() << " loc1 = " << 
							loc.y() << std::endl;
					}

					// NEW
					std::array<Acts::BoundIndices, 1> indices = {Acts::eBoundLoc0};

					Acts::ActsVector<1> loca;
					loca[Acts::eBoundLoc0] = drift;

					double sigma = 0.15_mm;
					Acts::ActsSquareMatrix<1> cova = Acts::ActsSquareMatrix<1>::Identity();
					cova(0, 0) = sigma * sigma; // 8.003_mm; // 
					
					ActsSourceLink::Index index = new_measurements.size();
					ActsSourceLink sl(surface->geometryId(), index);
					Acts::Measurement<Acts::BoundIndices, 1> meas(Acts::SourceLink{sl}, indices, loca, cova);
					SourceLinkVec.push_back(Acts::SourceLink{sl});
					new_measurements.push_back(meas);
				} // surface
			} // layer
		} // station
		std::cout << "Size of measurements is: " << new_measurements.size() << std::endl;

		auto iPos = Acts::Vector3{initialPos->X() - 1080._cm, initialPos->Y(), 
			initialPos->Z() - 333.45_cm};
		auto iMom = Acts::Vector3{initialMom->X(), initialMom->Y(), initialMom->Z()};
		result_vector.push_back(result{new_measurements, SourceLinkVec, iPos, iMom});
	} // event loop

	linkedfile.close();

	return result_vector;
}

// Reads through the events to create sourcelinks and measurements
// Disk measurements in 2D
std::vector<result> aCreateDiskLinks(
	std::map<std::string, std::vector<std::shared_ptr<const Acts::Surface>>> tracker_layout,
	std::shared_ptr<const Acts::TrackingGeometry> tGeometry,
	const Acts::GeometryContext& gctx, bool debug) {

	int NofStations  = 8;            // Number of stations
	double radius    = 0.49_cm;      // StrawTubeOuterRadius

	// get the tree and iterate over the branches
	auto *f = TFile::Open("detector_oa_all_t.root"); //  propagated_meas.root
	if (!f) {
		std::cout << "Couldn't open measurement file\n";
	}

	auto *tree = (TTree*)f->Get("detections");
	TVector3 *initialPos = 0; 
	TVector3 *initialMom = 0;
	std::vector<ROOT::Math::XYZVector> *spOfEvent = 0;
	tree->SetBranchAddress("spOfEvent", &spOfEvent);
	tree->SetBranchAddress("initialPos", &initialPos);	
	tree->SetBranchAddress("initialMom", &initialMom);

	std::vector<result> result_vector;

	std::fstream linkedfile;
	linkedfile.open("linked_meas.txt", std::ios::out);
	for (size_t ev = 0; ev < 20000; ++ev) {
	//for (size_t ev = 0; ev < tree->GetEntries(); ++ev) {
		
		std::vector<Acts::SourceLink> SourceLinkVec;
		std::vector<Acts::BoundVariantMeasurement> new_measurements;
		double ix, iy, iz;

		tree->GetEntry(ev);
		std::map<std::string, std::vector<Acts::Vector3>> hit_layout = 
			aReadDiskLayout(*spOfEvent, tGeometry, gctx, debug);

	    // Average the hits around each surface to get a single point in the surface
	    for (int i = 0; i < NofStations; ++i) {
	      	auto station = "Station_" + std::to_string(i);
		    for (auto& surface : tracker_layout[station]) {
		      	Acts::Vector3 surf_center = surface->center(gctx);
			    double posZ = 0; 
			    double posY = 0, posX = 0;
			    int k = 0;
				// Value vector for covariance 
				std::vector<Acts::Vector3> chi;
				chi.clear();
			    for (auto hit : hit_layout[station]) {
			      	if (abs(hit.z() - surf_center.z()) <= radius) {
			      		posX += hit.x();
			      		posY += hit.y();
			    		posZ  = surf_center.z();
			  			k    += 1; 
						chi.push_back(Acts::Vector3{hit.x(), hit.y(), hit.z()});
		      		}	
		      	}
		      		
		      	if (posY == 0 & posX == 0) {continue;} // do not add sp without hits
			      	
		    	Acts::Vector2 loc = surface->globalToLocal(gctx, 
					Acts::Vector3( posX/k, posY/k, posZ ),
					Acts::Vector3::Zero()).value();
					
		  		auto geoID        = surface->geometryId();
				
				Acts::Vector3 glob = surface->localToGlobal(gctx, loc, Acts::Vector3::Zero());
				linkedfile << posX/k << " " << posY/k << " " << posZ << "\n";
				// Caclulate covariace explicitly
				Acts::Vector3 cov = Acts::Vector3::Zero();  
				for (int j = 0; j < k; ++j) {
					cov[0] += std::pow(chi[j].x() - posX/k, 2);  
					cov[1] += std::pow(chi[j].y() - posY/k, 2);
				}
				cov[2] = posZ;

				Acts::Vector2 locCov = surface->globalToLocal(gctx, cov,
					Acts::Vector3::Zero()).value(); 
			    if (debug){
					//bool isit = surface->Acts::Surface::isOnSurface(*gctx,
						//Acts::Vector3( posX, posY/k, posZ/k ),
						//Acts::Vector3::Zero()); 
			      	std::cout << geoID << " x = "   << 
			 	  		posX/k  << " y = "    << 
				  		posY/k  << " z = "    << 
				  		posZ    << " loc0 = " <<
				  		loc.x() << " loc1 = " << 
				  		loc.y() << std::endl;
			    }
				
				// NEW
				std::array<Acts::BoundIndices, 2> indices;
				indices[0] = Acts::BoundIndices::eBoundLoc0;
				indices[1] = Acts::BoundIndices::eBoundLoc1;

				Acts::ActsVector<2> loca;
				loca[Acts::eBoundLoc0] = loc.x();
				loca[Acts::eBoundLoc1] = loc.y();

				Acts::ActsSquareMatrix<2> cova = Acts::ActsSquareMatrix<2>::Zero();
				cova(Acts::eBoundLoc0, Acts::eBoundLoc0) = 0.15_mm; 
                cova(Acts::eBoundLoc1, Acts::eBoundLoc1) = 0.15_mm;
				
				ActsSourceLink::Index index = new_measurements.size();
				ActsSourceLink sl(surface->geometryId(), index);
				Acts::Measurement<Acts::BoundIndices, 2> meas(Acts::SourceLink{sl}, indices, loca, cova);
				SourceLinkVec.push_back(Acts::SourceLink{sl});
				new_measurements.push_back(meas);
		    	} // surface
		} // station
		std::cout << "Size of measurements is: " << new_measurements.size() << std::endl;

		auto iPos = Acts::Vector3{initialPos->X() - 1080._cm, initialPos->Y(), 
			initialPos->Z() - 333.45_cm};
		auto iMom = Acts::Vector3{initialMom->X(), initialMom->Y(), initialMom->Z()};
		result_vector.push_back(result{new_measurements, SourceLinkVec, iPos, iMom});
	} // event loop

	linkedfile.close();

	return result_vector;
}
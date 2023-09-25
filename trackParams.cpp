// ---------------------------------------------------------------- //
//                Creator: Amaia Razquin Lizarraga                  //
//                Contact: amaiafisika@gmail.com                    //
//                Date:    7 September 2023                         //
// ---------------------------------------------------------------- //
/* 
    This program uses the track parameter estimator with all the available
    triplets in a trajectory to get the momentum estimation. The estimations
    for every layer combination are plotted in the same histogram and outputed
    as pEstim_combi.pdf. This program is set to perform in a single station
    with the stationID variable. The results are then printed in estim_stat.txt,
    if done with every station (and by changing the name of estimFile) all the
    estimations can be compared using the estim_plot.c program in root.
    The modified estimator is also used and compared (from Estimator.cpp).

    Command:
        g++ seedSelection.cpp Geometry.cpp -o params.x
        -L/path/to/acts_install/lib -lActsPluginTGeo -lActsPluginJson -lActsCore 
        -I/path/to/acts_install/include/ -I$(root-config --incdir) 
        -I$(root-config --evelibs) $(root-config --libs) 
        -lGeom -I/usr/include/eigen3
*/

#include <fstream>
#include <iostream>
#include <vector>

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/Units.hpp>

#include <Acts/Seeding/BinnedSPGroup.hpp>
#include <Acts/Seeding/InternalSeed.hpp>
#include <Acts/Seeding/InternalSpacePoint.hpp>
#include <Acts/Seeding/Seed.hpp>
#include <Acts/Seeding/SeedFinder.hpp>
#include <Acts/Seeding/SeedFilter.hpp>
#include <Acts/Seeding/EstimateTrackParamsFromSeed.hpp>

// ROOT
#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <Math/Vector3D.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TRandom3.h>
#include <TH1.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TPaveLabel.h>

// AMAIA:
#include "base/SpacePoint.hpp"
#include "Geometry.hpp"
#include "MagField.hpp"
#include "Estimator.hpp"

using namespace Acts::UnitLiterals;

extern std::map<std::string, std::vector<std::shared_ptr<const Acts::Surface>>> tracker_layout;
std::map<int, std::pair<int, Acts::Vector3>> event_map; // event number, seeds
extern std::fstream estimatorFile;

std::unique_ptr<SpacePoint> makeSpacePoint(const Acts::GeometryContext& gctx,
    Acts::Vector3 globalPos, 
    double time, 
    std::shared_ptr<const Acts::Surface> surf,
    const int event, const Acts::Vector3 p) {

    bool debug = true;

    float x = globalPos.x() - 1080._cm;
    float y = globalPos.y();
    float z = globalPos.z() - 333.45_cm;
    float r = std::sqrt(x * x + y * y);

    Acts::Vector3 mom(1, 1, 1); 
    Acts::SquareMatrix2 localCov = Acts::SquareMatrix2::Zero();
    localCov(0, 0) = 1._mm;
    localCov(1, 1) = 1._mm;

    /// The space point requires only the variance of the transverse and
    /// longitudinal position. Reduce computations by transforming the
    /// covariance directly from local to r/z.
    ///
    /// compute Jacobian from global coordinates to r/z
    ///
    ///         r = sqrt(x² + y²)
    /// dr/d{x,y} = (1 / sqrt(x² + y²)) * 2 * {x,y}
    ///             = 2 * {x,y} / r
    ///       dz/dz = 1 
    Acts::RotationMatrix3 rotLocalToGlobal = surf->referenceFrame(gctx, 
        Acts::Vector3{x, y, z}, mom);
    auto scale = 2.0 / std::hypot(x, y);
    Acts::ActsMatrix<2, 3> jacXyzToRhoZ = Acts::ActsMatrix<2, 3>::Zero();
    jacXyzToRhoZ(0, Acts::ePos0) = scale * x;
    jacXyzToRhoZ(0, Acts::ePos1) = scale * y;
    jacXyzToRhoZ(0, Acts::ePos2) = 1;
    // Compute Jacobian from local coordinates to rho/z
    Acts::ActsMatrix<2, 2> jac = jacXyzToRhoZ * rotLocalToGlobal.block<3, 2>(
        Acts::ePos0, Acts::ePos0);
    // Compute rho/z variance
    Acts::ActsVector<2> var = (jac * localCov * jac.transpose()).diagonal();

    /*
    * From Acts v17 to v19 the scattering uncertainty value allowed was changed.
    * This led to a decrease in efficiency. To offset this, we scale the 
    * uncertainties by a tuned factor that gives the v17 performance
    * Track reconstruction is an art as much as it is a science...
    */
    double uncfactor = 3.175; // from sPHENIX
    std::unique_ptr<SpacePoint> spPtr(new SpacePoint{x, y, z, r, time, surf->geometryId(), 
        var[0] * uncfactor, var[1] * uncfactor, event, p});

    /*if (debug) {
        std::cout << "Space point has " << x << " " << y << " " << z  
        << "\nwith rphi/z variances " << localCov(0, 0) 
        << " " << localCov(1, 1) << " and rotated variances " << var[0] << " " << var[1]
        << "\nand a geometry ID " << surf->geometryId() << std::endl;
    }*/

    return spPtr;
}

std::vector<const SpacePoint*> getDiskSpacePoints(const Acts::GeometryContext& gctx,
    Acts::Extent& rRangeSPExtent, int i) {

    bool debug = true;
    std::vector<const SpacePoint*> spVector;
    unsigned int nHits = 0;
    double radius    = 0.49_cm;      // StrawTubeOuterRadius

	// get the tree and iterate over the branches
	auto *f = TFile::Open("detector_oa_all_t.root");
	if (!f) {
		std::cout << "Couldn't open measurement file\n";
	}

	auto *tree = (TTree*)f->Get("detections");
	std::vector<ROOT::Math::XYZVector> *spOfEvent = 0;
    std::vector<ROOT::Math::XYZVector> *momOfEvent = 0;
	tree->SetBranchAddress("spOfEvent", &spOfEvent);
    tree->SetBranchAddress("momOfEvent", &momOfEvent);

    auto& second_stat = tracker_layout["Station_1"][0]; // First layer of second station
    //std::cout << "Number of entries: " << tree->GetEntries() << std::endl;

    int stationID = 0; double StationD = 39.0_cm; 
    //for (size_t ev = 0; ev < tree->GetEntries(); ++ev) {
    for (size_t ev = i; ev < (i + 1); ++ev) {
        std::vector<Acts::Vector3> posVec = {Acts::Vector3::Zero(), Acts::Vector3::Zero(), 
            Acts::Vector3::Zero(), Acts::Vector3::Zero()};
        
        std::vector<int> k = {0, 0, 0, 0};
        tree->GetEntry(ev);
        for (size_t it = 0; it < spOfEvent->size(); ++it) {
            double positionZ = (*spOfEvent)[it].Z() - 333.45_cm;
            if ((positionZ >= (-1380.25 + StationD * (stationID+1) - 1_cm)) && 
                (positionZ <= (-1380.25 + StationD * (stationID-1) + 1_cm))){
                continue;}

            if (positionZ <= -1380.25 + StationD * stationID + radius) {
                posVec[0].x() += (*spOfEvent)[it].X();
                posVec[0].y() += (*spOfEvent)[it].Y();
                posVec[0].z() += (*spOfEvent)[it].Z();
                k[0] += 1;
            }
            else if (positionZ <= -1370.9 + StationD * stationID + radius) {
                posVec[1].x() += (*spOfEvent)[it].X();
                posVec[1].y() += (*spOfEvent)[it].Y();
                posVec[1].z() += (*spOfEvent)[it].Z();
                k[1] += 1;
            }
            else if (positionZ <= -1359.1 + StationD * stationID + radius) {
                posVec[2].x() += (*spOfEvent)[it].X();
                posVec[2].y() += (*spOfEvent)[it].Y();
                posVec[2].z() += (*spOfEvent)[it].Z();
                k[2] += 1;
            }
            else if (positionZ <= -1349.75 + StationD * stationID + radius) {
                posVec[3].x() += (*spOfEvent)[it].X();
                posVec[3].y() += (*spOfEvent)[it].Y();
                posVec[3].z() += (*spOfEvent)[it].Z();
                k[3] += 1;
            }
        }

        auto p = Acts::Vector3{(*momOfEvent)[0].X(), (*momOfEvent)[0].Y(), (*momOfEvent)[0].Z()};

        int kount = 0;
        for (auto pos : posVec) {
            for (auto& surface : tracker_layout["Station_" + std::to_string(stationID)]) {
                Acts::Vector3 surf_center = surface->center(gctx);
                if (k[kount] == 0) {continue;}
                if (abs(pos.z()/k[kount] - 333.45_cm - surf_center.z()) < radius) {
                    // TODO: Handle the radial distance in z pos.z()/k[kount]
                    auto spacePoint = makeSpacePoint(gctx, 
                        Acts::Vector3{pos.x()/k[kount], pos.y()/k[kount], surf_center.z() + 333.45_cm}, 0, surface, ev, p).release();
                    spVector.push_back(spacePoint);
                    rRangeSPExtent.extend({spacePoint->x(), spacePoint->y(), spacePoint->z()});
                    nHits++;
                    //std::cout << pos.x()/k[kount]-10800 << " " << pos.y()/k[kount] << " " << pos.z()/k[kount]-3334.5 << std::endl;
                }
            } // surfaces in the first station
            kount++;
        } // spacepoints
        if (spVector.size() > 2) {
            event_map[ev] = std::make_pair(0, p); // initialize values for the seed event map
        }
        else {continue;} // not enough hits
    } // loop over events
    
    f->Close();

    return spVector;
}


int main() {
    std::cout << "Hello world: I am alive." << std::endl;
    estimatorFile.open("estimations.txt", std::ios::out);
    // Initialize geometry context
	Acts::GeometryContext* gctx;
    bool debug = false;
    // Build geometry
    std::shared_ptr<const Acts::TrackingGeometry> tGeometry{nullptr};
    //tGeometry = aBuildStrawGeometry(*gctx, debug);
    tGeometry = aBuildDiskGeometry(*gctx, debug);
	if (!tGeometry) {
    		throw std::runtime_error("ACTS::BuilGeometry: tGeometry is null");
    	}
	std::cout << "ACTS:: Tracking geometry has been created." << std::endl;
 	//-----------------------------------------------------------------------------------

    // Initialize magnetic field context
	Acts::MagneticFieldContext magctx;

	// Get magnetic field
	auto mapBins = [](std::array<size_t, 3> bins, std::array<size_t, 3> sizes) {
        	return (bins[0] * (sizes[1] * sizes[2]) + bins[1] * sizes[2] + bins[2]);
      	};

	auto BFieldmap = aGetMagneticField(debug, magctx, std::move(mapBins));
	auto bField = std::dynamic_pointer_cast<const Acts::InterpolatedMagneticField>(BFieldmap);
    	if (!bField) {
    		throw std::runtime_error("ACTS::GetMagneticField: bField is null");
    	}  
    auto bCache = BFieldmap->makeCache(magctx);
	std::cout << "ACTS:: Magnetic field has been interpolated from root file." << std::endl;
	//-----------------------------------------------------------------------------------

     // Initialize all the histograms
    TH1F *h_pEstim = new TH1F("pEstim", "With any layer combination", 100, -150, 200);
    TH1F *h_pEstim_l1 = new TH1F("pEstim_l1", "With layers 0 1 2", 100, -150, 200);
    TH1F *h_pEstim_l2 = new TH1F("pEstim_l2", "With layers 0 1 3", 100, -150, 200);
    TH1F *h_pEstim_l3 = new TH1F("pEstim_l3", "With layers 0 2 3", 100, -150, 200);
    TH1F *h_pEstim_l4 = new TH1F("pEstim_l4", "With layers 1 2 3", 100, -150, 200);
    std::fstream estimFile;
    estimFile.open("estim_stat.txt", std::ios::out);

    auto logger = Acts::getDefaultLogger("ParamEstimator", Acts::Logging::VERBOSE);
    std::vector<std::vector<int>> combi{{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};
    for (int i = 0; i < 100; i++) {  // 7514 31013
        // Extent used to store r range for middle spacepoint
        Acts::Extent rRangeSPExtent;
        std::cout << "---------------- Ev " << i << " -------------------" << std::endl;
        // Input space points
        std::vector<const SpacePoint*> inputSP = getDiskSpacePoints(*gctx, rRangeSPExtent, i);
        int j = 0;
        if (inputSP.size() < 3) {continue;}
        for (const auto& combo : combi) {
            if ((inputSP.size() < 4) && (j > 0)) {continue;}
            std::vector<const SpacePoint*> seedSP{inputSP[combo[0]], inputSP[combo[1]], inputSP[combo[2]]};
            
            std::cout << seedSP[0]->z() << " " << seedSP[1]->z() << " " << seedSP[2]->z() << std::endl;
            // Sort SP in z-axis
            std::sort(seedSP.begin(), seedSP.end(), [](const SpacePoint* lhs, 
                    const SpacePoint* rhs) {return (lhs->z() < rhs->z());});

            std::cout << seedSP[0]->z() << " " << seedSP[1]->z() << " " << seedSP[2]->z() << std::endl;
            
            const auto rSurface = tGeometry->findSurface(seedSP[0]->geoID());
            auto field_res = BFieldmap->getField(Acts::Vector3{seedSP[0]->x(), 
                seedSP[0]->y(), seedSP[0]->z()}, bCache);
            auto field = field_res.value();
            // Minimum magnetic field required to trigger the q/pt estimation (field>min)
            Acts::ActsScalar bFieldMin = 0;
            Acts::ActsScalar eMass = 0.511_MeV;
            auto trackParams = (Acts::estimateTrackParamsFromSeed(*gctx, seedSP.begin(), 
                seedSP.end(), *rSurface, field, bFieldMin, *logger, eMass)).value();

            auto trackParams2 = (estimateModified(*gctx, seedSP.begin(), 
                seedSP.end(), *rSurface, field, bFieldMin, *logger, eMass)).value();

            double p_estim = std::abs(- 1 / trackParams[Acts::eBoundQOverP] * 1000);
            double p_real  = (seedSP[0]->p()).norm();
            
            if (j == 0) {h_pEstim_l1->Fill(p_estim);}
            if (j == 1) {h_pEstim_l2->Fill(p_estim);}
            if (j == 2) {h_pEstim_l3->Fill(p_estim);}
            if (j == 3) {h_pEstim_l4->Fill(p_estim);}
            h_pEstim->Fill(p_estim);
            estimFile << i << " " << std::abs(p_estim) << " " << std::abs(p_estim - p_real) << "\n";
            j++;
        }
    }
    h_pEstim->GetXaxis()->SetTitle("p = p_{estim}-p_{real} [MeV/c]");
	h_pEstim->GetYaxis()->SetTitle("Entries");
    h_pEstim->SetStats(0);
    
    TCanvas *c0 = new TCanvas("c0", "c0", 1500, 1000);
	h_pEstim->Draw();
    h_pEstim->SetLineColor(1);
    h_pEstim_l1->Draw("SAME HIST");
    h_pEstim_l1->SetLineColor(2);
    h_pEstim_l2->Draw("SAME HIST");
    h_pEstim_l2->SetLineColor(4);
    h_pEstim_l3->Draw("SAME HIST");
    h_pEstim_l3->SetLineColor(6);
    h_pEstim_l4->Draw("SAME HIST");
    h_pEstim_l4->SetLineColor(8);
    gPad->BuildLegend(0.53,0.67,0.88,0.88, "Estimation of Momentum");
    //TPaveLabel *t = new TPaveLabel(0.15, 0.9, 0.85, 1.0, "Resolution of estimated momentum", "brNDC"); 
    //t->SetBorderSize(0);
    //t->SetFillColor(gStyle->GetTitleFillColor());
    //t->Draw();
    h_pEstim->SetTitle(0);
    c0->Print("histograms/pEstim_combi.pdf");
    estimFile.close();
    return 0;
}
// ---------------------------------------------------------------- //
//                Creator: Amaia Razquin Lizarraga                  //
//                Contact: amaiafisika@gmail.com                    //
//                Date:    7 September 2023                         //
// ---------------------------------------------------------------- //
/* 
    This program runs the seeding algorithm and explicitly performs the 
    conformal mapping within the initial parameter estimator in order to 
    visualize the differences in the mapping between seeds from the same
    trajectory. The mapping points are outputed in seedMap.root and can
    be plotted with mapSeeds.c in root. The time difference between SPs 
    of the same track is also calculated and the maximum is returned, as 
    well as a plot.
    Command:
        g++ seedMapping.cpp Geometry.cpp MagField.cpp Seeder.cpp -o mapping.x
        -L/path/to/acts_install/lib -lActsPluginTGeo -lActsPluginJson -lActsCore 
        -I/path/to/acts_install/include/ -I$(root-config --incdir) 
        -I$(root-config --evelibs) $(root-config --libs) 
        -lGeom -I/usr/include/eigen3
*/

#include <iostream>
#include <vector>

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/Units.hpp>

#include <Acts/EventData/TrackParameters.hpp>

#include <Acts/Surfaces/PerigeeSurface.hpp>

#include <Acts/EventData/VectorTrackContainer.hpp>

#include <Acts/Seeding/BinnedSPGroup.hpp>
#include <Acts/Seeding/InternalSeed.hpp>
#include <Acts/Seeding/InternalSpacePoint.hpp>
#include <Acts/Seeding/Seed.hpp>
#include <Acts/Seeding/SeedFinder.hpp>
#include <Acts/Seeding/SeedFilter.hpp>

#include <Acts/TrackFitting/GainMatrixSmoother.hpp>
#include <Acts/TrackFitting/GainMatrixUpdater.hpp>

#include <Acts/TrackFinding/CombinatorialKalmanFilter.hpp>
#include <Acts/TrackFinding/MeasurementSelector.hpp>

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

// Amaia:
#include "base/SpacePoint.hpp"
#include "base/ActsSourceLink.hpp"
#include "base/Calibrator.hpp"
#include "Geometry.hpp"
#include "MagField.hpp"
#include "Seeder.hpp"

using namespace Acts::UnitLiterals;

TFile* out_file;
TTree* out_tree;

extern std::map<std::string, std::vector<std::shared_ptr<const Acts::Surface>>> tracker_layout;

int main() {
    std::cout << "Hello world: I am alive." << std::endl;

    // ----------- Geometry -------------------------------------------------------
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
    // ----------------------------------------------------------------------------

    // ----------- Magnetic field -------------------------------------------------
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
    // ----------------------------------------------------------------------------

    // ----------- Input measurements ---------------------------------------------
    auto *f = TFile::Open("detector_oa_all_td.root");
	if (!f) {
		std::cout << "Couldn't open measurement file\n";
	}

	auto *tree = (TTree*)f->Get("detections");
	std::vector<ROOT::Math::XYZVector> *spOfEvent = 0;
    std::vector<ROOT::Math::XYZVector> *momOfEvent = 0;
	std::vector<double> *timeOfEvent = 0;
    std::vector<double> *eDeposition = 0;
	tree->SetBranchAddress("spOfEvent", &spOfEvent);
    tree->SetBranchAddress("momOfEvent", &momOfEvent);
    tree->SetBranchAddress("timeOfEvent", &timeOfEvent);
    tree->SetBranchAddress("eDeposition", &eDeposition);
    // ----------------------------------------------------------------------------

    // Output file ----------------------------------------------------------------
    out_file = TFile::Open("seedMap.root", "recreate");
    out_tree = new TTree("seeds", "seeds");

    std::vector<ROOT::Math::XYZVector> originalPos;
    std::vector<ROOT::Math::XYZVector> newPos;
    std::vector<ROOT::Math::XYZVector> comformalPos;
    std::vector<int> layer;
    int event;

    out_tree->Branch("originalPos", &originalPos);
    out_tree->Branch("newPos", &newPos);
    out_tree->Branch("comformalPos", &comformalPos);
    out_tree->Branch("layer", &layer);
    out_tree->Branch("event", &event);

    // ----------- Seeding tools --------------------------------------------------
    // Parameters
    double radLperSeed  = 5.381265264757346;
    double rMax         = 690.1954797042055;
    double deltaRMin    = 0.2113334751419763;
    double deltaRMax    = 104.36756509015152;
    double deltaZMax    = 29.233848440510886;
    double binSizeR     = 0.6922886787481911;
    double deltaRMinMid = 0.15801016669789986;
    double deltaRMaxMid = 0.25805540458126863;
    double collMin      = -44478.60870408864;
    double collMax      = 3395.7357360193287;
    double cotThetaMax  = 1052.032750551505;
    double minPt        = 128.04907527478224;
    double impactMax    = 4914.642127968815;
    double sigmaScat    = 2.624471976102091;
    double maxPtScat    = 389.1264111786207;
    // Grid
    double gminPt       = 119.09021636183944;
    double grMax        = 570.8036490784622;
    double gdeltaRMax   = 141.52768431029594;
    double gimpactMax   = 98.48661116503268;
    // Filter
    double fdeltaInv    = 0.0005257695602828333;
    double fdeltaRMin   = 1.5583724519088409;
    int fimpactWF       = 6;
    int fzOriginWF      = 8;
    int fcompatSW       = 37;
    int fcompatSL       = 4;
    int fseedWI         = 3;

    // Seed Finder Configuration
    Acts::SeedFinderConfig<SpacePoint> seedingConfig;
    // Limiting location of measurements (e.g. detector contraints)
    seedingConfig.rMin                     = 0._mm;
    seedingConfig.rMax                     = rMax * 1._mm;
    seedingConfig.zMin                     = -1385.15_mm; //-2800._mm;
    seedingConfig.zMax                     = 1385.15_mm; //2800._mm;
    // Min/max distance between two measurements in one seed
    seedingConfig.deltaRMin                = deltaRMin * 1._mm; //1._mm;
    seedingConfig.deltaRMax                = deltaRMax * 1._mm; //16._mm;
    seedingConfig.deltaZMax                = deltaZMax * 1._mm;
    seedingConfig.deltaRMinTopSP           = seedingConfig.deltaRMin;
    seedingConfig.deltaRMinBottomSP        = seedingConfig.deltaRMin;
    seedingConfig.deltaRMaxTopSP           = seedingConfig.deltaRMax;
    seedingConfig.deltaRMaxBottomSP        = seedingConfig.deltaRMax;
    seedingConfig.binSizeR                 = binSizeR * 1._mm; //R
    // Variable range based on the SP radius
    seedingConfig.useVariableMiddleSPRange = true;
    seedingConfig.deltaRMiddleMinSPRange   = deltaRMinMid * 1._mm; // 1_mm;
    seedingConfig.deltaRMiddleMaxSPRange   = deltaRMaxMid * 1._mm; //600_mm;
    // Limiting collision region in z
    seedingConfig.collisionRegionMin       = collMin * 1._mm;
    seedingConfig.collisionRegionMax       = collMax * 1._mm;
    seedingConfig.sigmaScattering          = sigmaScat; // was 1, default is 5
    seedingConfig.maxPtScattering          = maxPtScat * 1._MeV; // 10 GeV default
    seedingConfig.cotThetaMax              = cotThetaMax * 1._mm; //7.40627;
    seedingConfig.minPt                    = minPt * 1._MeV; //25._MeV;
    // How many seeds a given hit can be the middle hit of the seed
    seedingConfig.maxSeedsPerSpM           = 2;
    // Maximum impact parameter must be smaller than rMin
    seedingConfig.impactMax                = impactMax * 1._mm; //20._mm;
    // Average radiation length traversed per seed
    seedingConfig.radLengthPerSeed         = radLperSeed; // it is a percentage
    // Custom bins in z
    seedingConfig.zBinEdges = std::vector<float>{-1385.15, -1366, -1344.55};
    seedingConfig.zBinsCustomLooping = {1, 2, 3};

    // Configure seed filter
    Acts::SeedFilterConfig seedFilterConfig;
    seedFilterConfig.maxSeedsPerSpM        = 2;
    seedFilterConfig.deltaInvHelixDiameter = fdeltaInv * 1._mm;
    seedFilterConfig.deltaRMin             = fdeltaRMin * 1._mm;
    seedFilterConfig.impactWeightFactor    = fimpactWF;
    seedFilterConfig.zOriginWeightFactor   = fzOriginWF;
    seedFilterConfig.compatSeedWeight      = fcompatSW;
    seedFilterConfig.compatSeedLimit       = fcompatSL;
    seedFilterConfig.seedWeightIncrement   = fseedWI;
    seedFilterConfig = seedFilterConfig.toInternalUnits();
    seedingConfig.seedFilter = std::make_unique<Acts::SeedFilter<SpacePoint>>(
        Acts::SeedFilter<SpacePoint>(seedFilterConfig));

    seedingConfig = seedingConfig.toInternalUnits();
    auto seedFinder = Acts::SeedFinder<SpacePoint>(seedingConfig);
    // Seed Finder Options
    Acts::SeedFinderOptions seedingOptions;
    seedingOptions.beamPos          = {0._mm, 0._mm};
    seedingOptions.bFieldInZ        = 1.0_T;

    seedingOptions = seedingOptions.toInternalUnits();

    // Covariance converter function needed by the finder
    auto covConverter = 
        [=](const SpacePoint& sp, float, float, float)
        -> std::pair<Acts::Vector3, Acts::Vector2> {
            Acts::Vector3 position{sp.x(), sp.y(), sp.z()};
            Acts::Vector2 covariance{sp.varianceR(), sp.varianceZ()};
            return std::make_pair(position, covariance); 
    };

    // Vector containing the map of z bins in the top and bottom layers
    std::vector<std::pair<int, int>> zBinNeighborsTop;
    std::vector<std::pair<int, int>> zBinNeighborsBottom;
    // Number of phiBin neighbors at each side of the current bin that will be used
    // to search for SPs
    int numPhiNeighbors = 1;

    // Setup Space Point Grid configuration
    // minPt, rMax and deltaRMax set the number of bins in phi (-pi,pi)
    // zBinEdges sets explicitly the grid in z
    Acts::SpacePointGridConfig gridConfig;
    gridConfig.minPt       = gminPt * 1.0_MeV;
    gridConfig.rMax        = grMax * 1.0_mm;
    gridConfig.deltaRMax   = gdeltaRMax * 1.0_mm;
    gridConfig.impactMax   = gimpactMax * 1.0_mm;
    gridConfig.zBinEdges   = seedingConfig.zBinEdges;
    gridConfig = gridConfig.toInternalUnits();
        
    // Setup Space Point Grid options
    Acts::SpacePointGridOptions gridOptions;
    gridOptions.bFieldInZ  = 1.0_T;
    gridOptions = gridOptions.toInternalUnits();
        
    // Create grid with bin sizes according to the configured geometry
    std::unique_ptr<Acts::SpacePointGrid<SpacePoint>> grid = 
        Acts::SpacePointGridCreator::createGrid<SpacePoint>(gridConfig, gridOptions);

	//-----------------------------------------------------------------------------

    TH1F *h_time = new TH1F("h_time", "Timing difference between bottom and top SP; t [ns]; Entries", 25, 0, 0.3);
    
    // ----------- Event loop -----------------------------------------------------   
    std::vector<Acts::Vector3> magList;
    double tMax = 0;
    //for (int eventID = 0; eventID < tree->GetEntries(); eventID++) {
    for (int eventID = 0; eventID < 10; eventID++) {
        tree->GetEntry(eventID);
        //std::cout << "Size of event: " << (*spOfEvent).size() << std::endl;
        if ((*spOfEvent).size() < 28) {std::cout << "Ignoring event" << std::endl; continue;}

        // Randomly add time of creation to match BG arrival time
        TRandom3* rndm = new TRandom3();
        rndm->SetSeed(0);
        double creationT = rndm->Uniform(620, 1000); //rndm->Exp(864)
        
        for (int i = 0; i < 1; i++) {
            Acts::Extent rRangeSPExtent;
            auto bottomBinFinder = std::make_shared<Acts::BinFinder<SpacePoint>>(
                Acts::BinFinder<SpacePoint>(zBinNeighborsBottom, numPhiNeighbors));
            auto topBinFinder    = std::make_shared<Acts::BinFinder<SpacePoint>>(
                Acts::BinFinder<SpacePoint>(zBinNeighborsTop, numPhiNeighbors));
            std::unique_ptr<Acts::SpacePointGrid<SpacePoint>> grid = 
                Acts::SpacePointGridCreator::createGrid<SpacePoint>(gridConfig, gridOptions);
            int stationID = i;
            auto spVector = getSpacePoints(*gctx, *spOfEvent, *momOfEvent, *eDeposition, *timeOfEvent, 
                creationT, tracker_layout, eventID, stationID, rRangeSPExtent);

            if (spVector.size() < 3) {
                std::cout << "Station " << stationID << ": Not enough hits." << std::endl;
                continue;
            }
            
            auto spGroup = Acts::BinnedSPGroup<SpacePoint>(
                spVector.begin(), spVector.end(), covConverter, 
                bottomBinFinder, topBinFinder, std::move(grid), rRangeSPExtent, 
                seedingConfig, seedingOptions);
                
            // variable middle SP radial region of interest
            const Acts::Range1D<float> rMiddleSPRange(std::floor(
                rRangeSPExtent.min(Acts::binR) * 0.5) * 2 + seedingConfig.deltaRMiddleMinSPRange,
                std::floor(rRangeSPExtent.max(Acts::binR) * 0.5) * 2 
                - seedingConfig.deltaRMiddleMaxSPRange);

            // Seeding
            SeedContainer seeds;
            seeds.clear();

            decltype(seedFinder)::SeedingState state;
            state.spacePointData.resize(spVector.size());

            auto start = std::chrono::system_clock::now();
            
            // Run seeding
            for (const auto [bottom, middle, top] : spGroup) {
                seedFinder.createSeedsForGroup(seedingOptions, state, spGroup.grid(), 
                    std::back_inserter(seeds), bottom, middle, top, rMiddleSPRange);
            }

            auto end = std::chrono::system_clock::now();
            std::chrono::duration<double> elapsed_seconds = end - start;

            // Estimate initial track parameters and get spacepoints to track find from
            for (size_t j = 0; j < seeds.size(); j++) {
                originalPos.clear();
                newPos.clear();
                comformalPos.clear();
                layer.clear();
                // Create vector and estimate parameters
                std::vector<const SpacePoint*> seed{seeds[j].sp()[0], 
                    seeds[j].sp()[1], seeds[j].sp()[2]};

                // Check maximum delta t between sp in the same seed
                double deltaT = std::abs(seed[0]->t() - seed[2]->t());
                //std::cout << seed[0]->t() << " " << seed[1]->t() << " " << seed[2]->t() << " " << deltaT << std::endl;
                h_time->Fill(deltaT);
                if (deltaT > tMax) {tMax = deltaT;}

                auto paramsA = estimateParams(*gctx, tGeometry, magctx, BFieldmap, seed);
                double p_estimA = - 1 / paramsA[Acts::eBoundQOverP] * 1000;

                // Get the mapping
                // Global position helper
                std::array<Acts::Vector3, 3> spGlobalPositions = {Acts::Vector3::Zero(), Acts::Vector3::Zero(),
                    Acts::Vector3::Zero()};
                auto bCache = BFieldmap->makeCache(magctx);
                auto field_res = BFieldmap->getField(Acts::Vector3{seed[0]->x(), 
                    seed[0]->y(), seed[0]->z()}, bCache);
                auto bField = field_res.value();
                
                // Get original positions
                int i = 0;
                for (auto sp : seed) {
                    originalPos.push_back(ROOT::Math::XYZVector{sp->x(), sp->y(), sp->z()});
                    layer.push_back(sp->geoID().layer());
                    spGlobalPositions[i] = Acts::Vector3(sp->x(), sp->y(), sp->z());
                    i++;
                }
                event = eventID;

                // Transform SpacePoints to new coordinate system ---------------------------------
                Acts::Vector3 relVec   = spGlobalPositions[1] - spGlobalPositions[0];
                Acts::Vector3 newZAxis = bField.normalized();
                Acts::Vector3 newYAxis = newZAxis.cross(relVec).normalized();
                Acts::Vector3 newXAxis = newYAxis.cross(newZAxis);
                Acts::RotationMatrix3 rotation;
                rotation.col(0) = newXAxis;
                rotation.col(1) = newYAxis;
                rotation.col(2) = newZAxis;
                // The center of the new frame is at the bottom space point
                Acts::Translation3 trans(spGlobalPositions[0]);
                // The transform which constructs the new frame
                Acts::Transform3 transform(trans * rotation);
                // The coordinate of the middle and top space point in the new frame
                Acts::Vector3 local1 = transform.inverse() * spGlobalPositions[1]; // middle
                Acts::Vector3 local2 = transform.inverse() * spGlobalPositions[2]; // top

                // Get new positions
                ROOT::Math::XYZVector bottom = {0, 0, 0};
                ROOT::Math::XYZVector middle = {local1.x(), local1.y(), local1.z()};
                ROOT::Math::XYZVector top    = {local2.x(), local2.y(), local2.z()};
                newPos.push_back(bottom); newPos.push_back(middle); newPos.push_back(top);

                // Transform spacepoints to comformal space ---------------------------------------
                auto uvTransform = [](const Acts::Vector3& local) -> Acts::Vector2 {
                    Acts::Vector2 uv;
                    Acts::ActsScalar denominator = local.x() * local.x() + local.y() * local.y();
                    uv.x() = local.x() / denominator;
                    uv.y() = local.y() / denominator;
                    return uv;
                };
                // The uv1.y() should be zero
                Acts::Vector2 uv1 = uvTransform(local1);
                Acts::Vector2 uv2 = uvTransform(local2);
                
                // Get comformal positions
                ROOT::Math::XYZVector cBottom = {0, 0, 0};
                ROOT::Math::XYZVector cMiddle = {uv1.x(), uv1.y(), 0};
                ROOT::Math::XYZVector cTop    = {uv2.x(), uv2.y(), 0};
                comformalPos.push_back(cBottom); comformalPos.push_back(cMiddle); comformalPos.push_back(cTop);

                out_tree->Fill();
            } // seed loop
        } // Station loop
    } // event loop 
    std::cout << "MAXIMUM TIME DIFFERENCE " << tMax << std::endl;
    TCanvas* c0 = new TCanvas();
    h_time->Draw();
    c0->Print("histograms/timingDif.png");
    
    f->Delete();
    out_file->Write();
    out_file->Print();
    out_file->Close();
}
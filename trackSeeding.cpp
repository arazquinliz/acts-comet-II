// ---------------------------------------------------------------- //
//                Creator: Amaia Razquin Lizarraga                  //
//                Contact: amaiafisika@gmail.com                    //
//                Date:    7 September 2023                         //
// ---------------------------------------------------------------- //
/* 
    This program performs track seeding with a set of input measurements.
    It loops over all of the events in a file and returns the efficiency
    values of the seeding. It also applies additional cuts and outputs 
    the seeded spacepoints in seeded_meas.txt.

    Command:
        g++ trackSeeding.cpp Geometry.cpp MagField.cpp Seeder.cpp 
        KalmanFitting.cpp -o seeding.x
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

#include <Acts/Propagator/Propagator.hpp>
#include <Acts/Propagator/EigenStepper.hpp>
#include <Acts/Propagator/Navigator.hpp>
#include <Acts/Propagator/MaterialInteractor.hpp>
#include <Acts/Propagator/StandardAborters.hpp>
#include <Acts/Propagator/detail/SteppingLogger.hpp>

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
#include "KalmanFitting.hpp"

using namespace Acts::UnitLiterals;

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

    // ----------- Event loop -----------------------------------------------------   
    // Score function parameters
    int totalSignal        = 0;
    int totalBG            = 0;
    int totalMatched       = 0;
    int totalEvents        = 0;

    std::map<int, std::vector<SpacePoint>> spLayout;

    //for (int eventID = 0; eventID < tree->GetEntries(); eventID++) {
    for (int eventID = 0; eventID < 10; eventID++) {
        tree->GetEntry(eventID);
        std::cout << "Size of event: " << (*spOfEvent).size() << std::endl;
        if ((*spOfEvent).size() < 28) {std::cout << "Ignoring event" << std::endl; continue;}
        totalEvents++;

        std::vector<double> momList; // helper vector for initial track parameters
        std::vector<SpacePoint> spList; // helper vector for spacepoints
        std::vector<Acts::CurvilinearTrackParameters> initialParameters;
        SourceLinkContainer sourceLinks;
        std::vector<Acts::BoundVariantMeasurement> measurements;

        int bg     = 0;
        int signal = 0;
        std::fstream seededfile;
        seededfile.open("seeded_meas.txt", std::ios::out); 
        for (int i = 0; i < 8; i++) {
            Acts::Extent rRangeSPExtent;
            auto bottomBinFinder = std::make_shared<Acts::BinFinder<SpacePoint>>(
                Acts::BinFinder<SpacePoint>(zBinNeighborsBottom, numPhiNeighbors));
            auto topBinFinder    = std::make_shared<Acts::BinFinder<SpacePoint>>(
                Acts::BinFinder<SpacePoint>(zBinNeighborsTop, numPhiNeighbors));
            std::unique_ptr<Acts::SpacePointGrid<SpacePoint>> grid = 
                Acts::SpacePointGridCreator::createGrid<SpacePoint>(gridConfig, gridOptions);
            int stationID = i;
            auto spVector = getSpacePoints(*gctx, *spOfEvent, *momOfEvent, *eDeposition, *timeOfEvent, 0, tracker_layout, 
                eventID, stationID, rRangeSPExtent);

            if (spVector.size() < 3) {
                std::cout << "Station " << stationID << ": Not enough hits." << std::endl;
                continue;
            }
            totalEvents++;
            
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
            std::cout << "Station " << stationID << ": " << spVector.size() << 
                " hits, number of seeds: " << seeds.size() << std::endl;

            // Estimate initial track parameters and get spacepoints to track find from
            for (size_t j = 0; j < seeds.size(); j++) {
                // Create vector and estimate parameters
                std::vector<const SpacePoint*> seedSPA{seeds[j].sp()[0], 
                    seeds[j].sp()[1], seeds[j].sp()[2]};
                auto paramsA = estimateParams(*gctx, tGeometry, magctx, BFieldmap, seedSPA);
                double p_estimA = - 1 / paramsA[Acts::eBoundQOverP] * 1000;
                
               

                for (size_t k = (j + 1); k < seeds.size(); k++) {
                    // Create vector and estimate parameters
                    std::vector<const SpacePoint*> seedSPB{seeds[k].sp()[0], 
                    seeds[k].sp()[1], seeds[k].sp()[2]};
                    auto paramsB = estimateParams(*gctx, tGeometry, magctx, BFieldmap, seedSPB);
                    double p_estimB = - 1 / paramsB[Acts::eBoundQOverP] * 1000;
                    // Compare seeds
                    auto seedComb = compareSeeds(seeds[j], seeds[k]);

                    // If coincidence was found in the first station, get initial parameters
                    // and combine seeds
                    // TODO: consider doing this in a separate function
                    if ((i == 0) && (seedComb.size() == 4)) {
                        // Sort seeds in z-axis
                        std::sort(seedComb.begin(), seedComb.end(), [](const SpacePoint lhs, 
                            const SpacePoint rhs) {return (lhs.z() < rhs.z());});
                        
                        // If one layer contains two sp, ignore combination
                        bool repetition = false;
                        double prev_z   = 0.;
                        for (int spIt = 0; spIt < seedComb.size(); spIt++) {
                            if (seedComb[spIt].z() == prev_z) {repetition = true; continue;}
                            prev_z = seedComb[spIt].z();
                        }
                        if (repetition) {continue;}

                        // Time progression cut
                        if ((seedComb[3].t() - seedComb[2].t() < - 1.5) || 
                            (seedComb[2].t() - seedComb[1].t() < - 1.5) ||
                            (seedComb[1].t() - seedComb[0].t() < - 1.5)) {continue;}

                        auto purity = checkPurity(seedComb);
                        std::cout << "Purity is: " << purity.second << " for particle " << 
                            purity.first << " " << p_estimA << " " << p_estimB << std::endl;

                        // Initial parameters
                        Acts::Vector3 pos = Acts::Vector3{seedComb[0].x(), seedComb[0].y(), seedComb[0].z()};
                        double time  = 0.;
                        double dX    = seedComb[1].x() - seedComb[0].x();
                        double dY    = seedComb[1].y() - seedComb[0].y();
                        double dZ    = seedComb[1].z() - seedComb[0].z();
                        // This angles should be of the estimated direction
                        double phi   = dY/abs(dY) * atan(dX / std::sqrt(dX * dX + dY * dY));
                        double theta = acos(dZ / std::sqrt(dX * dX + dY * dY + dZ * dZ));
                        int    help  = p_estimA < 0 ? 2 : 0;
                        // Idealy this should me estimated in the first layer
                        double p     = seeds[j].sp()[help]->z() < seeds[k].sp()[help]->z() ? 
                            abs(p_estimA) : abs(p_estimB);
                        int    q     = -1;
                        auto cov     = setDefaultCovariance();
                        auto curvPar = Acts::CurvilinearTrackParameters(
                            Acts::VectorHelpers::makeVector4(pos, time), phi, theta, p*1._MeV, q, cov);
                        
                        // Only add new initial parameters
                        if (std::find(momList.begin(), momList.end(), p) == momList.end()) {
                            momList.push_back(p);
                            initialParameters.push_back(curvPar);
                            //std::cout << pos.x() << " " << pos.y() << " " << pos.z() << " " << phi << " " 
                            //   << theta << " " << p << std::endl;
                        }
                    }
                    
                    if (seedComb.size() == 4) {
                        // Sort seeds in z-axis
                        std::sort(seedComb.begin(), seedComb.end(), [](const SpacePoint lhs, 
                            const SpacePoint rhs) {return (lhs.z() < rhs.z());});
                            
                        // If one layer contains two sp, ignore combination
                        bool repetition = false;
                        double prev_z   = 0.;
                        for (int spIt = 0; spIt < seedComb.size(); spIt++) {
                            if (seedComb[spIt].z() == prev_z) {repetition = true; continue;}
                            prev_z = seedComb[spIt].z();
                        }
                        if (repetition) {continue;}

                        // Time progression cut
                        if ((seedComb[3].t() - seedComb[2].t() < - 1.5) || 
                            (seedComb[2].t() - seedComb[1].t() < - 1.5) ||
                            (seedComb[1].t() - seedComb[0].t() < - 1.5)) {continue;}

                        auto purity = checkPurity(seedComb);
                        std::cout << "Purity is: " << purity.second << " for particle " << 
                            purity.first << std::endl;

                        for (auto& sp : seedComb) {
                            // Get sourcelinks and measurements only for matching seeds
                            if (std::find(spList.begin(), spList.end(), sp) == spList.end()) {
                                ActsSourceLink::Index index = measurements.size();
                                ActsSourceLink sl(sp.geoID(), index);
                                sourceLinks.emplace(sl.geometryId(), std::move(sl));
                                spList.push_back(sp);
                                if (sp.event() > tree->GetEntries()) {
                                    bg++;
                                }
                                else {
                                    seededfile << sp.x() << " " << sp.y() << " " << sp.z() << "\n"; 
                                    signal++;
                                }
                                spLayout[sp.event()].push_back(sp);
                                std::cout << sp.x() << " " << sp.y() << " " << sp.z() << " " 
                                    << sp.t() << " " << sp.event() << std::endl;

                                std::array<Acts::BoundIndices, 2> indices;
                                indices[0] = Acts::BoundIndices::eBoundLoc0;
                                indices[1] = Acts::BoundIndices::eBoundLoc1;

                                const auto rSurface    = tGeometry->findSurface(sp.geoID());
                                Acts::Vector2 loc = rSurface->globalToLocal(*gctx, 
                                    Acts::Vector3( sp.x(), sp.y(), sp.z() ),
                                    Acts::Vector3::Zero()).value();

                                Acts::ActsVector<2> loca;
                                loca[Acts::eBoundLoc0] = loc.x();
                                loca[Acts::eBoundLoc1] = loc.y();


                                Acts::ActsSquareMatrix<2> cova = Acts::ActsSquareMatrix<2>::Zero();
                                cova(Acts::eBoundLoc0, Acts::eBoundLoc0) = 0.0_mm; 
                                cova(Acts::eBoundLoc1, Acts::eBoundLoc1) = 0.0_mm;
                
                                Acts::Measurement<Acts::BoundIndices, 2> meas(Acts::SourceLink{sl}, indices, loca, cova);
                                measurements.push_back(meas);

                            }  
                        } // spacepoint
                    } // if size == 4
                } // seedB
            } // seedA
        } // station
        totalBG      += bg;
        totalSignal  += signal;
        totalMatched += bg + signal;
        seededfile.close();
        std::cout << "Amount of signal hits: " << signal << ", amount of BG hits: " << bg << std::endl;
        //-----------------------------------------------------------------------------------
    } // event loop 
    f->Delete();
    std::cout << "Number of events considered:        " << totalEvents << "\n"
              << "Number of matched triplets:         " << totalMatched << "\n"
              << "Number of matched signal seeds:     " << totalSignal << "\n"
              << "Number of matched background seeds: " << totalBG << std::endl;

}
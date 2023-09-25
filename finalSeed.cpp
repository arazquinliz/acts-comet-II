// ---------------------------------------------------------------- //
//                Creator: Amaia Razquin Lizarraga                  //
//                Contact: amaiafisika@gmail.com                    //
//                Date:    7 September 2023                         //
// ---------------------------------------------------------------- //
/* 
    This program performs a realistic track seeding procedure. It outputs 
    (a) seedEfficiency.txt with the event number, the accuracy of the found seed
    and the estimated information, alongside the real one,
    (b) seedMap.txt with the spacepoint of each seed and event, and 
    (c) pEstim_seed.pdf with the estimated initial momentum. (a) can be visualised
    with efficiency_plot.c and (b) with mapSeeds.c.
    of the found seeds.

    Command:
        g++ finalSeed.cpp Geometry.cpp MagField.cpp Seeder.cpp KalmanFitting.cpp 
        -o finalSeeding.x
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
#include <Acts/Seeding/EstimateTrackParamsFromSeed.hpp>
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
std::map<int, std::pair<int, Acts::Vector3>> event_map; // event number, seeds

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
    auto *f = TFile::Open("detector_oa_all_95.root");
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

    // ----------- Fitting essentials ---------------------------------------------
    // The step length logger for testing & end of world aborter
    using MaterialInteractor = Acts::MaterialInteractor;
    using SteppingLogger = Acts::detail::SteppingLogger;
    using EndOfWorld = Acts::EndOfWorldReached;

    // Action list and abort list
    using ActionList = Acts::ActionList<SteppingLogger, MaterialInteractor>;
    using AbortList = Acts::AbortList<EndOfWorld>;
    using PropagatorOptions = Acts::DenseStepperPropagatorOptions<ActionList, AbortList>;

    using KalmanUpdater = Acts::GainMatrixUpdater;
    using KalmanSmoother = Acts::GainMatrixSmoother;

    KalmanUpdater kfUpdater;
    KalmanSmoother kfSmoother;
    // ----------------------------------------------------------------------------

    // ----------- Seeding tools --------------------------------------------------
    // Parameters
    double radLperSeed  = 2.8221273651251786; //1.82039471328576; //5.381265264757346;
    double rMax         = 618.4430331494052; //414.44127364553236; //690.1954797042055;
    double deltaRMin    = 0.2552817188931512; //0.7124756581567081; //0.2113334751419763;
    double deltaRMax    = 18.93532091460369; //205.5927616870714; //104.36756509015152;
    double deltaZMax    = 23.31336978526198; //17.84223560583328; //29.233848440510886;
    double binSizeR     = 0.6733079493913607; //0.5836506989813482; //0.6922886787481911;
    double deltaRMinMid = 0.0621793950931239; //0.01683702549520169; //0.15801016669789986;
    double deltaRMaxMid = 0.1487536496743246; //0.03526500801415552; //0.25805540458126863;
    double collMin      = -40601.59963597792; //-30300.88965230594; //-44478.60870408864;
    double collMax      = 4562.1278276630055; //1671.2037499575285; //3395.7357360193287;
    double cotThetaMax  = 3391.796214595956; //1727.4763850279494; //1052.032750551505;
    double minPt        = 140.67336385344095; //135.08640935511912; //128.04907527478224;
    double impactMax    = 1555.7500336299754; //4169.118798266818; //4914.642127968815;
    double sigmaScat    = 5.124755191149501; //4.92640165617069; //2.624471976102091;
    double maxPtScat    = 216.26260179148767; //386.8669313920498; //389.1264111786207;
    // Grid
    double gminPt       = 511.1400453935276; //495.23596016030933; //119.09021636183944;
    double grMax        = 276.689575326454; //269.72204202405516; //570.8036490784622;
    double gdeltaRMax   = 95.92672446228377; //80.76058831602866; //141.52768431029594;
    double gimpactMax   = 95.6080445233896; //24.18607493795143; //98.48661116503268;
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
    std::fstream seedfile, efficiency;
    seedfile.open("seedMap.txt", std::ios::out); 
    efficiency.open("seedEfficiency.txt", std::ios::out);
    int totalFake = 0;
    int totalMatched = 0;
    int eventMatch = 0;
    int totalEvents = 0;
    int bgMatches = 0;
    int total_seeds = 0;

    // Initialize all the histograms
    TH1F *h_pEstim = new TH1F("pEstim", "True signal seeds", 200, 0, 150);
    TH1F *h_pEstim_f = new TH1F("pEstim_l1", "Fake signal seeds", 200, 0, 150);
    TH1F *h_pEstim_bg = new TH1F("pEstim_l1", "True background seeds", 200, 0, 150);
    auto logger = Acts::getDefaultLogger("ParamEstimator", Acts::Logging::VERBOSE);
    std::cout << tree->GetEntries() << std::endl;
    //for (int eventID = 0; eventID < tree->GetEntries(); eventID++) {
    for (int eventID = 0; eventID < 16762; eventID++) {
        tree->GetEntry(eventID);
        std::cout << "------------------- EVENT: " << eventID << " --------------------------" << std::endl;
        std::cout << "Size of event " << eventID << " is " << (*spOfEvent).size() << std::endl;
        if ((*spOfEvent).size() < 28) {std::cout << "Ignoring event" << std::endl; continue;}
        
        // Randomly add time of creation to match BG arrival time
        TRandom3* rndm = new TRandom3();
        rndm->SetSeed(0);
        double creationT = rndm->Uniform(620, 1000); //rndm->Exp(864)

        // ----------- Fitting tools --------------------------------------------------
        // Calibrator
        Acts::CalibrationContext calctx;

        auto propagator = aBuildPropagator(tGeometry, bField);
        std::cout << "ACTS::BuildPropagator: Propagator has been created." << std::endl;

        // Build propagator options
        PropagatorOptions pOptions(*gctx, magctx);
        pOptions.loopProtection = false; // Very important to allow helical trajectories
        pOptions.maxStepSize    = 3_m;
        // Default options -change as needed-
        //pOptions.pathLimit               = std::numeric_limits<double>::max();
        pOptions.absPdgCode              = Acts::PdgParticle::eElectron; // 211 pion by default
        //pOptions.direction               = Acts::NavigationDirection::Forward; 
        pOptions.mass                    = 0.511_MeV;
        //pOptions.maxRungeKuttaStepTrials = 10000;
        //pOptions.targetTolerance         = 0.; // Tolerance to reach target
        //pOptions.loopFraction            = 0.5; // DON'T use
        //pOptions.tolerance               = 1e-4; // For the error integration on the stepper
        //pOptions.stepSizeCutOff          = 0.;
        
        // Switch the material interaction on
        auto& mInteractor = pOptions.actionList.get<MaterialInteractor>();
        mInteractor.multipleScattering = true;
        mInteractor.energyLoss         = true;
        mInteractor.recordInteractions = true;

        //Switch the logger to sterile, e.g. for timing checks 
        auto& sLogger   = pOptions.actionList.get<SteppingLogger>();
        sLogger.sterile = false;

        // Perigee surface as target surface
        auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
            Acts::Vector3{0., 0., -1380.25_mm});

        //-----------------------------------------------------------------------------

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
                std::cout << "Station : " << stationID << " Not enough hits." << std::endl;
                continue;
            }
            else {event_map[eventID * 10 + i] = std::make_pair(0, spVector[0]->p());}
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
            int matched = 0;
            int si = 0;
            int bg = 0;
            for (size_t j = 0; j < seeds.size(); j++) {
                // Create vector and estimate parameters
                std::vector<const SpacePoint*> seedSPA{seeds[j].sp()[0], 
                    seeds[j].sp()[1], seeds[j].sp()[2]};

                std::sort(seedSPA.begin(), seedSPA.end(), [](const SpacePoint* lhs, 
                    const SpacePoint* rhs) {return (lhs->z() < rhs->z());});
                const auto rSurface = tGeometry->findSurface(seedSPA[0]->geoID());
                auto field_res = BFieldmap->getField(Acts::Vector3{seedSPA[0]->x(), 
                    seedSPA[0]->y(), seedSPA[0]->z()}, bCache);
                auto field = field_res.value();
                //auto paramsA = estimateParams(*gctx, tGeometry, magctx, BFieldmap, seedSPA);
                Acts::ActsScalar bFieldMin = 0;
                Acts::ActsScalar eMass = 0.511_MeV;
                // Sort seeds in z-axis
                auto paramsA = (Acts::estimateTrackParamsFromSeed(*gctx, seedSPA.begin(), 
                    seedSPA.end(), *rSurface, field, bFieldMin, *logger, eMass)).value();
                double p_estimA = - 1 / paramsA[Acts::eBoundQOverP] * 1000;
                
                // Sort seeds in z-axis
                std::vector<SpacePoint> seed{*seeds[j].sp()[0], 
                    *seeds[j].sp()[1], *seeds[j].sp()[2]};
                std::sort(seed.begin(), seed.end(), [](const SpacePoint lhs, 
                            const SpacePoint rhs) {return (lhs.z() < rhs.z());});
                
                // Time progression cut
                if ((seed[2].t() - seed[1].t() < - 101.5) ||
                   (seed[1].t() - seed[0].t() < - 101.5)) {continue;}

                int coincidence = 0;
                // Compare seeds
                for (size_t k = (j + 1); k < seeds.size(); k++) {
                    std::vector<const SpacePoint*> seedSPB{seeds[k].sp()[0], 
                        seeds[k].sp()[1], seeds[k].sp()[2]};
                    auto seedComb = compareSeeds(seeds[j], seeds[k]);
                    if (seedComb.size() == 4) {
                        coincidence++;
                    }
                }
                double minC, maxC, minN, maxN;
                if (i == 0) {minC = 60; maxC = 1e6; minN = 40; maxN = 1e6;}
                else if (i == 1 || i == 2 || i == 3) {minC = 80; maxC = 1e6; minN = 70; maxN = 1e6;}
                else if (i == 4) {minC = 20; maxC = 1e6; minN = 0; maxN = 1e6;}
                else if (i == 5) {minC = 10; maxC = 1e6; minN = 5; maxN = 1e6;}
                else if (i == 6) {minC = 0; maxC = 180; minN = 0; maxN = 200;}
                else if (i ==7) {minC = 0; maxC = 20; minN = 0; maxN = 40;}
                if (coincidence < 3) {
                    // Apply momentum window cut
                    if ((abs(p_estimA) < minC) || (abs(p_estimA) > maxC)) {continue;}
                }
                else {
                    // Apply momentum window cut
                    if ((abs(p_estimA) < minN) || (abs(p_estimA) > maxN)) {continue;}
                }

                total_seeds++;
                auto mom = seeds[j].sp()[0]->p();
                double pT = std::sqrt(mom.x() * mom.x() + mom.y() * mom.y());
                double mu = std::atanh(mom.z() / mom.norm());

                // ADD INFORMATION 
                if ((seeds[j].sp()[0]->event() == seeds[j].sp()[1]->event()) && 
                    (seeds[j].sp()[0]->event() == seeds[j].sp()[2]->event())) {
                    matched++;
                    if (seeds[j].sp()[0]->event() == eventID) {
                        si++;
                        h_pEstim->Fill(p_estimA);
                        event_map[eventID*10 + i].first += 1;
                        efficiency << eventID << " " << 0 << " " << 1 << " " << mom.norm() << " "  << pT << " " << mu << "\n";
                    }
                    else {
                        bg++; 
                        h_pEstim_bg->Fill(p_estimA);
                        efficiency << eventID << " " << 1 << " " << 1 << " " << mom.norm() << " "  << pT << " " << mu << "\n";
                    }
                }
                else {
                    totalFake++; 
                    h_pEstim_f->Fill(p_estimA);
                    efficiency << eventID << " " << 2 << " " << 0 << " " << mom.norm() << " "  << pT << " " << mu << "\n";
                }
            } // seedA
            seedfile << eventID << " " << i << " " << si << " " << spVector.size() << " " << seeds.size() << "\n";
            totalMatched += matched;
            bgMatches += bg;
            if (si > 0) {eventMatch++;}  
        } // station
    } // event loop 

    std::cout << event_map.size() << std::endl;
    for (auto const& event : event_map) {
        if (event.second.first == 0) {
            auto mom = event.second.second;
            double pT = std::sqrt(mom.x() * mom.x() + mom.y() * mom.y());
            double mu = std::atanh(mom.z() / mom.norm());
            efficiency << event.first << " " << 0 << " " << 0 << " " << mom.norm() << " " << pT << " " << mu << "\n";
        }
    }
    std::cout << "TOTAL NUMBER OF SEEDS:  " << total_seeds  << std::endl;
    std::cout << "TOTAL NUMBER OF EVENTS: " << totalEvents  << std::endl;
    std::cout << "TOTAL MATCHED EVENTS:   " << eventMatch   << std::endl;
    std::cout << "TOTAL MATCHED SEEDS:    " << totalMatched << std::endl;
    std::cout << "TOTAL BG MATCHED SEEDS: " << bgMatches    << std::endl;
    std::cout << "NUMBER OF FAKE SEEDS:   " << totalFake    << std::endl;
    std::cout << "My job here is done" << std::endl;

    TCanvas *c5 = new TCanvas();
    
    h_pEstim_f->Scale(0.05);
    h_pEstim_bg->Scale(0.05);

	h_pEstim_bg->Draw("HIST");
    h_pEstim_bg->SetLineColor(4);
    h_pEstim->Draw("HIST SAME");
    h_pEstim->SetLineColor(2);
    h_pEstim_f->Draw("HIST SAME");
    h_pEstim_f->SetLineColor(6);
    
    gPad->BuildLegend(0.53,0.67,0.88,0.88, "Estimation of Momentum");
    /*TPaveLabel *t = new TPaveLabel(0.15, 0.9, 0.85, 1.0, "Estimated momentum for found seeds", "brNDC"); 
    t->SetBorderSize(0);
    t->SetFillColor(gStyle->GetTitleFillColor());
    t->Draw();*/
    h_pEstim_bg->SetTitle(0);
    h_pEstim_bg->SetStats(0);
    h_pEstim_bg->GetYaxis()->SetTitle("Entries");
    h_pEstim_bg->GetXaxis()->SetTitle("p_{estim} [MeV/c]");
    c5->Print("histograms/pEstim_seed.pdf");
    f->Delete();
}
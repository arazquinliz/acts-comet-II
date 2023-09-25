// ---------------------------------------------------------------- //
//                Creator: Amaia Razquin Lizarraga                  //
//                Contact: amaiafisika@gmail.com                    //
//                Date:    7 September 2023                         //
// ---------------------------------------------------------------- //
/* 
    This program performs a realistic track finding procedure (CKF). 
    The fitted momentum, original momentum and fit chi2 and ndf are outputed
	in fit_check.txt, which later on can be analysed with fitCheck.c in root.
    A recon_mom.txt file is also outputed that can later be analysed with 
    recon_stats.c. Finally, the initial position, momentum and spacepoints are
    stored in acts_output.root.

    Command:
        g++ finalFind.cpp Geometry.cpp MagField.cpp Measurements.cpp 
        KalmanFitting.cpp -o finalFind.x
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

TFile* out_file;
TTree* out_tree;

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
    auto f = TFile::Open("detector_oa_all_t.root");
	if (!f) {
		std::cout << "Couldn't open measurement file\n";
	}

	auto tree = (TTree*)f->Get("detections");
	std::vector<ROOT::Math::XYZVector> *spOfEvent = 0;
    std::vector<ROOT::Math::XYZVector> *momOfEvent = 0;
	std::vector<double> *timeOfEvent = 0;
    std::vector<double> *eDeposition = 0;
	tree->SetBranchAddress("spOfEvent", &spOfEvent);
    tree->SetBranchAddress("momOfEvent", &momOfEvent);
    tree->SetBranchAddress("timeOfEvent", &timeOfEvent);
    tree->SetBranchAddress("eDeposition", &eDeposition);
    // ----------------------------------------------------------------------------

    // ----- Output ---------------------------------------------------------------
	out_file = TFile::Open("acts_output.root", "recreate");
	out_tree = new TTree("tracks", "tracks");

    TVector3 initMomentum;
	TVector3 initPosition;
	std::vector<ROOT::Math::XYZVector> spacePoints;
	//std::vector<double> time;
	out_tree->Branch("initMomentum", &initMomentum);
	out_tree->Branch("initPosition", &initPosition);
	out_tree->Branch("spacePoints", &spacePoints);
	//out_tree->SetBranchAddress("time", &time); 
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

    // ----------- Event loop -----------------------------------------------------   
    std::fstream reconfile;
    reconfile.open("recon_mom.txt", std::ios::out); 

    // Score function parameters
    int totalReconstructed = 0;
    int totalDuplicated    = 0;
    int totalMatched       = 0;
    int totalFake          = 0;
    int totalEvents        = 0;

    std::fstream qualityfile;
	qualityfile.open("findCheck.txt", std::ios::out); 

    //for (int eventID = 0; eventID < tree->GetEntries(); eventID++) {
    for (int eventID = 0; eventID < 20000; eventID++) {
        tree->GetEntry(eventID);
        std::cout << "Size of event: " << (*spOfEvent).size() << std::endl;
        if ((*spOfEvent).size() < 28) {std::cout << "Ignoring event" << std::endl; continue;}
        totalEvents++;
        // Randomly add time of creation to match BG arrival time
        TRandom3* rndm = new TRandom3();
        rndm->SetSeed(0);
        double creationT = rndm->Uniform(620, 1000);

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

        std::vector<double> momList; // helper vector for initial track parameters
        std::map<int, std::vector<SpacePoint>> spMap; // helper vector for spacepoints
        std::vector<Acts::CurvilinearTrackParameters> initialParameters;
        SourceLinkContainer sourceLinks;
        std::vector<Acts::BoundVariantMeasurement> measurements;
        std::vector<Acts::Vector3> bgSP;
        std::vector<Acts::Vector3> signalSP;

        int bg     = 0;
        int signal = 0;
        
        std::fstream seededfile;
        seededfile.open("seeded_meas.txt", std::ios::out); 
        for (int i = 0; i < 8; i++) {
            Acts::Extent rRangeSPExtent;
            int stationID = i;
            auto spVector = getSpacePoints(*gctx, *spOfEvent, *momOfEvent, *eDeposition, *timeOfEvent, 
                creationT, tracker_layout, eventID, stationID, rRangeSPExtent);

            /*if (spVector.size() < 3) {
                std::cout << "Station : " << stationID << " Not enough hits." << std::endl;
                continue;
            }*/
            std::vector<SpacePoint> seed;
            for (auto& sp : spVector) {
                seed.push_back(*sp);
            }
            if (i == 0) {
                // Create vector and estimate parameters
                if (spVector.size() < 3) {break;}

                // The spacepoint within the seed should be ordered for this to be ok
                std::vector<const SpacePoint*> seedSPA{spVector[0], 
                    spVector[1], spVector[2]};
                auto params = estimateParams(*gctx, tGeometry, magctx, BFieldmap, seedSPA);
                double p_estim = - 1 / params[Acts::eBoundQOverP] * 1000;

                // Initial parameters
                Acts::Vector3 pos = Acts::Vector3{seed[0].x(), seed[0].y(), seed[0].z()};
                double time  = 0.;
                Acts::Vector3 mom = spVector[0]->p();
                int sign = mom.y() > 0 ? 1 : -1;
                double p     = mom.norm() * 0.001;  //std::abs(p_estim)* 0.001; //   
                double phi   = sign * acos(mom.x() / sqrt(mom.x() * mom.x() + mom.y() * mom.y())); // params[Acts::eBoundPhi]; //   
                double theta = acos(mom.z() / sqrt(mom.x() * mom.x() + mom.y() * mom.y() + mom.z() * mom.z())); // params[Acts::eBoundTheta]; // 
                int    q     = -1;
                auto cov     = setDefaultCovariance();
                auto curvPar = Acts::CurvilinearTrackParameters(
                    Acts::VectorHelpers::makeVector4(pos, time), phi, theta, p*1._GeV, q, cov);
                            
                momList.push_back(p);
                initialParameters.push_back(curvPar);
                    std::cout << pos.x() << " " << pos.y() << " " << pos.z() << " " << phi << " " 
                        << theta << " " << p << std::endl;
            }
            
            for (auto& sp : seed) {
                // Get sourcelinks and measurements only for matching seeds
                if (std::find(spMap[sp.event()].begin(), spMap[sp.event()].end(), sp) 
                    == spMap[sp.event()].end()) {
                    ActsSourceLink::Index index = measurements.size();
                    ActsSourceLink sl(sp.geoID(), index);
                    sourceLinks.emplace(sl.geometryId(), std::move(sl));
                    spMap[sp.event()].push_back(sp);
                    if (sp.event() > tree->GetEntries()) {
                        bgSP.push_back(Acts::Vector3{sp.x(), sp.y(), sp.z()});
                        bg++;
                    }
                    else {
                        seededfile << sp.x() << " " << sp.y() << " " << sp.z() << "\n"; 
                        signalSP.push_back(Acts::Vector3{sp.x(), sp.y(), sp.z()});
                        signal++;
                    }

                    // TODO: 
                    // eBoundTime should be added to the measurements
                    // if we consider realistic straws, there should only be
                    // eBoundLoc0 (with time) representing the drift distance
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
                    cova(Acts::eBoundLoc0, Acts::eBoundLoc0) = 0.15_mm; 
                    cova(Acts::eBoundLoc1, Acts::eBoundLoc1) = 0.15_mm;
                
                    Acts::Measurement<Acts::BoundIndices, 2> meas(Acts::SourceLink{sl}, indices, loca, cova);
                    measurements.push_back(meas);

                    //std::cout << sp.x() << " " << sp.y() << " " << sp.z() << " " << sp.event() << std::endl; 
                } // spacepoint
            } // seedA
        } // station
        seededfile.close();
        std::cout << "Amount of signal hits: " << signal << ", amount of BG hits: " << bg << std::endl;
        //-----------------------------------------------------------------------------------

        //----------- Find Tracsks ----------------------------------------------------------
        Calibrator calibrator{measurements};
        
        // Measurement Selector
        // One measurement per surface, no chi2 cut
        Acts::MeasurementSelector::Config measurementSelectorConfig = {{Acts::GeometryIdentifier(),
            {{}, {std::numeric_limits<double>::max()}, {1u}}},};
        // Alternative selector where
        std::vector<double> chi2Max = {700}; // Maximum chi2
        std::vector<double> etaBins = {2}; // bins in |eta| to specify variable selections 
        std::vector<size_t> nMax = {1}; // Maximum number of measurement candidates on a surface
        //Acts::MeasurementSelector::Config measurementSelectorConfig = {{Acts::GeometryIdentifier(),
            //{etaBins, chi2Max, {nMax.begin(), nMax.end()}}},};
        Acts::MeasurementSelector measSel{measurementSelectorConfig};
        
        // Combinatorial Kalman Filter 
        // Extensions
        Acts::CombinatorialKalmanFilterExtensions<Acts::VectorMultiTrajectory> extensions;
        extensions.calibrator.connect<&Calibrator::calibrate>(&calibrator);
        extensions.updater.connect<&KalmanUpdater::operator()<Acts::VectorMultiTrajectory>>(
            &kfUpdater);
        extensions.smoother.connect<&KalmanSmoother::operator()<Acts::VectorMultiTrajectory>>(
            &kfSmoother);
        extensions.measurementSelector.connect<&Acts::MeasurementSelector::select<Acts::VectorMultiTrajectory>>(
            &measSel);
        // Sourcelink Accessor
        SourceLinkAccessor slAccessor;
        slAccessor.container = &sourceLinks;
        Acts::SourceLinkAccessorDelegate<SourceLinkAccessor::Iterator> slAccessorDelegate;
        slAccessorDelegate.connect<&SourceLinkAccessor::range>(&slAccessor);
        // Options
        auto ckfOptions = Acts::CombinatorialKalmanFilterOptions<SourceLinkAccessor::Iterator, 
            Acts::VectorMultiTrajectory>(*gctx, magctx, calctx, slAccessorDelegate, 
            extensions, pOptions, &(*pSurface));
        // Finder
        auto ckFinder = Acts::CombinatorialKalmanFilter<Acts::Propagator<Acts::EigenStepper<>, 
            Acts::Navigator>, Acts::VectorMultiTrajectory>(std::move(propagator),
            Acts::getDefaultLogger("CKF", Acts::Logging::INFO));
        
        //-----------------------------------------------------------------------------------
        // Track containers
        auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
        auto trackStateContainer = std::make_shared<Acts::VectorMultiTrajectory>();
        Acts::TrackContainer tracks(trackContainer, trackStateContainer);
        tracks.addColumn<unsigned int>("trackGroup");
        Acts::TrackAccessor<unsigned int> seedNumber("trackGroup");
        std::vector<std::pair<std::vector<double>, std::vector<Acts::Vector3>>> tracksInEvent;
        // p, 1. * ndf, 1. * chi2, init_mom, numHits, purity, mainEv positionVec
        
        // Find tracks
        unsigned int nSeed = 0;
        int nTotalSeeds    = 0;
        int nFailedSeeds   = 0;

        std::fstream predictedfile;
        predictedfile.open("predicted_meas.txt", std::ios::out); 
        for (std::size_t iSeed = 0; iSeed < initialParameters.size(); ++iSeed) {
            //auto parsee = initialParameters.at(iSeed);
            //std::cout << parsee[Acts::eBound QOverP] 
            std::cout << "We are in the seed: " << iSeed << std::endl;
            auto result = ckFinder.findTracks(initialParameters.at(iSeed), ckfOptions, tracks);
            nTotalSeeds++;
            nSeed++;
            
            auto init = initialParameters.at(iSeed);
            double init_mom = init.absoluteMomentum() * 1000;
            if (!result.ok()) {
                nFailedSeeds++;
                std::cout << "Track finding failed for seed " << iSeed << " with error" << 
                    result.error() << std::endl;
                qualityfile << eventID << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << "\n";
                continue;
            }

            auto& tracksForSeed = result.value();
            for (auto& track : tracksForSeed) {
                std::cout << "Looping over finished tracks" << std::endl;
                int ndf, chi2;
                double p;
                spacePoints.clear();
                if (track.hasReferenceSurface()) {
                    const auto& params     = track.parameters().transpose();
                    const auto& rSurface   = track.referenceSurface();
                    const auto& covariance = track.covariance();
                    Acts::Vector3 globalPos = rSurface.localToGlobal(*gctx, 
                    Acts::Vector2(params[0], params[1]), Acts::Vector3::Zero());
                    p = - 1 / params[4] * 1000;
                    double phi = params[Acts::eBoundPhi];
                    double theta = params[Acts::eBoundTheta];
                    
                            std::cout << "Fitted parameters for track\n" <<
                    "Position: " << globalPos.x() << " " << globalPos.y() << " " << globalPos.z() << "\n" << 
                    "Direction: " << "Phi: " << params[2] << " Theta: " << params[3] << "\n" <<
                    "Momentum: " << - 1 / params[4] * 1000 << " MeV\n" << 
                    "For TrackTip == " << track.tipIndex() << std::endl;

                    initPosition = TVector3(globalPos.x(), globalPos.y(), globalPos.z());
                    initMomentum = TVector3(p * cos(theta) * sin(phi), p * sin(theta) * sin(phi), p * cos(phi));
                }
                seedNumber(track) = nSeed;

                std::vector<Acts::Vector3> positionVec;
                size_t numHits = 0u;
                for (const auto trackState : track.trackStates()) {
                    numHits += 1u;
                    const auto& traj = trackState.trajectory();

                    // Check MultiTrajectoryTests.cpp for all the INFOrmation in trajState
                    const auto& trajState = traj.getTrackState(trackState.index());
                    const auto& params    = trajState.smoothed().transpose(); // predicted, filtered or smoothed
                    const auto& rSurface  = trajState.referenceSurface();
                    Acts::Vector3 globalPos = rSurface.localToGlobal(*gctx, 
                        Acts::Vector2(params[0], params[1]), Acts::Vector3::Zero());

                    predictedfile << globalPos.x() << " " << globalPos.y() << " " << globalPos.z() << "\n";
                    ROOT::Math::XYZVector sps{globalPos.x(), globalPos.y(), globalPos.z()};
                    spacePoints.push_back(sps);
                    std::cout << sps.X() << " " << sps.Y() << " " << sps.Z() << std::endl;

                    if (trajState.hasCalibrated()) {
                        const auto& meas      = trajState.calibrated<2>().transpose();
                        Acts::Vector3 measPos = rSurface.localToGlobal(*gctx, 
                            Acts::Vector2(meas[0], meas[1]), Acts::Vector3::Zero());
                        positionVec.push_back(measPos);
                    }
                    
                    if (numHits == 1u) {
                        auto lastTrajState = Acts::MultiTrajectoryHelpers::trajectoryState(traj, trackState.index());
                        std::cout << "Fitting NDF is: " << lastTrajState.NDF  << 
                            " with Chi squared: " << lastTrajState.chi2Sum << std::endl;
                        ndf  = lastTrajState.NDF;
                        chi2 = lastTrajState.chi2Sum;
                    }
                } // globalPos
                std::cout << numHits << std::endl;
                if (numHits > 15) {
                    std::vector<double> vec = {p, 1. * ndf, 1. * chi2, init_mom, numHits};
                    tracksInEvent.emplace_back(vec, positionVec);
                    out_tree->Fill();
                }
                else {std::cout << "Not enough hits in reconstructed track" << std::endl;}
            } // tracks for seed
        } // seeds
        predictedfile.close();
        std::cout << "Number of initial seeds: " << nTotalSeeds << ", number of failed seeds: " << nFailedSeeds << std::endl; 
    
        // Score function parameters
        double trackPurity;
        int nDuplicates = 0;
        int nMatched    = 0;
        int nFake       = 0;
        totalReconstructed += tracksInEvent.size();
        for (int it = 0; it < tracksInEvent.size(); it++) {
            std::map<int, int> purityMap;
            int nHits = 0;
            int nBG   = 0;
            int nSignal = 0;
            if (tracksInEvent[it].second.size() == 0) {continue;}
            // Make purity map
            for (auto sp : tracksInEvent[it].second) {
                nHits++;
                for (auto const& [event, spVec] : spMap) {
                    if (isIn(spVec, sp)) {
                        purityMap[event] += 1;
                    }
                }
            }
            // Get purity from map
            int mainEv = 0;
            int mainN = 0;
            for (auto const& [event, N] : purityMap) {
                if (N > mainN) {
                    mainEv = event;
                    mainN  = N;
                }
            }
            trackPurity = mainN * 1.0 / nHits;
            tracksInEvent[it].first.push_back(trackPurity);
            tracksInEvent[it].first.push_back(mainEv);
            std::cout << "Track purity is: " << trackPurity << " for particle " << mainEv << std::endl;
        }

        for (auto track : tracksInEvent) {
            if (track.first[4] < 0.5) {nFake++;}
        }

        // Iteratively choose the track with most measurements
        for (int itA = 0; itA < tracksInEvent.size(); itA++) {
            int coin = 0;
            for (int itB = 0; itB < tracksInEvent.size(); itB++) {
                if (itA == itB) {continue;}
                for (auto posA : tracksInEvent[itA].second) {
                    for (auto posB : tracksInEvent[itB].second) {
                        if (posA == posB) {coin++;}
                    }
                }
                if (coin > 1) {
                    int remo = tracksInEvent[itA].second.size() > tracksInEvent[itB].second.size() ? itB : itA;
                    tracksInEvent.erase(tracksInEvent.begin() + remo);
                    nDuplicates++;
                    break;
                }
            }
        }
        
        for (auto track : tracksInEvent) {
            reconfile << track.first[0] << " " << track.first[1] << " " << 
            track.first[2] << " " << track.first[3] << " " << track.first[4] << 
            " "  << track.first[5] << " " << track.first[6] << "\n";
            qualityfile << track.first[6] << " " << 1 << " " << track.first[1] << " " << track.first[2] << " " << track.first[0] << " " << track.first[3] << " " << track.first[5] << "\n";
            if ((track.first[6] == eventID) && (track.first[4] >= 0.5)) {nMatched++;}
        }
        std::cout << "MATCHED EVENTS: " << nMatched << std::endl;
        totalDuplicated += nDuplicates;
        totalMatched    += nMatched;
        totalFake       += nFake;
        
    } // event loop 
    reconfile.close();
    out_file->Write();
    out_tree->Print();
    out_file->Close();
    f->Delete();

    int k = 50;
    double efficiency     = totalMatched * 1.0 / totalEvents;
    double fakeRate       = totalFake * 1.0 / totalReconstructed;
    double duplicatedRate = totalDuplicated * 1.0 / (totalReconstructed - totalFake);
    double scoreFunction  = efficiency - (fakeRate);

    std::cout << scoreFunction << " " << efficiency << " " << fakeRate << " " << duplicatedRate << std::endl;
}
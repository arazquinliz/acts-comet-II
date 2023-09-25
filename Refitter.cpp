#include <fstream>

// Definitions & utilities headers
#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/Utilities/CalibrationContext.hpp>

// Geometry headers
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/Geometry/GeometryHierarchyMap.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>

// Magnetic field headers
#include <Acts/MagneticField/ConstantBField.hpp>

// Event data headers
#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/EventData/SourceLink.hpp>
#include <Acts/EventData/VectorMultiTrajectory.hpp>
#include <Acts/EventData/VectorTrackContainer.hpp>

// Propagator headers
#include <Acts/Propagator/Propagator.hpp>
#include <Acts/Propagator/EigenStepper.hpp>
#include <Acts/Propagator/Navigator.hpp>
#include <Acts/Propagator/MaterialInteractor.hpp>
#include <Acts/Propagator/StandardAborters.hpp>
#include <Acts/Propagator/detail/SteppingLogger.hpp>

// Track fitting headers
#include <Acts/TrackFitting/GainMatrixSmoother.hpp>
#include <Acts/TrackFitting/GainMatrixUpdater.hpp>
#include <Acts/TrackFitting/KalmanFitter.hpp>
#include <Acts/TrackFinding/CombinatorialKalmanFilter.hpp>
#include <Acts/TrackFinding/MeasurementSelector.hpp>

#include "Seeder.hpp"
#include "Geometry.hpp"
#include "MagField.hpp"
#include "KalmanFitting.hpp"

using namespace Acts::UnitLiterals;

std::vector<std::pair<std::vector<double>, std::vector<Acts::Vector3>>> refitting(
    std::vector<double> chi2List,
    std::vector<Acts::BoundVariantMeasurement> measurements,
    std::vector<Acts::CurvilinearTrackParameters> initialParameters) {

    int it = 0;
    std::vector<Acts::Vector2> spVec;
    std::vector<Acts::BoundVariantMeasurement> measures;
    SourceLinkContainer sourceLinks;

    std::array<Acts::BoundIndices, 2> indices;
    indices[0] = Acts::BoundIndices::eBoundLoc0;
    indices[1] = Acts::BoundIndices::eBoundLoc1;

    Acts::ActsSymMatrix<2> cova = Acts::ActsSymMatrix<2>::Zero();
    cova(Acts::eBoundLoc0, Acts::eBoundLoc0) = 0.15_mm; 
    cova(Acts::eBoundLoc1, Acts::eBoundLoc1) = 0.15_mm;
    for (double chi2 : chi2List) {
        if (chi2 < 100) {
            const auto& meas = std::get<Acts::Measurement<Acts::BoundIndices, 2>>(measurements[it]);
            auto sp = meas.parameters();
            auto sl = meas.sourceLink();
            if (!(isIn(spVec, sp))) {
                ActsSourceLink::Index index = measures.size();
				ActsSourceLink source(sl.geometryId(), index);
                Acts::Measurement<Acts::BoundIndices, 2> meas_og(Acts::SourceLink{source}, 
                            indices, sp, cova);
                spVec.push_back(sp);
                sourceLinks.emplace(sl.geometryId(), std::move(source));
                measures.push_back(meas_og);
            }
        }
        it++;
    }
    std::cout << "NEW MEASUREMENT LIST: " << measures.size() << std::endl;

    // Initialize geometry context
	Acts::GeometryContext* gctx;
    
    // Build geometry
    std::shared_ptr<const Acts::TrackingGeometry> tGeometry{nullptr};
    tGeometry = aBuildDiskGeometry_no(*gctx, false);
	if (!tGeometry) {
    		throw std::runtime_error("ACTS::BuilGeometry: tGeometry is null");
    }

    // Initialize magnetic field context
	Acts::MagneticFieldContext magctx;

	// Get magnetic field
	auto mapBins = [](std::array<size_t, 3> bins, std::array<size_t, 3> sizes) {
        	return (bins[0] * (sizes[1] * sizes[2]) + bins[1] * sizes[2] + bins[2]);
      	};

	auto BFieldmap = aGetMagneticField(false, magctx, std::move(mapBins));
	auto bField = std::dynamic_pointer_cast<const Acts::InterpolatedMagneticField>(BFieldmap);
    if (!bField) {
    	throw std::runtime_error("ACTS::GetMagneticField: bField is null");
    }  

    // ----------- Fitting tools --------------------------------------------------
    // Calibrator
    Acts::CalibrationContext calctx;
    // Fitting essentials
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

	auto propagator = aBuildPropagator(tGeometry, bField);
	std::cout << "ACTS::BuildPropagator: Propagator has been created." << std::endl;

    	// Build propagator options
	PropagatorOptions pOptions(*gctx, magctx);
	pOptions.loopProtection = false; // Very important to allow helical trajectories
	pOptions.maxStepSize    = 3_m;
	// Default options -change as needed-
	//pOptions.pathLimit               = std::numeric_limits<double>::max();
	pOptions.absPdgCode              = 11; // 211 pion by default
	//pOptions.direction               = Acts::NavigationDirection::Forward; 
	pOptions.mass                    = 0.510_MeV; // 1000 by default
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

    //----------- Find Tracsks ----------------------------------------------------------
    Calibrator calibrator{measures};
        
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
    for (int iSeed = 0; iSeed < initialParameters.size(); iSeed++) {
        //auto fitRes = kFitter.fit(wrappedSLs.begin(), wrappedSLs.end(), startParameters, 
            //kfOptions, tracks);
        auto result = ckFinder.findTracks(initialParameters.at(iSeed), ckfOptions, tracks);
        if (!result.ok()) {
            std::cout << "Track finding failed for seed " << iSeed << " with error" << 
            result.error() << std::endl;
            continue;
        }
        auto& tracksForSeed = result.value();
        double p, ndf, chi2;
        auto init = initialParameters.at(iSeed);
        double init_mom   = init.absoluteMomentum() * 1000;
        for (auto& track : tracksForSeed) {
            if (track.hasReferenceSurface()) {
                const auto& params     = track.parameters().transpose();
                p = - 1 / params[4] * 1000;
                const auto& rSurface   = track.referenceSurface();
                const auto& covariance = track.covariance();
                Acts::Vector3 globalPos = rSurface.localToGlobal(*gctx, 
                    Acts::Vector2(params[0], params[1]), Acts::Vector3::Zero());
                double phi = params[Acts::eBoundPhi];
                double theta = params[Acts::eBoundTheta];
                //initPosition = TVector3(globalPos.x(), globalPos.y(), globalPos.z());
                //initMomentum = TVector3(p * cos(theta) * sin(phi), p * sin(theta) * sin(phi), p * cos(phi));
            }
                
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
            } // track state
            if (numHits > 15) {
                    std::vector<double> vec = {p, 1. * ndf, 1. * chi2, init_mom, numHits};
                    tracksInEvent.emplace_back(vec, positionVec);
            }
        } // tracks for seed
    }
    return tracksInEvent;
}

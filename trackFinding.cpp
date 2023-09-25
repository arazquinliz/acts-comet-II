// ---------------------------------------------------------------- //
//                Creator: Amaia Razquin Lizarraga                  //
//                Contact: amaiafisika@gmail.com                    //
//                Date:    7 September 2023                         //
// ---------------------------------------------------------------- //
/* 
    This program performs track finding with a set of input measurements
	coming from Measurements.cpp. It contains the standard basic structure
	for running a CKF. 

    Command:
        g++ trackFinding.cpp Geometry.cpp MagField.cpp Measurements.cpp 
        KalmanFitting.cpp -o finding.x
        -L/path/to/acts_install/lib -lActsPluginTGeo -lActsPluginJson -lActsCore 
        -I/path/to/acts_install/include/ -I$(root-config --incdir) 
        -I$(root-config --evelibs) $(root-config --libs) 
        -lGeom -I/usr/include/eigen3
*/

#include <iostream>
#include <ostream>
#include <vector>
#include <map>

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/Units.hpp>

#include <Acts/Utilities/CalibrationContext.hpp>

#include <Acts/EventData/VectorTrackContainer.hpp>

#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/Geometry/GeometryContext.hpp>

#include <Acts/Surfaces/PerigeeSurface.hpp>

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

// Amaia:
#include "base/ActsSourceLink.hpp"
#include "base/Calibrator.hpp"
#include "Geometry.hpp"
#include "MagField.hpp"
#include "Measurements.hpp"
#include "KalmanFitting.hpp"

using namespace Acts::UnitLiterals;

extern std::map<std::string, std::vector<std::shared_ptr<const Acts::Surface>>> tracker_layout;
int main() {
    std::cout << "ACTS:: Hello world: I am alive." << std::endl;

	bool debug = false; // Extensive output or not
	
	// Initialize geometry context
	Acts::GeometryContext* gctx;
    
    // Build geometry
    std::shared_ptr<const Acts::TrackingGeometry> tGeometry{nullptr};
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
	std::cout << "ACTS::GetMagneticField: Magnetic field has been interpolated from root file." << std::endl;
	//auto bField = std::make_shared<Acts::ConstantBField>(Acts::Vector3{0, 0, 1_T});

	//-----------------------------------------------------------------------------------
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
	//pOptions.absPdgCode              = Acts::PdgParticle::eElectron; // 211 pion by default
	//pOptions.direction               = Acts::NavigationDirection::Forward; 
	//pOptions.mass                    = 1000;
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

	//-----------------------------------------------------------------------------------
	// Get input data
	// Get measurements and sourcelinks
	auto event_vector = aCreateDiskLinks(
		tracker_layout, tGeometry, *gctx, debug);

	std::vector<Acts::CurvilinearTrackParameters> initialParameters;
	SourceLinkContainer sourceLinks;
	std::vector<Acts::BoundVariantMeasurement> measurements;

	for (int i = 0; i < 4; i++) {
		std::vector<Acts::BoundVariantMeasurement>& measurements_ev = event_vector[i].measurements;
		std::vector<Acts::SourceLink>& sourceLinkVec = event_vector[i].SourceLinks;
		Acts::Vector3 pos = event_vector[i].initialPos;
		Acts::Vector3 mom = event_vector[i].initialMom;

		if (measurements_ev.size() < 1) {
			continue;
		} 
		std::cout << ".............Event " << i << "........................\n"	
			<< "Measurments size is: " << measurements_ev.size() << "\n"
			<< "SourceLink vector size is: " << sourceLinkVec.size() << std::endl;
				
		for (int it = 0; it < measurements_ev.size(); it++) {
			const auto& meas = std::get<Acts::Measurement<Acts::BoundIndices, 2>>(measurements_ev[it]);
			const auto& sl = meas.sourceLink();
			auto geomID = sl.geometryId();
			auto params = meas.parameters();
			auto cov    = meas.covariance();
			std::cout << geomID << " " << params[0] << " " << params[1] << std::endl;
		}

		for (auto& sl : sourceLinkVec) {
			sourceLinks.emplace(sl.geometryId(), std::move(sl));
		}

		for (auto& meas : measurements_ev) {
			measurements.push_back(meas);
		}

		// Initial parameters (from seeds to track parameters)
		double time  = 0.;
		double p     = mom.norm() * 0.001; 
		double phi   = atan(mom.y() / mom.x());
		double theta = atan(mom.z() / mom.x());
		int    q     = -1;
		auto cov     = setDefaultCovariance();

		std::cout << pos.x() << " " << pos.y() << " " << pos.z() << " " << phi << " " << theta <<
			" " << p << std::endl;

		initialParameters.push_back(Acts::CurvilinearTrackParameters(
			Acts::VectorHelpers::makeVector4(pos, time), phi, theta, p, q, cov));
	}

    // Perigee surface as target surface
	auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
      Acts::Vector3{0., 0., -1380.25_mm});

	//-----------------------------------------------------------------------------------
    // Calibrator
	Acts::CalibrationContext calctx;
	Calibrator calibrator{measurements};

    // Measurement Selector
	// One measurement per surface, no chi2 cut
	Acts::MeasurementSelector::Config measurementSelectorConfig = {{Acts::GeometryIdentifier(),
		{{}, {std::numeric_limits<double>::max()}, {1u}}},};
	// Alternative selector where
	// chi2Max = Maximum chi2 
	// nMax = Maximum number of measurement candidates on a surface
	// etaBins = bins in |eta| to specify variable selections
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

	// Find tracks
	unsigned int nSeed = 0;
	int nTotalSeeds    = 0;
	int nFailedSeeds   = 0;
	
	for (std::size_t iSeed = 0; iSeed < initialParameters.size(); ++iSeed) {
		auto result = ckFinder.findTracks(initialParameters.at(iSeed), ckfOptions, tracks);
		nTotalSeeds++;
		nSeed++;

		if (!result.ok()) {
			nFailedSeeds++;
			std::cout << "Track finding failed for seed " << iSeed << " with error" << 
				result.error() << std::endl;
			continue;
		}

		auto& tracksForSeed = result.value();
		for (auto& track : tracksForSeed) {
			if (track.hasReferenceSurface()) {
				const auto& params     = track.parameters().transpose();
				const auto& rSurface   = track.referenceSurface();
				const auto& covariance = track.covariance();
				Acts::Vector3 globalPos = rSurface.localToGlobal(*gctx, 
				Acts::Vector2(params[0], params[1]), Acts::Vector3::Zero());
						std::cout << "Fitted parameters for track\n" <<
				"Position: " << globalPos.x() << " " << globalPos.y() << " " << globalPos.z() << "\n" << 
				"Direction: " << "Phi: " << params[2] << " Theta: " << params[3] << "\n" <<
				"Momentum: " << - 1 / params[4] * 1000 << " MeV\n" << 
				"For TrackTip == " << track.tipIndex() << std::endl;
			}
			seedNumber(track) = nSeed;
		}
	}

    // Statistics etc

	return 0;
}
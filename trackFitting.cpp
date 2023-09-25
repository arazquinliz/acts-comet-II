// ---------------------------------------------------------------- //
//                Creator: Amaia Razquin Lizarraga                  //
//                Contact: amaiafisika@gmail.com                    //
//                Date:    7 September 2023                         //
// ---------------------------------------------------------------- //
/* 
    This program performs track fitting with a set of input measurements
	coming from Measurements.cpp. It contains the standard basic structure
	for running a KF with known initial seeds and in a single trajectory.
	The fitted momentum, original momentum and fit chi2 and ndf are outputed
	in fit_check.txt, which later on can be analysed with fitCheck.c in root.
	The predicted measurements are printed in predicted_meas.txt so that they 
	can be compared with the real measurements in plot_fitting.c.

    Command:
        g++ trackFitting.cpp Geometry.cpp MagField.cpp Measurements.cpp 
        KalmanFitting.cpp -o fitting.x
        -L/path/to/acts_install/lib -lActsPluginTGeo -lActsPluginJson -lActsCore 
        -I/path/to/acts_install/include/ -I$(root-config --incdir) 
        -I$(root-config --evelibs) $(root-config --libs) 
        -lGeom -I/usr/include/eigen3
*/

#include <iostream>
#include <ostream>
#include <vector>
#include <map>

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

// Visualization
#include <Acts/Visualization/ObjVisualization3D.hpp>
#include <Acts/Visualization/PlyVisualization3D.hpp>
#include <Acts/Visualization/GeometryView3D.hpp>
#include <Acts/Visualization/ViewConfig.hpp>

// SVG
#include <Acts/Plugins/ActSVG/LayerSvgConverter.hpp>
#include <Acts/Plugins/ActSVG/SvgUtils.hpp>
#include <Acts/Plugins/ActSVG/TrackingGeometrySvgConverter.hpp>

// Amaia
#include "Geometry.hpp"
#include "MagField.hpp"
#include "Measurements.hpp"
#include "Logger.hpp"
#include "KalmanFitting.hpp"
#include "Seeder.hpp"
// Base classes
#include "base/ActsSourceLink.hpp"
#include "base/Calibrator.hpp"
#include "base/SpacePoint.hpp"

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

	aLogGeometry(tGeometry, *gctx);
	if (debug){
		aSurfaceStream(tGeometry, *gctx);
	}

	// Visualization	
	double outputScalor = 1.0;
  	size_t outputPrecision = 6;
  	Acts::ObjVisualization3D objVis(outputPrecision, outputScalor);
	Acts::PlyVisualization3D plyVis;

    Acts::ViewConfig containerView = Acts::ViewConfig({220, 220, 220});
  	Acts::ViewConfig volumeView = Acts::ViewConfig({220, 220, 0});
  	Acts::ViewConfig sensitiveView = Acts::ViewConfig({0, 180, 240});
  	Acts::ViewConfig passiveView = Acts::ViewConfig({240, 280, 0});
  	Acts::ViewConfig gridView = Acts::ViewConfig({220, 0, 0});

  	Acts::GeometryView3D::drawTrackingVolume(
      	objVis, *(tGeometry->highestTrackingVolume()), *gctx, containerView,
      	volumeView, passiveView, sensitiveView, gridView, true, ".obj", "./geometry_files");
	std::cout << "ACTS:: 3D OBJ Visualization of geometry has been dumped out." << std::endl;

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
	std::cout << "ACTS:: Magnetic field has been interpolated from root file." << std::endl;
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
	pOptions.absPdgCode              = Acts::PdgParticle::eElectron; // 211 pion by default
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

	// Build Kalman Fitter extensions
    	Acts::KalmanFitterExtensions<Acts::VectorMultiTrajectory> extensions;
	extensions.updater.connect<&KalmanUpdater::operator()<Acts::VectorMultiTrajectory>>(
		&kfUpdater);
	extensions.smoother.connect<&KalmanSmoother::operator()<Acts::VectorMultiTrajectory>>(
		&kfSmoother);
	//-----------------------------------------------------------------------------------
	// Get measurements and sourcelinks
	auto event_vector = aCreateDiskLinks(
		tracker_layout, tGeometry, *gctx, debug);

	std::fstream qualityfile;
	qualityfile.open("fit_check.txt", std::ios::out); 
	for (size_t i = 0; i < event_vector.size(); i++) {
		std::vector<Acts::BoundVariantMeasurement> measurements = event_vector[i].measurements;
		std::vector<Acts::SourceLink> sourceLinks = event_vector[i].SourceLinks;
		Acts::Vector3 pos = event_vector[i].initialPos;
		Acts::Vector3 mom = event_vector[i].initialMom;

		if (measurements.size() < 1) {
	    		throw std::runtime_error("ACTS::CreateLinks: Measurements are null");
	    	} 
		if (measurements.size() < 28) {continue;}
		std::cout << ".............Event " << i << "........................\n"	
			<< "Measurments size is: " << measurements.size() << "\n"
			<< "SourceLink vector size is: " << sourceLinks.size() << std::endl;
			
		std::vector<const SpacePoint*> seed;
		for (int it = 0; it < measurements.size(); it++) {
			const auto& meas = std::get<Acts::Measurement<Acts::BoundIndices, 2>>(measurements[it]);
			const auto& sl = meas.sourceLink();
			auto geomID = sl.geometryId();
			auto params = meas.parameters();
			auto cov    = meas.covariance();
			std::cout << geomID << " " << params[0] << " " << params[1] << std::endl;
			if (it == 0) {
				auto surface    = tGeometry->findSurface(geomID);
				auto surf = new std::shared_ptr<const Acts::Surface> (tGeometry->findSurface(geomID));
				auto posV = surface->localToGlobal(*gctx, 
					Acts::Vector2(params[0], params[1]), Acts::Vector3::Zero());
				auto sp = makeSpacePoint(*gctx,
    				Acts::Vector3{posV.x() + 1080_cm, posV.y(), posV.z() + 333.45_cm}, 0, 
    				*surf, i, Acts::Vector3{0, 0, 0});
				std::cout << sp->x() << " " << sp->y() << " " << sp->z() << std::endl;
				const SpacePoint* sp0 = new SpacePoint{sp->x(), sp->y(), sp->z(), sp->r(), 0, 
					geomID, sp->varianceR(), sp->varianceZ(), i, Acts::Vector3{0, 0, 0}};
				seed.push_back(sp0);
			}
			if (it == 1) {
				auto surface    = tGeometry->findSurface(geomID);
				auto surf = new std::shared_ptr<const Acts::Surface> (tGeometry->findSurface(geomID));
				auto posV = surface->localToGlobal(*gctx, 
					Acts::Vector2(params[0], params[1]), Acts::Vector3::Zero());
				auto sp = makeSpacePoint(*gctx,
    				Acts::Vector3{posV.x() + 1080_cm, posV.y(), posV.z() + 333.45_cm}, 0, 
    				*surf, i, Acts::Vector3{0, 0, 0});
				std::cout << sp->x() << " " << sp->y() << " " << sp->z() << std::endl;
				const SpacePoint* sp0 = new SpacePoint{sp->x(), sp->y(), sp->z(), sp->r(), 0, 
					geomID, sp->varianceR(), sp->varianceZ(), i, Acts::Vector3{0, 0, 0}};
				seed.push_back(sp0);
			}
			if (it == 2) {
				auto surface    = tGeometry->findSurface(geomID);
				auto surf = new std::shared_ptr<const Acts::Surface> (tGeometry->findSurface(geomID));
				auto posV = surface->localToGlobal(*gctx, 
					Acts::Vector2(params[0], params[1]), Acts::Vector3::Zero());
				auto sp = makeSpacePoint(*gctx,
    				Acts::Vector3{posV.x() + 1080_cm, posV.y(), posV.z() + 333.45_cm}, 0, 
    				*surf, i, Acts::Vector3{0, 0, 0});
				std::cout << sp->x() << " " << sp->y() << " " << sp->z() << std::endl;
				const SpacePoint* sp0 = new SpacePoint{sp->x(), sp->y(), sp->z(), sp->r(), 0, 
					geomID, sp->varianceR(), sp->varianceZ(), i, Acts::Vector3{0, 0, 0}};
				seed.push_back(sp0);
			}
		}
		
		// Sort based on r for bottom, middle, top
		std::sort(seed.begin(), seed.end(), [](const SpacePoint* lhs, 
                            const SpacePoint* rhs) {return (lhs->r() < rhs->r());});
	
		// Get initial parameters internally to check biases
        auto params = estimateParams(*gctx, tGeometry, magctx, BFieldmap, seed);
        double p_estim = - 1 / params[Acts::eBoundQOverP] * 1000;
		std::cout << "Estimated initial momentum is: " << p_estim << std::endl;

		// Acts requires a wrapped vector to get the memory access correct
		std::vector<std::reference_wrapper<const Acts::SourceLink>> wrappedSLs;
		for (const auto& sl : sourceLinks) {wrappedSLs.push_back(std::cref(sl));}
		//-----------------------------------------------------------------------------------
		// Define start parameters 
		//Acts::Vector3 pos(-254.9, -168.013, -1379.02);
		//Acts::Vector3 pos(-255, -167.949, -1378.28);
		double time  = 0.;
		//double p     = 105_MeV;
		//double phi  = 0._degree, theta = 90._degree;
		//double phi   = atan(-2.559 / 5.6);
		//double theta = atan(9.35 / 5.6);
		//double px    = p * cos(theta) * cos(phi);
		//double py    = p * cos(theta) * sin(phi);
		//double pz    = p * sin(theta);
		//double px    = 51.7647_MeV;
		//double py    = -25.4294_MeV;
		//double pz    = 88.3467_MeV;
		//Acts::Vector3 mom(px, py, pz);
		int sign = mom.y() > 0 ? 1 : -1;
		double p     = mom.norm() * 0.001;  // std::abs(p_estim)* 0.001; // 
		double phi   = acos(mom.z() / sqrt(mom.x() * mom.x() + mom.y() * mom.y() + mom.z() * mom.z())); // params[Acts::eBoundPhi]; // 
		double theta = sign * acos(mom.x() / sqrt(mom.x() * mom.x() + mom.y() * mom.y()));//0.142523 0.775937; params[Acts::eBoundTheta]; // 
		int    q     = -1;
		
		auto cov = setDefaultCovariance();
		Acts::CurvilinearTrackParameters startParameters(
			Acts::VectorHelpers::makeVector4(pos, time), phi, theta, p, q, cov);

		std::cout << "Initial position: " << pos.x() << " " << pos.y() << " " << pos.z() <<
			" Initial momentum: " << mom.x() << " " << mom.y() << " " << mom.z() << std::endl; 

		// Construct a perigee surface as the target surface
	    	// position is a vector3  for the moment it should be the first position
	    	auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(pos);
		//-----------------------------------------------------------------------------------
		Acts::CalibrationContext calctx;
		Calibrator calibrator{measurements};
	    	extensions.calibrator.connect<&Calibrator::calibrate>(
		&calibrator);

	    	Acts::KalmanFitterOptions<Acts::VectorMultiTrajectory> kfOptions(*gctx, magctx, calctx,
			extensions, pOptions, &(*pSurface));
	    	// Default options -change as needed-
	    	//kfOptions.multipleScattering                 = true;
	    	//kfOptions.energyLoss                         = true;
	    	kfOptions.reversedFiltering                  = false;
	    	//kfOptions.reversedFilteringCovarianceScaling = 1.0;
	    	//kfOptions.freeToBoundCorrection              = Acts::FreeToBoundCorrection(false);

		auto kFitter = Acts::KalmanFitter<Acts::Propagator<Acts::EigenStepper<>, 
			Acts::Navigator>, Acts::VectorMultiTrajectory>(propagator, 
			Acts::getDefaultLogger("KalmanFitter", Acts::Logging::VERBOSE));
	 
		Acts::TrackContainer tracks{Acts::VectorTrackContainer{},
	      				    Acts::VectorMultiTrajectory{}};

		//tracks.addColumn<bool>("reversed");
	    	//tracks.addColumn<bool>("smoothed");

	    	auto fitRes = kFitter.fit(wrappedSLs.begin(), wrappedSLs.end(), startParameters, 
			kfOptions, tracks);
		//-----------------------------------------------------------------------------------
		// Check that the track fit did not turn and error
		std::unordered_map<Acts::MultiTrajectoryTraits::IndexType, 
			Acts::BoundTrackParameters> indexedParams;
		if (fitRes.ok()) {
			std::cout << "Fitting result is OK." << std::endl;
			const auto& fitOutput  = fitRes.value();
			double momentum;
			int ndf;
			double chi2;
			std::fstream predictedfile;
			predictedfile.open("predicted_meas.txt", std::ios::out); 
			if (fitOutput.hasReferenceSurface()) {
				const auto& params     = fitOutput.parameters().transpose();
				const auto& rSurface   = fitOutput.referenceSurface();
				const auto& covariance = fitOutput.covariance();
				Acts::Vector3 globalPos = rSurface.localToGlobal(*gctx, 
					Acts::Vector2(params[0], params[1]), Acts::Vector3::Zero());
				std::cout << "Fitted parameters for track\n" <<
					"Position: " << globalPos.x() << " " << globalPos.y() << " " << globalPos.z() << "\n" << 
					"Direction: " << "Phi: " << params[2] << " " << phi << " Theta: " << params[3] << " " << theta << "\n" <<
					"Momentum: " << - 1 / params[4] * 1000 << " MeV\n" << 
					"For TrackTip == " << fitOutput.tipIndex() << std::endl;

				momentum = - 1 / params[4] * 1000;
				//predictedfile << globalPos.x() << " " << globalPos.y() << " " << globalPos.z() << "\n";

				indexedParams.emplace(std::pair(fitOutput.tipIndex(), 
					Acts::BoundTrackParameters{rSurface.getSharedPtr(), params, covariance}));
			}

			auto trackplease = tracks.getTrack(0);
			std::cout << "Number of track states: " << trackplease.nTrackStates() << std::endl;
			
			size_t numHits = 0u;
			for (const auto trackState : trackplease.trackStates()) {
				numHits += 1u;
				const auto& traj = trackState.trajectory();

				// Check MultiTrajectoryTests.cpp for all the VERBOSErmation in trajState
				const auto& trajState = traj.getTrackState(trackState.index());
				const auto& params    = trajState.smoothed().transpose(); // predicted, filtered or smoothed
				const auto& rSurface  = trajState.referenceSurface();
				Acts::Vector3 globalPos = rSurface.localToGlobal(*gctx, 
					Acts::Vector2(params[0], params[1]), Acts::Vector3::Zero());

				predictedfile << globalPos.x() << " " << globalPos.y() << " " << globalPos.z() << "\n";
			
				if (numHits == 1u) {
					auto lastTrajState = Acts::MultiTrajectoryHelpers::trajectoryState(traj, trackState.index());
					std::cout << "Fitting NDF is: " << lastTrajState.NDF  << 
						" with Chi squared: " << lastTrajState.chi2Sum << std::endl;
					ndf = lastTrajState.NDF;
					chi2 = lastTrajState.chi2Sum;
				}
			}
			qualityfile << i << " " << 1 << " " << ndf << " " << chi2 << " " << momentum << " " << mom.norm() << " " << p << "\n"; // ok ndf chi2 p
			predictedfile.close();
		}
		else {qualityfile << i << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << "\n";}

	} // event loop
	qualityfile.close();
	
    std::cout << "ACTS:: My job here is done." << std::endl;
	return 0;
}
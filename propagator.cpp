// ---------------------------------------------------------------- //
//                Creator: Amaia Razquin Lizarraga                  //
//                Contact: amaiafisika@gmail.com                    //
//                Date:    7 September 2023                         //
// ---------------------------------------------------------------- //
/* 
    This program serves to test the geometry with random partile propagation. 
    It uses the Writer.cpp module to write out the propagation information and
    the material interactions. It is essentially a combination of PropagationAlgorithm.cpp 
    & PropagatorInterface.hpp already contained in ACTS.
    The output information is stored in propagation_steps.txt, PropagationOutput.root
    and MaterialOutput.root. The propagation steps can be plotted out by using the
    plot_propagation.c program within root.
    Command:
        g++ propagator.cpp Geometry.cpp MagField.cpp Logger.cpp Writer.cpp -o propagator.x 
        -L/path/to/acts_install/lib -lActsPluginTGeo -lActsPluginJson -lActsCore 
        -I/path/to/acts_install/include/ -I$(root-config --incdir) 
        -I$(root-config --evelibs) $(root-config --libs) 
        -lGeom -I/usr/include/eigen3
*/

#include <cmath>
#include <limits>
#include <memory>
#include <optional>
#include <iostream>
#include <ostream>
#include <vector>
#include <map>
#include <random>

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/Utilities/Helpers.hpp>

#include "Acts/EventData/GenericBoundTrackParameters.hpp"
#include <Acts/EventData/TrackParameters.hpp>

#include <Acts/Propagator/AbortList.hpp>
#include <Acts/Propagator/ActionList.hpp>
#include <Acts/Propagator/DenseEnvironmentExtension.hpp>
#include <Acts/Propagator/MaterialInteractor.hpp>
#include <Acts/Propagator/Navigator.hpp>
#include <Acts/Propagator/Propagator.hpp>
#include <Acts/Propagator/StandardAborters.hpp>
#include <Acts/Propagator/detail/SteppingLogger.hpp>
#include <Acts/Propagator/EigenStepper.hpp>

#include <Acts/Surfaces/PerigeeSurface.hpp>

#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/Geometry/GeometryContext.hpp>

#include <Acts/MagneticField/MagneticFieldContext.hpp>
#include <Acts/MagneticField/ConstantBField.hpp>

// Visualization
#include <Acts/Visualization/ObjVisualization3D.hpp>
#include <Acts/Visualization/GeometryView3D.hpp>
#include <Acts/Visualization/ViewConfig.hpp>

#include "Writer.hpp"
#include "Geometry.hpp"
#include "MagField.hpp"
#include "Logger.hpp"

using RecordedMaterial = Acts::MaterialInteractor::result_type;
// start position, momentum and recorded material
using RecordedMaterialTrack = 
    std::pair<std::pair<Acts::Vector3, Acts::Vector3>, RecordedMaterial>;
// Output of the test
using PropagationOutput = std::pair<std::vector<Acts::detail::Step>, 
    RecordedMaterial>;

std::default_random_engine rng(0);

using namespace Acts::UnitLiterals;

std::optional<Acts::BoundSquareMatrix> generateCovariance(
    std::normal_distribution<double>& gauss) {
    
    // Free inputs
    bool covarianceTransport          = false;
    Acts::BoundSquareMatrix correlations = Acts::BoundSquareMatrix::Identity();
    Acts::BoundVector covariances     = Acts::BoundVector::Zero();
    
    if (covarianceTransport) {
        // We start from the correlation matrix
        Acts::BoundSquareMatrix newCov(correlations);
        // Then we draw errors according to the error values
        Acts::BoundVector covs_smeared = covariances;
        for (size_t k = 0; k < size_t(covs_smeared.size()); ++k) {
            covs_smeared[k] *= gauss(rng);
        }

        // and apply a double loop
        for (size_t i = 0; i < size_t(newCov.rows()); ++i) {
            for (size_t j = 0; j < size_t(newCov.cols()); ++j) {
                (newCov)(i, j) *= covs_smeared[i];
                (newCov)(i, j) *= covs_smeared[j];
            }
        }
        return newCov;
    }
    return std::nullopt;
}

int main() {

    bool debug = false;

    // The step length logger for testing & end of world aborter
    using MaterialInteractor = Acts::MaterialInteractor;
    using SteppingLogger = Acts::detail::SteppingLogger;
    using EndOfWorld = Acts::EndOfWorldReached;

    // Action list and abort list
    using ActionList = Acts::ActionList<SteppingLogger, MaterialInteractor>;
    using AbortList = Acts::AbortList<EndOfWorld>;
    using PropagatorOptions =
        Acts::DenseStepperPropagatorOptions<ActionList, AbortList>;

    std::normal_distribution<double> gauss(0., 1.);
    
    std::uniform_real_distribution<double> phiDist(0, M_PI*0.5);
    std::uniform_real_distribution<double> thetaDist(0, M_PI*0.5);
    //std::uniform_real_distribution<double> etaDist(-4, 4);
    //std::uniform_real_distribution<double> ptDist(100_MeV, 0.11_GeV);
    std::uniform_real_distribution<double> pDist(100_MeV, 110_MeV);
    std::uniform_real_distribution<double> qDist(0., 1.);

    // Some free inputs
    int nTests     = 10;
    double d0Sigma = 150_um;
    double z0Sigma = 300_um;
    double tSigma  = 1_ns;

    std::shared_ptr<const Acts::PerigeeSurface> surface = 
        Acts::Surface::makeShared<Acts::PerigeeSurface>(
            Acts::Vector3(0., 0., -1380.));
    
    // Propagation steps to output
    std::vector<std::vector<Acts::detail::Step>> propagationSteps;
    propagationSteps.reserve(nTests);

    // Recorded material to output
    std::unordered_map<size_t, Acts::RecordedMaterialTrack> recordedMaterial;

    // Get geometry
    // Initialize geometry context
	Acts::GeometryContext* gctx;
    
    // Build geometry
    std::shared_ptr<const Acts::TrackingGeometry> tGeometry{nullptr};
    tGeometry = aBuildDiskGeometry(*gctx, debug);

    aLogGeometry(tGeometry, *gctx);
	if (debug){
		aSurfaceStream(tGeometry, *gctx);
	}

    // Set up propagator
    Acts::Navigator::Config ncfg{tGeometry};
    ncfg.resolveMaterial  = true;
    ncfg.resolvePassive   = false;
    ncfg.resolveSensitive = true;
    Acts::Navigator navigator(ncfg, Acts::getDefaultLogger("Navigator",
      	Acts::Logging::Level::INFO));

    // Initialize magnetic field context
	Acts::MagneticFieldContext mctx;
    // Get magnetic field
	auto mapBins = [](std::array<size_t, 3> bins, std::array<size_t, 3> sizes) {
        	return (bins[0] * (sizes[1] * sizes[2]) + bins[1] * sizes[2] + bins[2]);
      	};

	auto BFieldmap = aGetMagneticField(debug, mctx, std::move(mapBins));
	auto bField = std::dynamic_pointer_cast<const Acts::InterpolatedMagneticField>(BFieldmap);
    if (!bField) {
    		throw std::runtime_error("ACTS:: GetMagneticField: bField is null");
    } 
	std::cout << "ACTS:: Magnetic field has been interpolated from root file." << std::endl;
	//auto bField = std::make_shared<Acts::ConstantBField>(Acts::Vector3{0, 0, 1._T});
	Acts::EigenStepper<> stepper(bField);
	Acts::Propagator<Acts::EigenStepper<>, Acts::Navigator> propagator(std::move(stepper), 
		std::move(navigator), Acts::getDefaultLogger("Propagator", Acts::Logging::INFO));

    // Loop over number of particles
    int mat = 0;
    for (size_t it = 0; it < nTests; ++it) {
        double d0     = d0Sigma * gauss(rng);
        double z0     = z0Sigma * gauss(rng);
        double phi    = phiDist(rng);
        double theta  = thetaDist(rng);
        //double eta = etaDist(rng);
        //double theta = 2 * atan(exp(-eta));
        //double pt     = ptDist(rng);
        //double p      = pt / sin(theta);
        double p      = pDist(rng);
        double charge = qDist(rng) > 0.5 ? 1. : -1.;
        double qop    = charge / p;
        double t      = 0; //tSigma * gauss(rng);
        
        // Parameters
        Acts::BoundVector pars;
        pars << d0, z0, phi, theta, qop, t;

        Acts::Vector3 sPosition(0., 0., 0.);
        Acts::Vector3 sMomentum(0., 0., 0.);

        auto cov = generateCovariance(gauss);

        // More free variables
        double pathLength               = std::numeric_limits<double>::max();
        double maxStepSize              = 3_m;
        double ptLoopers                = 500_MeV;
        bool energyLoss                 = true;
        bool multipleScattering         = true;
        bool recordMaterialInteractions = true;
        bool sterileLogger              = false;

    
        PropagationOutput pOutput;

        if (charge != 0.0) {
            // Charged extrapolation - with hit recording
            Acts::BoundTrackParameters startParameters(surface, pars, std::move(cov));
            sPosition = startParameters.position(*gctx);
            sMomentum = startParameters.momentum();
            std::cout << "ACTS:: Test propagation with charged particle starts." << std::endl;
            std::cout << "With particle at position: " << sPosition.x() << " " << sPosition.y()
                << " " << sPosition.z() << std::endl;
            std::cout << "Momentum: " << p << " Charge: " << charge << " Phi: " << phi << " Theta: " << theta << std::endl;

            PropagatorOptions options(*gctx, mctx);
            options.pathLimit = pathLength;
            // Activate loop protection at some pt value
            options.loopProtection = false;
                //(startParameters.transverseMomentum() < ptLoopers);
            // Set a maximum step size
            options.maxStepSize = maxStepSize;

            // Switch the material interaction on/off & eventually into logging mode
            auto& mInteractor = options.actionList.get<MaterialInteractor>();
            mInteractor.multipleScattering = multipleScattering;
            mInteractor.energyLoss = energyLoss;
            mInteractor.recordInteractions = recordMaterialInteractions;

            // Switch the logger to sterile, e.g. for timing checks
            auto& sLogger = options.actionList.get<SteppingLogger>();
            sLogger.sterile = sterileLogger;

            // Propagate using the propagator
            auto result = propagator.propagate(startParameters, options);
            if (result.ok()) {
                const auto& resultValue = result.value();
                auto steppingResults =
                    resultValue.template get<SteppingLogger::result_type>();

                // Set the stepping result
                pOutput.first = std::move(steppingResults.steps);
                // Also set the material recording result - if configured
                if (recordMaterialInteractions) {
                    auto materialResult =
                        resultValue.template get<MaterialInteractor::result_type>();
                    pOutput.second = std::move(materialResult);
                    mat++;
                }
            }
        }
        else {
            // Execute the test for neutral particles
            Acts::NeutralBoundTrackParameters neutralParameters(surface, pars, std::move(cov));
            sPosition = neutralParameters.position(*gctx);
            sMomentum = neutralParameters.momentum();
            std::cout << "ACTS:: Test propagation with neutral particle starts." << std::endl;

            PropagatorOptions options(*gctx, mctx);
            options.pathLimit = pathLength;
            // Activate loop protection at some pt value
            options.loopProtection =
                (neutralParameters.transverseMomentum() < ptLoopers);
            // Set a maximum step size
            options.maxStepSize = maxStepSize;

            // Switch the material interaction on/off & eventually into logging mode
            auto& mInteractor = options.actionList.get<MaterialInteractor>();
            mInteractor.multipleScattering = multipleScattering;
            mInteractor.energyLoss = energyLoss;
            mInteractor.recordInteractions = recordMaterialInteractions;

            // Switch the logger to sterile, e.g. for timing checks
            auto& sLogger = options.actionList.get<SteppingLogger>();
            sLogger.sterile = sterileLogger;

            // Propagate using the propagator
            auto result = propagator.propagate(neutralParameters, options);
            if (result.ok()) {
                const auto& resultValue = result.value();
                auto steppingResults =
                    resultValue.template get<SteppingLogger::result_type>();

                // Set the stepping result
                pOutput.first = std::move(steppingResults.steps);
                // Also set the material recording result - if configured
                if (recordMaterialInteractions) {
                    auto materialResult =
                        resultValue.template get<MaterialInteractor::result_type>();
                    pOutput.second = std::move(materialResult);
                    mat++;
                }
            }
        }
        //std::cout << mat << std::endl;
        // Record the propagation steps
        propagationSteps.push_back(std::move(pOutput.first));

        if (recordMaterialInteractions) {
            // Create a recorded material track
            RecordedMaterialTrack rmTrack;
            // Start position
            rmTrack.first.first = std::move(sPosition);
            // Start momentum
            rmTrack.first.second = std::move(sMomentum);
            // The material
            rmTrack.second = std::move(pOutput.second);
            // push it it
            recordedMaterial[it] = (std::move(rmTrack));
        }
    }

    std::cout << "ALL GOOD HERE, LENGTH OF PROPAGATION IS " << propagationSteps.size() << std::endl; 
    aPropagationStepsWriter(propagationSteps);
    aMaterialTrackWriter(recordedMaterial, *gctx);

    // Visualization	
	double outputScalor = 1.0;
  	size_t outputPrecision = 6;
  	Acts::ObjVisualization3D objVis(outputPrecision, outputScalor);

    Acts::ViewConfig containerView = Acts::ViewConfig({220, 220, 220});
  	Acts::ViewConfig volumeView = Acts::ViewConfig({220, 220, 0});
  	Acts::ViewConfig sensitiveView = Acts::ViewConfig({0, 180, 240});
  	Acts::ViewConfig passiveView = Acts::ViewConfig({240, 280, 0});
  	Acts::ViewConfig gridView = Acts::ViewConfig({220, 0, 0});

  	Acts::GeometryView3D::drawTrackingVolume(
      	objVis, *(tGeometry->highestTrackingVolume()), *gctx, containerView,
      	volumeView, passiveView, sensitiveView, gridView, true, "", "./geometry_files/");
	std::cout << "ACTS:: 3D Visualization of geometry has been dumped out." << std::endl;

    return 0;
}

// ---------------------------------------------------------------- //
//                Creator: Amaia Razquin Lizarraga                  //
//                Contact: amaiafisika@gmail.com                    //
//                Date:    7 September 2023                         //
// ---------------------------------------------------------------- //
/* This module contains helper functions for Kalman Fitting, it builds
   the propagator  and set default covariance (based on sPHENIX alg.)*/

#include <Acts/Definitions/Units.hpp>

// Track fitting headers
#include <Acts/TrackFitting/GainMatrixSmoother.hpp>
#include <Acts/TrackFitting/GainMatrixUpdater.hpp>

// Propagator headers
#include <Acts/Propagator/Propagator.hpp>
#include <Acts/Propagator/EigenStepper.hpp>
#include <Acts/Propagator/Navigator.hpp>
#include <Acts/Propagator/StraightLineStepper.hpp>
#include <Acts/Propagator/AbortList.hpp>
#include <Acts/Propagator/ActionList.hpp>
#include <Acts/Propagator/DenseEnvironmentExtension.hpp>
#include <Acts/Propagator/MaterialInteractor.hpp>
#include <Acts/Propagator/StandardAborters.hpp>
#include <Acts/Propagator/detail/SteppingLogger.hpp>

#include <Acts/EventData/VectorTrackContainer.hpp>

#include "KalmanFitting.hpp"

using namespace Acts::UnitLiterals;

// Code found in sPHENIX TrkFitter 
Acts::BoundSquareMatrix setDefaultCovariance() {
    Acts::BoundSquareMatrix cov = Acts::BoundSquareMatrix::Zero();

    double sigmaD0    = 1_mm;
    double sigmaZ0    = 1_mm;
    double sigmaPhi   = 0.1; // rad
    double sigmaTheta = 0.1; // rad
    double sigmaT     = 1_ns;

    cov(Acts::eBoundLoc0, Acts::eBoundLoc0)     = sigmaD0 * sigmaD0;
    cov(Acts::eBoundLoc1, Acts::eBoundLoc1)     = sigmaZ0 * sigmaZ0;
    cov(Acts::eBoundTime, Acts::eBoundTime)     = sigmaT * sigmaT;
    cov(Acts::eBoundPhi, Acts::eBoundPhi)       = sigmaPhi * sigmaPhi;
    cov(Acts::eBoundTheta, Acts::eBoundTheta)   = sigmaTheta * sigmaTheta;
    // Acts takes this value very seriously - tuned to be in a "sweet spot"
    // The bigger the number the smaller the confidence
    cov(Acts::eBoundQOverP, Acts::eBoundQOverP) = 0.001; //0.0001; tunned for COMET Phase-II

    return cov;
}

Acts::Propagator<Acts::EigenStepper<>, Acts::Navigator> aBuildPropagator( 
    std::shared_ptr<const Acts::TrackingGeometry> tGeometry, 
    std::shared_ptr<const Acts::InterpolatedMagneticField> bField) {
    //std::shared_ptr<const Acts::ConstantBField> bField) {
    
    // Create stepper using the magnetic field
	Acts::EigenStepper<> stepper(bField);

    // Create navigator using the tracking geometry
    Acts::Navigator::Config ncfg{tGeometry};
    ncfg.resolveMaterial             = true;
    ncfg.resolvePassive              = false;
    ncfg.resolveSensitive            = true;
	ncfg.boundaryCheckLayerResolving = true;
    Acts::Navigator navigator(ncfg, Acts::getDefaultLogger("Navigator",
      	Acts::Logging::Level::INFO));
	
    // Create propagator 
	Acts::Propagator<Acts::EigenStepper<>, Acts::Navigator> propagator(std::move(stepper), 
		std::move(navigator), Acts::getDefaultLogger("Propagator", Acts::Logging::INFO));

    return propagator;
}

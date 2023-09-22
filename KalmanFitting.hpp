#ifndef KALMANFITTING_HPP
#define KALMANFITTING_HPP

// Context
#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/MagneticField/MagneticFieldContext.hpp>

#include <Acts/MagneticField/MagneticFieldProvider.hpp>
#include <Acts/MagneticField/InterpolatedBFieldMap.hpp>
#include <Acts/MagneticField/ConstantBField.hpp>

#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/EventData/Measurement.hpp>

#include <Acts/TrackFitting/KalmanFitter.hpp>

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/EventData/TrackParameters.hpp>

#include "base/ActsSourceLink.hpp"
#include "base/Calibrator.hpp"

// The step length logger for testing & end of world aborter
using MaterialInteractor = Acts::MaterialInteractor;
using SteppingLogger = Acts::detail::SteppingLogger;
using EndOfWorld = Acts::EndOfWorldReached;

// Action list and abort list
using ActionList = Acts::ActionList<SteppingLogger, MaterialInteractor>;
using AbortList = Acts::AbortList<EndOfWorld>;

Acts::BoundSquareMatrix setDefaultCovariance();

Acts::Propagator<Acts::EigenStepper<>, Acts::Navigator> aBuildPropagator(
    std::shared_ptr<const Acts::TrackingGeometry> tGeometry, 
    std::shared_ptr<const Acts::InterpolatedMagneticField> bField);
    //std::shared_ptr<const Acts::ConstantBField> bField);

#endif

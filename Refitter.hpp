#ifndef ESTIMATOR_HPP
#define ESTIMATOR_HPP

#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Utilities/Logger.hpp>

#include "base/SpacePoint.hpp"

// Return spacepoints, final parameters
std::vector<std::pair<std::vector<double>, std::vector<Acts::Vector3>>> refitting(
    std::vector<double> chi2List,
    std::vector<Acts::BoundVariantMeasurement> measurements,
    std::vector<Acts::CurvilinearTrackParameters> initialParams);

#endif
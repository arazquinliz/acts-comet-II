#ifndef LOGGER_HPP
#define LOGGER_HPP

// Geometry headers
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/Geometry/TrackingGeometryBuilder.hpp>

void aLogGeometry(std::shared_ptr<const Acts::TrackingGeometry> tGeometry, 
    const Acts::GeometryContext& gctx);

void aSurfaceStream (std::shared_ptr<const Acts::TrackingGeometry> tGeometry, 
    const Acts::GeometryContext& gctx);

#endif
#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/Geometry/GeometryContext.hpp>

std::unique_ptr<const Acts::TrackingGeometry> 
    aBuildDiskGeometry(const Acts::GeometryContext& gctx, bool debug);
    
std::unique_ptr<const Acts::TrackingGeometry> 
    aBuildStrawGeometry(const Acts::GeometryContext& gctx, bool debug);

std::unique_ptr<const Acts::TrackingGeometry> 
    aBuildSimpleGeometry(const Acts::GeometryContext& gctx, bool debug);

#endif
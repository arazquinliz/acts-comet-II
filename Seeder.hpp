#ifndef SEEDER_HPP
#define SEEDER_HPP
 
#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/MagneticField/MagneticFieldContext.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>

#include <TVector3.h>
#include <Math/Vector3D.h>

#include "base/SpacePoint.hpp"

std::vector<SpacePoint> compareSeeds(Acts::Seed<SpacePoint> seedA, Acts::Seed<SpacePoint> seedB);

bool isIn(std::vector<Acts::Vector3> vector, Acts::Vector3 point);
bool isIn(std::vector<SpacePoint> vector, Acts::Vector3 point);
bool isIn(std::vector<Acts::Vector2> vector, Acts::Vector2 point);

std::pair<int, double> checkPurity(std::vector<SpacePoint> seed);

Acts::BoundVector estimateParams(const Acts::GeometryContext& gctx, 
    std::shared_ptr<const Acts::TrackingGeometry> tGeometry,
    Acts::MagneticFieldContext magctx, 
    std::shared_ptr<Acts::MagneticFieldProvider> BFieldmap,
    std::vector<const SpacePoint*> seedSP);

std::unique_ptr<SpacePoint> makeSpacePoint(
    const Acts::GeometryContext& gctx,
    Acts::Vector3 globalPos, 
    double time,
    std::shared_ptr<const Acts::Surface> surf,
    const int event, const Acts::Vector3 p);

void getRandomBG(std::vector<const SpacePoint*>& spVector, 
    const Acts::GeometryContext& gctx, 
    Acts::Extent& rRangeSPExtent,
    std::vector<std::shared_ptr<const Acts::Surface>> station);

void getRealBG(std::vector<const SpacePoint*>& spVector, 
    const Acts::GeometryContext& gctx, Acts::Extent& rRangeSPExtent,
    std::map<std::string, std::vector<std::shared_ptr<const Acts::Surface>>> tracker_layout,
    int stationID);

std::vector<const SpacePoint*> getSpacePoints(
    const Acts::GeometryContext& gctx,
    std::vector<ROOT::Math::XYZVector> spOfEvent,
    std::vector<ROOT::Math::XYZVector> momOfEvent,
    std::vector<double> eDeposition,
    std::vector<double> timeOfEvent, double creationT,
    std::map<std::string, std::vector<std::shared_ptr<const Acts::Surface>>> tracker_layout,
    int eventID, int stationID, Acts::Extent& rRangeSPExtent);

std::vector<const SpacePoint*> getSpacePointsTC(
    const Acts::GeometryContext& gctx,
    std::vector<ROOT::Math::XYZVector> spOfEvent,
    std::vector<ROOT::Math::XYZVector> momOfEvent,
    std::vector<double> eDeposition,
    std::vector<double> timeOfEvent, double creationT,
    std::map<std::string, std::vector<std::shared_ptr<const Acts::Surface>>> tracker_layout,
    int eventID, int stationID, Acts::Extent& rRangeSPExtent);

#endif
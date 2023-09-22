#ifndef MEASUREMENTS_HPP
#define MEASUREMENTS_HPP

#include <map>
#include <vector>

// ACTS
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/EventData/Measurement.hpp>
#include <Acts/Definitions/Algebra.hpp>

// ROOT
#include <Math/Vector3D.h>

// AMAIA
#include "base/ActsSourceLink.hpp"

struct result {std::vector<Acts::BoundVariantMeasurement> measurements; 
	std::vector<Acts::SourceLink> SourceLinks;
	Acts::Vector3 initialPos;
	Acts::Vector3 initialMom;};

std::map<std::string, std::vector<Acts::Vector3>> aReadDiskLayout(
	std::vector<ROOT::Math::XYZVector> spOfEvent,
	std::shared_ptr<const Acts::TrackingGeometry> tGeometry, 
	const Acts::GeometryContext& gctx, bool debug);

std::map<std::string, std::vector<std::pair<Acts::Vector3, Acts::Vector3>>> aReadStrawLayout(
	std::vector<ROOT::Math::XYZVector> spOfEvent,
	std::vector<ROOT::Math::XYZVector> momOfEvent,
	std::shared_ptr<const Acts::TrackingGeometry> tGeometry, 
	const Acts::GeometryContext& gctx, bool debug);

std::map<std::string, std::vector<std::tuple<Acts::Vector3, Acts::Vector3, double>>> aReadStrawLayoutNew(
	std::vector<ROOT::Math::XYZVector> spOfEvent,
	std::vector<ROOT::Math::XYZVector> momOfEvent,
	std::vector<double> driftD,
	std::shared_ptr<const Acts::TrackingGeometry> tGeometry, 
	const Acts::GeometryContext& gctx, bool debug);

std::vector<result> aCreateStrawLinks(
	std::map<std::string, std::vector<std::shared_ptr<const Acts::Surface>>> tracker_layout,
	std::shared_ptr<const Acts::TrackingGeometry> tGeometry,
	const Acts::GeometryContext& gctx, bool debug);

std::vector<result> aCreateStrawLinksNew(
	std::map<std::string, std::vector<std::shared_ptr<const Acts::Surface>>> tracker_layout,
	std::shared_ptr<const Acts::TrackingGeometry> tGeometry,
	const Acts::GeometryContext& gctx, bool debug);

std::vector<result> aCreateDiskLinks(
	std::map<std::string, std::vector<std::shared_ptr<const Acts::Surface>>> tracker_layout,
	std::shared_ptr<const Acts::TrackingGeometry> tGeometry,
	const Acts::GeometryContext& gctx, bool debug);

#endif
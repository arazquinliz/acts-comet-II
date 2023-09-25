#ifndef ESTIMATOR_HPP
#define ESTIMATOR_HPP

#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Utilities/Logger.hpp>

#include "base/SpacePoint.hpp"

struct Circle{
    double m_cx;
    double m_cy;
    double m_r;
    double m_sigma; // Root mean square error 
    double m_gradient; 
    int m_iteration; // Total outer iterations of parameter updates
    int m_inner; // Inner iterations of lambda adjustment

    double x() const { return m_cx; }
    double y() const { return m_cy; }
    double r() const { return m_r; }
    double sigma() const { return m_sigma; }
    double grad() const { return m_gradient; }
    int iter() const { return m_iteration; }
    int inner() const { return m_inner; }
};

// Estimate sigma = square root of data divided by N
double sigmaCalc(std::vector<TVector3> points, Circle Circ);

// Hyper fit by Al-Sharadqah and Chernov
Circle fitHyper(std::vector<TVector3> points);

// Taubin algebraic circle fit
Circle fitTaubin(std::vector<TVector3> points);

// Levenberg-Marquardt geometrical circle fit
double fitLevenberg(std::vector<TVector3> points);

std::optional<Acts::BoundVector> estimateModified(
    const Acts::GeometryContext& gctx, std::vector<const SpacePoint*>::iterator spBegin,
    std::vector<const SpacePoint*>::iterator spEnd, const Acts::Surface& surface, const Acts::Vector3& bField,
    Acts::ActsScalar bFieldMin, const Acts::Logger& logger,
    Acts::ActsScalar mass);

#endif
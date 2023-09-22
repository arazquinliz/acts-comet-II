// ---------------------------------------------------------------- //
//                Creator: Amaia Razquin Lizarraga                  //
//                Contact: amaiafisika@gmail.com                    //
//                Date:    7 September 2023                         //
// ---------------------------------------------------------------- //

#ifndef BASE_CalibratorNew_HPP
#define BASE_CalibratorNew_HPP

#include <Acts/EventData/Measurement.hpp>
#include <Acts/EventData/MultiTrajectory.hpp>
#include <Acts/EventData/VectorMultiTrajectory.hpp>
#include <Acts/EventData/SourceLink.hpp>

#include "ActsSourceLink.hpp"

// Calibrator found in Measurement.hpp of ActsExamples, also used in sPHENIX with modifications
// Different from Calibrator because the measurements are 1D here, used in StrawNew
class CalibratorNew
{
    using Measurement = Acts::BoundVariantMeasurement;
    // Container of measurements
    // Inn contrast to the source links, the measurements themselves must not be orderable.
    // The source links stored in the measurements are treated as opaque here and no ordering
    // is enforced on the stored measurements.
    using MeasurementContainer = std::vector<Measurement>;

  public:
    // Construct on invalid CalibratorNew. Required to arrow copying.
    CalibratorNew() = default;
    // Construct using a user-provided container to chose measurements from.
    CalibratorNew(const MeasurementContainer& measurements)
        : m_measurements(&measurements) {}

    // Find the measurement corresponding to the source link.
    // @tparam parameters_t Track parameters type
    // @param gctx The geometry context (unused)
    // @param trackState The track state to calibrate
    void calibrate(const Acts::GeometryContext& gctx,
        Acts::MultiTrajectory<Acts::VectorMultiTrajectory>::TrackStateProxy trackState)

const 
    {
        const auto& sourceLink = trackState.getUncalibratedSourceLink().get<ActsSourceLink>();
        assert(m_measurements and "Undefined measurement container in DigitizedCalibratorNew");
        assert((sourceLink.index() < m_measurements->size()) and
            "Source link index is outside the container boundes");
        std::visit(
            [&](const auto& uncalibmeas) {
                std::array<Acts::BoundIndices, 1> indices;
                indices[0] = Acts::BoundIndices::eBoundLoc0;

                Acts::ActsVector<1> loc;
                loc(0) = uncalibmeas.parameters()[Acts::eBoundLoc0];

                auto cov = uncalibmeas.covariance();
                const double misalignmentFactor = 1.0;

                Acts::ActsSymMatrix<1> expandedCov = Acts::ActsSymMatrix<1>::Identity();

                expandedCov(0, 0) = cov(0, 0)*misalignmentFactor;

                Acts::Measurement<Acts::BoundIndices, 1> meas(uncalibmeas.sourceLink(), indices,
                    loc, expandedCov);
                // back to normal
                trackState.allocateCalibrated(meas.size());
                trackState.setCalibrated(meas);
            },

            (*m_measurements)[sourceLink.index()]);
    }

  private:
    // Use pointer so the CalibratorNew is copyable and default constructible.
    const MeasurementContainer* m_measurements = nullptr;
};

#endif
// ---------------------------------------------------------------- //
//                Creator: Amaia Razquin Lizarraga                  //
//                Contact: amaiafisika@gmail.com                    //
//                Date:    7 September 2023                         //
// ---------------------------------------------------------------- //

#ifndef BASE_CALIBRATOR_HPP
#define BASE_CALIBRATOR_HPP

#include <Acts/EventData/Measurement.hpp>
#include <Acts/EventData/MultiTrajectory.hpp>
#include <Acts/EventData/VectorMultiTrajectory.hpp>
#include <Acts/EventData/SourceLink.hpp>

#include "ActsSourceLink.hpp"

// Calibrator found in Measurement.hpp of ActsExamples, also used in sPHENIX with modifications
class Calibrator
{
    using Measurement = Acts::BoundVariantMeasurement;
    // Container of measurements
    // Inn contrast to the source links, the measurements themselves must not be orderable.
    // The source links stored in the measurements are treated as opaque here and no ordering
    // is enforced on the stored measurements.
    using MeasurementContainer = std::vector<Measurement>;

  public:
    // Construct on invalid calibrator. Required to arrow copying.
    Calibrator() = default;
    // Construct using a user-provided container to chose measurements from.
    Calibrator(const MeasurementContainer& measurements)
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
        assert(m_measurements and "Undefined measurement container in DigitizedCalibrator");
        assert((sourceLink.index() < m_measurements->size()) and
            "Source link index is outside the container boundes");
        std::visit(
            [&](const auto& uncalibmeas) {
                // This part is exclusive from sPHENIX
                std::array<Acts::BoundIndices, 2> indices;
                indices[0] = Acts::BoundIndices::eBoundLoc0;
                indices[1] = Acts::BoundIndices::eBoundLoc1;

                Acts::ActsVector<2> loc;
                loc(0) = uncalibmeas.parameters()[Acts::eBoundLoc0];
                loc(1) = uncalibmeas.parameters()[Acts::eBoundLoc1];

                auto cov = uncalibmeas.covariance();
                const double misalignmentFactor = 1.0;

                Acts::ActsSquareMatrix<2> expandedCov = Acts::ActsSquareMatrix<2>::Zero();

                for (int i = 0; i < cov.rows(); i++) {
                    for (int j = 0; j < cov.cols(); j++) {
                            expandedCov(i, j) = cov(i, j)*misalignmentFactor;
                    }
                }

                Acts::Measurement<Acts::BoundIndices, 2> meas(uncalibmeas.sourceLink(), indices,
                    loc, expandedCov);
                // back to normal
                trackState.allocateCalibrated(meas.size());
                trackState.setCalibrated(meas);
            },

            (*m_measurements)[sourceLink.index()]);
    }

  private:
    // Use pointer so the calibrator is copyable and default constructible.
    const MeasurementContainer* m_measurements = nullptr;
};

#endif

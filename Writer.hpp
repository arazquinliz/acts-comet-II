#ifndef WRITER_HPP
#define WRITER_HPP

#include "Acts/Propagator/detail/SteppingLogger.hpp"

void aPropagationStepsWriter(
    std::vector<std::vector<Acts::detail::Step>>& stepCollection);

void aMaterialTrackWriter(
    const std::unordered_map<size_t, Acts::RecordedMaterialTrack>&
    materialTracks, const Acts::GeometryContext& gctx);

#endif
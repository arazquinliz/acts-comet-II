#ifndef MAGFIELD_HPP
#define MAGFIELD_HPP

// Definitions & utilities
#include <Acts/Definitions/Algebra.hpp>

// Magnetic field headers
#include <Acts/MagneticField/MagneticFieldContext.hpp>
#include <Acts/MagneticField/BFieldMapUtils.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>
#include <Acts/MagneticField/InterpolatedBFieldMap.hpp>

using InterpolatedMagneticField3 =
    Acts::InterpolatedBFieldMap<Acts::detail::Grid< 
        Acts::Vector3, Acts::detail::EquidistantAxis,
        Acts::detail::EquidistantAxis, Acts::detail::EquidistantAxis>>;     

InterpolatedMagneticField3 get_mf(
	const std::function<size_t(std::array<size_t, 3> binsXYZ,
    std::array<size_t, 3> nBinsXYZ)>& localToGlobalBin);

void writeout_mf(InterpolatedMagneticField3 BFieldmap);

std::shared_ptr<Acts::MagneticFieldProvider> aGetMagneticField(bool debug,
	Acts::MagneticFieldContext mctx,
	const std::function<size_t(std::array<size_t, 3> binsXYZ,
    std::array<size_t, 3> nBinsXYZ)>& localToGlobalBin);

#endif
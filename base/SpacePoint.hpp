// ---------------------------------------------------------------- //
//                Creator: Amaia Razquin Lizarraga                  //
//                Contact: amaiafisika@gmail.com                    //
//                Date:    7 September 2023                         //
// ---------------------------------------------------------------- //

#ifndef BASE_SPACEPOINT_HPP
#define BASE_SPACEPOINT_HPP

#include <memory>

#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <Acts/Seeding/Seed.hpp>

/** From sPHENIX SpacePoint.h and SpacePoint.hpp in UnitTests
 * A struct for Acts to take cluster information for seeding
 */
struct SpacePoint {
  double m_x;
  double m_y;
  double m_z;
  double m_r;
  double m_t;
  Acts::GeometryIdentifier m_geoId;
  double m_varianceR;
  double m_varianceZ;
  int m_eventNum;
  Acts::Vector3 m_momentum;

  double x() const { return m_x; }
  double y() const { return m_y; }
  double z() const { return m_z; }
  double r() const { return m_r; }
  double t() const { return m_t; }
  double varianceR() const { return m_varianceR; }
  double varianceZ() const { return m_varianceZ; }
  Acts::GeometryIdentifier geoID() const { return m_geoId; }
  int event() const {return m_eventNum;}
  Acts::Vector3 p() const { return m_momentum; }

};

/// This is needed by the Acts seedfinder 
inline bool operator==(SpacePoint a, SpacePoint b) {
    if (fabs(a.m_x - b.m_x) < 1e-6 && fabs(a.m_y - b.m_y) < 1e-6 &&
        fabs(a.m_z - b.m_z) < 1e-6) {
            return true;
    } else {
        return false;
    }
}

//using SpacePointPtr = std::unique_ptr<SpacePoint>;
using SeedContainer = std::vector<Acts::Seed<SpacePoint>>;

#endif
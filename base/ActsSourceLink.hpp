// ---------------------------------------------------------------- //
//                Creator: Amaia Razquin Lizarraga                  //
//                Contact: amaiafisika@gmail.com                    //
//                Date:    7 September 2023                         //
// ---------------------------------------------------------------- //

#ifndef BASE_ACTSSOURCELINK_HPP
#define BASE_ACTSSOURCELINK_HPP

#include <Acts/EventData/SourceLink.hpp>
#include <Acts/Surfaces/Surface.hpp>

#include <cassert>
#include <iostream>

// Class similar to ActsSourceLink in sPHENIX and almost identical to IndexSourceLink except 
// container and lookup methods are not added YET (potentially TODO in the future)

/// A source link that stores just an index.
/// Using an index instead of e.g. a pointer, means source link and
/// measurement are decoupled and the measurement represenation can be
/// easily changed without having to also change the source link.
class ActsSourceLink final // : public Acts::SourceLink cannot derive from SL because it is already final
{
  public:
    using Index = uint8_t;

    // Construct from geometry identifier and index.
    constexpr ActsSourceLink(Acts::GeometryIdentifier gid, Index idx)
        : m_geometryId(gid)
        , m_index(idx) {}

    // Construct an invalid source link. Must be default constructible to 
    // satisfy SourceLinkConcept.
    ActsSourceLink()
        : m_geometryId{Acts::GeometryIdentifier{}}
        , m_index(UINT8_MAX) {}

    ActsSourceLink(const ActsSourceLink&)            = default;
    ActsSourceLink(ActsSourceLink&&)                 = default;
    ActsSourceLink& operator=(const ActsSourceLink&) = default;
    ActsSourceLink& operator=(ActsSourceLink&&)      = default;

    // Access the index
    constexpr Index index() const { return m_index; }

    constexpr Acts::GeometryIdentifier geometryId() const { return m_geometryId; }

  private:
    Index m_index;
    Acts::GeometryIdentifier m_geometryId;

    friend constexpr bool operator==(const ActsSourceLink& lhs, 
                                     const ActsSourceLink& rhs) {
        return (lhs.geometryId() == rhs.geometryId()) and
               (lhs.m_index == rhs.m_index);
    }
    friend constexpr bool operator!=(const ActsSourceLink& lhs, 
                                     const ActsSourceLink& rhs){
        return not(lhs == rhs);
    }
};

/// The map(-like) container accessor
template <typename container_t>
struct ActsContainerAccessor {
  using Container = container_t;
  using Key = typename container_t::key_type;
  using Value = typename container_t::mapped_type;

  /// This iterator adapter is needed to have the deref operator return a single
  /// source link instead of the map pair <GeometryIdentifier,SourceLink>
  struct Iterator {
    using BaseIterator = typename container_t::const_iterator;

    using iterator_category = typename BaseIterator::iterator_category;
    using value_type = typename BaseIterator::value_type;
    using difference_type = typename BaseIterator::difference_type;
    using pointer = typename BaseIterator::pointer;
    using reference = typename BaseIterator::reference;

    Iterator& operator++() {
      ++m_iterator;
      return *this;
    }

    bool operator==(const Iterator& other) const {
      return m_iterator == other.m_iterator;
    }

    bool operator!=(const Iterator& other) const { return !(*this == other); }

    Acts::SourceLink operator*() const {
      const auto& sl = m_iterator->second;
      return Acts::SourceLink{sl};
    }

    BaseIterator m_iterator;
  };

  // pointer to the container
  const Container* container = nullptr;

  // get the range of elements with requested key
  std::pair<Iterator, Iterator> range(const Acts::Surface& surface) const {
    assert(container != nullptr);
    auto [begin, end] = container->equal_range(surface.geometryId());
    return {Iterator{begin}, Iterator{end}};
  }
};

using SourceLinkContainer = std::unordered_multimap<Acts::GeometryIdentifier, Acts::SourceLink>;

using SourceLinkAccessor = ActsContainerAccessor<SourceLinkContainer>;

#endif
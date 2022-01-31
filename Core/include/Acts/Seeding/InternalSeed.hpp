// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/InternalSpacePoint.hpp"
#include "Acts/Seeding/Seed.hpp"

#include <memory>

namespace Acts {
template <typename SpacePoint>
class InternalSeed {
  /////////////////////////////////////////////////////////////////////////////////
  // Public methods:
  /////////////////////////////////////////////////////////////////////////////////

 public:
  InternalSeed(const InternalSpacePoint<SpacePoint>& s0,
               const InternalSpacePoint<SpacePoint>& s1,
               const InternalSpacePoint<SpacePoint>& s2, float z);
  InternalSeed& operator=(const InternalSeed& seed);

  const std::array<const InternalSpacePoint<SpacePoint>*, 3> sp;
  float z() const { return m_z; }

 protected:
  float m_z;
};

/// @cond

/////////////////////////////////////////////////////////////////////////////////
// Inline methods
/////////////////////////////////////////////////////////////////////////////////

template <typename SpacePoint>
inline InternalSeed<SpacePoint>& InternalSeed<SpacePoint>::operator=(
    const InternalSeed<SpacePoint>& seed) {
  m_z = seed.m_z;
  sp = seed.sp;
  return (*this);
}

template <typename SpacePoint>
inline InternalSeed<SpacePoint>::InternalSeed(
    const InternalSpacePoint<SpacePoint>& s0,
    const InternalSpacePoint<SpacePoint>& s1,
    const InternalSpacePoint<SpacePoint>& s2, float z)
    : sp({&s0, &s1, &s2}) {
  m_z = z;
}

// template <typename SpacePoint>
// inline bool InternalSeed<SpacePoint>::setQuality(float q) {
//	if(q > sp[0]->quality() and q > sp[1]->quality() and q >
// sp[2]->quality()) return false; 	sp[0]->setQuality(q);
// sp[1]->setQuality(q); 	sp[2]->setQuality(q); 	return true;
//}

/// @endcond

}  // namespace Acts

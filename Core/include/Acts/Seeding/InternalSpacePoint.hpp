// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
//#include "ActsExamples/EventData/SimSpacePoint.hpp"

#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

namespace Acts {

struct DoubleMeasurementDetails {
  float topHalfStripLength;
  float bottomHalfStripLength;
  Acts::Vector3 topStripDirection;
  Acts::Vector3 bottomStripDirection;
  Acts::Vector3 stripCenterDistance;
  Acts::Vector3 bottomStripCenterPosition;
  bool validDoubleMeasurementDetails = false;
};

template <typename SpacePoint>
class InternalSpacePoint {
  /////////////////////////////////////////////////////////////////////////////////
  // Public methods:
  /////////////////////////////////////////////////////////////////////////////////

 public:
  InternalSpacePoint() = delete;
  InternalSpacePoint(const SpacePoint& sp, const Acts::Vector3& globalPos,
                     const Acts::Vector2& offsetXY,
                     const Acts::Vector2& variance);

  InternalSpacePoint(const InternalSpacePoint<SpacePoint>& sp);
  ~InternalSpacePoint() = default;

  InternalSpacePoint<SpacePoint>& operator=(
      const InternalSpacePoint<SpacePoint>&);

  const float& x() const { return m_x; }
  const float& y() const { return m_y; }
  const float& z() const { return m_z; }
  const float& radius() const { return m_r; }
  float phi() const { return atan2f(m_y, m_x); }
  // float phi() const { float value = atan2f(m_y, m_x); return value>0 ? value
  // : value+2*M_PI;}
  const float& varianceR() const { return m_varianceR; }
  const float& varianceZ() const { return m_varianceZ; }
  float& cotTheta() const { return m_cotTheta; }
  void setCotTheta(float& cotTheta) const { m_cotTheta = cotTheta; }
  float& deltaR() const { return m_deltaR; }
  void setDeltaR(const float& deltaR) const { m_deltaR = deltaR; }

  float topHalfStripLength() const {
    return m_doubleMeasurement.topHalfStripLength;
  };
  float bottomHalfStripLength() const {
    return m_doubleMeasurement.bottomHalfStripLength;
  };
  const Acts::Vector3& topStripDirection() const {
    return m_doubleMeasurement.topStripDirection;
  }
  const Acts::Vector3& bottomStripDirection() const {
    return m_doubleMeasurement.bottomStripDirection;
  }
  const Acts::Vector3& stripCenterDistance() const {
    return m_doubleMeasurement.stripCenterDistance;
  }
  const Acts::Vector3& bottomStripCenterPosition() const {
    return m_doubleMeasurement.bottomStripCenterPosition;
  }
  bool validDoubleMeasurementDetails() const {
    return m_doubleMeasurement.validDoubleMeasurementDetails;
  }

  const float& curvature() const { return m_cotTheta; }
  void curvature(float& cotTheta) const { m_cotTheta = cotTheta; }

  const SpacePoint& sp() const { return m_sp; }
  float& quality() const { return m_quality; }
  void setQuality(float& quality) const {
    if (quality >= m_quality) {
      std::cout << "|setQuality| Old quality :" << m_quality
                << " New quality: " << quality << std::endl;
      m_quality = quality;
    }
  }

 protected:
  float m_x;                 // x-coordinate in beam system coordinates
  float m_y;                 // y-coordinate in beam system coordinates
  float m_z;                 // z-coordinate in beam system coordinetes
  float m_r;                 // radius       in beam system coordinates
  float m_varianceR;         //
  float m_varianceZ;         //
  mutable float m_cotTheta;  //
  mutable float m_deltaR;    //
  const SpacePoint& m_sp;    // external space point
  mutable float m_quality;   // quality of the best seed containing this SP

  DoubleMeasurementDetails m_doubleMeasurement;
};

/////////////////////////////////////////////////////////////////////////////////
// Inline methods
/////////////////////////////////////////////////////////////////////////////////

template <typename SpacePoint>
inline InternalSpacePoint<SpacePoint>::InternalSpacePoint(
    const SpacePoint& sp, const Acts::Vector3& globalPos,
    const Acts::Vector2& offsetXY, const Acts::Vector2& variance)
    : m_sp(sp) {
  m_x = globalPos.x() - offsetXY.x();
  m_y = globalPos.y() - offsetXY.y();
  m_z = globalPos.z();
  m_r = std::sqrt(m_x * m_x + m_y * m_y);
  m_varianceR = variance.x();
  m_varianceZ = variance.y();
  m_quality = -100000.;

  m_doubleMeasurement.topHalfStripLength = sp.topHalfStripLength();
  m_doubleMeasurement.bottomHalfStripLength = sp.bottomHalfStripLength();
  m_doubleMeasurement.topStripDirection = sp.topStripDirection();
  m_doubleMeasurement.bottomStripDirection = sp.bottomStripDirection();
  m_doubleMeasurement.stripCenterDistance = sp.stripCenterDistance();
  m_doubleMeasurement.bottomStripCenterPosition =
      sp.bottomStripCenterPosition();
  m_doubleMeasurement.validDoubleMeasurementDetails =
      sp.validDoubleMeasurementDetails();
}

/////////////////////////////////////////////////////////////////////////////////
// Copy constructor
/////////////////////////////////////////////////////////////////////////////////

template <typename SpacePoint>
inline InternalSpacePoint<SpacePoint>::InternalSpacePoint(
    const InternalSpacePoint<SpacePoint>& sp)
    : m_sp(sp.sp()) {
  m_x = sp.m_x;
  m_y = sp.m_y;
  m_z = sp.m_z;
  m_r = sp.m_r;
  m_varianceR = sp.m_varianceR;
  m_varianceZ = sp.m_varianceZ;
  m_cotTheta = sp.m_cotTheta;
  m_deltaR = sp.deltaR;
  m_quality = sp.m_quality;

  m_doubleMeasurement.topHalfStripLength = sp.topHalfStripLength();
  m_doubleMeasurement.bottomHalfStripLength = sp.bottomHalfStripLength();
  m_doubleMeasurement.topStripDirection = sp.topStripDirection();
  m_doubleMeasurement.bottomStripDirection = sp.bottomStripDirection();
  m_doubleMeasurement.stripCenterDistance = sp.stripCenterDistance();
  m_doubleMeasurement.bottomStripCenterPosition =
      sp.bottomStripCenterPosition();
  m_doubleMeasurement.validDoubleMeasurementDetails =
      sp.validDoubleMeasurementDetails();
}

}  // end of namespace Acts

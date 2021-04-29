// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/detail/Axis.hpp"

#include <memory>

#include <iostream>

using AxisScalar = Acts::Vector3::Scalar;

template <typename SpacePoint>
std::unique_ptr<Acts::SpacePointGrid<SpacePoint>>
Acts::SpacePointGridCreator::createGrid(
    const Acts::SpacePointGridConfig& config) {
  // calculate circle intersections of helix and max detector radius
  float minHelixRadius = config.minPt / (300. * config.bFieldInZ);  // in mm  
  float maxR2 = config.rMax * config.rMax;
  float xOuter = maxR2 / (2 * minHelixRadius);
  float yOuter = std::sqrt(maxR2 - xOuter * xOuter);
  float outerAngle = std::atan(xOuter / yOuter);
  // intersection of helix and max detector radius minus maximum R distance from
  // middle SP to top SP
  float innerAngle = 0;
  float Rmin = config.rMax;
  if (config.rMax > config.deltaRMax) {
    Rmin = config.rMax - config.deltaRMax;
    float innerCircleR2 =
        (config.rMax - config.deltaRMax) * (config.rMax - config.deltaRMax);
    float xInner = innerCircleR2 / (2 * minHelixRadius);
    float yInner = std::sqrt(innerCircleR2 - xInner * xInner);
    innerAngle = std::atan(xInner / yInner);
  }
  
  std::cout << "==============================" << std::endl;
  std::cout << "Dumping configuration: " << std::endl;
  std::cout << "config.rMax = " << config.rMax << std::endl;
  std::cout << "config.rMin = " << (config.rMax-config.deltaRMax)  << std::endl;
  std::cout << "config.impactMax = " << config.impactMax << std::endl;
  std::cout << "config.numberOfPhiBins = " << config.numberOfPhiBins << std::endl;
  std::cout << "==============================" << std::endl;
  std::cout << "minHelixRadius = " << minHelixRadius << std::endl;
  std::cout << "sI = " << outerAngle - innerAngle << std::endl;
  
  // evaluating the azimutal deflection including the maximum impact parameter
  float deltaAngleWithMaxD0 = std::abs(std::asin(config.impactMax/(Rmin)) - std::asin(config.impactMax/config.rMax));
  
  std::cout << "sF = " << deltaAngleWithMaxD0 << std::endl;
  
  float deltaPhi = (outerAngle - innerAngle + deltaAngleWithMaxD0)/float(config.numberOfPhiBins);
  
  std::cout << "deltaPhi = " << deltaPhi << std::endl;
  
  // divide 2pi by angle delta to get number of phi-bins
  // size is always 2pi even for regions of interest
  int phiBins = std::ceil(2 * M_PI / deltaPhi);
  // need to scale the number of phi bins accordingly to the number of 
  // consecutive phi bins in the seed making step. 
  // Each individual bin should be approximately a fraction (depending on this number) 
  // of the maximum expected azimutal deflection.
  Acts::detail::Axis<detail::AxisType::Equidistant,
                     detail::AxisBoundaryType::Closed>
      phiAxis(-M_PI, M_PI, phiBins);
      
  std::cout << "Number of bins in phi = " << phiBins << std::endl;
  
  std::vector<AxisScalar> zValues;
  
  if (config.zBinEdges.empty()) {
    // TODO: can probably be optimized using smaller z bins
    // and returning (multiple) neighbors only in one z-direction for forward
    // seeds
    // FIXME: zBinSize must include scattering
    float zBinSize = config.cotThetaMax * config.deltaRMax;
    int zBins = std::floor((config.zMax - config.zMin) / zBinSize);
    for (int bin=0; bin<=zBins; bin++){
      AxisScalar edge = config.zMin + bin*zBinSize;
      zValues.push_back(edge);
    }
  } else {
    for (auto& bin : config.zBinEdges) {
      std::cout << "z bins --> " << bin << ", ";
      AxisScalar edge = bin;
      zValues.push_back(edge);
    }
    std::cout << std::endl;
  }
  std::cout << "==============================" << std::endl;
  
  detail::Axis<detail::AxisType::Variable, detail::AxisBoundaryType::Bound>
      zAxis(zValues);
  return std::make_unique<Acts::SpacePointGrid<SpacePoint>>(
      std::make_tuple(phiAxis, zAxis));  
}

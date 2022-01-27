// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/SeedFilter.hpp"

#include <cmath>
#include <iostream>
#include <numeric>
#include <type_traits>

namespace Acts {

template <typename external_spacepoint_t, typename platform_t>
Seedfinder<external_spacepoint_t, platform_t>::Seedfinder(
    Acts::SeedfinderConfig<external_spacepoint_t> config)
    : m_config(config.toInternalUnits()) {
  // calculation of scattering using the highland formula
  // convert pT to p once theta angle is known
  m_config.highland = 13.6 * std::sqrt(m_config.radLengthPerSeed) *
                      (1 + 0.038 * std::log(m_config.radLengthPerSeed));
  float maxScatteringAngle = m_config.highland / m_config.minPt;
  m_config.maxScatteringAngle2 = maxScatteringAngle * maxScatteringAngle;
  // helix radius in homogeneous magnetic field. Units are Kilotesla, MeV and
  // millimeter
  // TODO: change using ACTS units
  m_config.pTPerHelixRadius = 300. * m_config.bFieldInZ;
  m_config.minHelixDiameter2 =
      std::pow(m_config.minPt * 2 / m_config.pTPerHelixRadius, 2);
  m_config.pT2perRadius =
      std::pow(m_config.highland / m_config.pTPerHelixRadius, 2);
<<<<<<< Updated upstream
  m_config.sigmapT2perRadius =
      m_config.pT2perRadius * std::pow(2 * m_config.sigmaScattering, 2);
=======
  m_config.sigmapT2perRadius =
      m_config.pT2perRadius * std::pow(2 * m_config.sigmaScattering, 2);
>>>>>>> Stashed changes
}

template <typename external_spacepoint_t, typename platform_t>
template <template <typename...> typename container_t, typename sp_range_t>
void Seedfinder<external_spacepoint_t, platform_t>::createSeedsForGroup(
    State& state,
    std::back_insert_iterator<container_t<Seed<external_spacepoint_t>>> outIt,
<<<<<<< Updated upstream
    sp_range_t bottomSPs, sp_range_t middleSPs, sp_range_t topSPs,
    Extent rRangeSPExtent) const {
=======
    sp_range_t bottomSPs, sp_range_t middleSPs, sp_range_t topSPs,
    Extent rRangeSPExtent) const {
  //	std::sort(middleSPs.begin(), middleSPs.end(), [] (auto* a, auto* b) ->
  //bool { return (a->radius() < b->radius()); } );

  int i_n = 0;
>>>>>>> Stashed changes
  for (auto spM : middleSPs) {
    // std::cout << std::endl;
    // std::cout << " rRangeSPExtent min | " << rRangeSPExtent.min(Acts::binR)
    // << std::endl; std::cout << " rRangeSPExtent max | " <<
    // rRangeSPExtent.max(Acts::binR)  << std::endl;

    float rM = spM->radius();
    float zM = spM->z();
    float varianceRM = spM->varianceR();
    float varianceZM = spM->varianceZ();
    float rMin = m_config.rMax - m_config.deltaRMax;

    // std::cout << std::endl;
    // std::cout << "|Seeds| rM: " << rM << " deltaRMin: " << m_config.deltaRMin
    // << std::endl;

<<<<<<< Updated upstream
    /// check if spM is outside our radial region of interest
    if (m_config.useVariableMiddleSPRange) {
      float rMinMiddleSP = std::floor(rRangeSPExtent.min(Acts::binR) / 2) * 2 +
                           m_config.deltaRMiddleSPRange;
      float rMaxMiddleSP = std::floor(rRangeSPExtent.max(Acts::binR) / 2) * 2 -
                           m_config.deltaRMiddleSPRange;
      if (rM < rMinMiddleSP || rM > rMaxMiddleSP) {
        continue;
      }
    } else if (not m_config.rRangeMiddleSP.empty()) {
      /// get zBin position of the middle SP
      auto pVal = std::lower_bound(m_config.zBinEdges.begin(),
                                   m_config.zBinEdges.end(), zM);
      int zBin = std::distance(m_config.zBinEdges.begin(), pVal);
      /// protects against zM at the limit of zBinEdges
      zBin == 0 ? zBin : --zBin;
      if (rM < m_config.rRangeMiddleSP[zBin][0] ||
          rM > m_config.rRangeMiddleSP[zBin][1]) {
        continue;
      }
    }

    size_t nTopSeedConf = 0;
    if (m_config.seedConfirmation == true) {
      // check if middle SP is in the central or forward region
      SeedConfirmationRange seedConfRange =
          (zM > m_config.centralSeedConfirmationRange.zMaxSeedConf ||
           zM < m_config.centralSeedConfirmationRange.zMinSeedConf)
              ? m_config.forwardSeedConfirmationRange
              : m_config.centralSeedConfirmationRange;
      // set the minimum number of top SP depending on whether the middle SP is
      // in the central or forward region
      nTopSeedConf = rM > seedConfRange.rMaxSeedConf
                         ? seedConfRange.nTopForLargeR
                         : seedConfRange.nTopForSmallR;
    }

    state.compatTopSP.clear();

    for (auto topSP : topSPs) {
      float rT = topSP->radius();
      float deltaR = rT - rM;
      // if r-distance is too small, try next SP in bin
      if (deltaR < m_config.deltaRMinTopSP) {
        continue;
      }
      // if r-distance is too big, try next SP in bin
      if (deltaR > m_config.deltaRMaxTopSP) {
        continue;
      }
      // ratio Z/R (forward angle) of space point duplet
      float cotTheta = (topSP->z() - zM) / deltaR;
      if (std::fabs(cotTheta) > m_config.cotThetaMax) {
        continue;
      }
      // check if duplet origin on z axis within collision region
      float zOrigin = zM - rM * cotTheta;
      if (zOrigin < m_config.collisionRegionMin ||
          zOrigin > m_config.collisionRegionMax) {
        continue;
      }
      // cut on the max curvature between top SP and interaction point
      // first transform the space point coordinates into a frame such that the
      // central space point SPm is in the origin of the frame and the x axis
      // points away from the interaction point in addition to a translation
      // transformation we also perform a rotation in order to keep the
      // curvature of the circle tangent to the x axis
      float xVal = (topSP->x() - spM->x()) * (spM->x() / rM) +
                   (topSP->y() - spM->y()) * (spM->y() / rM);
      float yVal = (topSP->y() - spM->y()) * (spM->x() / rM) -
                   (topSP->x() - spM->x()) * (spM->y() / rM);
      if (std::abs(rM * yVal) > m_config.impactMax * xVal) {
        // conformal transformation u=x/(x²+y²) v=y/(x²+y²) transform the circle
        // into straight lines in the u/v plane the line equation can be
        // described in terms of aCoef and bCoef, where v = aCoef * u + bCoef
        float uT = xVal / (xVal * xVal + yVal * yVal);
        float vT = yVal / (xVal * xVal + yVal * yVal);
        // in the rotated frame the interaction point is positioned at x = -rM
        // and y ~= impactParam
        float uIP = -1. / rM;
        float vIP = m_config.impactMax / (rM * rM);
        if (yVal > 0.)
          vIP = -vIP;
        // we can obtain aCoef as the slope dv/du of the linear function,
        // estimated using du and dv between the two SP bCoef is obtained by
        // inserting aCoef into the linear equation
        float aCoef = (vT - vIP) / (uT - uIP);
        float bCoef = vIP - aCoef * uIP;
        // the distance of the straight line from the origin (radius of the
        // circle) is related to aCoef and bCoef by d^2 = bCoef^2 / (1 +
        // aCoef^2) = 1 / (radius^2) and we can apply the cut on the curvature
        if ((bCoef * bCoef) >
            (1 + aCoef * aCoef) / m_config.minHelixDiameter2) {
          continue;
        }
      }
      state.compatTopSP.push_back(topSP);
    }
    if (state.compatTopSP.empty()) {
      continue;
    }
    // apply cut on the number of top SP if seedConfirmation is true
    if (m_config.seedConfirmation == true &&
        state.compatTopSP.size() < nTopSeedConf) {
      continue;
    }

    state.compatBottomSP.clear();

    for (auto bottomSP : bottomSPs) {
      float rB = bottomSP->radius();
      float deltaR = rM - rB;
      // this condition is the opposite of the condition for top SP
      if (deltaR > m_config.deltaRMaxBottomSP) {
        continue;
      }
      if (deltaR < m_config.deltaRMinBottomSP) {
        continue;
      }
      // ratio Z/R (forward angle) of space point duplet
      float cotTheta = (zM - bottomSP->z()) / deltaR;
=======
    //		//std::cout << "|Seeds| rMinMiddleSP: " << m_config.rRanMiddleSP[0] <<
    //" rMaxMiddleSP: " << m_config.rRanMiddleSP[1] << std::endl;

    // delete this --> just for printing values:
    size_t nSeeds = 0;
    size_t nSeeds_filter1 = 0;
    size_t nSeeds_filter2 = 0;
    size_t nSeeds_test1 = 0;
    size_t nSeeds_test2 = 0;
    size_t nSeeds_test3 = 0;
    auto pVal = std::lower_bound(m_config.zBinEdges.begin(),
                                 m_config.zBinEdges.end(), zM);
    int zBin = std::distance(m_config.zBinEdges.begin(), pVal);
    zBin == 0 ? zBin : --zBin;

    /// check if spM is outside our radial region of interest
    if (m_config.useVariableMiddleSPRange) {
      float rMinMiddleSP = std::floor(rRangeSPExtent.min(Acts::binR) / 2) * 2 +
                           m_config.deltaRMiddleSPRange;
      float rMaxMiddleSP = std::floor(rRangeSPExtent.max(Acts::binR) / 2) * 2 -
                           m_config.deltaRMiddleSPRange;
      if (rM < rMinMiddleSP || rM > rMaxMiddleSP) {
        continue;
      }
    } else if (not m_config.rRangeMiddleSP.empty()) {
      /// get zBin position of the middle SP
      auto pVal = std::lower_bound(m_config.zBinEdges.begin(),
                                   m_config.zBinEdges.end(), zM);
      int zBin = std::distance(m_config.zBinEdges.begin(), pVal);
      /// protects against zM at the limit of zBinEdges
      zBin == 0 ? zBin : --zBin;
      if (rM < m_config.rRangeMiddleSP[zBin][0] ||
          rM > m_config.rRangeMiddleSP[zBin][1]) {
        // std::cout << "|Seeds| !!! rM < rRangeMiddleSP || rM > rRangeMiddleSP
        // == TRUE !!!" << std::endl;
        continue;
      }
    }

    /// check if spM is on the last disk
    //    if (zM > m_config.zMax || zM < 0) {
    //      //std::cout << "|Seeds| !!! zM > m_config.zMax || zM < m_config.zMin
    //      == TRUE !!!" << std::endl; continue;
    //    }

    // check if middle SP is in the central or forward region
    SeedConfirmationRange seedConfRange;
    if (m_config.seedConfirmation == true) {
      seedConfRange =
          (zM > m_config.centralSeedConfirmationRange.zMaxSeedConf ||
           zM < m_config.centralSeedConfirmationRange.zMinSeedConf)
              ? m_config.forwardSeedConfirmationRange
              : m_config.centralSeedConfirmationRange;
      seedConfRange.nTopSeedConf = rM > seedConfRange.rMaxSeedConf
                                       ? seedConfRange.nTopForLargeR
                                       : seedConfRange.nTopForSmallR;
    }

    state.compatTopSP.clear();

    // std::cout << i_n << ") TOP SP" << std::endl;

    for (auto topSP : topSPs) {
      float rT = topSP->radius();
      float deltaR = rT - rM;
      // std::cout << std::endl;
      // std::cout << "|Seeds| --> rT: " << rT << " deltaR: " << deltaR <<
      // std::endl; std::cout << "|Seeds| deltaRMax: " << m_config.deltaRMaxTopSP
      // << " deltaRMin: " << m_config.deltaRMinTopSP << std::endl;

      // **** CHANGE COMENTS ****
      // this condition is the opposite of the condition for bottom SP
      if (deltaR < m_config.deltaRMinTopSP) {
        // std::cout << "|Seeds| !!! deltaR < deltaRMin == TRUE !!!" <<
        // std::endl;
        continue;
      }
      if (deltaR > m_config.deltaRMaxTopSP) {
        // std::cout << "|Seeds| !!! deltaR > deltaRMax == TRUE !!!" <<
        // std::endl;
        continue;
      }
      float cotTheta = (topSP->z() - zM) / deltaR;
>>>>>>> Stashed changes
      if (std::fabs(cotTheta) > m_config.cotThetaMax) {
        // std::cout << "|Seeds| !!! fabs(cotTheta) > cotThetaMax == TRUE !!!"
        // << std::endl;
        continue;
      }
      // check if duplet origin on z axis within collision region
      float zOrigin = zM - rM * cotTheta;
      // std::cout << "|Seeds| zM - rM * cotTheta " << zM << " - " << rM << " *
      // " << cotTheta << std::endl; std::cout << "|Seeds| zOrigin: " << zOrigin
      // << " collisionRegionMin: " << m_config.collisionRegionMin << std::endl;
      if (zOrigin < m_config.collisionRegionMin ||
          zOrigin > m_config.collisionRegionMax) {
        // std::cout << "|Seeds| !!! zOrigin < collisionRegionMin || zOrigin >
        // collisionRegionMax == TRUE !!!" << std::endl;
        continue;
      }
<<<<<<< Updated upstream
      // cut on the max curvature between bottom SP and interaction point
      // first transform the space point coordinates into a frame such that the
      // central space point SPm is in the origin of the frame and the x axis
      // points away from the interaction point in addition to a translation
      // transformation we also perform a rotation in order to keep the
      // curvature of the circle tangent to the x axis
      float xVal = (bottomSP->x() - spM->x()) * (spM->x() / rM) +
                   (bottomSP->y() - spM->y()) * (spM->y() / rM);
      float yVal = (bottomSP->y() - spM->y()) * (spM->x() / rM) -
                   (bottomSP->x() - spM->x()) * (spM->y() / rM);
      if (std::abs(rM * yVal) > -m_config.impactMax * xVal) {
        // conformal transformation u=x/(x²+y²) v=y/(x²+y²) transform the circle
        // into straight lines in the u/v plane the line equation can be
        // described in terms of aCoef and bCoef, where v = aCoef * u + bCoef
        float uB = xVal / (xVal * xVal + yVal * yVal);
        float vB = yVal / (xVal * xVal + yVal * yVal);
        // in the rotated frame the interaction point is positioned at x = -rM
        // and y ~= impactParam
        float uIP = -1. / rM;
        float vIP = m_config.impactMax / (rM * rM);
        if (yVal < 0.)
          vIP = -vIP;
        // we can obtain aCoef as the slope dv/du of the linear function,
        // estimated using du and dv between the two SP bCoef is obtained by
        // inserting aCoef into the linear equation
        float aCoef = (vB - vIP) / (uB - uIP);
        float bCoef = vIP - aCoef * uIP;
        // the distance of the straight line from the origin (radius of the
        // circle) is related to aCoef and bCoef by d^2 = bCoef^2 / (1 +
        // aCoef^2) = 1 / (radius^2) and we can apply the cut on the curvature
        if ((bCoef * bCoef) >
            (1 + aCoef * aCoef) / m_config.minHelixDiameter2) {
          continue;
        }
      }
      state.compatBottomSP.push_back(bottomSP);
    }
    // no bottom SP found -> try next spM
    if (state.compatBottomSP.empty()) {
=======
      // cut on the max curvature between top SP and interaction point
      // first transform the space point coordinates into a frame such that the
      // central space point SPm is in the origin of the frame and the x axis
      // points away from the interaction point in addition to a translation
      // transformation we also perform a rotation in order to keep the
      // curvature of the circle tangent to the x axis
      float xVal = (topSP->x() - spM->x()) * (spM->x() / rM) +
                   (topSP->y() - spM->y()) * (spM->y() / rM);
      float yVal = (topSP->y() - spM->y()) * (spM->x() / rM) -
                   (topSP->x() - spM->x()) * (spM->y() / rM);
      // std::cout << std::abs(rM) << " * " << std::abs(yVal) << " > " <<
      // -m_config.impactMax << " * " << xVal << std::endl;
      if (std::abs(rM * yVal) > m_config.impactMax * xVal) {
        // conformal transformation u=x/(x²+y²) v=y/(x²+y²) transform the circle
        // into straight lines in the u/v plane the line equation can be
        // described in terms of aCoef and bCoef, where v = aCoef * u + bCoef
        float uT = xVal / (xVal * xVal + yVal * yVal);
        float vT = yVal / (xVal * xVal + yVal * yVal);
        // in the rotated frame the interaction point is positioned at x = -rM
        // and y ~= impactParam
        float uIP = -1. / rM;
        float vIP = m_config.impactMax / (rM * rM);
        if (yVal > 0.)
          vIP = -vIP;
        // we can obtain aCoef as the slope dv/du of the linear function,
        // estimated using du and dv between the two SP bCoef is obtained by
        // inserting aCoef into the linear equation
        float aCoef = (vT - vIP) / (uT - uIP);
        float bCoef = vIP - aCoef * uIP;
        // std::cout << bCoef << " * " << bCoef << " > 1+" << aCoef << " * " <<
        // aCoef << " / " << m_config.minHelixDiameter2 << std::endl;
        // the distance of the straight line from the origin is given by d^2 =
        // bCoef^2 / (1 + aCoef^2) and we can apply the cut on the curvature
        if ((bCoef * bCoef) >
            (1 + aCoef * aCoef) / m_config.minHelixDiameter2) {
          // std::cout << "|Seeds| !!!  impact parameter cut == TRUE !!!" <<
          // std::endl;
          continue;
        }
      }
      //			if ( std::abs((topSP->z() - zM) / std::sqrt(xVal
      //* xVal + yVal * yVal)) > m_config.cotThetaMax) {
      //				//std::cout << "|Seeds| !!!
      //fabs(cotTheta) > cotThetaMax == TRUE !!!" << std::endl; 				continue;
      //			}
      state.compatTopSP.push_back(topSP);
      std::cout << "|Seeds| # Fill SP" << std::endl;
    }
    if (state.compatTopSP.empty()) {
      // std::cout << "|Seeds| !!! top SP empty == TRUE !!!" << std::endl;
>>>>>>> Stashed changes
      continue;
    }
    if (m_config.seedConfirmation == true &&
        state.compatTopSP.size() < seedConfRange.nTopSeedConf) {
      // std::cout << "|Seeds| !!! seedConfirmation top SP cut == TRUE !!!" <<
      // std::endl;
      continue;
    }

    state.compatBottomSP.clear();

    // std::cout << i_n << ") BOTTOM SP" << std::endl;
    i_n++;

    for (auto bottomSP : bottomSPs) {
      float rB = bottomSP->radius();
      float deltaR = rM - rB;
      // std::cout << std::endl;
      // std::cout << "|Seeds| --> rB: " << rB << " deltaR: " << deltaR <<
      // std::endl;
      //      //std::cout << "|Seeds| deltaRMax: " << m_config.deltaRMaxBottomSP
      //      << " deltaRMin: " << m_config.deltaRMinBottomSP << std::endl;
      // if r-distance is too big, try next SP in bin
      if (deltaR > m_config.deltaRMaxBottomSP) {
        // std::cout << "|Seeds| !!! deltaR > deltaRMax == TRUE !!!" <<
        // std::endl;
        continue;
      }
      // if r-distance is too small, continue because bins are NOT r-sorted
      if (deltaR < m_config.deltaRMinBottomSP) {
        // std::cout << "|Seeds| !!! deltaR < deltaRMin == TRUE !!!" <<
        // std::endl;
        continue;
      }
      // ratio Z/R (forward angle) of space point duplet
      float cotTheta = (zM - bottomSP->z()) / deltaR;
      if (std::fabs(cotTheta) > m_config.cotThetaMax) {
        //				//std::cout << "|Seeds| cotTheta: " <<
        //cotTheta << " cotThetaMax: " << m_config.cotThetaMax << std::endl;
        // std::cout << "|Seeds| !!! fabs(cotTheta) > cotThetaMax == TRUE !!!"
        // << std::endl;
        continue;
      }
      // check if duplet origin on z axis within collision region
      float zOrigin = zM - rM * cotTheta;
      //			//std::cout << "|Seeds| zM - rM * cotTheta " <<
      //zM << " - " << rM << " * " << cotTheta << std::endl;
      //        //std::cout << "|Seeds| zOrigin: " << zOrigin << "
      //        collisionRegionMin: " << m_config.collisionRegionMin <<
      //        std::endl;
      if (zOrigin < m_config.collisionRegionMin ||
          zOrigin > m_config.collisionRegionMax) {
        // std::cout << "|Seeds| !!! zOrigin < collisionRegionMin || zOrigin >
        // collisionRegionMax == TRUE !!!" << std::endl;
        continue;
      }
      // cut on the max curvature between bottom SP and interaction point
      // first transform the space point coordinates into a frame such that the
      // central space point SPm is in the origin of the frame and the x axis
      // points away from the interaction point in addition to a translation
      // transformation we also perform a rotation in order to keep the
      // curvature of the circle tangent to the x axis
      float xVal = (bottomSP->x() - spM->x()) * (spM->x() / rM) +
                   (bottomSP->y() - spM->y()) * (spM->y() / rM);
      float yVal = (bottomSP->y() - spM->y()) * (spM->x() / rM) -
                   (bottomSP->x() - spM->x()) * (spM->y() / rM);
      // std::cout << std::abs(rM) << " * " << std::abs(yVal) << " > " <<
      // -m_config.impactMax << " * " << xVal << std::endl;
      if (std::abs(rM * yVal) > -m_config.impactMax * xVal) {
        // conformal transformation u=x/(x²+y²) v=y/(x²+y²) transform the circle
        // into straight lines in the u/v plane the line equation can be
        // described in terms of aCoef and bCoef, where v = aCoef * u + bCoef
        float uB = xVal / (xVal * xVal + yVal * yVal);
        float vB = yVal / (xVal * xVal + yVal * yVal);
        // in the rotated frame the interaction point is positioned at x = -rM
        // and y ~= impactParam
        float uIP = -1. / rM;
        float vIP = m_config.impactMax / (rM * rM);
        if (yVal < 0.)
          vIP = -vIP;
        // we can obtain aCoef as the slope dv/du of the linear function,
        // estimated using du and dv between the two SP bCoef is obtained by
        // inserting aCoef into the linear equation
        float aCoef = (vB - vIP) / (uB - uIP);
        float bCoef = vIP - aCoef * uIP;
        // the distance of the straight line from the origin (radius of the
        // circle) is related to aCoef and bCoef by d^2 = bCoef^2 / (1 +
        // aCoef^2) = 1 / (radius^2) and we can apply the cut on the curvature
        if ((bCoef * bCoef) >
            (1 + aCoef * aCoef) / m_config.minHelixDiameter2) {
          // std::cout << bCoef << "^2 > (1 - " << aCoef << "^2) /" <<
          // m_config.minHelixDiameter2 << std::endl; std::cout << "|Seeds| !!!
          // impact parameter cut == TRUE !!!" << std::endl;
          continue;
        }
      }
      state.compatBottomSP.push_back(bottomSP);
      std::cout << "|Seeds| # Fill SP" << std::endl;
    }
    // no bottom SP found -> try next spM
    if (state.compatBottomSP.empty()) {
      // std::cout << "|Seeds| !!! bottom SP empty == TRUE !!!" << std::endl;
      continue;
    }

    state.linCircleBottom.clear();
    state.linCircleTop.clear();

    transformCoordinates(state.compatBottomSP, *spM, true,
                         m_config.enableCutsForSortedSP, state.linCircleBottom);
    transformCoordinates(state.compatTopSP, *spM, false,
                         m_config.enableCutsForSortedSP, state.linCircleTop);

    state.topSpVec.clear();
    state.curvatures.clear();
    state.impactParameters.clear();
    state.seedsPerSpM.clear();

    size_t numBotSP = state.compatBottomSP.size();
    size_t numTopSP = state.compatTopSP.size();
<<<<<<< Updated upstream

    size_t t0 = 0;

=======

    int m_nOneSeedsQ = 0;

    // std::cout << std::endl;
    // std::cout << "					------- Filled SPs -------" <<
    // std::endl;

    for (size_t b = 0; b < numBotSP; b++) {
      for (size_t t = 0; t < numTopSP; t++) {
        auto lb = state.linCircleBottom[b];
        auto lt = state.linCircleTop[t];
        // std::cout << std::endl;
        // std::cout << "----> rM: " << rM << " rB: " << lb.r << " rT: " << lt.r
        // << std::endl; std::cout << " xB: " << lb.x << " yB: " << lb.y << " xB:
        // " << lb.y << " xT: " << lt.x << " yT: " << lt.y << " xT: " << lt.y <<
        // std::endl;
      }
    }

    //		//std::cout << std::endl;
    //		//std::cout << "|Seeds| --- compare SP ---"<< std::endl;
    // std::cout << std::endl;
    // std::cout << "					------- Compare SPs -------" <<
    // std::endl;

    size_t t0 = 0;
>>>>>>> Stashed changes
    for (size_t b = 0; b < numBotSP; b++) {
      auto lb = state.linCircleBottom[b];
      float Zob = lb.Zo;
      float cotThetaB = lb.cotTheta;
      float Vb = lb.V;
      float Ub = lb.U;
      float ErB = lb.Er;
      float iDeltaRB = lb.iDeltaR;

      // 1+(cot^2(theta)) = 1/sin^2(theta)
      float iSinTheta2 = (1. + cotThetaB * cotThetaB);
      // calculate max scattering for min momentum at the seed's theta angle
      // scaling scatteringAngle^2 by sin^2(theta) to convert pT^2 to p^2
      // accurate would be taking 1/atan(thetaBottom)-1/atan(thetaTop) <
      // scattering
      // but to avoid trig functions we approximate cot by scaling by
      // 1/sin^4(theta)
      // resolving with pT to p scaling --> only divide by sin^2(theta)
      // max approximation error for allowed scattering angles of 0.04 rad at
      // eta=infinity: ~8.5%
      float scatteringInRegion2 = m_config.maxScatteringAngle2 * iSinTheta2;
      // multiply the squared sigma onto the squared scattering
      scatteringInRegion2 *=
          m_config.sigmaScattering * m_config.sigmaScattering;

      // clear all vectors used in each inner for loop
      state.topSpVec.clear();
      state.curvatures.clear();
      state.impactParameters.clear();
      for (size_t t = t0; t < numTopSP; t++) {
        auto lt = state.linCircleTop[t];

        std::cout << std::endl;
        std::cout << "--> rM: " << rM << " rB: " << lb.r << " rT: " << lt.r
                  << std::endl;
        std::cout << "xB: " << lb.x << " yB: " << lb.y << " zB: " << lb.z
                  << " xT: " << lt.x << " yT: " << lt.y << " zT: " << lt.z
                  << std::endl;
        std::cout << "|Seeds Map| SP: rM, rB, rT, zM, zB, zT: " << rM << " "
                  << lb.r << " " << lt.r << " " << zM << " " << lb.z << " "
                  << lt.z << std::endl;

        nSeeds_test1 += 1;

        // add errors of spB-spM and spM-spT pairs and add the correlation term
        // for errors on spM
        float error2 = lt.Er + ErB +
                       2 * (cotThetaB * lt.cotTheta * varianceRM + varianceZM) *
                           iDeltaRB * lt.iDeltaR;

        float deltaCotTheta = cotThetaB - lt.cotTheta;
        float deltaCotTheta2 = deltaCotTheta * deltaCotTheta;
        float error;
        float dCotThetaMinusError2;
<<<<<<< Updated upstream
        if (m_config.enableCutsForSortedSP) {
          // if the error is larger than the difference in theta, no need to
          // compare with scattering
          if (deltaCotTheta2 - error2 > scatteringInRegion2) {
            // break if cotThetaB < lt.cotTheta because the SP are sorted by
            // cotTheta
            if (cotThetaB - lt.cotTheta < 0) {
              break;
            }
            // since cotThetaB > lt.cotTheta and the SP are sorted by cotTheta,
            // the next bottom SP is expected to have cotThetaB > lt.cotTheta as
            // well and deltaCotTheta2 - error2 > sigmaSquaredScatteringMinPt
            t0 = t + 1;
            continue;
          }
        } else {
          if (deltaCotTheta2 - error2 > 0) {
            deltaCotTheta = std::abs(deltaCotTheta);
            // if deltaTheta larger than the scattering for the lower pT cut,
            // skip
            error = std::sqrt(error2);
            dCotThetaMinusError2 =
                deltaCotTheta2 + error2 - 2 * deltaCotTheta * error;
            // avoid taking root of scatteringInRegion
            // if left side of ">" is positive, both sides of unequality can be
            // squared
            // (scattering is always positive)
            if (dCotThetaMinusError2 > scatteringInRegion2) {
              continue;
            }
          }
        }

=======
        // if the error is larger than the difference in theta, no need to
        // compare with scattering
        // std::cout << "|Seeds| dT: " << deltaCotTheta2 << " - " << lt.Er + ErB
        // << " - " << 2 * (cotThetaB * lt.cotTheta * varianceRM + varianceZM) *
        //				iDeltaRB * lt.iDeltaR << std::endl;
        // std::cout << "|Seeds| Ert + Erb: " << lt.Er << " + " << ErB <<
        // std::endl; std::cout << "|Seeds| 3ºterm: " << 2 << " * (" << cotThetaB
        // << " * " << lt.cotTheta << " * " << varianceRM << " + " << varianceZM
        // << ") * " <<
        //				iDeltaRB << " * " << lt.iDeltaR <<
        //std::endl; std::cout << "|Seeds| ICSA: " << iSinTheta2 << " * " <<
        // 134*.05*9.*(1./std::abs(m_config.minPt))*(1./std::abs(m_config.minPt))
        // << std::endl; std::cout << "|Seeds| cotThetaB: " << cotThetaB << "
        // cotThetaT: " << lt.cotTheta << " DT: " << deltaCotTheta << std::endl;

        if (deltaCotTheta2 - error2 > scatteringInRegion2) {
          deltaCotTheta = std::abs(deltaCotTheta);
          // std::cout << "deltaCotTheta2 - error2 >
          // sigmaSquaredScatteringMinPt" << std::endl;
          // if deltaTheta larger than the scattering for the lower pT cut, skip
          error = std::sqrt(error2);
          dCotThetaMinusError2 =
              deltaCotTheta2 + error2 - 2 * deltaCotTheta * error;
          // avoid taking root of scatteringInRegion
          // if left side of ">" is positive, both sides of unequality can be
          // squared
          // (scattering is always positive)
          //          if (dCotThetaMinusError2 > scatteringInRegion2) {
          //            //std::cout << "** Continue: dCotThetaMinusError2 >
          //            scatteringInRegion2" << std::endl; continue;
          //          }
          // since the SP are sorted by cotTheta, break if cotThetaB <
          // lt.cotTheta
          if (cotThetaB - lt.cotTheta < 0) {
            // std::cout << "** Continue: cotThetaB - lt.cotTheta < 0" <<
            // std::endl;
            break;
          }
          // since cotThetaB > lt.cotTheta and the SP are sorted by cotTheta,
          // the next bottom SP is expected to have cotThetaB > lt.cotTheta as
          // well and deltaCotTheta2 - error2 > sigmaSquaredScatteringMinPt
          if (m_config.cotSeedSort) {
            t0 = t + 1;
            continue;
          }
        }

        nSeeds_test2 += 1;
				
>>>>>>> Stashed changes
        // protects against division by 0
        float dU = lt.U - Ub;
        if (dU == 0.) {
          continue;
        }
        // A and B are evaluated as a function of the circumference parameters
        // x_0 and y_0
        float A = (lt.V - Vb) / dU;
        float S2 = 1. + A * A;
        float B = Vb - A * Ub;
        float B2 = B * B;
        // sqrt(S2)/B = 2 * helixradius
        // calculated radius must not be smaller than minimum radius
        // std::cout << "|Seeds| S2: " << S2 << " B2: " << B2 << "
        // minHelixDiameter2: " << m_config.minHelixDiameter2 << std::endl;
        // std::cout << "|Seeds| Vb - A * Ub: " << Vb << " - " << A << " * " <<
        // Ub << std::endl; std::cout << "|Seeds| dT*S2" << (deltaCotTheta2 -
        // error2) * S2 << " CSA: " << iSinTheta2 *
        // 134*.05*9*2*2/(m_config.pTPerHelixRadius*m_config.pTPerHelixRadius) <<
        // " m_COFK " << 134*.05*9*1000000/(300*300)  << " iSinTheta2: " <<
        // iSinTheta2 << std::endl;
        if (S2 < B2 * m_config.minHelixDiameter2) {
          // std::cout << "** Continue:  S2 < B2 * m_config.minHelixDiameter2"
          // << std::endl;
          continue;
        }
<<<<<<< Updated upstream

        // refinement of the cut on the compatibility between the r-z slope of
        // the two seed segments using a scattering term scaled by the actual
        // measured pT
        float iHelixDiameter2 = B2 / S2;
        // calculate scattering for p(T) calculated from seed curvature
        float pT2scatterSigma = iHelixDiameter2 * m_config.sigmapT2perRadius;
        // if pT > maxPtScattering, calculate allowed scattering angle using
        // maxPtScattering instead of pt.
        float pT = m_config.pTPerHelixRadius * std::sqrt(S2 / B2) / 2.;
        if (pT > m_config.maxPtScattering) {
          float pTscatterSigma =
              (m_config.highland / m_config.maxPtScattering) *
              m_config.sigmaScattering;
          pT2scatterSigma = pTscatterSigma * pTscatterSigma;
        }
        // convert p(T) to p scaling by sin^2(theta) AND scale by 1/sin^4(theta)
        // from rad to deltaCotTheta
        float p2scatterSigma = pT2scatterSigma * iSinTheta2;
        // if deltaTheta larger than allowed scattering for calculated pT, skip
        if (m_config.enableCutsForSortedSP) {
          if (deltaCotTheta2 - error2 > p2scatterSigma) {
            if (cotThetaB - lt.cotTheta < 0) {
              break;
            }
            t0 = t;
            continue;
          }
        } else {
          if ((deltaCotTheta2 - error2 > 0) &&
              (dCotThetaMinusError2 > p2scatterSigma)) {
            continue;
          }
        }

=======
        //				if ((deltaCotTheta2 - error2) * S2 > B2 *
        //iSinTheta2 * m_config.sigmapT2perRadius) {
        //					//std::cout << "** Continue:
        //(deltaCotTheta2 - error2) * S2 > B2 * iSinTheta2 *
        //134*.05*9*1000000/(300*300)" << std::endl;
        // 					//
        // 134*0.05*9*(2/(300*2))*(2/(300*2))*1000000 ???? (2/(300*B)) --->
        // check this
        //					if (cotThetaB - lt.cotTheta < 0)
        //{
        //						//std::cout << "** Continue:
        //cotThetaB - lt.cotTheta < 0" << std::endl; 						break;
        //					}
        //					if (m_config.cotSeedSort) {
        //						t0=t;
        //						continue;
        //					}
        //				}

        nSeeds_test3 += 1;

        // refinement of the cut on the compatibility between the r-z slope of
        // the two seed segments using a scattering term scaled by the actual
        // measured pT
        float iHelixDiameter2 = B2 / S2;
        // calculate scattering for p(T) calculated from seed curvature
        float pT2scatterSigma = iHelixDiameter2 * m_config.sigmapT2perRadius;
        // if pT > maxPtScattering, calculate allowed scattering angle using
        // maxPtScattering instead of pt.
        float pT = m_config.pTPerHelixRadius * std::sqrt(S2 / B2) / 2.;
        if (pT > m_config.maxPtScattering) {
          float pTscatterSigma =
              (m_config.highland / m_config.maxPtScattering) *
              m_config.sigmaScattering;
          pT2scatterSigma = pTscatterSigma * pTscatterSigma;
        }
        // convert p(T) to p scaling by sin^2(theta) AND scale by 1/sin^4(theta)
        // from rad to deltaCotTheta
        float p2scatterSigma = pT2scatterSigma * iSinTheta2;
        // if deltaTheta larger than allowed scattering for calculated pT, skip
        if (deltaCotTheta2 - error2 > p2scatterSigma) {
          if (cotThetaB - lt.cotTheta < 0) {
            break;
          }
          t0 = t;
          continue;
        }
				
>>>>>>> Stashed changes
        // A and B allow calculation of impact params in U/V plane with linear
        // function
        // (in contrast to having to solve a quadratic function in x/y plane)
        float Im = std::abs((A - B * rM) * rM);

        //        //std::cout << "|Seeds| A: " << A << " B: " << B << std::endl;
        // std::cout << "|Seeds| Im: " << Im << " impactMax: " <<
        // m_config.impactMax << std::endl;
        if (Im <= m_config.impactMax) {
          /// **** evaluate distance the two closest-by SP in this seed
          /// candidate
          float dr = lb.iDeltaR;
          if (lt.iDeltaR < lb.iDeltaR)
            dr = lt.iDeltaR;
          /// obtain a quality score - start from the d0 estimate, and add a
          /// penalty term corresponding to how far the seed segments deviate
          /// from a straight line in r-z
          double score =
              std::abs(deltaCotTheta / (dr * std::sqrt(1. + deltaCotTheta2)));
          state.dScores.push_back(score);
          //          //std::cout << "SCORE: " << score << std::endl;

          // evaluate eta and pT of the seed
          float theta = std::atan(1. / std::sqrt(cotThetaB * lt.cotTheta));
          //					if (1. / std::sqrt(cotThetaB * lt.cotTheta) <
          //1e-8) { 						theta = std::atan(1e-8);
          //					}
          float eta = -std::log(std::tan(0.5 * theta));
          //					eta = std::asinh(std::sqrt(std::abs(cotThetaB *
          //lt.cotTheta)));
          state.etaVec.push_back(eta);
          state.ptVec.push_back(pT);

          nSeeds += 1;

          std::cout << std::endl;
          std::cout << "|Seeds Map| Seeds: rM, rB, rT, zM, zB, zT: " << rM
                    << " " << lb.r << " " << lt.r << " " << zM << " " << lb.z
                    << " " << lt.z << std::endl;

          state.topSpVec.push_back(state.compatTopSP[t]);
          // inverse diameter is signed depending if the curvature is
          // positive/negative in phi
          state.curvatures.push_back(B / std::sqrt(S2));
          state.impactParameters.push_back(Im);
<<<<<<< Updated upstream

          // evaluate eta and pT of the seed
          float theta = std::atan(1. / std::sqrt(cotThetaB * lt.cotTheta));
          float eta = -std::log(std::tan(0.5 * theta));
          state.etaVec.push_back(eta);
          state.ptVec.push_back(pT);
=======
          // std::cout << "ACEPTED" << std::endl;
          std::cout << "|Seeds Map| pT, eta, dScore, curvature, Im: "
                    << std::setprecision(10) << pT / 1000 << " " << eta << " "
                    << 0 << " " << B / std::sqrt(S2) << " " << Im << std::endl;
>>>>>>> Stashed changes
        }
      }

      //			nSeeds, zBin, phiBin: 15 1 38
      std::cout << "|size1| " << (state.seedsPerSpM).size() << std::endl;

      if (!state.topSpVec.empty()) {
        m_config.seedFilter->filterSeeds_2SpFixed(
            *state.compatBottomSP[b], *spM, state.topSpVec, state.curvatures,
            state.impactParameters, Zob, m_nOneSeedsQ,
            std::back_inserter(state.seedsPerSpM));
      }
    }

    std::cout << "|size2| " << (state.seedsPerSpM).size() << std::endl;

    nSeeds_filter1 = (state.seedsPerSpM).size();

    m_config.seedFilter->filterSeeds_1SpFixed(state.seedsPerSpM, outIt);

    std::cout << "|size3| " << (state.seedsPerSpM).size() << std::endl;

    nSeeds_filter2 = (state.seedsPerSpM).size();

    // std::cout << std::endl;
    std::cout << "|Seeds Map| nSeeds, zBin, phiBin: " << nSeeds << " " << zBin
              << " "
              << std::abs(std::ceil(spM->phi() * 1 / (2 * 3.14159265359 / 138)))
              << std::endl;
    std::cout << "|Seeds Map| nSeeds_test: " << nSeeds_test1 << " "
              << nSeeds_test2 << " " << nSeeds_test3 << std::endl;
    //		std::cout << "|Seeds Map| nSeeds_filter: " << nSeeds_filter1 << " " <<
    //nSeeds_filter2 << std::endl;

    // std::cout << "|Seeds Map| nSeeds, zBin, phiBin: " <<
    // state.seedsPerSpM.size() << " " << zBin << " " <<
    // std::ceil(spM->phi()*1/(2*3.14159265359/138)) << std::endl;
  }
}
<<<<<<< Updated upstream
=======

template <typename external_spacepoint_t, typename platform_t>
void Seedfinder<external_spacepoint_t, platform_t>::transformCoordinates(
    std::vector<const InternalSpacePoint<external_spacepoint_t>*>& vec,
    const InternalSpacePoint<external_spacepoint_t>& spM, bool bottom,
    std::vector<LinCircle>& linCircleVec) const {
  float xM = spM.x();
  float yM = spM.y();
  float zM = spM.z();
  float rM = spM.radius();

  float varianceZM = spM.varianceZ();
  float varianceRM = spM.varianceR();
  float cosPhiM = xM / rM;
  float sinPhiM = yM / rM;
  for (auto sp : vec) {
    float deltaX = sp->x() - xM;
    float deltaY = sp->y() - yM;
    float deltaZ = sp->z() - zM;
    // calculate projection fraction of spM->sp vector pointing in same
    // direction as
    // vector origin->spM (x) and projection fraction of spM->sp vector pointing
    // orthogonal to origin->spM (y)
    float x = deltaX * cosPhiM + deltaY * sinPhiM;
    float y = deltaY * cosPhiM - deltaX * sinPhiM;
    // 1/(length of M -> SP)
    float iDeltaR2 = 1. / (deltaX * deltaX + deltaY * deltaY);
    float iDeltaR = std::sqrt(iDeltaR2);
    //
    int bottomFactor = 1 * (int(!bottom)) - 1 * (int(bottom));
    // cot_theta = (deltaZ/deltaR)
    float cot_theta = deltaZ * iDeltaR * bottomFactor;
    // VERY frequent (SP^3) access
    LinCircle l;
    l.cotTheta = cot_theta;
    // location on z-axis of this SP-duplet
    l.Zo = zM - rM * cot_theta;
    l.iDeltaR = iDeltaR;
    // transformation of circle equation (x,y) into linear equation (u,v)
    // x^2 + y^2 - 2x_0*x - 2y_0*y = 0
    // is transformed into
    // 1 - 2x_0*u - 2y_0*v = 0
    // using the following m_U and m_V
    // (u = A + B*v); A and B are created later on
    l.U = x * iDeltaR2;
    l.V = y * iDeltaR2;
    // error term for sp-pair without correlation of middle space point
    l.Er = ((varianceZM + sp->varianceZ()) +
            (cot_theta * cot_theta) * (varianceRM + sp->varianceR())) *
           iDeltaR2;

    l.x = sp->x();
    l.y = sp->y();
    l.z = sp->z();
    l.r = sp->radius();

    linCircleVec.push_back(l);
    sp->setCotTheta(cot_theta);
  }
  if (m_config.cotSeedSort) {
    std::sort(vec.begin(), vec.end(),
              [](const InternalSpacePoint<external_spacepoint_t>* a,
                 const InternalSpacePoint<external_spacepoint_t>* b) -> bool {
                return (a->cotTheta() < b->cotTheta());
              });

    std::sort(linCircleVec.begin(), linCircleVec.end(),
              [](const LinCircle& a, const LinCircle& b) -> bool {
                return (a.cotTheta < b.cotTheta);
              });
  }
}

>>>>>>> Stashed changes
}  // namespace Acts

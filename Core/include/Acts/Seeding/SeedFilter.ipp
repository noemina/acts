// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <algorithm>
#include <numeric>
#include <utility>

namespace Acts {
// constructor
template <typename external_spacepoint_t>
SeedFilter<external_spacepoint_t>::SeedFilter(
    SeedFilterConfig config,
    IExperimentCuts<external_spacepoint_t>* expCuts /* = 0*/)
    : m_cfg(config.toInternalUnits()), m_experimentCuts(expCuts) {}

// function to filter seeds based on all seeds with same bottom- and
// middle-spacepoint.
// return vector must contain weight of each seed
template <typename external_spacepoint_t>
void SeedFilter<external_spacepoint_t>::filterSeeds_2SpFixed(
    const InternalSpacePoint<external_spacepoint_t>& bottomSP,
    const InternalSpacePoint<external_spacepoint_t>& middleSP,
    std::vector<const InternalSpacePoint<external_spacepoint_t>*>& topSpVec,
    std::vector<float>& invHelixDiameterVec,
    std::vector<float>& impactParametersVec, std::vector<float>& cotThetaVec,
    float zOrigin, int& nOneSeedsQ,
    std::back_insert_iterator<std::vector<std::pair<
        float, std::unique_ptr<const InternalSeed<external_spacepoint_t>>>>>
        outIt) const {
  std::cout << " tedt QQ" << nOneSeedsQ << std::endl;

  int nTopSeedConf;
  if (m_cfg.seedConfirmation) {
    float rMaxSeedConfirmation =
        std::abs(bottomSP.z()) < m_cfg.centralSeedConfirmationRange.zMaxSeedConf
            ? m_cfg.centralSeedConfirmationRange.rMaxSeedConf
            : m_cfg.forwardSeedConfirmationRange.rMaxSeedConf;
    nTopSeedConf = 2;
    if (bottomSP.radius() > rMaxSeedConfirmation)
      nTopSeedConf = 1;
  }

  size_t minWeightSeedIndex = 0;
  bool minWeightSeed = false;
  float weightMin = 1.e20;
	
	int i_n = -1; //delete
	int j_n = -1; //delete

  // TODO: pass this from the Seedfinder, adding a new config called
  // enableSortingInFilter
  bool enableSorting = m_cfg.curvatureSortingInFilter;

  // initialize original index locations
  std::vector<size_t> idx(topSpVec.size());
  std::iota(idx.begin(), idx.end(), 0);

  if (enableSorting) {
    // sort indexes based on comparing values in cotThetaVec
    std::sort(idx.begin(), idx.end(),
              [&invHelixDiameterVec](size_t i1, size_t i2) {
                return invHelixDiameterVec[i1] < invHelixDiameterVec[i2];
              });
  }

  for (auto& i : idx) {
		i_n += 1;
		
    //  for (size_t i = 0; i < topSpVec.size(); i++) {
    //		std::cout << "|seed filter| i: " << i << std::endl;
    // if two compatible seeds with high distance in r are found, compatible
    // seeds span 5 layers
    // -> very good seed
    std::vector<float> compatibleSeedR;

    float invHelixDiameter = invHelixDiameterVec[i];
    float lowerLimitCurv = invHelixDiameter - m_cfg.deltaInvHelixDiameter;
    float upperLimitCurv = invHelixDiameter + m_cfg.deltaInvHelixDiameter;
    float currentTop_r = topSpVec[i]->radius();
//    currentTop_r = std::sqrt(std::pow((topSpVec[i]->x() - bottomSP.x()), 2) +
//                             std::pow((topSpVec[i]->y() - bottomSP.y()), 2) +
//                             std::pow((topSpVec[i]->z() - bottomSP.z()), 2));

    //		**** test ****
    float deltaXt = topSpVec[i]->x() - middleSP.x();
    float deltaYt = topSpVec[i]->y() - middleSP.y();
    float deltaZt = topSpVec[i]->z() - middleSP.z();
    float cosPhiMt = middleSP.x() / middleSP.radius();
    float sinPhiMt = middleSP.y() / middleSP.radius();
    float xt = deltaXt * cosPhiMt + deltaYt * sinPhiMt;
    float yt = deltaYt * cosPhiMt - deltaXt * sinPhiMt;
    currentTop_r = std::sqrt((xt * xt) + (yt * yt) + (deltaZt * deltaZt)); // ***
    //		**** **** ****

    float impact = impactParametersVec[i];

    float a = 0;
		
		std::vector<size_t> jn;

    float weight = -(impact * m_cfg.impactWeightFactor);
		for (auto& j : idx) {
			j_n += 1;
			
      //    for (size_t j = 0; j < topSpVec.size(); j++) {
      //			std::cout << "|seed filter| j: " << j <<
      // std::endl;
      // skip it if we are looking at the same SP
      if (i == j) {
        continue;
      }
			
//			if(std::find(jn.begin(), jn.end(), j) != jn.end()) {
//				break;
//			}

      std::cout << std::endl;

      std::cout << "----------------" << std::endl;
      std::cout << "(i, j): " << i_n << ", " << j_n << std::endl;

      // compared top SP should have at least deltaRMin distance
      float otherTop_r = topSpVec[j]->radius();
//      otherTop_r = std::sqrt(std::pow((topSpVec[j]->x() - bottomSP.x()), 2) +
//                             std::pow((topSpVec[j]->y() - bottomSP.y()), 2) +
//                             std::pow((topSpVec[j]->z() - bottomSP.z()), 2));

      //		**** test ****
      float deltaXot = topSpVec[j]->x() - middleSP.x();
      float deltaYot = topSpVec[j]->y() - middleSP.y();
      float deltaZot = topSpVec[j]->z() - middleSP.z();
      float cosPhiMot = middleSP.x() / middleSP.radius();
      float sinPhiMot = middleSP.y() / middleSP.radius();
      float xot = deltaXot * cosPhiMot + deltaYot * sinPhiMot;
      float yot = deltaYot * cosPhiMot - deltaXot * sinPhiMot;
      otherTop_r = std::sqrt((xot * xot) + (yot * yot) + (deltaZot * deltaZot)); // ***
      //		**** **** ****

      float deltaR = currentTop_r - otherTop_r;

      std::cout
          << "deltaR: " << deltaR << " " << otherTop_r << " " << currentTop_r
          << std::endl;  // ---> no athena o R do top é o dR entre o top e
                         // middle, preciso colocar isso no experimental cuts
                         //			std::cout << "t: " << deltaXt <<
                         //"
                         //"
                         //<< deltaYt
      //<< " " << deltaZt << " " << xt << " " << yt << std::endl;
      // std::cout << "ot: " << deltaXot << " " << deltaYot << " " << deltaZot
      // << " " << xot
      //<< " " << yot << std::endl;
      std::cout << "invHelixDiameterVec: " << invHelixDiameterVec[j] << " "
                << lowerLimitCurv << " " << upperLimitCurv << std::endl;
      std::cout << "invHelixDiameterVec: " << invHelixDiameter << " +/- "
                << m_cfg.deltaInvHelixDiameter << std::endl;
      std::cout << "tx: " << topSpVec[i]->x() << " otx " << topSpVec[j]->x()
                << " tz " << topSpVec[i]->z() << " otz " << topSpVec[j]->z()
                << std::endl;
      std::cout << "rB: " << bottomSP.radius() << " rM " << middleSP.radius()
                << " tR " << topSpVec[i]->radius() << " otR "
                << topSpVec[j]->radius() << std::endl;

      std::cout << "cotTheta: " << topSpVec[i]->cotTheta() << " "
                << topSpVec[j]->cotTheta() << std::endl;

      // curvature difference within limits?
      // TODO: how much slower than sorting all vectors by curvature
      // and breaking out of loop? i.e. is vector size large (e.g. in jets?)
      if (invHelixDiameterVec[j] < lowerLimitCurv) {
        std::cout
            << "|seed filter 1| invHelixDiameterVec[j] < lowerLimitCurv !!!"
            << std::endl;
				jn.push_back(j);
        continue;
      }
      if (invHelixDiameterVec[j] > upperLimitCurv) {
        std::cout
            << "|seed filter 1| invHelixDiameterVec[j] > upperLimitCurv !!!"
            << std::endl;
        break;
      }
      bool newCompSeed = true;
      for (float previousDiameter : compatibleSeedR) {
        // original ATLAS code uses higher min distance for 2nd found compatible
        // seed (20mm instead of 5mm)
        // add new compatible seed only if distance larger than rmin to all
        // other compatible seeds
        if (std::abs(previousDiameter - otherTop_r) < m_cfg.deltaRMin) {
          newCompSeed = false;
          std::cout << "previousDiameter -  otherTop_r " << previousDiameter
                    << " - " << otherTop_r << std::endl;
          std::cout << "|seed filter 1| std::abs(previousDiameter - "
                       "otherTop_r) < m_cfg.deltaRMin !!!"
                    << std::endl;
          break;
        }
      }
      if (std::abs(deltaR) < m_cfg.deltaRMin) {
        newCompSeed = false;
        std::cout << "previousDiameter -  otherTop_r " << currentTop_r << " - "
                  << otherTop_r << std::endl;
        std::cout << "|seed filter 1| std::abs(deltaR) < m_cfg.deltaRMin !!!"
                  << std::endl;
        continue;
      }
      if (newCompSeed) {
        compatibleSeedR.push_back(otherTop_r);
        weight += m_cfg.compatSeedWeight;
        a += m_cfg.compatSeedWeight;
        std::cout << "a = " << a << std::endl;
      }
      if (compatibleSeedR.size() >= m_cfg.compatSeedLimit) {
        std::cout << "compatibleSeedR.size() " << compatibleSeedR.size() << "  "
                  << m_cfg.compatSeedLimit + 1 << std::endl;
        std::cout << "|seed filter 1| compatibleSeedR.size() >= "
                     "m_cfg.compatSeedLimit !!!"
                  << std::endl;
        break;
      }
    }
		j_n = -1;

    int deltaSeedConf;
    if (m_cfg.seedConfirmation) {
      deltaSeedConf = compatibleSeedR.size() + 1 - nTopSeedConf;
      std::cout << "compatibleSeedR.size(), m_nOneSeedsQ, dN, NTc "
                << compatibleSeedR.size() + 1 << "  " << nOneSeedsQ << "  "
                << deltaSeedConf << "  " << nTopSeedConf << std::endl;
      if (deltaSeedConf < 0 || (nOneSeedsQ and !deltaSeedConf)) {
        std::cout << "|seed filter 1| (dN < 0 || (m_nOneSeedsQ and !dN)) !!!"
                  << std::endl;
        continue;
      }
			bool Qm = bottomSP.radius() < 60. || std::abs(zOrigin) > 150.; // ***
      if (Qm and !deltaSeedConf and impact > 1.) {
        std::cout << "Qm, dN, impact " << Qm << "  " << deltaSeedConf << "  " << impact
                  << std::endl;
        std::cout << "|seed filter 1| (Qm and !dN and impact > 1.) !!!"
                  << std::endl;
        continue;
      }
    }
		
		

    if (m_experimentCuts != nullptr) {
      // add detector specific considerations on the seed weight
      weight += m_experimentCuts->seedWeight(bottomSP, middleSP, *topSpVec[i]);
      // discard seeds according to detector specific cuts (e.g.: weight)
      if (!m_experimentCuts->singleSeedCut(weight, bottomSP, middleSP,
                                           *topSpVec[i])) {
        std::cout << "|seed filter 1| m_experimentCuts !!!" << std::endl;
        continue;
      }
    }

    weight += std::abs(zOrigin) + m_cfg.compatSeedWeight;
    //		weight = -weight;
    std::cout << "|seed filter 1| Q: " << weight << std::endl;
    std::cout << "|seed filter 1| Q = "
              << " -1* " << m_cfg.impactWeightFactor << " * " << impact << " + "
              << std::abs(zOrigin) << " + " << a - 100 << std::endl;
    std::cout << "|seed filter 1| SP quality: " << bottomSP.quality() << " "
              << middleSP.quality() << " " << topSpVec[i]->quality() << " "
              << std::endl;

    if (weight > bottomSP.quality() and weight > middleSP.quality() and
				weight > topSpVec[i]->quality()) {
			std::cout << "|seed filter 1| quality " << std::endl;
      continue;  // ***
		}

    if (deltaSeedConf) { // ***
      ++nOneSeedsQ;
      std::cout << "|seed filter 1| dN = true " << nOneSeedsQ << " " << minWeightSeedIndex
                << std::endl;
			outIt = std::make_pair(weight, std::make_unique<const InternalSeed<external_spacepoint_t>>(bottomSP, middleSP, *topSpVec[i], zOrigin));

			std::cout << "|seed filter 1| newOneSeedQ" << std::endl;
    } else if (weight < weightMin) {
      weightMin = weight;
      minWeightSeedIndex = i;
      minWeightSeed = true;
      std::cout << "|seed filter 1| weight < weightMin " << nOneSeedsQ << " "
                << minWeightSeedIndex << " " << weightMin << std::endl;
    }
		//    if (!deltaSeedConf)
		//      continue;
	}

  if (minWeightSeed and !nOneSeedsQ) {  // weightMin < 1.e20
    outIt = std::make_pair(
        weightMin, std::make_unique<const InternalSeed<external_spacepoint_t>>(
                       bottomSP, middleSP, *topSpVec[minWeightSeedIndex], zOrigin));
    std::cout << "|seed filter 1| newOneSeed test " << nOneSeedsQ << " " << minWeightSeedIndex
              << " " << weightMin << std::endl;
  }
}

// after creating all seeds with a common middle space point, filter again
template <typename external_spacepoint_t>
void SeedFilter<external_spacepoint_t>::filterSeeds_1SpFixed(
    std::vector<std::pair<
        float, std::unique_ptr<const InternalSeed<external_spacepoint_t>>>>&
        seedsPerSpM,
    std::back_insert_iterator<std::vector<Seed<external_spacepoint_t>>> outIt)
    const {
  // sort by weight and iterate only up to configured max number of seeds per
  // middle SP
  std::sort((seedsPerSpM.begin()), (seedsPerSpM.end()),
            [](const std::pair<float, std::unique_ptr<const Acts::InternalSeed<
                                          external_spacepoint_t>>>& i1,
               const std::pair<float, std::unique_ptr<const Acts::InternalSeed<
                                          external_spacepoint_t>>>& i2) {
              if (i1.first != i2.first) {
                //                return i1.first > i2.first;
                return i1.first < i2.first;  // ********* change this **********
              } else {
                // This is for the case when the weights from different seeds
                // are same. This makes cpu & cuda results same
                float seed1_sum = 0;
                float seed2_sum = 0;
                for (int i = 0; i < 3; i++) {
                  seed1_sum += pow(i1.second->sp[i]->sp().y(), 2) +
                               pow(i1.second->sp[i]->sp().z(), 2);
                  seed2_sum += pow(i2.second->sp[i]->sp().y(), 2) +
                               pow(i2.second->sp[i]->sp().z(), 2);
                }
                return seed1_sum > seed2_sum;
              }
            });
  if (m_experimentCuts != nullptr) {
    seedsPerSpM = m_experimentCuts->cutPerMiddleSP(std::move(seedsPerSpM));
  }
  unsigned int maxSeeds = seedsPerSpM.size();
  if (maxSeeds > m_cfg.maxSeedsPerSpM) {
    maxSeeds = m_cfg.maxSeedsPerSpM + 1;
  }
  auto itBegin = seedsPerSpM.begin();
  auto it = seedsPerSpM.begin();
  // default filter removes the last seeds if maximum amount exceeded
  // ordering by weight by filterSeeds_2SpFixed means these are the lowest
  // weight seeds

  int nSeeds = 0;

  for (; it < itBegin + maxSeeds; ++it) {
    float bestSeedQuality = (*it).first;

    std::cout << "|set quality| w = " << bestSeedQuality << " "
              << (*it).second->sp[0]->x() << " " << (*it).second->sp[1]->x()
              << " " << (*it).second->sp[2]->x() << std::endl;

//		bool testQ = (bestSeedQuality < (*it).second->sp[0]->quality() ||  bestSeedQuality < (*it).second->sp[1]->quality() || bestSeedQuality < (*it).second->sp[2]->quality());
//
//		if (!testQ) continue;
		
//    if (bestSeedQuality >= (*it).second->sp[0]->quality() and
//        bestSeedQuality >= (*it).second->sp[1]->quality() and
//        bestSeedQuality >= (*it).second->sp[2]->quality()) {
//      std::cout << "|set quality| " << (*it).second->sp[0]->quality() << " "
//                << (*it).second->sp[1]->quality() << " "
//                << (*it).second->sp[2]->quality() << std::endl;
//
//			(*it).second->sp[0]->setQuality(bestSeedQuality);
//			(*it).second->sp[1]->setQuality(bestSeedQuality);
//			(*it).second->sp[2]->setQuality(bestSeedQuality);
//
//      continue;
//    }

    std::cout << "|set quality| " << (*it).second->sp[0]->quality() << " "
              << (*it).second->sp[1]->quality() << " "
              << (*it).second->sp[2]->quality() << std::endl;

    std::cout << "|set quality| Acepted" << std::endl;

    (*it).second->sp[0]->setQuality(bestSeedQuality);
    (*it).second->sp[1]->setQuality(bestSeedQuality);
    (*it).second->sp[2]->setQuality(bestSeedQuality);

    outIt = Seed<external_spacepoint_t>{
        (*it).second->sp[0]->sp(), (*it).second->sp[1]->sp(),
        (*it).second->sp[2]->sp(), (*it).second->z()};
    nSeeds += 1;
  }

  std::cout << "|Seeds Map| nSeeds_filter: " << nSeeds << " " << 0 << std::endl;
}

}  // namespace Acts

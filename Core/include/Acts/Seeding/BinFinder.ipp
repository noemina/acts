// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

template <typename external_spacepoint_t>
//<<<<<<< Updated upstream
Acts::BinFinder<external_spacepoint_t>::BinFinder(
    const std::vector<std::pair<int, int> >&& zBinNeighbors,
    const int&& numPhiNeighbors)
    : m_zBinNeighbors(std::move(zBinNeighbors)),
      m_numPhiNeighbors(std::move(numPhiNeighbors)) {}

template <typename external_spacepoint_t>
Acts::BinFinder<external_spacepoint_t>::BinFinder(
    const std::vector<std::pair<int, int> >& zBinNeighbors,
    const int& numPhiNeighbors)
    : m_zBinNeighbors(zBinNeighbors), m_numPhiNeighbors(numPhiNeighbors) {}
//=======
//Acts::BinFinder<external_spacepoint_t>::BinFinder(
//    const std::vector<std::pair<int, int> >&& zBinNeighbors,
//    const int&& numPhiNeighbors)
//    : m_zBinNeighbors(std::move(zBinNeighbors)),
//      m_numPhiNeighbors(std::move(numPhiNeighbors)) {}
//
//template <typename external_spacepoint_t>
//Acts::BinFinder<external_spacepoint_t>::BinFinder(
//    const std::vector<std::pair<int, int> >& zBinNeighbors,
//    const int& numPhiNeighbors)
//    : m_zBinNeighbors(zBinNeighbors), m_numPhiNeighbors(numPhiNeighbors) {}
//>>>>>>> Stashed changes

template <typename external_spacepoint_t>
boost::container::small_vector<size_t, 10>
Acts::BinFinder<external_spacepoint_t>::findBins(
    size_t phiBin, size_t zBin,
    const Acts::SpacePointGrid<external_spacepoint_t>* binnedSP) {
//<<<<<<< Updated upstream
  boost::container::small_vector<size_t, 9> indices;
  // if zBinNeighbors is not defined, get the indices using
  // neighborHoodIndices
  if (m_zBinNeighbors.empty()) {
    indices = binnedSP->neighborHoodIndices({phiBin, zBin}).collect();
  }
  // if the zBinNeighbors is defined, get the indices from there
  else {
    std::array<std::pair<int, int>, 2> sizePerAxis;
    sizePerAxis.at(0) = std::make_pair(-m_numPhiNeighbors, m_numPhiNeighbors);
    sizePerAxis.at(1) = m_zBinNeighbors[zBin - 1];
    indices =
        binnedSP->neighborHoodIndices({phiBin, zBin}, sizePerAxis).collect();
  }
//=======
  //    std::cout << " === dumping bin information (phi,z): " << phiBin << ", "
  //    << zBin << std::endl;

//  boost::container::small_vector<size_t, 9> indices;
//  // if the map is not defined get the indices using neighborHoodIndices
//  if (m_zBinNeighbors.empty()) {
//    indices = binnedSP->neighborHoodIndices({phiBin, zBin}).collect();
//  }
//  // if the map is defined get the indices from the map
//  else {
//    std::array<std::pair<int, int>, 2> sizePerAxis;
//    sizePerAxis.at(0) = std::make_pair(-m_numPhiNeighbors, m_numPhiNeighbors);
//    sizePerAxis.at(1) = m_zBinNeighbors[zBin - 1];
//    indices =
//        binnedSP->neighborHoodIndices({phiBin, zBin}, sizePerAxis).collect();
    ////      std::cout << "phi,z" << std::endl;
    //      // loop over the phi range defined by m_numPhiNeighbors
    //      int phiNeighborRange = m_numPhiNeighbors;
    //      for (int phiBinIndex = -phiNeighborRange; phiBinIndex <=
    //      phiNeighborRange; phiBinIndex++) {
    //        // loop over the z bins inside zBinNeighbors
    //        for (size_t zBinIndex = 0; zBinIndex < m_zBinNeighbors[zBin -
    //        1].size(); zBinIndex++) {
    //          // get z bin local index from zBinNeighbors
    //          auto zBinLocalIndex = m_zBinNeighbors[zBin - 1][zBinIndex];
    //          // get phi bin local index
    //          int maxPhiBin = (binnedSP->numLocalBins())[0];
    //          // wrap around phi
    //          size_t phiBinLocalIndex = 1 + (phiBin + phiBinIndex +
    //          (maxPhiBin-1)) % (maxPhiBin); const std::array<size_t, 2>
    //          localIndexArray = {phiBinLocalIndex,
    //                                                            zBinLocalIndex};
    //          // get the global bin index from local bin index
    //          auto globalIndex =
    //              binnedSP->globalBinFromLocalBins(localIndexArray);
    //          indices.push_back(globalIndex);
    //        }
    //      }
//  }

  //    std::cout << "indices --> (phi,z)" << std::endl;
  //    for (unsigned int i_bin = 0; i_bin < indices.size(); i_bin++){
  //        auto position = binnedSP->localBinsFromGlobalBin(indices[i_bin]);
  //        std::cout << indices[i_bin] << " --> " << position[0] << "," <<
  //        position[1] << std::endl;
  //    }
//>>>>>>> Stashed changes

  return {indices.begin(), indices.end()};
}

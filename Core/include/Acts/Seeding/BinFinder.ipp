// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// LUT for building neighbors
/* Guide for the following: 
 * z == 5: central z region, |z|<250mm
 * [-2700, -2500., -1400., -925., -450., -250., 250., 450, 925., 1400., 2500, 2700]
 *       1       2       3      4      5      6     7    8     9     10     11      z bin index
 * ---------------------------------------------------------------------------------> Z[mm]
 * Z=-2700                                  IP,Z=0                           Z=+2700
 */

// std::map < int, std::vector<int> > indices_top = { {1, {1}},
//                                                    {2, {1, 2}},
//                                                    {3, {2, 3}},
//                                                    {4, {3, 4}},
//                                                    {5, {4, 5, 6}},
//                                                    {6, {6, 7}},
//                                                    {7, {7, 8}},
//                                                    {8, {8, 9}},
//                                                    {9, {9, 10}},
//                                                    {10, {10, 11}},
//                                                    {11, {11}},
// };
// 
// std::map < int, std::vector<int> > indices_bottom = { {1, {1, 2}},
//                                                       {2, {2, 3}},
//                                                       {3, {3, 4}},
//                                                       {4, {4, 5}},
//                                                       {5, {5}},
//                                                       {6, {5, 6}},
//                                                       {7, {6, 7}},
//                                                       {8, {7, 8}},
//                                                       {9, {8, 9}},
//                                                       {10, {9, 10}},
//                                                       {11, {10, 11}},
// };



template <typename external_spacepoint_t>
std::vector<size_t> Acts::BinFinder<external_spacepoint_t>::findBins(
    size_t phiBin, size_t zBin,
    const Acts::SpacePointGrid<external_spacepoint_t>* binnedSP) {
  std::cout << " === dumping bin information: " <<  phiBin << ", " << zBin << std::endl;
  return binnedSP->neighborHoodIndices({phiBin, zBin}).collect();
}

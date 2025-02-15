// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Utilities/ProtoTracksToTracks.hpp"

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/EventDataTransforms.hpp"

#include <algorithm>

namespace ActsExamples {

PrototracksToTracks::PrototracksToTracks(Config cfg, Acts::Logging::Level lvl)
    : IAlgorithm("PrototracksToTracks", lvl), m_cfg(std::move(cfg)) {
  m_outputTracks.initialize(m_cfg.outputTracks);
  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
  m_inputProtoTracks.initialize(m_cfg.inputProtoTracks);
}

ProcessCode PrototracksToTracks::execute(const AlgorithmContext& ctx) const {
  const auto& measurements = m_inputMeasurements(ctx);

  auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
  auto mtj = std::make_shared<Acts::VectorMultiTrajectory>();
  TrackContainer tracks(trackContainer, mtj);

  const auto& prototracks = m_inputProtoTracks(ctx);
  ACTS_DEBUG("Received " << prototracks.size() << " prototracks");

  float avgSize = 0;
  std::size_t minSize = std::numeric_limits<std::size_t>::max();
  std::size_t maxSize = 0;

  for (const auto& protoTrack : prototracks) {
    if (protoTrack.empty()) {
      continue;
    }

    avgSize += static_cast<float>(protoTrack.size());
    minSize = std::min(minSize, protoTrack.size());
    maxSize = std::max(maxSize, protoTrack.size());

    auto track = tracks.makeTrack();
    for (auto measIndex : protoTrack) {
      ConstVariableBoundMeasurementProxy measurement =
          measurements.getMeasurement(measIndex);
      IndexSourceLink sourceLink(measurement.geometryId(), measIndex);

      auto trackStateProxy =
          track.appendTrackState(Acts::TrackStatePropMask::None);
      trackStateProxy.typeFlags().set(Acts::TrackStateFlag::MeasurementFlag);
      trackStateProxy.setUncalibratedSourceLink(Acts::SourceLink(sourceLink));
    }

    track.nMeasurements() = static_cast<std::uint32_t>(protoTrack.size());
    track.nHoles() = 0;
    track.nOutliers() = 0;
  }

  ConstTrackContainer constTracks{
      std::make_shared<Acts::ConstVectorTrackContainer>(
          std::move(*trackContainer)),
      std::make_shared<Acts::ConstVectorMultiTrajectory>(std::move(*mtj))};

  ACTS_DEBUG("Produced " << constTracks.size() << " tracks");
  ACTS_DEBUG("Avg track size: " << avgSize / constTracks.size());
  ACTS_DEBUG("Min track size: " << minSize << ", max track size " << maxSize);

  m_outputTracks(ctx, std::move(constTracks));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples

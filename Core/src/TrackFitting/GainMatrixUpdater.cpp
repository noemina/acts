// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/TrackFitting/GainMatrixUpdater.hpp"

#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/Utilities/Helpers.hpp"

namespace Acts {

std::tuple<double, std::error_code> GainMatrixUpdater::visitMeasurement(
    InternalTrackState trackState, NavigationDirection direction,
    LoggerWrapper logger) const {
  // default-constructed error represents success, i.e. an invalid error code
  std::error_code error;
  double chi2 = 0;

  visit_measurement(trackState.calibratedSize, [&](auto N) -> bool {
    constexpr size_t kMeasurementSize = decltype(N)::value;
    using ParametersVector = ActsVector<kMeasurementSize>;
    using CovarianceMatrix = ActsSymMatrix<kMeasurementSize>;

    typename TrackStateTraits<kMeasurementSize, true>::Measurement calibrated{
        trackState.calibrated};
    typename TrackStateTraits<kMeasurementSize, true>::MeasurementCovariance
        calibratedCovariance{trackState.calibratedCovariance};

    ACTS_INFO("Measurement dimension: " << kMeasurementSize);
    ACTS_INFO("Calibrated measurement: " << calibrated.transpose());
    ACTS_INFO("Calibrated measurement covariance:\n"
                 << calibratedCovariance);

    const auto H = trackState.projector
                       .template topLeftCorner<kMeasurementSize, eBoundSize>()
                       .eval();

    ACTS_INFO("Measurement projector H:\n" << H);

    const auto K = (trackState.predictedCovariance * H.transpose() *
                    (H * trackState.predictedCovariance * H.transpose() +
                     calibratedCovariance)
                        .inverse())
                       .eval();

    ACTS_INFO("Gain Matrix K:\n" << K);

    if (K.hasNaN()) {
      error = (direction == NavigationDirection::Forward)
                  ? KalmanFitterError::ForwardUpdateFailed
                  : KalmanFitterError::BackwardUpdateFailed;  // set to error
      return false;                                           // abort execution
    }

    trackState.filtered =
        trackState.predicted + K * (calibrated - H * trackState.predicted);
    trackState.filteredCovariance =
        (BoundSymMatrix::Identity() - K * H) * trackState.predictedCovariance;
    ACTS_INFO("Filtered parameters: " << trackState.filtered.transpose());
    ACTS_INFO("Filtered covariance:\n" << trackState.filteredCovariance);

    // calculate filtered residual
    //
    // FIXME: Without separate residual construction and assignment, we
    //        currently take a +0.7GB build memory consumption hit in the
    //        EventDataView unit tests. Revisit this once Measurement
    //        overhead problems (Acts issue #350) are sorted out.
    //
    ParametersVector residual;
    residual = calibrated - H * trackState.filtered;
    ACTS_INFO("Residual: " << residual.transpose());

    chi2 = (residual.transpose() *
            ((CovarianceMatrix::Identity() - H * K) * calibratedCovariance)
                .inverse() *
            residual)
               .value();

    ACTS_INFO("Chi2: " << chi2);
    return true;  // continue execution
  });

  return {chi2, error};
}

}  // namespace Acts

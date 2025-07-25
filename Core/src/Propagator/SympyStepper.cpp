// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/SympyStepper.hpp"

#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/Material/IVolumeMaterial.hpp"
#include "Acts/Material/Interactions.hpp"
#include "Acts/Propagator/EigenStepperError.hpp"
#include "Acts/Propagator/detail/SympyCovarianceEngine.hpp"
#include "Acts/Propagator/detail/SympyJacobianEngine.hpp"

#include <cmath>

#include "codegen/sympy_stepper_math.hpp"

namespace Acts {

SympyStepper::SympyStepper(std::shared_ptr<const MagneticFieldProvider> bField)
    : m_bField(std::move(bField)) {}

SympyStepper::SympyStepper(const Config& config) : m_bField(config.bField) {}

SympyStepper::State SympyStepper::makeState(const Options& options) const {
  State state{options, m_bField->makeCache(options.magFieldContext)};
  return state;
}

void SympyStepper::initialize(State& state,
                              const BoundTrackParameters& par) const {
  return initialize(state, par.parameters(), par.covariance(),
                    par.particleHypothesis(), par.referenceSurface());
}

void SympyStepper::initialize(State& state, const BoundVector& boundParams,
                              const std::optional<BoundMatrix>& cov,
                              ParticleHypothesis particleHypothesis,
                              const Surface& surface) const {
  FreeVector freeParams = transformBoundToFreeParameters(
      surface, state.options.geoContext, boundParams);

  state.particleHypothesis = particleHypothesis;

  state.pathAccumulated = 0;
  state.nSteps = 0;
  state.nStepTrials = 0;
  state.stepSize = ConstrainedStep();
  state.stepSize.setAccuracy(state.options.initialStepSize);
  state.stepSize.setUser(state.options.maxStepSize);
  state.previousStepSize = 0;
  state.statistics = StepperStatistics();

  state.pars = freeParams;

  // Init the jacobian matrix if needed
  state.covTransport = cov.has_value();
  if (state.covTransport) {
    // set the covariance transport flag to true and copy
    state.cov = *cov;
    state.jacToGlobal = surface.boundToFreeJacobian(
        state.options.geoContext, freeParams.segment<3>(eFreePos0),
        freeParams.segment<3>(eFreeDir0));
    state.jacobian = BoundMatrix::Identity();
    state.jacTransport = FreeMatrix::Identity();
    state.derivative = FreeVector::Zero();
  }
}

Result<std::tuple<BoundTrackParameters, BoundMatrix, double>>
SympyStepper::boundState(
    State& state, const Surface& surface, bool transportCov,
    const FreeToBoundCorrection& freeToBoundCorrection) const {
  std::optional<FreeMatrix> additionalFreeCovariance =
      state.materialEffectsAccumulator.computeAdditionalFreeCovariance(
          direction(state));
  state.materialEffectsAccumulator.reset();
  return detail::sympy::boundState(
      state.options.geoContext, surface, state.cov, state.jacobian,
      state.jacTransport, state.derivative, state.jacToGlobal,
      additionalFreeCovariance, state.pars, state.particleHypothesis,
      state.covTransport && transportCov, state.pathAccumulated,
      freeToBoundCorrection);
}

bool SympyStepper::prepareCurvilinearState(State& state) const {
  // TODO implement like in EigenStepper
  (void)state;
  return true;
}

std::tuple<BoundTrackParameters, BoundMatrix, double>
SympyStepper::curvilinearState(State& state, bool transportCov) const {
  std::optional<FreeMatrix> additionalFreeCovariance =
      state.materialEffectsAccumulator.computeAdditionalFreeCovariance(
          direction(state));
  state.materialEffectsAccumulator.reset();
  return detail::sympy::curvilinearState(
      state.cov, state.jacobian, state.jacTransport, state.derivative,
      state.jacToGlobal, additionalFreeCovariance, state.pars,
      state.particleHypothesis, state.covTransport && transportCov,
      state.pathAccumulated);
}

void SympyStepper::update(State& state, const FreeVector& freeParams,
                          const BoundVector& /*boundParams*/,
                          const Covariance& covariance,
                          const Surface& surface) const {
  state.pars = freeParams;
  state.cov = covariance;
  state.jacToGlobal = surface.boundToFreeJacobian(
      state.options.geoContext, freeParams.template segment<3>(eFreePos0),
      freeParams.template segment<3>(eFreeDir0));
}

void SympyStepper::update(State& state, const Vector3& uposition,
                          const Vector3& udirection, double qOverP,
                          double time) const {
  state.pars.template segment<3>(eFreePos0) = uposition;
  state.pars.template segment<3>(eFreeDir0) = udirection;
  state.pars[eFreeTime] = time;
  state.pars[eFreeQOverP] = qOverP;
}

void SympyStepper::transportCovarianceToCurvilinear(State& state) const {
  detail::sympy::transportCovarianceToCurvilinear(
      state.cov, state.jacobian, state.jacTransport, state.derivative,
      state.jacToGlobal, std::nullopt,
      state.pars.template segment<3>(eFreeDir0));
}

void SympyStepper::transportCovarianceToBound(
    State& state, const Surface& surface,
    const FreeToBoundCorrection& freeToBoundCorrection) const {
  detail::sympy::transportCovarianceToBound(
      state.options.geoContext, surface, state.cov, state.jacobian,
      state.jacTransport, state.derivative, state.jacToGlobal, std::nullopt,
      state.pars, freeToBoundCorrection);
}

Result<double> SympyStepper::step(State& state, Direction propDir,
                                  const IVolumeMaterial* material) const {
  double h = state.stepSize.value() * propDir;

  const double initialH = h;
  const Direction timeDirection = Direction::fromScalarZeroAsPositive(h);

  const Vector3 pos = position(state);
  const Vector3 dir = direction(state);
  const double t = time(state);
  const double qop = qOverP(state);
  const double pabs = absoluteMomentum(state);
  const double m = particleHypothesis(state).mass();
  const PdgParticle absPdg = particleHypothesis(state).absolutePdg();
  const double q = charge(state);
  const double absQ = std::abs(q);

  if (state.options.doDense && material != nullptr &&
      pabs < state.options.dense.momentumCutOff) {
    return EigenStepperError::StepInvalid;
  }

  const auto getB = [&](const double* p) -> Result<Vector3> {
    return getField(state, {p[0], p[1], p[2]});
  };

  const auto getG = [&](const double* p, double l) -> double {
    double newPabs = particleHypothesis(state).extractMomentum(l);
    if (newPabs < state.options.dense.momentumCutOff) {
      return 0.;
    }

    if (state.options.dense.meanEnergyLoss) {
      return timeDirection *
             computeEnergyLossMean(
                 MaterialSlab(material->material({p[0], p[1], p[2]}),
                              1.0f * UnitConstants::mm),
                 absPdg, m, l, absQ);
    } else {
      return timeDirection *
             computeEnergyLossMode(
                 MaterialSlab(material->material({p[0], p[1], p[2]}),
                              1.0f * UnitConstants::mm),
                 absPdg, m, l, absQ);
    }
  };

  const auto calcStepSizeScaling = [&](const double errorEstimate_) -> double {
    // For details about these values see ATL-SOFT-PUB-2009-001
    constexpr double lower = 0.25;
    constexpr double upper = 4.0;
    // This is given by the order of the Runge-Kutta method
    constexpr double exponent = 0.25;

    double x = state.options.stepTolerance / errorEstimate_;

    if constexpr (exponent == 0.25) {
      // This is 3x faster than std::pow
      x = std::sqrt(std::sqrt(x));
    } else {
      x = std::pow(x, exponent);
    }

    return std::clamp(x, lower, upper);
  };

  std::size_t nStepTrials = 0;
  double errorEstimate = 0.;

  while (true) {
    ++nStepTrials;
    ++state.statistics.nAttemptedSteps;

    // For details about the factor 4 see ATL-SOFT-PUB-2009-001
    Result<bool> res = Result<bool>::success(false);
    if (!state.options.doDense || material == nullptr) {
      res =
          rk4_vacuum(pos.data(), dir.data(), t, h, qop, m, pabs, getB,
                     &errorEstimate, 4 * state.options.stepTolerance,
                     state.pars.template segment<3>(eFreePos0).data(),
                     state.pars.template segment<1>(eFreeTime).data(),
                     state.pars.template segment<3>(eFreeDir0).data(),
                     state.derivative.data(),
                     state.covTransport ? state.jacTransport.data() : nullptr);
    } else {
      res = rk4_dense(pos.data(), dir.data(), t, h, qop, m, q, pabs, getB, getG,
                      &errorEstimate, 4 * state.options.stepTolerance,
                      state.pars.template segment<3>(eFreePos0).data(),
                      state.pars.template segment<1>(eFreeTime).data(),
                      state.pars.template segment<3>(eFreeDir0).data(),
                      state.pars.template segment<1>(eFreeQOverP).data(),
                      state.derivative.data(),
                      state.covTransport ? state.jacTransport.data() : nullptr);
    }
    if (!res.ok()) {
      return res.error();
    }
    // Protect against division by zero
    errorEstimate = std::max(1e-20, errorEstimate);

    if (*res) {
      break;
    }

    ++state.statistics.nRejectedSteps;

    const double stepSizeScaling = calcStepSizeScaling(errorEstimate);
    h *= stepSizeScaling;

    // If step size becomes too small the particle remains at the initial
    // place
    if (std::abs(h) < std::abs(state.options.stepSizeCutOff)) {
      // Not moving due to too low momentum needs an aborter
      return EigenStepperError::StepSizeStalled;
    }

    // If the parameter is off track too much or given stepSize is not
    // appropriate
    if (nStepTrials > state.options.maxRungeKuttaStepTrials) {
      // Too many trials, have to abort
      return EigenStepperError::StepSizeAdjustmentFailed;
    }
  }

  state.pathAccumulated += h;
  ++state.nSteps;
  state.nStepTrials += nStepTrials;

  ++state.statistics.nSuccessfulSteps;
  if (propDir != Direction::fromScalarZeroAsPositive(initialH)) {
    ++state.statistics.nReverseSteps;
  }
  state.statistics.pathLength += h;
  state.statistics.absolutePathLength += std::abs(h);

  const double stepSizeScaling = calcStepSizeScaling(errorEstimate);
  const double nextAccuracy = std::abs(h * stepSizeScaling);
  const double previousAccuracy = std::abs(state.stepSize.accuracy());
  const double initialStepLength = std::abs(initialH);
  if (nextAccuracy < initialStepLength || nextAccuracy > previousAccuracy) {
    state.stepSize.setAccuracy(nextAccuracy);
  }

  if (state.options.doDense &&
      (material != nullptr || !state.materialEffectsAccumulator.isVacuum())) {
    if (state.materialEffectsAccumulator.isVacuum()) {
      state.materialEffectsAccumulator.initialize(
          state.options.maxXOverX0Step, particleHypothesis(state), pabs);
    }

    Material mat =
        material != nullptr ? material->material(pos) : Material::Vacuum();

    state.materialEffectsAccumulator.accumulate(mat, propDir * h, qop,
                                                qOverP(state));
  }

  return h;
}

void SympyStepper::setIdentityJacobian(State& state) const {
  state.jacobian = BoundMatrix::Identity();
}

}  // namespace Acts

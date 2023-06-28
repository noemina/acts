// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/CylindricalContainerBuilder.hpp"
#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/DetectorBuilder.hpp"
#include "Acts/Detector/DetectorComponents.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/PortalGenerators.hpp"
#include "Acts/Detector/interface/IDetectorComponentBuilder.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/SurfaceCandidatesUpdators.hpp"
#include "Acts/Plugins/ActSVG/DetectorSvgConverter.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <fstream>
#include <memory>
#include <vector>

using namespace Acts;
using namespace Acts::Experimental;

GeometryContext tContext;

/// @brief A mockup volume builder, it generates volumes with
/// a single surface filled in in order to use the CylindricalContainerBuilder
/// infrastructure.
template <typename surface_type, typename surface_bounds_type>
class CylindricalVolumeBuilder : public IDetectorComponentBuilder {
 public:
  CylindricalVolumeBuilder(const Transform3& transform,
                           const CylinderVolumeBounds& vBounds,
                           const surface_bounds_type& sBounds,
                           const std::string& vName)
      : IDetectorComponentBuilder(),
        m_transform(transform),
        m_volumeBounds(vBounds),
        m_surfaceBounds(sBounds),
        m_name(vName) {}

  DetectorComponent construct(
      [[maybe_unused]] const GeometryContext& gctx) const final {
    // The outgoing root volumes
    std::vector<std::shared_ptr<DetectorVolume>> rootVolumes;

    // Ingredients
    auto surface = Surface::makeShared<surface_type>(
        (m_transform), std::make_shared<surface_bounds_type>(m_surfaceBounds));

    auto bounds = std::make_unique<CylinderVolumeBounds>(m_volumeBounds);
    auto portalGenerator = defaultPortalGenerator();
    auto volume = DetectorVolumeFactory::construct(
        portalGenerator, tContext, m_name, m_transform, std::move(bounds),
        {surface}, {}, tryNoVolumes(), tryAllPortalsAndSurfaces());

    // Add to the roots
    rootVolumes.push_back(volume);

    DetectorComponent::PortalContainer dContainer;
    for (auto [ip, p] : enumerate(volume->portalPtrs())) {
      dContainer[ip] = p;
    }
    return DetectorComponent{
        {volume},
        dContainer,
        RootDetectorVolumes{rootVolumes, tryRootVolumes()}};
  }

 private:
  Transform3 m_transform;
  CylinderVolumeBounds m_volumeBounds;
  surface_bounds_type m_surfaceBounds;
  std::string m_name;
};

BOOST_AUTO_TEST_SUITE(ActSvg)

BOOST_AUTO_TEST_CASE(CylindricalDetector) {
  // Declare a barrel sub builder
  auto beampipe = std::make_shared<
      CylindricalVolumeBuilder<CylinderSurface, CylinderBounds>>(
      Transform3::Identity(), CylinderVolumeBounds(0., 50., 1000.),
      CylinderBounds(25., 380.), "BeamPipe");

  // Declare a negative endcap builder
  std::vector<Acts::ActsScalar> ecNs = {-900, -700, -500, -300};
  std::vector<std::shared_ptr<const IDetectorComponentBuilder>>
      negativeEndcapsBuilders;

  for (auto [iz, nz] : Acts::enumerate(ecNs)) {
    Transform3 negZ = Transform3::Identity();
    negZ.pretranslate(Vector3(0., 0., nz));
    negativeEndcapsBuilders.push_back(
        std::make_shared<CylindricalVolumeBuilder<DiscSurface, RadialBounds>>(
            negZ, CylinderVolumeBounds(50., 140., 100.),
            RadialBounds(60., 120.), "NegativeEndcap" + std::to_string(iz)));
  }
  // Create the barrel container builder
  CylindricalContainerBuilder::Config negativeEndcapCfg;
  negativeEndcapCfg.builders = negativeEndcapsBuilders;
  negativeEndcapCfg.binning = {binZ};

  auto negativeEndcap = std::make_shared<CylindricalContainerBuilder>(
      negativeEndcapCfg, getDefaultLogger("NegativeEndcapZ", Logging::VERBOSE));

  // Declare a barrel sub builder
  auto barrel0 = std::make_shared<
      CylindricalVolumeBuilder<CylinderSurface, CylinderBounds>>(
      Transform3::Identity(), CylinderVolumeBounds(50., 80., 200.),
      CylinderBounds(65., 180.), "Barrel0");

  // Declare a barrel sub builder
  auto barrel1 = std::make_shared<
      CylindricalVolumeBuilder<CylinderSurface, CylinderBounds>>(
      Transform3::Identity(), CylinderVolumeBounds(80., 110., 200.),
      CylinderBounds(95., 180.), "Barrel1");

  // Declare a barrel sub builder
  auto barrel2 = std::make_shared<
      CylindricalVolumeBuilder<CylinderSurface, CylinderBounds>>(
      Transform3::Identity(), CylinderVolumeBounds(110., 140., 200.),
      CylinderBounds(125., 180.), "Barrel2");

  // Create the barrel container builder
  CylindricalContainerBuilder::Config barrelRCfg;
  barrelRCfg.builders = {barrel0, barrel1, barrel2};
  barrelRCfg.binning = {binR};

  auto barrel = std::make_shared<CylindricalContainerBuilder>(
      barrelRCfg, getDefaultLogger("BarrelBuilderR", Logging::VERBOSE));

  // Declare a positive endcap builder
  std::vector<Acts::ActsScalar> ecPs = {300, 500, 700, 900};
  std::vector<std::shared_ptr<const IDetectorComponentBuilder>>
      positiveEndcapsBuilders;

  for (auto [iz, pz] : Acts::enumerate(ecPs)) {
    Transform3 posZ = Transform3::Identity();
    posZ.pretranslate(Vector3(0., 0., pz));
    positiveEndcapsBuilders.push_back(
        std::make_shared<CylindricalVolumeBuilder<DiscSurface, RadialBounds>>(
            posZ, CylinderVolumeBounds(50., 140., 100.),
            RadialBounds(60., 120.), "PositiveEndcap" + std::to_string(iz)));
  }
  // Create the barrel container builder
  CylindricalContainerBuilder::Config positiveEndcapCfg;
  positiveEndcapCfg.builders = positiveEndcapsBuilders;
  positiveEndcapCfg.binning = {binZ};

  auto positiveEndcap = std::make_shared<CylindricalContainerBuilder>(
      positiveEndcapCfg, getDefaultLogger("PositiveEndcapZ", Logging::VERBOSE));

  // Create the barrel container builder
  CylindricalContainerBuilder::Config barrelEndcapCfg;
  barrelEndcapCfg.builders = {negativeEndcap, barrel, positiveEndcap};
  barrelEndcapCfg.binning = {binZ};

  auto barrelEndcap = std::make_shared<CylindricalContainerBuilder>(
      barrelEndcapCfg,
      getDefaultLogger("BarrelEndcapBuilder", Logging::VERBOSE));

  auto support = std::make_shared<
      CylindricalVolumeBuilder<CylinderSurface, CylinderBounds>>(
      Transform3::Identity(), CylinderVolumeBounds(140., 200., 1000.),
      CylinderBounds(180., 980.), "Support");

  // Create the barrel container builder
  CylindricalContainerBuilder::Config detectorCfg;
  detectorCfg.builders = {beampipe, barrelEndcap, support};
  detectorCfg.binning = {binR};

  auto containerBuilder = std::make_shared<CylindricalContainerBuilder>(
      detectorCfg, getDefaultLogger("DetectorBuilder", Logging::VERBOSE));

  // Detector builder
  Acts::Experimental::DetectorBuilder::Config dCfg;
  dCfg.auxiliary = "*** Test : Cylindrical Detector ***";
  dCfg.name = "CylindricalDetector";
  dCfg.builder = containerBuilder;

  auto detector = DetectorBuilder(dCfg).construct(tContext);

  Acts::Svg::DetectorConverter::Options detectorOptions;
  auto pDetector = Acts::Svg::DetectorConverter::convert(tContext, *detector,
                                                         detectorOptions);
  pDetector._name = detector->name();

  // Colorize in blue
  actsvg::style::color gray({{155, 155, 155}});
  actsvg::style::color pink({{255, 153, 255}});
  actsvg::style::color brown({{153, 102, 51}});
  actsvg::style::color red({{255, 0, 0}});
  actsvg::style::color green({{0, 255, 0}});
  actsvg::style::color blue({{0, 0, 255}});
  actsvg::style::color yellow({{246, 250, 5}});
  actsvg::style::color purple({{132, 5, 250}});
  actsvg::style::color turqoise({{2, 247, 170}});
  actsvg::style::color magenta({{247, 2, 178}});
  actsvg::style::color marine({{31, 1, 84}});
  actsvg::style::color black({{0, 0, 0}});
  actsvg::style::color orange({{250, 114, 2}});

  std::vector<actsvg::style::color> colors = {
      gray,   pink,     brown,   red,    green, blue,  yellow,
      purple, turqoise, magenta, marine, black, orange};
  for (auto& c : colors) {
    c._opacity = 0.1;
  }

  pDetector.colorize(colors);

  // As sheet
  auto dv_zr = Acts::Svg::View::zr(pDetector, pDetector._name);
  Acts::Svg::toFile({dv_zr}, pDetector._name + "_zr.svg");
}

BOOST_AUTO_TEST_SUITE_END()
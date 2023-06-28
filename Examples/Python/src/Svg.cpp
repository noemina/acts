// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Plugins/ActSVG/DetectorSvgConverter.hpp"
#include "Acts/Plugins/ActSVG/DetectorVolumeSvgConverter.hpp"
#include "Acts/Plugins/ActSVG/LayerSvgConverter.hpp"
#include "Acts/Plugins/ActSVG/SurfaceSvgConverter.hpp"
#include "Acts/Plugins/ActSVG/SvgUtils.hpp"
#include "Acts/Plugins/ActSVG/TrackingGeometrySvgConverter.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Io/Svg/SvgPointWriter.hpp"
#include "ActsExamples/Io/Svg/SvgTrackingGeometryWriter.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;
using namespace ActsExamples;

namespace Acts::Python {
void addSvg(Context& ctx) {
  auto [m, mex] = ctx.get("main", "examples");

  {
    auto c = py::class_<Svg::Style>(m, "SvgStyle").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Svg::Style);
    ACTS_PYTHON_MEMBER(fillColor);
    ACTS_PYTHON_MEMBER(fillOpacity);
    ACTS_PYTHON_MEMBER(highlightColor);
    ACTS_PYTHON_MEMBER(highlights);
    ACTS_PYTHON_MEMBER(strokeWidth);
    ACTS_PYTHON_MEMBER(strokeColor);
    ACTS_PYTHON_MEMBER(nSegments);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    auto c = py::class_<Svg::SurfaceConverter::Options>(m, "SvgSurfaceOptions")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Svg::SurfaceConverter::Options);
    ACTS_PYTHON_MEMBER(style);
    ACTS_PYTHON_MEMBER(templateSurface);
    ACTS_PYTHON_MEMBER(logLevel);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using DefinedStyle = std::pair<GeometryIdentifier, Svg::Style>;
    using DefinedStyleSet = std::vector<DefinedStyle>;

    auto sm = py::class_<GeometryHierarchyMap<Svg::Style>>(m, "SvgStyleMap")
                  .def(py::init<DefinedStyleSet>(), py::arg("elements"));

    auto c = py::class_<Svg::LayerConverter::Options>(m, "SvgLayerOptions")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Svg::LayerConverter::Options);
    ACTS_PYTHON_MEMBER(name);
    ACTS_PYTHON_MEMBER(surfaceStyles);
    ACTS_PYTHON_MEMBER(zRange);
    ACTS_PYTHON_MEMBER(phiRange);
    ACTS_PYTHON_MEMBER(gridInfo);
    ACTS_PYTHON_MEMBER(moduleInfo);
    ACTS_PYTHON_MEMBER(projectionInfo);
    ACTS_PYTHON_MEMBER(labelProjection);
    ACTS_PYTHON_MEMBER(labelGauge);
    ACTS_PYTHON_MEMBER(logLevel);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using DefinedLayerOptions =
        std::pair<GeometryIdentifier, Svg::LayerConverter::Options>;
    using DefinedLayerOptionsSet = std::vector<DefinedLayerOptions>;

    auto lom =
        py::class_<GeometryHierarchyMap<Svg::LayerConverter::Options>>(
            m, "SvgLayerOptionMap")
            .def(py::init<DefinedLayerOptionsSet>(), py::arg("elements"));

    auto c = py::class_<Svg::TrackingGeometryConverter::Options>(
                 m, "SvgTrackingGeometryOptions")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Svg::TrackingGeometryConverter::Options);
    ACTS_PYTHON_MEMBER(prefix);
    ACTS_PYTHON_MEMBER(layerOptions);
    ACTS_PYTHON_MEMBER(logLevel);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Writer = ActsExamples::SvgTrackingGeometryWriter;
    auto w = py::class_<Writer, std::shared_ptr<Writer>>(
                 mex, "SvgTrackingGeometryWriter")
                 .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                      py::arg("config"), py::arg("level"))
                 .def("write", py::overload_cast<const AlgorithmContext&,
                                                 const Acts::TrackingGeometry&>(
                                   &Writer::write));

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Writer::Config);
    ACTS_PYTHON_MEMBER(outputDir);
    ACTS_PYTHON_MEMBER(converterOptions);
    ACTS_PYTHON_STRUCT_END();
  }
  {
    using Writer = ActsExamples::SvgPointWriter<ActsExamples::SimSpacePoint>;
    auto w =
        py::class_<Writer, ActsExamples::IWriter, std::shared_ptr<Writer>>(
            mex, "SvgSimSpacePointWriter")
            .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"))
            .def("write",
                 py::overload_cast<const AlgorithmContext&>(&Writer::write));

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Writer::Config);
    ACTS_PYTHON_MEMBER(writerName);
    ACTS_PYTHON_MEMBER(trackingGeometry);
    ACTS_PYTHON_MEMBER(inputCollection);
    ACTS_PYTHON_MEMBER(infoBoxTitle);
    ACTS_PYTHON_MEMBER(outputDir);
    ACTS_PYTHON_STRUCT_END();
  }
  {
    using Writer =
        ActsExamples::SvgPointWriter<ActsExamples::SimHit,
                                     ActsExamples::AccessorPositionXYZ>;
    auto w =
        py::class_<Writer, ActsExamples::IWriter, std::shared_ptr<Writer>>(
            mex, "SvgSimHitWriter")
            .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"))
            .def("write",
                 py::overload_cast<const AlgorithmContext&>(&Writer::write));

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Writer::Config);
    ACTS_PYTHON_MEMBER(writerName);
    ACTS_PYTHON_MEMBER(trackingGeometry);
    ACTS_PYTHON_MEMBER(inputCollection);
    ACTS_PYTHON_MEMBER(infoBoxTitle);
    ACTS_PYTHON_MEMBER(outputDir);
    ACTS_PYTHON_STRUCT_END();
  }
  {
    mex.def("drawDetectorSvg", [](const Acts::Experimental::Detector& detector,
                                  const std::vector<std::string>& modes = {}) {
      Acts::Svg::DetectorConverter::Options detectorOptions;
      auto pDetector = Acts::Svg::DetectorConverter::convert(
          Acts::GeometryContext(), detector, detectorOptions);
      pDetector._name = detector.name();

      for (const auto& mode : modes) {
        if (mode == "zr") {
          auto dZr = Acts::Svg::View::zr(pDetector, pDetector._name);
          Acts::Svg::toFile({dZr}, pDetector._name + "_zr.svg");
        } else if (mode == "xy") {
          auto dXy = Acts::Svg::View::xy(pDetector, pDetector._name);
          Acts::Svg::toFile({dXy}, pDetector._name + "_xy.svg");
        }
      }
    });
  }

  {
    mex.def("drawDetectorVolumesSvg",
            [](const Acts::Experimental::Detector& detector,
               const std::vector<std::string>& modes = {}) {
              Acts::Svg::DetectorVolumeConverter::Options volumeOptions;

              Acts::Svg::SurfaceConverter::Options surfaceOptions;
              surfaceOptions.style.fillColor = {50, 121, 168};
              surfaceOptions.style.fillOpacity = 0.5;
              volumeOptions.surfaceOptions = surfaceOptions;

              for (const auto& volume : detector.volumes()) {
                auto [pVolume, pGrid] =
                    Acts::Svg::DetectorVolumeConverter::convert(
                        Acts::GeometryContext(), *volume, volumeOptions);

                std::string pvname =  pVolume._name;
                

                for (const auto& mode : modes) {
                  // x-y view
                  if (mode == "xy") {
                    auto volumeXY = Acts::Svg::View::xy(pVolume, pVolume._name);
                    Acts::Svg::toFile({volumeXY}, pvname + "_xy.svg");
                  }

                  // z-r view
                  if (mode == "zr") {
                    auto volumeZR = Acts::Svg::View::zr(pVolume, pVolume._name);
                    Acts::Svg::toFile({volumeZR}, pvname + "_zr.svg");
                  }

                  // x-y grid view
                  if (mode == "grid_xy" and not volume->surfaces().empty()) {
                    auto gridXY =
                        Acts::Svg::View::xy(pGrid, pVolume._name + "_grid_xy");
                    Acts::Svg::toFile({gridXY}, pvname + "_grid_xy.svg");
                  }

                  // x-y grid view
                  if (mode == "grid_zphi" and not volume->surfaces().empty()) {
                    auto gridZPhi = Acts::Svg::View::zphi(
                        pGrid, pVolume._name + "grid_zphi");
                    Acts::Svg::toFile({gridZPhi}, pvname + "_grid_zphi.svg");
                  }
                }
              }
            });
  }
}
}  // namespace Acts::Python

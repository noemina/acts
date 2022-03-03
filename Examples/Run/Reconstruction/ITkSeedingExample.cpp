// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "ActsExamples/Detector/IBaseDetector.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Geometry/CommonGeometry.hpp"
#include "ActsExamples/Io/Csv/CsvOptionsReader.hpp"
#include "ActsExamples/Io/Csv/CsvParticleReader.hpp"
#include "ActsExamples/Io/Csv/CsvSpacePointReader.hpp"
#include "ActsExamples/Io/Csv/CsvSimHitReader.hpp"
#include "ActsExamples/Io/Performance/SeedingPerformanceWriter.hpp"
#include "ActsExamples/Io/Performance/TrackFinderPerformanceWriter.hpp"
#include "ActsExamples/Io/Root/RootTrackParameterWriter.hpp"
#include "ActsExamples/MagneticField/MagneticFieldOptions.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/TrackFinding/SeedingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/SpacePointMaker.hpp"
#include "ActsExamples/TrackFinding/TrackParamsEstimationAlgorithm.hpp"
#include "ActsExamples/TruthTracking/TruthSeedSelector.hpp"
#include "ActsExamples/Utilities/Options.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <memory>

#include <boost/program_options.hpp>

#include "RecInput.hpp"

//#include "../../Tests/UnitTests/Core/Seeding/ATLASCuts.hpp"
#include "Acts/Seeding/SeedFilter.hpp"

using namespace Acts::UnitLiterals;
using namespace ActsExamples;

int main(int argc, char* argv[]) {
  
  // Setup and parse options
  auto desc = Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addOutputOptions(desc, OutputFormat::DirectoryOnly);
  Options::addInputOptions(desc);
  Options::addMagneticFieldOptions(desc);
  
  // Add specific options for this geometry
  auto vm = Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }
  
  for (const auto& it : vm) {
    std::cout << it.first.c_str() << " ";
    auto& value = it.second.value();
    if (auto v = boost::any_cast<std::string>(&value))
      std::cout << *v;
    else if (auto a = boost::any_cast<float>(&value))
      std::cout << *a;
    else if (auto b = boost::any_cast<double>(&value))
      std::cout << *b;
    else if (auto c = boost::any_cast<int>(&value))
      std::cout << *c;
    else if (auto d = boost::any_cast<bool>(&value))
      std::cout << *d;
    else if (auto e = boost::any_cast<size_t>(&value))
      std::cout << *e;
    else 
      std::cout << "wrong cast...";
    std::cout << std::endl;
  }
  
  Sequencer sequencer(Options::readSequencerConfig(vm));  
  
  // Now read the standard options
  auto logLevel = Options::readLogLevel(vm);
  auto outputDir = ensureWritableDirectory(vm["output-dir"].as<std::string>());
 
  // Setup the magnetic field
  Options::setupMagneticFieldServices(vm, sequencer);
  auto magneticField = Options::readMagneticField(vm);

  //    Read the space points and build the container
  auto spReaderCfg = setupSpacePointReading(vm, sequencer, "pixel");
  //spReaderCfg.outputSpacePoints = {"PixelSpacePoints", "StripSpacePoints", "OverlapSpacePoints"};
  spReaderCfg.outputSpacePoints = "PixelSpacePoints";


  // Seeding addAlgorithm
  SeedingAlgorithm::Config seedingCfg;
  seedingCfg.inputSpacePoints = {"PixelSpacePoints"};
  seedingCfg.outputSeeds = "PixelSeeds";
  seedingCfg.outputProtoTracks = "prototracks";
	
//	Acts::ATLASCuts<SimSpacePoint> atlasCuts = Acts::ATLASCuts<SimSpacePoint>();
//	seedingCfg.expCuts = atlasCuts;

  seedingCfg.gridConfig.rMax = 320._mm;
  seedingCfg.seedFinderConfig.rMax = seedingCfg.gridConfig.rMax;

  seedingCfg.seedFilterConfig.deltaRMin = 20_mm;
  seedingCfg.seedFinderConfig.deltaRMin = seedingCfg.seedFilterConfig.deltaRMin;
	
	seedingCfg.seedFinderConfig.deltaRMinTopSP = 6_mm;
	seedingCfg.seedFinderConfig.deltaRMinBottomSP = 6_mm;

  seedingCfg.gridConfig.deltaRMax = 280._mm;
  seedingCfg.seedFinderConfig.deltaRMax = seedingCfg.gridConfig.deltaRMax;
	
	seedingCfg.seedFinderConfig.deltaRMaxTopSP = 280_mm;
	seedingCfg.seedFinderConfig.deltaRMaxBottomSP = 120_mm;

  seedingCfg.seedFinderConfig.collisionRegionMin = -200_mm;
  seedingCfg.seedFinderConfig.collisionRegionMax = 200._mm;

  seedingCfg.gridConfig.zMin = -3000._mm;
  seedingCfg.gridConfig.zMax = 3000._mm;
  seedingCfg.seedFinderConfig.zMin = seedingCfg.gridConfig.zMin;
  seedingCfg.seedFinderConfig.zMax = seedingCfg.gridConfig.zMax;

  seedingCfg.seedFilterConfig.maxSeedsPerSpM = 4;
  seedingCfg.seedFinderConfig.maxSeedsPerSpM =
      seedingCfg.seedFilterConfig.maxSeedsPerSpM;

	seedingCfg.gridConfig.cotThetaMax = 27.2899; // 4.0 eta ---> 27.2899 = 1/np.tan(2*np.arctan(np.exp(-4)))
  seedingCfg.seedFinderConfig.cotThetaMax = seedingCfg.gridConfig.cotThetaMax;

	// number of standard deviations of Coulomb scattering angle that should be considered
	seedingCfg.seedFinderConfig.sigmaScattering = 2;
	// radiation length in Highland equation
	seedingCfg.seedFinderConfig.radLengthPerSeed = 0.09804522341059585;

  seedingCfg.gridConfig.minPt = 900._MeV;
  seedingCfg.seedFinderConfig.minPt = seedingCfg.gridConfig.minPt;

	seedingCfg.gridConfig.bFieldInZ =  1.997244311_T; // 1.997244192_T; //
  seedingCfg.seedFinderConfig.bFieldInZ = seedingCfg.gridConfig.bFieldInZ;

  seedingCfg.seedFinderConfig.beamPos = {0_mm, 0_mm};

  seedingCfg.gridConfig.impactMax = 2._mm;
  seedingCfg.seedFinderConfig.impactMax = seedingCfg.gridConfig.impactMax;
	
	seedingCfg.seedFinderConfig.maxPtScattering = 1000000._GeV;
  
	//
  seedingCfg.gridConfig.zBinEdges = {-3000., -2500., -1400., -925., -450., -250., 
                                            250., 450., 925., 1400., 2500., 3000.};
  seedingCfg.seedFinderConfig.zBinEdges = seedingCfg.gridConfig.zBinEdges;

	seedingCfg.seedFinderConfig.zBinsCustomLooping = {1, 2, 3, 4, 11, 10, 9, 8, 6, 5, 7};
	
	seedingCfg.seedFinderConfig.enableCutsForSortedSP = true;
	
  // LUT for building neighbors
  /* Guide for the following:
   * z == 6: central z region, |z|<250mm
   * [-3000, -2500., -1400., -925., -450., -250.,  250.,  450.,  925.,  1400.,  2500.,  3000]
   *       1       2       3      4      5      6      7      8      9       10      11        z bin index
   * --------------------------------------------------------------------------------------------> Z[mm]
   * Z=-3000                                  IP,Z=0                                  Z=+3000
   */
	// allows to specify the number of neighbors desired for each bin
	// {-1,1} means one neighbor on the left and one on the right
  // vector containing the map of z bins for the top SP
  seedingCfg.zBinNeighborsTop = {{0,0},{-1,0},{-1,0},{-1,0},{-1,0},{-1,1},{0,1},{0,1},{0,1},{0,1},{0,0}};
  // vector containing the map of z bins for the bottom SP
  seedingCfg.zBinNeighborsBottom = {{0,1},{0,1},{0,1},{0,1},{0,1},{0,0},{-1,0},{-1,0},{-1,0},{-1,0},{-1,0}};
	// numPhiNeighbors for the Grid means "how many phiBin neighbors (plus the current bin) should cover the full deflection of a minimum pT particle".
	// numPhiNeighbors for the BinFinder sets how many neighboring phi bins at each side of the current bin are returned
	seedingCfg.gridConfig.numPhiNeighbors = 1;
	
	
  // radial range for middle SP
  seedingCfg.seedFinderConfig.rRangeMiddleSP = {{40., 90.},{40., 200.},{46., 200.},{46., 200.},{46., 250.},{46., 250.},{46., 250.},{46., 200.},{46., 200.},{40., 200.},{40., 90.}};
	seedingCfg.seedFinderConfig.useVariableMiddleSPRange = true;
	seedingCfg.seedFinderConfig.deltaRMiddleSPRange = 10.;

	// seed confirmation
	seedingCfg.seedFinderConfig.seedConfirmation = true;
	seedingCfg.seedFilterConfig.seedConfirmation = seedingCfg.seedFinderConfig.seedConfirmation;
	
	// enable curvature sorting in SeedFilter
	seedingCfg.seedFilterConfig.curvatureSortingInFilter = true;

	// contains parameters for central seed confirmation
	seedingCfg.seedFinderConfig.centralSeedConfirmationRange = Acts::SeedConfirmationRange(250., -250., 140., 1, 2);
	seedingCfg.seedFilterConfig.centralSeedConfirmationRange = seedingCfg.seedFinderConfig.centralSeedConfirmationRange;
	// contains parameters for forward seed confirmation
	seedingCfg.seedFinderConfig.forwardSeedConfirmationRange = Acts::SeedConfirmationRange(3000., -3000., 140., 1, 2);
	seedingCfg.seedFilterConfig.forwardSeedConfirmationRange = seedingCfg.seedFinderConfig.forwardSeedConfirmationRange;
	
	seedingCfg.seedFilterConfig.impactWeightFactor = 100.;
	seedingCfg.seedFilterConfig.compatSeedWeight = 100.;
	seedingCfg.seedFilterConfig.compatSeedLimit = 3;
	
  sequencer.addAlgorithm(
      std::make_shared<SeedingAlgorithm>(seedingCfg, logLevel));

  return sequencer.run(); 
  
//   // Algorithm estimating track parameter from seed
//   TrackParamsEstimationAlgorithm::Config paramsEstimationCfg;
//   paramsEstimationCfg.inputSeeds = seedingCfg.outputSeeds;
//   paramsEstimationCfg.inputSourceLinks = hitSmearingCfg.outputSourceLinks;
//   paramsEstimationCfg.outputTrackParameters = "estimatedparameters";
//   paramsEstimationCfg.outputTrackParametersSeedMap = "estimatedparams_seed_map";
//   paramsEstimationCfg.trackingGeometry = tGeometry;
//   paramsEstimationCfg.magneticField = magneticField;
//   sequencer.addAlgorithm(std::make_shared<TrackParamsEstimationAlgorithm>(
//       paramsEstimationCfg, logLevel));
// 
//   // Seeding performance Writers
//   TrackFinderPerformanceWriter::Config tfPerfCfg;
//   tfPerfCfg.inputProtoTracks = seedingCfg.outputProtoTracks;
//   tfPerfCfg.inputParticles = inputParticles;
//   tfPerfCfg.inputMeasurementParticlesMap =
//       hitSmearingCfg.outputMeasurementParticlesMap;
//   tfPerfCfg.outputDir = outputDir;
//   tfPerfCfg.outputFilename = "performance_seeding_trees.root";
//   sequencer.addWriter(
//       std::make_shared<TrackFinderPerformanceWriter>(tfPerfCfg, logLevel));
// 
//   SeedingPerformanceWriter::Config seedPerfCfg;
//   seedPerfCfg.inputSeeds = seedingCfg.outputSeeds;
//   seedPerfCfg.inputParticles = inputParticles;
//   seedPerfCfg.inputMeasurementParticlesMap =
//       hitSmearingCfg.outputMeasurementParticlesMap;
//   seedPerfCfg.outputDir = outputDir;
//   seedPerfCfg.outputFilename = "performance_seeding_hists.root";
//   sequencer.addWriter(
//       std::make_shared<SeedingPerformanceWriter>(seedPerfCfg, logLevel));
// 
//   // The track parameters estimation writer
//   RootTrackParameterWriter::Config trackParamsWriterCfg;
//   trackParamsWriterCfg.inputSeeds = seedingCfg.outputSeeds;
//   trackParamsWriterCfg.inputTrackParameters =
//       paramsEstimationCfg.outputTrackParameters;
//   trackParamsWriterCfg.inputTrackParametersSeedMap =
//       paramsEstimationCfg.outputTrackParametersSeedMap;
//   trackParamsWriterCfg.inputParticles = particleReader.outputParticles;
//   trackParamsWriterCfg.inputSimHits = simHitReaderCfg.outputSimHits;
//   trackParamsWriterCfg.inputMeasurementParticlesMap =
//       hitSmearingCfg.outputMeasurementParticlesMap;
//   trackParamsWriterCfg.inputMeasurementSimHitsMap =
//       hitSmearingCfg.outputMeasurementSimHitsMap;
//   trackParamsWriterCfg.outputDir = outputDir;
//   trackParamsWriterCfg.outputFilename = "estimatedparams.root";
//   trackParamsWriterCfg.outputTreename = "estimatedparams";
//   sequencer.addWriter(std::make_shared<RootTrackParameterWriter>(
//       trackParamsWriterCfg, logLevel));
// 
//   return sequencer.run();
}

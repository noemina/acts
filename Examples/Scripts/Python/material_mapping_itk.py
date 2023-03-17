#!/usr/bin/env python3
import os
import argparse
from pathlib import Path


from acts.examples import (
    Sequencer,
    WhiteBoard,
    AlgorithmContext,
    ProcessCode,
    RootMaterialTrackReader,
    RootMaterialTrackWriter,
    MaterialMapping,
    JsonMaterialWriter,
    JsonFormat,
)

import acts
from acts import (
    Vector4,
    UnitConstants as u,
    SurfaceMaterialMapper,
    VolumeMaterialMapper,
    Navigator,
    Propagator,
    StraightLineStepper,
    MaterialMapJsonConverter,
)
from common import getOpenDataDetectorDirectory
from acts.examples.odd import getOpenDataDetector


def runMaterialMapping(
    trackingGeometry,
    decorators,
    outputDir,
    inputDir,
    mapName="material-map",
    mapSurface=True,
    readCachedSurfaceInformation=False,
    dumpMaterialTracks=False,
    s=None,
):
    s = s or Sequencer(numThreads=1)

    for decorator in decorators:
        s.addContextDecorator(decorator)

    wb = WhiteBoard(acts.logging.INFO)

    context = AlgorithmContext(0, 0, wb)

    for decorator in decorators:
        assert decorator.decorate(context) == ProcessCode.SUCCESS

    # Read material step information from a ROOT TTRee
    s.addReader(
        RootMaterialTrackReader(
            level=acts.logging.INFO,
            collection="material-tracks",
            fileList=[
                os.path.join(
                    inputDir,
                    mapName + "_tracks.root"
                    if readCachedSurfaceInformation
                    else "MaterialTracks_mapping-ATLAS-P2-RUN4-01-00-00.root",
                )
            ],
            readCachedSurfaceInformation=readCachedSurfaceInformation,
        )
    )

    stepper = StraightLineStepper()

    mmAlgCfg = MaterialMapping.Config(context.geoContext, context.magFieldContext)
    mmAlgCfg.trackingGeometry = trackingGeometry
    mmAlgCfg.collection = "material-tracks"

    if mapSurface:
        navigator = Navigator(
            trackingGeometry=trackingGeometry,
            resolveSensitive=False,
            resolveMaterial=True,
            resolvePassive=True
        )
        propagator = Propagator(stepper, navigator)
        mapper = SurfaceMaterialMapper(level=acts.logging.INFO, propagator=propagator)
        mmAlgCfg.materialSurfaceMapper = mapper

    jmConverterCfg = MaterialMapJsonConverter.Config(
        processSensitives=False,
        processNonMaterial=False,
        processApproaches=True,
        processRepresenting=True,
        processBoundaries=True,
        processVolumes=True,
        context=context.geoContext,
    )

    jmw = JsonMaterialWriter(
        level=acts.logging.VERBOSE,
        converterCfg=jmConverterCfg,
        fileName=os.path.join(outputDir, mapName),
        writeFormat=JsonFormat.Json,
    )

    if dumpMaterialTracks:
        s.addWriter(
            RootMaterialTrackWriter(
                level=acts.logging.INFO,
                collection=mmAlgCfg.collection,
                filePath=os.path.join(
                    outputDir,
                    mapName + "_tracks.root",
                ),
                storeSurface=True,
                storeVolume=True,
            )
        )

    mmAlgCfg.materialWriters = [jmw]

    s.addAlgorithm(MaterialMapping(level=acts.logging.INFO, config=mmAlgCfg))

    return s


if "__main__" == __name__:
    p = argparse.ArgumentParser(
        description="Script to run material mapping on ITk geometry"
    )
    p.add_argument(
        "geo_dir",
        help="Input directory containing the ITk standalone geometry. Get in touch if you don't have this.",
    )
    p.add_argument(
        "--material", type=str, default = "", help="Geometry file to define layers used in material mapping"
    )

    args = p.parse_args()

    geo_example_dir = Path(args.geo_dir)
    assert geo_example_dir.exists(), "Detector example input directory missing"
    assert os.path.exists(args.material), "Invalid file path/name in --material. Please check your input!"

    from acts.examples.itk import buildITkGeometry

    detector, trackingGeometry, decorators = buildITkGeometry(
        geo_example_dir,
        customMaterialFile=args.material
    )

    runMaterialMapping(
        trackingGeometry,
        decorators,
        outputDir=os.getcwd(),
        inputDir=os.getcwd(),
        readCachedSurfaceInformation=False).run()

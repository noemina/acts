import acts
import math
import acts.examples


from acts import logging, Binning, ProtoBinning, LayerStructureBuilder, VolumeStructureBuilder, DetectorVolumeBuilder, DetectorBuilder, CylindricalContainerBuilder, Transform3

from detector import *

jsOptions = acts.examples.SurfaceJsonOptions()
jsOptions.inputFile = '/home/noemi/Work/ACTS/detector/itk-geometry-map.json'

# Where to pick the surfaces from

surfacesHierarchyMap = acts.examples.readSurfaceFromJson(jsOptions)
svlMap = acts.examples.extractVolumeLayerSurfaces(surfacesHierarchyMap, True)

# Create the detector structure
detectorRmin = 0.
detectorRmax = 1148.0
detectorZmin = -3100
detectorZmax = -detectorZmin
beamPipeRmax = 27.


# Beam pipe section
beamPipe = EmptyVolume(
    'BeamPipe', [detectorRmin, beamPipeRmax, -detectorZmax, detectorZmax])

# # IPT
# #
# iptRmax = 29.5
# ipt = EmptyVolume('IPT', [beamPipeRmax, iptRmax, -detectorZmax, detectorZmax])
#
# # Pixel section
# #
# innerPixelRmax = 124.
# innerPixelZmid = 253.
#
# # Inner Pixel negative/positive endcap
# # binning r / binning phi
# coupledRingsBinningR = acts.ProtoBinning(
#     acts.Binning.r, acts.Binning.bound, 30., 124., 2, 1)
# coupledRingsBinningPhi = acts.ProtoBinning(
#     acts.Binning.phi, acts.Binning.closed, -math.pi, math.pi, 20, 1)
#
# innerRingsBinningR = acts.ProtoBinning(
#     acts.Binning.r, acts.Binning.bound, 55., 80., 1, 1)
# innerRingsBinningPhi = acts.ProtoBinning(
#     acts.Binning.phi, acts.Binning.closed, -math.pi, math.pi, 30, 1)
#
# outerRingsBinningR = acts.ProtoBinning(
#     acts.Binning.r, acts.Binning.bound, 75., 124., 1, 1)
# outerRingsBinningPhi = acts.ProtoBinning(
#     acts.Binning.phi, acts.Binning.closed, -math.pi, math.pi, 20, 1)
#
# coupledRingsBinning = [coupledRingsBinningR, coupledRingsBinningPhi]
# innerRingsBinning = [innerRingsBinningR, innerRingsBinningPhi]
# outerRingsBinning = [outerRingsBinningR, outerRingsBinningPhi]
#
# innerPixelEndcapPositions = [263., 291., 322., 357., 396., 437., 486., 543., 604., 675., 749.,
#                              835., 925., 1026., 1103., 1142., 1229., 1272., 1359., 1403., 1503.,
#                              1553., 1665., 1721., 1846., 1909., 2120., 2357., 2621.]
# endcapBinnings = 14*[coupledRingsBinning] + [innerRingsBinning] + [coupledRingsBinning] \
#     + 5*[innerRingsBinning,outerRingsBinning] +  3*[outerRingsBinning]
# innerPixelEndcapVolumeIds = [8, 10]
# innerPixelEndcapLayerIdOffsets = [2, 2]
#
# innerPixelEndcaps = createEndcaps('InnerPixel::Endcap',
#                                   [iptRmax, innerPixelRmax, innerPixelZmid, detectorZmax],
#                                   svlMap, innerPixelEndcapVolumeIds, innerPixelEndcapLayerIdOffsets,
#                                   innerPixelEndcapPositions, 10.,
#                                   endcapBinnings)
#
# # Inner Pixel barrel
# innerPixelBarrelPositions = [34., 99]
# innerPixelBarrelBinningZ = acts.ProtoBinning(
#     acts.Binning.z, acts.Binning.bound, -245., 245., 12, 1)
# innerPixelBarrel = createBarrel('InnerPixel::Barrel', [iptRmax, innerPixelRmax, -innerPixelZmid, innerPixelZmid],
#                                 svlMap, 9, innerPixelBarrelPositions, 10., innerPixelBarrelBinningZ, [12, 20])
#
# innerPixelBarrel.printDimensions()
#
# # Inner Pixel overall
# innerPixelBuilder = ContainerStructure('InnerPixel',
#                                        [iptRmax, innerPixelRmax, detectorZmin, detectorZmax],
#                                        [innerPixelEndcaps[0], innerPixelBarrel, innerPixelEndcaps[1]], None)
# # IST
# #
# istRmax = 142.5
# ist = EmptyVolume('IST', [innerPixelRmax, istRmax, -detectorZmax, detectorZmax])
#
# # Outer Pixel negative/positive endcap
# # binning r / binning phi
#
# outerPixelZmid = 376.
# #Ring 0
# outerPixelRmaxRing0 = 205.
# outerPixelEndcapBinningR = acts.ProtoBinning(
#     acts.Binning.r, acts.Binning.bound, 150., 200., 1, 1)
# outerPixelEndcapBinningPhi = acts.ProtoBinning(
#     acts.Binning.phi, acts.Binning.closed, -math.pi, math.pi, 32, 1)
# outerPixelEndcapPositions = [404.85, 482.6, 579.6, 700.5, 851.2, 1039.14, 1145.5, 1249.0, 1365.,
#                              1492.0, 1633.0, 1789.0, 1961.0, 2151.0, 2361.0, 2593.0, 2850.]
# outerPixelEndcapVolumeIds = [13, 18]
# outerPixelEndcapLayerIdOffsets = [2, 2]
#
# outerPixelEndcapsRing0 = createEndcaps('OuterPixel::EndcapRing0',
#                                        [istRmax, outerPixelRmaxRing0, outerPixelZmid, detectorZmax],
#                                        svlMap, outerPixelEndcapVolumeIds, outerPixelEndcapLayerIdOffsets,
#                                        outerPixelEndcapPositions, 15.,
#                                        len(outerPixelEndcapPositions)*[[outerPixelEndcapBinningR, outerPixelEndcapBinningPhi]])
#
# outerPixelRmaxRing1 = 265.
# outerPixelEndcapBinningR = acts.ProtoBinning(
#     acts.Binning.r, acts.Binning.bound, 210., 260., 1, 1)
# outerPixelEndcapBinningPhi = acts.ProtoBinning(
#     acts.Binning.phi, acts.Binning.closed, -math.pi, math.pi, 44, 1)
# outerPixelEndcapPositions = [400.6, 457.0, 522.1, 597.4, 684.5, 785.2, 901.7, 1036.4, 1145.5,
#                              1297.0, 1473.0, 1676.0, 1911.0, 2180.0, 2491.0, 2850.0]
# outerPixelEndcapVolumeIds = [14, 19]
# outerPixelEndcapLayerIdOffsets = [2, 2]
#
# outerPixelEndcapsRing1 = createEndcaps('OuterPixel::EndcapRing1',
#                                        [outerPixelRmaxRing0, outerPixelRmaxRing1, outerPixelZmid, detectorZmax],
#                                        svlMap, outerPixelEndcapVolumeIds, outerPixelEndcapLayerIdOffsets,
#                                        outerPixelEndcapPositions, 15.,
#                                        len(outerPixelEndcapPositions)*[[outerPixelEndcapBinningR, outerPixelEndcapBinningPhi]])
#
# outerPixelRmax = 320.
# outerPixelEndcapBinningR = acts.ProtoBinning(
#     acts.Binning.r, acts.Binning.bound, 270., 320., 1, 1)
#
# outerPixelInclinedEndcapBinningPhi = acts.ProtoBinning(
#     acts.Binning.phi, acts.Binning.closed, -math.pi, math.pi, 56, 1)
# outerPixelRingEndcapBinningPhi = acts.ProtoBinning(
#     acts.Binning.phi, acts.Binning.closed, -math.pi, math.pi, 52, 1)
#
# outerPixelInclinedEndcapBinning = [outerPixelEndcapBinningR, outerPixelInclinedEndcapBinningPhi]
# outerPixelRingEndcapBinning = [outerPixelEndcapBinningR, outerPixelRingEndcapBinningPhi]
#
# outerPixelEndcapPositions = [397.5, 449.9, 508.5, 573.9, 646.9, 728.6, 819.7, 921.6, 1035.4,
#                              1145.5, 1277.0, 1427.0, 1597.0, 1789.0, 2007.0, 2253.0, 2533.0, 2850.0]
# outerPixelEndcapVolumeIds = [15, 20]
# outerPixelEndcapLayerIdOffsets = [2, 2]
#
# outerPixelEndcapsRing2 = createEndcaps('OuterPixel::EndcapRing2',
#                                        [outerPixelRmaxRing1, outerPixelRmax, outerPixelZmid, detectorZmax],
#                                        svlMap, outerPixelEndcapVolumeIds, outerPixelEndcapLayerIdOffsets,
#                                        outerPixelEndcapPositions, 15.,
#                                        9*[outerPixelInclinedEndcapBinning] + 9*[outerPixelRingEndcapBinning])
#
# # Outer Pixel barrel
# outerPixelBarrelPositions = [160., 228., 291.]
# outerPixelBarrelBinningZ = acts.ProtoBinning(
#     acts.Binning.z, acts.Binning.bound, -376., 376., 18, 1)
# outerPixelBarrel = createBarrel('OuterPixel::Barrel', [istRmax, outerPixelRmax, -outerPixelZmid, outerPixelZmid],
#                                 svlMap, 16, outerPixelBarrelPositions, 10., outerPixelBarrelBinningZ, [32, 44, 56])
#
# outerPixelBarrel.printDimensions()
#
# # Outer Pixel overall
# outerPixelBuilder = ContainerStructure('OuterPixel',
#                                        [istRmax, outerPixelRmax, detectorZmin, detectorZmax],
#                                        [outerPixelEndcapsRing0[0], outerPixelEndcapsRing1[0], outerPixelEndcapsRing2[0], outerPixelBarrel, outerPixelEndcapsRing0[1], outerPixelEndcapsRing1[1], outerPixelEndcapsRing2[1]], None)
#
#
# # PST
# #
pstRmax = 346.7
# pst = EmptyVolume('PST', [outerPixelRmax, pstRmax, -detectorZmax, detectorZmax])

# Strip negative/positive endcap
# binning r / binning phi
# stripEndcapBinningR = acts.ProtoBinning(
#     acts.Binning.r, acts.Binning.bound, 380., 970., 6, 1)
# stripEndcapBinningPhi = acts.ProtoBinning(
#     acts.Binning.phi, acts.Binning.closed, -math.pi, math.pi, 32, 1)
# stripEndcapPositions = [1512., 1702., 1952., 2237., 2532., 2850.]
# stripEndcapVolumeIds = [22, 24]
# stripEndcapLayerIdOffsets = [2, 2]
#
stripZmid=1375.
#
# stripEndcaps = createEndcaps('Strip::Endcap',
#                              [pstRmax, detectorRmax, stripZmid, detectorZmax],
#                              svlMap, stripEndcapVolumeIds, stripEndcapLayerIdOffsets,
#                              stripEndcapPositions, 25.,
#                              len(stripEndcapPositions)*[[stripEndcapBinningR, stripEndcapBinningPhi]])

# # Strip barrel
# stripBarrelPositions = [399., 562., 762., 1000.]
# stripBarrelBinningZ = acts.ProtoBinning(
#     acts.Binning.z, acts.Binning.bound, -1370., 1370., 28, 1)
# stripBarrel = createBarrel('Strip::Barrel', [pstRmax, detectorRmax, -stripZmid, stripZmid],
#                            svlMap, 23, stripBarrelPositions, 25., stripBarrelBinningZ,
#                            [28, 40, 56, 72])
#
# stripBarrel.printDimensions()

# Strip overall
# stripBuilder = ContainerStructure('Strip',
#                                   [pstRmax, detectorRmax, detectorZmin, detectorZmax],
#                                   [stripEndcaps[0], stripBarrel, stripEndcaps[1]], None)

stripBarrelPositions = [399.]
stripBarrelBinningZ = acts.ProtoBinning(
    acts.Binning.z, acts.Binning.bound, -1370., 1370., 28, 1)
stripBarrel = createBarrel('Strip::Barrel', [pstRmax, detectorRmax, -stripZmid, stripZmid],
                           svlMap, 23, stripBarrelPositions, 25., stripBarrelBinningZ,
                           [28])

stripBarrel.printDimensions()

stripBuilder = ContainerStructure('Strip',
                                  [pstRmax, detectorRmax, detectorZmin, detectorZmax],
                                  [stripBarrel], None)


detectorBuilder = ContainerStructure('Detector',
                                     [detectorRmin, detectorRmax, detectorZmin, detectorZmax],
                                     # [beamPipe, ipt, innerPixelBuilder, ist, outerPixelBuilder, pst, stripBuilder],
                                     [beamPipe, stripBuilder],
                                     None)

# detectorBuilder = ContainerStructure('Detector',
#                                      [detectorRmin, detectorRmax, detectorZmin, detectorZmax],
#                                      [beamPipe, ipt, innerPixelBuilder],
#                                      None)


detectorBuilderConf = acts.DetectorBuilder.Config()
detectorBuilderConf.name = 'ITk'
detectorBuilderConf.builder = detectorBuilder.getBuilder()

detectorBuilder = acts.DetectorBuilder(
    detectorBuilderConf, 'ITk Detector Builder', acts.logging.VERBOSE)

detector = detectorBuilder.construct(acts.GeometryContext())

acts.examples.drawDetectorSvg(detector,  ['zr', 'xy'])

acts.examples.drawDetectorVolumesSvg(detector, ['grid_zphi', 'grid_xy'])



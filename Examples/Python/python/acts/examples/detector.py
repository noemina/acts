import acts
import math
import acts.examples


from acts import logging, Binning, ProtoBinning, LayerStructureBuilder, VolumeStructureBuilder, DetectorVolumeBuilder, DetectorBuilder, CylindricalContainerBuilder, Transform3


class LayerVolume:

    def __init__(self, name, dimensions, surfaces, binnings, supports):
        """Creates a layer structure with the given dimensions, surfaces, binnings and supports.

        Parameters
        ----------
        dimensions : list of float, rmin, rmax, zmin, zmax
           The dimensions of the layer structure.
        surfaces : list of Surface
            The surfaces of the layer structure.
        binnings : list of int
            The binnings of the layer structure.
        supports : list of Surface
            The supports of the layer structure.

        Returns
        -------
        LayerStructure
            The created layer structure.
        """
        self.__name = name
        self.__dimensions = dimensions
        self.__surfaces = surfaces
        self.__binnings = binnings
        self.__supports = supports

    def print(self):
        print("LayerVolume: ", self.__name)
        print("  dimensions: ", self.__dimensions)
        print("  surfaces: ", self.__surfaces)
        print("  binnings: ", self.__binnings)
        print("  supports: ", self.__supports)

    def getDimensions(self):
        """Returns the dimensions of the layer structure.

        Returns
        -------
        list of float
            The dimensions of the layer structure.
        """
        return self.__dimensions

    def setDimensions(self, dimensions):
        """Sets the dimensions of the layer structure.

        Parameters
        ----------
        dimensions : list of float, rmin, rmax, zmin, zmax
           The dimensions of the layer structure.
        """
        self.__dimensions = dimensions

    def getSurfaces(self):
        """Returns the surfaces of the layer structure.

        Returns
        -------
        list of Surface
            The surfaces of the layer structure.
        """
        return self.__surfaces

    def getBinnings(self):
        """Returns the binnings of the layer structure.

        Returns
        -------
        list of int
            The binnings of the layer structure.
        """
        return self.__binnings

    def getSupports(self):
        """Returns the supports of the layer structure.

        Returns
        -------
        list of Surface
            The supports of the layer structure.
        """
        return self.__supports

    def setName(self, name):
        self.__name = name

    def getName(self):
        return self.__name

    def getBuilder(self):
        layerStructureConfig = acts.LayerStructureBuilder.Config()
        layerStructureConfig.surfacesProvider = self.__surfaces
        layerStructureConfig.binnings = self.__binnings

        layerStructureBuilder = acts.LayerStructureBuilder(
            layerStructureConfig, self.__name + 'LayerBuilder', acts.logging.VERBOSE)

        # Volume externals
        volumeStructureConfig = acts.VolumeStructureBuilder.Config()
        volumeStructureConfig.boundsType = acts.VolumeBoundsType.Cylinder
        volumeStructureConfig.boundValues = [self.__dimensions[0], self.__dimensions[1], 0.5*(
            self.__dimensions[3]-self.__dimensions[2]), math.pi, 0.]
        zPos = 0.5*(self.__dimensions[2] + self.__dimensions[3])
        if zPos != 0.:
            volumeStructureConfig.transform = acts.Transform3([0, 0, zPos])

        volumeStructureBuilder = acts.VolumeStructureBuilder(
            volumeStructureConfig, self.__name + 'ShapeBuilder', acts.logging.VERBOSE)

        # Volume builder
        volumeBuilderConfig = acts.DetectorVolumeBuilder.Config()
        volumeBuilderConfig.name = self.__name
        volumeBuilderConfig.internalsBuilder = layerStructureBuilder
        volumeBuilderConfig.externalsBuilder = volumeStructureBuilder

        volumeBuilder = acts.DetectorVolumeBuilder(
            volumeBuilderConfig, self.__name + 'VolumeBuilder', acts.logging.VERBOSE)

        return volumeBuilder


class EmptyVolume:
    def __init__(self, name, dimensions):
        """Creates an empty volume with the given dimensions.

        Parameters
        ----------
        dimensions : list of float, rmin, rmax, zmin, zmax
           The dimensions of the empty volume.

        Returns
        -------
        EmptyVolume
            The created empty volume.
        """
        self.__dimensions = dimensions
        self.__name = name

    def getDimensions(self):
        """Returns the dimensions of the empty volume.

        Returns
        -------
        list of float
            The dimensions of the empty volume.
        """
        return self.__dimensions

    def getName(self):
        return self.__name

    def getBuilder(self):
        # Volume externals
        volumeStructureConfig = acts.VolumeStructureBuilder.Config()
        volumeStructureConfig.boundsType = acts.VolumeBoundsType.Cylinder
        volumeStructureConfig.boundValues = [self.__dimensions[0], self.__dimensions[1], 0.5*(
            self.__dimensions[3]-self.__dimensions[2]), math.pi, 0.]
        zPos = 0.5*(self.__dimensions[2] + self.__dimensions[3])
        if zPos != 0.:
            volumeStructureConfig.transform = acts.Transform3([0, 0, zPos])

        volumeStructureBuilder = acts.VolumeStructureBuilder(
            volumeStructureConfig, self.__name + 'ShapeBuilder', acts.logging.VERBOSE)

        # Volume builder
        volumeBuilderConfig = acts.DetectorVolumeBuilder.Config()
        volumeBuilderConfig.name = self.__name
        volumeBuilderConfig.externalsBuilder = volumeStructureBuilder

        volumeBuilder = acts.DetectorVolumeBuilder(
            volumeBuilderConfig, self.__name + 'VolumeBuilder', acts.logging.VERBOSE)

        return volumeBuilder

    def print(self):
        print("EmptyVolume")
        print("  dimensions: ", self.__dimensions)


class ContainerStructure:
    def __init__(self, name, dimensions, layers, ordering):
        """Creates a barrel structure with the given dimensions and layers.

        Parameters
        ----------
        dimensions : list of float, rmin, rmax, zmin, zmax
           The dimensions of the barrel structure.
        layers : list of LayerStructure
            The layers of the barrel structure.
        ordering : is an acts.BinningValue

        Returns
        -------
        BarrelStructure
            The created barrel structure.
        """
        self.__name = name
        self.__dimensions = dimensions
        self.__layers = layers
        self.__volumes = []

        if (ordering is acts.Binning.z):
            self.harmonizeInZ()

        elif (ordering is acts.Binning.r):
            self.harmonizeInR()

    def harmonizeInR(self):
        """Harmonizes the internal layer structures to fit into the dimension.
        """
        """ Gap volumes are introduced to fit the overal volume dimension
        """
        volumeRmin = self.__dimensions[0]
        volumeRmax = self.__dimensions[1]
        volumeZmin = self.__dimensions[2]
        volumeZmax = self.__dimensions[3]
        lastRmax = volumeRmin

        ginx = 0
        for layer in self.__layers:
            lDim = layer.getDimensions()
            if lDim[0] > lastRmax:
                self.__volumes += [EmptyVolume(self.__name + '::Gap'+str(ginx), [lastRmax,
                                                                                 lDim[0], volumeZmin, volumeZmax])]
            layer.setDimensions([lDim[0], lDim[1], volumeZmin, volumeZmax])
            layer.setName(self.__name + '::' + layer.getName())
            self.__volumes += [layer]
            lastRmax = lDim[1]
            ginx = ginx + 1

        if lastRmax < volumeRmax:
            self.__volumes += [EmptyVolume(self.__name + '::Gap'+str(ginx), [lastRmax,
                                                                             volumeRmax, volumeZmin, volumeZmax])]
        return self

    def harmonizeInZ(self):
        """Harmonizes the internal layer structures to fit into the dimension.
        """
        """ Gap volumes are introduced to fit the overal volume dimension
        """
        volumeRmin = self.__dimensions[0]
        volumeRmax = self.__dimensions[1]
        volumeZmin = self.__dimensions[2]
        volumeZmax = self.__dimensions[3]

        lastZmax = volumeZmin
        ginx = 0
        for layer in self.__layers:
            lDim = layer.getDimensions()
            if lDim[2] > lastZmax:
                self.__volumes += [EmptyVolume(self.__name + '::Gap'+str(ginx), [volumeRmin,
                                                                                 volumeRmax, lastZmax, lDim[2]])]
            layer.setDimensions([volumeRmin, volumeRmax, lDim[2], lDim[3]])
            layer.setName(self.__name + '::' + layer.getName())
            self.__volumes += [layer]
            lastZmax = lDim[3]
            ginx = ginx + 1

        if lastZmax < volumeZmax:
            self.__volumes += [EmptyVolume(self.__name + '::Gap'+str(ginx), [volumeRmin,
                                                                             volumeRmax, lastZmax, volumeZmax])]
        return self

    def getName(self):
        return self.__name

    def getBuilder(self):

        volumeBuilders = []
        if len(self.__volumes) == 0:
            self.__volumes = self.__layers
        for volume in self.__volumes:
            volumeBuilders += [volume.getBuilder()]

        containerConfig = acts.CylindricalContainerBuilder.Config()
        containerConfig.builders = volumeBuilders

        containerBuilder = acts.CylindricalContainerBuilder(
            containerConfig, self.__name + 'Builder', acts.logging.VERBOSE)
        return containerBuilder

    def printDimensions(self):
        print("ContainerStructure: ", self.__name)
        for volume in self.__volumes:
            print(" * ", volume.getName(), " - ", volume.getDimensions())

# Pixel barrel


def createBarrel(name, dimensions, surfacesMap, volumeID, layerPositions, layerHalfR, binningInZ, modulesInPhi):
    # The barrel detector
    barrels = []
    for i, r in enumerate(layerPositions):
        rMin = r - layerHalfR
        rMax = r + layerHalfR
        if rMin < dimensions[0]:
            rMin = dimensions[0]
        if rMax > dimensions[1]:
            rMax = dimensions[1]

        binningInR = acts.ProtoBinning(
            acts.Binning.phi, acts.Binning.closed, -math.pi, math.pi, modulesInPhi[i], 1)

        barrelLayers = LayerVolume('L'+str(i),
                                   [rMin, rMax], acts.LayerStructureBuilder.SurfacesHolder(
                                       surfacesMap[volumeID][(i+1)*2]),
                                   [binningInZ, binningInR], [])
        barrels += [barrelLayers]

    barrel = ContainerStructure(name, dimensions, barrels, acts.Binning.r)
    barrel.printDimensions()
    return barrel

# Convencience function


def createEndcaps(name, dimensions, surfacesMap, endcapVolumeIds, endcapLayerIdOffsets, endcapPositions, layerHalfZ, endcapBinnings):
    """Creates a barrel structure with the given dimensions and layers.

    """
    endcapRmin, endcapRmax, endcapZmid, endcapZmax = dimensions
    endcaps = []
    # Helpers
    sign = 1
    tags = ['N', 'P']

    for it, t in enumerate(tags):
        endcapLayers = []
        endcapPositions.reverse()
        endcapBinnings.reverse()
        sign = -sign
        endcapDimensionsZ = [sign*endcapZmid, sign*endcapZmax]
        endcapDimensionsZ.sort()
        offset = 0
        if sign < 0:
            offset = len(endcapPositions)-1
        # build the layer loops
        for i in range(0, len(endcapPositions)):
            #  calculate z min / z max
            zMin = sign*endcapPositions[i] - layerHalfZ
            zMax = sign*endcapPositions[i] + layerHalfZ
            # create the layer
            volID = endcapVolumeIds[it]
            layerID = endcapLayerIdOffsets[it] + i * 2

            endcapLayers += [LayerVolume(t+str(offset+sign*i),
                                         [endcapRmin, endcapRmax, zMin, zMax],
                                         acts.LayerStructureBuilder.SurfacesHolder(
                surfacesMap[volID][layerID]),
                endcapBinnings[i], [])]

        endcap = ContainerStructure(name+t,
                                    [endcapRmin, endcapRmax, endcapDimensionsZ[0], endcapDimensionsZ[1]], endcapLayers, acts.Binning.z)
        endcap.printDimensions()
        endcaps += [endcap]
    return endcaps


# # # jsOptions = acts.examples.SurfaceJsonOptions()
# # # jsOptions.inputFile = 'odd-geometry-map.json'
# # #
# # # # Where to pick the surfaces from
# # #
# # # surfacesHierarchyMap = acts.examples.readSurfaceFromJson(jsOptions)
# # # svlMap = acts.examples.extractVolumeLayerSurfaces(surfacesHierarchyMap, True)
# # #
# # # # Create the detector structure
# # # detectorRmin = 0.
# # # detectorRmax = 1200
# # # detectorZmin = -3100
# # # detectorZmax = -detectorZmin
# # # beamPipeRmax = 27.
# # #
# # #
# # # # Beam pipe section
# # # beamPipe = EmptyVolume(
# # #     'BeamPipe', [detectorRmin, beamPipeRmax, -detectorZmax, detectorZmax])
# # #
# # #
# # # # Pixel section
# # # #
# # # pixelRmax = 200.
# # # pixelZmid = 580.
# # #
# # # # Pixel negative/positive endcap
# # # # binning r / binning phi
# # # pixelEndcapBinningR = acts.ProtoBinning(
# # #     acts.Binning.r, acts.Binning.bound, 40, 175, 2, 1)
# # # pixelEndcapBinningPhi = acts.ProtoBinning(
# # #     acts.Binning.phi, acts.Binning.closed, -math.pi, math.pi, 56, 1)
# # # pixelEndcapPositions = [620., 720., 840., 980., 1120., 1320., 1520.]
# # # pixelEndcapVolumeIds = [16, 18]
# # # pixelEndcapLayerIdOffsets = [4, 2]
# # #
# # # pixelEndcaps = createEndcaps('Pixel::Endcap',
# # #                              [beamPipeRmax, pixelRmax, pixelZmid, detectorZmax],
# # #                              svlMap, pixelEndcapVolumeIds, pixelEndcapLayerIdOffsets,
# # #                              pixelEndcapPositions, 20.,
# # #                              [pixelEndcapBinningR, pixelEndcapBinningPhi])
# # #
# # # # Pixel barrel
# # # pixelBarrelPositions = [34., 70., 116., 172.]
# # # pixelBarrelBinningZ = acts.ProtoBinning(
# # #     acts.Binning.z, acts.Binning.bound, -500, 500, 14, 1)
# # # pixelBarrel = createBarrel('Pixel::Barrel', [beamPipeRmax, pixelRmax, -pixelZmid, pixelZmid],
# # #                            svlMap, 17, pixelBarrelPositions, 10., pixelBarrelBinningZ, [16, 32, 52, 78])
# # #
# # # pixelBarrel.printDimensions()
# # #
# # # # Pxiel overall
# # # pixelBuilder = ContainerStructure('Pixel',
# # #                                   [beamPipeRmax, pixelRmax, detectorZmin, detectorZmax], [pixelEndcaps[0], pixelBarrel, pixelEndcaps[1]], None)
# # # # PST
# # # #
# # # pstRmax = 220.
# # # pst = EmptyVolume('PST', [pixelRmax, pstRmax, -detectorZmax, detectorZmax])
# # #
# # # # SStrip section
# # # #
# # # sstripRmax = 720.
# # # sstripZmid = 1250.
# # # # SStrip negative/positive endcap
# # # sstripEndcapBinningR = acts.ProtoBinning(
# # #     acts.Binning.r, acts.Binning.bound, 230, 710, 3, 1)
# # # sstripEndcapBinningPhi = acts.ProtoBinning(
# # #     acts.Binning.phi, acts.Binning.closed, -math.pi, math.pi, 40, 1)
# # # sstripEndcapPositions = [1300., 1550., 1850., 2200, 2550., 2950.]
# # # sstripEndcapVolumeIds = [23, 25]
# # # sstripEndcapLayerIdOffsets = [2, 2]
# # #
# # # sstripEndcaps = createEndcaps('ShortStrips::Endcap',
# # #                               [pstRmax, sstripRmax, sstripZmid, detectorZmax],
# # #                               svlMap, sstripEndcapVolumeIds, sstripEndcapLayerIdOffsets,
# # #                               sstripEndcapPositions, 50.,
# # #                               [sstripEndcapBinningR, sstripEndcapBinningPhi])
# # #
# # # # SStrip barrel
# # # sstripBarrelPositions = [260., 360., 500., 660.]
# # # sstripBarrelBinningZ = acts.ProtoBinning(
# # #     acts.Binning.z, acts.Binning.bound, -1100, 1100, 21, 1)
# # # sstripBarrel = createBarrel('ShortStrips::Barrel', [pstRmax, sstripRmax, -sstripZmid, sstripZmid],
# # #                             svlMap, 24, sstripBarrelPositions, 20., sstripBarrelBinningZ, [40, 56, 78, 102])
# # #
# # # sstripBarrel.printDimensions()
# # #
# # # sstripBuilder = ContainerStructure('ShortStrips',
# # #                                    [pstRmax, sstripRmax, detectorZmin, detectorZmax], [sstripEndcaps[0], sstripBarrel, sstripEndcaps[1]], None)
# # #
# # #
# # # # LStrip section
# # # lstripRmax = 1100.
# # # lstripZmid = sstripZmid
# # #
# # # # LStrip negative/positive endcap
# # # lstripEndcapBinningR = acts.ProtoBinning(
# # #     acts.Binning.r, acts.Binning.bound, 720, 1020, 2, 1)
# # # lstripEndcapBinningPhi = acts.ProtoBinning(
# # #     acts.Binning.phi, acts.Binning.closed, -math.pi, math.pi, 40, 1)
# # # lstripEndcapPositions = [1300., 1600., 1900., 2250, 2600., 3000.]
# # # lstripEndcapVolumeIds = [28, 30]
# # # lstripEndcapLayerIdOffsets = [2, 2]
# # #
# # # lstripEndcaps = createEndcaps('LongStrips::Endcap',
# # #                               [sstripRmax, lstripRmax, lstripZmid, detectorZmax],
# # #                               svlMap, lstripEndcapVolumeIds, lstripEndcapLayerIdOffsets,
# # #                               lstripEndcapPositions, 50.,
# # #                               [lstripEndcapBinningR, lstripEndcapBinningPhi])
# # #
# # # # SStrip barrel
# # # lstripBarrelPositions = [830, 1030]
# # # lstripBarrelBinningZ = acts.ProtoBinning(
# # #     acts.Binning.z, acts.Binning.bound, -1100, 1100, 21, 1)
# # # lstripBarrel = createBarrel('LongStrips::Barrel', [sstripRmax, lstripRmax, -lstripZmid, lstripZmid],
# # #                             svlMap, 29, lstripBarrelPositions, 30., lstripBarrelBinningZ, [40, 56])
# # #
# # # lstripBarrel.printDimensions()
# # #
# # # lstripBuilder = ContainerStructure('ShortStrips',
# # #                                    [sstripRmax, lstripRmax, detectorZmin, detectorZmax], [lstripEndcaps[0], lstripBarrel, lstripEndcaps[1]], None)
# # #
# # #
# # # detectorBuilder = ContainerStructure('Detector',
# # #                                      [detectorRmin, lstripRmax, detectorZmin, detectorZmax], [beamPipe, pixelBuilder, pst, sstripBuilder, lstripBuilder], None)
# # #
# # #
# # # detectorBuilderConf = acts.DetectorBuilder.Config()
# # # detectorBuilderConf.name = 'ODD'
# # # detectorBuilderConf.builder = detectorBuilder.getBuilder()
# # #
# # # detectorBuilder = acts.DetectorBuilder(
# # #     detectorBuilderConf, 'ODD Detector Builder', acts.logging.VERBOSE)
# # #
# # # detector = detectorBuilder.construct(acts.GeometryContext())
# # #
# # #
# # # acts.examples.drawDetectorSvg(detector,  ['zr', 'xy'])
# # # acts.examples.drawDetectorVolumesSvg(detector, ['grid_zphi', 'grid_xy'])
# # #
# # #

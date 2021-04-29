// This file is part of the Acts project.
//
// Copyright (C) 2017-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/GdmlDetectorConstruction.hpp"

#include <G4GDMLParser.hh>
#include <G4VisAttributes.hh>
#include <G4LogicalVolume.hh>
#include <G4UnitsTable.hh>
#include <G4PVPlacement.hh>

#include "GeoModelKernel/GeoBox.h"
#include "GeoModelKernel/GeoPhysVol.h"
#include "GeoModelKernel/GeoVGeometryPlugin.h"
#include "GeoModelKernel/GeoGeometryPluginLoader.h"
#include "GeoModel2G4/ExtParameterisedVolumeBuilder.h"
#include "GeoMaterial2G4/GeoExtendedMaterial.h"

using namespace ActsExamples;

GdmlDetectorConstruction::GdmlDetectorConstruction(std::string path)
    : G4VUserDetectorConstruction(), m_path(std::move(path)) {}

G4VPhysicalVolume* GdmlDetectorConstruction::Construct() {
//   G4GDMLParser parser;
//   // TODO how to handle errors
//   parser.Read(m_path,false);
//   G4VPhysicalVolume * fWorld = (G4VPhysicalVolume *)parser.GetWorldVolume();
//   fWorld->GetLogicalVolume()->SetVisAttributes(G4VisAttributes::Invisible);
// //   parser.Write("test.gdml", fWorld->GetLogicalVolume());
//   return (G4VPhysicalVolume *)fWorld;
  
  std::cout<< "Bulding the detector from a plugin: "<<m_path<<std::endl;
  GeoGeometryPluginLoader loader;
  GeoVGeometryPlugin *factory=loader.load(m_path);
  if (!factory) {
    std::cout<<"Error!Cannot load geometry from factory. Exiting!"<<std::endl;
    exit(0);
  }

  GeoPhysVol* world = CreateTheWorld(nullptr);
  factory->create(world);

  std::cout << "ReadGeoModel::buildGeoModel() done." << std::endl;
  std::cout << "First step done. GeoModelTree built from the gdml file." << std::endl;
  
  // build the Geant4 geometry and get an hanlde to the world' volume
  ExtParameterisedVolumeBuilder* builder = new ExtParameterisedVolumeBuilder("ATLAS");
  
  std::cout << "Building G4 geometry."<<std::endl;
  G4LogicalVolume* envelope = builder->Build(world);
  
  G4VPhysicalVolume* physWorld= new G4PVPlacement(0,G4ThreeVector(),envelope,envelope->GetName(),0,false,0,false);
  
  G4VPhysicalVolume* fWorld = physWorld;
  fWorld->GetLogicalVolume()->SetVisAttributes(G4VisAttributes::Invisible);
  
  std::cout << "Second step done. Geant4 geometry created from GeoModeltree "<< std::endl;
  return fWorld;  
}

GeoPhysVol*  GdmlDetectorConstruction::CreateTheWorld(GeoPhysVol* world)
{
  if (world == nullptr)
  {
    // Setup the 'World' volume from which everything else will be suspended
    // Get the materials that we shall use.
    // -------------------------------------//
    GeoMaterial* Air=new GeoMaterial("Air", 1.290*SYSTEM_OF_UNITS::mg/SYSTEM_OF_UNITS::cm3);
    GeoElement* Oxigen = new GeoElement("Oxygen",  "O", 8.0, 16.0*SYSTEM_OF_UNITS::g/SYSTEM_OF_UNITS::mole);
    GeoElement* Nitrogen = new GeoElement("Nitrogen", "N", 7., 14.0067*SYSTEM_OF_UNITS::g/SYSTEM_OF_UNITS::mole);
    
    Air->add(Nitrogen, .8);
    Air->add(Oxigen, .2);
    Air->lock();
    const GeoBox* worldBox = new GeoBox(2000*SYSTEM_OF_UNITS::cm, 2000*SYSTEM_OF_UNITS::cm, 4000*SYSTEM_OF_UNITS::cm);
    const GeoLogVol* worldLog = new GeoLogVol("WorldLog", worldBox, Air);
    //const GeoLogVol* worldLog = new GeoLogVol("WorldLog", worldBox, Ether);
    world = new GeoPhysVol(worldLog);
  }
  return world;
}

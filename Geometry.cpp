// ---------------------------------------------------------------- //
//                Creator: Amaia Razquin Lizarraga                  //
//                Contact: amaiafisika@gmail.com                    //
//                Date:    7 September 2023                         //
// ---------------------------------------------------------------- //
/* This module creates a realistic straw detector geometry, 
   a simplified disk geometry with disk detectors instead straws, 
   and a simple square detector geometry for testing */

#include <iostream>

// Definitions & utilities
#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/Utilities/BinUtility.hpp>
#include <Acts/Utilities/BinnedArrayXD.hpp>

// Surface headers
#include <Acts/Surfaces/DiscSurface.hpp>
#include <Acts/Surfaces/StrawSurface.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Acts/Surfaces/RectangleBounds.hpp>
#include <Acts/Surfaces/RadialBounds.hpp>

// Geometry headers
#include <Acts/Geometry/CuboidVolumeBuilder.hpp>
#include <Acts/Geometry/CuboidVolumeBounds.hpp>
#include <Acts/Geometry/CylinderVolumeBounds.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/Geometry/TrackingGeometryBuilder.hpp>
#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/Geometry/Layer.hpp>
#include <Acts/Geometry/DiscLayer.hpp>
#include <Acts/Geometry/LayerArrayCreator.hpp>
#include <Acts/Geometry/LayerCreator.hpp>
#include <Acts/Geometry/SurfaceArrayCreator.hpp>
#include <Acts/Geometry/NavigationLayer.hpp>

// Material headers
#include <Acts/Material/Material.hpp>
#include <Acts/Material/HomogeneousSurfaceMaterial.hpp>
#include <Acts/Material/HomogeneousVolumeMaterial.hpp>

// Visualization
#include <Acts/Visualization/ObjVisualization3D.hpp>
#include <Acts/Visualization/GeometryView3D.hpp>
#include <Acts/Visualization/ViewConfig.hpp>

#include "Geometry.hpp"

using namespace Acts::UnitLiterals;

// Create layour so that it can be use in every necessary module
std::map<std::string, std::vector<std::shared_ptr<const Acts::Surface>>> tracker_layout;

/* Note that the material in this functions is manually inputted from known data.
However, it can be automatically set by using the following code:
    #include <Acts/Plugins/TGeo/TGeoParser.hpp> 
    #include <Acts/Plugins/TGeo/TGeoSurfaceConverter.hpp>

    #include <TGeoManager.h>
    #include <TFile.h>
    #include <TTree.h>

    // Import geometry (based on TGeoParserTest)
	new TGeoManager("Geometry", "Comet geometry");
	TGeoManager *geometry = TGeoManager::Import("comet-II.root");
	
	// Initialize geometry context
	Acts::GeometryContext tgContext = Acts::GeometryContext();

	// If geometry exists, get sensitive nodes
	if (gGeoManager != nullptr) {
		std::string volumeName = "comet";
		
		Acts::TGeoParser::Options tgpOptions;
		tgpOptions.volumeNames = {volumeName};
		// Identify the sensor(s)/target(s) by name
		tgpOptions.targetNames = {"StrawTube_239"};
		
		Acts::TGeoParser::State tgpState;
		tgpState.volume = gGeoManager->GetTopVolume();
		
		// Parse the full ones
    	Acts::TGeoParser::select(tgpState, tgpOptions);
		
		// Material properties
		double radlen, intlen, ar, z, massRho;
		double STradlen, STintlen, STar, STz, STmassRho;
		
		// Crete Acts::Surface of the ECAL
		int i = 0;
		for (auto& snode : tgpState.selectedNodes) {
			i++;
			
			// Material
			const auto& mat = *(snode.node->GetVolume()->GetMaterial());
			
			radlen  = mat.GetRadLen() * Acts::UnitConstants::cm;
			intlen  = mat.GetIntLen() * Acts::UnitConstants::cm;
			ar      = mat.GetA();
			z       = mat.GetZ(); 
			massRho = mat.GetDensity(); 	
        }
    }
*/

// Disk geometry where instead of straw layers, a single disk detector is placed
std::unique_ptr<const Acts::TrackingGeometry> 
    aBuildDiskGeometry(const Acts::GeometryContext& gctx, bool debug){

    std::cout << "ACTS:: Building tracking geometry." << std::endl;
	// Visualization	
	double outputScalor = 1.0;
  	size_t outputPrecision = 6;
  	Acts::ObjVisualization3D objVis(outputPrecision, outputScalor);

	// Material (manually set, see line 47)
	double radlen  = 283.182_mm; // in mm
	double intlen  = 558.747_mm; // in mm
	double ar      = 12.972;
	double za      = 6.50013;
	double massRho = 1.40453; // g/cm3
	if (debug){
		std::cout << "Straw tube material is:\n" <<
					"Radiation length:      " << radlen << "\n"  <<
					"Interaction length:    " << intlen << "\n"  << 
					"Relative atomic mass:  " << ar << "\n"      << 
					"Nuclear charge number: " << za << "\n"      <<
					"Mass density:          " << massRho << "\n" << 
					"-----------------------------------------" << std::endl;
	}	
	Acts::Material almy = Acts::Material::fromMassDensity(radlen, intlen,
				ar, za, massRho); // aluminium + mylar	
	double straw_thickness = 0.00125_cm*2;   // StrawTubeThickness
	Acts::MaterialSlab straw_wall(almy, straw_thickness);
  	
  	// Create StrawTube surface using StrawTrkPhase-II.macro 
    // Should be automatically inputted using the TGeometry plugin
	double radius     = 0.49_cm;      // StrawTubeOuterRadius
	double StationD   = 39.0_cm;      // Distance between stations
	double tubeD      = 1.08_cm;      // Distance between tubes    
	double layerD     = 0.935_cm;     // Distance between layers   
    double separation = 1.0575_cm;    // Separation to build geometry
	double manifoldD  = 1.08_cm;      // Distance between mnifolds 
	double rMin       = 0.0_cm;       // No hole inside disc
	double rMax       = 160_cm / 2;   // Longest straw
	double halfPhiSector(M_PI);       // Full circle 
	int NofStations   = 8;            // Number of stations
    double clearance  = 0.25;         // For the protolayer envelope
	int mani, lay;

    // Vector of volumes to glue
    std::vector<std::shared_ptr<Acts::TrackingVolume>> volumes;

  	// Loop to create straws in each of the 8 stations
	for (int i = 0; i < NofStations; i++) {
		auto station = "Station_" + std::to_string(i);
        // Build one array per station
        // Volume information (station dimension and filling)
        auto name                   =  station;
        Acts::Vector3 tVol_position = {0., 0., -136.5_cm + i * StationD};
        Acts::Vector3 tVol_length   = {rMax * 2, rMax * 2, 2 * separation + layerD *3.0/2}; 
        Acts::Material subdet_mat   = Acts::Material(); // Vacuum 
        auto tVol_material          =
      	    std::make_shared<Acts::HomogeneousVolumeMaterial>(subdet_mat);

        Acts::LayerVector layVec;
        // Binning range for layer array
        std::pair<double, double> minMax = 
            std::make_pair(std::numeric_limits<double>::max(), -std::numeric_limits<double>::max());

        // Build four discs per layer
        int manifold = 0;
        while (manifold < 2) {
            if (manifold == 0) {mani = -1;}
            else if (manifold == 1) {mani = 1;}
		    int layer = 0;
		    while (layer < 2) {
                if (layer == 0) {lay = -1;}
                else if (layer == 1) {lay = 1;}

                // One Acts layer corresponds to one StrawTube layer
                Acts::LayerCreator::Config lCfg;
                lCfg.surfaceArrayCreator = std::make_shared<const Acts::SurfaceArrayCreator>(Acts::getDefaultLogger("SurfaceArrayCreator", 
      		        Acts::Logging::INFO));
                Acts::LayerCreator layerCreator(lCfg);
                layerCreator.setLogger(Acts::getDefaultLogger("LayerCreator", 
      		        Acts::Logging::INFO));
                Acts::RotationMatrix3 rot{{1, 0, 0}, {1, 0, 0}, {0, 0, 1}}; 
                Acts::Translation3 trans{0, 0, -136.5_cm + i * StationD + mani*separation + lay*layerD/2}; //  + 7.1175_cm
                Acts::Transform3 trafo = Acts::Transform3(trans * rot);
                
                // The center is in the middle of the Detector Solenoid (1080, 0, 333.45) cm wrt COMET
                Acts::Translation3 translation{0., 0., -136.5_cm + i * StationD + mani*separation + lay*layerD/2}; 
                auto pTransform = Acts::Transform3(translation);
                
                // Create discs 
                std::vector<std::shared_ptr<const Acts::Surface>> discs;
                if (debug) {
                    std::cout << "Disc in " << station << std::endl;
                    std::cout << "THE SENSOR TRANSFORM - TRANSLATION" << std::endl;
                    std::cout << pTransform.translation()(0)   << 
                        " " << pTransform.translation()(1) << 
                        " " << pTransform.translation()(2) << std::endl;
                    std::cout << "THE SENSOR TRNASFORM - ROTATION"    << std::endl;
                    std::cout << pTransform.rotation() << std::endl;
                }
                
                auto strawDisc = Acts::Surface::makeShared<Acts::DiscSurface>(pTransform, rMin, rMax, halfPhiSector);
                strawDisc->assignSurfaceMaterial(
                        std::make_shared<Acts::HomogeneousSurfaceMaterial>(
                        straw_wall));
                tracker_layout[station].push_back(strawDisc);
                discs.push_back(strawDisc);
                layer++;
                Acts::GeometryView3D::drawSurface(objVis, *strawDisc, gctx);

                // Create proto layer 
                Acts::ProtoLayer pl{gctx, discs};
                pl.envelope[Acts::binX] = {0, 0};
                pl.envelope[Acts::binY] = {0, 0};
                pl.envelope[Acts::binZ] = {straw_thickness * 0.5 + clearance,
      		        straw_thickness * 0.5 + clearance};
                size_t bins = 1;
                std::shared_ptr<const Acts::Layer> layer_txiki = layerCreator.planeLayer(
                    gctx, discs, bins, 
                    bins, Acts::BinningValue::binZ, pl, trafo); 

                /* Alternative proto layer with radial binning
                pl.extent.range(Acts::binZ).set(trans.z() - (straw_thickness * 0.5 + clearance),
      		        trans.z() + straw_thickness * 0.5 + clearance);
                pl.envelope[Acts::binR] = {0., 0.};
                pl.envelope[Acts::binZ] = {straw_thickness * 0.5 + clearance,
      		        straw_thickness * 0.5 + clearance};
                int bins = 1;
                
                std::cout << "Creating disc layer" << std::endl;
                std::shared_ptr<Acts::Layer> layer_txiki = layerCreator.discLayer( //Acts::BinningType::arbitrary
                    gctx, std::move(discs), 3, 
                    15, pl);*/
                
                // Binning range
                double surfacePosMin = pl.min(Acts::binZ);
                double surfacePosMax = pl.max(Acts::binZ);

                // Test if new extreme is found and set it
                if (surfacePosMin < minMax.first) {
                    minMax.first = surfacePosMin;
                }
                if (surfacePosMax > minMax.second) {
                    minMax.second = surfacePosMax;
                }  

                layVec.push_back(layer_txiki);
      		} // layer in manifold
            manifold++;
        } // manifold in station
        // Build layer array (one array per station)
        Acts::Vector3 minVolumeBoundaries = tVol_position - 0.5 * tVol_length;
        Acts::Vector3 maxVolumeBoundaries = tVol_position + 0.5 * tVol_length;

        minMax.first  = std::min(minMax.first, minVolumeBoundaries(Acts::binZ));
        minMax.second = std::max(minMax.second, maxVolumeBoundaries(Acts::binZ));
    
        Acts::LayerArrayCreator::Config lacCnf;
        Acts::LayerArrayCreator layArrCreator(
            lacCnf, Acts::getDefaultLogger("LayerArrayCreator", Acts::Logging::INFO));
        std::unique_ptr<const Acts::LayerArray> layArr(
            layArrCreator.layerArray(gctx, layVec, minMax.first, minMax.second, 
        Acts::BinningType::arbitrary, Acts::BinningValue::binZ));
    
        // Tracking volume (one station)
        std::shared_ptr<Acts::TrackingVolume> trackVolume;
        Acts::Transform3 trafo(Acts::Transform3::Identity());
        trafo.translation() = tVol_position;
        double rmin  = 0.0; 
        auto zbounds = std::make_shared<const Acts::CylinderVolumeBounds>(
            rmin, tVol_length.x() * 0.5, tVol_length.z() * 0.5);
        auto qbounds = std::make_shared<const Acts::CuboidVolumeBounds>(
            tVol_length.x() * 0.5, tVol_length.y() * 0.5, tVol_length.z() * 0.5);
    
        trackVolume  = Acts::TrackingVolume::create(trafo, zbounds, tVol_material, 
            std::move(layArr), nullptr, {}, name);

        volumes.push_back(trackVolume);

    } // stations

    // Loop to create 7 gap volumes in between stations
    for (int i = 0; i < NofStations; i++) {
        using LayerOrderPosition = std::pair<std::shared_ptr<const Acts::Layer>, Acts::Vector3>;
        auto name                   =  "Gap_" + std::to_string(i);
		Acts::Vector3 tVol_position = {0., 0., -136.5_cm + i * StationD - StationD * 0.5}; 
        Acts::Vector3 tVol_length   = {rMax * 2, rMax * 2, (StationD - (2 * separation + layerD * 3./2))}; 
        Acts::Material subdet_mat   = Acts::Material(); // Vacuum 
        auto tVol_material          =
      	    std::make_shared<Acts::HomogeneousVolumeMaterial>(subdet_mat);

        Acts::LayerVector layVec;
        // Binning range for layer array
       	std::pair<double, double> minMax = 
            std::make_pair(std::numeric_limits<double>::max(), -std::numeric_limits<double>::max());

        std::vector<LayerOrderPosition> layerOrderVector;

        Acts::MaterialSlab nav_mat(subdet_mat, straw_thickness);

        // Build one navigation layer inside the volume
        Acts::LayerCreator::Config lCfg;
        lCfg.surfaceArrayCreator = std::make_shared<const Acts::SurfaceArrayCreator>(Acts::getDefaultLogger("SurfaceArrayCreator", 
      		Acts::Logging::INFO));
        Acts::LayerCreator layerCreator(lCfg);
        layerCreator.setLogger(Acts::getDefaultLogger("LayerCreator", 
      		Acts::Logging::INFO));
        Acts::RotationMatrix3 rot{{1, 0, 0}, {1, 0, 0}, {0, 0, 1}}; 
        Acts::Translation3 trans(tVol_position); 

        auto pTransform = Acts::Transform3(trans);
            
        std::vector<std::shared_ptr<const Acts::Surface>> discs;
        if (debug) {
            std::cout << "Disc in " << name << std::endl;
        	std::cout << "THE SENSOR TRANSFORM - TRANSLATION" << std::endl;
            std::cout << pTransform.translation()(0)   << 
            " " << pTransform.translation()(1) << 
            " " << pTransform.translation()(2) << std::endl;
            std::cout << "THE SENSOR TRNASFORM - ROTATION"    << std::endl;
            std::cout << pTransform.rotation() << std::endl;
        }
                
        auto gapSurface = Acts::Surface::makeShared<Acts::DiscSurface>(pTransform, rMin, rMax, halfPhiSector);
        gapSurface->assignSurfaceMaterial(
            std::make_shared<Acts::HomogeneousSurfaceMaterial>(
            nav_mat));
        tracker_layout[name].push_back(gapSurface);
        discs.push_back(gapSurface);
        Acts::GeometryView3D::drawSurface(objVis, *gapSurface, gctx);
        std::shared_ptr<const Acts::Layer> navLayer = Acts::NavigationLayer::create(
            std::move(gapSurface));

        layerOrderVector.push_back(LayerOrderPosition(navLayer, 
            navLayer->binningPosition(gctx, Acts::BinningValue::binZ)));

        // Create proto layer (just to get min max)
        Acts::ProtoLayer pl{gctx, discs};
        pl.envelope[Acts::binX] = {0, 0};
        pl.envelope[Acts::binY] = {0, 0};
        pl.envelope[Acts::binZ] = {straw_thickness / 2. + clearance,
      	    straw_thickness / 2. + clearance};

        // Binning range
        double surfacePosMin = pl.min(Acts::binZ);
        double surfacePosMax = pl.max(Acts::binZ);

        // Test if new extreme is found and set it
        if (surfacePosMin < minMax.first) {
        	minMax.first = surfacePosMin;
        }
        if (surfacePosMax > minMax.second) {
    	    minMax.second = surfacePosMax;
        }  
        
        // Build layer array (one array per station)
        Acts::Vector3 minVolumeBoundaries = tVol_position - 0.5 * tVol_length;
        Acts::Vector3 maxVolumeBoundaries = tVol_position + 0.5 * tVol_length;

        minMax.first  = std::min(minMax.first, minVolumeBoundaries(Acts::binZ));
        minMax.second = std::max(minMax.second, maxVolumeBoundaries(Acts::binZ));

        // Create equidistanct bin utility
        std::unique_ptr<const Acts::BinUtility> binUtility = std::make_unique<const
    		Acts::BinUtility>(layerOrderVector.size(), minMax.first, minMax.second,
    		Acts::open, Acts::BinningValue::binZ);
        std::unique_ptr<const Acts::LayerArray> layArr = std::make_unique<const 
    		Acts::BinnedArrayXD<Acts::LayerPtr>>(layerOrderVector, std::move(binUtility));

        // Tracking volume (one station)
        std::shared_ptr<Acts::TrackingVolume> trackVolume;
        Acts::Transform3 trafo(Acts::Transform3::Identity());
        trafo.translation() = tVol_position;
        double rmin  = 0.0; 
        auto zbounds = std::make_shared<const Acts::CylinderVolumeBounds>(
    		rmin, tVol_length.x() * 0.5, tVol_length.z() * 0.5);
        auto qbounds = std::make_shared<const Acts::CuboidVolumeBounds>(
        	tVol_length.x() * 0.5, tVol_length.y() * 0.5, tVol_length.z() * 0.5);
    
        trackVolume  = Acts::TrackingVolume::create(trafo, zbounds, tVol_material, 
        	std::move(layArr), nullptr, {}, name);

        volumes.push_back(trackVolume);
    } // gap 
    objVis.write("StrECAL_disks");

    // Sort volumes according to the center location
    std::sort(volumes.begin(), volumes.end(),
        [](const Acts::TrackingVolumePtr& lhs, const Acts::TrackingVolumePtr& rhs) {
            return lhs->center().z() < rhs->center().z();
        });
    
    // Print out volume position
    //for (auto& volume : volumes) {
        //std::cout << volume->volumeName() << " volume is at: " << volume->center().z() << std::endl;
    //}

    // Glue volumes
    for (unsigned int i = 0; i < volumes.size() - 1; i++) {
        volumes[i + 1]->glueTrackingVolume(
            gctx, Acts::BoundarySurfaceFace::negativeFaceXY, volumes[i].get(),
            Acts::BoundarySurfaceFace::positiveFaceXY);
        volumes[i]->glueTrackingVolume(
            gctx, Acts::BoundarySurfaceFace::positiveFaceXY, volumes[i + 1].get(),
            Acts::BoundarySurfaceFace::negativeFaceXY);
    }

    Acts::Translation3 world_trans{0., 0., 0.}; 
    Acts::Transform3 world_trafo = Acts::Transform3(world_trans); 
    auto world_qsize = std::make_shared<const Acts::CuboidVolumeBounds>(
        2568_mm * 0.5, 2568_mm * 0.5, 3150_mm * 0.5); // half lengths
    auto world_zsize = std::make_shared<const Acts::CylinderVolumeBounds>(
            0., 2568_mm * 0.5, 3150_mm * 0.5);

    // Vector of confined volumes
    std::vector<std::pair<Acts::TrackingVolumePtr, Acts::Vector3>> tapVec;
    tapVec.reserve(volumes.size());
    for (auto& tVol : volumes) {
        tapVec.push_back(std::make_pair(tVol, tVol->center()));
    }

    // Set bin boundaries along binning
    std::vector<float> binBoundaries;
    binBoundaries.push_back(volumes[0]->center().z() - 
       (StationD - (2 * separation + layerD  * 3./2)) * 0.5);

    for (size_t i = 0; i < volumes.size(); i++) {
        if (i % 2 == 0) {
      		binBoundaries.push_back(volumes[i]->center().z() + 
                (StationD - (2 * separation + layerD  * 3./2)) * 0.5);
            //std::cout << "Volume boundary at: " << volumes[i]->center().z() + 
                //(StationD - (2 * separation + layerD + 2 * straw_thickness)) * 0.5 << std::endl;
      	}
      	else if (i % 2 != 0) {
      		binBoundaries.push_back(volumes[i]->center().z() + 
                (2 * separation + layerD  * 3./2) * 0.5);
            //std::cout << "Volume boundary at: " << volumes[i]->center().z() + 
                //(2 * separation + layerD + 2 * straw_thickness) * 0.5 << std::endl;
      	}
    }

    // Build binning
    Acts::BinningData binData(Acts::BinningOption::open, 
        Acts::BinningValue::binZ, binBoundaries);
    auto bu = std::make_unique<const Acts::BinUtility>(binData);

    // Build Tracking Volume Array
    std::shared_ptr<const Acts::TrackingVolumeArray> trVolArr(
        new Acts::BinnedArrayXD<Acts::TrackingVolumePtr>(tapVec, std::move(bu)));
    
    // Create world volume
    Acts::MutableTrackingVolumePtr mtvp(Acts::TrackingVolume::create(
        world_trafo, world_zsize, trVolArr, "World"));
    // I'm not sure what is the purpose of this hook
    std::shared_ptr<const Acts::GeometryIdentifierHook> geometryIdentifierHook =
        std::make_shared<Acts::GeometryIdentifierHook>();

    return std::make_unique<Acts::TrackingGeometry>(
        mtvp, nullptr, *geometryIdentifierHook);
}

// Straw geometry with the disk accurately created
std::unique_ptr<const Acts::TrackingGeometry> 
    aBuildStrawGeometry(const Acts::GeometryContext& gctx, bool debug){

    std::cout << "ACTS:: Building tracking geometry." << std::endl;
	// Visualization	
	double outputScalor = 1.0;
  	size_t outputPrecision = 6;
  	Acts::ObjVisualization3D objVis(outputPrecision, outputScalor);

	// Material (manually set, see line 47)
	double radlen  = 283.182_mm; // in mm
	double intlen  = 558.747_mm; // in mm
	double ar      = 12.972;
	double za      = 6.50013;
	double massRho = 1.40453;    // g/cm3
	if (debug){
		std::cout << "Straw tube material is:\n" <<
					"Radiation length:      " << radlen << "\n"  <<
					"Interaction length:    " << intlen << "\n"  << 
					"Relative atomic mass:  " << ar << "\n"      << 
					"Nuclear charge number: " << za << "\n"      <<
					"Mass density:          " << massRho << "\n" << 
					"-----------------------------------------" << std::endl;
	}	
	Acts::Material almy = Acts::Material::fromMassDensity(radlen, intlen,
				ar, za, massRho); // aluminium + mylar	
	double straw_thickness = 0.00125_cm;   // StrawTubeThickness
	Acts::MaterialSlab straw_wall(almy, straw_thickness);
  	
  	// Create StrawTube surface using StrawTrkPhase-II.macro 
	double radius     = 0.49_cm;      // StrawTubeOuterRadius
	double StationD   = 39.0_cm;      // Distance between stations
	double tubeD      = 1.08_cm;      // Distance between tubes    
	double layerD     = 0.935_cm;     // Distance between layers   
    double separation = 1.0575_cm;    // Separation to build geometry
	double manifoldD  = 1.08_cm;      // Distance between mnifolds 
	double rMin       = 0.0_cm;       // No hole inside disc
	double rMax       = 160_cm / 2;   // Longest straw
	double halfPhiSector(M_PI);       // Full circle 
	int NofStations   = 8;            // Number of stations
    int NofGroups     = 15;           // Number of groups in one layer
    int NofStrinGrup  = 8;            // Number of straws in group
    double clearance  = 0.03;         // For the protolayer envelope
    std::vector<double> length{790, 990, 1130, 1230, 1310, 1370, 1370, 1390, 1370,
			1370, 1310, 1230, 1130, 990, 790}; // length of the straws
	int mani, lay;
    
    std::cout << "ACTS:: Building layout for Straw Tracker." << std::endl;

    // Vector of volumes to glue
    std::vector<std::shared_ptr<Acts::TrackingVolume>> volumes;

  	// Loop to create straws in each of the 8 stations
	for (int i = 0; i < NofStations; i++) {
		auto station = "Station_" + std::to_string(i);
        Acts::Material subdet_mat   = Acts::Material(); // Vacuum 
        auto tVol_material          =
      	    std::make_shared<Acts::HomogeneousVolumeMaterial>(subdet_mat);

        // Build four discs per layer
        int manifold = 0;
        int layernumber = 2;
        while (manifold < 2) {
            if (manifold == 0) {mani = -1;}
            else if (manifold == 1) {mani = 1;}

            // Build two arrays per station, two volumes closely tight to the layers
            // Volume information (manifold dimension and filling)
            auto name                   =  station + "_" + std::to_string(manifold);
            Acts::Vector3 tVol_position = {0., 0., -136.5_cm + i * StationD + mani * separation - 2. * 0.5};
            Acts::Vector3 tVol_length   = {rMax * 2, rMax * 2, layerD + radius * 2 + 2.}; 
            Acts::LayerVector layVec;
            std::cout << tVol_position.z() - tVol_length.z() * 0.5 << " " << tVol_position.z() << " " << tVol_position.z() + tVol_length.z() * 0.5 << std::endl;
            // Binning range for layer array
            std::pair<double, double> minMax = 
                std::make_pair(std::numeric_limits<double>::max(), -std::numeric_limits<double>::max());

		    int layer = 0;
		    while (layer < 2) {
                if (layer == 0) {lay = -1;}
                else if (layer == 1) {lay = 1;}

                auto station_layer = station + "_" + std::to_string(layernumber);

                // One Acts layer corresponds to one StrawTube layer
                Acts::LayerCreator::Config lCfg;
                lCfg.surfaceArrayCreator = std::make_shared<const Acts::SurfaceArrayCreator>(Acts::getDefaultLogger("SurfaceArrayCreator", 
      		        Acts::Logging::INFO));
                Acts::LayerCreator layerCreator(lCfg);
                layerCreator.setLogger(Acts::getDefaultLogger("LayerCreator", 
      		        Acts::Logging::INFO));
                Acts::RotationMatrix3 rot{{1, 0, 0}, {1, 0, 0}, {0, 0, 1}}; 
                Acts::Translation3 trans{0, 0, -136.5_cm + i * StationD + mani*separation + lay*layerD/2}; //  + 7.1175_cm
                Acts::Transform3 trafo = Acts::Transform3(trans * rot);
                
                // Loop over the groups and the straws within each group
                std::vector<std::shared_ptr<const Acts::Surface>> straws;
                int count = 0;
                for (int j = 0; j < NofGroups; j++) {
                    double halfZ = length[j] * 0.5;
                    for (int k = 0; k < NofStrinGrup; k++) {
                        Acts::Transform3 pTransform;
                        if (mani == -1) {
                            Acts::Translation3 translation{
                                -64.53_cm + (count + k) * tubeD + tubeD * 0.5 * layer, 
                                0.0, // + radius*layer, 
                                -136.5_cm + i * StationD + mani*separation + lay*layerD/2};
                            Acts::RotationMatrix3 rot{{0, -1, 0}, {1, 0, 0}, {0, 0, 1}};
                            pTransform = Acts::Transform3(translation*rot);
                            pTransform.rotate(Eigen::AngleAxisd(1.5*M_PI, 
                                Acts::Vector3(0, 1, 0)));
                            //pTransform.rotate(Eigen::AngleAxisd(0.5*M_PI, Acts::Vector3(0, 0, 1)));
                        }
                        else if (mani == 1) {
                            Acts::Translation3 translation{
                                0.0, 
                                -64.53_cm + (count + k) * tubeD + tubeD * 0.5 * layer, // + radius*layer, 
                                -136.5_cm + i * StationD + mani*separation + lay*layerD/2};
                            pTransform = Acts::Transform3(translation);
                            pTransform.rotate(Eigen::AngleAxisd(1.5*M_PI, 
						        Acts::Vector3(0, 1, 0)));
                            //pTransform.rotate(Eigen::AngleAxisd(0.5*M_PI, Acts::Vector3(0, 0, 1)));
                        }
                        else {std::cout << "There's been a problem, check the manifolds\n";}

                        if (debug) {
                            std::cout << "Disc in " << station_layer << std::endl;
                            std::cout << "THE SENSOR TRANSFORM - TRANSLATION" << std::endl;
                            std::cout << pTransform.translation()(0)   << 
                                " " << pTransform.translation()(1) << 
                                " " << pTransform.translation()(2) << std::endl;
                            std::cout << "THE SENSOR TRNASFORM - ROTATION"    << std::endl;
                            std::cout << pTransform.rotation() << std::endl;
                        }

                        auto straw = Acts::Surface::makeShared<Acts::StrawSurface>(pTransform, radius, halfZ);
                        straw->assignSurfaceMaterial(
                            std::make_shared<Acts::HomogeneousSurfaceMaterial>(straw_wall));
                        straws.push_back(straw);
                        tracker_layout[station_layer].push_back(straw);
                        Acts::GeometryView3D::drawSurface(objVis, *straw, gctx);
                    } // Straws in group
                    count = count + NofStrinGrup;
                } // groups in layer

                // Create proto layer with binning depending on straw orientation
                size_t binsX, binsY;
                double xi, xf, yi, yf;
                double add = 44.8 * 0.5;
                if (layernumber == 2 || layernumber == 4) {
                    binsX = 1; binsY = 1;
                    if (layernumber == 2) {
                        xi    = 0 + add;   xf    = tubeD * 0.5 + add;
                    }
                    else if (layernumber == 4) {
                        xi    = tubeD * 0.5 + add;   xf    = 0 + add;
                    }
                    yi    = 0;   yf    = 0;
                }
                else if (layernumber == 6 || layernumber == 8) {
                    binsX = 1; binsY = 1;
                    xi    = 0;   xf    = 0;
                    if (layernumber == 6) {
                        yi    = 0 + add;   yf    = tubeD * 0.5 + add;
                    }
                    else if (layernumber == 8) {
                        yi    = tubeD * 0.5 + add;   yf    = 0 + add;
                    }
                }
                
                Acts::ProtoLayer pl{gctx, straws};
                pl.envelope[Acts::binX] = {xi, xf};
                pl.envelope[Acts::binY] = {yi, yf};
                pl.envelope[Acts::binZ] = {layerD * 0.5 - clearance,
      		        layerD * 0.5 - clearance};
                
                std::shared_ptr<const Acts::Layer> layer_txiki = layerCreator.planeLayer(
                    gctx, straws, binsX, binsY, 
                    Acts::BinningValue::binZ, pl, trafo);
                
                // Binning range
                double surfacePosMin = pl.min(Acts::binZ);
                double surfacePosMax = pl.max(Acts::binZ);

                // Test if new extreme is found and set it
                if (surfacePosMin < minMax.first) {
                    minMax.first = surfacePosMin;
                }
                if (surfacePosMax > minMax.second) {
                    minMax.second = surfacePosMax;
                }  

                layVec.push_back(layer_txiki);
                layer++;
                layernumber = layernumber + 2;
      		} // layer in manifold

            // -------------- Build volume -----------------------------------------------
            // Build layer array (two arrays per station)
            Acts::Vector3 minVolumeBoundaries = tVol_position - 0.5 * tVol_length;
            Acts::Vector3 maxVolumeBoundaries = tVol_position + 0.5 * tVol_length;

            minMax.first  = std::min(minMax.first, minVolumeBoundaries(Acts::binZ));
            minMax.second = std::max(minMax.second, maxVolumeBoundaries(Acts::binZ));

            Acts::LayerArrayCreator::Config lacCnf;
            Acts::LayerArrayCreator layArrCreator(
                lacCnf, Acts::getDefaultLogger("LayerArrayCreator", Acts::Logging::INFO));
            std::unique_ptr<const Acts::LayerArray> layArr(
                layArrCreator.layerArray(gctx, layVec, minMax.first, minMax.second, 
            Acts::BinningType::equidistant, Acts::BinningValue::binZ));
            
            // Tracking volume (one station)
            std::shared_ptr<Acts::TrackingVolume> trackVolume;
            Acts::Transform3 vol_trafo(Acts::Transform3::Identity());
            vol_trafo.translation() = tVol_position;
            double rmin  = 0.0; 
            auto zbounds = std::make_shared<const Acts::CylinderVolumeBounds>(
                rmin, tVol_length.x() * 0.5, tVol_length.z() * 0.5);
            auto qbounds = std::make_shared<const Acts::CuboidVolumeBounds>(
                tVol_length.x() * 0.5, tVol_length.y() * 0.5, tVol_length.z() * 0.5);
            
            trackVolume  = Acts::TrackingVolume::create(vol_trafo, zbounds, tVol_material, 
                std::move(layArr), nullptr, {}, name);

            volumes.push_back(trackVolume);

            manifold++;
        } // manifold in station
    } // stations

    // Loop to create 7 gap volumes in between stations
    for (int i = 0; i < NofStations; i++) {
        using LayerOrderPosition = std::pair<std::shared_ptr<const Acts::Layer>, Acts::Vector3>;
        	auto name                   =  "Gap_" + std::to_string(i);
		Acts::Vector3 tVol_position = {0., 0., -136.5_cm + i * StationD - StationD * 0.5 - 2. * 0.5}; 
        Acts::Vector3 tVol_length   = {rMax * 2, rMax * 2, (StationD - (2 * separation + layerD + 2 * radius + 2.))}; 
        Acts::Material subdet_mat   = Acts::Material(); // Vacuum 
        auto tVol_material          =
      		std::make_shared<Acts::HomogeneousVolumeMaterial>(subdet_mat);

        Acts::LayerVector layVec;
        // Binning range for layer array
       	std::pair<double, double> minMax = 
            std::make_pair(std::numeric_limits<double>::max(), -std::numeric_limits<double>::max());

        std::vector<LayerOrderPosition> layerOrderVector;

        Acts::MaterialSlab nav_mat(subdet_mat, straw_thickness);

        // Build one navigation layer inside the volume
        Acts::LayerCreator::Config lCfg;
        lCfg.surfaceArrayCreator = std::make_shared<const Acts::SurfaceArrayCreator>(Acts::getDefaultLogger("SurfaceArrayCreator", 
      		Acts::Logging::INFO));
        Acts::LayerCreator layerCreator(lCfg);
        layerCreator.setLogger(Acts::getDefaultLogger("LayerCreator", 
      	    Acts::Logging::INFO));
        Acts::RotationMatrix3 rot{{1, 0, 0}, {1, 0, 0}, {0, 0, 1}}; 
        Acts::Translation3 trans(tVol_position);

        auto pTransform = Acts::Transform3(trans);
            
        std::vector<std::shared_ptr<const Acts::Surface>> discs;
        if (debug) {
            std::cout << "Disc in " << name << std::endl;
            	std::cout << "THE SENSOR TRANSFORM - TRANSLATION" << std::endl;
            std::cout << pTransform.translation()(0)   << 
                " " << pTransform.translation()(1) << 
               	" " << pTransform.translation()(2) << std::endl;
            	std::cout << "THE SENSOR TRNASFORM - ROTATION"    << std::endl;
            	std::cout << pTransform.rotation() << std::endl;
        }
                
        auto gapSurface = Acts::Surface::makeShared<Acts::DiscSurface>(pTransform, rMin, rMax, halfPhiSector);
        gapSurface->assignSurfaceMaterial(
            std::make_shared<Acts::HomogeneousSurfaceMaterial>(
                nav_mat));
        tracker_layout[name].push_back(gapSurface);
        discs.push_back(gapSurface);
        Acts::GeometryView3D::drawSurface(objVis, *gapSurface, gctx);
        std::shared_ptr<const Acts::Layer> navLayer = Acts::NavigationLayer::create(
            std::move(gapSurface));

        layerOrderVector.push_back(LayerOrderPosition(navLayer, 
            navLayer->binningPosition(gctx, Acts::BinningValue::binZ)));

        // Create proto layer (just to get min max)
        Acts::ProtoLayer pl{gctx, discs};
        pl.envelope[Acts::binX] = {0, 0};
        pl.envelope[Acts::binY] = {0, 0};
        pl.envelope[Acts::binZ] = {straw_thickness / 2. + clearance,
      	    straw_thickness / 2. + clearance};

        // Binning range
        double surfacePosMin = pl.min(Acts::binZ);
        double surfacePosMax = pl.max(Acts::binZ);

        // Test if new extreme is found and set it
        if (surfacePosMin < minMax.first) {
            minMax.first = surfacePosMin;
        }
        if (surfacePosMax > minMax.second) {
            minMax.second = surfacePosMax;
        }  
        
        // Build layer array (one array per station)
        Acts::Vector3 minVolumeBoundaries = tVol_position - 0.5 * tVol_length;
        Acts::Vector3 maxVolumeBoundaries = tVol_position + 0.5 * tVol_length;

        minMax.first  = std::min(minMax.first, minVolumeBoundaries(Acts::binZ));
        minMax.second = std::max(minMax.second, maxVolumeBoundaries(Acts::binZ));

        // Create equidistanct bin utility
        std::unique_ptr<const Acts::BinUtility> binUtility = std::make_unique<const
            Acts::BinUtility>(layerOrderVector.size(), minMax.first, minMax.second,
            	Acts::open, Acts::BinningValue::binZ);
        std::unique_ptr<const Acts::LayerArray> layArr = std::make_unique<const 
            Acts::BinnedArrayXD<Acts::LayerPtr>>(layerOrderVector, std::move(binUtility));

        // Tracking volume (one station)
        std::shared_ptr<Acts::TrackingVolume> trackVolume;
        Acts::Transform3 trafo(Acts::Transform3::Identity());
        trafo.translation() = tVol_position;
        double rmin  = 0.0; 
        auto zbounds = std::make_shared<const Acts::CylinderVolumeBounds>(
            rmin, tVol_length.x() * 0.5, tVol_length.z() * 0.5);
        auto qbounds = std::make_shared<const Acts::CuboidVolumeBounds>(
            tVol_length.x() * 0.5, tVol_length.y() * 0.5, tVol_length.z() * 0.5);
    
        trackVolume  = Acts::TrackingVolume::create(trafo, zbounds, tVol_material, 
            std::move(layArr), nullptr, {}, name);

        volumes.push_back(trackVolume);
    } // gap 

    objVis.write("StrECAL_straws");

    std::cout << "AMOUNT OF VOLUMES: " << volumes.size() << std::endl;

    // Sort volumes according to the center location
    std::sort(volumes.begin(), volumes.end(),
        [](const Acts::TrackingVolumePtr& lhs, const Acts::TrackingVolumePtr& rhs) {
            return lhs->center().z() < rhs->center().z();
        });
    
    // Print out volume position
    for (auto& volume : volumes) {
        std::cout << volume->volumeName() << " volume is at: " << volume->center().z() << std::endl;
    }

    // Glue volumes
    for (unsigned int i = 0; i < volumes.size() - 1; i++) {
        volumes[i + 1]->glueTrackingVolume(
            gctx, Acts::BoundarySurfaceFace::negativeFaceXY, volumes[i].get(),
            Acts::BoundarySurfaceFace::positiveFaceXY);
        volumes[i]->glueTrackingVolume(
            gctx, Acts::BoundarySurfaceFace::positiveFaceXY, volumes[i + 1].get(),
            Acts::BoundarySurfaceFace::negativeFaceXY);
    }

    Acts::Translation3 world_trans{0., 0., 0.}; 
    Acts::Transform3 world_trafo = Acts::Transform3(world_trans); 
    auto world_qsize = std::make_shared<const Acts::CuboidVolumeBounds>(
        2568_mm * 0.5, 2568_mm * 0.5, 3150_mm * 0.5); // half lengths
    auto world_zsize = std::make_shared<const Acts::CylinderVolumeBounds>(
            0., 2568_mm * 0.5, 3150_mm * 0.5);

    // Vector of confined volumes
    std::vector<std::pair<Acts::TrackingVolumePtr, Acts::Vector3>> tapVec;
    tapVec.reserve(volumes.size());
    for (auto& tVol : volumes) {
        tapVec.push_back(std::make_pair(tVol, tVol->center()));
    }

    // Set bin boundaries along binning
    std::set<int> manifoldX {1, 4, 7, 10, 13, 16, 19, 22};
    std::set<int> manifoldY {2, 5, 8, 11, 14, 17, 20, 23};
    std::vector<float> binBoundaries;
    binBoundaries.push_back(volumes[0]->center().z() - 
       (StationD - (2 * separation + layerD + 2 * radius +  2.)) * 0.5);
    for (size_t i = 0; i < volumes.size(); i++) {
        if (i % 4 == 0) {
            binBoundaries.push_back(volumes[i]->center().z() + 
                (StationD - (2 * separation + layerD + 2 * radius +  2.)) * 0.5);
            std::cout << "Gap volume " << i << " found at " << volumes[i]->center().z() << ", boundary set at " << volumes[i]->center().z() + 
                (StationD - (2 * separation + layerD + 2 * radius +  2.)) * 0.5 << std::endl;
        }
        else if (manifoldX.find(i) != manifoldX.end()) {
            binBoundaries.push_back(volumes[i]->center().z() + 
                (layerD + 2 * radius + 2.) * 0.5);
            std::cout << "ManifoldX volume " << i << " found " << volumes[i]->center().z() << ", boundary set at " << volumes[i]->center().z() + 
               (layerD + 2 * radius + 2.) * 0.5 << std::endl;
        }
        else if (manifoldY.find(i) != manifoldY.end()) {
            binBoundaries.push_back(volumes[i]->center().z() + 
                (layerD + 2 * radius + 2.) * 0.5);
            std::cout << "ManifoldY volume " << i << " found " << volumes[i]->center().z() << ", boundary set at " << volumes[i]->center().z() + 
               (layerD + 2 * radius + 2.) * 0.5 << std::endl;
        }
    }

    // Build binning
    Acts::BinningData binData(Acts::BinningOption::open, 
        Acts::BinningValue::binZ, binBoundaries);
    auto bu = std::make_unique<const Acts::BinUtility>(binData);

    // Build Tracking Volume Array
    std::shared_ptr<const Acts::TrackingVolumeArray> trVolArr(
        new Acts::BinnedArrayXD<Acts::TrackingVolumePtr>(tapVec, std::move(bu)));
    
    // Create world volume
    Acts::MutableTrackingVolumePtr mtvp(Acts::TrackingVolume::create(
        world_trafo, world_zsize, trVolArr, "World"));
    // I'm not sure what is the purpose of this hook
    std::shared_ptr<const Acts::GeometryIdentifierHook> geometryIdentifierHook =
        std::make_shared<Acts::GeometryIdentifierHook>();

    return std::make_unique<Acts::TrackingGeometry>(
        mtvp, nullptr, *geometryIdentifierHook);
}

// Simple test geometry with square detectors
std::unique_ptr<const Acts::TrackingGeometry> 
    aBuildSimpleGeometry(const Acts::GeometryContext& gctx, bool debug){

    std::cout << "ACTS:: Building tracking geometry." << std::endl;
	// Visualization	
	double outputScalor = 1.0;
  	size_t outputPrecision = 6;
  	Acts::ObjVisualization3D objVis(outputPrecision, outputScalor);

	// Material
	double radlen  = 283.182_mm; // in mm
	double intlen  = 558.747_mm; // in mm
	double ar      = 12.972;
	double za      = 6.50013;
	double massRho = 1.40453; // g/cm3
	if (debug){
		std::cout << "Straw tube material is:\n" <<
					"Radiation length:      " << radlen << "\n"  <<
					"Interaction length:    " << intlen << "\n"  << 
					"Relative atomic mass:  " << ar << "\n"      << 
					"Nuclear charge number: " << za << "\n"      <<
					"Mass density:          " << massRho << "\n" << 
					"-----------------------------------------" << std::endl;
	}	
	Acts::Material almy = Acts::Material::fromMassDensity(radlen, intlen,
				ar, za, massRho); // aluminium + mylar	
	double straw_thickness = 0.00125_cm*2;   // StrawTubeThickness
	Acts::MaterialSlab straw_wall(almy, straw_thickness);
  	
    // Dimensions
    double squareL  = 60_cm;
    double StationD = 40_cm;
    int NofStations = 4;
    
    std::cout << "ACTS:: Building layout for Straw Tracker." << std::endl;

    // Vector of volumes to glue
    std::vector<std::shared_ptr<Acts::TrackingVolume>> volumes;

  	// Loop to create straws in each of the 8 stations
	for (int i = 0; i < NofStations; i++) {
		auto station = "Station_" + std::to_string(i);
        // Build one array per station
        // Volume INFOrmation (station dimension and filling)
        auto name                   =  station;
        Acts::Vector3 tVol_position = {0., 0., -60_cm + i * StationD};
        Acts::Vector3 tVol_length   = {squareL, squareL, 10_cm}; 
        Acts::Material subdet_mat   = Acts::Material(); // Vacuum 
        auto tVol_material          =
      	    std::make_shared<Acts::HomogeneousVolumeMaterial>(subdet_mat);

        Acts::LayerVector layVec;
        // Binning range for layer array
        std::pair<double, double> minMax = 
            std::make_pair(std::numeric_limits<double>::max(), -std::numeric_limits<double>::max());

        // Layer
        Acts::LayerCreator::Config lCfg;
        lCfg.surfaceArrayCreator = std::make_shared<const Acts::SurfaceArrayCreator>(Acts::getDefaultLogger("SurfaceArrayCreator", 
    	    Acts::Logging::INFO));
        Acts::LayerCreator layerCreator(lCfg);
        layerCreator.setLogger(Acts::getDefaultLogger("LayerCreator", 
      		Acts::Logging::INFO));
        Acts::RotationMatrix3 rot{{1, 0, 0}, {1, 0, 0}, {0, 0, 1}}; 
        Acts::Translation3 trans{0, 0, -60_cm + i * StationD}; 
        Acts::Transform3 trafo = Acts::Transform3(trans * rot);
                
        // Surface
        Acts::Translation3 translation{0, 0, -60_cm + i * StationD}; 
        auto pTransform = Acts::Transform3(translation);
                
        std::vector<std::shared_ptr<const Acts::Surface>> surfaces;
        if (debug) {
            std::cout << "Disc in " << station << std::endl;
            std::cout << "THE SENSOR TRANSFORM - TRANSLATION" << std::endl;
            std::cout << pTransform.translation()(0)   << 
                " " << pTransform.translation()(1) << 
                " " << pTransform.translation()(2) << std::endl;
            std::cout << "THE SENSOR TRNASFORM - ROTATION"    << std::endl;
            std::cout << pTransform.rotation() << std::endl;
        }
        
        std::shared_ptr<const Acts::RectangleBounds> bounds = std::make_shared<const Acts::RectangleBounds>(Acts::RectangleBounds(60_cm*0.5, 60_cm*0.5));
        auto square = Acts::Surface::makeShared<Acts::PlaneSurface>(pTransform, bounds);
        square->assignSurfaceMaterial(
            std::make_shared<Acts::HomogeneousSurfaceMaterial>(straw_wall));
        surfaces.push_back(square); // this should be inside another for loop for straws
    
        Acts::GeometryView3D::drawSurface(objVis, *square, gctx);

        // Create proto layer (in the future it will contain all the straws in a layer)
        Acts::ProtoLayer pl{gctx, surfaces};
        pl.envelope[Acts::binX] = {0., 0.};
        pl.envelope[Acts::binY] = {0., 0.};
        pl.envelope[Acts::binZ] = {1, 1};
        size_t bins = 120;
        std::shared_ptr<const Acts::Layer> layer_txiki = layerCreator.planeLayer(
            gctx, surfaces, Acts::arbitrary, Acts::arbitrary, Acts::BinningValue::binY, pl, trafo);

        // Binning range
        double surfacePosMin = pl.min(Acts::binZ);
        double surfacePosMax = pl.max(Acts::binZ);

        // Test if new extreme is found and set it
        if (surfacePosMin < minMax.first) {
            minMax.first = surfacePosMin;
        }
        if (surfacePosMax > minMax.second) {
            minMax.second = surfacePosMax;
        }  

        layVec.push_back(layer_txiki);
        
        // Build layer array (one array per station)
        Acts::Vector3 minVolumeBoundaries = tVol_position - 0.5 * tVol_length;
        Acts::Vector3 maxVolumeBoundaries = tVol_position + 0.5 * tVol_length;

        minMax.first  = std::min(minMax.first, minVolumeBoundaries(Acts::binZ));
        minMax.second = std::max(minMax.second, maxVolumeBoundaries(Acts::binZ));
    
        Acts::LayerArrayCreator::Config lacCnf;
        Acts::LayerArrayCreator layArrCreator(
            lacCnf, Acts::getDefaultLogger("LayerArrayCreator", Acts::Logging::INFO));
        std::unique_ptr<const Acts::LayerArray> layArr(
            layArrCreator.layerArray(gctx, layVec, minMax.first, minMax.second, 
        Acts::BinningType::arbitrary, Acts::BinningValue::binZ));
    
        // Tracking volume (one station)
        std::shared_ptr<Acts::TrackingVolume> trackVolume;
        Acts::Transform3 trafoo(Acts::Transform3::Identity());
        trafoo.translation() = tVol_position;
        double rmin  = 0.0; 
        auto zbounds = std::make_shared<const Acts::CylinderVolumeBounds>(
            rmin, tVol_length.x() * 0.5, tVol_length.z() * 0.5);
        auto qbounds = std::make_shared<const Acts::CuboidVolumeBounds>(
            tVol_length.x() * 0.5, tVol_length.y() * 0.5, tVol_length.z() * 0.5);
    
        trackVolume  = Acts::TrackingVolume::create(trafoo, qbounds, tVol_material, 
            std::move(layArr), nullptr, {}, name);

        volumes.push_back(trackVolume);

    } // stations

    // Loop to create 7 gap volumes in between stations
    for (int i = 0; i < NofStations; i++) {
        using LayerOrderPosition = std::pair<std::shared_ptr<const Acts::Layer>, Acts::Vector3>;
        auto name                   =  "Gap_" + std::to_string(i);
		Acts::Vector3 tVol_position = {0., 0., -60.0_cm - StationD * 0.5 + StationD * i}; 
        Acts::Vector3 tVol_length   = {squareL, squareL, (StationD - (2 * 10_cm) * 0.5)}; 
        Acts::Material subdet_mat   = Acts::Material(); // Vacuum 
        auto tVol_material          =
      	    std::make_shared<Acts::HomogeneousVolumeMaterial>(subdet_mat);

        Acts::LayerVector layVec;
        // Binning range for layer array
       	std::pair<double, double> minMax = 
            std::make_pair(std::numeric_limits<double>::max(), -std::numeric_limits<double>::max());

        std::vector<LayerOrderPosition> layerOrderVector;

        Acts::MaterialSlab nav_mat(subdet_mat, straw_thickness);

        // Build one navigation layer inside the volume
        Acts::LayerCreator::Config lCfg;
        lCfg.surfaceArrayCreator = std::make_shared<const Acts::SurfaceArrayCreator>(Acts::getDefaultLogger("SurfaceArrayCreator", 
      		Acts::Logging::INFO));
        Acts::LayerCreator layerCreator(lCfg);
        layerCreator.setLogger(Acts::getDefaultLogger("LayerCreator", 
      		Acts::Logging::INFO));
        Acts::RotationMatrix3 rot{{1, 0, 0}, {1, 0, 0}, {0, 0, 1}}; 
        Acts::Translation3 trans(tVol_position); 

        auto pTransform = Acts::Transform3(trans);
                
        // There should be a for loop here for the individual straws
        std::vector<std::shared_ptr<const Acts::Surface>> discs;
        if (debug) {
            std::cout << "Disc in " << name << std::endl;
        	std::cout << "THE SENSOR TRANSFORM - TRANSLATION" << std::endl;
            std::cout << pTransform.translation()(0)   << 
            " " << pTransform.translation()(1) << 
            " " << pTransform.translation()(2) << std::endl;
            std::cout << "THE SENSOR TRNASFORM - ROTATION"    << std::endl;
            std::cout << pTransform.rotation() << std::endl;
        }

        std::shared_ptr<const Acts::RectangleBounds> bounds = std::make_shared<const Acts::RectangleBounds>(Acts::RectangleBounds(60_cm*0.5, 60_cm*0.5));
        auto gapSurface =  Acts::Surface::makeShared<Acts::PlaneSurface>(pTransform, bounds);
        gapSurface->assignSurfaceMaterial(
            std::make_shared<Acts::HomogeneousSurfaceMaterial>(
            nav_mat));
        discs.push_back(gapSurface);
        Acts::GeometryView3D::drawSurface(objVis, *gapSurface, gctx);
        std::shared_ptr<const Acts::Layer> navLayer = Acts::NavigationLayer::create(
            std::move(gapSurface));

        layerOrderVector.push_back(LayerOrderPosition(navLayer, 
            navLayer->binningPosition(gctx, Acts::BinningValue::binZ)));

        // Create proto layer (just to get min max)
        Acts::ProtoLayer pl{gctx, discs};
        pl.envelope[Acts::binX] = {0, 0};
        pl.envelope[Acts::binY] = {0, 0};
        pl.envelope[Acts::binZ] = {1, 1};

        // Binning range
        double surfacePosMin = pl.min(Acts::binZ);
        double surfacePosMax = pl.max(Acts::binZ);

        // Test if new extreme is found and set it
        if (surfacePosMin < minMax.first) {
        	minMax.first = surfacePosMin;
        }
        if (surfacePosMax > minMax.second) {
    	    minMax.second = surfacePosMax;
        }  
        
        // Build layer array (one array per station)
        Acts::Vector3 minVolumeBoundaries = tVol_position - 0.5 * tVol_length;
        Acts::Vector3 maxVolumeBoundaries = tVol_position + 0.5 * tVol_length;

        minMax.first  = std::min(minMax.first, minVolumeBoundaries(Acts::binZ));
        minMax.second = std::max(minMax.second, maxVolumeBoundaries(Acts::binZ));

        // Create equidistanct bin utility
        std::unique_ptr<const Acts::BinUtility> binUtility = std::make_unique<const
    		Acts::BinUtility>(layerOrderVector.size(), minMax.first, minMax.second,
    		Acts::open, Acts::BinningValue::binZ);
        std::unique_ptr<const Acts::LayerArray> layArr = std::make_unique<const 
    		Acts::BinnedArrayXD<Acts::LayerPtr>>(layerOrderVector, std::move(binUtility));

        // Tracking volume (one station)
        std::shared_ptr<Acts::TrackingVolume> trackVolume;
        Acts::Transform3 trafo(Acts::Transform3::Identity());
        trafo.translation() = tVol_position;
        double rmin  = 0.0; 
        auto zbounds = std::make_shared<const Acts::CylinderVolumeBounds>(
    		rmin, tVol_length.x() * 0.5, tVol_length.z() * 0.5);
        auto qbounds = std::make_shared<const Acts::CuboidVolumeBounds>(
        	tVol_length.x() * 0.5, tVol_length.y() * 0.5, tVol_length.z() * 0.5);
    
        trackVolume  = Acts::TrackingVolume::create(trafo, qbounds, tVol_material, 
        	std::move(layArr), nullptr, {}, name);

        volumes.push_back(trackVolume);
    } // gap */
    objVis.write("StrECAL");

    std::cout << "AMOUNT OF VOLUMES: " << volumes.size() << std::endl;

    // Sort volumes according to the center location
    std::sort(volumes.begin(), volumes.end(),
        [](const Acts::TrackingVolumePtr& lhs, const Acts::TrackingVolumePtr& rhs) {
            return lhs->center().z() < rhs->center().z();
        });
    
    // Print out volume position
    for (auto& volume : volumes) {
        std::cout << volume->volumeName() << " volume is at: " << volume->center().z() << std::endl;
    }

    // Glue volumes
    for (unsigned int i = 0; i < volumes.size() - 1; i++) {
        volumes[i + 1]->glueTrackingVolume(
            gctx, Acts::BoundarySurfaceFace::negativeFaceXY, volumes[i].get(),
            Acts::BoundarySurfaceFace::positiveFaceXY);
        volumes[i]->glueTrackingVolume(
            gctx, Acts::BoundarySurfaceFace::positiveFaceXY, volumes[i + 1].get(),
            Acts::BoundarySurfaceFace::negativeFaceXY);
    }

    Acts::Translation3 world_trans{0., 0., 0.}; 
    Acts::Transform3 world_trafo = Acts::Transform3(world_trans); 
    auto world_qsize = std::make_shared<const Acts::CuboidVolumeBounds>(
        70_cm * 0.5, 70_cm * 0.5, 150_cm * 0.5); // half lengths
    auto world_zsize = std::make_shared<const Acts::CylinderVolumeBounds>(
            0., 2568_mm * 0.5, 3150_mm * 0.5);

    // Vector of confined volumes
    std::vector<std::pair<Acts::TrackingVolumePtr, Acts::Vector3>> tapVec;
    tapVec.reserve(volumes.size());
    for (auto& tVol : volumes) {
        tapVec.push_back(std::make_pair(tVol, tVol->center()));
    }

    // Set bin boundaries along binning
    std::vector<float> binBoundaries;
    binBoundaries.push_back(volumes[0]->center().z() - 
       (StationD - (2 * 10_cm)) * 0.5);
    for (size_t i = 0; i < volumes.size(); i++) {
        if (i % 2 == 0) {
      		binBoundaries.push_back(volumes[i]->center().z() + 
                StationD * 0.5);
      	}
      	else if (i % 2 != 0) {
      		binBoundaries.push_back(volumes[i]->center().z() + 
                (StationD - (2 * 10_cm)) * 0.5);
      	}
    }

    // Build binning
    Acts::BinningData binData(Acts::BinningOption::open, 
        Acts::BinningValue::binZ, binBoundaries);
    auto bu = std::make_unique<const Acts::BinUtility>(binData);

    // Build Tracking Volume Array
    std::shared_ptr<const Acts::TrackingVolumeArray> trVolArr(
        new Acts::BinnedArrayXD<Acts::TrackingVolumePtr>(tapVec, std::move(bu)));
    
    // Create world volume
    Acts::MutableTrackingVolumePtr mtvp(Acts::TrackingVolume::create(
        world_trafo, world_qsize, trVolArr, "World"));
    // I'm not sure what is the purpose of this hook
    std::shared_ptr<const Acts::GeometryIdentifierHook> geometryIdentifierHook =
        std::make_shared<Acts::GeometryIdentifierHook>();

    return std::make_unique<Acts::TrackingGeometry>(
        mtvp, nullptr, *geometryIdentifierHook);
}
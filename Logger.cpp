// ---------------------------------------------------------------- //
//                Creator: Amaia Razquin Lizarraga                  //
//                Contact: amaiafisika@gmail.com                    //
//                Date:    7 September 2023                         //
// ---------------------------------------------------------------- //
/* This module takes the tracking geometry and outputs the information
   in a txt file or in the terminal to be able to check the geometry
   formation. */

#include <iostream>
#include <vector>

#include <Acts/Material/ISurfaceMaterial.hpp>

#include "Logger.hpp" 

// Outputs geometry information into a Geometry_log.txt
void aLogGeometry(std::shared_ptr<const Acts::TrackingGeometry> tGeometry, 
    const Acts::GeometryContext& gctx){

	std::fstream geomfile;
	geomfile.open("Geometry_log.txt", std::ios::out); 

  	const Acts::TrackingVolume* tVolume = tGeometry->highestTrackingVolume();
	geomfile << "Tracking VOLUME with GeometryID: " << tVolume->geometryId() << " and name: " 
		<< tVolume->volumeName() << "\n";
	geomfile << "VOLUME position: " << tVolume->center().x() << " " << tVolume->center().y()
		<< " " << tVolume->center().z() << "\n";
	geomfile << "VOLUME material: " << tVolume->volumeMaterial() << "\n";
  	if (tVolume->confinedVolumes()) {
    	for (auto volume : tVolume->confinedVolumes()->arrayObjects()) {
			geomfile << "Confined VOLUME with GeometryID: " << volume->geometryId() << " "
				<< volume->volumeName() << "\n";
			geomfile << "VOLUME position: " << volume->center().x() << " " << volume->center().y()
				<< " " << volume->center().z() << "\n";
			geomfile << "VOLUME material: " << volume->volumeMaterial() << "\n";
			auto vol_mat = volume->volumeMaterial()->material(volume->center());
			geomfile << "Volume material is:\n" <<
					"Relative atomic mass:            " << vol_mat.Ar() << "\n"  <<
					"Nuclear interaction length:      " << vol_mat.L0() << "\n"  << 
					"Mass density:                    " << vol_mat.massDensity() << "\n" <<
					"Mean electron excitation energy: " << vol_mat.meanExcitationEnergy() << "\n" << 
					"Radiation length:                " << vol_mat.X0() << "\n" <<
					"Nuclear charge number:           " << vol_mat.Z() << "\n";
      		if (volume->confinedLayers()) {
        		for (const auto& layer : volume->confinedLayers()->arrayObjects()) {
					geomfile << "Found a LAYER with GeometryID: " << layer->geometryId() << 
						" of type " << layer->layerType() << " and thickness: " 
						<< layer->thickness() << "\n";
					if (layer->layerType() == Acts::navigation) {continue;}
					for (auto surface : layer->approachDescriptor()->containedSurfaces()){
						geomfile << "Found an approach SURFACE with GeometryID: " <<
						surface->geometryId() << " at " << surface->center(gctx).x() << " " 
							<< surface->center(gctx).y() << " " << surface->center(gctx).z() << "\n";
						geomfile << "with bounds " << surface->bounds() << "\n" << 
							"and normal " << surface->normal(gctx).x() << " " << 
							surface->normal(gctx).y() << " " << surface->normal(gctx).z() << "\n";
					}
          			for (auto surface : layer->surfaceArray()->surfaces()) {
						geomfile << "Found a sensitive SURFACE with GeometryID:: " << surface->geometryId() 
							<< " of type " << surface->type() << " at " << surface->center(gctx).x() << " " 
							<< surface->center(gctx).y() << " " << surface->center(gctx).z() << "\n";
						geomfile << "with bounds " << surface->bounds() << "\n" << 
							"and normal " << surface->normal(gctx).x() << " " << 
							surface->normal(gctx).y() << " " << surface->normal(gctx).z() << "\n";
						geomfile << "SURFACE material " << surface->surfaceMaterial() << "\n";
						auto& mat_slab = surface->surfaceMaterial()->materialSlab(surface->center(gctx));
						auto surf_mat  = mat_slab.material();
						geomfile << "Surface material is:\n" <<
							"Relative atomic mass:            " << surf_mat.Ar() << "\n"  <<
							"Nuclear interaction length:      " << surf_mat.L0() << "\n"  << 
							"Mass density:                    " << surf_mat.massDensity() << "\n" <<
							"Mean electron excitation energy: " << surf_mat.meanExcitationEnergy() << "\n" << 
							"Radiation length:                " << surf_mat.X0() << "\n" <<
							"Nuclear charge number:           " << surf_mat.Z() << "\n";
          			}    // surfaces
        		}      // layers objects
      		}        // confined layers
    	}          // volumes objects
  	}            // confined volumes
	geomfile.close();
	std::cout << "ACTS::aLogGeometry: Tracking geometry log file has been written out." << std::endl;
}

// Loops over tracking geometry to output to terminal all the surfaces and their position
void aSurfaceStream(std::shared_ptr<const Acts::TrackingGeometry> tGeometry, 
    const Acts::GeometryContext& gctx) {

    std::vector<const Acts::Surface*> surfaces;
	int i = 0;
  	const Acts::TrackingVolume* tVolume = tGeometry->highestTrackingVolume();
  	if (tVolume->confinedVolumes()) {
    	for (auto volume : tVolume->confinedVolumes()->arrayObjects()) {
      		if (volume->confinedLayers()) {
        		for (const auto& layer : volume->confinedLayers()->arrayObjects()) {
					if (layer->layerType() == Acts::navigation) {continue;}
          			for (auto surface : layer->surfaceArray()->surfaces()) {
              			std::cout<< i << " GeometryID:: "<<surface->geometryId() 
				            << " " << surface->center(gctx).z() << std::endl;
          			}    // surfaces
					i++;
        		}      // layers objects
      		}        // confined layers
    	}          // volumes objects
  	}            // confined volumes
}
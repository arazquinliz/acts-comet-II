// ---------------------------------------------------------------- //
//                Creator: Amaia Razquin Lizarraga                  //
//                Contact: amaiafisika@gmail.com                    //
//                Date:    7 September 2023                         //
// ---------------------------------------------------------------- //
/* 
    This program calls ACTS and all of the modules created by the author
    It also executes some of the esential functions (geometry, magnetic field,...) 
    Command:
        g++ helloworld.cpp Geometry.cpp MagField.cpp Measurements.cpp Seeder.cpp 
        Logger.cpp KalmanFitting.cpp Writer.cpp Estimator.cpp -o helloworld.x 
        -L/path/to/acts_install/lib -lActsPluginTGeo -lActsPluginJson -lActsCore 
        -I/path/to/acts_install/include/ -I$(root-config --incdir) 
        -I$(root-config --evelibs) $(root-config --libs) 
        -lGeom -I/usr/include/eigen3
*/

// STD
#include <iostream>
#include <ostream>
#include <vector>
#include <map>

// ACTS
// Definitions & utilities headers
#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/Utilities/CalibrationContext.hpp>

// Geometry headers
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/Geometry/GeometryHierarchyMap.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>

// Magnetic field headers
#include <Acts/MagneticField/ConstantBField.hpp>

// Event data headers
#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/EventData/SourceLink.hpp>
#include <Acts/EventData/VectorMultiTrajectory.hpp>
#include <Acts/EventData/VectorTrackContainer.hpp>

// Propagator headers
#include <Acts/Propagator/Propagator.hpp>
#include <Acts/Propagator/EigenStepper.hpp>
#include <Acts/Propagator/Navigator.hpp>
#include <Acts/Propagator/MaterialInteractor.hpp>
#include <Acts/Propagator/StandardAborters.hpp>
#include <Acts/Propagator/detail/SteppingLogger.hpp>

// Track fitting headers
#include <Acts/TrackFitting/GainMatrixSmoother.hpp>
#include <Acts/TrackFitting/GainMatrixUpdater.hpp>
#include <Acts/TrackFitting/KalmanFitter.hpp>

// Visualization
#include <Acts/Visualization/ObjVisualization3D.hpp>
#include <Acts/Visualization/PlyVisualization3D.hpp>
#include <Acts/Visualization/GeometryView3D.hpp>
#include <Acts/Visualization/ViewConfig.hpp>

// SVG
#include <Acts/Plugins/ActSVG/LayerSvgConverter.hpp>
#include <Acts/Plugins/ActSVG/SvgUtils.hpp>
#include <Acts/Plugins/ActSVG/TrackingGeometrySvgConverter.hpp>

// Amaia:
#include "base/SpacePoint.hpp"
#include "base/ActsSourceLink.hpp"
#include "base/Calibrator.hpp"

#include "Geometry.hpp"
#include "Measurements.hpp"
#include "MagField.hpp"
#include "Seeder.hpp"
#include "Logger.hpp"
#include "KalmanFitting.hpp"
#include "Writer.hpp"
#include "Estimator.hpp"

int main() {
    std::cout << "Hello world: I am alive." << std::endl;

    // Initialize geometry context
    Acts::GeometryContext* gctx;
    bool debug = true; // Get some extra information

    std::shared_ptr<const Acts::TrackingGeometry> tGeometry;
    tGeometry = aBuildDiskGeometry(*gctx, debug);
    // tGeometry = aBuildStrawGeometry(*gctx, debug);

    if (!tGeometry) {
    		throw std::runtime_error("ACTS::BuilGeometry: tGeometry is null");
    }

    // Initialize magnetic field context
	Acts::MagneticFieldContext magctx;

	// Create the object for the fieldmap binning
	auto mapBins = [](std::array<size_t, 3> bins, std::array<size_t, 3> sizes) {
        	return (bins[0] * (sizes[1] * sizes[2]) + bins[1] * sizes[2] + bins[2]);
    };
    // Get magnetic field
	auto BFieldmap = aGetMagneticField(debug, magctx, std::move(mapBins));
	auto bField = std::dynamic_pointer_cast<const Acts::InterpolatedMagneticField>(BFieldmap);
    if (!bField) {
    	throw std::runtime_error("ACTS::GetMagneticField: bField is null");
    }  
    auto bCache = BFieldmap->makeCache(magctx);

    std::cout << "----------------------------------------------------------------------\n"
        << "-------------- Everything seems fine, congratulations! ---------------\n"
        << "----------------------------------------------------------------------" << std::endl;
}

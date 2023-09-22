// ---------------------------------------------------------------- //
//                Creator: Amaia Razquin Lizarraga                  //
//                Contact: amaiafisika@gmail.com                    //
//                Date:    7 September 2023                         //
// ---------------------------------------------------------------- //
/* This module gets the magnetic field vector in each position from a 
   fieldmap root file (in get_mf()). The resulting map is printed out 
   via writeout_mf() as magnetic_field.root, then it can be visualized
   with printBField.c. The magnetic field is introduced in the provider
   in aGetMagneticField, the main function. */

#include <map>

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Seeding/EstimateTrackParamsFromSeed.hpp>
 
#include <TRandom3.h>
#include <TFile.h>
#include <TTree.h>

#include "base/SpacePoint.hpp"
#include "Seeder.hpp"

using namespace Acts::UnitLiterals;

// Compare if two seeds share two spacepoints and return the vector with
// a vector with the shared spacepoints 
std::vector<SpacePoint> compareSeeds(Acts::Seed<SpacePoint> seedA, 
    Acts::Seed<SpacePoint> seedB) {
    
    int coin = 0; // coincidence (max = 2)
    std::vector<int> all = {0, 1, 2};
    std::vector<int> inA;
    std::vector<int> inB;
    std::vector<SpacePoint> seed;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            if (seedA.sp()[i] == seedB.sp()[j]) {
                seed.push_back(*(seedA.sp()[i]));
                coin++;
                inA.push_back(i);
                inB.push_back(j);
            }   
        }
    }
    if (coin == 2) {
        for (auto& a : all) {
            if (std::find(inA.begin(), inA.end(), a) == inA.end()) {
                seed.push_back(*(seedA.sp()[a]));
            }
            if (std::find(inB.begin(), inB.end(), a) == inB.end()) {
                seed.push_back(*(seedB.sp()[a]));
            }
        }
        return seed;
    }
    else {
        std::vector<SpacePoint> emptySeed;
        return emptySeed;
    }
}

// Checks if a point is within a list
bool isIn(std::vector<Acts::Vector3> vector, Acts::Vector3 point) {
    bool inside = false;
    for (auto p : vector) {
        if ((abs(p.x() - point.x()) < 0.0001) && (abs(p.y() - point.y()) < 0.0001) && 
            (abs(p.z() - point.z()) < 0.0001)) {
            inside = true;
            break;
        }
    }
    return inside;
} // Overload function

bool isIn(std::vector<SpacePoint> vector, Acts::Vector3 point) {
    bool inside = false;
    for (auto p : vector) {
        if ((abs(p.x() - point.x()) < 0.0001) && (abs(p.y() - point.y()) < 0.0001) && 
            (abs(p.z() - point.z()) < 0.0001)) {
            inside = true;
            break;
        }
    }
    return inside;
}

bool isIn(std::vector<Acts::Vector2> vector, Acts::Vector2 point) {
    bool inside = false;
    for (auto p : vector) {
        if ((abs(p.x() - point.x()) < 0.0001) && (abs(p.y() - point.y()) < 0.0001)) {
            inside = true;
            break;
        }
    }
    return inside;
}

// Gets the event with biggest presence in the vector and returns the purity
// i.e. number of SP of that event w.r.t the vector size
std::pair<int, double> checkPurity(std::vector<SpacePoint> seed) {

    int nHits = 0;
    std::map<int, int> purityMap; // event, N sp
    for (auto& sp : seed) {
        purityMap[sp.event()] += 1;
        nHits++;
    }

    int mainEv;
    int mainN = 0;
    for (auto const& [event, N] : purityMap)
    {
        if (N > mainN) {
            mainEv = event;
            mainN  = N;
        }
    }

    return std::make_pair(mainEv, mainN * 1.0 / nHits);
}

// Helper function to use the track parameter estimator within ACTS
Acts::BoundVector estimateParams(const Acts::GeometryContext& gctx, 
    std::shared_ptr<const Acts::TrackingGeometry> tGeometry,
    Acts::MagneticFieldContext magctx, 
    std::shared_ptr<Acts::MagneticFieldProvider> BFieldmap,
    std::vector<const SpacePoint*> seedSP) {

    auto logger = Acts::getDefaultLogger("ParamEstimator", Acts::Logging::VERBOSE);

    //std::sort(seedSP.begin(), seedSP.end(), [](const SpacePoint* lhs, 
        //const SpacePoint* rhs) {return (lhs->z() < rhs->z());});
    
    const auto rSurface        = tGeometry->findSurface(seedSP[0]->geoID());
    auto bCache    = BFieldmap->makeCache(magctx);
    auto field_res = BFieldmap->getField(Acts::Vector3{seedSP[0]->x(), 
                    seedSP[0]->y(), seedSP[0]->z()}, bCache);
    auto field     = field_res.value();
    Acts::ActsScalar bFieldMin = 0;
    Acts::ActsScalar eMass     = 0.511_MeV;

    auto trackParams = (Acts::estimateTrackParamsFromSeed(gctx, seedSP.begin(), 
        seedSP.end(), *rSurface, field, bFieldMin, *logger, eMass)).value();
    //double p_estim = - 1 / trackParams[Acts::eBoundQOverP] * 1000;

    return trackParams;
}

// Takes a position an creates a SpacePoint class object
// Uses code from sPHENIX 
std::unique_ptr<SpacePoint> makeSpacePoint(const Acts::GeometryContext& gctx,
    Acts::Vector3 globalPos, 
    double time, 
    std::shared_ptr<const Acts::Surface> surf,
    const int event, const Acts::Vector3 p) {

    bool debug = true;

    float x = globalPos.x() - 1080._cm;
    float y = globalPos.y();
    float z = globalPos.z() - 333.45_cm;
    float r = std::sqrt(x * x + y * y);

    Acts::Vector3 mom(1, 1, 1); 
    Acts::SquareMatrix2 localCov = Acts::SquareMatrix2::Zero();
    localCov(0, 0) = 1._mm;
    localCov(1, 1) = 1._mm;

    /// The space point requires only the variance of the transverse and
    /// longitudinal position. Reduce computations by transforming the
    /// covariance directly from local to r/z.
    ///
    /// compute Jacobian from global coordinates to r/z
    ///
    ///         r = sqrt(x² + y²)
    /// dr/d{x,y} = (1 / sqrt(x² + y²)) * 2 * {x,y}
    ///             = 2 * {x,y} / r
    ///       dz/dz = 1 
    Acts::RotationMatrix3 rotLocalToGlobal = surf->referenceFrame(gctx, 
        Acts::Vector3{x, y, z}, mom);
    auto scale = 2.0 / std::hypot(x, y);
    Acts::ActsMatrix<2, 3> jacXyzToRhoZ = Acts::ActsMatrix<2, 3>::Zero();
    jacXyzToRhoZ(0, Acts::ePos0) = scale * x;
    jacXyzToRhoZ(0, Acts::ePos1) = scale * y;
    jacXyzToRhoZ(0, Acts::ePos2) = 1;
    // Compute Jacobian from local coordinates to rho/z
    Acts::ActsMatrix<2, 2> jac   = jacXyzToRhoZ * rotLocalToGlobal.block<3, 2>(
        Acts::ePos0, Acts::ePos0);
    // Compute rho/z variance
    Acts::ActsVector<2> var      = (jac * localCov * jac.transpose()).diagonal();

    /*
    * From Acts v17 to v19 the scattering uncertainty value allowed was changed.
    * This led to a decrease in efficiency. To offset this, we scale the 
    * uncertainties by a tuned factor that gives the v17 performance
    * Track reconstruction is an art as much as it is a science...
    */
    double uncfactor = 3.175; // from sPHENIX
    std::unique_ptr<SpacePoint> spPtr(new SpacePoint{x, y, z, r, time, surf->geometryId(), 
        var[0] * uncfactor, var[1] * uncfactor, event, p});

    return spPtr;
}

// Generate random backgrounds per layer
void getRandomBG(std::vector<const SpacePoint*>& spVector, 
    const Acts::GeometryContext& gctx, Acts::Extent& rRangeSPExtent,
    std::vector<std::shared_ptr<const Acts::Surface>> station) {
    
    TRandom3* rndm = new TRandom3();
    rndm->SetSeed(0);

    double disk_rad = 1390 * 0.5;
    auto zPositions = {-1380.25, -1370.9, -1359.1, -1349.75};
    int surf_num    = 0;
    int nHits       = 0;
    double time     = 0;
    for (auto& surface : station) {
        for (int i = 0; i < 10; i++) {
            double x = rndm->Uniform(-1*disk_rad, disk_rad);
            double y = rndm->Uniform(-1*disk_rad, disk_rad);
            double z = surface->center(gctx).z(); //  + rndm->Uniform(-1*radius, radius)
            
            auto spacePoint = makeSpacePoint(gctx, 
                Acts::Vector3{x + 1080._cm, y, z + 333.45_cm}, time, surface, 800000, Acts::Vector3::Zero()).release();
            spVector.push_back(spacePoint);
            rRangeSPExtent.extend({spacePoint->x(), spacePoint->y(), spacePoint->z()});
            nHits++;
        }
        surf_num++;
    }
}

// Get real backgrounds from a truth MC background file
void getRealBG(std::vector<const SpacePoint*>& spVector, 
    const Acts::GeometryContext& gctx, Acts::Extent& rRangeSPExtent,
    std::map<std::string, std::vector<std::shared_ptr<const Acts::Surface>>> tracker_layout,
    int stationID) {
    //std::cout << "Getting real BG" << std::endl;
    
    TRandom3* rndm = new TRandom3();
    rndm->SetSeed(0);

    // get the tree and iterate over the branches
	auto *f_BG = TFile::Open("detector_truth_all_MC6.root");
	if (!f_BG) {
		std::cout << "Couldn't open measurement file\n";
	}

	auto *tree_BG = (TTree*)f_BG->Get("detections");
	std::vector<ROOT::Math::XYZVector> *spOfEvent_BG = 0;
    std::vector<ROOT::Math::XYZVector> *momOfEvent_BG = 0;
    std::vector<double> *timeOfEvent = 0;
	std::vector<double> *eDeposition = 0;
	std::vector<std::string> *particleID = 0;
	tree_BG->SetBranchAddress("spOfEvent", &spOfEvent_BG);
    tree_BG->SetBranchAddress("momOfEvent", &momOfEvent_BG);
    tree_BG->SetBranchAddress("timeOfEvent", &timeOfEvent);
    tree_BG->SetBranchAddress("eDeposition", &eDeposition);
    tree_BG->SetBranchAddress("particleID", &particleID);

    unsigned int nHits = 0;
    int NofStations    = 8; 
    double radius      = 0.49_cm;      // StrawTubeOuterRadius
	double StationD    = 39.0_cm;      // Distance between stations
	double layerD      = 0.935_cm;     // Distance between layers   
    double separation  = 1.0575_cm;    // Separation to build geometry
    double thickness   = 2 * separation + layerD + radius * 2;
    auto StationName   = "Station_" + std::to_string(stationID);
    
    for (size_t ev = 0; ev < tree_BG->GetEntries(); ++ev) {
            tree_BG->GetEntry(ev);

            std::vector<Acts::Vector3> posVec = {Acts::Vector3::Zero(), Acts::Vector3::Zero(), 
                Acts::Vector3::Zero(), Acts::Vector3::Zero()};
            std::vector<double> timeVec = {0, 0 , 0, 0};
            std::vector<int> k = {0, 0, 0, 0};
            
            if ((*spOfEvent_BG).size() == 0) {continue;}
            for (size_t it = 0; it < (*spOfEvent_BG).size(); ++it) {

                double positionZ = (*spOfEvent_BG)[it].Z() - 333.45_cm;
                auto& surface_1 = tracker_layout[StationName][0];
                auto& surface_2 = tracker_layout[StationName][1];
                auto& surface_3 = tracker_layout[StationName][2];
                auto& surface_4 = tracker_layout[StationName][3];
                
                double time = (*timeOfEvent)[it];
                // Time window cut
                if ((time < 650) || (time > 1200)) {continue;}

                // Energy deposition cut
                if (((*eDeposition)[it] < 0) || ((*eDeposition)[it] > 6)) {continue;}

                if (positionZ >= (-1380.25 + StationD * stationID + thickness)) {continue;}
                if (std::abs(positionZ-surface_1->center(gctx).z()) <= radius) {
                    posVec[0].x() += (*spOfEvent_BG)[it].X();
                    posVec[0].y() += (*spOfEvent_BG)[it].Y();
                    posVec[0].z() += (*spOfEvent_BG)[it].Z();
                    timeVec[0] += (*timeOfEvent)[it];
                    k[0] += 1;
                }
                else if (std::abs(positionZ-surface_2->center(gctx).z()) <= radius) {
                    posVec[1].x() += (*spOfEvent_BG)[it].X();
                    posVec[1].y() += (*spOfEvent_BG)[it].Y();
                    posVec[1].z() += (*spOfEvent_BG)[it].Z();
                    timeVec[1] += (*timeOfEvent)[it];
                    k[1] += 1;
                }
                else if (std::abs(positionZ-surface_3->center(gctx).z()) <= radius) {
                    posVec[2].x() += (*spOfEvent_BG)[it].X();
                    posVec[2].y() += (*spOfEvent_BG)[it].Y();
                    posVec[2].z() += (*spOfEvent_BG)[it].Z();
                    timeVec[2] += (*timeOfEvent)[it];
                    k[2] += 1;
                }
                else if (std::abs(positionZ-surface_4->center(gctx).z()) <= radius) {
                    posVec[3].x() += (*spOfEvent_BG)[it].X();
                    posVec[3].y() += (*spOfEvent_BG)[it].Y();
                    posVec[3].z() += (*spOfEvent_BG)[it].Z();
                    timeVec[3] += (*timeOfEvent)[it];
                    k[3] += 1;
                }
            }

            auto p = Acts::Vector3{(*momOfEvent_BG)[0].X(), (*momOfEvent_BG)[0].Y(), (*momOfEvent_BG)[0].Z()};
            int kount = 0;
            for (auto pos : posVec) {
                for (auto& surface : tracker_layout[StationName]) {
                    Acts::Vector3 surf_center = surface->center(gctx);

                    if (k[kount] == 0) {continue;}

                    if (abs(pos.z()/k[kount] - 333.45_cm - surf_center.z()) < radius) {

                        // 20% acceptance (currently deactivated)
                        double acc = rndm->Uniform(0, 1);
                        //if (acc > 0.2) {continue;}

                        auto spacePoint = makeSpacePoint(gctx, 
                            Acts::Vector3{pos.x()/k[kount], pos.y()/k[kount], surf_center.z() + 333.45_cm}, 
                            timeVec[kount]/k[kount], surface, std::stol((*particleID)[0]), p).release();
                        //std::cout << pos.x()/k[kount]-10800 << " " << pos.y()/k[kount] << " " << 
                            //pos.z()/k[kount]-3334.5 << " " << timeVec[kount]/k[kount] << std::endl;
                        spVector.push_back(spacePoint);
                        rRangeSPExtent.extend({spacePoint->x(), spacePoint->y(), spacePoint->z()});
                        nHits++;
                    }
                } // surfaces in the first station
                kount++;
            } // spacepoints
    } // loop over events
    f_BG->Close();
}

// Get spacepoint from data vectors from root file
std::vector<const SpacePoint*> getSpacePoints(
    const Acts::GeometryContext& gctx,
    std::vector<ROOT::Math::XYZVector> spOfEvent,
    std::vector<ROOT::Math::XYZVector> momOfEvent,
    std::vector<double> eDeposition,
    std::vector<double> timeOfEvent, double creationT,
    std::map<std::string, std::vector<std::shared_ptr<const Acts::Surface>>> tracker_layout,
    int eventID, int stationID, Acts::Extent& rRangeSPExtent) {
    
    int NofStations   = 8; 
    double radius     = 0.49_cm;      // StrawTubeOuterRadius
	double StationD   = 39.0_cm;      // Distance between stations
	double layerD     = 0.935_cm;     // Distance between layers   
    double separation = 1.0575_cm;    // Separation to build geometry
    double thickness  = 2 * separation + layerD + radius * 2;
    auto StationName  = "Station_" + std::to_string(stationID);

    std::vector<const SpacePoint*> spVector;

    std::vector<Acts::Vector3> posVec = {Acts::Vector3::Zero(), Acts::Vector3::Zero(), 
            Acts::Vector3::Zero(), Acts::Vector3::Zero()};
    std::vector<double> timeVec = {0, 0 , 0, 0};
    std::vector<int> k = {0, 0, 0, 0};

    for (size_t it = 0; it < spOfEvent.size(); ++it) {
        
        double positionZ = spOfEvent[it].Z() - 333.45_cm;
        auto& surface_1 = tracker_layout[StationName][0];
        auto& surface_2 = tracker_layout[StationName][1];
        auto& surface_3 = tracker_layout[StationName][2];
        auto& surface_4 = tracker_layout[StationName][3];
        
        if (positionZ >= (-1380.25 + StationD * stationID + thickness)) {continue;}
        // Energy deposition cut
        if ((eDeposition[it] * 1000 < 0) || (eDeposition[it] * 1000 > 6)) {continue;}

        if (std::abs(positionZ-surface_1->center(gctx).z()) <= radius) {
            posVec[0].x() += spOfEvent[it].X();
            posVec[0].y() += spOfEvent[it].Y();
            posVec[0].z() += spOfEvent[it].Z();
            timeVec[0] += timeOfEvent[it];
            k[0] += 1;
        }
        else if (std::abs(positionZ-surface_2->center(gctx).z()) <= radius) {
            posVec[1].x() += spOfEvent[it].X();
            posVec[1].y() += spOfEvent[it].Y();
            posVec[1].z() += spOfEvent[it].Z();
            timeVec[1] += timeOfEvent[it];
            k[1] += 1;
        }
        else if (std::abs(positionZ-surface_3->center(gctx).z()) <= radius) {
            posVec[2].x() += spOfEvent[it].X();
            posVec[2].y() += spOfEvent[it].Y();
            posVec[2].z() += spOfEvent[it].Z();
            timeVec[2] += timeOfEvent[it];
            k[2] += 1;
        }
        else if (std::abs(positionZ-surface_4->center(gctx).z()) <= radius) {
            posVec[3].x() += spOfEvent[it].X();
            posVec[3].y() += spOfEvent[it].Y();
            posVec[3].z() += spOfEvent[it].Z();
            timeVec[3] += timeOfEvent[it];
            k[3] += 1;
        }
    }

    auto p = Acts::Vector3{momOfEvent[0].X(), momOfEvent[0].Y(), momOfEvent[0].Z()};

    int kount = 0;
    for (auto pos : posVec) {
        for (auto& surface : tracker_layout[StationName]) {
            Acts::Vector3 surf_center = surface->center(gctx);
            if (k[kount] == 0) {continue;}
            if (abs(pos.z()/k[kount] - 333.45_cm - surf_center.z()) < radius) {
                auto spacePoint = makeSpacePoint(gctx, 
                    Acts::Vector3{pos.x()/k[kount], pos.y()/k[kount], surf_center.z() + 333.45_cm}, 
                        timeVec[kount]/k[kount] + creationT, surface, eventID, p).release();
                double time = spacePoint->t();
                if ((time < 650) || (time > 1200)) {std::cout << "AHHHHHHHH SOMETHING IS WRONG " << time << std::endl; continue;}
                spVector.push_back(spacePoint);
                rRangeSPExtent.extend({spacePoint->x(), spacePoint->y(), spacePoint->z()});
            }
        } // surfaces in the first station
        kount++;
    } // spacepoints
    
    const auto station = tracker_layout[StationName];
    //getRandomBG(spVector, gctx, rRangeSPExtent, station);
    getRealBG(spVector, gctx, rRangeSPExtent, tracker_layout, stationID);

    return spVector;
} 

// Get spacepoint from data vectors from root file
// Function overload for energy deposition cut test (if not enough deposition)
std::vector<const SpacePoint*> getSpacePointsTC(
    const Acts::GeometryContext& gctx,
    std::vector<ROOT::Math::XYZVector> spOfEvent,
    std::vector<ROOT::Math::XYZVector> momOfEvent,
    std::vector<double> eDeposition,
    std::vector<double> timeOfEvent, double creationT,
    std::map<std::string, std::vector<std::shared_ptr<const Acts::Surface>>> tracker_layout,
    int eventID, int stationID, Acts::Extent& rRangeSPExtent) {
    
    int NofStations   = 8; 
    double radius     = 0.49_cm;      // StrawTubeOuterRadius
	double StationD   = 39.0_cm;      // Distance between stations
	double layerD     = 0.935_cm;     // Distance between layers   
    double separation = 1.0575_cm;    // Separation to build geometry
    double thickness  = 2 * separation + layerD + radius * 2;
    auto StationName  = "Station_" + std::to_string(stationID);

    std::vector<const SpacePoint*> spVector;

    std::vector<Acts::Vector3> posVec = {Acts::Vector3::Zero(), Acts::Vector3::Zero(), 
            Acts::Vector3::Zero(), Acts::Vector3::Zero()};
    std::vector<double> timeVec = {0, 0 , 0, 0};
    std::vector<int> k = {0, 0, 0, 0};

    for (size_t it = 0; it < spOfEvent.size(); ++it) {
        
        double positionZ = spOfEvent[it].Z() - 333.45_cm;
        auto& surface_1 = tracker_layout[StationName][0];
        auto& surface_2 = tracker_layout[StationName][1];
        auto& surface_3 = tracker_layout[StationName][2];
        auto& surface_4 = tracker_layout[StationName][3]; 
        
        if (positionZ >= (-1380.25 + StationD * stationID + thickness)) {continue;}
        if ((eDeposition[it] * 1000 < 0) || (eDeposition[it] * 1000 > 6)) {continue;}

        if (std::abs(positionZ-surface_1->center(gctx).z()) <= radius) {
            posVec[0].x() += spOfEvent[it].X();
            posVec[0].y() += spOfEvent[it].Y();
            posVec[0].z() += spOfEvent[it].Z();
            timeVec[0] += timeOfEvent[it];
            k[0] += 1;
        }
        else if (std::abs(positionZ-surface_2->center(gctx).z()) <= radius) {
            posVec[1].x() += spOfEvent[it].X();
            posVec[1].y() += spOfEvent[it].Y();
            posVec[1].z() += spOfEvent[it].Z();
            timeVec[1] += timeOfEvent[it];
            k[1] += 1;
        }
        else if (std::abs(positionZ-surface_3->center(gctx).z()) <= radius) {
            posVec[2].x() += spOfEvent[it].X();
            posVec[2].y() += spOfEvent[it].Y();
            posVec[2].z() += spOfEvent[it].Z();
            timeVec[2] += timeOfEvent[it];
            k[2] += 1;
        }
        else if (std::abs(positionZ-surface_4->center(gctx).z()) <= radius) {
            posVec[3].x() += spOfEvent[it].X();
            posVec[3].y() += spOfEvent[it].Y();
            posVec[3].z() += spOfEvent[it].Z();
            timeVec[3] += timeOfEvent[it];
            k[3] += 1;
        }
    }

    auto p = Acts::Vector3{momOfEvent[0].X(), momOfEvent[0].Y(), momOfEvent[0].Z()};

    int kount = 0;
    for (auto pos : posVec) {
        for (auto& surface : tracker_layout[StationName]) {
            Acts::Vector3 surf_center = surface->center(gctx);
            if (k[kount] == 0) {continue;}
            if (abs(pos.z()/k[kount] - 333.45_cm - surf_center.z()) < radius) {
                auto spacePoint = makeSpacePoint(gctx, 
                    Acts::Vector3{pos.x()/k[kount], pos.y()/k[kount], surf_center.z() + 333.45_cm}, 
                        timeVec[kount]/k[kount] + creationT, surface, eventID, p).release();
                spVector.push_back(spacePoint);
                rRangeSPExtent.extend({spacePoint->x(), spacePoint->y(), spacePoint->z()});
            }
        } // surfaces in the first station
        kount++;
    } // spacepoints
    
    const auto station = tracker_layout[StationName];
    //getRandomBG(spVector, gctx, rRangeSPExtent, station);
    getRealBG(spVector, gctx, rRangeSPExtent, tracker_layout, stationID);

    return spVector;
} 
// ---------------------------------------------------------------- //
//                Creator: Amaia Razquin Lizarraga                  //
//                Contact: amaiafisika@gmail.com                    //
//                Date:    7 September 2023                         //
// ---------------------------------------------------------------- //
/* This module serves the Propagator program to write out the information
   of the propagation steps (aPropagationStepsWriter) and the material
   interactions in the detector regions (aMaterialTrackWriter). Both 
   functions are a copy of the existing writers in ACTS*/

#include <ios>
#include <iostream>
#include <ostream>
#include <stdexcept>

#include "Acts/Propagator/detail/SteppingLogger.hpp"
#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <Acts/Geometry/TrackingVolume.hpp>
#include <Acts/Propagator/ConstrainedStep.hpp>
#include <Acts/Surfaces/Surface.hpp>

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include <Acts/Propagator/MaterialInteractor.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <TFile.h>
#include <TTree.h>

#include "Writer.hpp"

using RecordedMaterial = Acts::MaterialInteractor::result_type;
// start position, momentum and recorded material
using RecordedMaterialTrack = 
    std::pair<std::pair<Acts::Vector3, Acts::Vector3>, RecordedMaterial>;

void aPropagationStepsWriter(
    std::vector<std::vector<Acts::detail::Step>>& stepCollection){
    
    std::fstream measfile;
    measfile.open("propagation_steps.txt", std::ios::out); 
    // All input options
    int m_eventNr = 0;               ///< the event number of
    std::vector<int> m_volumeID;     ///< volume identifier
    std::vector<int> m_boundaryID;   ///< boundary identifier
    std::vector<int> m_layerID;      ///< layer identifier if
    std::vector<int> m_approachID;   ///< surface identifier
    std::vector<int> m_sensitiveID;  ///< surface identifier
    std::vector<int> m_material;     ///< flag material if present
    std::vector<float> m_x;          ///< global x
    std::vector<float> m_y;          ///< global y
    std::vector<float> m_z;          ///< global z
    std::vector<float> m_r;          ///< global r
    std::vector<float> m_dx;         ///< global direction x    
    std::vector<float> m_dy;         ///< global direction y
    std::vector<float> m_dz;         ///< global direction z
    std::vector<int> m_step_type;    ///< step type
    std::vector<float> m_step_acc;   ///< accuracy
    std::vector<float> m_step_act;   ///< actor check
    std::vector<float> m_step_abt;   ///< aborter
    std::vector<float> m_step_usr;   ///< user
    std::vector<size_t> 
         m_nStepTrials;  ///< Number of iterations needed by the stepsize
                      ///  finder (e.g. Runge-Kutta) of the stepper.
    // Setup Root I/O
    TFile *m_outputFile = new TFile("PropagationOutput.root", "Recreate");
    TTree *m_outputTree = new TTree("propagation_steps", "TTree from PropagationStepsWriter");

    // Set the branches
    m_outputTree->Branch("event_nr", &m_eventNr);
    m_outputTree->Branch("volume_id", &m_volumeID);
    m_outputTree->Branch("boundary_id", &m_boundaryID);
    m_outputTree->Branch("layer_id", &m_layerID);
    m_outputTree->Branch("approach_id", &m_approachID);
    m_outputTree->Branch("sensitive_id", &m_sensitiveID);
    m_outputTree->Branch("material", &m_material);
    m_outputTree->Branch("g_x", &m_x);
    m_outputTree->Branch("g_y", &m_y);
    m_outputTree->Branch("g_z", &m_z);
    m_outputTree->Branch("g_r", &m_r);
    m_outputTree->Branch("d_x", &m_dx);
    m_outputTree->Branch("d_y", &m_dy);
    m_outputTree->Branch("d_z", &m_dz);
    m_outputTree->Branch("type", &m_step_type);
    m_outputTree->Branch("step_acc", &m_step_acc);
    m_outputTree->Branch("step_act", &m_step_act);
    m_outputTree->Branch("step_abt", &m_step_abt);
    m_outputTree->Branch("step_usr", &m_step_usr);
    m_outputTree->Branch("nStepTrials", &m_nStepTrials);

    // Get the event number from input
    //m_eventNr = N;
    // loop over the step vector of each test propagation in this
    int n_steps = 0;
    for (auto& steps : stepCollection) {
        int n_step = 0;
        // clear the vectors for each collection
        m_volumeID.clear();
        m_boundaryID.clear();
        m_layerID.clear();
        m_approachID.clear();
        m_sensitiveID.clear();
        m_material.clear();
        m_x.clear();
        m_y.clear();
        m_z.clear();
        m_r.clear();
        m_dx.clear();
        m_dy.clear();
        m_dz.clear();
        m_step_type.clear();
        m_step_acc.clear();
        m_step_act.clear();
        m_step_abt.clear();
        m_step_usr.clear();
        m_nStepTrials.clear();
        // loop over single steps
        for (auto& step : steps) {
            // the identification of the step
            Acts::GeometryIdentifier::Value volumeID = 0;
            Acts::GeometryIdentifier::Value boundaryID = 0;
            Acts::GeometryIdentifier::Value layerID = 0;
            Acts::GeometryIdentifier::Value approachID = 0;
            Acts::GeometryIdentifier::Value sensitiveID = 0;
            int material = 0;
            // get the identification from the surface first
            if (step.surface) {
                auto geoID  = step.surface->geometryId();
                volumeID    = geoID.volume();
                boundaryID  = geoID.boundary();
                layerID     = geoID.layer();
                approachID  = geoID.approach();
                sensitiveID = geoID.sensitive();
                if (step.surface->surfaceMaterial() != nullptr) {
                    material = 1;
                }
            }
            // a current volume overwrites the surface tagged one
            if (step.volume != nullptr) {
                volumeID = step.volume->geometryId().volume();
            }
            // now fill
            m_sensitiveID.push_back(sensitiveID);
            m_approachID.push_back(approachID);
            m_layerID.push_back(layerID);
            m_boundaryID.push_back(boundaryID);
            m_volumeID.push_back(volumeID);
            m_material.push_back(material);

            // kinematic information
            m_x.push_back(step.position.x());
            m_y.push_back(step.position.y());
            m_z.push_back(step.position.z());
            m_r.push_back(sqrt(step.position.x() * step.position.x() + 
                step.position.y() + step.position.y()));
            measfile << step.position.x() << " " << step.position.y() << " " << step.position.z()  << " " << 
                step.momentum.x() << " " << step.momentum.y() << " " << step.momentum.z() << 
                //sqrt(step.position.x() * step.position.x() + step.position.y() + step.position.y()) << 
                " " << volumeID << " " << layerID << " " << sensitiveID << "\n";
            auto direction = step.momentum.normalized();
            m_dx.push_back(direction.x());
            m_dy.push_back(direction.y());
            m_dz.push_back(direction.z());

            double accuracy = step.stepSize.value(Acts::ConstrainedStep::accuracy);
            double actor = step.stepSize.value(Acts::ConstrainedStep::actor);
            double aborter = step.stepSize.value(Acts::ConstrainedStep::aborter);
            double user = step.stepSize.value(Acts::ConstrainedStep::user);
            double actAbs = std::abs(actor);
            double accAbs = std::abs(accuracy);
            double aboAbs = std::abs(aborter);
            double usrAbs = std::abs(user);

            // todo - fold with direction
            if (actAbs < accAbs && actAbs < aboAbs && actAbs < usrAbs) {
                m_step_type.push_back(0);
            } else if (accAbs < aboAbs && accAbs < usrAbs) {
                m_step_type.push_back(1);
            } else if (aboAbs < usrAbs) {
                m_step_type.push_back(2);
            } else {
                m_step_type.push_back(3);
            }

            // step size information
            m_step_acc.push_back(Acts::clampValue<float>(accuracy));
            m_step_act.push_back(Acts::clampValue<float>(actor));
            m_step_abt.push_back(Acts::clampValue<float>(aborter));
            m_step_usr.push_back(Acts::clampValue<float>(user));

            // stepper efficiency
            m_nStepTrials.push_back(step.stepSize.nStepTrials);

            n_step++;
        }
        std::cout << "Steps = " << n_steps << " with number of Step = " << n_step << std::endl;
        n_steps++;
        m_outputTree->Fill();
    }
    m_outputTree->Write();
    m_outputTree->Print();
    m_outputFile->Close();
    measfile.close();
}

void aMaterialTrackWriter(
    const std::unordered_map<size_t, Acts::RecordedMaterialTrack>&
    materialTracks, const Acts::GeometryContext& gctx){
    
    // Configuration options
    /// Re-calculate total values from individual steps (for cross-checks)
    bool recalculateTotals = false;
    /// Write aut pre and post step (for G4), otherwise central step position
    bool prePostStep       = false;
    /// Write the surface to which the material step correpond
    bool storeSurface      = false;
    /// Write the volume to which the material step correpond
    bool storeVolume       = false;
    // Variables
    uint32_t m_eventId = 0;

    float m_v_x = 0;    ///< start global x
    float m_v_y = 0;    ///< start global y
    float m_v_z = 0;    ///< start global z
    float m_v_px = 0;   ///< start global momentum x
    float m_v_py = 0;   ///< start global momentum y
    float m_v_pz = 0;   ///< start global momentum z
    float m_v_phi = 0;  ///< start phi direction
    float m_v_eta = 0;  ///< start eta direction
    float m_tX0 = 0;    ///< thickness in X0/L0
    float m_tL0 = 0;    ///< thickness in X0/L0

    std::vector<float> m_step_sx;      ///< step x (start) position (optional)
    std::vector<float> m_step_sy;      ///< step y (start) position (optional)
    std::vector<float> m_step_sz;      ///< step z (start) position (optional)
    std::vector<float> m_step_x;       ///< step x position
    std::vector<float> m_step_y;       ///< step y position
    std::vector<float> m_step_z;       ///< step z position
    std::vector<float> m_step_ex;      ///< step x (end) position (optional)
    std::vector<float> m_step_ey;      ///< step y (end) position (optional)
    std::vector<float> m_step_ez;      ///< step z (end) position (optional)
    std::vector<float> m_step_dx;      ///< step x direction
    std::vector<float> m_step_dy;      ///< step y direction
    std::vector<float> m_step_dz;      ///< step z direction
    std::vector<float> m_step_length;  ///< step length
    std::vector<float> m_step_X0;      ///< step material x0
    std::vector<float> m_step_L0;      ///< step material l0
    std::vector<float> m_step_A;       ///< step material A
    std::vector<float> m_step_Z;       ///< step material Z
    std::vector<float> m_step_rho;     ///< step material rho

    std::vector<std::uint64_t>
        m_sur_id;  ///< ID of the suface associated with the step
    std::vector<int32_t>
        m_sur_type;              ///< Type of the suface associated with the step
    std::vector<float> m_sur_x;  ///< x position of the center of the suface
                                ///< associated with the step
    std::vector<float> m_sur_y;  ///< y position of the center of the suface
                                ///< associated with the step
    std::vector<float> m_sur_z;  ///< z position of the center of the suface
                                ///< associated with the step
    std::vector<float>
        m_sur_pathCorrection;  ///< path correction when associating
                                ///< material to the given surface
    std::vector<float>
        m_sur_range_min;  ///< Min range of the suface associated with the step
    std::vector<float>
        m_sur_range_max;  ///< Max range of the suface associated with the step

    std::vector<std::uint64_t>
        m_vol_id;  ///< ID of the volume associated with the step

    // Setup Root I/O
    TFile *m_outputFile = new TFile("MaterialOutput.root", "Recreate");
    TTree *m_outputTree = new TTree("material-tracks", "TTree from MaterialTrackWriter");

    // Set the branches
    m_outputTree->Branch("event_id", &m_eventId);
    m_outputTree->Branch("v_x", &m_v_x);
    m_outputTree->Branch("v_y", &m_v_y);
    m_outputTree->Branch("v_z", &m_v_z);
    m_outputTree->Branch("v_px", &m_v_px);
    m_outputTree->Branch("v_py", &m_v_py);
    m_outputTree->Branch("v_pz", &m_v_pz);
    m_outputTree->Branch("v_phi", &m_v_phi);
    m_outputTree->Branch("v_eta", &m_v_eta);
    m_outputTree->Branch("t_X0", &m_tX0);
    m_outputTree->Branch("t_L0", &m_tL0);
    m_outputTree->Branch("mat_x", &m_step_x);
    m_outputTree->Branch("mat_y", &m_step_y);
    m_outputTree->Branch("mat_z", &m_step_z);
    m_outputTree->Branch("mat_dx", &m_step_dx);
    m_outputTree->Branch("mat_dy", &m_step_dy);
    m_outputTree->Branch("mat_dz", &m_step_dz);
    m_outputTree->Branch("mat_step_length", &m_step_length);
    m_outputTree->Branch("mat_X0", &m_step_X0);
    m_outputTree->Branch("mat_L0", &m_step_L0);
    m_outputTree->Branch("mat_A", &m_step_A);
    m_outputTree->Branch("mat_Z", &m_step_Z);
    m_outputTree->Branch("mat_rho", &m_step_rho);

    if (prePostStep) {
        m_outputTree->Branch("mat_sx", &m_step_sx);
        m_outputTree->Branch("mat_sy", &m_step_sy);
        m_outputTree->Branch("mat_sz", &m_step_sz);
        m_outputTree->Branch("mat_ex", &m_step_ex);
        m_outputTree->Branch("mat_ey", &m_step_ey);
        m_outputTree->Branch("mat_ez", &m_step_ez);
    }
    if (storeSurface) {
        m_outputTree->Branch("sur_id", &m_sur_id);
        m_outputTree->Branch("sur_type", &m_sur_type);
        m_outputTree->Branch("sur_x", &m_sur_x);
        m_outputTree->Branch("sur_y", &m_sur_y);
        m_outputTree->Branch("sur_z", &m_sur_z);
        m_outputTree->Branch("sur_pathCorrection", &m_sur_pathCorrection);
        m_outputTree->Branch("sur_range_min", &m_sur_range_min);
        m_outputTree->Branch("sur_range_max", &m_sur_range_max);
    }
    if (storeVolume) {
        m_outputTree->Branch("vol_id", &m_vol_id);
    }

    // Loop over the material tracks and write them out
    for (auto& [idTrack, mtrack] : materialTracks) {
        std::cout << "Material interaction in track" << std::endl;
        // Clearing the vector first
        m_step_sx.clear();
        m_step_sy.clear();
        m_step_sz.clear();
        m_step_x.clear();
        m_step_y.clear();
        m_step_z.clear();
        m_step_ex.clear();
        m_step_ey.clear();
        m_step_ez.clear();
        m_step_dx.clear();
        m_step_dy.clear();
        m_step_dz.clear();
        m_step_length.clear();
        m_step_X0.clear();
        m_step_L0.clear();
        m_step_A.clear();
        m_step_Z.clear();
        m_step_rho.clear();

        m_sur_id.clear();
        m_sur_type.clear();
        m_sur_x.clear();
        m_sur_y.clear();
        m_sur_z.clear();
        m_sur_pathCorrection.clear();
        m_sur_range_min.clear();
        m_sur_range_max.clear();

        m_vol_id.clear();

        // Reserve the vector then
        size_t mints = mtrack.second.materialInteractions.size();
        m_step_sx.reserve(mints);
        m_step_sy.reserve(mints);
        m_step_sz.reserve(mints);
        m_step_x.reserve(mints);
        m_step_y.reserve(mints);
        m_step_z.reserve(mints);
        m_step_ex.reserve(mints);
        m_step_ey.reserve(mints);
        m_step_ez.reserve(mints);
        m_step_dx.reserve(mints);
        m_step_dy.reserve(mints);
        m_step_dz.reserve(mints);
        m_step_length.reserve(mints);
        m_step_X0.reserve(mints);
        m_step_L0.reserve(mints);
        m_step_A.reserve(mints);
        m_step_Z.reserve(mints);
        m_step_rho.reserve(mints);

        m_sur_id.reserve(mints);
        m_sur_type.reserve(mints);
        m_sur_x.reserve(mints);
        m_sur_y.reserve(mints);
        m_sur_z.reserve(mints);
        m_sur_pathCorrection.reserve(mints);
        m_sur_range_min.reserve(mints);
        m_sur_range_max.reserve(mints);

        m_vol_id.reserve(mints);

        // reset the global counter
        if (recalculateTotals) {
            m_tX0 = 0.;
            m_tL0 = 0.;
        } else {
            m_tX0 = mtrack.second.materialInX0;
            m_tL0 = mtrack.second.materialInL0;
        }

        // set the track information at vertex
        m_v_x = mtrack.first.first.x();
        m_v_y = mtrack.first.first.y();
        m_v_z = mtrack.first.first.z();
        m_v_px = mtrack.first.second.x();
        m_v_py = mtrack.first.second.y();
        m_v_pz = mtrack.first.second.z();
        m_v_phi = Acts::VectorHelpers::phi(mtrack.first.second);
        m_v_eta = Acts::VectorHelpers::eta(mtrack.first.second);

        // an now loop over the material
        for (auto& mint : mtrack.second.materialInteractions) {
            auto direction = mint.direction.normalized();
            
            // Save it somewhere so that it can be compared and added to the output file
            std::cout << mint.position.x() << " " << mint.position.y() << " " 
                << mint.position.z() << " " << mint.time << std::endl;

            // The material step position information
            m_step_x.push_back(mint.position.x());
            m_step_y.push_back(mint.position.y());
            m_step_z.push_back(mint.position.z());
            m_step_dx.push_back(direction.x());
            m_step_dy.push_back(direction.y());
            m_step_dz.push_back(direction.z());

            if (prePostStep) {
                Acts::Vector3 prePos =
                    mint.position - 0.5 * mint.pathCorrection * direction;
                Acts::Vector3 posPos =
                    mint.position + 0.5 * mint.pathCorrection * direction;

                m_step_sx.push_back(prePos.x());
                m_step_sy.push_back(prePos.y());
                m_step_sz.push_back(prePos.z());
                m_step_ex.push_back(posPos.x());
                m_step_ey.push_back(posPos.y());
                m_step_ez.push_back(posPos.z());
            }

            // Store surface information
            if (storeSurface) {
                const Acts::Surface* surface = mint.surface;
                if (mint.intersectionID.value() != 0) {
                    m_sur_id.push_back(mint.intersectionID.value());
                    m_sur_pathCorrection.push_back(mint.pathCorrection);
                    m_sur_x.push_back(mint.intersection.x());
                    m_sur_y.push_back(mint.intersection.y());
                    m_sur_z.push_back(mint.intersection.z());
                } else if (surface != nullptr) {
                    auto sfIntersection = surface->intersect(
                        gctx, mint.position, mint.direction, true);
                    m_sur_id.push_back(surface->geometryId().value());
                    m_sur_pathCorrection.push_back(1.0);
                    m_sur_x.push_back(sfIntersection.intersection.position.x());
                    m_sur_y.push_back(sfIntersection.intersection.position.y());
                    m_sur_z.push_back(sfIntersection.intersection.position.z());
                } else {
                    m_sur_id.push_back(Acts::GeometryIdentifier().value());
                    m_sur_x.push_back(0);
                    m_sur_y.push_back(0);
                    m_sur_z.push_back(0);
                    m_sur_pathCorrection.push_back(1.0);
                }
                if (surface != nullptr) {
                    m_sur_type.push_back(surface->type());
                    const Acts::SurfaceBounds& surfaceBounds = surface->bounds();
                    const Acts::RadialBounds* radialBounds =
                        dynamic_cast<const Acts::RadialBounds*>(&surfaceBounds);
                    const Acts::CylinderBounds* cylinderBounds =
                        dynamic_cast<const Acts::CylinderBounds*>(&surfaceBounds);
                    if (radialBounds != nullptr) {
                        m_sur_range_min.push_back(radialBounds->rMin());
                        m_sur_range_max.push_back(radialBounds->rMax());
                    } else if (cylinderBounds != nullptr) {
                        m_sur_range_min.push_back(
                            -cylinderBounds->get(Acts::CylinderBounds::eHalfLengthZ));
                        m_sur_range_max.push_back(
                            cylinderBounds->get(Acts::CylinderBounds::eHalfLengthZ));
                    } else {
                        m_sur_range_min.push_back(0);
                        m_sur_range_max.push_back(0);
                    }
                } else {
                    m_sur_type.push_back(-1);
                    m_sur_range_min.push_back(0);
                    m_sur_range_max.push_back(0);
                }
            }

            // store volume information
            if (storeVolume) {
                const Acts::Volume* volume = mint.volume;
                Acts::GeometryIdentifier vlayerID;
                if (volume != nullptr) {
                    vlayerID = volume->geometryId();
                    m_vol_id.push_back(vlayerID.value());
                } else {
                    vlayerID.setVolume(0);
                    vlayerID.setBoundary(0);
                    vlayerID.setLayer(0);
                    vlayerID.setApproach(0);
                    vlayerID.setSensitive(0);
                    m_vol_id.push_back(vlayerID.value());
                }
            }

            // the material information
            const auto& mprops = mint.materialSlab;
            m_step_length.push_back(mprops.thickness());
            m_step_X0.push_back(mprops.material().X0());
            m_step_L0.push_back(mprops.material().L0());
            m_step_A.push_back(mprops.material().Ar());
            m_step_Z.push_back(mprops.material().Z());
            m_step_rho.push_back(mprops.material().massDensity());
            // re-calculate if defined to do so
            if (recalculateTotals) {
                m_tX0 += mprops.thicknessInX0();
                m_tL0 += mprops.thicknessInL0();
            }
        }
    // write to
    m_outputTree->Fill();
  }
  m_outputTree->Write();
  m_outputTree->Print();
  m_outputFile->Close();
}
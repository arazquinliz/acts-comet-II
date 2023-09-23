// ---------------------------------------------------------------- //
//                Creator: Amaia Razquin Lizarraga                  //
//                Contact: amaiafisika@gmail.com                    //
//                Date:    7 September 2023                         //
// ---------------------------------------------------------------- //
/* 
    This program replicates the selection processes of the seeding algorithm
    in order to get the values of the selection criterias (cuts) for all the
    triplets in a trajectory. A plot for each selected cut is then printed out
    in order to tune the cuts. A better fine-tunning step should be done with 
    a hyperparameter optimiser, this program only offers a first step tuning.
    All plots are saved in a histogram folder. Note it runs under the assumption
    that the maximum and minimum r values are on the edges between middle and 
    top/bottom SP means between to points of a duplet.
    Command:
        g++ seedSelection.cpp Geometry.cpp -o selection.x
        -L/path/to/acts_install/lib -lActsPluginTGeo -lActsPluginJson -lActsCore 
        -I/path/to/acts_install/include/ -I$(root-config --incdir) 
        -I$(root-config --evelibs) $(root-config --libs) 
        -lGeom -I/usr/include/eigen3
*/

// Std
#include <iostream>
#include <vector>

// Acts
#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/Seeding/BinnedSPGroup.hpp>
#include <Acts/Seeding/InternalSeed.hpp>
#include <Acts/Seeding/InternalSpacePoint.hpp>
#include <Acts/Seeding/Seed.hpp>
#include <Acts/Seeding/SeedFinder.hpp>
#include <Acts/Seeding/SeedFilter.hpp>

// ROOT
#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <Math/Vector3D.h>
#include <TH1.h>
#include <TCanvas.h>

// Amaia
#include "Geometry.hpp"

using namespace Acts::UnitLiterals;

int main() {

    // get the tree and iterate over the branches
	auto *f = TFile::Open("detector_oa_all_td.root");
	if (!f) {
		std::cout << "Couldn't open measurement file\n";
	}

    double radius    = 0.49_cm;

	auto *tree = (TTree*)f->Get("detections");
	std::vector<ROOT::Math::XYZVector> *spOfEvent = 0;
    std::vector<ROOT::Math::XYZVector> *momOfEvent = 0;
	tree->SetBranchAddress("spOfEvent", &spOfEvent);
    tree->SetBranchAddress("momOfEvent", &momOfEvent);

    std::cout << "Number of entries: " << tree->GetEntries() << std::endl;

    // Initialize all the histograms
    TH1F *h_deltaRMidMax = new TH1F("deltaRMidMax", "deltaRMiddleMaxSPRange", 100, 0, 15); // Min
    TH1F *h_deltaRMidMin = new TH1F("deltaRMidMin", "deltaRMiddleMinSPRange", 100, 0, 20); // Min
    TH1F *h_deltaRMaxB   = new TH1F("deltaRMaxB", "deltaRMaxBottomSP", 100, 0, 20); // Max
    TH1F *h_deltaRMinB   = new TH1F("deltaRMinB", "deltaRMinBottomSP", 100, 0, 20); // Min
    TH1F *h_deltaRMaxT   = new TH1F("deltaRMaxT", "deltaRMaxTopSP", 100, 0, 20); // Max
    TH1F *h_deltaRMinT   = new TH1F("deltaRMinT", "deltaRMinTopSP", 100, 0, 15); // Min
    TH1F *h_deltaZMax    = new TH1F("deltaZMax", "deltaZMax", 100, 5, 30); // Max
    TH1F *h_cotThetaMax  = new TH1F("cotThetaMax", "cotThetaMax", 100, 0, 10); // Max
    TH1F *h_collMax      = new TH1F("collMax", "collisionRegionMax", 10000, -4000, 1400); // Max
    TH1F *h_collMin      = new TH1F("collMin", "collisionRegionMin", 10000, -4000, 1400); // Min
    TH1F *h_impactMax    = new TH1F("impactMax", "impactMax", 1000, 0, 1500); // Max
    TH1F *h_cotThetaT    = new TH1F("cotThetaT", "cotThetaT", 100, 0, 10); // Comparative (unused) 
    TH1F *h_cotThetaB    = new TH1F("cotThetaB", "cotThetaB", 100, 0, 10); // Comparative (unused) 
    TH1F *h_minHelix     = new TH1F("minHelix", "minHelixRadius", 1000, 0, 600); // Min (unused)
    TH1F *h_deltaCot     = new TH1F("deltaCot", "deltaCotTheta2", 100, 0, 0.6); // Max (unused) error2 + scatteringInRegion2
    TH1F *h_minPt        = new TH1F("minPt", "minPt", 100, 0, 105); // Min 
    TH1F *h_zOrigin      = new TH1F("zOrigin", "zOrigin", 4000, -10000, -1000);

    std::vector<int> unseeded = {15, 25, 31, 48, 53, 75, 84, 90, 97, 111, 122, 133, 134, 178, 183, 184, 190, 215, 229, 238, 256, 276, 279, 359, 375, 394, 400, 411, 413, 417, 419, 425, 445, 450, 518, 522, 526, 543, 581, 592, 610, 612, 645, 663, 666, 692, 727, 751, 853, 888, 899, 910, 917, 939, 960, 977, 988, 994, 1015, 1037, 1059, 1064, 1069, 1104, 1112, 1118, 1125, 1128, 1131, 1132, 1133, 1144, 1166, 1167, 1180, 1218, 1222, 1233, 1239, 1244, 1253, 1269, 1270, 1275, 1278, 1287, 1293, 1311, 1336, 1369, 1387, 1395, 1401, 1411, 1417, 1435, 1436, 1459, 1473, 1513, 1524, 1542, 1546, 1563, 1567, 1589, 1602, 1604, 1613, 1617, 1626, 1627, 1628, 1633, 1667, 1672, 1680, 1684, 1694, 1727, 1734, 1745, 1747, 1772, 1781, 1784, 1789, 1805, 1834, 1835, 1841, 1870, 1887, 1899, 1908, 1921, 1926, 1928, 1963, 1969, 1973, 1975, 1982, 2021, 2022, 2026, 2032, 2035, 2067, 2085, 2086, 2126, 2149, 2151, 2162, 2232, 2239, 2293, 2331, 2349, 2350, 2357, 2362, 2372, 2391, 2397, 2402, 2437, 2445, 2447, 2459, 2495, 2499, 2505, 2506, 2518, 2532, 2534, 2537, 2539, 2553, 2562, 2573, 2577, 2594, 2601, 2611, 2613, 2637, 2657, 2658, 2661, 2718, 2725, 2732, 2785, 2789, 2791, 2817, 2832, 2837, 2844, 2853, 2890, 2903, 2906, 2922, 2928, 2931, 2942, 2949, 2952, 2987, 3002, 3017, 3037, 3040, 3048, 3088, 3091, 3119, 3126, 3128, 3181, 3191, 3211, 3260, 3281, 3283, 3293, 3314, 3322, 3341, 3346, 3351, 3353, 3358, 3360, 3363, 3365, 3378, 3379, 3380, 3384, 3386, 3399, 3409, 3439, 3444, 3476, 3485, 3504, 3506, 3514, 3542, 3586, 3641, 3646, 3650, 3656, 3664, 3665, 3668, 3669, 3679, 3782, 3787, 3790, 3803, 3807, 3819, 3866, 3869, 3879, 3889, 3894, 3900, 3904, 3906, 3911, 3947, 3950, 3954, 3972, 4018, 4027, 4034, 4043, 4046, 4056, 4062, 4075, 4099, 4121, 4136, 4141, 4156, 4180, 4185, 4208, 4211, 4240, 4243, 4245, 4255, 4271, 4314, 4330, 4343, 4357, 4358, 4360, 4394, 4396, 4402, 4416, 4420, 4434, 4442, 4488, 4505, 4506, 4510, 4517, 4528, 4533, 4544, 4578, 4611, 4615, 4661, 4675, 4693, 4733, 4752, 4763, 4770, 4772, 4782, 4799, 4841, 4862, 4865, 4870, 4871, 4872, 4878, 4887, 4893, 4900, 4934, 4953, 4956, 4960, 4961, 4981, 5003, 5017, 5032, 5041, 5051, 5064, 5075, 5080, 5116, 5125, 5142, 5148, 5159, 5160, 5173, 5177, 5191, 5244, 5249, 5262, 5291, 5294, 5306, 5331, 5346, 5349, 5370, 5431, 5438, 5442, 5445, 5456, 5466, 5486, 5489, 5496, 5498, 5507, 5516, 5564, 5565, 5567, 5596, 5598, 5630, 5647, 5651, 5654, 5665, 5677, 5687, 5694, 5722, 5728, 5734, 5775, 5780, 5782, 5784, 5801, 5817, 5820, 5872, 5876, 5895, 5896, 5900, 5928, 5944, 5955, 5960, 5961, 5964, 5969, 5978, 5990, 5992, 5997, 6002, 6008, 6036, 6037, 6040, 6048, 6062, 6065, 6068, 6116, 6120, 6121, 6123, 6144, 6154, 6158, 6160, 6169, 6171, 6198, 6216, 6222, 6231, 6234, 6243, 6246, 6256, 6268, 6290, 6296, 6299, 6311, 6337, 6353, 6357, 6360, 6376, 6377, 6383, 6428, 6430, 6459, 6464, 6506, 6515, 6557, 6576, 6620, 6627, 6638, 6651, 6654, 6670, 6690, 6694, 6712, 6724, 6736, 6741, 6750, 6768, 6777, 6811, 6837, 6841, 6917, 6928, 6947, 6957, 6970, 7027, 7033, 7044, 7045, 7046, 7049, 7057, 7063, 7064, 7082, 7100, 7110, 7125, 7131, 7138, 7142, 7180, 7184, 7205, 7216, 7235, 7256, 7266, 7272, 7295, 7311, 7342, 7360, 7389, 7413, 7416, 7429, 7451, 7457, 7464, 7475, 7478, 7481, 7484, 7485};

    //std::vector<int> unseeded = {50, 53, 69, 75, 143, 175, 184, 300, 315, 332, 336, 479, 628, 722, 732, 773, 882, 892, 899, 962, 981, 1018, 1054, 1141, 1156, 1218, 1257, 1312, 1395, 1411, 1435, 1436, 1450, 1543, 1587, 1604, 1627, 1667, 1687, 1728, 1799, 1871, 1872, 1888, 1928, 1957, 2030, 2085, 2211, 2244, 2293, 2340, 2367, 2583, 2598, 2608, 2637, 2656, 2657, 2661, 2677, 2736, 2757, 2810, 2813, 2873, 2880, 2945, 3008, 3146, 3268, 3278, 3353, 3358, 3360, 3373, 3409, 3542, 3546, 3671, 3702, 3745, 3807, 3832, 3904, 3911, 3922, 4097, 4105, 4145, 4159, 4167, 4208, 4267, 4314, 4360, 4370, 4488, 4509, 4518, 4619, 4722, 4732, 4757, 4758, 4763, 4766, 4836, 4840, 4872, 4884, 4893, 4953, 4956, 4970, 4981, 4993, 5008, 5021, 5042, 5100, 5254, 5262, 5291, 5305, 5395, 5432, 5440, 5471, 5476, 5489, 5498, 5516, 5526, 5536, 5575, 5598, 5640, 5670, 5728, 5827, 5895, 5918, 5997, 5998, 6002, 6026, 6062, 6065, 6116, 6140, 6171, 6198, 6231, 6237, 6309, 6311, 6357, 6388, 6419, 6469, 6566, 6578, 6590, 6638, 6648, 6654, 6685, 6691, 6716, 6736, 6739, 6743, 6750, 6841, 6885, 6895, 7004, 7006, 7045, 7049, 7159, 7216, 7239, 7260, 7277, 7291, 7314, 7344};
    //std::vector<int> unseeded = {1921, 3017, 4314};

    // Get limits directly
    double dRMidMax     = 1400.; // Min
    double dRMidMin     = 1400.; // Min
    double dRMaxB       = -1400.; // Max
    double dRMinB       = 1400.; // Min
    double dRMaxT       = -1400.; // Max
    double dRMinT       = 1400.; // Min
    double dZMax        = -1400.; // Max
    double ctThetaMax   = -1400.; // Max
    double collisionMax = -1400000.; // Max
    double collisionMin = 1400.; // Min
    double impactMax    = -1400.; // Max
    //double cotThetaT    = 0.; // Comparative (unused) 
    //double cotThetaB    = 0.; // Comparative (unused) 
    //double minHelix     = 0.; // Min (unused)
    //double deltaCot     = 0.; // Max (unused) error2 + scatteringInRegion2
    double mnPt        = 1400.; // Min   

    std::fstream spacepointfile;
	spacepointfile.open("spacepoints.txt", std::ios::out); 

    for (size_t ev = 0; ev < tree->GetEntries(); ++ev) {
    //for (size_t ev = 0; ev < 10; ++ev) {
    //for (auto ev : unseeded) {
        std::vector<Acts::Vector3> posVec = {Acts::Vector3::Zero(), Acts::Vector3::Zero(), 
            Acts::Vector3::Zero(), Acts::Vector3::Zero()};
        std::vector<int> k = {0, 0, 0, 0};
        tree->GetEntry(ev);
        
        for (size_t it = 0; it < spOfEvent->size(); ++it) {
            double positionZ = (*spOfEvent)[it].Z() - 333.45_cm;
            if (positionZ >= -1200) {
                continue;}

            if (positionZ <= -1380.25 + radius) {
                posVec[0].x() += (*spOfEvent)[it].X();
                posVec[0].y() += (*spOfEvent)[it].Y();
                posVec[0].z() += (*spOfEvent)[it].Z();
                k[0] += 1;
            }
            else if (positionZ <= -1370.9 + radius) {
                posVec[1].x() += (*spOfEvent)[it].X();
                posVec[1].y() += (*spOfEvent)[it].Y();
                posVec[1].z() += (*spOfEvent)[it].Z();
                k[1] += 1;
            }
            else if (positionZ <= -1359.1 + radius) {
                posVec[2].x() += (*spOfEvent)[it].X();
                posVec[2].y() += (*spOfEvent)[it].Y();
                posVec[2].z() += (*spOfEvent)[it].Z();
                k[2] += 1;
            }
            else if (positionZ <= -1349.75 + radius) {
                posVec[3].x() += (*spOfEvent)[it].X();
                posVec[3].y() += (*spOfEvent)[it].Y();
                posVec[3].z() += (*spOfEvent)[it].Z();
                k[3] += 1;
            }
        }
        
        int kount = 0;
        std::vector<Acts::Vector3> spacePoints;
        for (auto pos : posVec) {
            if (k[kount] > 0) { // Only save if there are measurements in layer
                spacePoints.push_back(Acts::Vector3{pos.x()/k[kount] - 1080._cm, 
                    pos.y()/k[kount], pos.z()/k[kount] - 333.45_cm});
                spacepointfile << pos.x()/k[kount] - 1080._cm << " " <<
                    pos.y()/k[kount] << " " << pos.z()/k[kount] - 333.45_cm << "\n";
            }
            kount++;
        } // spacepoints

        if (spacePoints.size() < 3) {continue;} // Cannot find triplets, too few measurements
        
        // Middle SpacePoint cuts
        // Get variableMiddelSPRange
        double rMax = 0., rMin = 600.;
        double zMax = -1400., zMin = 0.;
        for (auto spM : spacePoints) {
            float rM = sqrt(spM.x() * spM.x() + spM.y() * spM.y());
            float zM = spM.z(); 
            if (rM > rMax) {rMax = std::floor(rM * 0.5) * 2;}
            if (rM < rMin) {rMin = std::floor(rM * 0.5) * 2;}
            if (zM > zMax) {zMax = zM;}
            if (zM < zMin) {zMin = zM;}
        }
        double deltaRMidMin = 600;
        double deltaRMidMax = 600;
        for (int i = 1; i < spacePoints.size() - 1; i++) {
            float rM = sqrt(spacePoints[i].x() * spacePoints[i].x() + 
                spacePoints[i].y() * spacePoints[i].y());
            if ((rMax - rM) < deltaRMidMax) {deltaRMidMax = rMax - rM;}
            if ((rM - rMin) < deltaRMidMin) {deltaRMidMin = rM - rMin;}
        } // rM have to be WITHIN (rMin + deltaRMin, rMax - deltaRMax)

        if (deltaRMidMax < dRMidMax) {dRMidMax = deltaRMidMax;}
        if (deltaRMidMin < dRMidMin) {dRMidMin = deltaRMidMin;}
        // Compatible doublets (top/bottom) cuts
        double deltaRMinT  = 600;
        double deltaRMaxT  = 0;
        double deltaRMinB  = 600;
        double deltaRMaxB  = 0;
        double deltaZMax   = 0;
        double cotThetaMax = 0;
        double collMin     = 1400;
        double collMax     = -4000;
        double zOrigin     = 0;
        for (int i = 1; i < spacePoints.size() - 1; i++) {
            float rM = sqrt(spacePoints[i].x() * spacePoints[i].x() + 
                spacePoints[i].y() * spacePoints[i].y());
            float zM = spacePoints[i].z(); 

            for (int j = 0; j < spacePoints.size(); j++) {
                if (i == j) {continue;}
                float r0 = sqrt(spacePoints[j].x() * spacePoints[j].x() + 
                    spacePoints[j].y() * spacePoints[j].y());
                float z0 = spacePoints[j].z(); 

                // Get min/max r distance between middle and top/bottom SP
                float deltaR = r0 - rM;
                if (deltaR > 0) {
                    if (deltaR > deltaRMaxT) {deltaRMaxT = deltaR;}
                    if (deltaR < deltaRMinT) {deltaRMinT = deltaR;}
                }
                else if (deltaR <= 0) {
                    if (abs(deltaR) > deltaRMaxB) {deltaRMaxB = abs(deltaR);}
                    if (abs(deltaR) < deltaRMinB) {deltaRMinB = abs(deltaR);}
                }

                // Get  min/max z distance between middle and top/bottom SP 
                // Only for consistency --> It is already limited by the station
                float deltaZ = abs(z0 - zM);
                if (deltaZ > deltaZMax) {deltaZMax = deltaZ;}

                // Get cotTheta (forward angle, z/R ratio) between middle and top/bottom SP
                float cotTheta = deltaZ / abs(deltaR);
                if (cotTheta > cotThetaMax) {cotThetaMax = cotTheta;} 

                // Get douplet origin for collision region
                zOrigin = zM - rM * cotTheta; // - to account for T/B dynamics
                if (zOrigin < collMin) {collMin = zOrigin;}
                if (zOrigin > collMax) {collMax = zOrigin;}
            }   
        }
        if ((deltaRMaxB > 30) | (deltaRMaxT > 30)) {continue;}
        if (deltaRMaxB > dRMaxB) {dRMaxB = deltaRMaxB;}
        if (deltaRMinB < dRMinB) {dRMinB = deltaRMinB;}
        if (deltaRMaxT > dRMaxT) {dRMaxT = deltaRMaxT;}
        if (deltaRMinT < dRMinT) {dRMinT = deltaRMinT;}
        if (deltaZMax > dZMax) {dZMax = deltaZMax;}
        //if (cotThetaMax > 10) {std::cout << ev << ", "; continue;}
        if (cotThetaMax > ctThetaMax) {ctThetaMax = cotThetaMax;}
        if (collMax > collisionMax) {collisionMax = collMax;}
        if (collMin < collisionMin) {collisionMin = collMin;}

        // Candidate filter cuts
        float cotThetaT, cotThetaB;
        float minHelixRadius;
        float minPt;
        float deltaCotTheta2;

        double pT = std::sqrt((*momOfEvent)[0].X() * (*momOfEvent)[0].X() + 
            (*momOfEvent)[0].Y() * (*momOfEvent)[0].Y());

        float Im = 0;
        for (int i = 1; i < spacePoints.size() - 1; i++) {
            float rM = sqrt(spacePoints[i].x() * spacePoints[i].x() + 
                spacePoints[i].y() * spacePoints[i].y());
            float zM = spacePoints[i].z();

            double r0 = sqrt(spacePoints[0].x() * spacePoints[0].x() + 
                spacePoints[0].y() * spacePoints[0].y());
            double r3 = sqrt(spacePoints[3].x() * spacePoints[3].x() + 
                spacePoints[3].y() * spacePoints[3].y());

            Acts::Vector3 top, bottom;
            if (r0 > r3) {
                top    = spacePoints[0];
                bottom = spacePoints[3];
            }
            else if (r0 < r3) {
                top    = spacePoints[3];
                bottom = spacePoints[0];
            }

            float rT = std::sqrt(top.x() * top.x() + top.y() * top.y());
            float zT = top.z();
            float rB = std::sqrt(bottom.x() * bottom.x() + bottom.y() * bottom.y());
            float zB = bottom.z();

            // conformal transformation u=x/(x²+y²) v=y/(x²+y²)
            float Vb = bottom.y() / (bottom.x() * bottom.x() + bottom.y() * bottom.y());
            float Ub = bottom.x() / (bottom.x() * bottom.x() + bottom.y() * bottom.y());
            float Vt = top.y() / (top.x() * top.x() + top.y() * top.y());
            float Ut = top.x() / (top.x() * top.x() + top.y() * top.y());
            float dU = Ut - Ub;
            float A  = (Vt - Vb) / dU;
            float B  = Vb - A * Ub;
            float S2 = 1. + A * A;
            float B2 = B * B;

            cotThetaT = abs(zT - zM) / abs(rT - rM);
            cotThetaB = abs(zB - zM) / abs(rB - rM);

            float minHelixDiameter2 = S2 / B2; // Has to be bigger (in options)
            minHelixRadius = std::sqrt(minHelixDiameter2) * 0.5;
            
            // Compatibility between the slopes of the two seed segments (not done yet)
            
            float deltaCotTheta = cotThetaB - cotThetaT;
            deltaCotTheta2 = deltaCotTheta * deltaCotTheta; 

            minPt = std::sqrt(S2 / B2) * 150 * 0.001; // for the proper calc of pT
            // Impact parameter
            Im = abs((A - B * rM) * rM);
        }
        //std::cout << minPt << std::endl;
        if (Im > impactMax) {impactMax = Im;}
        if (minPt < mnPt) {mnPt = minPt;}

        // Fill histograms only for events with enough hits
        if (spacePoints.size() > 2) {
            //std::cout << deltaRMaxB << " " << deltaRMinB << " " << deltaRMaxT << " " << deltaRMinT << std::endl;
            //std::cout << cotThetaMax << " " << collMax << " " << collMin << std::endl;
            //std::cout << Im << std::endl;
            h_deltaRMidMax->Fill(deltaRMidMax);
            h_deltaRMidMin->Fill(deltaRMidMin);
            h_deltaRMaxB->Fill(deltaRMaxB);
            h_deltaRMinB->Fill(deltaRMinB);
            h_deltaRMaxT->Fill(deltaRMaxT);
            h_deltaRMinT->Fill(deltaRMinT);
            h_deltaZMax->Fill(deltaZMax);
            h_cotThetaMax->Fill(cotThetaMax);
            h_collMax->Fill(collMax);
            h_collMin->Fill(collMin);
            h_impactMax->Fill(Im);
            h_cotThetaT->Fill(cotThetaT);
            h_cotThetaB->Fill(cotThetaB);
            h_minHelix->Fill(minHelixRadius);
            h_deltaCot->Fill(deltaCotTheta2);
            h_minPt->Fill(minPt);
            h_zOrigin->Fill(zOrigin);
        }
    } // loop over events
    spacepointfile.close(); 
    // Global limits
    std::cout << "deltaRMidMax: " << dRMidMax << " deltaRMidMin: " << dRMidMin << 
        "\ndeltaRMaxB: " << dRMaxB << " deltaRMinB: " << dRMinB << " deltaRMinT: " << dRMinT << " deltaRMaxT: " << dRMaxT <<
        "\ndeltaZMax: " << dZMax << " cotThetaMax: " << ctThetaMax << " collMax: " << collisionMax << " collMin: " << collisionMin << 
        "\nimpactMax: " << impactMax << " minPt: " << mnPt << std::endl;

    // Draw histograms in different files inside the histograms folder
    // In the future just output all of them as ONE single vroot file
    TCanvas *c0 = new TCanvas();
	h_deltaRMidMax->Draw();
    c0->Print("histograms/deltaRMidMax.png");

    TCanvas *c1 = new TCanvas();
	h_deltaRMidMin->Draw();
    c1->Print("histograms/deltaRMidMin.png");

    TCanvas *c2 = new TCanvas();
	h_deltaRMinB->Draw();
    c2->Print("histograms/deltaRMinB.png");
    
    TCanvas *c3 = new TCanvas();
	h_deltaRMaxB->Draw();
    c3->Print("histograms/deltaRMaxB.png");

    TCanvas *c4 = new TCanvas();
	h_deltaRMinT->Draw();
    c4->Print("histograms/deltaRMinT.png");

    TCanvas *c5 = new TCanvas();
	h_deltaRMaxT->Draw();
    c5->Print("histograms/deltaRMaxT.png");

    TCanvas *c6 = new TCanvas();
    h_deltaZMax->Draw();
    c6->Print("histograms/deltaZMax.png");

    TCanvas *c7 = new TCanvas();
    h_cotThetaMax->Draw();
    c7->Print("histograms/cotThetaMax.png");

    TCanvas *c8 = new TCanvas();
    h_collMax->Draw();
    c8->Print("histograms/collMax.png");

    TCanvas *c9 = new TCanvas();
    h_collMin->Draw();
    c9->Print("histograms/collMin.png");

    TCanvas *c10 = new TCanvas();
    h_impactMax->Draw();
    c10->Print("histograms/impactMax.png");

    TCanvas *c11 = new TCanvas();
    h_cotThetaT->Draw();
    c11->Print("histograms/cotThetaT.png");

    TCanvas *c12 = new TCanvas();
    h_cotThetaB->Draw();
    c12->Print("histograms/cotThetaB.png");

    TCanvas *c13 = new TCanvas();
    h_minHelix->Draw();
    c13->Print("histograms/minHelix.png");

    TCanvas *c14 = new TCanvas();
    h_deltaCot->Draw();
    c14->Print("histograms/deltaCot.png");

    TCanvas *c15 = new TCanvas();
    h_minPt->Draw();
    c15->Print("histograms/minPt.png");

    TCanvas *c16 = new TCanvas();
    h_zOrigin->Draw();
    c16->Print("histograms/zOrigin.png");

    return 0;
}

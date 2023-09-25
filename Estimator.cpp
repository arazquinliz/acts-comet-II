// ---------------------------------------------------------------- //
//                Creator: Amaia Razquin Lizarraga                  //
//                Contact: amaiafisika@gmail.com                    //
//                Date:    7 September 2023                         //
// ---------------------------------------------------------------- //
/* This module contains three circle fits: Levenberg-Marquardt (geometrical),
   Taubin (analytical), Hyper fit (analytical). The three are then used within
   estimateModified. This function is an update of the track parameter estimator
   in ACTS, here the new circle fits can be used and also a line-search minimizer
   is applied. More details can be found in the comments. The original algorithm
   for the circle fit functions are by Nicolai Chernov (2012), and can be found at: 
   https://people.cas.uab.edu/~mosya/cl/CPPcircle.html */

#include <fstream>

#include <Acts/Definitions/TrackParametrization.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/Seeding/Seed.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Utilities/Logger.hpp>

// ROOT
#include <TVector3.h>
#include <Math/Vector3D.h> 
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <Math/RootFinder.h>
#include <TF1.h>

#include "Estimator.hpp"

using namespace Acts::UnitLiterals;
std::fstream estimatorFile;

// Estimate sigma = square root of data divided by N
double sigmaCalc(std::vector<TVector3> points, Circle Circ) {
    double sum = 0;
    int size = points.size();
    for (int i = 0; i < points.size(); i++) {
        double dx  = points[i].X() - Circ.x();
        double dy  = points[i].Y() - Circ.y();
        sum += std::pow(std::sqrt(dx * dx + dy * dy) - Circ.r(), 2);
    }
    return std::sqrt(sum / size);
}

// Hyper fit by Al-Sharadqah and Chernov
Circle fitHyper(std::vector<TVector3> points){
    int iter, iterMax = 1e4;
    double Xi, Yi, Zi;
    double Mz, Mxy, Mxx, Myy, Mxz, Myz, Mzz, Cov_xy, Var_z;
    double A0, A1, A2, A22;
    double Dy, xnew, x, ynew, y;
    double DET, Xcenter, Ycenter;

    int size = points.size();

    int mnX = 0, mnY = 0;
    for (auto point : points) {
        mnX += point.X();
        mnY += point.Y();
    }
    double meanX = mnX / size;
    double meanY = mnY / size;

    // Compute moments
    Mxx = Myy = Mxy = Mxz = Myz = Mzz = 0.;
    for (auto point : points) {
        Xi = point.X() - meanX;
        Yi = point.Y() - meanY;
        Zi = Xi*Xi + Yi*Yi;
        
        Mxy += Xi*Yi;
        Mxx += Xi*Xi;
        Myy += Yi*Yi;
        Mxz += Xi*Zi;
        Myz += Yi*Zi;
        Mzz += Zi*Zi;
    }
    Mxx /= size;
    Myy /= size;
    Mxy /= size;
    Mxz /= size;
    Myz /= size;
    Mzz /= size;

    // Coefficients of the characteristic polynomial
    Mz = Mxx + Myy;
    Cov_xy = Mxx * Myy - Mxy * Mxy;
    Var_z = Mzz - Mz * Mz;

    A2 = 4. * Cov_xy - 3. * Mz * Mz - Mzz;
    A1 = Var_z * Mz + 4.*Cov_xy*Mz - Mxz * Mxz - Myz * Myz;
    A0 = Mxz * (Mxz * Myy - Myz * Mxy) + Myz * (Myz * Mxx - Mxz * Mxy) - Var_z * Cov_xy;
    A22 = A2 + A2;

    // Find root of polynomial with Newton's method at x = 0
    // it will find the first positive root
    for (x = 0., y = A0, iter = 0; iter < iterMax; iter++)  // usually, 4-6 iterations are enough
    {
        Dy = A1 + x * (A22 + 16. * x * x);
        xnew = x - y / Dy;
        if ((xnew == x)||(!isfinite(xnew))) {break;}
        ynew = A0 + xnew * (A1 + xnew*(A2 + 4. * xnew * xnew));
        if (abs(ynew)>=abs(y))  {break;}
        x = xnew;  
        y = ynew;
    }

    // Parameters of circle
    DET = x * x - x * Mz + Cov_xy;
    Xcenter = (Mxz * (Myy - x) - Myz * Mxy) / DET / 2.;
    Ycenter = (Myz * (Mxx - x) - Mxz * Mxy) / DET / 2.;

    double cx    = Xcenter + meanX;
    double cy    = Ycenter + meanY;
    double rad   = sqrt(Xcenter * Xcenter + Ycenter * Ycenter + Mz - x - x);
    double sigma = sigmaCalc(points, Circle{cx, cy, rad, 0, 0, 0});
    
    return Circle{cx, cy, rad, sigma, 0, 0, iter};
}

// Taubin algebraic circle fit
Circle fitTaubin(std::vector<TVector3> points) {
    int iter, iterMax = 1e4;
    double Xi, Yi, Zi;
    double Mz, Mxy, Mxx, Myy, Mxz, Myz, Mzz, Cov_xy, Var_z;
    double A0, A1, A2, A22, A3, A33;
    double Dy, xnew, x, ynew, y;
    double DET, Xcenter, Ycenter;

    int size = points.size();

    int mnX = 0, mnY = 0;
    for (auto point : points) {
        mnX += point.X();
        mnY += point.Y();
    }
    double meanX = mnX / size;
    double meanY = mnY / size;

    // Compute moments
    Mxx = Myy = Mxy = Mxz = Myz = Mzz = 0.;
    for (auto point : points) {
        Xi = point.X() - meanX;
        Yi = point.Y() - meanY;
        Zi = Xi*Xi + Yi*Yi;
        
        Mxy += Xi*Yi;
        Mxx += Xi*Xi;
        Myy += Yi*Yi;
        Mxz += Xi*Zi;
        Myz += Yi*Zi;
        Mzz += Zi*Zi;
    }
    Mxx /= size;
    Myy /= size;
    Mxy /= size;
    Mxz /= size;
    Myz /= size;
    Mzz /= size;

    // Coefficients of characteristic polynomial
    Mz = Mxx + Myy;
    Cov_xy = Mxx * Myy - Mxy * Mxy;
    Var_z = Mzz - Mz * Mz;
    A3 = 4. * Mz;
    A2 = -3. * Mz * Mz - Mzz;
    A1 = Var_z * Mz + 4. * Cov_xy * Mz - Mxz * Mxz - Myz * Myz;
    A0 = Mxz * (Mxz * Myy - Myz * Mxy) + Myz * (Myz * Mxx - Mxz * Mxy) - Var_z * Cov_xy;
    A22 = A2 + A2;
    A33 = A3 + A3 + A3;

    // Find root of polynomial with Newton's method at x = 0
    // it will find the first positive root
    for (x=0, y=A0, iter=  0; iter < iterMax; iter++)  // usually, 4-6 iterations are enough
    {
    	Dy = A1 + x * (A22 + A33 * x);
        xnew = x - y/Dy;
        if ((xnew == x)||(!isfinite(xnew))) {break;}
        ynew = A0 + xnew * (A1 + xnew * (A2 + xnew * A3));
        if (abs(ynew) >= abs(y)) {break;}
        x = xnew;  
        y = ynew;
    }

    // Parameters of the circle
    DET = x * x - x * Mz + Cov_xy;
    Xcenter = (Mxz * (Myy - x) - Myz * Mxy) / DET / 2.;
    Ycenter = (Myz * (Mxx - x) - Mxz * Mxy) / DET / 2.;

    double cx    = Xcenter + meanX;
    double cy    = Ycenter + meanY;
    double rad   = sqrt(Xcenter * Xcenter + Ycenter * Ycenter + Mz);
    double sigma = sigmaCalc(points, Circle{cx, cy, rad, 0, 0, 0});
    
    return Circle{cx, cy, rad, sigma, 0, 0, iter};
}

// Levenberg-Marquardt geometrical circle fit
double fitLevenberg(std::vector<TVector3> points) {
    double size = points.size();
    double final_x, final_y, final_r;

    // Get initial prefitted values with the Taubin fit
    //auto initCirc = fitTaubin(points);
    //std::cout << "Taubin fit: " << initCirc.x() << " " << initCirc.y() << " " << initCirc.r() << std::endl;
    // Get iniitial prefitted values with the Hyper fit
    auto initCirc = fitHyper(points);
    std::cout << "Hyper fit: " << initCirc.x() << " " << initCirc.y() << " " << initCirc.r() << std::endl;

    // Levenberg-Marquardt geometrical circle fit 
    double factorUp = 10., factorDown = 0.04, ParLimit = 695, epsilon = 3.e-8;
    double lambda; // Control parameter for the Levenberg-Marquardt procedure (small positive number)
    double Mu, Mv, Muu, Mvv, Muv, Mr, UUl, VVl, Nl, F1, F2, F3, dX, dY, dR;
    double G11, G22, G33, G12, G13, G23, D1, D2, D3;
    double mnX, mnY, meanY, meanX;

    double sigma = sigmaCalc(points, Circle{initCirc.x(), initCirc.y(), initCirc.r(), 0, 0, 0});

    Circle newCirc = Circle{initCirc.x(), initCirc.y(), initCirc.r(), sigma, 0, 0};

    int code, iter, inner; // verbose, inner iterations
    int IterMAX = 99;
    lambda = 0.008;
    inner = code = 0;

    // Start loop
    for (iter = 0; iter < IterMAX; iter++) {
        Circle oldCirc = newCirc;
        if (++iter > IterMAX) {
            code = 1; 
            break;
            std::cout << "CODE: " << code << " maximum iterations have been reached." << std::endl;
        }
        // Compute moments
        Mu = Mv = Muu = Mvv = Muv = Mr = mnY, mnX, meanY, meanX = 0.;
        for (int i = 0; i < points.size(); i++) {
            double dx   = points[i].X() - oldCirc.x();
            double dy   = points[i].Y() - oldCirc.y();
            double ri   = std::sqrt(dx * dx + dy * dy);
            double u    = dx/ri;
            double v    = dy/ri;
            Mu  += u;
            Mv  += v;
            Muu += u*u;
            Mvv += v*v;
            Muv += u*v;
            Mr  += ri;
            mnX += points[i].Y();
            mnY += points[i].Y();
        }
        Mu   /= size;
        Mv   /= size;
        Muu  /= size;
        Mvv  /= size;
        Muv  /= size;
        Mr   /= size;
        meanX = mnX / size;
        meanY = mnY / size;

        // Compute matrices
        F1 = oldCirc.x() + oldCirc.r() * Mu - meanX;
        F2 = oldCirc.y() + oldCirc.r() * Mv - meanY;
        F3 = oldCirc.r() - Mr;

        // Update gradient
        double newGrad = std::sqrt(F1 * F1 + F2 * F2 + F3 * F3);

        // TRY AGAIN IS HERE
        try_again:

        // Cholesly decomposition
        UUl = Muu + lambda;
        VVl = Mvv + lambda;
        Nl = 1. + lambda;

        G11 = std::sqrt(UUl);
        G12 = Muv / G11;
        G13 = Mu / G11;
        G22 = std::sqrt(VVl - G12 * G12);
        G23 = (Mv - G12 * G13) / G22;
        G33 = std::sqrt(Nl - G13 * G13 - G23 * G23);
        
        D1 = F1 / G11;
        D2 = (F2 - G12 * D1) / G22;
        D3 = (F3 - G13 * D1 - G23 * D2) / G33;
        dR = D3 / G33;
        dY = (D2 - G23 * dR) / G22;
        dX = (D1 - G12*dY - G13 * dR) / G11;

        // Check procedure termination with tolerance
        if ((std::abs(dR) + std::abs(dX) + std::abs(dY)) / (1. + oldCirc.r()) < epsilon) {
            std::cout << "Tolerance has been reched." << std::endl;
            final_x = oldCirc.x();
            final_y = oldCirc.y();
            final_r = oldCirc.r();
            break;
        }
    
        // Update parameters otherwise
        // Check if center is out of the detector
        if ((std::abs(oldCirc.x()-dX) > ParLimit) || (std::abs(oldCirc.y()-dY) > ParLimit)) {
            code = 3; 
            final_x = oldCirc.x();
            final_y = oldCirc.y();
            final_r = oldCirc.r();
            std::cout << "CODE: " << code << " x = " << oldCirc.x()-dX << 
                " y = " << oldCirc.y()-dY << " " << oldCirc.x() << " " << 
                dX << " " << oldCirc.y() << " " << dY <<std::endl;
            break;
        } 
        // Check if radius is imposible
        if (oldCirc.r()-dR <= 0) {
            lambda *= factorUp;
            if (++inner > IterMAX) {
                code = 2; 
                final_x = oldCirc.x();
                final_y = oldCirc.y();
                final_r = oldCirc.r();
                std::cout << "CODE: " << code << " r = " << oldCirc.r()-dR << std::endl;
                break;}
            // otherwise try again with the updated lambda
            goto try_again;
        }
        // Calculate new sigma
        double newSigma = sigmaCalc(points, Circle{oldCirc.x()-dX, oldCirc.y()-dY, oldCirc.r()-dR, 0, 0, 0});
        // If there is improvement fact lambda up
        if (newSigma < oldCirc.sigma()) {lambda *= factorDown;}
        else {
            if (++inner > IterMAX) {
                code = 2; 
                final_x = oldCirc.x();
                final_y = oldCirc.y();
                final_r = oldCirc.r();
                std::cout << "CODE: " << code << " maximum lambda iterations" << std::endl;
                break;
            }
            lambda *= factorUp;
            // otherwise try again with the updated lambda
            goto try_again;
        }
        
        final_x = oldCirc.x();
        final_y = oldCirc.y();
        final_r = oldCirc.r();

        newCirc = Circle{oldCirc.x()-dX, oldCirc.y()-dY, oldCirc.r()-dR, newSigma, newGrad, iter, inner};
    }

    std::cout << "FINAL VALUES: " << final_x << " " << final_y << " " << final_r << std::endl;
    return final_r;
}


std::optional<Acts::BoundVector> estimateModified(
    const Acts::GeometryContext& gctx, std::vector<const SpacePoint*>::iterator spBegin,
    std::vector<const SpacePoint*>::iterator spEnd, const Acts::Surface& surface, const Acts::Vector3& bField,
    Acts::ActsScalar bFieldMin, const Acts::Logger& logger,
    Acts::ActsScalar mass) {

    // Check the number of provided space points
    size_t numSP = std::distance(spBegin, spEnd);
    if (numSP != 3) {
        std::cout << "There should be exactly three space points provided." << std::endl;
        return std::nullopt;
    }

    // Convert bField to Tesla
    Acts::ActsScalar bFieldInTesla = bField.norm() / Acts::UnitConstants::T;
    Acts::ActsScalar bFieldMinInTesla = bFieldMin / Acts::UnitConstants::T;
    // Check if magnetic field is too small
    if (bFieldInTesla < bFieldMinInTesla) {
        // @todo shall we use straight-line estimation and use default q/pt in such
        // case?
        std::cout << "The magnetic field at the bottom space point: B = "
                 << bFieldInTesla << " T is smaller than |B|_min = "
                 << bFieldMinInTesla << " T. Estimation is not performed." << std::endl;
        return std::nullopt;
    }

    // The global positions of the bottom, middle and space points
    std::array<Acts::Vector3, 3> spGlobalPositions = {Acts::Vector3::Zero(), Acts::Vector3::Zero(),
                                              Acts::Vector3::Zero()};
    // The first, second and third space point are assumed to be bottom, middle
    // and top space point, respectively
    Acts::Vector3 p;
    for (size_t isp = 0; isp < 3; ++isp) {
    std::vector<const SpacePoint*>::iterator it = std::next(spBegin, isp);
        if (*it == nullptr) {
            std::cout << "Empty space point found. This should not happen." << std::endl;
            return std::nullopt;
        }
        const auto& sp = *it;
        spGlobalPositions[isp] = Acts::Vector3(sp->x(), sp->y(), sp->z());
        if (isp == 0) {p = sp->p();}
    }

    // Define a new coordinate frame with its origin at the bottom space point, z
    // axis long the magnetic field direction and y axis perpendicular to vector
    // from the bottom to middle space point. Hence, the projection of the middle
    // space point on the tranverse plane will be located at the x axis of the new
    // frame.
    Acts::Vector3 relVec = spGlobalPositions[1] - spGlobalPositions[0];
    Acts::Vector3 newZAxis = bField.normalized();
    Acts::Vector3 newYAxis = newZAxis.cross(relVec).normalized();
    Acts::Vector3 newXAxis = newYAxis.cross(newZAxis);
    Acts::RotationMatrix3 rotation;
    rotation.col(0) = newXAxis;
    rotation.col(1) = newYAxis;
    rotation.col(2) = newZAxis;
    // The center of the new frame is at the bottom space point
    Acts::Translation3 trans(spGlobalPositions[0]);
    // The transform which constructs the new frame
    Acts::Transform3 transform(trans * rotation);

    // The coordinate of the middle and top space point in the new frame
    Acts::Vector3 local1 = transform.inverse() * spGlobalPositions[1];
    Acts::Vector3 local2 = transform.inverse() * spGlobalPositions[2];

    // Lambda to transform the coordinates to the (u, v) space
    auto uvTransform = [](const Acts::Vector3& local) -> Acts::Vector2 {
        Acts::Vector2 uv;
        Acts::ActsScalar denominator = local.x() * local.x() + local.y() * local.y();
        uv.x() = local.x() / denominator;
        uv.y() = local.y() / denominator;
        return uv;
    };
    // The uv1.y() should be zero
    Acts::Vector2 uv1 = uvTransform(local1);
    Acts::Vector2 uv2 = uvTransform(local2);

    // A,B are slope and intercept of the straight line in the u,v plane
    // connecting the three points
    Acts::ActsScalar A = (uv2.y() - uv1.y()) / (uv2.x() - uv1.x());
    Acts::ActsScalar B = uv2.y() - A * uv2.x();
    // Curvature (with a sign) estimate
    Acts::ActsScalar rho = -2.0 * B / std::hypot(1., A);
    // The projection of the top space point on the transverse plane of the new
    // frame
    Acts::ActsScalar rn = local2.x() * local2.x() + local2.y() * local2.y();
    // The (1/tanTheta) of momentum in the new frame,
    Acts::ActsScalar invTanTheta =
        local2.z() * std::sqrt(1. / rn) / (1. + rho * rho * rn);
    // The momentum direction in the new frame (the center of the circle has the
    // coordinate (-1.*A/(2*B), 1./(2*B)))
    Acts::Vector3 transDirection(1., A, std::hypot(1, A) * invTanTheta);
    double r = 0.5 * std::sqrt((A*A+1)/(B*B));
    std::cout << "Parameters: " << -1.*A/(2*B) << ", " << 1./(2*B) << ", " << r << std::endl;
    // Transform it back to the original frame
    Acts::Vector3 direction = rotation * transDirection.normalized();

    // ------------------ AMAIA --------------------------------------------------------------
    // Vector of circle parameters (x0, y0, radius) (in mm)
    std::vector<double> circle = {-1.*A/(2*B), 1./(2*B), 0.5 * std::sqrt((A*A+1)/(B*B))};
    std::vector<double> cerror = {5, 5, 20};

    const Int_t nParams = 3;
    const std::string parNames[nParams] = {"x0", "y0", "radius"};

    std::vector<TVector3> points;
    double x, y, z;

    ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

    // Set tolerance, etc.
    minimizer->SetMaxFunctionCalls(1000);
    minimizer->SetTolerance(0.01);
    minimizer->SetPrintLevel(0);

    // Create function wrapper for minimzer
    auto circleFitFunction = [&](const double *par) {
        // par[0] : X center of circle
        // par[1] : Y center of circle
        // par[2] : radius of circle

        double chi2   = 0; // Chi2 for minimization
        Int_t nGoodHit  = 0; // Number of good hits

        double x0  = par[0];
        double y0  = par[1];
        double r   = par[2];

        // Chi2 calculation
        for (auto it = points.begin(); it != points.end(); ++it) {
            double circFunc = std::sqrt(std::abs((x - x0) * (x - x0) + (y - y0) * (y - y0) - r * r)); // this value has to be 0
            chi2 += pow(circFunc, 2);
            //std::cout << nGoodHit << " " << chi2 << std::endl;
            nGoodHit++;
        }
        if (!nGoodHit) { return 1e10; }
        else { return chi2 / ((double)(nGoodHit)); }
    };

    ROOT::Math::Functor func(circleFitFunction, nParams);
    minimizer->SetFunction(func);

    for (size_t isp = 0; isp < 3; ++isp) {
        std::vector<const SpacePoint*>::iterator it = std::next(spBegin, isp);
        const auto& sp = *it;
        x = sp->x();
        y = sp->y();
        z = sp->z();
        points.push_back(TVector3{sp->x(), sp->y(), sp->z()});
    }

    double chi2 = 0;
    for (auto it = points.begin(); it != points.end(); ++it) {
            double circFunc = std::sqrt(std::abs((x - circle.at(0)) * (x - circle.at(0)) + 
            (y - circle.at(1)) * (y - circle.at(1)) - circle.at(2) * circle.at(2))); // this value has to be 0
            chi2 += pow(circFunc, 2);
    }

    double acts_minVal = chi2 / 3.0;

    // Get orthogonal line to the triplet
    TVector3 v1  = TVector3{points[2].X()-points[0].X(), points[2].Y()-points[0].Y(), 0};
    TVector3 v2  = v1.Orthogonal();
    double malda = v2.Y() / v2.X();
    double cnst  = points[1].Y() - malda * points[1].X(); 

    // Make quadratic function to find the intersection points between the circle and the
    // orthogonal line
    // TODO: Consider dynamic range depending on r
    TF1 *qF = new TF1("Quadratic", "(1+[0]*[0])*x*x+(2*[0]*([1]-[4])-2*[3])*x+([1]-[4])*([1]-[4])+[3]*[3]-[2]*[2]", -695, 695);
    qF->SetParameter(0, malda);
    qF->SetParameter(1, cnst);
    qF->SetParameter(2, circle.at(2));  // r (original from ACTS)
    qF->SetParameter(3, points[1].X()); // x
    qF->SetParameter(4, points[1].Y()); // y

    // Find first root of function
    // We need to find the two roots and select the one "inside" the angle of the points
    // in order to get the center of the circle
    // If there is no root, or the root is 2r or in the original SP, we already have the solution
    ROOT::Math::RootFinder finder;
    finder.Solve(*qF, -695, 695); // First root 
    double sol_x_1 = finder.Root();
    double sol_y_1 = malda * sol_x_1 + cnst;
    std::cout << sol_x_1 << " " << sol_y_1 << std::endl;
    finder.Solve(*qF, sol_x_1 + 0.1, 695); // Second root
    double sol_x_2 = finder.Root();
    double sol_y_2 = malda * sol_x_2 + cnst;
    std::cout << sol_x_2 << " " << sol_y_1 << std::endl;

    // Having the two potential points, find the one INSIDE the angle sustented by the points
    // Thus, the point with the smaller angle
    double angle_1  = TVector3{points[0].X() - sol_x_1, points[0].Y() - sol_y_1, 0}.Angle(
        TVector3{points[2].X() - sol_x_1, points[2].Y() - sol_y_1, 0});
    double angle_2  = TVector3{points[0].X() - sol_x_2, points[0].Y() - sol_y_2, 0}.Angle(
        TVector3{points[2].X() - sol_x_2, points[2].Y() - sol_y_2, 0});
    std::cout << angle_1 << " " << angle_2 << std::endl;

    double sol_x, sol_y;
    if (angle_1 < angle_2) {
        sol_x = sol_x_1;
        sol_y = sol_y_1;
    }
    else if (angle_1 > angle_2) {
        sol_x = sol_x_2;
        sol_y = sol_y_2;
    }

    // Second iteration with new circle fit center
    minimizer->SetLimitedVariable(0, parNames[0], sol_x, 0.02 * cerror.at(0),
        sol_x - 5, sol_x + 5);
    minimizer->SetLimitedVariable(1, parNames[1], sol_y, 0.02 * cerror.at(1),
        sol_y - 5, sol_y + 5);
    minimizer->SetLimitedVariable(2, parNames[2], circle.at(2), 0.02 * cerror.at(2),
        circle.at(2) -50, circle.at(2) +50);

    minimizer->Minimize();
    minimizer->PrintResults();
    const double *xs = minimizer->X();
    double new_minVal = minimizer->MinValue();

    std::cout << "FINAL RESULTS x = " << xs[0] << " y = " << xs[1] << " r = " << xs[2] << std::endl;

    double rLevenberg = fitLevenberg(points);
    // Gather all values (assume Bz = 1 T)
    double c = 2.99792458e8; // speed of light in vacuum
    double real_pT = std::sqrt(p.x() * p.x() + p.y() * p.y());
    double acts_pT = r * c / 1e9; // MeV/c
    double new_pT  = xs[2] * c / 1e9; // MeV/c
    double leven_pT = rLevenberg * c / 1e9;
    int sign = p.y() > 0 ? 1 : -1;
    double real_phi = sign * acos(p.x() / sqrt(p.x() * p.x() + p.y() * p.y()));
    double acts_phi = Acts::VectorHelpers::phi(direction);
    double real_theta = acos(p.z() / sqrt(p.x() * p.x() + p.y() * p.y() + p.z() * p.z()));
    double acts_theta = Acts::VectorHelpers::theta(direction);
    //if (new_minVal < acts_minVal) {
        estimatorFile << real_pT << " " << acts_pT << " " << new_pT << " " << leven_pT << " " <<
        real_phi << " " << acts_phi << " " << real_theta << " " << acts_theta << "\n"; 
    //}
    //else if (new_minVal > acts_minVal) {
       // estimatorFile << real_pT << " " << acts_pT << " " << acts_pT << "\n";
    //}
    // ---------------------------------------------------------------------------------------

    // Initialize the bound parameters vector
    Acts::BoundVector params = Acts::BoundVector::Zero();

    // The estimated phi and theta
    params[Acts::eBoundPhi] = Acts::VectorHelpers::phi(direction);
    params[Acts::eBoundTheta] = Acts::VectorHelpers::theta(direction);

    // Transform the bottom space point to local coordinates of the provided
    // surface
    auto lpResult = surface.globalToLocal(gctx, spGlobalPositions[0], direction);
    if (not lpResult.ok()) {
        std::cout << 
            "Global to local transformation did not succeed. Please make sure the "
            "bottom space point lies on the provided surface." << std::endl;
        return std::nullopt;
    }
    Acts::Vector2 bottomLocalPos = lpResult.value();
    // The estimated loc0 and loc1
    params[Acts::eBoundLoc0] = bottomLocalPos.x();
    params[Acts::eBoundLoc1] = bottomLocalPos.y();

    // The estimated q/pt in [GeV/c]^-1 (note that the pt is the projection of
    // momentum on the transverse plane of the new frame)
    Acts::ActsScalar qOverPt = rho * (Acts::UnitConstants::m) / (0.3 * bFieldInTesla);
    // The estimated q/p in [GeV/c]^-1
    params[Acts::eBoundQOverP] = qOverPt / std::hypot(1., invTanTheta);
    std::cout << rho << " " << 1/qOverPt << " " << params[Acts::eBoundQOverP] << std::endl;

    // The estimated momentum, and its projection along the magnetic field
    // diretion
    Acts::ActsScalar pInGeV = std::abs(1.0 / params[Acts::eBoundQOverP]);
    Acts::ActsScalar pzInGeV = 1.0 / std::abs(qOverPt) * invTanTheta;
    Acts::ActsScalar massInGeV = mass / Acts::UnitConstants::GeV;
    // The estimated velocity, and its projection along the magnetic field
    // diretion
    Acts::ActsScalar v = pInGeV / std::hypot(pInGeV, massInGeV);
    Acts::ActsScalar vz = pzInGeV / std::hypot(pInGeV, massInGeV);
    // The z coordinate of the bottom space point along the magnetic field
    // direction
    Acts::ActsScalar pathz = spGlobalPositions[0].dot(bField) / bField.norm();
    // The estimated time (use path length along magnetic field only if it's not
    // zero)
    if (pathz != 0) {
        params[Acts::eBoundTime] = pathz / vz;
    } else {
        params[Acts::eBoundTime] = spGlobalPositions[0].norm() / v;
    }

    if (params.hasNaN()) {
        std::cout << 
            "The NaN value exists at the estimated track parameters from seed with"
            << "\nbottom sp: " << spGlobalPositions[0] << "\nmiddle sp: "
            << spGlobalPositions[1] << "\ntop sp: " << spGlobalPositions[2] << std::endl;
        return std::nullopt;
    }
    return params;
}
#include <Eigen/Core>
#include <unsupported/Eigen/NonLinearOptimization>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>
#include <matplot/matplot.h>

#include "data.hpp"
#include "fitting.hpp"
#include "linalg.hpp"
#include "material.hpp"
#include "simulation.hpp"


void Fitting::genInPlaneWavevector() {
  // Cumulative sum of thicknesses
  matstack.z0.resize(matstack.numLayers - 1);
  matstack.z0(0) = 0.0;
  std::partial_sum(mThickness.begin(), mThickness.end(), std::next(matstack.z0.begin()), std::plus<double>());
  matstack.z0 -= (matstack.z0(mDipoleLayer - 1) + mDipolePosition);

  //getting sim data
//QUADRUPLE CHECK THIS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  matstack.u = Eigen::real(Eigen::sqrt(matstack.epsilon(matstack.numLayers - 1)/matstack.epsilon(mDipoleLayer)*(1- pow(Eigen::cos(mIntensityData.col(0)), 2))));
  matstack.x.resize(matstack.u.size());
  matstack.x = matstack.u.acos();
  matstack.numKVectors = matstack.u.size();
  
  // Differences
  matstack.dU = matstack.u.segment(1, matstack.u.size() - 1) - matstack.u.segment(0, matstack.u.size() - 1);
  matstack.dX = matstack.x.segment(1, matstack.x.size() - 1) - matstack.x.segment(0, matstack.x.size() - 1);
}

void Fitting::genOutofPlaneWavevector() {
  // Out of plane wavevector
  matstack.k = 2 * M_PI / mWvl / 1e-9 * matstack.epsilon.sqrt();
  matstack.h.resize(matstack.numLayers, matstack.u.size());
  matstack.h = matstack.k(mDipoleLayer) *
               (((matstack.epsilon.replicate(1, matstack.x.size())) / matstack.epsilon(mDipoleLayer)).rowwise() -
                 matstack.u.pow(2).transpose()).sqrt();
}

void Fitting::discretize() {
  loadMaterialData();
  genInPlaneWavevector();
  genOutofPlaneWavevector();
}

Fitting::Fitting(const std::vector<Material>& materials,
  const std::vector<double>& thickness,
  const size_t dipoleLayer,
  const double dipolePosition,
  const double wavelength,
  const std::string& fittingFilePath) :
  BaseSolver(materials,
    thickness,
    dipoleLayer,
    dipolePosition,
    wavelength) // initialization must be performed this way due to const members 
  {
  if (materials.size() != thickness.size() + 2) {
    throw std::runtime_error("Invalid Input! Number of materials different than number of layers.");
  }
  if (dipoleLayer >= materials.size()) { throw std::runtime_error("Invalid Input! Dipole position is out of bounds."); }
  // Log initialization of Simulation
  std::cout << "\n\n\n"
            << "-----------------------------------------------------------------\n";
  std::cout << "              Initializing Fitting             \n";
  std::cout << "-----------------------------------------------------------------\n"
            << "\n\n";
  mIntensityData = Data::loadFromFile(fittingFilePath, 2);
  discretize();
  calculate();
  //setting up functor for fitting
  mResidual.intensities = mIntensityData.col(1);
  mResidual.powerGlass = calculateEmissionSubstrate();
  }

Matrix Fitting::calculateEmissionSubstrate() {
  Vector powerPerppPolGlass;
  Vector powerParapPolGlass;

  powerPerppPolGlass = ((Eigen::real(mPowerPerpUpPol(matstack.numLayers - 2, Eigen::all))) *
                    std::sqrt(std::real(matstack.epsilon(matstack.numLayers - 1) / matstack.epsilon(mDipoleLayer))));
  powerPerppPolGlass /= Eigen::sin(mIntensityData.col(0));

  powerParapPolGlass = ((Eigen::real(mPowerParaUpPol(matstack.numLayers - 2, Eigen::all))) *
                    std::sqrt(std::real(matstack.epsilon(matstack.numLayers - 1) / matstack.epsilon(mDipoleLayer))));
  powerParapPolGlass /= Eigen::sin(mIntensityData.col(0));


  Matrix powerGlass(2, powerPerppPolGlass.size());
  powerGlass.row(0) = powerPerppPolGlass;
  powerGlass.row(1) = powerParapPolGlass;
  return powerGlass;
}

int ResFunctor::operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const {
  // x here is vector of fitting params
  for (size_t i = 0; i < intensities.size(); ++i) {
    fvec(i) = intensities(i) - x(0) * (x(1)*powerGlass(0, i) + (1 - x(1))*powerGlass(1, i)); // residual of each sample
  }
  return 0;
} 

int ResFunctor::inputs() const {return 2;}

int ResFunctor::values() const {return intensities.size();}

std::pair<Eigen::VectorXd, Eigen::ArrayXd> Fitting::fitEmissionSubstrate() {
  //returns the vector of parameters and the fitted intensities as a std::pair

  std::vector<double> theta(matstack.x.rows()), yFit(matstack.x.rows()), yExp(mResidual.intensities.rows());

  Eigen::ArrayXd::Map(&theta[0], mIntensityData.rows()) = mIntensityData.col(0);
  Eigen::ArrayXd::Map(&yExp[0], mIntensityData.rows()) = mResidual.intensities;

  // Setup
  Eigen::VectorXd x(2);
  // Initial guess
  x(0) = 1.0;
  x(1) = 0.0;

  Eigen::LevenbergMarquardt<ResFunctorNumericalDiff> lm(mResidual);
  lm.parameters.maxfev = 2000;
  lm.parameters.xtol = 1e-10;
  lm.parameters.ftol = 1e-10;

  int status = lm.minimize(x);
  std::cout << "Number of iterations: " << lm.iter << '\n';
  std::cout << "Status: " << status << '\n';
  std::cout << "Fitting result: " << x << '\n' << '\n';

  // simulation results
  Eigen::ArrayXd optIntensities(matstack.x.rows());
  optIntensities = x(0) * (x(1)*mResidual.powerGlass.row(0) + (1 - x(1))*mResidual.powerGlass.row(1));
  Eigen::ArrayXd::Map(&yFit[0], matstack.x.rows()) = optIntensities;

  //plotting the results
  matplot::figure();
  
  matplot::scatter(theta, yExp);
  matplot::hold(matplot::on);
  matplot::plot(theta, yFit)->line_width(2).color("red");
  matplot::show();

  return std::pair<Eigen::VectorXd, Eigen::ArrayXd>(x, optIntensities);
};
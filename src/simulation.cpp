#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

#include "linalg.hpp"
#include "data.hpp"
#include "material.hpp"
#include "simulation.hpp"

void Simulation::genInPlaneWavevector()
{
  // Cumulative sum of thicknesses
  matstack.z0.resize(matstack.numLayers - 1);
  matstack.z0(0) = 0.0;
  std::vector<double> thicknesses;
  for (size_t i=1; i < mLayers.size()-1; ++i) {
    thicknesses.push_back(mLayers[i].getThickness());
  }
  std::partial_sum(thicknesses.begin(), thicknesses.end(), std::next(matstack.z0.begin()), std::plus<double>());
  matstack.z0 -= (matstack.z0(mDipoleLayer - 1) + mDipolePosition);

  // Discretization of in-plane wavevector
  CMPLX I(0.0, 1.0);
  double x_res = 5e-4;
  if (_mode == SimulationMode::AngleSweep) {
    matstack.x = arange<Vector>(_sweepStart * M_PI / 180 + x_res, _sweepStop * M_PI / 180 - x_res, x_res);
    matstack.u = Eigen::real(Eigen::sqrt(matstack.epsilon(matstack.numLayers - 1)/matstack.epsilon(mDipoleLayer)*(1- pow(Eigen::cos(matstack.x), 2))));
  }
  else if (_mode == SimulationMode::ModeDissipation) {
    Vector x_real = arange<Vector>(-std::acos(_sweepStart), -x_res, x_res);
    CVector x_imag = I * arange<Vector>(x_res, -std::acos(std::complex<double>(_sweepStop, 0.0)).imag(), x_res);
    matstack.x.resize(x_real.rows() + x_imag.rows());
    matstack.x.head(x_real.size()) = x_real.cast<CMPLX>();
    matstack.x.segment(x_real.size(), x_imag.size()) = x_imag;
    matstack.u = matstack.x.cos().real();
  }
  else {throw std::runtime_error("Invalid mode!");}

  matstack.numKVectors = matstack.u.size();

  // Differences
  matstack.dU = matstack.u.segment(1, matstack.u.size() - 1) - matstack.u.segment(0, matstack.u.size() - 1);
  matstack.dX = matstack.x.segment(1, matstack.x.size() - 1) - matstack.x.segment(0, matstack.x.size() - 1);
}

void Simulation::genOutofPlaneWavevector()
{
  // Out of plane wavevector
  matstack.k = 2 * M_PI / mWvl / 1e-9 * matstack.epsilon.sqrt();
  matstack.h.resize(matstack.numLayers, matstack.u.size());
  matstack.h = matstack.k(mDipoleLayer) *
               (((matstack.epsilon.replicate(1, matstack.u.size())) / matstack.epsilon(mDipoleLayer)).rowwise() -
                 matstack.u.pow(2).transpose())
                 .sqrt();
}

void Simulation::discretize()
{
  loadMaterialData();
  genInPlaneWavevector();
  genOutofPlaneWavevector();
}

void Simulation::init() {
  // Log initialization of Simulation
  std::cout << "\n\n\n"
            << "-----------------------------------------------------------------\n";
  std::cout << "              Initializing Simulation             \n";
  std::cout << "-----------------------------------------------------------------\n"
            << "\n\n";
  this->discretize();
}

Simulation::Simulation(SimulationMode mode,
  const std::vector<Layer>& layers,
  const double dipolePosition,
  const double wavelength,
  const double sweepStart,
  const double sweepStop) :
  BaseSolver(mode,
    layers,
    dipolePosition,
    wavelength,
    sweepStart,
    sweepStop)
{
  init();
}

Simulation::Simulation(SimulationMode mode,
      const std::vector<Layer>& layers,
      const double dipolePosition,
      const std::string& spectrumFile,
      const double sweepStart,
      const double sweepStop) :
      BaseSolver(mode,
                 layers,
                 dipolePosition,
                 spectrumFile,
                 sweepStart,
                 sweepStop)
{
  init();
}

Simulation::Simulation(SimulationMode mode,
      const std::vector<Layer>& layers,
      const double dipolePosition,
      const GaussianSpectrum& spectrum,
      const double sweepStart,
      const double sweepStop) :
      BaseSolver(mode,
                 layers,
                 dipolePosition,
                 spectrum,
                 sweepStart,
                 sweepStop)
{
  init();
}

Simulation::Simulation(SimulationMode mode,
      const std::vector<Layer>& layers,
      const DipoleDistribution& dipoleDist,
      const GaussianSpectrum& spectrum,
      const double sweepStart,
      const double sweepStop) :
      BaseSolver(mode,
                 layers,
                 dipoleDist,
                 spectrum,
                 sweepStart,
                 sweepStop)
{
  init();
}
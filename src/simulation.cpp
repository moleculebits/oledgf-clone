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
  std::partial_sum(mThickness.begin(), mThickness.end(), std::next(matstack.z0.begin()), std::plus<double>());
  matstack.z0 -= (matstack.z0(mDipoleLayer - 1) + mDipolePosition);

  // Discretization of in-plane wavevector
  CMPLX I(0.0, 1.0);
  double x_res = 5e-4;
  Vector x_real = arange<Vector>(-M_PI_2, -x_res, x_res);
  CVector x_imag = I * arange<Vector>(x_res, 1.8, x_res);
  matstack.x.resize(x_real.rows() + x_imag.rows());
  matstack.x.head(x_real.size()) = x_real.cast<CMPLX>();
  matstack.x.segment(x_real.size(), x_imag.size()) = x_imag;
  matstack.u = matstack.x.cos().real();
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
               (((matstack.epsilon.replicate(1, matstack.x.size())) / matstack.epsilon(mDipoleLayer)).rowwise() -
                 matstack.x.cos().pow(2).transpose())
                 .sqrt();
}

void Simulation::discretize()
{
  loadMaterialData();
  genInPlaneWavevector();
  genOutofPlaneWavevector();
}

Simulation::Simulation(const std::vector<Material>& materials,
  const std::vector<double>& thickness,
  const size_t dipoleLayer,
  const double dipolePosition,
  const double wavelength) :
  BaseSolver(materials,
    thickness,
    dipoleLayer,
    dipolePosition,
    wavelength) // initialization must be performed this way due to const members
{
  // Log initialization of Simulation
  std::cout << "\n\n\n"
            << "-----------------------------------------------------------------\n";
  std::cout << "              Initializing Simulation             \n";
  std::cout << "-----------------------------------------------------------------\n"
            << "\n\n";
  discretize();
}

GaussianSpectrum::GaussianSpectrum(double xmin,
                                   double xmax,
                                   double x0,
                                   double sigma)
{
  spectrum.resize(50 , 2);
  Eigen::ArrayXd x = Eigen::ArrayXd::LinSpaced(50, xmin, xmax);
  spectrum.col(0) = x;
  spectrum.col(1) = (1.0/sqrt(2 * M_PI * pow(sigma, 2))) * (-0.5 * ((x - x0) / sigma).pow(2)).exp();
}

SimulationSweep::SimulationSweep(const std::vector<Material>& materials,
      const std::vector<double>& thickness,
      const size_t dipoleLayer,
      const double dipolePosition,
      const std::string& spectrumFile) :
      Simulation(materials,
                 thickness,
                 dipoleLayer,
                 dipolePosition,
                 -1.0)
{
  _spectrum = Data::loadFromFile(spectrumFile, 2);
}

SimulationSweep::SimulationSweep(const std::vector<Material>& materials,
      const std::vector<double>& thickness,
      const size_t dipoleLayer,
      const double dipolePosition,
      const GaussianSpectrum& spectrum) :
      Simulation(materials,
                 thickness,
                 dipoleLayer,
                 dipolePosition,
                 -1.0)
{
  _spectrum = std::move(spectrum.spectrum);
}

void SimulationSweep::calculate() {
  CMatrix pPerpUpPol = CMatrix::Zero(matstack.numLayers - 1, matstack.u.size() - 1);
  CMatrix pParaUpPol = CMatrix::Zero(matstack.numLayers - 1, matstack.u.size() - 1);
  CMatrix pParaUsPol = CMatrix::Zero(matstack.numLayers - 1, matstack.u.size() - 1);
  double dX = _spectrum(1, 0) - _spectrum(0, 0);
  for (Eigen::Index i=0; i < _spectrum.rows(); ++i) {
    mWvl = _spectrum(i, 0);
    discretize();
    BaseSolver::calculate();
    // Integration
    if (i == 0 || i == _spectrum.rows() - 1) {
      mPowerPerpUpPol *= 0.5;
      mPowerParaUpPol *= 0.5;
      mPowerParaUsPol *= 0.5;
    }
    pPerpUpPol += (mPowerPerpUpPol * _spectrum(i, 1)) ;
    pParaUpPol += (mPowerParaUpPol * _spectrum(i, 1));
    pParaUsPol += (mPowerParaUsPol * _spectrum(i, 1));
  }
  mPowerPerpUpPol = pPerpUpPol * dX;
  mPowerParaUpPol = pParaUpPol * dX;
  mPowerParaUsPol = pParaUsPol * dX;
}
    
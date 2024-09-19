#pragma once

#include <vector>

#include <simulation.hpp>
#include <jsonsimplecpp/parser.hpp>

enum struct SolverMode {Simulation, Fitting};

struct Input {
    SolverMode mode;

    std::vector<Material> materials;

    std::vector<double> thicknesses;
    
    Eigen::Index emitterIndex;

    double emitterPosition;
    DipoleDistribution dipoleDist;
    bool hasDipoleDistribution;


    double wavelength;
    GaussianSpectrum spectrum;
    std::string spectrumFile;
    bool hasSpectrum;
};

class ConfigurationManager {
   Json::JsonParser jsonParser; 
   Input input;

   public:
    ConfigurationManager(const std::string& filename);

    void configure();
    const Input& getInput() const;

};

class SimulationManager {
    ConfigurationManager _configurator;

    public:
        SimulationManager(const std::string& filename);

        std::unique_ptr<BaseSolver> create();

};
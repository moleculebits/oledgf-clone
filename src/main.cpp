#include <filesystem>
#include <iostream>

#include <fmt/core.h>
#include <fmt/ranges.h>
#include <linalg.hpp>
#include <material.hpp>
#include <matplotlib.hpp>
#include <matrix.hpp>

int main()
{
  // Simple Test fmtlib
  fmt::print("Hello World!\n");
  // Test input data
  std::filesystem::path dataFile("../mat/alq3_literature.dat");

  Material alq3 = Material(dataFile, ',');
  std::vector<double> alq3Data = alq3.getRefIndex();

  // Test slicing. alq3Data contains (wvl, n, k) as | wvl | n | k | wvl | n | k | ...
  // Here we print the first 10 values of wavelength and n
  std::vector<double> wvl = slice(alq3Data, 0, 30, 3);
  std::vector<double> alq3N = slice(alq3Data, 1, 30, 3);
  // Test fmtlib to print containers
  fmt::print("{}\n", wvl);
  fmt::print("{}\n", alq3N);
  // Test matplotlib embedding
  // plot(wvl, alq3N);
  // save("./test.png");
  // show();

  // Test Matrix class
  Matrix<double> m(3, 3);
  for (size_t i = 0; i < m.size().first; ++i) {
    for (size_t j = 0; j < m.size().second; ++j) {
      if (i == j) m(i, j) = 1.0;
    }
  }

  for (size_t i = 0; i < m.size().first; ++i) {
    for (size_t j = 0; j < m.size().second; ++j) { std::cout << m(i, j) << ", "; }
    std::cout << '\n';
  }
}
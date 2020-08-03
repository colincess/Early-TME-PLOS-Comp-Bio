#include <array>
#include <string>
#include <vector>
#include "CellGrids.h"
#include "diffusibles.h"

#ifndef MTC_MACROPHAGE_H
#define MTC_MACROPHAGE_H

class Macrophage{
public:
    Macrophage(std::array<int, 2> loc, double ifngBase, double il4Base,
               std::vector<std::vector<std::vector<double>>> nnWeights,
               std::vector<std::vector<std::vector<double>>> nnBiases,
               int initial, int identity);
    void migrate(CellGrids &cg, Diffusibles &diff);
    void simulate(double tstep, CellGrids &cg, Diffusibles &diff, double depletion);
    int activate(std::vector<std::vector<double>> outside);
    int neuralNetwork(std::vector<std::vector<double>> outside);
    std::vector<std::vector<double>> sigmoid(std::vector<std::vector<double>> x);
    std::vector<std::vector<double>> matmul(std::vector<std::vector<double>> a, std::vector<std::vector<double>> b);
    std::string state;
    std::array<int, 2> location;
    double IFNG0;
    double IFNGR0;
    double IL40;
    double IL4R0;
    double SOCS10;
    double SOCS30;
    double k15;
    double arg1Prod;

    double engagement;
    double idx;

private:
    double life_span;  // hr, max age
    double age;           // hr, current age

    double dying;

    double activationThreshold;

    std::vector<std::vector<double>> w1;
    std::vector<std::vector<double>> wo;
    std::vector<std::vector<double>> b1;
    std::vector<std::vector<double>> bo;

    std::vector<std::vector<std::vector<double>>> weights;
    std::vector<std::vector<std::vector<double>>> biases;

    double reDiff;

};

#endif //MTC_MACROPHAGE_H

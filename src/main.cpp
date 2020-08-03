#include "environment.h"
#include <iostream>
#include <omp.h>
#include <fstream>
#include <random>
#include <thread>
#include <chrono>
#include "macrophage.h"
#include <sstream>

/*
 * WHAT TO CLEAN UP BEFORE FINAL RUNS
 * ----------------------------------
 * - remove cancer stem cells
 * - I'm pretty sure max cell divisions for T cells doesn't matter
 * - choose macrophage lifespan, apparently it's "months to years"
 */

int main(int argc, char **argv){
    std::string folder = argv[1];
    std::string set = argv[2];
    int N = std::stoi(argv[3]);

    double ifng = std::stod(argv[4]);
    double il4 = std::stod(argv[5]);

    double simTime = 200;

    int perturb = std::stoi(argv[6]);
    double perturbTime = std::stod(argv[7]);
    double perturbLvl = std::stod(argv[8]);

    double treatmentModulation = std::stod(argv[9]);
    double timeOn = std::stod(argv[10]);

    std::string str = "mkdir -p ./"+folder+"/set_" + set;
    const char *command = str.c_str();
    std::system(command);

    double start = omp_get_wtime();
    //std::cout << "Running...\n";
    for (int i = 0; i < N; i++) {
        std::cout << i << std::endl;
        Environment model(0.5, folder, std::stoi(set), ifng, il4, perturb, perturbTime, perturbLvl);
        model.simulate(simTime, treatmentModulation, timeOn);
    }
    double stop = omp_get_wtime();
    std::cout << "Duration: " << (stop-start)/(60) << std::endl;

    return 0;
}

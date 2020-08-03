#include "environment.h"
#include <iostream>
#include <random>
#include <fstream>
#include <sstream>
#include <algorithm>

/*
 * See if M1 should recruit more T cells
 *
 * Work On
 * -------
 * - clean up unused variables and shit
 *
 * Perturbations
 * -------------
 * 0 - Recruitment inhibition
 * 1 - Macrophage depletion
 * 2 - PI3K
 * 3 - CAR T Cells
 *
 * LHS Params
 * ----------
 * 0 - macProb: environment.cpp:92
 * 1 - tumor IL-4 secretion: environment.cpp:54
 * 2 - t cell IFNG secretion: environment.cpp:79
 * 3 - cancer proliferation rate: cancer.cpp:19
 * 4 - M2 IL-4 secretion: environment.cpp:1123
 * 5 - macrophage lifespan (days): macrophage.cpp:19
 */

std::random_device rd;
std::vector<double> paramsLHS;

int totalCancerCells;

Environment::Environment(double stepSize, std::string folder, int set, double ifng, double il4, int perturb, double perturbTime, double perturbLvl){
    // infg is the secretion rate for T cells

    // directory to save time courses to
    saveDir = "./"+folder+"/set_"+std::to_string(set);

    // load neural network weights
    loadNNParameters();
    std::cout << "Loaded neural network parameters...\n";

    // load other parameters
    loadParameters(perturb,perturbTime,perturbLvl);
    std::cout << "Loaded LHS parameters\n";
    std::cout << "Num LHS params: " << paramsLHS.size() << std::endl;

    treatmentOn = 0;
    treatmentOff = 0;

    // base cytokine production rates of tumor cells
    baseIFNG = 0;
    baseIL4 = il4; //paramsLHS[1];//il4;

    // base macrophage receptor expression
    // not changing this at the moment
    // because there isn't much of an affect by doing so
    baseIFNR = 106;
    baseIL4R = 16;

    // other parameters
    tstep = stepSize; // hours

    // diffusion parameters
    D = 3e-7; // Wells 2015
    dx = 0.0015; // from Hao 2018 it seems that the avg cancer cell diameter is around 15 um
    // The following secretions rates are in pg/(cell*sec)
    // taken from Wells 2015
    k_M1f = 1e-7;//6e-8;

    // based on Han 2010
    // sorta made up for now
    // molecules/second
    k_IFNG_T = ifng;//paramsLHS[2];//ifng;

    // T Cell recruitment parameters
    // Gong 2017
    tdelay = 5; // days
    twindow = 1; // days
    ka = 15; // dimensionless
    ki = 0.01; // dimensionless
    r1 = 6.0; // cells/hr. Originally, theirs was 1/time_step, with time_step = 10 min

    // i'm using macProb as a probability for now
    macProb = 2*1e-8;//paramsLHS[0];//1e-8; // cells/(lattice site * sec * vasculature) Wells 2015

    endTime = 0;
}

void Environment::loadNNParameters() {
    // loads neural network weights and biases
    // for the macrophage model

    int n_layers = 2;

    for(int i=1; i<n_layers+1; ++i){
        std::ifstream dataW("./nnParams/w"+std::to_string(i)+".csv");
        std::string line;
        std::vector<std::vector<double>> matrixW;
        while(std::getline(dataW, line)){
            std::stringstream lineStream(line);
            std::string cell;
            std::vector<double> parsedRow;
            while(std::getline(lineStream, cell, ',')){
                parsedRow.push_back(std::stod(cell));
            }
            matrixW.push_back(parsedRow);
        }

        macrophageWeights.push_back(matrixW);

        std::ifstream dataB("./nnParams/b"+std::to_string(i)+".csv");
        std::vector<std::vector<double>> matrixB;
        while(std::getline(dataB, line)){
            std::stringstream lineStream(line);
            std::string cell;
            std::vector<double> parsedRow;
            while(std::getline(lineStream, cell, ',')){
                parsedRow.push_back(std::stod(cell));
            }
            matrixB.push_back(parsedRow);
        }

        macrophageBiases.push_back(matrixB);
    }
}

void Environment::loadParameters(int perturb, double time, double level) {
    // loads parameters that are allowed to change from
    // simulation to simulation

    // load treatment parameters
    std::vector<double> params = {0,0,0,0,0,0,0,0};

    params[perturb*2] = time;
    params[perturb*2 + 1] = level;

    recTime = params[0]; // time (days) when treatment is introduced
    reclvl = params[1]; // percent to reduce k15
    depTime = params[2]; // time (days) when treatment is introduced
    deplvl = params[3]; // percent to reduce arg1Prod
    PI3KTime = params[4]; // time (days) when treatment is introduced
    PI3Klvl = params[5]; // percent to reduce macprob in macrophage recruitment

    // load LHS parameters
    /*std::ifstream dataLHS(saveDir+"/paramsLHS.csv");
    std::string line;
    while (std::getline(dataLHS, line)) {
        std::stringstream lineStream(line);
        std::string cell;
        std::vector<double> parsedRow;
        while (std::getline(lineStream, cell, ',')) {
            parsedRow.push_back(std::stod(cell));
        }
        paramsLHS.push_back(parsedRow[0]);
    }*/

}

void Environment::plot(CellGrids &cg, Diffusibles &diff, int s){
    // outputs final spatial state of the simulation

    std::string str = "mkdir -p "+saveDir+"/spatial/"+std::to_string(s);
    std::system(str.c_str());

    std::ofstream myfile;
    myfile.open(saveDir+"/spatial/"+std::to_string(s)+"/ccg.csv");
    for(int i=0; i<100; ++i){
        myfile << cg.ccg[i][0] + 2*cg.m0g[i][0] + 3*cg.m1g[i][0] + 4*cg.m2g[i][0] + 5*cg.c8g[i][0] + 6*cg.actT[i][0];
        for(int j=1; j<100; ++j){
            myfile << "," << cg.ccg[i][j] + 2*cg.m0g[i][j] + 3*cg.m1g[i][j] + 4*cg.m2g[i][j] + 5*cg.c8g[i][j] + 6*cg.actT[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();

    /*std::ofstream myfile;
    myfile.open(saveDir+"/spatial/"+std::to_string(s)+"/ccg.csv");
    for(int i=0; i<100; ++i){
        myfile << cg.ccg[i][0];
        for(int j=1; j<100; ++j){
            myfile << "," << cg.ccg[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();

    myfile.open(saveDir+"/spatial/"+std::to_string(s)+"/m0g.csv");
    for(int i=0; i<100; ++i){
        myfile << cg.m0g[i][0];
        for(int j=1; j<100; ++j){
            myfile << "," << cg.m0g[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();

    myfile.open(saveDir+"/spatial/"+std::to_string(s)+"/m1g.csv");
    for(int i=0; i<100; ++i){
        myfile << cg.m1g[i][0];
        for(int j=1; j<100; ++j){
            myfile << "," << cg.m1g[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();

    myfile.open(saveDir+"/spatial/"+std::to_string(s)+"/m2g.csv");
    for(int i=0; i<100; ++i){
        myfile << cg.m2g[i][0];
        for(int j=1; j<100; ++j){
            myfile << "," << cg.m2g[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();

    myfile.open(saveDir+"/spatial/"+std::to_string(s)+"/c8g.csv");
    for(int i=0; i<100; ++i){
        myfile << cg.c8g[i][0];
        for(int j=1; j<100; ++j){
            myfile << "," << cg.c8g[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();

    myfile.open(saveDir+"/spatial/"+std::to_string(s)+"/actT.csv");
    for(int i=0; i<100; ++i){
        myfile << cg.actT[i][0];
        for(int j=1; j<100; ++j){
            myfile << "," << cg.actT[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();*/

    myfile.open(saveDir+"/spatial/"+std::to_string(s)+"/IFNG.csv");
    for(int i=0; i<100; ++i){
        myfile << diff.IFNG[i][0];
        for(int j=1; j<100; ++j){
            myfile << "," << diff.IFNG[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();

    myfile.open(saveDir+"/spatial/"+std::to_string(s)+"/IL4.csv");
    for(int i=0; i<100; ++i){
        myfile << diff.IL4[i][0];
        for(int j=1; j<100; ++j){
            myfile << "," << diff.IL4[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();
}

void Environment::save(){
    // outputs time courses of the simulation

    std::fstream myfile;

    /*myfile.open(saveDir+"/treatmentOnOff.csv", std::ios_base::app);
    myfile << treatmentOnOff[0];
    for(int i=1; i<treatmentOnOff.size(); ++i){
        myfile << "," << treatmentOnOff[i];
        //std::cout << i;
    }
    myfile << std::endl;
    myfile.close();*/

    /*myfile.open(saveDir+"/times.csv", std::ios_base::app);
    myfile << times[0];
    for(int i=1; i<times.size(); i++){
        myfile << "," << times[i];
    }
    myfile << std::endl;
    myfile.close();*/

    myfile.open(saveDir+"/cc_num.csv", std::ios_base::app);
    myfile << cc_num[0];
    for(int i=1; i<cc_num.size(); i++){
        myfile << "," << cc_num[i];
    }
    myfile << std::endl;
    myfile.close();

    int maxCancer = 0;
    int maxActT = 0;
    int maxC8 = 0;
    int maxM1 = 0;
    int maxM2 = 0;
    for(int i=0; i<cc_num.size(); ++i){
        if(cc_num[i] > maxCancer){maxCancer = cc_num[i];}
        if(actT_num[i] > maxActT){maxActT = actT_num[i];}
        if(c8_num[i] > maxC8){maxC8 = c8_num[i];}
        if(m1_num[i] > maxM1){maxM1 = m1_num[i];}
        if(m2_num[i] > maxM2){maxM2 = m2_num[i];}
    }

    myfile.open(saveDir+"/assortedInfo.csv", std::ios_base::app);
    myfile << times[times.size()-1] << ",";
    myfile << cc_num[cc_num.size()-1] << ",";
    myfile << maxCancer << ",";
    myfile << maxM1 << ",";
    myfile << maxM2 << ",";
    myfile << maxC8 << ",";
    myfile << maxActT;
    myfile << std::endl;
    myfile.close();

    myfile.open(saveDir+"/m0_num.csv", std::ios_base::app);
    myfile << m0_num[0];
    for(int i=1; i<m0_num.size(); i++){
        myfile << "," << m0_num[i];
    }
    myfile << std::endl;
    myfile.close();

    myfile.open(saveDir+"/m1_num.csv", std::ios_base::app);
    myfile << m1_num[0];
    for(int i=1; i<m1_num.size(); i++){
        myfile << "," << m1_num[i];
    }
    myfile << std::endl;
    myfile.close();

    myfile.open(saveDir+"/m2_num.csv", std::ios_base::app);
    myfile << m2_num[0];
    for(int i=1; i<m2_num.size(); i++){
        myfile << "," << m2_num[i];
    }
    myfile << std::endl;
    myfile.close();

    myfile.open(saveDir+"/c8_num.csv", std::ios_base::app);
    myfile << c8_num[0];
    for(int i=1; i<c8_num.size(); i++){
        myfile << "," << c8_num[i];
    }
    myfile << std::endl;
    myfile.close();

    myfile.open(saveDir+"/actT_num.csv", std::ios_base::app);
    myfile << actT_num[0];
    for(int i=1; i<actT_num.size(); i++){
        myfile << "," << actT_num[i];
    }
    myfile << std::endl;
    myfile.close();

    myfile.open(saveDir+"/avgIL4.csv", std::ios_base::app);
    myfile << IL4_avg[0];
    for(int i=1; i<IL4_avg.size(); ++i){
	myfile << "," << IL4_avg[i];
    }
    myfile << std::endl;
    myfile.close();

    myfile.open(saveDir+"/maxIL4.csv", std::ios_base::app);
    myfile << IL4_max[0];
    for(int i=1; i<IL4_max.size(); ++i){
	myfile << "," << IL4_max[i];
    }
    myfile << std::endl;
    myfile.close();

    myfile.open(saveDir+"/avgIFNG.csv", std::ios_base::app);
    myfile << IFNG_avg[0];
    for(int i=1; i<IFNG_avg.size(); ++i){
	myfile <<  "," << IFNG_avg[i];
    }
    myfile << std::endl;
    myfile.close();

    myfile.open(saveDir+"/maxIFNG.csv", std::ios_base::app);
    myfile << IFNG_max[0];
    for(int i=1; i<IFNG_max.size(); ++i){
	myfile << "," << IFNG_max[i];
    }
    myfile << std::endl;
    myfile.close();

    /*myfile.open(saveDir+"/mIFNGR0.csv", std::ios_base::app);
    myfile << IFNGRExpression[0];
    for(int i=1; i<IFNGRExpression.size(); i++){
        myfile << "," << IFNGRExpression[i];
    }
    myfile << std::endl;
    myfile.close();

    myfile.open(saveDir+"/mIL4R0.csv", std::ios_base::app);
    myfile << IL4RExpression[0];
    for(int i=1; i<IL4RExpression.size(); i++){
        myfile << "," << IL4RExpression[i];
    }
    myfile << std::endl;
    myfile.close();

    myfile.open(saveDir+"/mSOCS10.csv", std::ios_base::app);
    myfile << SOCS1Expression[0];
    for(int i=1; i<SOCS1Expression.size(); i++){
        myfile << "," << SOCS1Expression[i];
    }
    myfile << std::endl;
    myfile.close();

    myfile.open(saveDir+"/mSOCS30.csv", std::ios_base::app);
    myfile << SOCS3Expression[0];
    for(int i=1; i<SOCS3Expression.size(); i++){
        myfile << "," << SOCS3Expression[i];
    }
    myfile << std::endl;
    myfile.close();*/
}

void Environment::clean(CellGrids &cg){
    // delete dead cells from lists

    for(int i=1; i<99; i++){
        for(int j=1; j<99; j++){
            for(int k=1; k<99; k++){
                cg.ccg[i][j] = 0;
                cg.ccid[i][j] = -1;
                cg.mid[i][j] = -1;
                cg.allCells[i][j] = 0;
                cg.mpg[i][j] = 0;
                cg.m0g[i][j] = 0;
                cg.m1g[i][j] = 0;
                cg.m2g[i][j] = 0;
                cg.c8g[i][j] = 0;
            }
        }
    }

    std::vector<Cancer> new_cancer;
    int z = 0;
    for(auto & cell : cc_list){
        if(cell.state!="dead" && cell.state!="newly_dead"){
            int i = cell.location[0];
            int j = cell.location[1];

            cell.idx = z;

            new_cancer.push_back(cell);

            cg.allCells[i][j] = 1;
            cg.ccg[i][j] = 1;
            cg.ccid[i][j] = z;
            z++;
        }
    }

    std::vector<CD8> new_c8;
    for(auto & cell : c8_list){
        if(cell.state!="dead" && cell.state!="newly_dead"){
            int i = cell.location[0];
            int j = cell.location[1];

            new_c8.push_back(cell);

            cg.allCells[i][j] = 1;
            cg.c8g[i][j] = 1;
            if(cell.state == "active"){
                cg.actT[i][j] = 1;
            }
        }
    }

    z = 0;
    std::vector<Macrophage> new_mp;
    for(auto & cell : mp_list){
        if(cell.state!="dead" && cell.state!="newly_dead"){
            int i = cell.location[0];
            int j = cell.location[1];

            cell.idx = z;

            new_mp.push_back(cell);

            cg.allCells[i][j] = 1;
            cg.mpg[i][j] = 1;
            cg.mid[i][j] = z;
            if(cell.state == "M0"){
                cg.m0g[i][j] = 1;
            }
            if(cell.state == "M1"){
                cg.m1g[i][j] = 1;
            }
            if(cell.state == "M2"){
                cg.m2g[i][j] = 1;
            }

            z++;
        }
    }

    cc_list = new_cancer;
    mp_list = new_mp;
    c8_list = new_c8;
}

void Environment::recruitMacrophage(CellGrids &cg, double rec) {
    // recruits macrophages to the simulation
    // randomly place them in available sites based on 
    // vasculature presence 

    std::mt19937 g(rd());
    std::uniform_real_distribution<double> dis(0.0,1.0);
    double p;

    for(int i=1; i<99; i++){
        for(int j=1; j<99; j++){
            p = dis(g);
            if(p<(rec*(macProb)*tstep*60*60) && cg.allCells[i][j]==0 && cg.vas[i][j]==1){
                cg.mpg[i][j] = 1;
                cg.allCells[i][j] = 1;
                cg.m0g[i][j] = 1;
                cg.mid[i][j] = mp_list.size();
                mp_list.push_back(Macrophage({i,j}, baseIFNR, baseIL4R, macrophageWeights, macrophageBiases, 0, mp_list.size()));
            }
        }
    }
}

void Environment::recruitTcells(int simStep, CellGrids &cg) {
    // taken from Gong 2017.
    // recruit T cells based on cancer cell deaths
    // delayed by tdelay days

    if(simStep*tstep/24 >= (tdelay + 0.5*twindow)){
        int Nc = 0;
        int delay = tdelay*24/tstep;

        for(int i=(-12*twindow/tstep); i<(12*twindow/tstep); i++){
            Nc += cancer_deaths[cancer_deaths.size() - 1 - delay - i];
        }

        double rate = r1*ka*Nc/((1/ki) + Nc);
        int num2rec = rate*tstep;
        if(rate>0 && num2rec==0){
            num2rec = 1;
        }

        std::vector<std::vector<int>> entry_points;

        int i, j;
        for(int q=0; q<num2rec; q++){
            i = 50;
            j = 50;

            std::uniform_real_distribution<> big(1.0,98.0);
            while(cg.vas[i][j] == 0){
                i = big(rd);
                j = big(rd);
            }
            entry_points.push_back({i,j});
        }

        for(auto & entry_point : entry_points){
            i = entry_point[0];
            j = entry_point[1];

            if(cg.allCells[i][j] == 0){
                cg.c8g[i][j] = 1;
                cg.allCells[i][j] = 1;
                c8_list.push_back(CD8({i,j}));
            }
        }
    }
}

void Environment::run_Macrophage(CellGrids &cg, Diffusibles &diff, double pi3k, double depletion) {
    // simulate macrophages present in simulation

    std::vector<int> run;
    for(int q=0; q<mp_list.size(); q++){
        run.push_back(q);
    }
    std::mt19937 g(rd());
    std::shuffle(run.begin(), run.end(), g);

    for(int l : run){
        mp_list[l].k15 = pi3k;
        mp_list[l].simulate(tstep, cg, diff, depletion);
    }
}

void Environment::run_CD8(CellGrids &cg) {
    // simulate T cells present in simulation

    std::vector<int> run;
    for(int q=0; q<c8_list.size(); q++){
        run.push_back(q);
    }
    std::mt19937 g(rd());
    std::shuffle(run.begin(), run.end(), g);

    for(int l : run){
        c8_list[l].simulate(tstep, cg, c8_list);
        if(c8_list[l].cc!=-1){attackedCancer.push_back(c8_list[l].cc);}
    }
}

void Environment::run_Cancer(CellGrids &cg){
    // simulate cancer cells present in simulation

    std::vector<int> run;
    for(int q=0; q<cc_list.size(); q++){
        run.push_back(q);
    }
    std::shuffle(run.begin(), run.end(), rd);

    int deaths = 0;

    for(int l : run){
        cc_list[l].simulate(attackedCancer, tstep, cg, cc_list);
        if(cc_list[l].state == "dead") {
            deaths++;
        }
    }
    cancer_deaths.push_back(deaths);
    attackedCancer.clear();
}

void Environment::removeDead(CellGrids &cg) {
    // remove any dead cells from the grids

    for(auto & cell : cc_list){
        if(cell.state=="dead"){
            int i = cell.location[0];
            int j = cell.location[1];

            cg.ccg[i][j] = 0;
            cg.allCells[i][j] = 0;
        }
    }

    for(auto & cell : c8_list){
        if(cell.state=="dead"){
            int i = cell.location[0];
            int j = cell.location[1];

            cg.c8g[i][j] = 0;
            cg.actT[i][j] = 0;
            cg.allCells[i][j] = 0;
        }
    }

    for(auto & cell : mp_list){
        if(cell.state=="dead"){
            int i = cell.location[0];
            int j = cell.location[1];

            cg.m0g[i][j] = 0;
            cg.m1g[i][j] = 0;
            cg.m2g[i][j] = 0;
            cg.allCells[i][j] = 0;
        }
    }
}

void Environment::initializeCells(CellGrids &cg) {
    // place initial cancer cells in center of environment
    // place initial macrophages randomly throughout the environment

    std::uniform_real_distribution<double> ilMod(-0.5,0.5);
    std::uniform_real_distribution<double> randType(0.0,1.0);
    int z = 0;
    totalCancerCells = 0;
    for(int i=48; i<53; i++){
        for(int j=48; j<53; j++){
            double kIL4c = baseIL4;//*pow(10.0, ilMod(rd));
            double kIFNGc = baseIFNG;//*pow(10.0, ilMod(rd));
            cg.ccg[i][j] = 1;
            cg.allCells[i][j] = 1;
            cg.ccid[i][j] = z;

            cc_list.push_back(Cancer({i,j}, z, kIL4c, kIFNGc, 1));
            z++;
            totalCancerCells++;
        }
    }

    std::mt19937 g(rd());
    std::uniform_real_distribution<> dis(0.0,1.0);

    z = 0;
    for(int i=1; i<99; i++){
        for(int j=1; j<99; j++){
            if(dis(g) < 2e-3 && cg.ccg[i][j]==0){
                cg.mpg[i][j] = 1;
                cg.allCells[i][j] = 1;
                cg.m0g[i][j] = 1;
                cg.mid[i][j] = z;
                mp_list.push_back(Macrophage({i,j}, baseIFNR, baseIL4R, macrophageWeights, macrophageBiases, 1, z));
                z++;
            }
        }
    }
}

void Environment::initializeTimeCourses(double days) {
    // set size of time course vectors

    for(int q=0; q<days*24/tstep; q++){
        times.push_back(0);
        cc_num.push_back(0);
        m0_num.push_back(0);
        m1_num.push_back(0);
        m2_num.push_back(0);
        c8_num.push_back(0);
        actT_num.push_back(0);
        M1f_max.push_back(0);
        M1f_avg.push_back(0);
        IL4_max.push_back(0);
        IL4_avg.push_back(0);
        TLS_max.push_back(0);
        TLS_avg.push_back(0);
        IFNG_max.push_back(0);
        IFNG_avg.push_back(0);
        volume.push_back(0);

        IFNGRExpression.push_back(0);
        IL4RExpression.push_back(0);
        SOCS1Expression.push_back(0);
        SOCS3Expression.push_back(0);
    }
}

void Environment::updateTimeCourses(int s, CellGrids &cg, Diffusibles &diff) {
    // for each time course, set the number of 
    // each cell type at step s

    int sumc = 0;
    int sum0 = 0;
    int sum1 = 0;
    int sum2 = 0;
    int sum8 = 0;

    double avgIL4 = 0;
    double avgAct = 0;
    double avgTLS = 0;
    double avgIFNG = 0;
    double maxIL4 = 0;
    double maxAct = 0;
    double maxTLS = 0;
    double maxIFNG = 0;

    for(int i=0;i<100;i++){
        for(int j=0;j<100;j++){
            sumc += cg.ccg[i][j];
            sum0 += cg.m0g[i][j];
            sum1 += cg.m1g[i][j];
            sum2 += cg.m2g[i][j];
            sum8 += cg.c8g[i][j];
            avgIL4 += diff.IL4[i][j];
            avgAct += diff.M1f[i][j];
            avgIFNG += diff.IFNG[i][j];

            if(diff.IL4[i][j] > maxIL4){maxIL4=diff.IL4[i][j];}
            if(diff.M1f[i][j] > maxAct){maxAct=diff.M1f[i][j];}
            if(diff.IFNG[i][j] > maxIFNG){maxIFNG=diff.IFNG[i][j];}

            if(cg.c8g[i][j]==0 && cg.actT[i][j]==1){
                std::cout << "-------\n";
                std::cout << "Step: " << s << std::endl;
                std::cout << "Num t cells: " << c8_list.size() << std::endl;
                std::cout << "i j: " << i << " " << j << std::endl;
                std::cout << "actT: " << cg.actT[i][j] << std::endl;
                std::cout << "c8g: " <<  cg.c8g[i][j] << std::endl;
                std::cout << cg.allCells[i][j] << std::endl;
                throw std::runtime_error("done fucked up");
            }
        }
    }

    avgIL4 = avgIL4/(100*100);
    avgAct = avgAct/(100*100);
    avgTLS = avgTLS/(100*100);
    avgIFNG = avgIFNG/(100*100);

    /*times[s] = (tstep*(s+1)/24);
    cc_num[s] = (sumc);
    m0_num[s] = (sum0);
    m1_num[s] = (sum1);
    m2_num[s] = (sum2);
    c8_num[s] = (sum8);

    M1f_max[s] = maxAct;
    M1f_avg[s] = avgAct;
    IL4_max[s] = maxIL4;
    IL4_avg[s] = avgIL4;
    TLS_max[s] = maxTLS;
    TLS_avg[s] = avgTLS;
    IFNG_avg[s] = avgIFNG;
    IFNG_max[s] = maxIFNG;*/

    times.push_back(tstep*(s+1)/24);
    cc_num.push_back(sumc);
    m0_num.push_back(sum0);
    m1_num.push_back(sum1);
    m2_num.push_back(sum2);
    c8_num.push_back(sum8);

    M1f_max.push_back(maxAct);
    M1f_avg.push_back(avgAct);
    IL4_max.push_back(maxIL4);
    IL4_avg.push_back(avgIL4);
    TLS_max.push_back(maxTLS);
    TLS_avg.push_back(avgTLS);
    IFNG_avg.push_back(avgIFNG);
    IFNG_max.push_back(maxIFNG);

    double tumorVolume = sumc*(dx*dx); // cm^2
    //volume[s] = tumorVolume;
    volume.push_back(tumorVolume);

    int activeTcells = 0;
    for(auto & cell : c8_list){
        if(cell.state == "active"){
            activeTcells++;
        }
    }

    //actT_num[s] = activeTcells;
    actT_num.push_back(activeTcells);

    //int totalDeaths=0;
    //for(auto & q:cancer_deaths){totalDeaths+=q;}

    endTime = (s+1)*(tstep/24);
}

void Environment::checkError(int s, CellGrids &cg) {
    // make sure that all of the cell counts
    // equal the correct values
    // and that no cells are missing from the grid
    // or occupying the same space

    int sumc = cc_num[s];
    int sum8 = c8_num[s];
    int sum0 = m0_num[s];
    int sum1 = m1_num[s];
    int sum2 = m2_num[s];

    int allSum = 0;
    for(int i=1; i<99; ++i){
        for(int j=1; j<99; ++j){
            if(cg.allCells[i][j] > 1){
                throw std::runtime_error("You have forgotten the face of your father.");
            }
            allSum += cg.allCells[i][j];
        }
    }


    int macrophages = 0;
    for(auto & cell : mp_list){
        if(cell.state=="M0" || cell.state=="M1" || cell.state=="M2"){
            macrophages++;
        }
    }

    int tcells = 0;
    for(auto & cell : c8_list){
        if(cell.state=="active" || cell.state=="inactive" || cell.state=="exhausted"){
            tcells++;
        }
    }

    int cancercels = 0;
    for(auto & cell : cc_list){
        if(cell.state == "alive"){
            cancercels++;
        }
    }

    if(cancercels != sumc){
        std::cout << cancercels << " " << sumc << " " << cc_list.size() << " " << allSum << std::endl;
        throw std::runtime_error("Cancer cells not equal");
    }
    if(tcells != sum8){
        throw std::runtime_error("T Cells not equal");
    }
    if(macrophages != (sum0 + sum1 + sum2)){
        std::cout << macrophages << " " << sum0 << " " << sum1 << " " << sum2 << std::endl;
        throw std::runtime_error("Macros not equal");
    }
    if(sumc+sum0+sum1+sum2+sum8 != allSum){
        std::cout << "allSum: " << allSum << std::endl;
        std::cout << "ccg: " << sumc << std::endl;
        std::cout << "m0g: " << sum0 << std::endl;
        std::cout << "macros: " << macrophages << std::endl;
        std::cout << s << std::endl;
        std::cout << std::endl;

        for(int i=1; i<99; ++i){
            for(int j=1; j<99; ++j){
                if(cg.allCells[i][j]==1 && cg.ccg[i][j]+cg.m0g[i][j]!=1){
                    std::cout << i << " " << j << " " << cg.ccg[i][j]+cg.m0g[i][j] <<  std::endl;
                    for(auto & cell : cc_list){
                        if(cell.location[0]==i && cell.location[1]==j){
                            std::cout << "cancer!\n";
                        }
                    }
                    for(auto & cell : mp_list){
                        if(cell.location[0]==i && cell.location[1]==j){
                            std::cout << "macro!\n";
                        }
                    }
                }
            }
        }
        throw std::runtime_error("Cry your pardon.");
    }
}

void Environment::printStep() {
    // print current state of the simulation
    int s = cc_num.size() - 1;

    std::cout<< "-------------------------\n"
             << "Sim time: " << (s+1)*tstep/24 << "\n"
             << "Num cancer cells: " << cc_num[s]  << " " << cc_list.size() << "\n"
             << "Num M0: " << m0_num[s] << "\n"
             << "Num M1: " << m1_num[s] << "\n"
             << "Num M2: " << m2_num[s] << "\n"
             << "Total M: " << m0_num[s] + m1_num[s] + m2_num[s] << " " << mp_list.size() << "\n"
             << "Num CD8: " << c8_num[s] << " " << c8_list.size() << "\n"
             << "Active CD8: " << actT_num[s] << "\n"
             << "Avg/Max Act: " << M1f_avg[s] << " | " << M1f_max[s] << "\n"
             << "Avg/Max TLS: " << TLS_avg[s] << " | " << TLS_max[s] << "\n"
             << "Avg/Max IFNG: " << IFNG_avg[s] << " | " << IFNG_max[s] << "\n"
             << "Avg/Max IL4: " << IL4_avg[s] << " | " << IL4_max[s] << "\n"
             << "Cancer Deaths: " << cancer_deaths[cancer_deaths.size() - 1] << "\n";
}

void Environment::clearGrids(CellGrids &cg, Diffusibles &diff) {
    // resets all the global variables
    // I probably shouldn't use those, but I'll deal with that some other time

    cc_list.clear();
    mp_list.clear();
    c8_list.clear();

    for(int i=0; i<100; i++){
        for(int j=0; j<100; j++){
            cg.ccg[i][j] = 0;
            cg.ccid[i][j] = -1;
            cg.mpg[i][j] = 0;
            cg.m0g[i][j] = 0;
            cg.m1g[i][j] = 0;
            cg.m2g[i][j] = 0;
            cg.c8g[i][j] = 0;
            cg.actT[i][j] = 0;
            cg.allCells[i][j] = 0;
            cg.vas[i][j] = 1;

            diff.M1f[i][j] = 0;
            diff.IL4[i][j] = 0;
            diff.IFNG[i][j] = 0;
        }
    }
}

void Environment::printFinalInfo(CellGrids &cg) {
    // print final cell counts

    int ccs = 0;
    int m1s = 0;
    int m2s = 0;
    int c8s = 0;
    int acts = 0;

    for(int i=1; i<99; i++){
        for(int j=1; j<99; j++){
            ccs += cg.ccg[i][j];
            m1s += cg.m1g[i][j];
            m2s += cg.m2g[i][j];
            c8s += cg.c8g[i][j];
            acts += cg.actT[i][j];
        }
    }

    std::cout << endTime  << " " << ccs << " " << m1s << " " << m2s << " " << c8s << " " << acts << std::endl;
}

void Environment::saveVideo(int s, CellGrids &cg) {
    std::ofstream myfile;
    myfile.open(saveDir+"/cells/cells_"+std::to_string(s)+".csv");
    for(int i=0; i<100; ++i){
        myfile << cg.ccg[i][0] + 2*cg.c8g[i][0] + 3*cg.m0g[i][0] + 4*cg.m1g[i][0] + 5*cg.m2g[i][0] + 6*cg.actT[i][0];
        for(int j=1; j<100; ++j){
            myfile << "," << cg.ccg[i][j] + 2*cg.c8g[i][j] + 3*cg.m0g[i][j] + 4*cg.m1g[i][j] + 5*cg.m2g[i][j] + 6*cg.actT[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();
}

void  Environment::macroInfo(int s) {
    double avgIFNGR = 0;
    double avgIL4R = 0;
    double avgSOCS1 = 0;
    double avgSOCS3 = 0;

    for(auto & m: mp_list){
        avgIFNGR += m.IFNGR0;
        avgIL4R += m.IL4R0;
        avgSOCS1 += m.SOCS10;
        avgSOCS3 += m.SOCS30;
    }

    avgIFNGR /= mp_list.size();
    avgIL4R /= mp_list.size();
    avgSOCS1 /= mp_list.size();
    avgSOCS3 /= mp_list.size();

    /*IFNGRExpression[s] = avgIFNGR;
    IL4RExpression[s] = avgIL4R;
    SOCS1Expression[s] = avgSOCS1;
    SOCS3Expression[s] = avgSOCS3;*/

    IFNGRExpression.push_back(avgIFNGR);
    IL4RExpression.push_back(avgIL4R);
    SOCS1Expression.push_back(avgSOCS1);
    SOCS3Expression.push_back(avgSOCS3);

    /*std::ofstream myfile;
    myfile.open(saveDir+"/mIFNGR0.csv", std::ios_base::app);
    for(auto & m : allMacro){
        myfile << m.IFNGR0 << ",";
    }
    myfile << std::endl;
    myfile.close();

    myfile.open(saveDir+"/mIL4R0.csv", std::ios_base::app);
    for(auto & m : allMacro){
        myfile << m.IL4R0 << ",";
    }
    myfile << std::endl;
    myfile.close();

    myfile.open(saveDir+"/mSOCS10.csv", std::ios_base::app);
    for(auto & m : allMacro){
        myfile << m.IL4R0 << ",";
    }
    myfile << std::endl;
    myfile.close();

    myfile.open(saveDir+"/mSOCS30.csv", std::ios_base::app);
    for(auto & m : allMacro){
        myfile << m.IL4R0 << ",";
    }
    myfile << std::endl;
    myfile.close();*/
}

double Environment::treatmentLevel(double treatmentModulation, double timeOn) {
    // if treatmentModulation = 0, assume always on
    if(treatmentModulation == 0){
        return 1;
    }

    if(timeOn >= treatmentModulation){
        return 1;
    }

    if(treatmentOn > 0){
        treatmentOn -= tstep;
        if(treatmentOn <= 0){treatmentOff = treatmentModulation - timeOn;}
        //std::cout << "On\n";
        return 1;
    }
    if(treatmentOff > 0){
        treatmentOff -= tstep;
        if(treatmentOff <= 0){treatmentOn = timeOn;}
        //std::cout << "Off\n";
        return 0;
    }

    if((treatmentOn<=0 && treatmentOff<=0)){
        throw std::runtime_error("Environment::treatmentLevel");
    }

    throw std::runtime_error("end of treatmentLevel");
}

void Environment::simulate(double days, double treatmentModulation, double timeOn){
    // run a simulation

    // time on is the percent of treatmentModulation that
    // treatment is on for
    //treatmentOn = treatmentModulation*24*timeOn;
    treatmentOn = timeOn*24;

    //create grids for cells and diffusibles
    CellGrids cg;
    Diffusibles diff(dx, D, baseIL4, 10, k_IFNG_T, k_M1f);
    //Diffusibles diff(dx, D, baseIL4, paramsLHS[4], k_IFNG_T, k_M1f);

    // initialize cell number lists
    //initializeTimeCourses(days);

    // place initial cells
    initializeCells(cg);

    // treatment targets
    double pi3k = 1.16;
    double depletion = 0;
    double rec = 1;

    for(int s=0; s<days*24/tstep; s++){
        if(s*tstep/24 >= PI3KTime && PI3Klvl > 0){pi3k = 1.16*(1 - PI3Klvl*treatmentLevel(treatmentModulation*24, 24*timeOn));}
        if(s*tstep/24 >= depTime && deplvl > 0){depletion = deplvl*treatmentLevel(treatmentModulation*24, 24*timeOn);}
        if(s*tstep/24 >= recTime && reclvl > 0){rec = 1*(1 - reclvl*treatmentLevel(treatmentModulation*24, 24*timeOn));}

        diff.diffusion(cg, tstep);
        recruitMacrophage(cg, rec);
        run_Macrophage(cg, diff, pi3k, depletion);
        run_CD8(cg);
        recruitTcells(s, cg);
        run_Cancer(cg);
        removeDead(cg);
        clean(cg);

        updateTimeCourses(s, cg, diff);
        checkError(s, cg);

        //printStep();
        //if(fmod(s*tstep, 24) == 0 && (s*tstep/24)>=100){
        //    plot(cg, diff, s);
        //}

        if(cc_num[s] == 0 || cc_num[s] > 5000){break;}
    }
    //plot();
    save();
    //printFinalInfo(cg);
    clearGrids(cg, diff);
    paramsLHS.clear();
    //std::cout <<"Unique Cancer Cells: " << totalCancerCells <<  std::endl;
}

/*
 * REFERENCES
 * ----------
 * Han. "Multidimensional analysis of the frequencies and rates of cytokine secretion from single cells by quantitative microengraving". Lab Chip. 2010
 *
 * Hao. "Size-based separation methods of circulating tumor cells". Advanced Drug Delivery Reviews. 2018.
 *
 * Wells. "Spatial and Functional Heterogeneities Shape Collective Behavior of Tumor-Immune Networks". PLoS Comp Bio. 2015.
 *
 */

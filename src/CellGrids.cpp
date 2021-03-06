//
// Created by Colin Cess on 3/11/20.
//

#include "CellGrids.h"

CellGrids::CellGrids() {
    for(int i=0; i<100; ++i){
        for(int j=0; j<100; ++j){
            ccg[i][j] = 0;
            ccid[i][j] = -1;
            m0g[i][j] = 0;
            m1g[i][j] = 0;
            m2g[i][j] = 0;
            mid[i][j] = -1;
            c8g[i][j] = 0;
            actT[i][j] = 0;
            allCells[i][j] = 0;
            vas[i][j] = 1;
        }
    }

    // set boundries
    for(int i=0;i<100;i++){
        for(int j=0;j<100;j++){
            allCells[i][0] = 1;
            allCells[i][99] = 1;
            allCells[0][j] = 1;
            allCells[99][j] = 1;
        }
    }

    // set vasculature (entry points for cells)
    for(int i=30;i<70;i++){
        for(int j=30;j<70;j++){
            for(int k=30;k<70;k++){
                vas[i][j] = 0;
            }
        }
    }
}
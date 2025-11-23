#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <algorithm>
#include <array>
#include <cstdlib> 
#include <chrono>
#include <iomanip>

enum ConservedVariables { rho_p, rhovx, rhovy, rhovz, E, Bx_p, By_p, Bz_p };
enum PrimitiveVariables { rho, vx, vy, vz, p, Bx, By, Bz } ;
enum Parameters {xStart, xEnd, yStart, yEnd, tStart, tEnd, xCells, yCells, CFL, g };

//Constants and parameters storage function
std::array<double, 10> getPars(){
    std::array<double, 10> parameters;
    parameters[xStart] = 0;
    parameters[xEnd] = 1;
    parameters[yStart] = 0;
    parameters[yEnd] = 1;    
    parameters[tStart] = 0;
    parameters[tEnd] = 1;
    parameters[xCells] = 256;
    parameters[yCells] = 256;
    parameters[CFL] = 0.8;
    parameters[g] = 5.0/3.0;
    return parameters;
}

// Primitive to conservative and conservative to primitive functions

std::array<double, 8> PrimToCons(std::array<double, 8> prim){
    //prim and cons here are arrays which store the variables at a specific cell i
    std::array<double, 8> cons;
    double gamma = getPars()[g];

    double rho = prim[0];
    double vx = prim[1];
    double vy = prim[2];
    double vz = prim[3];
    double p = prim[4];
    double Bx = prim[5];
    double By = prim[6];
    double Bz = prim[7];
    double v_squared = vx*vx + vy*vy + vz*vz;
    double B_squared = Bx*Bx + By*By + Bz*Bz;

    cons[0] = rho; 
    cons[1] = rho*vx;
    cons[2] = rho*vy;
    cons[3] = rho*vz;
    cons[4] = p/(gamma-1) + 0.5*rho*v_squared + 0.5*B_squared; //Energy
    cons[5] = Bx;
    cons[6] = By;
    cons[7] = Bz;

    return cons;
}

std::array<double, 8> ConsToPrim(std::array<double, 8> cons){
    //prim and cons here are arrays which store the variables at a specific cell i
    std::array<double, 8> prim;
    double gamma = getPars()[g];

    double rho = cons[0];
    double rhoVx = cons[1];
    double rhoVy = cons[2];
    double rhoVz = cons[3];
    double E = cons[4];
    double Bx = cons[5];
    double By = cons[6];
    double Bz = cons[7];
    /*if (rho <= 0.0) {
        std::cout<<"Unphysical state detected: rho is non-positive."<<std::endl;
    }*/
    double vx = rhoVx/rho;
    double vy = rhoVy/rho;
    double vz = rhoVz/rho;
    double v_squared = vx*vx + vy*vy + vz*vz;
    double B_squared = Bx*Bx + By*By + Bz*Bz;
    double p = (E - 0.5*rho*v_squared - 0.5*B_squared) * (gamma-1) ;

    prim[0] = rho;
    prim[1] = vx;
    prim[2] = vy;
    prim[3] = vz;
    prim[4] = p;
    prim[5] = Bx;
    prim[6] = By;
    prim[7] = Bz;
    return prim;
}

// x-dimension functions

std::array<double, 8> getXFluxFunction(std::array<double,8> cons) {
    std::array<double, 8> fluxFunction;
    std::array<double, 8> prim = ConsToPrim(cons);
    double rho = prim[0];
    double vx = prim[1];
    double vy = prim[2];
    double vz = prim[3];
    double p = prim[4];
    double Bx = prim[5];
    double By = prim[6];
    double Bz = prim[7];
    double B_squared = Bx*Bx + By*By + Bz*Bz;
    double E = cons[4];

    /*if (rho <= 0.0) {
        std::cout<<"Unphysical state detected: rho is non-positive in getFluxFunction."<<std::endl;
    }
    if (p < 0.0) {
        std::cout<<"Unphysical state detected: p is non-positive in getFluxFunction."<<std::endl;
    }*/

    // Compute fluxes
    fluxFunction[0] = rho*vx;  
    fluxFunction[1] = rho*vx*vx + p + 0.5*B_squared - Bx*Bx;
    fluxFunction[2] = rho*vx*vy - Bx*By;
    fluxFunction[3] = rho*vx*vz - Bx*Bz;
    fluxFunction[4] = (E + p + 0.5*B_squared)*vx - Bx*(Bx*vx + By*vy + Bz*vz);
    fluxFunction[5] = 0.0;  //is Bxvx - vxBx
    fluxFunction[6] = By*vx - Bx*vy;
    fluxFunction[7] = Bz*vx - Bx*vz;

    return fluxFunction;
}

std::array<double, 8> getXFluxLF(std::array<double, 8> u_i, std::array<double, 8>  u_iPlus1, double dt){
    //NOTE: u_iPlus1 here is u(i+1,n), NOT u(i,n+1)
    int nxCells = getPars()[xCells];
    double x0 = getPars()[xStart];
    double x1 = getPars()[xEnd];
    double dx = (x1-x0)/nxCells;
    std::array<double, 8> f_iPlusHalf;
    std::array<double,8> flux_ui = getXFluxFunction(u_i);
    std::array<double,8> flux_uiPlus1 = getXFluxFunction(u_iPlus1);
    for (int v = 0; v < 8; v++){
        f_iPlusHalf[v] = 0.5*(dx/dt)*(u_i[v] - u_iPlus1[v]) + 0.5*(flux_ui[v] + flux_uiPlus1[v]);
    }
    return f_iPlusHalf;
}

std::array<double, 8> getXFluxRI(std::array<double,8> u_i, std::array<double,8> u_iPlus1, double dt){
    int nxCells = getPars()[xCells];
    double x0 = getPars()[xStart];
    double x1 = getPars()[xEnd];
    double dx = (x1-x0)/nxCells;
    std::array<double, 8> u_iPlusHalf;
    std::array<double, 8> f_iPlusHalf;
    std::array<double, 8> flux_ui = getXFluxFunction(u_i);
    std::array<double, 8> flux_uiPlus1 = getXFluxFunction(u_iPlus1);   

    for (int v = 0; v < 8; v++){
        u_iPlusHalf[v] = 0.5*(u_i[v] + u_iPlus1[v]) - 0.5*(dt/dx)*(flux_uiPlus1[v] - flux_ui[v]);
    }

    std::array<double, 8> flux_u_iPlusHalf = getXFluxFunction(u_iPlusHalf);

    for (int v = 0; v < 8; v++){
        f_iPlusHalf[v] = flux_u_iPlusHalf[v];
    }

    return f_iPlusHalf;
}

std::array<double, 8> getXfluxFORCE(std::array<double,8> u_i, std::array<double,8> u_iPlus1, double dt){
    std::array<double, 8> f_iPlusHalf_FORCE;
    std::array<double, 8> fluxLF = getXFluxLF(u_i, u_iPlus1, dt);
    std::array<double, 8> fluxRI = getXFluxRI(u_i, u_iPlus1, dt);    
    for (int i = 0; i < 8; i++){
        f_iPlusHalf_FORCE[i] = 0.5*fluxLF[i] + 0.5*fluxRI[i];
    }
    return f_iPlusHalf_FORCE;
}

std::array<double, 8> getXuHLL(std::array<double,8> uL, std::array<double,8> uR){
    std::array<double,8> uHLL;
    std::array<double,8> uprimL = ConsToPrim(uL);
    std::array<double,8> uprimR = ConsToPrim(uR);
    double gamma = getPars()[g];
    //Calculate cfL and cfR 
    double csL = sqrt(gamma*uprimL[p]/uprimL[rho_p]);
    double BxL = uprimL[Bx];
    double BL_squared = pow(uprimL[Bx], 2) + pow(uprimL[By], 2) + pow(uprimL[Bz], 2);
    double rhoL = uprimL[rho_p];
    double cfL = sqrt(0.5*(csL*csL + BL_squared/rhoL) + sqrt((pow(csL*csL + BL_squared/rhoL, 2) - 4*csL*csL*BxL*BxL/rhoL)));

    double csR = sqrt(gamma*uprimR[p]/uprimR[rho_p]);
    double BxR = uprimR[Bx];
    double BR_squared = pow(uprimR[Bx], 2) + pow(uprimR[By], 2) + pow(uprimR[Bz], 2);
    double rhoR = uprimR[rho_p];
    double cfR = sqrt(0.5*((csR*csR + BR_squared/rhoR) + sqrt(pow((csR*csR + BR_squared/rhoR), 2) - 4*csR*csR*BxR*BxR/rhoR)));
    //Use cfL and cfR to get sL and sR
    double sL = std::min(uprimL[vx],uprimR[vx]) - std::max(cfL, cfR);
    double sR = std::max(uprimL[vx],uprimR[vx]) + std::max(cfL, cfR);   
    for (int i = 0; i < 8; i++){
        uHLL[i] = (sR*uR[i] - sL*uL[i] + getXFluxFunction(uL)[i] - getXFluxFunction(uR)[i])/(sR-sL);
    }
    return uHLL;
}

std::array<double, 8> getXfluxHLL(std::array<double,8> uL, std::array<double,8> uR){
    std::array<double,8> uHLL;
    std::array<double,8> fHLL;
    std::array<double,8> f_iPlusHalf;
    std::array<double,8> uprimL = ConsToPrim(uL);
    std::array<double,8> uprimR = ConsToPrim(uR);
    double gamma = getPars()[g];

    double csL = sqrt(gamma*uprimL[p]/uprimL[rho_p]);
    double BxL = uprimL[Bx];
    double BL_squared = pow(uprimL[Bx], 2) + pow(uprimL[By], 2) + pow(uprimL[Bz], 2);
    double rhoL = uprimL[rho_p];
    double cfL = sqrt(0.5*(csL*csL + BL_squared/rhoL) + sqrt((pow(csL*csL + BL_squared/rhoL, 2) - 4*csL*csL*BxL*BxL/rhoL)));

    double csR = sqrt(gamma*uprimR[p]/uprimR[rho_p]);
    double BxR = uprimR[Bx];
    double BR_squared = pow(uprimR[Bx], 2) + pow(uprimR[By], 2) + pow(uprimR[Bz], 2);
    double rhoR = uprimR[rho_p];
    double cfR = sqrt(0.5*((csR*csR + BR_squared/rhoR) + sqrt(pow((csR*csR + BR_squared/rhoR), 2) - 4*csR*csR*BxR*BxR/rhoR)));
    double sL = std::min(uprimL[vx],uprimR[vx]) - std::max(cfL, cfR);
    double sR = std::max(uprimL[vx], uprimR[vx]) + std::max(cfL, cfR);
    //std::cout << "X HLL: " << "sL = " << sL << ", sR = " << sR << std::endl;
    for (int i = 0; i < uHLL.size(); i++){
        uHLL[i] = (sR*uR[i] - sL*uL[i] + getXFluxFunction(uL)[i] - getXFluxFunction(uR)[i])/(sR-sL);
        fHLL[i] = (sR*getXFluxFunction(uL)[i] - sL*getXFluxFunction(uR)[i] + sL*sR*(uR[i] - uL[i]))/(sR-sL);
        if (0 <= sL){
            f_iPlusHalf[i] = getXFluxFunction(uL)[i];
        }
        if (sL < 0 && 0 < sR){
            f_iPlusHalf[i] = fHLL[i];
        }
        if (0 >= sR){
            f_iPlusHalf[i] = getXFluxFunction(uR)[i];
        }
    }
    return f_iPlusHalf;
}

std::array<double, 8> getXfluxHLLC(std::array<double,8> uL, std::array<double,8> uR){
    std::array<double,8> uHLL = getXuHLL(uL, uR);
    std::array<double,8> f_iPlusHalf_HLLC;
    std::array<double,8> uStarL;
    std::array<double,8> uStarR;
    std::array<double,8> uprimL = ConsToPrim(uL);
    std::array<double,8> uprimR = ConsToPrim(uR);
    double gamma = getPars()[g];
    
    double csL = sqrt(gamma*uprimL[p]/uprimL[rho_p]);
    double BxL = uprimL[Bx];
    double BL_squared = pow(uprimL[Bx], 2) + pow(uprimL[By], 2) + pow(uprimL[Bz], 2);
    double rhoL = uprimL[rho_p];
    double cfL = sqrt(0.5*(csL*csL + BL_squared/rhoL) + sqrt((pow(csL*csL + BL_squared/rhoL, 2) - 4*csL*csL*BxL*BxL/rhoL)));

    double csR = sqrt(gamma*uprimR[p]/uprimR[rho_p]);
    double BxR = uprimR[Bx];
    double BR_squared = pow(uprimR[Bx], 2) + pow(uprimR[By], 2) + pow(uprimR[Bz], 2);
    double rhoR = uprimR[rho_p];
    double cfR = sqrt(0.5*((csR*csR + BR_squared/rhoR) + sqrt(pow((csR*csR + BR_squared/rhoR), 2) - 4*csR*csR*BxR*BxR/rhoR)));
    double vxR = uprimR[vx];
    double vxL = uprimL[vx];
    double sL = std::min(uprimL[vx],uprimR[vx]) - std::max(cfL, cfR);
    double sR = std::max(uprimL[vx],uprimR[vx]) + std::max(cfL, cfR);
    //double sL = std::min(std::abs(uprimL[vx]),std::abs(uprimR[vx])) - std::max(std::abs(cfL), std::abs(cfR));
    //double sR = std::max(std::abs(uprimL[vx]),std::abs(uprimR[vx])) + std::max(std::abs(cfL), std::abs(cfR));
    double ptL = uprimL[p] + 0.5*BL_squared;
    double ptR = uprimR[p] + 0.5*BR_squared;
    double sStar = (rhoR*vxR*(sR - vxR) - rhoL*vxL*(sL - vxL) + ptL - ptR - BxL*BxL + BxR*BxR)\
    /(rhoR*(sR - vxR) - rhoL*(sL-vxL));
    //std::cout << "X HLL: " << "sL = " << sL << ", sR = " << sR << "s* = " << sStar <<std::endl;

    uStarL[rho] = rhoL*(sL - vxL)/(sL - sStar);
    uStarR[rho] = rhoR*(sR - vxR)/(sR - sStar);
    uStarL[rhovx] = uStarL[rho]*sStar;
    uStarR[rhovx] = uStarR[rho]*sStar;
    double ByL = uL[By];
    double ByR = uR[By];
    double BzL = uL[Bz];
    double BzR = uR[Bz];
    uStarL[rhovy] = uL[rhovy]*(sL - vxL)/(sL - sStar) - (uHLL[Bx]*uHLL[By] - BxL*ByL)/(sL - sStar);
    uStarR[rhovy] = uR[rhovy]*(sR - vxR)/(sR - sStar) - (uHLL[Bx]*uHLL[By] - BxR*ByR)/(sR - sStar);
    uStarL[rhovz] = uL[rhovz]*(sL - vxL)/(sL - sStar) - (uHLL[Bx]*uHLL[Bz] - BxL*BzL)/(sL - sStar);
    uStarR[rhovz] = uR[rhovz]*(sR - vxR)/(sR - sStar) - (uHLL[Bx]*uHLL[Bz] - BxR*BzR)/(sR - sStar);
    uStarR[Bx] = uHLL[Bx];
    uStarL[Bx] = uHLL[Bx];
    uStarR[By] = uHLL[By];
    uStarL[By] = uHLL[By];
    uStarR[Bz] = uHLL[Bz];
    uStarL[Bz] = uHLL[Bz];
    double vxHLL = uHLL[rhovx]/uHLL[rho];
    double vyHLL = uHLL[rhovy]/uHLL[rho];
    double vzHLL = uHLL[rhovz]/uHLL[rho];
    double Bv_HLL = uHLL[Bx]*vxHLL + uHLL[By]*vyHLL + uHLL[Bz]*vzHLL;
    double Bv_R = uprimR[Bx]*uprimR[vx] + uprimR[By]*uprimR[vy] + uprimR[Bz]*uprimR[vz];
    double Bv_L = uprimL[Bx]*uprimL[vx] + uprimL[By]*uprimL[vy] + uprimL[Bz]*uprimL[vz];
    double ptStarL = rhoL*(sL - vxL)*(sStar - vxL) + ptL - BxL*BxL + uHLL[Bx]*uHLL[Bx];
    double ptStarR = rhoR*(sR - vxR)*(sStar - vxR) + ptR - BxR*BxR + uHLL[Bx]*uHLL[Bx];
    uStarR[E] = uR[E]*(sR - vxR)/(sR - sStar) + (ptStarR*sStar - ptR*vxR - (uHLL[Bx]*(Bv_HLL)- BxR*Bv_R))/(sR - sStar);
    uStarL[E] = uL[E]*(sL - vxL)/(sL - sStar) + (ptStarL*sStar - ptL*vxL - (uHLL[Bx]*(Bv_HLL)- BxL*Bv_L))/(sL - sStar);
        for (int i = 0; i < 8; i++){
        if (0 <= sL){
            f_iPlusHalf_HLLC[i] = getXFluxFunction(uL)[i];
        }
        if (sL < 0 && 0 <= sStar){
            f_iPlusHalf_HLLC[i] = getXFluxFunction(uL)[i] + sL*(uStarL[i] - uL[i]);
        }
        if (sStar < 0 && 0 <= sR){
            f_iPlusHalf_HLLC[i] = getXFluxFunction(uR)[i] + sR*(uStarR[i] - uR[i]);
        }
        if (0 > sR){
            f_iPlusHalf_HLLC[i] = getXFluxFunction(uR)[i];
        }
    }
    return f_iPlusHalf_HLLC;

}

// y-dimension functions

std::array<double, 8> getYFluxFunction(std::array<double,8> cons) {
    std::array<double, 8> fluxFunction;
    std::array<double, 8> prim = ConsToPrim(cons);
    double rho = prim[0];
    double vx = prim[1];
    double vy = prim[2];
    double vz = prim[3];
    double p = prim[4];
    double Bx = prim[5];
    double By = prim[6];
    double Bz = prim[7];
    double B_squared = Bx*Bx + By*By + Bz*Bz;
    double E = cons[4];

    /*if (rho <= 0.0) {
        std::cout<<"Unphysical state detected: rho is non-positive in getFluxFunction."<<std::endl;
    }
    if (p < 0.0) {
        std::cout<<"Unphysical state detected: p is non-positive in getFluxFunction."<<std::endl;
    }*/

    // Compute fluxes
    fluxFunction[0] = rho*vy;  
    fluxFunction[1] = rho*vx*vy - By*Bx;
    fluxFunction[2] = rho*vy*vy + p + 0.5*B_squared - By*By ;
    fluxFunction[3] = rho*vz*vy - By*Bz;
    fluxFunction[4] = (E + p + 0.5*B_squared)*vy - By*(Bx*vx + By*vy + Bz*vz);
    fluxFunction[5] = Bx*vy - By*vx;
    fluxFunction[6] = 0 ; // is Byvy - vyBy
    fluxFunction[7] = Bz*vy - By*vz;

    return fluxFunction;
}

std::array<double, 8> getYFluxLF(std::array<double, 8> u_i, std::array<double, 8>  u_iPlus1, double dt){
    //NOTE: u_iPlus1 here is u(i+1,n), NOT u(i,n+1)
    int nyCells = getPars()[yCells];
    double y0 = getPars()[yStart];
    double y1 = getPars()[yEnd];
    double dy = (y1-y0)/nyCells;
    std::array<double, 8> f_iPlusHalf;
    std::array<double,8> flux_ui = getYFluxFunction(u_i);
    std::array<double,8> flux_uiPlus1 = getYFluxFunction(u_iPlus1);
    for (int v = 0; v < 8; v++){
        f_iPlusHalf[v] = 0.5*(dy/dt)*(u_i[v] - u_iPlus1[v]) + 0.5*(flux_ui[v] + flux_uiPlus1[v]);
    }
    return f_iPlusHalf;
}

std::array<double, 8> getYFluxRI(std::array<double,8> u_i, std::array<double,8> u_iPlus1, double dt){
    int nyCells = getPars()[yCells];
    double y0 = getPars()[yStart];
    double y1 = getPars()[yEnd];
    double dy = (y1-y0)/nyCells;
    std::array<double, 8> u_iPlusHalf;
    std::array<double, 8> f_iPlusHalf;
    std::array<double, 8> flux_ui = getYFluxFunction(u_i);
    std::array<double, 8> flux_uiPlus1 = getYFluxFunction(u_iPlus1);   

    for (int v = 0; v < 8; v++){
        u_iPlusHalf[v] = 0.5*(u_i[v] + u_iPlus1[v]) - 0.5*(dt/dy)*(flux_uiPlus1[v] - flux_ui[v]);
    }
    
    std::array<double, 8> flux_u_iPlusHalf = getYFluxFunction(u_iPlusHalf);

    for (int v = 0; v < 8; v++){
        //f_iPlusHalf[v] = getYFluxFunction(u_iPlusHalf)[v];
        f_iPlusHalf[v] = flux_u_iPlusHalf[v];
    }

    return f_iPlusHalf;
}

std::array<double, 8> getYfluxFORCE(std::array<double,8> u_i, std::array<double,8> u_iPlus1, double dt){
    std::array<double, 8> f_iPlusHalf_FORCE;
    std::array<double ,8> fluxLF = getYFluxLF(u_i, u_iPlus1, dt);
    std::array<double ,8> fluxRI = getYFluxRI(u_i, u_iPlus1, dt);    
    for (int var = 0; var < 8; var++){
        f_iPlusHalf_FORCE[var] = 0.5*fluxLF[var] + 0.5*fluxRI[var];
    }
    return f_iPlusHalf_FORCE;
}

std::array<double, 8> getYuHLL(std::array<double,8> uL, std::array<double,8> uR){
    std::array<double,8> uHLL;
    std::array<double,8> uprimL = ConsToPrim(uL);
    std::array<double,8> uprimR = ConsToPrim(uR);
    double gamma = getPars()[g];
    //Calculate cfL and cfR 
    double csL = sqrt(gamma*uprimL[p]/uprimL[rho_p]);
    double ByL = uprimL[By];
    double BL_squared = pow(uprimL[Bx], 2) + pow(uprimL[By], 2) + pow(uprimL[Bz], 2);
    double rhoL = uprimL[rho_p];
    double cfL = sqrt(0.5*(csL*csL + BL_squared/rhoL) + sqrt((pow(csL*csL + BL_squared/rhoL, 2) - 4*csL*csL*ByL*ByL/rhoL)));

    double csR = sqrt(gamma*uprimR[p]/uprimR[rho_p]);
    double ByR = uprimR[By];
    double BR_squared = pow(uprimR[Bx], 2) + pow(uprimR[By], 2) + pow(uprimR[Bz], 2);
    double rhoR = uprimR[rho_p];
    double cfR = sqrt(0.5*((csR*csR + BR_squared/rhoR) + sqrt(pow((csR*csR + BR_squared/rhoR), 2) - 4*csR*csR*ByR*ByR/rhoR)));
    //Use cfL and cfR to get sL and sR
    double sL = std::min(uprimL[vy],uprimR[vy]) - std::max(cfL, cfR);
    double sR = std::max(uprimL[vy], uprimR[vy]) + std::max(cfL, cfR);
    for (int i = 0; i < 8; i++){
        uHLL[i] = (sR*uR[i] - sL*uL[i] + getYFluxFunction(uL)[i] - getYFluxFunction(uR)[i])/(sR-sL);
    }
    return uHLL;
}

std::array<double, 8> getYfluxHLL(std::array<double,8> uL, std::array<double,8> uR){
    std::array<double,8> uHLL;
    std::array<double,8> fHLL;
    std::array<double,8> f_iPlusHalf;
    std::array<double,8> uprimL = ConsToPrim(uL);
    std::array<double,8> uprimR = ConsToPrim(uR);
    double gamma = getPars()[g];

    double csL = sqrt(gamma*uprimL[p]/uprimL[rho_p]);
    double ByL = uprimL[By];
    double BL_squared = pow(uprimL[Bx], 2) + pow(uprimL[By], 2) + pow(uprimL[Bz], 2);
    double rhoL = uprimL[rho_p];
    double cfL = sqrt(0.5*(csL*csL + BL_squared/rhoL) + sqrt((pow(csL*csL + BL_squared/rhoL, 2) - 4*csL*csL*ByL*ByL/rhoL)));

    double csR = sqrt(gamma*uprimR[p]/uprimR[rho_p]);
    double ByR = uprimR[By];
    double BR_squared = pow(uprimR[Bx], 2) + pow(uprimR[By], 2) + pow(uprimR[Bz], 2);
    double rhoR = uprimR[rho_p];
    double cfR = sqrt(0.5*((csR*csR + BR_squared/rhoR) + sqrt(pow((csR*csR + BR_squared/rhoR), 2) - 4*csR*csR*ByR*ByR/rhoR)));
    double sL = std::min(uprimL[vy],uprimR[vy]) - std::max(cfL, cfR);
    double sR = std::max(uprimL[vy], uprimR[vy]) + std::max(cfL, cfR);
    //std::cout << "Y HLL: " << "sL = " << sL << ", sR = " << sR << std::endl;
    for (int i = 0; i < uHLL.size(); i++){
        uHLL[i] = (sR*uR[i] - sL*uL[i] + getYFluxFunction(uL)[i] - getYFluxFunction(uR)[i])/(sR-sL);
        fHLL[i] = (sR*getYFluxFunction(uL)[i] - sL*getYFluxFunction(uR)[i] + sL*sR*(uR[i] - uL[i]))/(sR-sL);
        if (0 <= sL){
            f_iPlusHalf[i] = getYFluxFunction(uL)[i];
        }
        if (sL < 0 && 0 < sR){
            f_iPlusHalf[i] = fHLL[i];
        }
        if (0 >= sR){
            f_iPlusHalf[i] = getYFluxFunction(uR)[i];
        }
    }
    return f_iPlusHalf;
}

std::array<double, 8> getYfluxHLLC(std::array<double,8> uL, std::array<double,8> uR){
    std::array<double,8> uHLL = getYuHLL(uL, uR);
    std::array<double,8> f_iPlusHalf_HLLC;
    std::array<double,8> uStarL;
    std::array<double,8> uStarR;
    std::array<double,8> uprimL = ConsToPrim(uL);
    std::array<double,8> uprimR = ConsToPrim(uR);
    double gamma = getPars()[g];
    
    double csL = sqrt(gamma*uprimL[p]/uprimL[rho_p]);
    double ByL = uprimL[By];
    double BL_squared = pow(uprimL[Bx], 2) + pow(uprimL[By], 2) + pow(uprimL[Bz], 2);
    double rhoL = uprimL[rho_p];
    double cfL = sqrt(0.5*(csL*csL + BL_squared/rhoL) + sqrt((pow(csL*csL + BL_squared/rhoL, 2) - 4*csL*csL*ByL*ByL/rhoL)));

    double csR = sqrt(gamma*uprimR[p]/uprimR[rho_p]);
    double ByR = uprimR[By];
    double BR_squared = pow(uprimR[Bx], 2) + pow(uprimR[By], 2) + pow(uprimR[Bz], 2);
    double rhoR = uprimR[rho_p];
    double cfR = sqrt(0.5*((csR*csR + BR_squared/rhoR) + sqrt(pow((csR*csR + BR_squared/rhoR), 2) - 4*csR*csR*ByR*ByR/rhoR)));
    //double sL = std::min(std::abs(uprimL[vy]),std::abs(uprimR[vy])) - std::max(std::abs(cfL), std::abs(cfR));
    //double sR = std::max(std::abs(uprimL[vy]),std::abs(uprimR[vy])) + std::max(std::abs(cfL), std::abs(cfR));
    double sL = std::min(uprimL[vy],uprimR[vy]) - std::max(cfL, cfR);
    double sR = std::max(uprimL[vy],uprimR[vy]) + std::max(cfL, cfR);
    double vyR = uprimR[vy];
    double vyL = uprimL[vy];
    //double sL = std::min(vyL - cfL, vyR - cfR);
    //double sR = std::max(vyL + cfL, vyR + cfR);

    double ptL = uprimL[p] + 0.5*BL_squared;
    double ptR = uprimR[p] + 0.5*BR_squared;
    double sStar = (rhoR*vyR*(sR - vyR) - rhoL*vyL*(sL - vyL) + ptL - ptR - ByL*ByL + ByR*ByR)\
    /(rhoR*(sR - vyR) - rhoL*(sL-vyL));
    //std::cout << "Y HLL: " << "sL = " << sL << ", sR = " << sR << "s* = " << sStar <<std::endl;
    uStarL[rho] = rhoL*(sL - vyL)/(sL - sStar);
    uStarR[rho] = rhoR*(sR - vyR)/(sR - sStar);
    uStarL[rhovy] = uStarL[rho]*sStar;
    uStarR[rhovy] = uStarR[rho]*sStar;
    double BxL = uL[Bx];
    double BxR = uR[Bx];
    double BzL = uL[Bz];
    double BzR = uR[Bz];
    uStarL[rhovx] = uL[rhovx]*(sL - vyL)/(sL - sStar) - (uHLL[By]*uHLL[Bx] - ByL*BxL)/(sL - sStar);
    uStarR[rhovx] = uR[rhovx]*(sR - vyR)/(sR - sStar) - (uHLL[By]*uHLL[Bx] - ByR*BxR)/(sR - sStar);
    uStarL[rhovz] = uL[rhovz]*(sL - vyL)/(sL - sStar) - (uHLL[By]*uHLL[Bz] - ByL*BzL)/(sL - sStar);
    uStarR[rhovz] = uR[rhovz]*(sR - vyR)/(sR - sStar) - (uHLL[By]*uHLL[Bz] - ByR*BzR)/(sR - sStar);
    uStarR[Bx] = uHLL[Bx];
    uStarL[Bx] = uHLL[Bx];
    uStarR[By] = uHLL[By];
    uStarL[By] = uHLL[By];
    uStarR[Bz] = uHLL[Bz];
    uStarL[Bz] = uHLL[Bz];
    double vxHLL = uHLL[rhovx]/uHLL[rho];
    double vyHLL = uHLL[rhovy]/uHLL[rho];
    double vzHLL = uHLL[rhovz]/uHLL[rho];
    double Bv_HLL = uHLL[Bx]*vxHLL + uHLL[By]*vyHLL + uHLL[Bz]*vzHLL;
    double Bv_R = uprimR[Bx]*uprimR[vx] + uprimR[By]*uprimR[vy] + uprimR[Bz]*uprimR[vz];
    double Bv_L = uprimL[Bx]*uprimL[vx] + uprimL[By]*uprimL[vy] + uprimL[Bz]*uprimL[vz];
    double ptStarL = rhoL*(sL - vyL)*(sStar - vyL) + ptL - ByL*ByL + uHLL[By]*uHLL[By];
    double ptStarR = rhoR*(sR - vyR)*(sStar - vyR) + ptR - ByR*ByR + uHLL[By]*uHLL[By];
    uStarR[E] = uR[E]*(sR - vyR)/(sR - sStar) + (ptStarR*sStar - ptR*vyR - (uHLL[By]*(Bv_HLL)- ByR*Bv_R))/(sR - sStar);
    uStarL[E] = uL[E]*(sL - vyL)/(sL - sStar) + (ptStarL*sStar - ptL*vyL - (uHLL[By]*(Bv_HLL)- ByL*Bv_L))/(sL - sStar);
        for (int i = 0; i < 8; i++){
        if (0 <= sL){
            f_iPlusHalf_HLLC[i] = getYFluxFunction(uL)[i];
        }
        if (sL < 0 && 0 <= sStar){
            f_iPlusHalf_HLLC[i] = getYFluxFunction(uL)[i] + sL*(uStarL[i] - uL[i]);
        }
        if (sStar < 0 && 0 <= sR){
            f_iPlusHalf_HLLC[i] = getYFluxFunction(uR)[i] + sR*(uStarR[i] - uR[i]);
        }
        if (0 > sR){
            f_iPlusHalf_HLLC[i] = getYFluxFunction(uR)[i];
        }
    }
    for (int i = 0; i < 8; i++){
        if (0 <= sL){
            f_iPlusHalf_HLLC[i] = getYFluxFunction(uL)[i];
        }
        if (sL < 0 && 0 <= sStar){
            f_iPlusHalf_HLLC[i] = getYFluxFunction(uL)[i] + sL*(uStarL[i] - uL[i]);
        }
        if (sStar < 0 && 0 <= sR){
            f_iPlusHalf_HLLC[i] = getYFluxFunction(uR)[i] + sR*(uStarR[i] - uR[i]);
        }
        if (0 > sR){
            f_iPlusHalf_HLLC[i] = getYFluxFunction(uR)[i];
        }
    }

   return f_iPlusHalf_HLLC;
}

// Compute time step function
double computeTimeStep(const std::vector<std::vector<std::array<double, 8>>>& cons){
    // prim and cons here are the entire u vector which stores all variables at all i's
    int nxCells = getPars()[xCells];
    int nyCells = getPars()[yCells];
    std::vector<std::vector<std::array<double, 8>>> prim;
    prim.resize(nxCells+4, std::vector<std::array<double, 8>>(nyCells+4));
    if (cons.size() < nxCells + 4 || cons[0].size() < nyCells + 4) {
        throw std::runtime_error("Input vector 'cons' is too small.");
    }
    //std::cout << "cons size: " << cons.size() << ", " << cons[0].size() << std::endl;
    for (int i = 0; i < nyCells+4; i++){
        for (int j = 0; j < nxCells+4; j++){
            prim[i][j] = ConsToPrim(cons[i][j]);
        }
    }
    double gamma = getPars()[g];
    //double vx, vy, Bx, By, Bz, rho, B_squared, p, cs, termx, termy, term,  cfx, cfy, v, cf;
    std::vector<std::vector<double>> amax(nxCells + 4, std::vector<double>(nyCells + 4));
    for (int i = 0; i < nxCells+4; i++){
        for (int j = 0; j < nyCells+4; j++){
            double vx_value = prim[i][j][vx];
            double vy_value = prim[i][j][vy];
            double Bx_value = prim[i][j][Bx];
            double By_value = prim[i][j][By];
            double Bz_value = prim[i][j][Bz];
            double rho_value = prim[i][j][rho];
            double B_squared = Bx_value*Bx_value + By_value*By_value + Bz_value*Bz_value ;
            double p_value = prim[i][j][p];
            /*if (rho <= 0) {
                std::cout<<"Unphysical state detected: rho is non-positive in computeTimeStep"<<std::endl;;
            }
            if (p <= 0) {
                std::cout<<"Unphysical state detected: p is non-positive in computeTimeStep"<<std::endl;;
            }*/
            double cs = sqrt(gamma*p_value/rho_value);
            double term = (cs*cs + B_squared/rho_value)*(cs*cs + B_squared/rho_value) - 4*(cs*cs*B_squared/rho_value);
            double v = sqrt(vx_value*vx_value + vy_value*vy_value);
            double cf = sqrt(0.5*((cs*cs + (B_squared)/rho_value) + sqrt(term)));
            amax[i][j] = std::abs(v) + cf;
        }
    }
    double max_wave_speed = 0.0;
    for (int i = 0; i < amax.size(); i++){
        for (int j = 0; j < amax[i].size(); j++){
            if (amax[i][j] > max_wave_speed){
                max_wave_speed = amax[i][j];
            }
        }
    }

    double C = getPars()[CFL];
    double x0 = getPars()[xStart];
    double x1 = getPars()[xEnd];
    double dx = (x1-x0)/nxCells;
    double dy = (getPars()[yEnd] - getPars()[yStart])/nyCells;
    //if (max_wave_speed < 1e-12) throw std::runtime_error("Max wave speed is too small!");
    double dt = C*(std::min(dx,dy)/max_wave_speed);
    if (dt < 1e-20) {
        std::cout<<" !!!!!!! dt too small !!!!!!!!"<<std::endl;
    }
    //below for debugging only
    //std::cout<<"dt is "<<dt<<std::endl;
    return dt;
}

// Limiter related functions

double minbee(double r){
    const double eps = 1e-12;  // Small tolerance for floating-point comparisons

    // Check for NaN or Inf
    if (std::isnan(r) || std::isinf(r)) {
        return 0.0;  // Fallback to first-order (no limiter)
    }

    // Handle negative or zero r
    if (r <= 0) {
        return 0.0;
    }

    // Handle r in (0, 1]
    if (r <= 1) {
        return r;
    }

    // Handle r > 1
    double denominator = 1.0 + r;
    if (std::abs(denominator) < eps) {
        return 0.0;  // Avoid division by zero
    }
    return std::min(1.0, 2.0 / denominator);
}  
double vanLeer(double r) {
    if (r <= 0) return 0.0;
    return (2 * r) / (1 + r);
}
double superbee(double r) {
    return std::max(0.0, std::min(2 * r, 1.0));
}
double minmod(double r1, double r2) {
    if (r1 * r2 <= 0.0) {
        return 0.0;  // If the slopes have opposite signs, return zero
    }
    return (std::fabs(r1) < std::fabs(r2)) ? r1 : r2;  // Return the smaller magnitude of the two slopes
}

double calculateSlopeRatio(double numerator, double denominator) {
    const double eps = 1e-12;  // tolerance

    // Check for NaN or Inf in inputs
    if (std::isnan(numerator) || std::isinf(numerator) ||
        std::isnan(denominator) || std::isinf(denominator)) {
        return 0.0;  // Fallback to first-order (no limiter)
    }

    // Avoid division by zero
    if (std::abs(denominator) < eps) {
        return 0.0;
    }

    return numerator / denominator;
}

//Save to file function
void saveToDatFile(const std::vector<std::vector<std::array<double, 8>>>& u, const std::string& filename, int nxCells, int nyCells) {
   
    double x0 = getPars()[xStart]; 
    double x1 = getPars()[xEnd]; 
    double y0 = getPars()[yStart]; 
    double y1 = getPars()[yEnd];
    double dx = (x1 - x0) / nxCells;
    double dy = (y1 - y0) / nyCells;
    double gamma = getPars()[g];

    std::ofstream outFile(filename);

    outFile << "# i, j, Density, x-Velocity, y-Velocity, Pressure, Bx, By, Bz\n";

    for (int i = 1; i < nxCells + 2; i++) {
        double x = x0 + (i-1)*dx;
        for (int j = 1; j < nyCells + 2; j++) {
            double y = y0 + (j-1)*dy;
            double specific_internal_energy = u[i][j][p]/(u[i][j][rho_p]*(gamma - 1));
            outFile << x << " " << y << " ";  

            for (int var = 0; var < 8; var++) {
                outFile << u[i][j][var] << " ";  
            }
            outFile << specific_internal_energy << " ";
            outFile << "\n";  
        }
    }

    outFile.close();
}

void saveToDatFileDiagonal(const std::vector<std::vector<std::array<double, 8>>>& u, const std::string& filename, int nxCells, int nyCells) {
   
    double x0 = getPars()[xStart]; 
    double x1 = getPars()[xEnd]; 
    double y0 = getPars()[yStart]; 
    double y1 = getPars()[yEnd];
    double dx = (x1 - x0) / nxCells;
    double dy = (y1 - y0) / nyCells;
    double gamma = getPars()[g];
    //double theta = 3M_PI / 4;
    double tolerance = dx;

    std::ofstream outFile(filename);

    outFile << "# x, y, Diag. dist., Density, x-Velocity, y-Velocity, Pressure, Bx, By, Bz\n";

    double cosTheta = sqrt(2)/2.0;
    double sinTheta = sqrt(2)/2.0;
    
    for (int i = 2; i < nxCells+2; i++) {
        for (int j = 2; j < nyCells+2; j++) {
            double x = x0 + (i - 1.5) * dx;
            double y = y0 + (j - 1.5) * dy;
            //double dist_to_diag = sqrt(x*x + y*y); // Distance normal to diagonal

            if (x == y ) {  
                double diag_dist = sqrt(x*x + y*y);  
                double vt = u[i][j][vx]*sinTheta + u[i][j][vy]*cosTheta;
                double vp = u[i][j][vy]*sinTheta - u[i][j][vx]*cosTheta;
                double Bp = -u[i][j][Bx]*sinTheta + u[i][j][By]*cosTheta;
                double Bt = u[i][j][Bx]*sinTheta + u[i][j][By]*cosTheta;
                double specific_internal_energy = u[i][j][p]/(u[i][j][rho_p]*(gamma-1));
                outFile << x << " " << y << " " << diag_dist << " " << u[i][j][rho_p] << " " << vt << " " << vp << \
                " " << u[i][j][3] << " " << u[i][j][4] << " " << Bt << " " << Bp << \
                " " << u[i][j][7] << " " << specific_internal_energy << "\n";  
            }
        }
    }

    outFile.close();
}

//Initial data functions
std::vector<std::vector<std::array<double,8>>> setInitialDataCylindricalExplosion(){
    int nxCells = getPars()[xCells]; 
    int nyCells = getPars()[yCells];
    double x0 = getPars()[xStart]; 
    double x1 = getPars()[xEnd]; 
    double y0 = getPars()[yStart]; 
    double y1 = getPars()[yEnd]; 
    double dx = (x1 - x0) / nxCells;
    double dy = (y1 - y0) / nyCells;
    double Lx = x1-x0;
    double Ly = y1-y0;
    std::vector<std::vector<std::array<double, 8>>> u;
    u.resize(nxCells+4, std::vector<std::array<double, 8> >(nyCells + 4));
    for(int i = 0; i < nxCells+4; i++){
        double x = x0 + (i-0.5) * dx;
        for (int j = 0; j < nyCells+4; j++){
            std::array<double, 8> prim;
            double y = y0 + (j-0.5) * dy;
            // Compute distance from the center of the circle (1,1)
            double r = sqrt((x - 1.0) * (x - 1.0) + (y - 1.0) * (y - 1.0));
            double R = 0.4 ;
            if(r <= R){
                prim = {1, 0, 0, 0, 1, 0, 0, 0};
            }
            else {
                prim = {0.125, 0, 0, 0, 0.1, 0, 0, 0};
            }
            u[i][j] = PrimToCons(prim);
        }
    }
    return u;
}

std::vector<std::vector<std::array<double,8>>> setInitialDataTest1_X(){
    int nxCells = getPars()[xCells]; 
    int nyCells = getPars()[yCells];
    double x0 = getPars()[xStart]; 
    double x1 = getPars()[xEnd]; 
    double y0 = getPars()[yStart]; 
    double y1 = getPars()[yEnd]; 
    double dx = (x1 - x0) / nxCells;
    double dy = (y1 - y0) / nyCells;
    double Lx = x1-x0;
    double Ly = y1-y0;
    std::vector<std::vector<std::array<double, 8>>> u;
    u.resize(nxCells+4, std::vector<std::array<double, 8> >(nyCells + 4));
    for(int i = 0; i < nxCells+4; i++){
        double x = x0 + (i-1.5) * dx;
        for (int j = 0; j < nyCells+4; j++){
            std::array<double, 8> prim;
            double y = y0 + (j-1.5) * dy;
            if(x <= 0.5*(x1-x0)){
                prim = {1, 0, 0, 0, 1, 0.75, 1, 0};
            }
            else if (x > 0.5*(x1-x0)){
                prim = {0.125, 0, 0, 0, 0.1, 0.75, -1, 0};
            }
            u[i][j] = PrimToCons(prim);
        }
    } 
    return u;
}

std::vector<std::vector<std::array<double, 8>>> setInitialDataDiagonal() {
    int nxCells = getPars()[xCells]; 
    int nyCells = getPars()[yCells];
    double x0 = getPars()[xStart]; 
    double x1 = getPars()[xEnd]; 
    double y0 = getPars()[yStart]; 
    double y1 = getPars()[yEnd]; 
    double dx = (x1 - x0) / nxCells;
    double dy = (y1 - y0) / nyCells;
    double theta = M_PI / 4;  // 45-degree diagonal line

    std::vector<std::vector<std::array<double, 8>>> u;
    u.resize(nxCells + 4, std::vector<std::array<double, 8>>(nyCells + 4));

    double cosTheta = std::cos(theta);
    double sinTheta = std::sin(theta);
    double sqrt2 = sqrt(2.0);
    for (int i = 0; i < nxCells + 4; i++) {
        double x = x0 + (i - 0.5) * dx;
        for (int j = 0; j < nyCells + 4; j++) {
            std::array<double, 8> prim;
            double y = y0 + (j - 0.5) * dy;

            // Define diagonal interface: y = tan(theta) * x
            if  (x + y < 800){   
                prim = {1.0, 0.0, 0.0, 0.0, 1.0, -0.177 , 1.237, 0.0};
            } else {   
                prim = {0.125, 0.0, 0.0, 0.0, 0.1, 1.237 , -0.177, 0.0};
            }
            u[i][j] = PrimToCons(prim);  // Convert to conservative variables
        }
    }
    return u;
}

std::vector<std::vector<std::array<double,8>>> setInitialDataTest1_Y() {
    int nxCells = getPars()[xCells]; 
    int nyCells = getPars()[yCells];
    double x0 = getPars()[xStart]; 
    double x1 = getPars()[xEnd]; 
    double y0 = getPars()[yStart]; 
    double y1 = getPars()[yEnd]; 
    double dx = (x1 - x0) / nxCells;
    double dy = (y1 - y0) / nyCells;
    double Lx = x1 - x0;
    double Ly = y1 - y0;

    std::vector<std::vector<std::array<double, 8>>> u;
    u.resize(nxCells + 4, std::vector<std::array<double, 8>>(nyCells + 4));

    for(int i = 0; i < nxCells + 4; i++) {
        double x = x0 + (i - 0.5) * dx;  // Compute x-coordinate (not used in condition)
        for (int j = 0; j < nyCells + 4; j++) {
            std::array<double, 8> prim;
            double y = y0 + (j - 0.5) * dy;  // Compute y-coordinate
            
            if (y <= 0.5 * (y1 - y0)) {  // Left state (bottom half)
                prim = {1, 0, 0, 0, 1, 1, 0.75, 0};
            } 
            else {  // Right state (top half)
                prim = {0.125, 0, 0, 0, 0.1, -1, 0.75, 0};
            }
            u[i][j] = PrimToCons(prim);  // Convert primitive variables to conservative form
        }
    } 
    return u;
}

std::vector<std::vector<std::array<double, 8>>> setInitialDataOrszagTang() {
    int nxCells = getPars()[xCells]; 
    int nyCells = getPars()[yCells];
    double x0 = getPars()[xStart]; 
    double x1 = getPars()[xEnd]; 
    double y0 = getPars()[yStart]; 
    double y1 = getPars()[yEnd]; 
    double dx = (x1 - x0) / nxCells;
    double dy = (y1 - y0) / nyCells;
    double Lx = x1 - x0;
    double Ly = y1 - y0;
    double gamma = getPars()[g];
    
    std::vector<std::vector<std::array<double, 8>>> u(nxCells + 4, std::vector<std::array<double, 8>>(nyCells + 4));

    for (int i = 2; i < nxCells + 4; i++) {
        double x = x0 + (i - 1.5) * dx;  
        for (int j = 2; j < nyCells + 4; j++) {
            double y = y0 + (j -1.5) * dy;  
            
            std::array<double, 8> prim;
            prim[rho] = gamma * gamma;               
            prim[vx] = -sin(2.0 * M_PI * y);        
            prim[vy] = sin(2.0 * M_PI * x);         
            prim[vz] = 0.0;   
            prim[p] = gamma;                         
            prim[Bx] = -sin(2.0 * M_PI * y);        
            prim[By] = sin(4.0 * M_PI * x);         
            prim[Bz] = 0.0;                         

            u[i][j] = PrimToCons(prim);
        }
    }
    return u;
}

std::vector<std::vector<std::array<double, 8>>> setInitialDataMHDRotorPaper() {
    int nxCells = getPars()[xCells]; 
    int nyCells = getPars()[yCells];
    double x0 = getPars()[xStart]; 
    double x1 = getPars()[xEnd]; 
    double y0 = getPars()[yStart]; 
    double y1 = getPars()[yEnd]; 
    double dx = (x1 - x0) / nxCells;
    double dy = (y1 - y0) / nyCells;

    std::vector<std::vector<std::array<double, 8>>> u(nxCells + 4, std::vector<std::array<double, 8>>(nyCells + 4));


    for (int i = 2; i < nxCells + 4; i++) {
        double x = x0 + (i - 1.5) * dx;  
        for (int j = 2; j < nyCells + 4; j++) {
            double y = y0 + (j - 1.5) * dy;  
            double r = sqrt((x-0.5) *(x-0.5) + (y-0.5)*(y-0.5));  
            double r0 = 0.1;
            double r1 = 0.115;
            double f = (r1-r)*(r1-r0);
            double v0 = 2;
            std::array<double, 8> prim;
            prim[p] = 0.5;   // Pressure
            prim[Bx] = 2.5 / sqrt(4.0 * M_PI);  // Magnetic field Bx
            prim[By] = 0.0; // Magnetic field By
            prim[Bz] = 0.0; // Magnetic field Bz
            prim[vz] = 0.0;  // Velocity vz (always 0)

            if (r < 0.1) {
                prim[rho_p] = 10.0;  // Density
                prim[vx] = - v0 * (y - 0.5) / r0 ;    // vx
                prim[vy] = v0 * (x - 0.5) / r0;    // vy
            } else if (r >= 0.1 && r < 0.115) {
                prim[rho_p] = 1.0 + 9.0 * f;   // Density
                prim[vx] = -f * v0 * (y - 0.5) / r ;          // vx
                prim[vy] = f * v0 * (x - 0.5) / r ;    // vy
            } else {
                prim[rho_p] = 1.0;  // Density
                prim[vx] = 0.0;  // vx
                prim[vy] = 0.0;  // vy
            }

            u[i][j] = PrimToCons(prim);
        }
    }
    return u;
}

std::vector<std::vector<std::array<double, 8>>> setInitialDataMHDRotorAssignment() {
    int nxCells = getPars()[xCells]; 
    int nyCells = getPars()[yCells];
    double x0 = getPars()[xStart]; 
    double x1 = getPars()[xEnd]; 
    double y0 = getPars()[yStart]; 
    double y1 = getPars()[yEnd]; 
    double dx = (x1 - x0) / nxCells;
    double dy = (y1 - y0) / nyCells;
    double gamma = getPars()[g];

    std::vector<std::vector<std::array<double, 8>>> u(nxCells + 4, std::vector<std::array<double, 8>>(nyCells + 4));

    for (int i = 2; i < nxCells + 4; i++) {
        double x = x0 + (i - 1.5) * dx;  // Cell-centered x-coordinate
        for (int j = 2; j < nyCells + 4; j++) {
            double y = y0 + (j - 1.5) * dy;  // Cell-centered y-coordinate
            double r = sqrt((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5));  // Distance from center
            double ur = 5.0 - 10.0*y;            
            double vr = 10.0*x - 5.0;
            double r0 = 0.1;
            double r1 = 0.115;
            double fr = (23-200*r)/27;
            double f = (r1-r)*(r1-r0);
            double v0 = 1;

            std::array<double, 8> prim;

            if (r >= 0 && r < r0) {
                prim[rho] = 10.0;
                prim[vx] = ur;
                prim[vy] =vr;
            } else if (r >= r0 && r < r1) {
                prim[rho] = 1.0 + 9.0 * fr;
                prim[vx] = ur*fr;
                prim[vy] = vr*fr;
            } else if (r > r1) {
                prim[rho] = 1.0;
                prim[vx] = 0.0;
                prim[vy] = 0.0;
            }

            prim[vz] = 0.0;
            prim[p] = 1;
            prim[Bx] = 2.5 / sqrt(4 * M_PI);
            prim[By] = 0;
            prim[Bz] = 0;
            //prim[psi] = 0.0;  // If needed for divergence cleaning

            // Convert primitive variables to conservative form if needed
            u[i][j] = PrimToCons(prim);
        }
    }
    return u;
}
//Boundary conditions functions
std::vector<std::vector<std::array<double, 8>>> applyTransmissiveBoundaries(std::vector<std::vector<std::array<double, 8>>> u, int nxCells, int nyCells) {

    // Apply transmissive boundary conditions in the x-direction
    for (int j = 0; j < nyCells + 4; j++) {
        u[1][j] = u[2][j];      
        u[0][j] = u[1][j];      

        u[nxCells+2][j] = u[nxCells+1][j]; 
        u[nxCells+3][j] = u[nxCells+2][j];  
    }

    // Apply transmissive boundary conditions in the y-direction
    for (int i = 0; i < nxCells + 4; i++) {
        u[i][1] = u[i][2];      
        u[i][0] = u[i][1];      

        u[i][nyCells+2] = u[i][nyCells+1];  
        u[i][nyCells+3] = u[i][nyCells+2];  
    } 

    return u;  // Return modified u
}

std::vector<std::vector<std::array<double, 8>>> applyPeriodicBoundaries(std::vector<std::vector<std::array<double, 8>>> u, int nxCells, int nyCells) {

    // Apply periodic boundary conditions in the x-direction
    for (int j = 0; j < nyCells + 4; j++) {
        u[0][j] = u[nxCells][j];      
        u[1][j] = u[nxCells+1][j];     

        u[nxCells+2][j] = u[2][j];  
        u[nxCells+3][j] = u[3][j];  
    }

    // Apply periodic boundary conditions in the y-direction
    for (int i = 0; i < nxCells + 4; i++) {
        u[i][0] = u[i][nyCells];      
        u[i][1] = u[i][nyCells+1];    

        u[i][nyCells+2] = u[i][2];  
        u[i][nyCells+3] = u[i][3];  
    } 

    return u;  // Return modified u
}

/*std::vector<std::vector<std::array<double, 8>>> applyTransmissiveBCDiagonal(std::vector<std::vector<std::array<double, 8>>>& u, int nxCells, int nyCells) {
    double cosTheta = std::cos(M_PI / 4.0);
    double sinTheta = std::sin(M_PI / 4.0);

    //Top left corner
    for (int j = 0; j < nyCells+4; j++) {
            u[0][j] = u[2][j];  // Copy from the first interior cells
            u[1][j] = u[2][j];
    }

    // Bottom-right corner (outflow side)
    for (int j = 0; j < nyCells+4 ; j++) {
            u[nyCells+2][j] = u[nyCells + 1][j];  // Extrapolate from last interior cells
            u[nyCells+3][j] = u
    }

    // Left and Right Boundaries (aligned with diagonal direction)
    for (int i = 0; i < nxCells+4; i++) {
        for (int ghost = 0; ghost < 2; ghost++) {
            u[i][ghost] = u[i][2];
            u[i][nyCells + ghost] = u[i][nyCells + 1];
        }
    }
    return u;

}*/


double computeMaxDivB(const std::vector<std::vector<std::array<double, 8>>>& u, double dx, double dy, int nxCells, int nyCells){
    std::vector<std::vector<double>> divB(nxCells+4, std::vector<double>(nyCells+4, 0.0));
    double max_divB = 0;
    for (int i = 2; i < nxCells+2; i++) {
        for (int j = 2; j < nyCells+2; j++) {
            double dBx_dx = (u[i+1][j][Bx] - u[i-1][j][Bx]) / (2.0 * dx);
            double dBy_dy = (u[i][j+1][By] - u[i][j-1][By]) / (2.0 * dy);
            double divB = dBx_dx + dBy_dy;
            max_divB = std::max(max_divB, std::abs(divB));
        }   
    }
    /*double max_divB = 0;
    for (int i = 2; i < nxCells+2; i++) {
        for (int j = 2; j < nyCells+2; j++) {
            if (divB[i][j] > max_divB){
                max_divB = divB[i][j];
            }
        }
    } */
    //divB = applyTransmissiveBoundaries(divB, nxCells, nyCells);
    return max_divB;
}

double computeMaxDivB2(const std::vector<std::vector<std::array<double, 8>>>& u,
    double dx, double dy, int nx, int ny) {
    double maxDivB = 0.0;
    double divB_total = 0.0;
    for (int i = 2; i < nx + 2; ++i) {
            for (int j = 2; j < ny + 2; ++j) {
                double dBx_dx = (u[i+1][j][Bx] - u[i-1][j][Bx]) / (2.0 * dx);
                double dBy_dy = (u[i][j+1][By] - u[i][j-1][By]) / (2.0 * dy);
                double divB = std::abs(dBx_dx + dBy_dy);
                maxDivB = std::max(maxDivB, divB);
                divB_total += divB;
            }
    }
    return maxDivB;
}

//Time update function
int SolveMHD2D(){
    int nxCells = getPars()[xCells]; 
    int nyCells = getPars()[yCells];
    double x0 = getPars()[xStart]; 
    double x1 = getPars()[xEnd]; 
    double y0 = getPars()[yStart]; 
    double y1 = getPars()[yEnd]; 
    double t0 = getPars()[tStart]; 
    double tStop = getPars()[tEnd];
    double dx = (x1 - x0) / nxCells;
    double dy = (y1 - y0) / nyCells;
    double dt;

    std::vector<std::vector<std::array<double, 8>>> u;
    u.resize(nxCells+4, std::vector<std::array<double, 8> >(nyCells + 4));
    
    std::vector<std::vector<std::array<double, 8>>> xflux;
    xflux.resize(nxCells+4, std::vector<std::array<double, 8> >(nyCells + 4));

    std::vector<std::vector<std::array<double, 8>>> yflux;
    yflux.resize(nxCells+4, std::vector<std::array<double, 8> >(nyCells + 4));
 
    std::vector<std::vector<std::array<double, 8>>> uPlus1;
    uPlus1.resize(nxCells+4, std::vector<std::array<double, 8> >(nyCells + 4));

    u = setInitialDataOrszagTang();

    int ConstantData = 0; //set to 1 for constant data in all time steps
    double t = t0;
    int dt_counter = 0;
    int disableLimiter = 0; //set to 1 to set all limiters to 0 


    //Declare vector to store the limited reconstrcuted states
    std::vector<std::vector<std::array<double, 8>>> uBarL;
    uBarL.resize(nxCells+4, std::vector<std::array<double, 8> >(nyCells + 4));
    std::vector<std::vector<std::array<double, 8>>> uBarR;
    uBarR.resize(nxCells+4, std::vector<std::array<double, 8> >(nyCells + 4));
    std::vector<std::vector<std::array<double, 8>>> uBarT;
    uBarT.resize(nxCells+4, std::vector<std::array<double, 8> >(nyCells + 4));
    std::vector<std::vector<std::array<double, 8>>> uBarB;
    uBarB.resize(nxCells+4, std::vector<std::array<double, 8> >(nyCells + 4));

    //Declare vectors for half-timestep updated reconstructed states
    std::vector<std::vector<std::array<double, 8>>> uBarLUpdate;
    uBarLUpdate.resize(nxCells+4, std::vector<std::array<double, 8> >(nyCells + 4));
    std::vector<std::vector<std::array<double, 8>>> uBarRUpdate;
    uBarRUpdate.resize(nxCells+4, std::vector<std::array<double, 8> >(nyCells + 4));
    std::vector<std::vector<std::array<double, 8>>> uBarTUpdate;
    uBarTUpdate.resize(nxCells+4, std::vector<std::array<double, 8> >(nyCells + 4));
    std::vector<std::vector<std::array<double, 8>>> uBarBUpdate;
    uBarBUpdate.resize(nxCells+4, std::vector<std::array<double, 8> >(nyCells + 4));

    //Declare vecotr that stores the flux function f(u) for each limited variable to use at half time step update
    std::vector<std::array<double,8>> xfluxFunction;
    xfluxFunction.resize(nxCells+4);
    std::vector<std::array<double,8>> yfluxFunction;
    yfluxFunction.resize(nyCells+4);
    std::ofstream outFile("Max_DivB_Orszag_noCleaning.dat");
    outFile << "# t, max_DivB\n";
    do{
        dt = computeTimeStep(u);
        t = t + dt;
        dt_counter++;
        std::cout<< "t= " << t<< " Timestep: " << dt_counter << " dt: "<<dt<<std::endl;            
        u = applyPeriodicBoundaries(u, nxCells, nyCells);
        //Debugging step of initial data////////////////////////////////////
            for (int i = 0; i < nxCells+4; i++) {
                for (int j = 0; j < nyCells+4; j++) {
                    for (int var = 0; var < 8; var++) {
                        if (std::isnan(u[i][j][var])) {
                            std::cout << "NaN detected in u BEFORE flux calculation at i=" << i << ", j=" << j << ", var=" << var << " timestep: " << dt_counter<< std::endl;
                            //return -1;
                        }
                    }
                }
            }

         //Apply limiter to get uBarL for cells i to nCells
        for (int i = 1; i < nxCells+3; i++){
            for (int j = 1; j < nxCells+3; j++){
                for (int var = 0; var < 8; var++){
                double rx_num = u[i][j][var] - u[i-1][j][var];
                double rx_den = u[i+1][j][var] - u[i][j][var];
                double rx = calculateSlopeRatio(rx_num,rx_den);
                //double dix = 0.5*(u[i+1][j][var]-u[i-1][j][var]);
                double dix = minmod(u[i+1][j][var] - u[i][j][var], u[i][j][var] - u[i-1][j][var]);
                if (std::isnan(u[i+1][j][var]) || std::isnan(u[i-1][j][var])) {
                    std::cout << "NaN detected in u array at i=" << i << ", j=" << j << " at timestep: " << dt_counter << std::endl;
                    std::cout << "u[i+1][j][var]: " << u[i+1][j][var] << ", u[i-1][j][var]: " << u[i-1][j][var] << std::endl;
                    //return -1;
                }
                double limx;
                if (disableLimiter == 0){
                    limx = minbee(rx);
                }
                if (disableLimiter == 1){
                    limx = 0;
                }
                uBarL[i][j][var] = u[i][j][var] - 0.5*limx*dix; 
                uBarR[i][j][var] = u[i][j][var] + 0.5*limx*dix;   
                }       
            }
        }
        //Apply half time step update to get half timestep updates
        for (int i = 1; i < nxCells+3; i++){
            for (int j = 1;  j < nyCells+3; j++){
                for (int var = 0; var < 8; var++){
                    uBarLUpdate[i][j][var] = uBarL[i][j][var] - \
                    0.5*(dt/dx)*(getXFluxFunction(uBarR[i][j])[var] - getXFluxFunction(uBarL[i][j])[var]);
                    uBarRUpdate[i][j][var] = uBarR[i][j][var] -\
                    0.5*(dt/dx)*(getXFluxFunction(uBarR[i][j])[var] - getXFluxFunction(uBarL[i][j])[var]);
                }
            }
        }
        //Fill flux arrays
        for (int i = 1; i < nxCells+2; i++){
            for (int j = 1; j < nyCells+2; j++){
                xflux[i][j] = getXfluxFORCE(uBarRUpdate[i][j],uBarLUpdate[i+1][j],dt);
            }
        }
        if (ConstantData == 1){ //this is for tests 1 and 2 only
            for (int i = 2; i < nxCells+1; i++){
                for (int j = 2; j < nyCells+1; j++){
                    for (int var = 0; var < 8; var++){
                        uPlus1[i][j][var] = u[i][j][var];
                    }
                }
            }
        }
        if (ConstantData == 0){ //real update
            for (int i = 2; i < nxCells+2; i++){
                for (int j = 2; j < nyCells+2; j++){
                    for (int var = 0; var < 8; var++){
                        uPlus1[i][j][var] = u[i][j][var] - (dt/dx)*(xflux[i][j][var] - xflux[i-1][j][var]);
                    }
                }
            }
        }

        uPlus1 = applyPeriodicBoundaries(uPlus1, nxCells, nyCells);

        for (int i = 1; i < nxCells+3; i++){
            for (int j = 1; j < nyCells+3; j++){
                for (int var = 0; var < 8; var++){
                double ry_num = uPlus1[i][j][var] - uPlus1[i][j-1][var];
                double ry_den = uPlus1[i][j+1][var] - uPlus1[i][j][var];
                double ry = calculateSlopeRatio(ry_num,ry_den);
                //double diy = 0.5*(uPlus1[i][j+1][var]-uPlus1[i][j-1][var]);
                double diy = minmod(uPlus1[i][j+1][var] - uPlus1[i][j][var], uPlus1[i][j][var] - uPlus1[i][j-1][var]);
                if (std::isnan(uPlus1[i][j+1][var]) || std::isnan(uPlus1[i][j-1][var])) {
                    std::cout << "NaN detected in u array at i=" << i << ", j=" << j << " at timestep: " << dt_counter << std::endl;
                    std::cout << "u[i][j+1][var]: " << uPlus1[i][j+1][var] << ", u[i][j-1][var]: " << uPlus1[i][j-1][var] << std::endl;
                    //return -1;
                }
                double limy;
                if (disableLimiter == 0){
                    limy = minbee(ry);

                }
                if (disableLimiter == 1){
                    limy = 0;
                }
                uBarB[i][j][var] = uPlus1[i][j][var] - 0.5*limy*diy; 
                uBarT[i][j][var] = uPlus1[i][j][var] + 0.5*limy*diy;
                }
            }
        }

        //Apply half time step update to get half timestep updates
        for (int i = 1; i < nxCells+3; i++){
            for (int j = 1;  j < nyCells+3; j++){
                for (int var = 0; var < 8; var++){
                    uBarTUpdate[i][j][var] = uBarT[i][j][var] - \
                    0.5*(dt/dy)*(getYFluxFunction(uBarT[i][j])[var] - getYFluxFunction(uBarB[i][j])[var]); 
                    uBarBUpdate[i][j][var] = uBarB[i][j][var] -\
                    0.5*(dt/dy)*(getYFluxFunction(uBarT[i][j])[var] - getYFluxFunction(uBarB[i][j])[var]);
                }
            }
        }

        for (int i = 1; i < nxCells+2; i++){
            for (int j = 1; j < nyCells+2; j++){
                yflux[i][j] = getYfluxFORCE(uBarTUpdate[i][j],uBarBUpdate[i][j+1],dt);
            }
        }

        for (int i = 2; i < nxCells+2; i++){
            for (int j = 2; j < nyCells+2; j++){
                for (int var = 0; var < 8; var++){
                    uPlus1[i][j][var] = uPlus1[i][j][var] - (dt/dy)*(yflux[i][j][var] - yflux[i][j-1][var]) ;
                }
            }
        }

        uPlus1 = applyPeriodicBoundaries(uPlus1, nxCells, nyCells);

        u = uPlus1;
        double maxDivB = computeMaxDivB2(u, dx, dy, nxCells, nyCells);
        outFile << t << " " << maxDivB << "\n";

    }while(t < tStop);
    outFile.close();

    std::cout << "number of timesteps: " << dt_counter << std::endl;

    //Convert back to primitive variable before outputting to file
    for (int i = 0; i < nxCells+4; i++){
        for (int j = 0; j < nyCells+4; j++){
            u[i][j] = ConsToPrim(u[i][j]);
        }
    }
    
    saveToDatFile(u, "MHD2D_Orszag_256.dat", nxCells, nyCells);
    return 0;
}

int main(){
    auto start = std::chrono::high_resolution_clock::now();
    SolveMHD2D();
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << "Execution time: " << duration.count() << " ms" << std::endl;
    return 0;
}
#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <algorithm>
#include <array>
#include <cstdlib>
#include <chrono>
#include <iomanip>

enum ConservedVariables {rho, rhoV, E, rhoL};
enum PrimitiveVariables {rho_p, v, p, l};
enum Parameters {xStart, xEnd, tStart, tEnd, xCells, CFL, g, heat, eA, dt_s};

std::array<double, 10> getPars(){
    std::array<double, 10> parameters;
    parameters[xStart] = 0;
    parameters[xEnd] = 1;
    parameters[tStart] = 0;
    parameters[tEnd] = 1.4;
    parameters[xCells] = 10000;
    parameters[CFL] = 0.8;
    parameters[g] = 1.4;
    parameters[heat] = 2.0;
    parameters[eA] = 14.99;
    parameters[dt_s] = 1e-5;
    return parameters;
}

std::array<double, 4> PrimToCons(std::array<double, 4> prim){
    std::array<double, 4> cons;
    double gamma = getPars()[g];
    cons[rho] = prim[rho_p];  // ρ = ρ
    cons[rhoV] = prim[rho_p]*prim[v];  // ρv = ρ*v
    cons[E] = prim[p]/(gamma*(gamma-1)) + 0.5*prim[rho_p]*prim[v]*prim[v]; // E = p/(γ-1) + 1/2*ρv^2
    cons[rhoL] = prim[rho_p]*prim[l];
    return cons;
}

std::array<double, 4> ConsToPrim(std::array<double, 4> cons){
    std::array<double, 4> prim;
    double gamma = getPars()[g];
    prim[rho_p] = cons[rho]; // ρ = ρ
    prim[v] = cons[rhoV]/cons[rho]; // v = ρ*v/ρ
    prim[p] = (cons[E] - 0.5*cons[rhoV]*cons[rhoV]/cons[rho])*(gamma*(gamma-1)); // p = [E - 1/2*(ρv)^2/ρ]*(γ-1)
    prim[l] = cons[rhoL]/cons[rho];
    return prim;
}

std::array<double, 4> getFluxFunction(std::array<double,4> u_i){
    std::array<double, 4> fluxFunction;
    std::array<double, 4> prim = ConsToPrim(u_i);
    double gamma = getPars()[g];
    fluxFunction[rho] = u_i[rhoV]; // rho flux: rho*v
    fluxFunction[rhoV] = u_i[rhoV]*prim[v] + prim[p]/gamma; // rho*v flux: rho*v*v + p
    fluxFunction[E] = (u_i[E] + prim[p]/gamma)*prim[v]; // E flux: (E+p)v
    fluxFunction[rhoL] = u_i[rho]*prim[l]*prim[v];
    return fluxFunction;
}

double computeHyperbolicTimeStep(const std::vector<std::array<double, 4>>& u_cons){
    int nCells = getPars()[xCells];
    std::vector<std::array<double, 4>> u_prim;
    u_prim.resize(nCells+3);
    for (int i = 0; i < nCells+3; i++){
        u_prim[i] = ConsToPrim(u_cons[i]);
    }
    double gamma = getPars()[g];
    std::vector<double> amax;
    amax.resize(nCells+3);
    for (int i = 0; i < nCells+3; i++){
         amax[i] = std::abs(u_prim[i][v]) + sqrt((gamma*u_prim[i][p])/u_prim[i][rho_p]); // amax = abs(v) + cs
    }
    double max_wave_speed = *std::max_element(amax.begin(), amax.end());
    double C = getPars()[CFL];
    double x0 = getPars()[xStart];
    double x1 = getPars()[xEnd];
    double dx = (x1-x0)/nCells;
    double dt = C*(dx/max_wave_speed);

    return dt;
}

std::array<double, 4> getFluxLF(std::array<double, 4> u_i, std::array<double, 4>  u_iPlus1, double dt){
    //NOTE: u_iPlus1 here is u(i+1,n), NOT u(i,n+1)
    int nCells = getPars()[xCells];
    double x0 = getPars()[xStart];
    double x1 = getPars()[xEnd];
    double dx = (x1-x0)/nCells;
    std::array<double, 4> f_iPlusHalf;
    std::array<double,4> flux_ui = getFluxFunction(u_i);
    std::array<double,4> flux_uiPlus1 = getFluxFunction(u_iPlus1);
    for (int var = 0; var < 4; var++){
        f_iPlusHalf[var] = 0.5*(dx/dt)*(u_i[var] - u_iPlus1[var]) + 0.5*(flux_ui[var] + flux_uiPlus1[var]);
    }

    return f_iPlusHalf;
}

std::array<double, 4> getFluxRI(std::array<double,4> u_i, std::array<double,4> u_iPlus1, double dt){
    int nCells = getPars()[xCells];
    double x0 = getPars()[xStart];
    double x1 = getPars()[xEnd];
    double dx = (x1-x0)/nCells;
    std::array<double, 4> u_iPlusHalf;
    std::array<double, 4> f_iPlusHalf;
    std::array<double, 4> flux_ui = getFluxFunction(u_i);
    std::array<double, 4> flux_uiPlus1 = getFluxFunction(u_iPlus1);   

    for (int var = 0; var < 4; var++){
        u_iPlusHalf[var] = 0.5*(u_i[var] + u_iPlus1[var]) - 0.5*(dt/dx)*(flux_uiPlus1[var] - flux_ui[var]);
    }

    std::array<double, 4> flux_u_iPlusHalf = getFluxFunction(u_iPlusHalf);

    for (int var = 0; var < 4; var++){
        f_iPlusHalf[var] = flux_u_iPlusHalf[var];
    }

    return f_iPlusHalf;
}

std::array<double, 4> getfluxFORCE(std::array<double,4> u_i, std::array<double,4> u_iPlus1, double dt){
    int nCells = getPars()[xCells];
    double x0 = getPars()[xStart];
    double x1 = getPars()[xEnd];
    double dx = (x1-x0)/nCells;
    std::array<double, 4> f_iPlusHalf_FORCE;
    for (int var = 0; var < 4; var++){
        f_iPlusHalf_FORCE[var] = 0.5*getFluxLF(u_i, u_iPlus1, dt)[var] + 0.5*getFluxRI(u_i, u_iPlus1, dt)[var];
    }
    return f_iPlusHalf_FORCE;
}

std::array<double, 4> getfluxHLL(std::array<double,4> uL, std::array<double,4> uR){
    std::array<double,4> uHLL;
    std::array<double,4> fHLL;
    std::array<double,4> f_iPlusHalf;
    std::array<double,4> primL = ConsToPrim(uL);
    std::array<double,4> primR = ConsToPrim(uR);
    double gamma = getPars()[g];
    double rhoL = uL[0];
    double rhoR = uR[0];
    double rhoVxL = uL[1];
    double rhoVxR = uR[1];
    double rhoL_L = uL[3];
    double rhoL_R = uR[3];
    double EL = uL[2];
    double ER = uR[2];

    double vxL = primL[1];
    double vxR = primR[1];
    double lamda_L = primL[3];
    double lamda_R = primR[3];
    double pL = primL[2];
    double pR = primR[2];

    // sound speeds
    double csL = sqrt(gamma*pL/rhoL);
    double csR = sqrt(gamma*pR/rhoR);

    // sL, sR, and s* 
    double sCross = std::max(std::abs(vxL) +  csL, std::abs(vxR) + csR);
    double sL = -sCross;
    double sR = sCross;

    for (int i = 0; i < uHLL.size(); i++){
        uHLL[i] = (sR*uR[i] - sL*uL[i] + getFluxFunction(uL)[i] - getFluxFunction(uR)[i])/(sR-sL);
        fHLL[i] = (sR*getFluxFunction(uL)[i] - sL*getFluxFunction(uR)[i] + sL*sR*(uR[i] - uL[i]))/(sR-sL);
        if (0 <= sL){
            f_iPlusHalf[i] = getFluxFunction(uL)[i];
        }
        else if (sL < 0 && 0 < sR){
            f_iPlusHalf[i] = fHLL[i];
        }
        else if (0 >= sR){
            f_iPlusHalf[i] = getFluxFunction(uR)[i];
        }
    }
    return f_iPlusHalf;
}

std::array<double, 4> getfluxHLLC(std::array<double, 4> uL, std::array<double, 4> uR){
    std::array<double, 4> f_iPlusHalf_HLLC;
    //Create vectors to hold intermediate state varibal values
    std::array<double, 4> uStarL;
    std::array<double, 4> uStarR;
    double gamma = getPars()[g];
    double rhoL = uL[0];
    double rhoR = uR[0];
    double rhoVL = uL[1];
    double rhoVR = uR[1];
    double EL = uL[2];
    double ER = uR[2];
    double pL = ConsToPrim(uL)[2];
    double pR = ConsToPrim(uR)[2];
    double vL = ConsToPrim(uL)[1];
    double vR = ConsToPrim(uR)[1];
    double rhoL_L = uL[3];
    double rhoL_R = uR[3];
    double lamdaL = ConsToPrim(uL)[3];
    double lamdaR = ConsToPrim(uR)[3];
    
    // sound speeds
    double csL = sqrt(gamma*pL/rhoL);
    double csR = sqrt(gamma*pR/rhoR);

    // sL, sR, and s* 
    double sCross = std::max(std::abs(vL) +  csL, std::abs(vR) + csR);
    double sL = -sCross;
    double sR = sCross;
    double sStar = (pR - pL + rhoVL*(sL - vL) - rhoVR*(sR - vR))/(rhoL*(sL - vL) - rhoR*(sR - vR));

    //Intermediate states
    double rhoStarL = rhoL*(sL - vL)/(sL - sStar);
    double rhoStarR = rhoR*(sR - vR)/(sR - sStar);
    double rhoVStarL = rhoStarL*sStar;
    double rhoVStarR = rhoStarR*sStar;
    double pStarL = sStar*rhoL*(sL- vL) - rhoVL*(sL - vL) + pL;
    double pStarR = sStar*rhoR*(sR - vR) - rhoVR*(sR - vR) + pR;
    double pStar  = 0.5*(pStarL + pStarR); 
    double EstarL = (sL*EL - (EL + pL)*vL + sStar*pStar)/(sL - sStar);
    double EstarR = (sR*ER - (ER + pR)*vR + sStar*pStar)/(sR - sStar);
    double rhoL_starL = rhoStarL*lamdaL;
    double rhoL_starR = rhoStarR*lamdaR;

    uStarL[0] = rhoStarL;
    uStarL[1] = rhoVStarL;
    uStarL[2] = EstarL;
    uStarR[0] = rhoStarR;
    uStarR[1] = rhoVStarR;
    uStarR[2] = EstarR;

    uStarL[3] = rhoL_starL;
    uStarR[3] = rhoL_starR;


    //Fluxes
    for (int i = 0; i < 4; i++){
        if (0 <= sL){
            f_iPlusHalf_HLLC[i] = getFluxFunction(uL)[i];
        }
        if (sL < 0 && 0 <= sStar){
            f_iPlusHalf_HLLC[i] = getFluxFunction(uL)[i] + sL*(uStarL[i] - uL[i]);
        }
        if (sStar < 0 && 0 <= sR){
            f_iPlusHalf_HLLC[i] = getFluxFunction(uR)[i] + sR*(uStarR[i] - uR[i]);
        }
        if (0 > sR){
            f_iPlusHalf_HLLC[i] = getFluxFunction(uR)[i];
        }
    }
    return f_iPlusHalf_HLLC;
}

double minMod(double r1, double r2) {
    if (r1 * r2 <= 0.0) {
        return 0.0;  // If the slopes have opposite signs, return zero
    }
    return (std::fabs(r1) < std::fabs(r2)) ? r1 : r2;  // Return the smaller magnitude of the two slopes
}

double minbeeLim(double r){
    const double eps = 1e-12;  // Small tolerance for floating-point comparisons

    // Check for NaN or Inf
    if (std::isnan(r) /*|| std::isinf(r)*/) {
        return 0.0;  // Fallback to first-order (no limiter)
    }
    if (std::isinf(r)){
        return 1;
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
        return 1.0;  // Avoid division by zero
    }
    return std::min(1.0, 2.0 / denominator);
}

double calculateSlopeRatio(double numerator, double denominator) {
    const double eps = 1e-12;  
    double r = numerator/denominator;
    // Avoid division by zero
    if (std::abs(denominator) < eps) {
        return 0.0;
    }
    // Check for NaN or Inf in inputs 
    if (std::isnan(r) || std::isinf(r)) {
        return 10; 
    }
    return r;
}

std::array<double,4> computeSourceTerm(const std::array<double, 4>& u) {
    std::array<double,4> sourceTerms;
    std::array<double,4> prim = ConsToPrim(u);
    double gamma = getPars()[g];
    double epsilon = 1/getPars()[eA];
    double Q = getPars()[heat];
    double theta = prim[p]/prim[rho_p];
    sourceTerms[rho] = 0;
    sourceTerms[rhoV] = 0;
    sourceTerms[E] = 1/(gamma-1) * u[rhoL] * epsilon * exp(1/epsilon * (1 - 1/theta));
    sourceTerms[rhoL] = -u[rhoL]/Q * epsilon * exp(1/epsilon * (1 - 1/theta));
    return sourceTerms;
}

std::array<double,4> updateRK4(std::array<double,4>& u, double dt) {
    std::array<double, 4> uUpdated = u;
    // k1 calculatiom
    std::array<double,4> k1;
    for (int var = 0; var < 4; var++){
        k1[var] = dt*computeSourceTerm(u)[var];
    }

    //std::cout << "k1[E]: " << k1[E] << std::endl;

    // k2 calculation
    std::array<double,4> u2;
    for (int var = 0; var < 4; var++){
        u2[var] = u[var] + 0.5 * k1[var];
    }
    std::array<double,4> k2;
    for (int var = 0; var < 4; var++){
        k2[var] = dt*computeSourceTerm(u2)[var];
    }

    //std::cout << "k2[E]: " << k2[E] << std::endl;

    // k3 calculation
    std::array<double,4> u3;
    for (int var = 0; var < 4; var++){
        u3[var] = u[var] + 0.5 * k2[var];
    }
    std::array<double,4> k3;
    for (int var = 0; var < 4; var++){
        k3[var] = dt*computeSourceTerm(u3)[var];
    }

    //std::cout << "k3[E]: " << k3[E] << std::endl;

    // k4 calculation
    std::array<double,4> u4;
    for (int var = 0; var < 4; var++){
        u4[var] = u[var] + k3[var];
    }
    std::array<double,4> k4;
    for (int var = 0; var < 4; var++){
        k4[var] = dt*computeSourceTerm(u4)[var];
    }
        
    //std::cout << "k4[E]: " << k4[E] << std::endl;

    // Final RK4 update
    for (int var = 0; var < 4; var++){
        uUpdated[var] +=  (1.0 / 6.0) * (k1[var] + 2.0 * k2[var] + 2.0 * k3[var] + k4[var]);
    }

    //std::cout << "E before RK4: " << u[E] << ", E after RK4: " << uUpdated[E] << std::endl;

    return uUpdated;
}

std::array<double, 4> updateODE(){
    std::array<double, 4> prim_u = {1.0, 0.0, 1.0, 1.0};
    std::array<double, 4> u = PrimToCons(prim_u);

    double dt = 1e-5; 
    double t = getPars()[tStart];
    double tStop = getPars()[tEnd];
    double gamma = getPars()[g];

    std::ofstream output("IgnitionFull.dat");
    output << "t rho v p lambda internal_energy theta\n";
    int dtcounter = 0;

    do{
        u = updateRK4(u, dt);
        t += dt;
        dtcounter++;
        double pressure = ConsToPrim(u)[p];
        double theta = pressure / u[rho];
        double internal_energy = u[E] - 0.5 * u[rhoV] * u[rhoV] / u[rho];
        double velocity = ConsToPrim(u)[v];
        double lambda = ConsToPrim(u)[l];
        output << t << " " << u[rho] << " " << velocity << " " << pressure << " " << lambda << " " 
                << internal_energy << " " << theta << "\n";
  
    }while(t < tStop);
    return u;
}

std::vector<std::array<double,4>> updateODEinSpace(std::vector<std::array<double,4>>&u, double dtCFL){
    int nCells = getPars()[xCells]; 
    double x0 = getPars()[xStart]; 
    double x1 = getPars()[xEnd];
    double t0 = 0;
    double tStop = dtCFL;
    double dx = (x1 - x0) / nCells;
    double gamma = getPars()[g];
    double Q = getPars()[heat];
    double e = 1/getPars()[eA];
    double dt = getPars()[dt_s];
    int dtcounter = 0;
    double t = t0;
    do{
        t += dt;
        dtcounter++;
        for (int i = 2; i < nCells+2; i++){
            u[i] = updateRK4(u[i], dt);
        }
    }while(t < dtCFL);

    return u;
}

void saveTimeframe(double t, std::vector<std::array<double,4>>u) {
    double nCells = getPars()[xCells];
    double x0 = getPars()[xStart]; 
    double x1 = getPars()[xEnd];
    double dx = (x1 - x0) / nCells;
    
    std::ostringstream filename;
    filename << "results_t" << std::fixed << std::setprecision(3) << t << "_10000cellsdt1e-5_Q2.dat";

    std::ofstream file(filename.str());
    file << "# t x density velocity pressure Lamda int.energy T\n";

    for (int i = 1; i < nCells+2; i++) {
        double x = x0 + (i-1)*dx;
        std::array<double, 4> prim = ConsToPrim(u[i]);
        double internalEnergy = PrimToCons(prim)[E] - 0.5*prim[rho_p]*prim[v]*prim[v];
        double temp = prim[p]/prim[rho_p];

        file << t << " " << x << " " << prim[0] << " " << prim[1] << " " \
        << prim[2] << " " << prim[3] << " " << internalEnergy << " " << temp <<  "\n";
    }
    file.close();
}



int updatePDESubcycling(){
    int nCells = getPars()[xCells]; 
    double x0 = getPars()[xStart]; 
    double x1 = getPars()[xEnd];
    double tStart = getPars()[tStart];
    double tStop = getPars()[tEnd];
    double dx = (x1 - x0) / nCells;
    double gamma = getPars()[g];
    double Q = getPars()[heat];
    double e = 1/getPars()[eA];
    double dt; // hyperbolic time step

    std::vector<std::array<double, 4>> u; // u vector
    u.resize(nCells+4);
    
    std::vector<std::array<double, 4>> flux; // Fluxes array: 3 fluxes at each i
    flux.resize(nCells+3);
 
    std::vector<std::array<double, 4>> uPlus1; //time-updated u vector
    uPlus1.resize(nCells+4);

    //Initial Data
    for(int i = 0; i < u.size(); i++){
        double x = x0 + (i-0.5) * dx;
        std::array<double, 4> prim;
        prim[0] = 0.259259; // Density
        prim[1] = -1.35769; // Velocity
        prim[2] = 0.0967742; // Pressure
        prim[3] = 1; // Lamda
        u[i] = PrimToCons(prim);
    }
    //Debugging print of initial data
    std::cout << "Initial data:" << std::endl;
    for (int i = 0; i < u.size(); i++){
        std::cout << "Cell " << i << ": "<< ConsToPrim(u[i])[0] << " " << ConsToPrim(u[i])[1] << " " \
        << ConsToPrim(u[i])[2] << " " << ConsToPrim(u[i])[3] << std::endl;
    }

    int debuggingInitialData = 0; //set to 1 for constant data in all time steps
    double t = tStart;
    int dt_counter = 0;
    int debuggingLimiter = 0; //set to 1 to set all limiters to 0 ie revert to force

    double dt_save;
    double next_save_time;

    //Declare vectors for reconstrcuted states
    std::vector<std::array<double, 4>> uL;
    uL.resize(nCells+4);
    std::vector<std::array<double, 4>> uR;
    uR.resize(nCells+4);

    //Declare vector to store the limited reconstrcuted states
    std::vector<std::array<double,4>> uBarL;
    uBarL.resize(nCells+4);
    std::vector<std::array<double,4>> uBarR;
    uBarR.resize(nCells+4);

    //Declare vectors for half-timestep updated reconstructed states
    std::vector<std::array<double,4>> uBarLUpdate;
    uBarLUpdate.resize(nCells+4);
    std::vector<std::array<double,4>> uBarRUpdate;
    uBarRUpdate.resize(nCells+4);

    //Declare vecotr that stores the flux function f(u) for each limited variable to use at half time step update
    std::vector<std::array<double,4>> fluxFunction;
    fluxFunction.resize(nCells+4);
    do{
        dt = computeHyperbolicTimeStep(u);
        t = t + dt;
        dt_counter++;
        std::cout<<"t = " << t << ", dt = "  << dt << " timestep no: " << dt_counter <<std::endl;

        if (t >= 1.02 && next_save_time > 1.02) {
            dt_save = 0.01;
            next_save_time = 1.02;  // Reset save time 
        } 
        else if (t >= 1.13 && next_save_time > 1.13) {
            dt_save = 0.02;
            next_save_time = 1.13;  // Reset save time 
        }

        u[1][rhoV] = -u[2][rhoV];
        u[0][rhoV] = u[1][rhoV];
        u[nCells+2][rhoV] = u[nCells+1][rhoV];
        u[nCells+3][rhoV] = u[nCells+2][rhoV];

        u[1][rho]= u[2][rho];
        u[0][rho] = u[1][rho];
        u[nCells+2][rho] = u[nCells+1][rho];
        u[nCells+3][rho] = u[nCells+2][rho];

        u[1][E]= u[2][E];
        u[0][E] = u[1][E];
        u[nCells+2][E] = u[nCells+1][E];
        u[nCells+3][E] = u[nCells+2][E];

        u[1][rhoL]= u[2][rhoL];
        u[0][rhoL] = u[1][rhoL];
        u[nCells+2][rhoL] = u[nCells+1][rhoL];
        u[nCells+3][rhoL] = u[nCells+2][rhoL];        


         //Apply limiter to get uBarL for cells i to nCells
        for (int i = 1; i < nCells+3; i++){
            for (int j = 0; j < 4; j++){
                double r_num = u[i][j] - u[i-1][j];
                double r_den = u[i+1][j] - u[i][j];
                double r = calculateSlopeRatio(r_num,r_den);
                //double di = minMod(u[i+1][j] - u[i][j], u[i][j] - u[i-1][j]);
                double di = 0.5*(u[i+1][j] - u[i-1][j]);
                double lim;
                if (debuggingLimiter == 0){
                    lim = minbeeLim(r);
                }
                if (debuggingLimiter == 1){
                    lim = 0;
                }
                uBarL[i][j] = u[i][j] - 0.5*lim*di; 
                uBarR[i][j] = u[i][j] + 0.5*lim*di;         
            }
        }
        //Apply half time step update to get uBarLUpdate
        for (int i = 1; i < nCells+3; i++){
            for (int j = 0;  j < 4; j++){
                uBarLUpdate[i][j] = uBarL[i][j] - 0.5*(dt/dx)*(getFluxFunction(uBarR[i])[j] - getFluxFunction(uBarL[i])[j]);
            }
        }
        //Apply half time step update to get uBarRUpdate
        for (int i = 0; i < nCells+3; i++){
            for (int j = 0;  j < 4; j++){
                uBarRUpdate[i][j] = uBarR[i][j] - 0.5*(dt/dx)*(getFluxFunction(uBarR[i])[j] - getFluxFunction(uBarL[i])[j]);
            }
        }
        //Fill flux array
        for (int i = 1; i < nCells+2; i++){
            flux[i] = getfluxHLLC(uBarRUpdate[i],uBarLUpdate[i+1]);
        }

        if (debuggingInitialData == 1){ //this is for tests 1 and 2 only
            for (int i = 2; i < nCells+1; i++){
                for (int j = 0; j < 3; j++){
                    uPlus1[i][j] = u[i][j];
                }
            }
        }
        if (debuggingInitialData == 0){ //real update
        // Hyperbolic update
        for (int i = 2; i < nCells+2; i++){
                for (int j = 0; j < 4; j++){
                    uPlus1[i][j] = u[i][j] - (dt/dx)*(flux[i][j] - flux[i-1][j]);
                }
            }
        }

        uPlus1 = updateODEinSpace(uPlus1, dt);

        uPlus1[1][rhoV] = -uPlus1[2][rhoV];
        uPlus1[0][rhoV] = uPlus1[1][rhoV];
        uPlus1[nCells+2][rhoV] = uPlus1[nCells+1][rhoV];
        uPlus1[nCells+3][rhoV] = uPlus1[nCells+2][rhoV];

        uPlus1[1][rho]= uPlus1[2][rho];
        uPlus1[0][rho] = uPlus1[1][rho];
        uPlus1[nCells+2][rho] = uPlus1[nCells+1][rho];
        uPlus1[nCells+3][rho] = uPlus1[nCells+2][rho];

        uPlus1[1][E]= uPlus1[2][E];
        uPlus1[0][E] = uPlus1[1][E];
        uPlus1[nCells+2][E] = uPlus1[nCells+1][E];
        uPlus1[nCells+3][E] = uPlus1[nCells+2][E];

        uPlus1[1][rhoL]= uPlus1[2][rhoL];
        uPlus1[0][rhoL] = uPlus1[1][rhoL];
        uPlus1[nCells+2][rhoL] = uPlus1[nCells+1][rhoL];
        uPlus1[nCells+3][rhoL] = uPlus1[nCells+2][rhoL]; 

        u = uPlus1;
        
        /*if (t >= next_save_time) {
            saveTimeframe(t, u);
            next_save_time += dt_save;
        }*/
        const double epsilon = 1e-10;
        if (t + epsilon >= next_save_time) {
            saveTimeframe(t, u);
            next_save_time += dt_save;
        }
    }while(t < tStop);

    std::cout << "number of timesteps: " << dt_counter << std::endl;

    //Convert back to primitive variable before outputting to file
    for (int i = 0; i < nCells+4; i++){
        u[i] = ConsToPrim(u[i]);
    }

    //Output the data
    std::ofstream output("IgnitionFull.dat");
    output << "x rho v p lambda e T\n";
    for(int i = 1; i < nCells+2; i++) {
        double x = x0 + (i-1)*dx;
        double internalEnergy = PrimToCons(u[i])[E] - 0.5*u[i][rho_p]*u[i][v]*u[i][v];
        double pressure = u[i][p];
        double temp = pressure / (u[i][rho_p]);
        output << x << " " << u[i][0] <<" "<< u[i][1] << " "<< u[i][2] << " " << u[i][3] << \
        " " << internalEnergy << " " << temp << std::endl;
    }
    output.close(); // Close the file

    return 0;
}

int main(){
    auto start = std::chrono::high_resolution_clock::now();
    updatePDESubcycling();
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << "Execution time: " << duration.count() << " ms" << std::endl;
    return 0;
}

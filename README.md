<pre>
# Mars_Ice_Sheets (unparallelized - working on creating parallel version)  
Numerical model simulating the deposition, compaction, and thermal evolution of Martian ice sheets  

The main file to run the script (and make edits in to tailor to your specific problem) is Mars_Ice_HiRes_Run.m  

In this file you can change several input parameters for your specific simulation:  

%% Input Variables  
output_t=5000;            % How often to output/save data (iterations/timesteps)  
tf=31557600*5000;         % Run Duration (s) this will be final time of run  
Ttop=Gale_Temp;           % Surface Temp (K) - Gale_Temp is a surface temperature matrix  
                          % that describes the surface temperature (K) for each grid cell of the simulation region  
Tbottom=230.0;            % Basal Temp - Initial Temperature (K) of the bottom of the ice sheet  
Geo_Flux=Gale_GeoFlux;    % Basal Heat Flux (W/m^2) - Gale_Geoflux is a matrix of the basal geothermal heat flux  
                          % for each grid cell of the simulation region  
k_i=2;                    % Ice Thermal Conductivity (W/m*K)  
k_br=0.6;                 % Water Thermal Conductivity (W/m*K)  
k_v=0.012;                % Void Space Thermal Conductivity (W/m*K)  
dt=31557600*1000;         % Time Step (s) (31557600s=1yr)  
dz=5;                     % Spatial Discretization (m)  
H_0=2*dz;                 % Initial Ice Thickness (m) - easiest to start simulation with thin ice sheet at the Tbottom temperature  
c_i=2000;                 % Specific Heat of Ice (J/kg*K)  
c_br=3985;                % Specific Heat of Water (J/kg*K)  
c_v=790;                  % Specific Heat of Voids in Ice (J/kg*K)  
L=334774;                 % Latent Heat of Fusion (J/kg)  
rho_i=917;                % Density of Ice (Kg/m^3)  
rho_br=1000;              % Density of Water (Kg/m^3)  
rho_v=1;                  % Density of Voids in Ice (Kg/m^3)  
Tm=273.15;                % Melting Temperature of Ice (K)  
TTol=0.01;                % Temperature Tolerance (K)  
PhiTol=0.01;              % Liquid Fraction Tolerance  
reservoir=2.8208*10^16;   % Total Available Water (kg) available to be deposited as ice  
res=463;                  % Grid resolution (m)  

% Deposition parameters following Cassanelli and Head (2015)  
b=Gale_Depo/31557600;     % Deposition Rate (m/s) - Gale_Depo is a matrix of snowfall deposition in (m/yr) for each  
                          % grid cell of the simulation - here it is divided by 31557600 to get (m/s - SI units)  
Q=45600;                  % Activation Energy (J/mol)  
R=8.314;                  % Gas Constant (J/mol*K)  
g=3.711;                  % Gravity (m/s^2)  
n1=3;                     % Power Law for Ice Deformation  
f_0=0.00000003/31557600;  % Constant Coefficient (couldn't find actual value in papers - this is from my own fit to  
                          %    Siple Dome data and b in m/s)  
rho_1=350;                % Freshly Deposited Ice Density (Kg/m^3)  

-------------------------------------------
The primary outputs of the code will be:  

melts3D         - A 3D matrix that catalogues the amount of basal melt produced (m) at each grid cell throughout the simulation,  
                  that is an Nx x Ny X Nz matrix where Nx and Ny are spatial coordinates and Nz is the temporal steps of the simulation  
                  at intervals of output_t  
H3D             - A 3D matrix that catalogues the thickness of the ice sheet (m) at each grid cell throughout the simulation, same Nx x Ny x Nz  
                  format as melts3D, but Nz is ice thickness  
res_mat         - Catalogues the remaining water reservoir available for deposition onto the ice sheet throughout the simulation - in temporal  
                  intervals of output_T  
Temperature     - Final vertical temperature profile throughout the ice sheet, Nx x Ny x Nz matrix where Nx and Ny are spatial coordinates and   
                  Nz is depth withing the ice sheet  
rho             - Final vertical density profile throughout the ice sheet, Nx x Ny x Nz matrix where Nx and Ny are spatial coordinates and   
                  Nz is depth withing the ice sheet  
Liquid_Fraction - Final vertical liquid fraction profile throughout the ice sheet, Nx x Ny x Nz matrix where Nx and Ny are spatial coordinates and   
                  Nz is depth withing the ice sheet  

-------------------------------------------
Accompanying paper can be found at: https://www.sciencedirect.com/science/article/pii/S0012821X22003351  
</pre>

%% Platform for running Mars Ice code
%  v3 has latitudinally and altidudinally varying surface temp, altidude
%  dependent ice deposition, and a supply limited ice deposition (Many of
%  these functions are closely related to the equations of Fastook and Head
%  (2012))

% Make sure the following input parameters are set right:
% tf - final time in sec
% Ttop and Tbottom - should be the desired surface temp in K
% Geoflux - make sure the right resolution nhf is used
% b - this is the deposition rate in m/s
% Htot - total ice thickness in m (param sweep in parfor loop)
% nhf_deg - nhf input
% nhf_res - resolution of nhf used (used to name saved matricies)
% b_rate - deposition rate in cm/yr
% heat_flux - mantle heat flux in mW

tf=31557600*10000000;       % Run Duration (s) this will be final time of 
                            % run (lower down will determine deposition vs.
                            % compaction

%% Input Variables
Ttop=Tsurf_Obliq25_5deg_alt;           % Surface Temp (K)
Tbottom=230.0;      % Basal Temp - Initial (K)
Geo_Flux=nhf_20mW_5deg/1000;     % Basal Heat Flux (W/m^2)
k_i=2;              % Ice Thermal Conductivity (W/m*K)
k_br=0.6;           % Water Thermal Conductivity (W/m*K)
k_v=0.012;              % Void Space Thermal Conductivity (W/m*K)
dt=31557600*1000;       % Time Step (s) (31557600s=1yr)
dz=5;               % Spatial Discretization (m)
H_0=2*dz;             % Initial Ice Thickness (m)
c_i=2000;           % Specific Heat of Ice (J/kg*K)
c_br=3985;          % Specific Heat of Water (J/kg*K)
c_v=790;              % Specific Heat of Voids in Ice (J/kg*K)
L=334774;           % Latent Heat of Fusion (J/kg)
rho_i=917;          % Density of Ice (Kg/m^3)
rho_br=1000;        % Density of Water (Kg/m^3)
rho_v=1;            % Density of Voids in Ice (Kg/m^3)
Tm=273.15;          % Melting Temperature of Ice (K)
TTol=0.01;          % Temperature Tolerance (K)
PhiTol=0.01;        % Liquid Fraction Tolerance
reservoir=10*(5*10^18);  % Total Available Water (kg)

% Deposition parameters following Cassanelli and Head (2015)
%b=[0.001,0.01,0.1,1];             % Deposition Rate (m/yr) - CONSTANT
%b=[0.001/31557600,0.01/31557600,0.1/31557600];      % Deposition Rate (m/s)
b=Depo_mat_5deg/31557600;        % Deposition Rate (m/s)
Q=45600;                % Activation Energy (J/mol)
R=8.314;                % Gas Constant (J/mol*K)
g=3.711;                % Gravity (m/s^2)
n1=3;                    % Power Law for Ice Deformation
f_0=0.00000003/31557600;   % Constant Coefficient (couldn't find actual value in papers - this is from my own fit to
                        %    Siple Dome data and b in m/s)
rho_1=350;              % Freshly Deposited Ice Density (Kg/m^3)

%% Initial Density, Liquid Fraction, Temperature, and Void Space Profiles
% Temperature=Ttop+((Tbottom-Ttop)/H_0)*[dz:dz:H_0]';
% rho=Ice_Density(H_0,b,Temperature,g,f_0,dz,rho_1);
% Liquid_Fraction=0*[dz:dz:H_0]';
% phi_v=1-(rho'./rho_i);
% 
% H=H_0;
% hold_snow=0;

%% Predefine Plotting Matrices for Speed Using Deposition Rate (b) and Run Time (tf)
% Tplot=Ttop+zeros(1500/dz,tf/dt);
% Phiplot=zeros(1500/dz,tf/dt);
% pc_plot=zeros(1500/dz,tf/dt);
% k_plot=zeros(1500/dz,tf/dt);
% rho_plot=zeros(1500/dz,tf/dt);

%% Sweeping through total ice thickness
%Htot=[1300,1400,1500,2000];
%Htot=2000;

%% Empty matrix for collecting melt thickness results
%melts=Geo_Flux;
melt_hold=zeros(size(Geo_Flux,1),size(Geo_Flux,2),tf/(31557600*100000));
h_hold=zeros(size(Geo_Flux,1),size(Geo_Flux,2),tf/(31557600*100000));

%% Resolution and stats for naming matricies during loops
nhf_deg=nhf_20mW_5deg;
nhf_res=5;
b_rate=1;
heat_flux=20;

%     parfor jj=1:length(Htot)
%         j=Htot(jj);
        
        %% start function for parallel loops
        [melts3D,H3D,res_mat]=Mars_parloop_v3(heat_flux,nhf_deg,nhf_res,Geo_Flux,H_0,Ttop,Tbottom,dz,b,g,f_0,rho_1,...
            rho_i,dt,tf,Q,R,n1,Tm,k_i,k_br,k_v,c_i,c_br,c_v,L,TTol,PhiTol,rho_br,rho_v,melt_hold,reservoir,h_hold);
        %% end of parallel loop
%     end

figure
subplot(2,2,1)
image(1000*Geo_Flux,'CDataMapping','scaled')
colorbar
title('Basal Heat Flux (mW)')
subplot(2,2,2)
image(melts3D(:,:,end),'CDataMapping','scaled')
colorbar
title('Total Basal Melt (m)')
subplot(2,2,3)
image(H3D(:,:,end),'CDataMapping','scaled')
colorbar
title('Ice Thickness (m)')
subplot(2,2,4)
plot((dt/31557600)*[1:length(res_mat)],res_mat)
title('Total H2O Reservoir (kg)')
xlabel('Time (years)')
ylabel('Available H2O (kg)')

% figure
% image(melts,'CDataMapping','scaled')
% colorbar
% caxis([0 max(max(melts3D(:,:,end)))])

% figure
% subplot(2,2,1)
% image(Tplot,'CDataMapping','scaled')
% title('Temperature (K)')
% xlabel(strcat('Time (',num2str(dt/31557600),'yr)'))
% ylabel(strcat('Depth (',num2str(dz),' m)'))
% colorbar
% subplot(2,2,2)
% image(Phiplot,'CDataMapping','scaled')
% title('Liquid Fraction')
% xlabel(strcat('Time (',num2str(dt/31557600),'yr)'))
% ylabel(strcat('Depth (',num2str(dz),' m)'))
% colorbar
% subplot(2,2,3)
% image(rho_plot,'CDataMapping','scaled')
% title('Density (kg/m^3)')
% xlabel(strcat('Time (',num2str(dt/31557600),'yr)'))
% ylabel(strcat('Depth (',num2str(dz),' m)'))
% colorbar
% subplot(2,2,4)
% image(k_plot./pc_plot,'CDataMapping','scaled')
% title('Thermal Diffusivity [k/rho*c]')
% xlabel(strcat('Time (',num2str(dt/31557600),'yr)'))
% ylabel(strcat('Depth (',num2str(dz),' m)'))
% colorbar
% 
% figure
% plot([dt:dt:tf],dz*sum(Phiplot))
% title('Water Depth vs. Time')
% xlabel('Time (s)')
% ylabel('Water Depth (m)')

figure
subplot(1,3,1)
hold on
hold1=zeros(1,size(H3D,3));
for i=1:size(H3D,1)
    for j=1:size(H3D,2)
        hold1(1:size(H3D,3))=H3D(i,j,:);
        plot(1:size(H3D,3),hold1)
    end
end
title('Ice Height (m)')
subplot(1,3,2)
hold on
hold1=zeros(1,size(H3D,3));
for i=1:size(H3D,1)
    for j=1:size(H3D,2)
        hold1(1:size(H3D,3))=melts3D(i,j,:);
        plot(1:size(H3D,3),hold1)
    end
end
title('Total Basal Melt (m)')
subplot(1,3,3)
hold1=zeros(1,size(H3D,3)-1);
for i=2:size(H3D,3)
    hold1(i-1)=sum(sum(melts3D(:,:,i)))-sum(sum(melts3D(:,:,i-1)));
end
plot(1:length(hold1),hold1)
title('Melt Rate (m/dt)')
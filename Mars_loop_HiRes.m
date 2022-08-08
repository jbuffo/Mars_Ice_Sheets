%% loop for MarsIce_Tsurf_vary that simulates latitudinally varying surface temperature
function [melts3D,H3D,res_mat,Temperature,rho,Liquid_Fraction]=Mars_loop_HiRes(heat_flux,save_name,Geo_Flux,H_0,...
Ttop,Tbottom,dz,b,g,f_0,rho_1,rho_i,dt,tf,Q,R,n1,Tm,k_i,k_br,k_v,c_i,c_br,c_v,L,TTol,PhiTol,rho_br,rho_v,melt_hold,...
reservoir,h_hold,res)

% Initial Conditions
t=0;
H=(b>0)*H_0;
hold_snow=zeros(size(Geo_Flux,1),size(Geo_Flux,2));
new_snow=zeros(size(Geo_Flux,1),size(Geo_Flux,2));
melt_tot=zeros(size(Geo_Flux,1),size(Geo_Flux,2));
Temperature=zeros(size(Geo_Flux,1),size(Geo_Flux,2),10000/dz);
rho=zeros(size(Geo_Flux,1),size(Geo_Flux,2),10000/dz);
Liquid_Fraction=zeros(size(Geo_Flux,1),size(Geo_Flux,2),10000/dz);
for m=1:size(Geo_Flux,1)
    for n=1:size(Geo_Flux,2)
        if b(m,n)>0
            Temperature(m,n,1:H(m,n)/dz)=Ttop(m,n)+((Tbottom-Ttop(m,n))/H(m,n))*[dz:dz:H(m,n)]';   % Conductive temp profile
            rho(m,n,1:H(m,n)/dz)=Ice_Density(H(m,n),b(m,n),Temperature(m,n,1:H(m,n)/dz),g,f_0,dz,rho_1);
        else
        end
    end
end
phi_v=1-(rho/rho_i);    % phi_v is void space (air/vacuum)

reservoir=reservoir-sum((dz*res^2).*rho,'all');
res_mat=zeros(1,tf/dt);

%% Discretization Matrices
up_avg_large=full(gallery('tridiag',10000/dz,1,1,0));

down_avg_large=full(gallery('tridiag',10000/dz,0,1,1));

up_disc_large=full(gallery('tridiag',10000/dz,1,-1,0));

down_disc_large=full(gallery('tridiag',10000/dz,0,1,-1));

%% Start loop over time
while t<tf
    t=t+dt;
    % initiate empty set to track deposition to subtract from reservoir
    res_sub=0;
    res_add=0;
    for m=1:size(Geo_Flux,1)
        for n=1:size(Geo_Flux,2)
            k=Geo_Flux(m,n);            

            if b(m,n)==0
                melts=0;
                ice_height=0;
            else

                % Melt=0 if no geothermal heat flux data
                if k==0
                    cells=0;
                else  

                %% Deposition of new snow onto ice sheet (or compaction) and calculation of density profile
                if reservoir>sum(rho_i*res^2*dt.*b,'all')
                    new_snow(m,n)=dt*b(m,n);
                    hold_snow(m,n)=hold_snow(m,n)+new_snow(m,n);
                    if hold_snow(m,n)>=dz
                        new_layers=floor(hold_snow(m,n)/dz);
                        hold_snow(m,n)=hold_snow(m,n)-dz*new_layers;
                        H(m,n)=H(m,n)+dz*new_layers;
                        res_sub=res_sub+dz*new_layers*rho_i*res^2;
                        holdT=[];
                        holdT(1:(H(m,n)/dz)-new_layers)=Temperature(m,n,1:(H(m,n)/dz)-new_layers);
                        Temperature(m,n,1:H(m,n)/dz)=[Ttop(m,n)*(1+0*[1:new_layers]) holdT];
                        holdrho=[];
                        holdrho(1:H(m,n)/dz-new_layers)=rho(m,n,1:H(m,n)/dz-new_layers);
                        rho(m,n,1:H(m,n)/dz)=[rho_1*(1+0*[1:new_layers]) compact(holdrho,dt,f_0,Q,R,holdT,rho_i,g,dz,n1)];
                        holdPhi=[];
                        holdPhi(1:(H(m,n)/dz)-new_layers)=Liquid_Fraction(m,n,1:(H(m,n)/dz)-new_layers);
                        Liquid_Fraction(m,n,1:H(m,n)/dz)=[0*[1:new_layers] holdPhi];
                        phi_v(m,n,1:H(m,n)/dz)=1-(rho(m,n,1:H(m,n)/dz)./rho_i);
                    else
                    end
                else
                    holdrho=[];
                    holdrho(1:H(m,n)/dz)=rho(m,n,1:H(m,n)/dz);
                    holdT=[];
                    holdT(1:H(m,n)/dz)=Temperature(m,n,1:H(m,n)/dz);
                    rho(m,n,1:H(m,n)/dz)=compact(holdrho,dt,f_0,Q,R,holdT,rho_i,g,dz,n1);
                    phi_v(m,n,1:H(m,n)/dz)=1-(rho(m,n,1:H(m,n)/dz)./rho_i);
                end

                %% Discretization Matrices
                up_avg=up_avg_large(H(m,n)/dz,H(m,n)/dz);
                up_avg(1,1)=2;

                down_avg=down_avg_large(H(m,n)/dz,H(m,n)/dz);
                down_avg(H(m,n)/dz,H(m,n)/dz)=2;

                up_disc=up_disc_large(H(m,n)/dz,H(m,n)/dz);

                down_disc=down_disc_large(H(m,n)/dz,H(m,n)/dz);
                down_disc(H(m,n)/dz,H(m,n)/dz)=0;

                %% Thermal Evolution of the Ice Sheet via Conduction
                hold_Phi=[];
                hold_Phi(1:(H(m,n)/dz))=Liquid_Fraction(m,n,1:(H(m,n)/dz));
                hold_T=[];
                hold_T(1:H(m,n)/dz)=Temperature(m,n,1:H(m,n)/dz);
                hold_Phiv=[];
                hold_Phiv(1:H(m,n)/dz)=phi_v(m,n,1:H(m,n)/dz);
                [NewTemp,NewPhi,pc_bar,k_bar,cells]=Mars_Ice_v2(...
                    k_i,k_br,k_v,dt,dz,H(m,n),c_i,c_br,c_v,L,TTol,PhiTol,hold_Phi',hold_T',...
                    hold_Phiv',rho_i,rho_br,rho_v,Tm,Ttop(m,n),k,up_avg,down_avg,up_disc,down_disc);

                melt_tot(m,n)=melt_tot(m,n)+cells*dz;
                res_add=res_add+dz*cells*rho_i*res^2;
                H(m,n)=H(m,n)-cells*dz;
                Temperature(m,n,1:H(m,n)/dz)=NewTemp(1:H(m,n)/dz);
                Liquid_Fraction(m,n,1:H(m,n)/dz)=NewPhi(1:H(m,n)/dz);

                %     Tplot(H_0/dz+ceil((b*(tf/31557600))/dz)-H/dz+1:H_0/dz+ceil((b*(tf/31557600))/dz),counter)=Temperature;
                %     Phiplot(H_0/dz+ceil((b*(tf/31557600))/dz)-H/dz+1:H_0/dz+ceil((b*(tf/31557600))/dz),counter)=Liquid_Fraction;
                %     pc_plot(H_0/dz+ceil((b*(tf/31557600))/dz)-H/dz+1:H_0/dz+ceil((b*(tf/31557600))/dz),counter)=pc_bar;
                %     k_plot(H_0/dz+ceil((b*(tf/31557600))/dz)-H/dz+1:H_0/dz+ceil((b*(tf/31557600))/dz),counter)=k_bar;
                %     rho_plot(H_0/dz+ceil((b*(tf/31557600))/dz)-H/dz+1:H_0/dz+ceil((b*(tf/31557600))/dz),counter)=rho;

                %     if ismember(t,[31557600*100000:31557600*100000:tf])==1
                %         melt_hold(m,n,t/(31557600*100000))=melt_tot(m,n);
                %         h_hold(m,n,t/(31557600*100000))=H(m,n);
                %     else
                %     end
                if ismember(t,[31557600*5000:31557600*5000:tf])==1
                    melt_hold(m,n,t/(31557600*5000))=melt_tot(m,n);
                    h_hold(m,n,t/(31557600*5000))=H(m,n);
                else
                end

                end
            end
        end
        % melts(m,n)=cells*dz;
    end

reservoir=reservoir-res_sub+res_add;
res_mat(t/dt)=reservoir;
end
        
melts3D=melt_hold;
H3D=h_hold;
save(strcat(save_name,num2str(heat_flux),'mW_',...
    '10My_by_100ky.mat'),'melts3D','H3D','res_mat','Temperature','rho','Liquid_Fraction')
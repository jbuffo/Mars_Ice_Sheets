function [T_np1_k_j,phi_np1_k_j,pc_bar,k_bar,cells]=Mars_Ice_v2(...
    k_i,k_br,k_v,dt,dz,H,c_i,c_br,c_v,L,TTol,PhiTol,phi_initial,T_initial,...
    phi_v,rho_i,rho_br,rho_v,Tm,Ttop,Geo_Flux,up_avg,down_avg,up_disc,down_disc)

%% Initial Condition Vectors for Phi, T, S, w
phi_initial=phi_initial(:);
T_initial=T_initial(:);
phi_np1_km2_j=phi_initial(:);
T_np1_km1_j=T_initial(:);

T_evolve=[]';
phi_evolve=[]';
TErr=TTol+1;
PhiErr=PhiTol+1;
counter=0;
matrix_dimension=H/dz;


%% Itterattions over k
while TErr>TTol || PhiErr>PhiTol;
     counter=counter+1;
%     if counter>25
%         disp(counter)
%     else
%     end
    
%   delta_T=1.853*(S_np1_km1_j/28);        %% Old Freezing point depression
%   due to salt (linear approx)

%     delta_T=Tm-(-(1.333489497*(10^-5)*S_np1_km1_j.^2)-0.01612951864*S_np1_km1_j+...
%         273.055175687);                       %% Liquidus point using FREZCHEM Europa ocean
    
%   delta_T=Tm-(-(9.1969758*(10^-5)*S_np1_km1_j.^2)-0.03942059*S_np1_km1_j+...
%        272.63617665);                        %% Liquidus point using FREZCHEM seawater

delta_T=0*[dz:dz:H];

    Hs=c_i*(Tm-delta_T);                     %% Enthalpy Calculation
    Hsmelt=c_i*Tm;                           %% Enthalpy to melt ice
    %% Calculating solid fraction at k-1 iteration for use in
    %% temperature, salinity, and brine velocity at k iteration
    
    %% Value of Phi(n+1,k-1,j)
    phi_np1_km1_j=phi_np1_km2_j;
    for i=1:matrix_dimension;
    %% Test to see if melting ice or freezing brine
    En=c_i*T_np1_km1_j(i)+L*phi_np1_km2_j(i);
        if En<Hs(i)
            phi_np1_km1_j(i)=0;
        elseif En>Hs(i)+L*(1-phi_v(i));
            phi_np1_km1_j(i)=1-phi_v(i);
        else
            phi_np1_km1_j(i)=((c_i*T_np1_km1_j(i)+phi_np1_km2_j(i)*L)-Hs(i))/...
            L;
        end;
    end;

    %% Reassigning Phi(n+1,k-2) as Phi(n+1,k-1) for next
    %% iteration
    phi_np1_km2_j=phi_np1_km1_j;
    
    %% Volume averaged quantities
    pc_bar=rho_br*c_br*phi_np1_km1_j+rho_v*c_v*phi_v+rho_i*c_i*(1-phi_np1_km1_j-phi_v);
    %k_bar=k_br*phi_np1_km1_j+k_v*phi_v+k_i*(1-phi_np1_km1_j-phi_v);
    k_bar=k_br*phi_np1_km1_j+k_v*phi_v+((651./T_np1_km1_j)-0.2).*(1-phi_np1_km1_j-phi_v); % Temperature dependent thermal 
                                                                                          % conductivity for ice
    %% Solve for Temperature profile
    T_mat=eye(H/dz)-(dt./(pc_bar.*dz^2)).*(((1/2)*up_avg*k_bar).*up_disc-((1/2)*down_avg*k_bar).*down_disc);

    a=diag(T_mat);
    b=diag(T_mat,1);
    c=diag(T_mat,-1);
    
    y=T_initial-(rho_i*L./pc_bar).*(phi_np1_km1_j-phi_initial);
    
    y(1)=y(1)+(dt./(pc_bar(1)*dz^2))*k_bar(1)*Ttop;
    y(end)=y(end)+(dt./(pc_bar(end)*dz))*Geo_Flux;

    % Solve for T
    T_np1_k_j=Thomas_Trid(a',b',c',y);
%     T_np1_k_j=T_mat\y;
    
    T_np1_km1_j=T_np1_k_j;
    
    %% Appending value to matrix to check for convergence
    T_evolve=[T_evolve T_np1_k_j];
    phi_evolve=[phi_evolve phi_np1_km1_j];
    
    %% Reassign T vector to use in next k iteration
    T_np1_km1_j=T_np1_k_j;
    phi_np1_k_j=phi_np1_km1_j;
    % tolerance tests
    if counter==1;
        TErr=TErr;
        PhiErr=PhiErr;
    else
        TErr=max(abs(T_evolve(:,counter)-T_evolve(:,counter-1)));
        PhiErr=max(abs(phi_evolve(:,counter)-phi_evolve(:,counter-1)));
    end;
end

%% Mixing fully melted cells (assume mixing time is much less than time step dt)
cells=0;
hold_temp=0;
for i=H/dz:-1:1
    if phi_np1_km1_j(i)==1-phi_v(i)
        cells=cells+1;
        hold_temp=hold_temp+T_np1_km1_j(i);
    else
        break
    end
end

if cells>0
    T_np1_k_j(H/dz-cells+1:H/dz)=hold_temp/cells;
else
end


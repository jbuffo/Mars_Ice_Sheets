%% Compaction function FDE from Cassanelli and Head eq. (1)

function rho=compact(rho_in,dt,f_0,Q,R,T,rho_i,g,dz,n)

del_rho=[];
del_rho(1)=0;

for i=2:length(rho_in)
    del_rho(i)=f_0*rho_in(i)*exp(-Q/(R*T(i)))*...
        (((rho_i/rho_in(i))-1)*((rho_i*g*dz*sum(rho_in(1:i-1)))/rho_in(i)))^n;
end

del_rho;

rho=rho_in+dt*del_rho;

for i=1:length(rho)
    if rho(i)>rho_i
        rho(i)=rho_i;
    else
    end
end

% figure
% plot(rho_in,[-1:-dz:-length(rho_in)],rho,[-1:-dz:-length(rho_in)])
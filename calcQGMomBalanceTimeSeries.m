%%%
%%% calcMomBalanceTimeSeries.m
%%%
%%% Computes time series of wind stress and form stresses under the quasi-geostrophic approximation.
%%%

%%% Load constant parameters
constants;
Ny = 128;
Nlay = 3;
tau_mean = [0.15];
tau_pert = 0.075;
tau_freq = t1year * 2^4;
% AABW_mean = 1.5;
AABW_mean = 0;
AABW_pert = 0;
AABW_freq = 0;
quad_drag = 2e-3;
lin_drag = 0e-4;  
topog_width = 150;
topog_height = 1000;
rough_topog = true;
n_E = 2;
run_name = constructRunName (false,Ny,Nlay, ...
                                  tau_mean,tau_pert,tau_freq, ...
                                  AABW_mean,AABW_pert,AABW_freq, ...
                                  quad_drag,lin_drag,topog_width,topog_height,rough_topog,n_E);

% run_name = 'ACC_AABW_Ny128_Nlay3_tauM0.15_tauP0_tauF0_wDiaM1.5_wDiaP0_wDiaF0_Cd2.000e-03_rb0.000e+00_E1_doublewind';

%%% Load parameters   
local_home_dir = '/Volumes/Stewart-RAID1-A/UCLA/Projects/AWSIM_WindAABW/runs_varywind';
prod_dir = fullfile(local_home_dir,'products');
loadParams;
dirpath = fullfile(local_home_dir,run_name);
gtild = reshape(cumsum(gg),[Nlay 1 1]);
rho0 = 1000;

%%% Reference layer thicknesses for QG approx
HH = [1000 1500 1500]';

%%% Reference Coriolis parameter and beta for QG approx
f0 = mean(mean(0.25 * (2*Omega_z(1:Nx,1:Ny) + 2*Omega_z(1:Nx,2:Ny+1) + 2*Omega_z(2:Nx+1,1:Ny) + 2*Omega_z(2:Nx+1,2:Ny+1))));
beta = mean(mean(0.5 * (2*Omega_z(1:Nx,2:Ny+1)+2*Omega_z(2:Nx+1,2:Ny+1)-2*Omega_z(1:Nx,1:Ny)-2*Omega_z(2:Nx+1,1:Ny))/dy));

%%% Max time at which to load transports (in case we want to look at a
%%% subset)
% tend = 0.6*t1year;
% tend = 22*t1year;
tend = 160*t1year;

%%% Set true to use time-averaged momentum budget diagnostics. This will
%%% only work if those diagnostics are available!
use_avg_diags = true;






%%% Precompute iteration numbers and output times
iters = n0+1:1:n0+Nframes-1;
tt = startTime + (iters-n0)*dt_s;
iters(tt>tend) = [];
tt(tt>tend) = [];
Niters = length(iters);




%%% Time average baroclinic mode amplitudes
Mpsig_tavg = zeros(Nx,Ny,Nlay);
for n=1:Niters

  disp(n)
  disp(Niters)
  
  %%% Load model output fields
  [pi,hh,eta] = readPiHEta (dirpath,n,Nx,Ny,Nlay,hhb,use_avg_diags);

  %%% Compute QG streamfunction and meridional velocity
  [psig,etag] = calcQGStreamfunction (pi,eta,HH,gg,f0);    

  %%% Decompose streamfunction into baroclinic modes
  [Mpsig,EE,eps] = calcBaroclinicModes (psig,HH,gg);

  %%% Add to time average
  Mpsig_tavg = Mpsig_tavg + Mpsig/Niters;  

end





%%% Compute form stress at each time step
cntr = 0;  
formStressQG = zeros(Ny,Nlay,Niters);
formStressQG_fixedBC = zeros(Ny,Nlay,Niters);
formStressQG_fixedBT = zeros(Ny,Nlay,Niters);
for n=1:Niters

  disp(n)
  disp(Niters)
  
  %%% Load model output fields
  [pi,hh,eta] = readPiHEta (dirpath,n,Nx,Ny,Nlay,hhb,use_avg_diags);

  %%% Compute QG streamfunction and meridional velocity
  [psig,etag] = calcQGStreamfunction (pi,eta,HH,gg,f0);    
  vg = (psig(1:Nx,1:Ny,:)-psig([Nx 1:Nx-1],1:Ny,:)) / dx; %%% N.B. defined on C-grid u-point

  %%% Compute layer thickness fluxes (equivalent to form stresses)
  vghg = 0.5.*(vg(1:Nx,:,:)+vg([2:Nx 1],:,:)) .* hh;  
  vghg(:,2:Ny,:) = 0.5*(vghg(:,1:Ny-1,:)+vghg(:,2:Ny,:)); %%% Interpolate to C-grid v-points
  vghg(:,1,:) = 0;

  %%% Form stress from thickness fluxes
  % vghg_int = squeeze(sum(vghg,1)*dx*rho0*f0);
  % formStressQG(:,:,n) = cumsum(vghg_int,2);  

  %%% Form stress from streamfunctions
  psig_mid = 0.5*(psig(1:Nx,:,:)+psig([Nx 1:Nx-1],:,:));
  dpsig_dx = (psig(1:Nx,:,:)-psig([Nx 1:Nx-1],:,:)) / dx;
  detab_dx = (hhb(1:Nx,:,:)-hhb([Nx 1:Nx-1],:,:)) / dx;
  for k=1:Nlay-1    
    formStressQG(:,k,n) = -rho0*f0/gg(k+1)*sum(psig_mid(:,:,k).*(dpsig_dx(:,:,k+1)));
  end
  formStressQG(:,Nlay,n) = -rho0*sum(psig_mid(:,:,Nlay).*(detab_dx));

  %%% Decompose streamfunction into baroclinic modes
  [Mpsig,EE,eps] = calcBaroclinicModes (psig,HH,gg);

  %%% Construct a modified representation of the modes in which the
  %%% barotropic mode varies but the baroclinic modes are held fixed
  Mpsig_fixedBC = Mpsig;
  Mpsig_fixedBC(:,:,2:end) = Mpsig_tavg(:,:,2:end);
  %%% Construct a modified representation of the modes in which the
  %%% baroclinic modes vary but the barotropic mode is held fixed
  Mpsig_fixedBT = Mpsig;
  Mpsig_fixedBT(:,:,1) = Mpsig_tavg(:,:,1);

  %%% Reconstruct modified QG streamfunction
  psig_fixedBC = zeros(Nx,Ny,Nlay);
  psig_fixedBT = zeros(Nx,Ny,Nlay);
  for i=1:Nx
    for j=1:Ny
      psig_fixedBC(i,j,:) = reshape(EE * squeeze(Mpsig_fixedBC(i,j,:)),[1 1 Nlay]);
      psig_fixedBT(i,j,:) = reshape(EE * squeeze(Mpsig_fixedBT(i,j,:)),[1 1 Nlay]);
    end
  end

  %%% Form stress from streamfunctions with fixed baroclinic components
  psig_mid = 0.5*(psig_fixedBC(1:Nx,:,:)+psig_fixedBC([Nx 1:Nx-1],:,:));
  dpsig_dx = (psig_fixedBC(1:Nx,:,:)-psig_fixedBC([Nx 1:Nx-1],:,:)) / dx;
  detab_dx = (hhb(1:Nx,:,:)-hhb([Nx 1:Nx-1],:,:)) / dx;
  for k=1:Nlay-1    
    formStressQG_fixedBC(:,k,n) = -rho0*f0/gg(k+1)*sum(psig_mid(:,:,k).*(dpsig_dx(:,:,k+1)));
  end
  formStressQG_fixedBC(:,Nlay,n) = -rho0*sum(psig_mid(:,:,Nlay).*(detab_dx));

  %%% Form stress from streamfunctions with fixed barotropic components
  psig_mid = 0.5*(psig_fixedBT(1:Nx,:,:)+psig_fixedBT([Nx 1:Nx-1],:,:));
  dpsig_dx = (psig_fixedBT(1:Nx,:,:)-psig_fixedBT([Nx 1:Nx-1],:,:)) / dx;
  detab_dx = (hhb(1:Nx,:,:)-hhb([Nx 1:Nx-1],:,:)) / dx;
  for k=1:Nlay-1    
    formStressQG_fixedBT(:,k,n) = -rho0*f0/gg(k+1)*sum(psig_mid(:,:,k).*(dpsig_dx(:,:,k+1)));
  end
  formStressQG_fixedBT(:,Nlay,n) = -rho0*sum(psig_mid(:,:,Nlay).*(detab_dx));

end

%%% Write to .mat file
save(fullfile(prod_dir,[run_name,'_QGMomBalance.mat']), ...  
  'XX_h','YY_h','hhb','gg', ...
  'tt','formStressQG','formStressQG_fixedBC','formStressQG_fixedBT', ...
  '-v7.3');


formStressQG_cavg = squeeze(sum(formStressQG*dy,1)/Lx/Ly);
formStressQG_fixedBC_cavg = squeeze(sum(formStressQG_fixedBC*dy,1)/Lx/Ly);
formStressQG_fixedBT_cavg = squeeze(sum(formStressQG_fixedBT*dy,1)/Lx/Ly);

%%% Sample plot
figure(10);
plot(tt/t1year,formStressQG_cavg(1,:));
hold on;
plot(tt/t1year,formStressQG_cavg(2,:));
plot(tt/t1year,formStressQG_cavg(3,:));
hold off;
legend('IFS_u_p_p_e_r','IFS_l_o_w_e_r','TFS');
title('Momentum balance, channel-averaged');
xlabel('Time (years)')
ylabel('N/m^2');

%%% Sample plot
figure(11);
plot(tt/t1year,formStressQG_fixedBC_cavg(1,:));
hold on;
plot(tt/t1year,formStressQG_fixedBC_cavg(2,:));
plot(tt/t1year,formStressQG_fixedBC_cavg(3,:));
hold off;
legend('IFS_u_p_p_e_r','IFS_l_o_w_e_r','TFS');
title('Momentum balance, channel-averaged');
xlabel('Time (years)')
ylabel('N/m^2');

%%% Sample plot
figure(9);
plot(tt/t1year,formStressQG_fixedBT_cavg(1,:));
hold on;
plot(tt/t1year,formStressQG_fixedBT_cavg(2,:));
plot(tt/t1year,formStressQG_fixedBT_cavg(3,:));
hold off;
legend('IFS_u_p_p_e_r','IFS_l_o_w_e_r','TFS');
title('Momentum balance, channel-averaged');
xlabel('Time (years)')
ylabel('N/m^2');


figure(12);
scatter(formStressQG_cavg(1,:),formStressQG_fixedBC_cavg(1,:));
corr(formStressQG_cavg(1,:)',formStressQG_fixedBC_cavg(1,:)')

figure(13);
scatter(formStressQG_cavg(2,:),formStressQG_fixedBC_cavg(2,:));
corr(formStressQG_cavg(2,:)',formStressQG_fixedBC_cavg(2,:)')

figure(14);
scatter(formStressQG_cavg(3,:),formStressQG_fixedBC_cavg(3,:));
corr(formStressQG_cavg(3,:)',formStressQG_fixedBC_cavg(3,:)')

figure(15);
scatter(formStressQG_cavg(1,:),formStressQG_fixedBT_cavg(1,:));
corr(formStressQG_cavg(1,:)',formStressQG_fixedBT_cavg(1,:)')

figure(16);
scatter(formStressQG_cavg(2,:),formStressQG_fixedBT_cavg(2,:));
corr(formStressQG_cavg(2,:)',formStressQG_fixedBT_cavg(2,:)')

figure(17);
scatter(formStressQG_cavg(3,:),formStressQG_fixedBT_cavg(3,:));
corr(formStressQG_cavg(3,:)',formStressQG_fixedBT_cavg(3,:)')



%%%
%%% Convenience function to read in needed fields from a given output
%%% snapshot time.
%%%
function [pi,hh,eta] = readPiHEta (dirpath,n,Nx,Ny,Nlay,etab,use_avg_diags)

  constants;

  %%% To store data read from model output files
  eta = zeros(Nx,Ny,Nlay+1);  
  hh = zeros(Nx,Ny,Nlay);

  if (use_avg_diags)
      
    %%% Load time-averaged model state
    for k=1:Nlay
      data_file = fullfile(dirpath,[OUTN_H_AVG,num2str(k-1),'_n=',num2str(n),'.dat']);
      hh(:,:,k) = readOutputFile(data_file,Nx,Ny);
    end  
    data_file = fullfile(dirpath,[OUTN_PI_AVG,'_n=',num2str(n),'.dat']);
    pi = readOutputFile(data_file,Nx,Ny);
    
  else
    
    %%% Load instantaneous model state
    for k=1:Nlay
      data_file = fullfile(dirpath,[OUTN_H,num2str(k-1),'_n=',num2str(n),'.dat']);
      hh(:,:,k) = readOutputFile(data_file,Nx,Ny);
    end  
    data_file = fullfile(dirpath,[OUTN_PI,'_n=',num2str(n),'.dat']);
    pi = readOutputFile(data_file,Nx,Ny);

  end

  %%% Calculate layer surface heights and layer mid-depths
  eta(:,:,Nlay+1) = etab;      
  for k=Nlay:-1:1               
    eta(:,:,k) = eta(:,:,k+1) + hh(:,:,k);
  end

end
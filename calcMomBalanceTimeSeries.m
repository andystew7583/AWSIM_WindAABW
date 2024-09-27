%%%
%%% calcMomBalanceTimeSeries.m
%%%
%%% Computes time series of wind stress and form stresses.
%%%

%%% Run to load
% run_name = 'ACC_AABW_Ny128_Nlay3_tauM0.15_tauP0.075_tauF45.625_wDiaM0_wDiaP0_wDiaF0_Cd2.000e-03_rb0.000e+00_E7_diags';
% run_name = 'ACC_AABW_Ny128_Nlay3_tauM0.15_tauP0_tauF0_wDiaM1.5_wDiaP0_wDiaF0_Cd2.000e-03_rb0.000e+00_E9_spinup';
run_name = 'ACC_AABW_Ny128_Nlay3_tauM0.15_tauP0_tauF0_wDiaM0_wDiaP0_wDiaF0_Cd2.000e-03_rb0.000e+00_E9_diags';


%%% Load parameters   
local_home_dir = '/Volumes/Stewart-RAID1-A/UCLA/Projects/AWSIM_WindAABW/runs_varywind';
prod_dir = fullfile(local_home_dir,'products');
loadParams;
dirpath = fullfile(local_home_dir,run_name);
gtild = reshape(cumsum(gg),[Nlay 1 1]);
rho0 = 1000;

%%% Max time at which to load transports
% tend = 0.6*t1year;
tend = 22*t1year;

%%% Set true to use time-averaged momentum budget diagnostics. This will
%%% only work if those diagnostics are available!
use_avg_diags = false;







%%% Time series of wind forcing
[tauNrecs tauNrecs_found] = readparam(params_file,'tauNrecs','%lf'); 
[tauPeriod tauPeriod_found] = readparam(params_file,'tauPeriod','%lf'); 
tauFile = 'taux.dat';
tauInt = zeros(1,tauNrecs);
tauTimes = (0:1:tauNrecs-1)/tauNrecs*tauPeriod;
fid = fopen(fullfile(dirpath,tauFile),'r','b');
if (fid == -1)
  error(['Could not open ',tauFile]);
end
taux = zeros(tauNrecs,Nx,Ny);
for j=1:Ny
  for i=1:Nx
    taux(:,i,j) = fread(fid,[tauNrecs 1],'real*8','ieee-le');
  end
end
for n=1:tauNrecs
  tauInt(n) = sum(sum(squeeze(taux(n,:,:)*dx*dy)*rho0));
end
fclose(fid);

%%% Time series of "true" AABW export
[wDiaNrecs wDiaNrecs_found] = readparam(params_file,'wDiaNrecs','%lf'); 
[wDiaPeriod wDiaPeriod_found] = readparam(params_file,'wDiaPeriod','%lf'); 
wDiaFile = 'wDiaFile.dat';
wDiaInt = zeros(1,wDiaNrecs);
wDiaTimes = (0:1:wDiaNrecs-1)/wDiaNrecs*wDiaPeriod;
fid = fopen(fullfile(dirpath,wDiaFile),'r','b');
if (fid == -1)
  error(['Could not open ',wDiaFile]);
end
wDia = zeros(wDiaNrecs,Nlay+1,Nx,Ny);
for j=1:Ny
  for i=1:Nx
    wDia(:,:,i,j) = fread(fid,[wDiaNrecs Nlay+1],'real*8','ieee-le');
  end
end
for n=1:wDiaNrecs
  wDiaInt(n) = sum(sum(squeeze(wDia(n,Nlay,:,:))*dx*dy));
end
fclose(fid);



%%% Precompute iteration numbers and output times
iters = n0+1:1:n0+Nframes-1;
tt = startTime + (iters-n0)*dt_s;
iters(tt>tend) = [];
tt(tt>tend) = [];
Niters = length(iters);

%%% At each time iteration...
cntr = 0;  
formStress = zeros(Nlay,Niters);
MOC = zeros(Nlay,Niters);
surfStress = zeros(1,Niters);
formStress_mid = zeros(Nlay,Niters);
MOC_mid = zeros(Nlay,Niters);
surfStress_mid = zeros(1,Niters);
h = zeros(Nlay,Nx,Ny);
v = zeros(Nlay,Nx,Ny);
u = zeros(Nlay,Nx,Ny);
M = zeros(Nlay,Nx,Ny);
hdMdx = zeros(Nlay,Nx,Ny);
taux = zeros(Nlay,Nx,Ny);
hv = zeros(Nlay,Nx,Ny);
eta = zeros(Nlay+1,Nx,Ny);
pi = zeros(Nx,Ny);  
for n=1:Niters

  disp(n)
  disp(Niters)
  
  if (use_avg_diags)
      
    %%% Load pressure gradient tendency directly
    for k=1:Nlay
      data_file = fullfile(dirpath,[OUTN_UMOM_GRADM,num2str(k-1),'_n=',num2str(n),'.dat']);
      hdMdx(k,:,:) = readOutputFile(data_file,Nx,Ny);
      data_file = fullfile(dirpath,[OUTN_UMOM_WIND,num2str(k-1),'_n=',num2str(n),'.dat']);
      taux(k,:,:) = readOutputFile(data_file,Nx,Ny);
      data_file = fullfile(dirpath,[OUTN_HV_AVG,num2str(k-1),'_n=',num2str(n),'.dat']);
      hv(k,:,:) = readOutputFile(data_file,Nx,Ny);
    end
    
    surfStress(n) = sum(sum(sum(taux*dx*dy*rho0)));
    surfStress_mid(n) = sum(sum(taux(:,:,Ny/2)*dx*rho0));
    MOC(:,n) = mean(sum(hv*dx,2),3);
    MOC_mid(:,n) = sum(hv(:,:,Ny/2)*dx,2);
    
  else
    
    %%% Load instantaneous model state
    for k=1:Nlay
      data_file = fullfile(dirpath,[OUTN_H,num2str(k-1),'_n=',num2str(n),'.dat']);
      h(k,:,:) = readOutputFile(data_file,Nx,Ny);
    end  
    for k=1:Nlay
      data_file = fullfile(dirpath,[OUTN_U,num2str(k-1),'_n=',num2str(n),'.dat']);
      u(k,:,:) = readOutputFile(data_file,Nx,Ny);
    end
    for k=1:Nlay
      data_file = fullfile(dirpath,[OUTN_V,num2str(k-1),'_n=',num2str(n),'.dat']);
      v(k,:,:) = readOutputFile(data_file,Nx,Ny);
    end       
    data_file = fullfile(dirpath,[OUTN_PI,'_n=',num2str(n),'.dat']);
    pi = readOutputFile(data_file,Nx,Ny);

    %%% Calculate layer surface heights and layer mid-depths
    eta(Nlay+1,:,:) = hhb;      
    for k=Nlay:-1:1               
      eta(k,:,:) = eta(k+1,:,:) + h(k,:,:);
    end

    %%% Calculate Montgomery potential       
    M(1,:,:) = pi;    
    for k=2:Nlay
      M(k,:,:) = M(k-1,:,:) + gg(k)*eta(k,:,:);          
    end     
    %%% Salmon layer adds nothing to integrated mom balance
    % for k=1:Nlay
    %   M(k,:,:) = M(k,:,:) - gsum(k) .* h0.^4 ./ hh(k,:,:).^3 ./ 3;
    % end   

    %%% Tendency due to Montgomery potential gradient
    h_w = 0.5*(h(:,1:Nx,:)+h(:,[Nx 1:Nx-1],:));
    hdMdx = h_w.*(M(:,1:Nx,:)-M(:,[Nx 1:Nx-1],:))/dx;
    
  end
    
  %%% Form stress
  hdMdx_int = squeeze(sum(sum(hdMdx,2),3)*dx*dy*rho0);
  formStress(:,n) = -cumsum(hdMdx_int,1);
  hdMdx_int = squeeze(sum(hdMdx(:,:,Ny/2),2)*dx*rho0);
  formStress_mid(:,n) = -cumsum(hdMdx_int,1);

end

%%% For cases with steady forcing
if (~use_avg_diags && tauPeriod == 0)
  tauTimes = tt;
  tauInt = tauInt(1)*ones(size(tt));
end

%%% Write to .mat file
save(fullfile(prod_dir,[run_name,'_MomBalance.mat']), ...
  'wDiaTimes','wDiaInt', ...
  'tauTimes','tauInt', ...
  'XX_h','YY_h','hhb','gg', ...
  'tt','formStress','surfStress','MOC', ...
  'tt','formStress_mid','surfStress_mid','MOC_mid', ...
  '-v7.3');

%%% Sample plot
figure(1);
plot(tt/t1year,surfStress/Lx/Ly);
hold on;
plot(tt/t1year,formStress(1,:)/Lx/Ly);
plot(tt/t1year,formStress(2,:)/Lx/Ly);
plot(tt/t1year,formStress(3,:)/Lx/Ly);
hold off;
legend('Surface stress','IFS_u_p_p_e_r','IFS_l_o_w_e_r','TFS');
title('Momentum balance, channel-averaged');

f0 = 2*mean(Omega_z(:));
%%% Sample plot
figure(2);
plot(tt/t1year,surfStress/Ly/rho0/abs(f0)/1e6);
hold on;
plot(tt/t1year,MOC(1,:)/1e6);
plot(tt/t1year,MOC(2,:)/1e6);
plot(tt/t1year,MOC(3,:)/1e6);
hold off;
legend('T_E_k_m_a_n','T_1','T_2','T_3');
title('Overturning, channel-averaged');

%%% Sample plot
figure(3);
plot(tt/t1year,surfStress_mid/Lx);
hold on;
plot(tt/t1year,formStress_mid(1,:)/Lx);
plot(tt/t1year,formStress_mid(2,:)/Lx);
plot(tt/t1year,formStress_mid(3,:)/Lx);
hold off;
legend('Surface stress','IFS_u_p_p_e_r','IFS_l_o_w_e_r','TFS');
title('Momentum balance in channel center');


%%% Sample plot
figure(4);
plot(tt/t1year,surfStress_mid/rho0/abs(f0)/1e6);
hold on;
plot(tt/t1year,MOC_mid(1,:)/1e6);
plot(tt/t1year,MOC_mid(2,:)/1e6);
plot(tt/t1year,MOC_mid(3,:)/1e6);
hold off;
legend('T_E_k_m_a_n','T_1','T_2','T_3');
title('Overturning in channel center');
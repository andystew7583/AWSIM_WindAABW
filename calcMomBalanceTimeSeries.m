%%%
%%% calcMomBalanceTimeSeries.m
%%%
%%% Computes time series of wind stress and form stresses.
%%%

%%% Load constant parameters
constants;
Ny = 128;
Nlay = 3;
tau_mean = [0.15];
% tau_pert = 0.075;
tau_pert = 0;
% tau_freq = t1year * 2^0;
tau_freq = 0;
% AABW_mean = 1.5;
AABW_mean = 0;
AABW_pert = 0;
AABW_freq = 0;
quad_drag = 2e-3;
lin_drag = 0e-4;  
topog_width = 150;
topog_height = 1000;
rough_topog = false;
rough_topog_height = 50;
double_wind = true;
double_ridge = true;
n_E_batch = 1:10;

init = true;
for n_E = n_E_batch
  run_name = constructRunName (false,Ny,Nlay, ...
                                    tau_mean,tau_pert,tau_freq, ...
                                    AABW_mean,AABW_pert,AABW_freq, ...
                                    quad_drag,lin_drag,...
                                    topog_width,topog_height,rough_topog,rough_topog_height,...
                                    double_ridge,double_wind,n_E);
  
  % run_name = 'ACC_AABW_Ny128_Nlay3_tauM0.15_tauP0_tauF0_wDiaM1.5_wDiaP0_wDiaF0_Cd2.000e-03_rb0.000e+00_E1_doublewind';
  
  %%% Load parameters   
  local_home_dir = '/Volumes/Stewart-RAID1-A/UCLA/Projects/AWSIM_WindAABW/runs_varywind';
  prod_dir = fullfile(local_home_dir,'products');
  loadParams;
  dirpath = fullfile(local_home_dir,run_name);
  gtild = reshape(cumsum(gg),[Nlay 1 1]);
  rho0 = 1000;
  
  %%% Max time at which to load transports
  % tend = 0.6*t1year;
  % tend = 22*t1year;
  tend = 160*t1year;
  
  %%% Set true to use time-averaged momentum budget diagnostics. This will
  %%% only work if those diagnostics are available!
  use_avg_diags = true;
  
  
  
  
  
  
  
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

  %%% Initialize storage arrays
  if (init)
    formStress = zeros(Nlay,Niters,length(n_E_batch));
    formStress_mean_tot = zeros(Nlay,Niters,length(n_E_batch));
    formStress_mean_ridge = zeros(Nlay,Niters,length(n_E_batch));
    formStress_mean_bumps = zeros(Nlay,Niters,length(n_E_batch));
    TFS_tot = zeros(Nlay,Niters,length(n_E_batch));
    TFS_ridge = zeros(Nlay,Niters,length(n_E_batch));
    TFS_bumps = zeros(Nlay,Niters,length(n_E_batch));
    MOC = zeros(Nlay,Niters,length(n_E_batch));
    surfStress = zeros(1,Niters,length(n_E_batch));
    formStress_mid = zeros(Nlay,Niters,length(n_E_batch));
    MOC_mid = zeros(Nlay,Niters,length(n_E_batch));
    surfStress_mid = zeros(1,Niters,length(n_E_batch));
    Tacc = zeros(1,Niters,length(n_E_batch));
    Tacc_bc = zeros(1,Niters,length(n_E_batch));
    Tacc_bt = zeros(1,Niters,length(n_E_batch));
    init = false;
  end
  
  %%% At each time iteration...
  cntr = 0;   
 
  h = zeros(Nlay,Nx,Ny);
  v = zeros(Nlay,Nx,Ny);
  u = zeros(Nlay,Nx,Ny);
  M = zeros(Nlay,Nx,Ny);
  hdMdx = zeros(Nlay,Nx,Ny);
  taux = zeros(Nlay,Nx,Ny);
  hv = zeros(Nlay,Nx,Ny);
  hu = zeros(Nlay,Nx,Ny);
  u = zeros(Nlay,Nx,Ny);
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
        data_file = fullfile(dirpath,[OUTN_HU_AVG,num2str(k-1),'_n=',num2str(n),'.dat']);
        hu(k,:,:) = readOutputFile(data_file,Nx,Ny);
        data_file = fullfile(dirpath,[OUTN_U_AVG,num2str(k-1),'_n=',num2str(n),'.dat']);
        u(k,:,:) = readOutputFile(data_file,Nx,Ny);
        data_file = fullfile(dirpath,[OUTN_H_AVG,num2str(k-1),'_n=',num2str(n),'.dat']);
        h(k,:,:) = readOutputFile(data_file,Nx,Ny);
        data_file = fullfile(dirpath,[OUTN_M_AVG,num2str(k-1),'_n=',num2str(n),'.dat']);
        M(k,:,:) = readOutputFile(data_file,Nx,Ny);
      end
      
      surfStress(1,n,n_E) = sum(sum(sum(taux*dx*dy*rho0)));
      surfStress_mid(1,n,n_E) = sum(sum(sum(taux(:,:,Ny/4:3*Ny/4)*dy*dx*rho0)));
      MOC(:,n,n_E) = mean(sum(hv*dx,2),3);
      MOC_mid(:,n,n_E) = mean(sum(hv(:,:,Ny/4:3*Ny/4)*dx,2),3);
      
      Tacc(1,n,n_E) = sum(sum(hu(:,1,:)*dy,3),1);
      Tacc_bt_tmp = sum(squeeze(u(Nlay,:,:)).*(-hhb)*dy,2);
      Tacc_bt(1,n,n_E) = mean(Tacc_bt_tmp);
      Tacc_bc(1,n,n_E) = Tacc(1,n,n_E) - Tacc_bt(1,n,n_E);
      
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
    formStress(:,n,n_E) = -cumsum(hdMdx_int,1);
    hdMdx_int = squeeze(sum(sum(hdMdx(:,:,Ny/4:3*Ny/4),2),3)*dx*dy*rho0);
    formStress_mid(:,n,n_E) = -cumsum(hdMdx_int,1);

    %%% Form stress decomposition
    h_yavg = mean(h,3);
    M_yavg = mean(M,3);
    hw_yavg = 0.5*(h_yavg(:,1:Nx)+h_yavg(:,[Nx 1:Nx-1]));
    dMdx_yavg = (M_yavg(:,1:Nx)-M_yavg(:,[Nx 1:Nx-1])) / dx;    
    h_w = 0.5*(h(:,1:Nx,:)+h(:,[Nx 1:Nx-1],:));
    hdMdx_mean = h_w.*(M(:,1:Nx,:)-M(:,[Nx 1:Nx-1],:))/dx;    
    hdMdx_mean_int = -squeeze(sum(sum(hdMdx_mean,2),3)*dx*dy*rho0);
    hdMdx_mean_ridge = -repmat(hw_yavg.*dMdx_yavg,[1 1 Ny]);
    hdMdx_mean_ridge_int = squeeze(sum(sum(hdMdx_mean_ridge,2),3)*dx*dy*rho0);
    formStress_mean_tot(:,n,n_E) = -cumsum(hdMdx_mean_int,1);
    formStress_mean_ridge(:,n,n_E) = -cumsum(hdMdx_mean_ridge_int,1);
    formStress_mean_bumps(:,n,n_E) = formStress_mean_tot(:,n,n_E) - formStress_mean_ridge(:,n,n_E);

    %%% Latitudinal ridge
    Wb = 150*1000;
    Xb = 1000*m1km;
    if (double_ridge)
      Xb = Ly/2;
      Xb2 = 3*Ly/2;
    end
    Hb = 1000;
    H = 4000;
    etab_ridge = Hb*exp(-((XX_h-Xb)/Wb).^2);    
    if (double_ridge)
      etab_ridge = etab_ridge+Hb*exp(-((XX_h-Xb2)/Wb).^2);
    end
    etab_ridge = etab_ridge - H;
    etab_bumps = hhb - etab_ridge;
    
    %%% Load instantaneous model state    
    for k=1:Nlay
      data_file = fullfile(dirpath,[OUTN_M_AVG,num2str(k-1),'_n=',num2str(n),'.dat']);
      M(k,:,:) = readOutputFile(data_file,Nx,Ny);
    end       
    detab_dx = (hhb(1:Nx,:) - hhb([Nx 1:Nx-1],:)) / dx;
    detab_ridge_dx = (etab_ridge(1:Nx,:) - etab_ridge([Nx 1:Nx-1],:)) / dx;
    detab_bumps_dx = (etab_bumps(1:Nx,:) - etab_bumps([Nx 1:Nx-1],:)) / dx;
    TFS_tot(1,n,n_E) = sum(sum(squeeze(0.5*(M(Nlay,1:Nx,:)+M(Nlay,[Nx 1:Nx-1],:))).*detab_dx*dx*dy*rho0,1),2);
    TFS_ridge(1,n,n_E) = sum(sum(squeeze(0.5*(M(Nlay,1:Nx,:)+M(Nlay,[Nx 1:Nx-1],:))).*detab_ridge_dx*dx*dy*rho0,1),2);
    TFS_bumps(1,n,n_E) = sum(sum(squeeze(0.5*(M(Nlay,1:Nx,:)+M(Nlay,[Nx 1:Nx-1],:))).*detab_bumps_dx*dx*dy*rho0,1),2);
 
  end
  
  %%% For cases with steady forcing
  if (~use_avg_diags && tauPeriod == 0)
    tauTimes = tt;
    tauInt = tauInt(1)*ones(size(tt));
  end

end

%%% Write to .mat file
save(fullfile(prod_dir,[run_name,'_MomBalance.mat']), ...
  'n_E_batch',...
  'wDiaTimes','wDiaInt', ...
  'tauTimes','tauInt', ...
  'XX_h','YY_h','hhb','gg', ...
  'tt','formStress','surfStress','MOC', ...
  'formStress_mid','surfStress_mid','MOC_mid', ...
  'Tacc','Tacc_bt','Tacc_bc', ...
  'formStress_mean_tot','formStress_mean_bumps','formStress_mean_ridge', ...
  'TFS_tot','TFS_bumps','TFS_ridge', ...
  '-v7.3');

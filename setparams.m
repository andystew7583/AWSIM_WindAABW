%%%
%%% setparams.m
%%%
%%% Sets parameters for AWSIM. This file configures an ACC-like
%%% channel.
%%%
%%% local_home_dir  Directory to hold simulation folder
%%% run_name        Name of simulation
%%% is_spinup       Set true if this is a spinup simulation, false if we're
%%%                 averaging diagnostics online
%%% grid_size       Number of meridional grid points
%%% num_layers      Number of isopycnal layers to use (should be 2 or 3)
%%% tau_mean        Time-mean wind stress in N/m^2
%%% tau_pert        Amplitude of wind stress fluctuations in N/m^2
%%% tau_freq        Period of wind stress fluctuations in s
%%% AABW_mean       Time-mean AABW formation in Sv
%%% AABW_pert       Amplitude of AABW formation fluctuations in Sv
%%% AABW_freq       Period of AABW formation fluctuations in s
%%% quad_drag       Quadratic drag coefficient (dimensionless)
%%% lin_drag        Linear drag coefficient (m/s)
%%% topog_width     Zonal width of topographic bump in km
%%% topog_height    Height of topographic bump in m
%%% restart_index   Index of the output file used to restart the run
%%% end_time        Simulation end time (s)
%%% 
function setparams (local_home_dir,run_name, ...
  is_spinup,grid_size,num_layers, ...
  tau_mean,tau_pert,tau_freq, ...
  AABW_mean,AABW_pert,AABW_freq, ...
  quad_drag, lin_drag, topog_width, ...
  topog_height,restart_index, end_time)

  %%% Set true to run with random forcing, rather than periodic forcing.
  random_forcing = false;

  %%% Load constant parameters 
  constants;
  
  %%% Run directory
  run_name = strtrim(run_name); 
  local_home_dir = strtrim(local_home_dir); 
  local_run_dir = fullfile(local_home_dir,run_name);
  mkdir(local_run_dir);
  pfname = fullfile(local_run_dir,[run_name,'_in']);   
  model_code_dir = fullfile('../AWSIM/',model_code_dir_name);
  
  %%% Cluster config
  %%% NOTE: You will need to edit matlab_common/createRunScript to add a
  %%% configuration for your cluster!
  use_cluster = true;
  use_intel = false;
  use_pbs = use_cluster;
  uname = 'astewart';
  cluster_addr = 'caolila.atmos.ucla.edu';
  cluster_home_dir = '/data2/astewart/AWSIM_WindAABW/runs';
  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%% PARAMETERS %%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  %%% To store parameters
  paramTypes;
  PARAMS = {};
  
  %%% Select vertical resolution
  Nlay = num_layers;
  
  %%% Physical parameters
  rho0 = 1000;                  %%% Reference density 
  f0 = -1e-4;                   %%% Coriolis parameter
  beta = 1.5e-11;               %%% Coriolis parameter gradient
  beta_t = 0e-11;               %%% Topographic beta
  Ly = 1600*m1km;               %%% Domain length.   
  if (Nlay == 2)                %%% Reduced gravities at layer interfaces
    geff = [g .5e-2];
  end
  if (Nlay == 3)
    geff = [g .5e-2 .2e-2];        
  end
  h0 = 0;                       %%% Salmon layer thickness
%   Xb = 1000*m1km;               %%% Zonal position of topography  
  Xb = 1600*m1km;               %%% Zonal position of topography  
  Wb = topog_width*m1km;                %%% Zonal width of topography  
  Hb = topog_height;                    %%% Height of topography  
  H = 4000;                     %%% Ocean depth  
  if (Nlay == 2)
    H0 = [1500 2500];        %%% Initial layer thicknesses - used for wave speed calculation  
  end
  if (Nlay == 3)
    H0 = [1000 1000 2000];        %%% Initial layer thicknesses - used for wave speed calculation  
  end
  E0 = 0.01;                    %%% Initial EKE density
  rb = lin_drag;                %%% Linear bottom drag (reference value 2e-4 m/s)
  Cd = quad_drag;               %%% Quadratic bottom drag (reference value 2e-3)
  deta2 = 1000;                 %%% Initial isopycnal depth change across the channel
%   eta_north = [-1350 -2350];    %%% Relaxation layer depths at northern boundary
  if (Nlay == 2)
    eta_north = [-1750];    %%% Relaxation layer depths at northern boundary
    eta_south = [-1150];     %%% Relaxation layer deptsh at southern boundary
  end
  if (Nlay == 3)
    eta_north = [-1250 -2750];    %%% Relaxation layer depths at northern boundary
%     eta_north = [-1000 -2750];    %%% Relaxation layer depths at northern boundary
    eta_south = [-650 -1650];     %%% Relaxation layer deptsh at southern boundary
  end
  tRelax = 7*t1day;            %%% Relaxation time scale
  Lrelax = 100*m1km;            %%% Width of buoyancy forcing zone    
  
  %%% Time-varying wind stress parameters
  tau0 = tau_mean;                  %%% Wind stress maximum  (reference value 0.15)
  dtau0 = tau_pert;                 %%% Amplitude of wind stress fluctuations (reference value 0.05)
  tauPeriod = tau_freq;           %%% Wind stress period (reference value 500 years)
  tauNrecs = 24;            %%% Number of temporal wind stress records to write  
%   tau0 = 0.15;                  %%% Wind stress maximum  
%   dtau0 = 0;                 %%% Amplitude of wind stress fluctuations
%   tauPeriod = 500*t1year;       %%% Wind stress period 
%   tauNrecs = 1;            %%% Number of temporal wind stress records to write  
  
  %%% Time-varying buoyancy forcing parameters
  Psi0_upper = 0e6;                   %%% Imposed AAIW formation and export (reference value 1.5 Sv)
  dPsi0_upper = 0.5e6;                  %%% Amplitude of AAIW export fluctuations (reference value 0.5 Sv)
  Psi0_lower = AABW_mean*1e6;                   %%% Imposed AABW formation and export (reference value 1.5 Sv)
  dPsi0_lower = AABW_pert*1e6;                  %%% Amplitude of AABW export fluctuations (reference value 0.5 Sv)
  wDiaPeriod = AABW_freq;               %%% Diapycnal velocity period (reference value (reference value 500 years)                
  wDiaNrecs = 24;                %%% Number of temporal diapycnal velocity records to write  
  

  if ((AABW_freq <= 0) && (tau_freq <= 0))
    
    %%% For runs with totally steady forcing, just use fixed output
    %%% frequencies
    if (is_spinup)
      tmax = end_time; 
      savefreq = 1*t1year;   
      savefreqEZ = 1*t1day;  
      savefreqAvg = -t1year;
      savefreqUMom = -t1year;
      savefreqVMom= -1;
      savefreqThic = -1;     
    else
      tmax = end_time; 
%       savefreq = 5*t1day;         
      savefreq = 1*t1year;   
      savefreqEZ = 1*t1day;  
      savefreqAvg = 5*t1year;
      savefreqUMom = 5*t1year;
      savefreqVMom= 5*t1year;
      savefreqThic = -1; 
    end
    
  else
    
    %%% Frequency of diagnostic output chosen to resolve the smallest forcing
    %%% period
    diag_freq = -1;
    if ((AABW_freq > 0) && (tau_freq > 0))
      diag_freq = min([AABW_freq tau_freq]);
    else
      if (AABW_freq > 0)
        diag_freq = AABW_freq;
      else
        if (tau_freq > 0)
          diag_freq = tau_freq;
        end
      end
    end

    %%% Temporal parameters  
    if (is_spinup)
      %%% In spinup simulations we don't need much output - just enough to
      %%% check the model state looks ok and check equilibration
      tmax = end_time; 
      savefreq = 1*t1year;   
      savefreqEZ = 1*t1day;  
      savefreqAvg = -t1year;
      savefreqUMom = -t1year;
      savefreqVMom= -1;
      savefreqThic = -1;     
    else
      %%% For diagnostic runs we need lots of output, sufficient to resolve
      %%% the frequency of the forcing
      tmax = 10*diag_freq; %%% Long enough to sample many forcing cycles
      savefreq = diag_freq  / 12; %%% Enough subdivisions to fully resolve the forcing cycle 
      savefreqEZ = min([1*t1day, diag_freq/12]);  
      savefreqAvg = diag_freq/12;
      savefreqUMom = diag_freq/12;
      savefreqVMom= -1;
      savefreqThic = -1;  
    end
    
  end  
  
  %%% Always restart from a previous run
  startIdx = restart_index;
  restart = startIdx >= 0;
  
      
  %%% Rigid lid-related parameters
  useRL = 1; %%% Set to 1 to use rigid lid, or 0 not to  
  use_MG = 1;
  tol = 1e-7;
  SOR_rp_max = 2.0;
  SOR_rp_min = 1.7;
  SOR_opt_freq = 1000; 
  
  %%% Grids  
  Ny = grid_size;
  Nx = 2*Ny;
  d = Ly/Ny;  
  Lx = Nx*d;
  xx_q = 0:d:Lx;
  yy_q = 0:d:Ly;  
  xx_h = d/2:d:Lx-d/2;
  yy_h = d/2:d:Ly-d/2;    
  [YY_q,XX_q] = meshgrid(yy_q,xx_q);
  [YY_h,XX_h] = meshgrid(yy_h,xx_h);
  [YY_u,XX_u] = meshgrid(yy_h,xx_q(1:Nx));
  [YY_v,XX_v] = meshgrid(yy_q(1:Ny),xx_h);    
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% BOTTOM TOPOGRAPHY %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
   
  %%% Latitudinal ridge
  etab = Hb*exp(-((XX_h-Xb)/Wb).^2);
  etab = etab - H;
  
  %%% Latitudinal bottom slope 
  etab = etab + (beta_t*H/f0)*(YY_h-(Ly/2));
  
  %%% Plot topography
  figure(10);
  surf(XX_h,YY_h,etab);
  shading interp;
  colorbar;
%   plot(xx_h,etab(:,1));
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% OTHER PARAMETERS %%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
   
  %%% Set timestep based on full gravity wave speed, assuming a that the
  %%% absolute velocity never exceeds the gravity wave speed. Then correct
  %%% dt so that tmax is an integer number of time steps.  
  c = calcWaveSpeed(H0,geff,useRL);
  Umax = 2;  
  dt = 0.25*d/(c+Umax);  
  Nt = ceil(tmax/dt)+1;       
  
  %%% Set viscosities; 
  A2 = 0;
%   A4 = 0.01*d^3*Umax;
%   tau4 = 100*t1day;
%   A4 = d^4/tau4; 
  A4 = 0;
  A4smag = 4;
         
  %%% Define parameters 
  PARAMS = addParameter(PARAMS,'Nlay',Nlay,PARM_INT);
  PARAMS = addParameter(PARAMS,'Nx',Nx,PARM_INT);
  PARAMS = addParameter(PARAMS,'Ny',Ny,PARM_INT);
  PARAMS = addParameter(PARAMS,'dt',dt,PARM_REALF);
  PARAMS = addParameter(PARAMS,'savefrequency',savefreq,PARM_REALF);  
  PARAMS = addParameter(PARAMS,'savefreqEZ',savefreqEZ,PARM_REALF);
  PARAMS = addParameter(PARAMS,'savefreqAvg',savefreqAvg,PARM_REALF);
  PARAMS = addParameter(PARAMS,'savefreqUMom',savefreqUMom,PARM_REALF);
  PARAMS = addParameter(PARAMS,'savefreqVMom',savefreqVMom,PARM_REALF);
  PARAMS = addParameter(PARAMS,'savefreqThic',savefreqThic,PARM_REALF);  
  PARAMS = addParameter(PARAMS,'tmax',tmax,PARM_REALF);
  PARAMS = addParameter(PARAMS,'restart',restart,PARM_INT);
  PARAMS = addParameter(PARAMS,'startIdx',startIdx,PARM_INT);
  PARAMS = addParameter(PARAMS,'Ly',Ly,PARM_REALF);
  PARAMS = addParameter(PARAMS,'h0',h0,PARM_REALF);
  PARAMS = addParameter(PARAMS,'A2',A2,PARM_REALF);
  PARAMS = addParameter(PARAMS,'A4',A4,PARM_REALF);  
  PARAMS = addParameter(PARAMS,'A4smag',A4smag,PARM_REALF);
  PARAMS = addParameter(PARAMS,'tauPeriod',tauPeriod,PARM_REALF);
  PARAMS = addParameter(PARAMS,'tauNrecs',tauNrecs,PARM_INT);
  PARAMS = addParameter(PARAMS,'wDiaPeriod',wDiaPeriod,PARM_REALF);
  PARAMS = addParameter(PARAMS,'wDiaNrecs',wDiaNrecs,PARM_INT);
  PARAMS = addParameter(PARAMS,'linDragCoeff',rb,PARM_REALF);
  PARAMS = addParameter(PARAMS,'quadDragCoeff',Cd,PARM_REALF);
  PARAMS = addParameter(PARAMS,'use_MG',use_MG,PARM_INT);
  PARAMS = addParameter(PARAMS,'useRL',useRL,PARM_INT);
  PARAMS = addParameter(PARAMS,'useWallEW',0,PARM_INT);
  PARAMS = addParameter(PARAMS,'useWallNS',1,PARM_INT);
  PARAMS = addParameter(PARAMS,'tol',tol,PARM_REALE);
  PARAMS = addParameter(PARAMS,'SOR_rp_max',SOR_rp_max,PARM_REALF);
  PARAMS = addParameter(PARAMS,'SOR_rp_min',SOR_rp_min,PARM_REALF);   
  PARAMS = addParameter(PARAMS,'SOR_opt_freq',SOR_opt_freq,PARM_INT);   
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% BACKGROUND ROTATION %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  
  %%% Background rotation rate
  Omegaz = 0.5* (f0*ones(Nx+1,Ny+1) + beta*(YY_q-(Ly/2)));
  


  

  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%% WIND STRESS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%
  
  %%% Time-varying wind stress
  taux = zeros(tauNrecs,Nx,Ny);
  deltaT_tau = tauPeriod/tauNrecs;
  tt_tau = [0:1:(tauNrecs-1)]*deltaT_tau;
  
  %%% Create wind stress time series
  if (random_forcing)
  
    %%% Random fluctuations
    freq = [0:1:ceil(tauNrecs-1)/2 -floor(tauNrecs/2):1:-1]/tauPeriod;
    phase = 2*pi*rand(1,tauNrecs);
    taufft = exp(1i*phase).*abs(freq).^(-.5);
    taufft(1) = 0;
    tau_amp = real(ifft(taufft));
    tau_amp = tau_amp * dtau0/std(tau_amp) + tau0;

  else
    
    %%% Periodic oscillation 
    if (tauPeriod > 0)        
      tau_amp = tau0 + dtau0*sin(2*pi*tt_tau/tauPeriod);                               
    else
      tau_amp = tau0*ones(size(tt_tau));
    end 
  
  end
  
  %%% Create wind stress input matrix
  for n=1:tauNrecs
    taux(n,:,:) = tau_amp(n) * sin(pi*YY_u/Ly).^2 / rho0;
  end
  
  %%% Visualize forcing  
  figure(40);
  plot(yy_h/1000,squeeze(mean(taux(:,1,:),1)));   
  
  figure(42)
  plot(tt_tau/t1year,tau_amp);
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% INITIAL CONDITIONS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    %%% Set sea surface height
  Rd = c/abs(f0)
  lambdaK = 4*Rd;  
  eta1 = (f0/g)*genRandIC(lambdaK,E0,Nx,Ny,Ly);
  eta1 = (1 - exp(-(YY_h./(lambdaK)).^2)) .* (1 - exp(-((YY_h-Ly)./(lambdaK)).^2)) .* eta1;  
  eta1 = eta1 + (geff(2)/g)*deta2*(YY_h-Ly/2)/Ly;  
  eta2 = -H0(1) - (g/geff(2))*eta1;
  h1 = -eta2;
  if (Nlay == 2)
    h2 = eta2-etab;
  end
  if (Nlay == 3)
    eta3 = -H0(1)-H0(2) - deta2*(YY_h-Ly/2)/Ly;  
    h2 = eta2-eta3;
    h3 = eta3-etab;
  end
 
  %%% Plot layer interfaces
  figure(1);
  contourf(XX_h,YY_h,eta1);
  colorbar;  
  figure(2);
  contourf(XX_h,YY_h,eta2);
  colorbar;
  if (Nlay == 3)
    figure(3);
    contourf(XX_h,YY_h,eta3);
    colorbar;
  end
  
  %%% Set geostrophic zonal velocity
  eta1_mat1 = circshift(eta1, [0,-1]);  
  eta1_mat2 = circshift(eta1, [1,-1]);
  eta1_mat3 = circshift(eta1, [1,0]);
  eta1_mat4 = circshift(eta1, [1,1]);
  eta1_mat5 = circshift(eta1, [0,1]); 
  u1 = -(g/f0*(1/d)) .* ((1/4).*(eta1_mat1+eta1_mat2+eta1_mat3+eta1) - (1/4).*(eta1+eta1_mat3+eta1_mat4+eta1_mat5));
  u2 = 0*u1;
  u3 = 0*u1;
  
  %%% Quick fix for near-boundary points
  u1(:,1) = u1(:,2);  
  u2(:,1) = u2(:,2);
  u3(:,1) = u3(:,2);
  u1(:,Ny) = u1(:,Ny-1);
  u2(:,Ny) = u2(:,Ny-1);
  u3(:,Ny) = u3(:,Ny-1);
  
  %%% Plot along-slope velocity
  figure(4);
  contourf(XX_u,YY_u,u1);
  colorbar;
  
  %%% Set geostrophic meridional velocity
  eta1_mat6 = circshift(eta1, [-1,1]);
  eta1_mat7 = circshift(eta1, [-1,0]);  
  v1 = -(g/f0*(1/d)) .* ((1/4).*(eta1_mat3+eta1_mat4+eta1_mat5+eta1) - (1/4).*(eta1+eta1_mat5+eta1_mat6+eta1_mat7));
  v2 = 0*v1;
  v3 = 0*v1;
  
  %%% Apply initial velocity profiles in each layer
  uu = zeros(Nlay,Nx,Ny);
  uu(1,:,:) = u1;
  uu(2,:,:) = u2;  
  vv = zeros(Nlay,Nx,Ny);
  vv(1,:,:) = v1;
  vv(2,:,:) = v2;  
  hh = zeros(Nlay,Nx,Ny);
  hh(1,:,:) = h1;
  hh(2,:,:) = h2;
  if (Nlay == 3)
    uu(3,:,:) = u3;
    vv(3,:,:) = v3;
    hh(3,:,:) = h3; 
  end
  
  %%% Plot initial layer thickness
  figure(8);
  plot(xx_h,squeeze(hh(1,:,round(Ny/2))));
  
  
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% RELAXATION SETUP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%% Relaxation targets for layer thicknesses
  h1Relax = 0*h1;
  h2Relax = 0*h2;   
  h1Relax(YY_h<Ly/2) = - eta_south(1);
  h1Relax(YY_h>=Ly/2) = - eta_north(1);  
  if (Nlay == 2)
    h2Relax(YY_h<Ly/2) = eta_south(1) - etab(YY_h<Ly/2);
    h2Relax(YY_h>=Ly/2) = eta_north(1) - etab(YY_h>=Ly/2);      
  end
  if (Nlay == 3)
    h3Relax = 0*h3; 
    h2Relax(YY_h<Ly/2) = eta_south(1) - eta_south(2);
    h2Relax(YY_h>=Ly/2) = eta_north(1) - eta_north(2);  
    h3Relax(YY_h<Ly/2) = eta_south(2) - etab(YY_h<Ly/2);
    h3Relax(YY_h>=Ly/2) = eta_north(2) - etab(YY_h>=Ly/2);
  end  
  %%% Relaxation time scale
  hTime = -ones(Nx,Ny);  
%   hTime(YY_h<Lrelax) = tRelax ./ (1-YY_h(YY_h<Lrelax)/Lrelax);
  hTime(YY_h>Ly-Lrelax) = tRelax ./ (1-(Ly-YY_h(YY_h>Ly-Lrelax))/Lrelax);
  
  figure(99)
  plot(YY_h(1,:),h1Relax(1,:));
  hold on;
  plot(YY_h(1,:),h2Relax(1,:));
  hold off;
  
  figure(100);
  plot(YY_h(1,:),hTime(1,:));
  
  %%% Create input matrices
  hRelax = zeros(Nlay,Nx,Ny);
  hRelax(1,:,:) = h1Relax;
  hRelax(2,:,:) = h2Relax;
  if (Nlay == 3)    
    hRelax(3,:,:) = h3Relax;
  end
  
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% DIAPYCNAL VELOCITY %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  wDia = zeros(wDiaNrecs,Nlay+1,Nx,Ny);
  deltaT_wDia = wDiaPeriod/wDiaNrecs;
  southIdx = find(yy_h<Lrelax);
  northIdx = find(yy_h>Ly-Lrelax);
  Asouth = length(southIdx)*Nx*d^2;
  Anorth = length(northIdx)*Nx*d^2; 
  tt_wDia = [0:1:(wDiaNrecs-1)]*deltaT_wDia;
  
  %%% Create diapycnal velocity time series
  if (random_forcing)
    
    %%% Random fluctuations in AABW formation and AAIW formation
    freq = [0:1:ceil(wDiaNrecs-1)/2 -floor(wDiaNrecs/2):1:-1]/wDiaPeriod;
    phase_upper = 2*pi*rand(1,wDiaNrecs);
    phase_lower = 2*pi*rand(1,wDiaNrecs);
    Psifft_upper = exp(1i*phase_upper).*abs(freq).^(-.5);
    Psifft_upper(1) = 0;
    Psifft_lower = exp(1i*phase_lower).*abs(freq).^(-.5);
    Psifft_lower(1) = 0;
    Psi_upper = real(ifft(Psifft_upper));
    Psi_upper = Psi_upper * dPsi0_upper/std(Psi_upper) + Psi0_upper;
    Psi_lower = real(ifft(Psifft_lower));
    Psi_lower = Psi_lower * dPsi0_lower/std(Psi_lower) + Psi0_lower;
    
  else
    
    %%% Periodic oscillation in AABW formation 
    if (wDiaPeriod > 0)        
      Psi_lower = Psi0_lower + dPsi0_lower*sin(2*pi*tt_wDia/wDiaPeriod);                               
      Psi_upper = Psi0_upper + dPsi0_upper*sin(2*pi*tt_wDia/wDiaPeriod);                               
    else
      Psi_lower = Psi0_lower*ones(size(tt_wDia));
      Psi_upper = Psi0_upper*ones(size(tt_wDia));
    end         
    
  end
    
  %%% Create diapycnal velocity matrix
  for n=1:wDiaNrecs
    wDia(n,Nlay,:,southIdx) = - Psi_lower(n)/Asouth; %%% AABW formation all occurs between deepest two layers
    if (Nlay > 2)
      wDia(n,Nlay-1,:,southIdx) = Psi_upper(n)/Asouth; %%% AAIW formation all occurs between lightest two layers
    end
  end
  
  %%% Visualizations
  figure(30);
  plot(yy_h/1000,squeeze(wDia(1,:,1,:)));
  legend('1','2','3','4');
%   
%   figure(31);
%   loglog(freq,abs(Psifft_upper));
%   hold on;
%   loglog(freq,abs(Psifft_lower));
%   hold off;
%   legend('upper cell','lower cell');

  figure(32);
  plot(tt_wDia/t1year,Psi_upper);  
  hold on;
  plot(tt_wDia/t1year,smooth(Psi_upper,24));
  plot(tt_wDia/t1year,Psi_lower);  
  plot(tt_wDia/t1year,smooth(Psi_lower,24));
  hold off;
  legend('Upper cell','Upper cell (2 year smoothing)','Lower cell','Lower cell (2 year smoothing)');
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% CREATE INPUT FILES %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
  
  %%% Initial h 
  hInitFile = 'hInit.dat';
  writeDataFile(fullfile(local_run_dir,hInitFile),hh);
  PARAMS = addParameter(PARAMS,'hInitFile',hInitFile,PARM_STR);  
  
  %%% Initial u  
  uInitFile = 'uInit.dat';
  writeDataFile(fullfile(local_run_dir,uInitFile),uu);
  PARAMS = addParameter(PARAMS,'uInitFile',uInitFile,PARM_STR); 
  
  %%% Initial v
  vInitFile = 'vInit.dat';
  writeDataFile(fullfile(local_run_dir,vInitFile),vv);
  PARAMS = addParameter(PARAMS,'vInitFile',vInitFile,PARM_STR);  
  
  %%% Bathymetry
  etabFile = 'etab.dat';          
  writeDataFile(fullfile(local_run_dir,etabFile),etab);
  PARAMS = addParameter(PARAMS,'hbFile',etabFile,PARM_STR);
  
  %%% Background rotation
  OmegazFile = 'Omegaz.dat';          
  writeDataFile(fullfile(local_run_dir,OmegazFile),Omegaz);
  PARAMS = addParameter(PARAMS,'OmegazFile',OmegazFile,PARM_STR);
  
  %%% Reduced gravity
  gFile = 'gg.dat';          
  writeDataFile(fullfile(local_run_dir,gFile),geff);
  PARAMS = addParameter(PARAMS,'gFile',gFile,PARM_STR);

  %%% Background rotation
  tauxFile = 'taux.dat';          
  writeDataFile(fullfile(local_run_dir,tauxFile),taux);
  PARAMS = addParameter(PARAMS,'tauxFile',tauxFile,PARM_STR);
  
  %%% Relaxation values for h 
  hRelaxFile = 'hRelax.dat';
  writeDataFile(fullfile(local_run_dir,hRelaxFile),hRelax);
  PARAMS = addParameter(PARAMS,'hRelaxFile',hRelaxFile,PARM_STR);  
  
  %%% Relaxation timescale for h  
  hTimeFile = 'hTime.dat';
  writeDataFile(fullfile(local_run_dir,hTimeFile),hTime);
  PARAMS = addParameter(PARAMS,'hTimeFile',hTimeFile,PARM_STR);  
  
  %%% Imposed diapycnal velocity
  wDiaFile = 'wDiaFile.dat';
  writeDataFile(fullfile(local_run_dir,wDiaFile),wDia);
  PARAMS = addParameter(PARAMS,'wDiaFile',wDiaFile,PARM_STR); 

  %%% Create the input parameter file
  writeParamFile(pfname,PARAMS);      

  %%% Create a run script
  createRunScript (  local_home_dir, ...
                     run_name, ...
                     model_code_dir, ...
                     exec_name, ...
                     use_intel, ...
                     use_pbs, ... 
                     use_cluster, ...
                     uname, ...
                     cluster_addr, ...
                     cluster_home_dir);

end
%%%
%%% plotDragSensitivity.m
%%%
%%% Plots sensitivity of ACC transport to wind for different bottom drag coefficients.
%%%

%%% Load static definitions
constants;

%%% Directory to store runs
local_home_dir = '/Volumes/Kilchoman/UCLA/Projects/AWSIM_WindAABW/runs';

%%% Spinup simulations are long and produce no diagnostic, diagnostic
%%% simulations output high-frequency diagnostics to resolve the forcing
%%% period
is_spinup = false;

%%% Grid resolution 
Ny = 128;
Nlay = 2;

%%% Averaging period
% tmin = 170.5*t1year;
% tmax = 200.5*t1year;
tmin = 0.5*t1year;
tmax = 30.5*t1year;

%%% Parameters defining the batch of runs to plot
tau_mean = 0.1;
tau_pert = 0;
tau_freq = 0;
% AABW_mean = [-1.5 -.75 0 .75 1.5];
% AABW_mean = [-1.5 0 1.5];
AABW_mean = 0;
AABW_pert = 0;
AABW_freq = 0;
quad_drag = [.5e-3:.5e-3:4e-3];
lin_drag = 0e-4;
topog_width = 150;
topog_height = 1000;
% quad_drag = 0e-3;
% lin_drag = 2e-4;
N_Cd = length(quad_drag);

%%% Loop over runs and compute transport
Ttot = zeros(1,N_Cd);
Tbt = zeros(1,N_Cd);
Tbc = zeros(1,N_Cd);
kap = zeros(1,N_Cd);
r_kap = zeros(1,N_Cd);
for n_Cd = 1:N_Cd
    
  [n_Cd]
  %%% Simulation name
  run_name = constructRunName (is_spinup,Ny,Nlay, ...
                          tau_mean,tau_pert,tau_freq, ...
                          AABW_mean,AABW_pert,AABW_freq, ...
                          quad_drag(n_Cd),lin_drag,topog_width,topog_height);
  loadParams;

  %%% Read time-mean zonal flux   
  hu_tavg = do_avg(dirpath,OUTN_HU_AVG,Nx,Ny,Nlay,n0_avg,N_avg,dt_avg,tmin,tmax,startTime);
  u_tavg = do_avg(dirpath,OUTN_U_AVG,Nx,Ny,Nlay,n0_avg,N_avg,dt_avg,tmin,tmax,startTime);

  %%% Compute transports
  Ttot(n_Cd) = mean(sum(sum(hu_tavg,3),2),1)*dy;
  Tbt(n_Cd) = mean(sum(u_tavg(:,:,end).*(-hhb).*dy,2),1);
  Tbc(n_Cd) = Ttot(n_Cd) - Tbt(n_Cd);
  
  %%% Compute diffusivity
  [kap_bulk,nu_bulk,r_kap_bulk,r_nu_bulk,EKE_zavg] = calcBulkEddyViscDiff(local_home_dir,run_name);
  kap(n_Cd) = kap_bulk;
  r_kap(n_Cd) = r_kap_bulk;
   
end

%%% Make figure
figure(1);
clf;
plot(quad_drag,Ttot/1e6,'o-');
xlabel('Quadratic drag coefficient');
ylabel('Total transport (Sv)');
grid on;
print('-dpng',fullfile('figures',['TotalTransportSensitivity_drag_',datestr(datenum(datetime('now')),'ddmmmyyyy'),'.png']));



%%% Make figure
figure(2);
clf;
plot(quad_drag,Tbt/1e6,'o-');
xlabel('Quadratic drag coefficient');
ylabel('Barotropic transport (Sv)');
grid on;
print('-dpng',fullfile('figures',['BTTransportSensitivity_drag_',datestr(datenum(datetime('now')),'ddmmmyyyy'),'.png']));



%%% Make figure
figure(3);
clf;
plot(quad_drag,Tbc/1e6,'o-');
xlabel('Quadratic drag coefficient');
ylabel('Baroclinic transport (Sv)');
grid on;
print('-dpng',fullfile('figures',['BCTransportSensitivity_drag_',datestr(datenum(datetime('now')),'ddmmmyyyy'),'.png']));



%%% Make figure
figure(4);
plot(quad_drag,kap,'o-');
xlabel('Quadratic drag coefficient');
ylabel('Transient eddy diffusivity (m^2/s)');
grid on;
print('-dpng',fullfile('figures',['DiffusivitySensitivity_drag_',datestr(datenum(datetime('now')),'ddmmmyyyy'),'.png']));



%%% Make figure
figure(5);
plot(quad_drag,r_kap,'o-');
xlabel('Quadratic drag coefficient');
ylabel('Correlation coefficient for transient eddy diffusivity');
grid on;
print('-dpng',fullfile('figures',['DiffusivityCorrelationSensitivity_drag_',datestr(datenum(datetime('now')),'ddmmmyyyy'),'.png']));


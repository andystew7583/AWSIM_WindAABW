%%%
%%% plotLinDragSensitivity.m
%%%
%%% Plots sensitivity of ACC transport to wind for different linear bottom drag coefficients.
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
quad_drag = 0;
lin_drag = [1e-4:1e-4:1e-3];
topog_width = 150;
topog_height = 1000;
% quad_drag = 0e-3;
% lin_drag = 2e-4;
N_rb = length(lin_drag);

%%% Loop over runs and compute transport
Ttot = zeros(1,N_rb);
Tbt = zeros(1,N_rb);
Tbc = zeros(1,N_rb);
kap = zeros(1,N_rb);
r_kap = zeros(1,N_rb);
for n_rb = 1:N_rb
    
  [n_rb]
  %%% Simulation name
  run_name = constructRunName (is_spinup,Ny,Nlay, ...
                          tau_mean,tau_pert,tau_freq, ...
                          AABW_mean,AABW_pert,AABW_freq, ...
                          quad_drag,lin_drag(n_rb),topog_width,topog_height);
  loadParams;

  %%% Read time-mean zonal flux   
  hu_tavg = do_avg(dirpath,OUTN_HU_AVG,Nx,Ny,Nlay,n0_avg,N_avg,dt_avg,tmin,tmax,startTime);
  u_tavg = do_avg(dirpath,OUTN_U_AVG,Nx,Ny,Nlay,n0_avg,N_avg,dt_avg,tmin,tmax,startTime);

  %%% Compute transports
  Ttot(n_rb) = mean(sum(sum(hu_tavg,3),2),1)*dy;
  Tbt(n_rb) = mean(sum(u_tavg(:,:,end).*(-hhb).*dy,2),1);
  Tbc(n_rb) = Ttot(n_rb) - Tbt(n_rb);
  
  %%% Compute diffusivity
  [kap_bulk,nu_bulk,r_kap_bulk,r_nu_bulk,EKE_zavg] = calcBulkEddyViscDiff(local_home_dir,run_name);
  kap(n_rb) = kap_bulk;
  r_kap(n_rb) = r_kap_bulk;
   
end

%%% Make figure
figure(1);
clf;
plot(lin_drag,Ttot/1e6,'o-');
xlabel('Linear drag coefficient');
ylabel('Total transport (Sv)');
grid on;
print('-dpng',fullfile('figures',['TotalTransportSensitivity_lindrag_',datestr(datenum(datetime('now')),'ddmmmyyyy'),'.png']));



%%% Make figure
figure(2);
clf;
plot(lin_drag,Tbt/1e6,'o-');
xlabel('Linear drag coefficient');
ylabel('Barotropic transport (Sv)');
grid on;
print('-dpng',fullfile('figures',['BTTransportSensitivity_lindrag_',datestr(datenum(datetime('now')),'ddmmmyyyy'),'.png']));



%%% Make figure
figure(3);
clf;
plot(lin_drag,Tbc/1e6,'o-');
xlabel('Linear drag coefficient');
ylabel('Baroclinic transport (Sv)');
grid on;
print('-dpng',fullfile('figures',['BCTransportSensitivity_lindrag_',datestr(datenum(datetime('now')),'ddmmmyyyy'),'.png']));


%%% Make figure
figure(4);
plot(lin_drag,kap,'o-');
xlabel('Linear drag coefficient');
ylabel('Transient eddy diffusivity (m^2/s)');
grid on;
print('-dpng',fullfile('figures',['DiffusivitySensitivity_lindrag_',datestr(datenum(datetime('now')),'ddmmmyyyy'),'.png']));



%%% Make figure
figure(5);
plot(lin_drag,r_kap,'o-');
xlabel('Linear drag coefficient');
ylabel('Correlation coefficient for transient eddy diffusivity');
grid on;
print('-dpng',fullfile('figures',['DiffusivityCorrelationSensitivity_lindrag_',datestr(datenum(datetime('now')),'ddmmmyyyy'),'.png']));



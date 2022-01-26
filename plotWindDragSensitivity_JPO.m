%%%
%%% plotWindDragSensitivity_JPO.m
%%%
%%% Plots sensitivity of ACC transport and eddy properties to wind and drag
%%% for our JPO paper.
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
% tau_mean = [0.01 0.017 0.03 0.05 0.1 0.17 0.3];
tau_mean = [0.01 0.013 0.017 0.022 0.03 0.039 0.05 0.07 0.1 0.13 0.17 0.22 0.3];
tau_pert = 0;
tau_freq = 0;
% AABW_mean = [-1.5 -.75 0 .75 1.5];
AABW_mean = 0;
% AABW_mean = [-1.5 0 1.5];
AABW_pert = 0;
AABW_freq = 0;
quad_drag = 2e-3;
lin_drag = 0e-4;
% quad_drag = 0e-3;
% lin_drag = 2e-4;
topog_width = 150;
topog_height = 1000;
N_tm = length(tau_mean);

%%% Loop over runs and compute transport
Ttot = zeros(N_tm,N_Cd);
Tbt = zeros(N_tm,N_Cd);
Tbc = zeros(N_tm,N_Cd);
kap = zeros(N_tm,N_Cd);
r_kap = zeros(N_tm,N_Cd);

for n_tm = 1:N_tm
  for n_Cd = 1:N_Cd
    
    [n_tm,n_Cd]
    %%% Simulation name
    run_name = constructRunName (is_spinup,Ny,Nlay, ...
                            tau_mean(n_tm),tau_pert,tau_freq, ...
                            AABW_mean,AABW_pert,AABW_freq, ...
                            quad_drag(n_Cd),lin_drag,topog_width,topog_height);
    loadParams;
    
    %%% Read time-mean zonal flux
    hu_tavg = do_avg(dirpath,OUTN_HU_AVG,Nx,Ny,Nlay,n0_avg,N_avg,dt_avg,tmin,tmax,startTime);
    u_tavg = do_avg(dirpath,OUTN_U_AVG,Nx,Ny,Nlay,n0_avg,N_avg,dt_avg,tmin,tmax,startTime);
    
    %%% Compute transports
    Ttot(n_am,n_tm) = mean(sum(sum(hu_tavg,3),2),1)*dy;
    Tbt(n_am,n_tm) = mean(sum(u_tavg(:,:,end).*(-hhb).*dy,2),1);
    Tbc(n_am,n_tm) = Ttot(n_am,n_tm) - Tbt(n_am,n_tm);
    
    %%% Compute diffusivity
    [kap_bulk,nu_bulk,r_kap_bulk,r_nu_bulk,EKE_zavg] = calcBulkEddyViscDiff(local_home_dir,run_name);
    kap(n_am,n_tm) = kap_bulk;
    r_kap(n_am,n_tm) = r_kap_bulk;
    
  end
end

%%% Make figure
figure(1);
for n_am=1:N_am
  semilogx(tau_mean,Ttot(n_am,:)/1e6,'o-');
  if (n_am == 1)
  hold on;
  end
end
hold off;
xlabel('Wind Stress (N/m^2)');
ylabel('Total transport (Sv)');
set(gca,'XTick',tau_mean);
set(gca,'XLim',[0 0.3]);
title(['C_d = ',num2str(quad_drag),', r_b = ',num2str(lin_drag),' m/s']);
grid on;
legstr = {};
for n_am=1:N_am
  legstr = {legstr{:},['Taabw=',num2str(AABW_mean(n_am)),' Sv']}; 
end
legend(legstr,'Location','SouthEast');
print('-dpng',fullfile('figures',['TotalTransportSensitivity_Cd=',num2str(quad_drag),'_rb=',num2str(lin_drag),'_',datestr(datenum(datetime('now')),'ddmmmyyyy'),'.png']));



%%% Make figure
figure(2);
for n_am=1:N_am
  plot(tau_mean,Tbt(n_am,:)/1e6,'o-');
  if (n_am == 1)
  hold on;
  end
end
hold off;
xlabel('Wind Stress (N/m^2)');
ylabel('Barotropic transport (Sv)');
set(gca,'XTick',tau_mean);
set(gca,'XLim',[0 0.3]);
title(['C_d = ',num2str(quad_drag),', r_b = ',num2str(lin_drag),' m/s']);
grid on;
legstr = {};
for n_am=1:N_am
  legstr = {legstr{:},['Taabw=',num2str(AABW_mean(n_am)),' Sv']}; 
end
legend(legstr,'Location','SouthEast');
print('-dpng',fullfile('figures',['BTTransportSensitivity_Cd=',num2str(quad_drag),'_rb=',num2str(lin_drag),'_',datestr(datenum(datetime('now')),'ddmmmyyyy'),'.png']));



%%% Make figure
figure(3);
for n_am=1:N_am
  semilogx(tau_mean,Tbc(n_am,:)/1e6,'o-');
  if (n_am == 1)
  hold on;
  end
end
hold off;
xlabel('Wind Stress (N/m^2)');
ylabel('Baroclinic transport (Sv)');
set(gca,'XTick',tau_mean);
set(gca,'XLim',[0 0.3]);
title(['C_d = ',num2str(quad_drag),', r_b = ',num2str(lin_drag),' m/s']);
grid on;
legstr = {};
for n_am=1:N_am
  legstr = {legstr{:},['Taabw=',num2str(AABW_mean(n_am)),' Sv']}; 
end
legend(legstr,'Location','SouthEast');
print('-dpng',fullfile('figures',['BCTransportSensitivity_Cd=',num2str(quad_drag),'_rb=',num2str(lin_drag),'_',datestr(datenum(datetime('now')),'ddmmmyyyy'),'.png']));



%%% Make figure
figure(4);
for n_am=1:N_am
%   semilogx(tau_mean,kap(n_am,:),'o-');
  plot(tau_mean,kap(n_am,:),'o-');
  if (n_am == 1)
  hold on;
  end
end
hold off;
xlabel('Wind Stress (N/m^2)');
ylabel('Transient eddy diffusivity (m^2/s)');
set(gca,'XTick',tau_mean);
set(gca,'XLim',[0 0.3]);
title(['C_d = ',num2str(quad_drag),', r_b = ',num2str(lin_drag),' m/s']);
grid on;
legstr = {};
for n_am=1:N_am
  legstr = {legstr{:},['Taabw=',num2str(AABW_mean(n_am)),' Sv']}; 
end
legend(legstr,'Location','SouthEast');
print('-dpng',fullfile('figures',['DiffusivitySensitivity_Cd=',num2str(quad_drag),'_rb=',num2str(lin_drag),'_',datestr(datenum(datetime('now')),'ddmmmyyyy'),'.png']));



%%% Make figure
figure(5);
for n_am=1:N_am
  semilogx(tau_mean,r_kap(n_am,:),'o-');
%   plot(tau_mean,kap(n_am,:),'o-');
  if (n_am == 1)
  hold on;
  end
end
hold off;
xlabel('Wind Stress (N/m^2)');
ylabel('Correlation coefficient for transient eddy diffusivity');
set(gca,'XTick',tau_mean);
set(gca,'XLim',[0 0.3]);
title(['C_d = ',num2str(quad_drag),', r_b = ',num2str(lin_drag),' m/s']);
grid on;
legstr = {};
for n_am=1:N_am
  legstr = {legstr{:},['Taabw=',num2str(AABW_mean(n_am)),' Sv']}; 
end
legend(legstr,'Location','SouthEast');
print('-dpng',fullfile('figures',['DiffusivityCorrelationSensitivity_Cd=',num2str(quad_drag),'_rb=',num2str(lin_drag),'_',datestr(datenum(datetime('now')),'ddmmmyyyy'),'.png']));

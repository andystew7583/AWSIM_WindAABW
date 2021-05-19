%%%
%%% plotWindSensitivity.m
%%%
%%% Plots sensitivity of ACC transport to wind in various experiments.
%%%

%%% Load static definitions
constants;

%%% Directory to store runs
local_home_dir = './runs';

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
tau_mean = [0.01 0.017 0.03 0.05 0.1 0.17 0.3];
tau_pert = 0;
tau_freq = 0;
AABW_mean = [-1.5 -.75 .75 1.5];
% AABW_mean = [-1.5 0 1.5];
AABW_pert = 0;
% AABW_freq = 0;
quad_drag = 2e-3;
lin_drag = 0e-4;
% quad_drag = 0e-3;
% lin_drag = 2e-4;
N_am = length(AABW_mean);
N_tm = length(tau_mean);

%%% Loop over runs and compute transport
Ttot = zeros(N_am,N_tm);
Tbt = zeros(N_am,N_tm);
Tbc = zeros(N_am,N_tm);
for n_am = 1:N_am
  for n_tm = 1:N_tm
    
    [n_am n_tm]
    %%% Simulation name
    run_name = constructRunName (is_spinup,Ny,Nlay, ...
                            tau_mean(n_tm),tau_pert,tau_freq, ...
                            AABW_mean(n_am),AABW_pert,AABW_freq, ...
                            quad_drag,lin_drag);
    loadParams;
    
    %%% Read time-mean zonal flux
    hu_tavg = do_avg(dirpath,OUTN_HU_AVG,Nx,Ny,Nlay,n0_avg,N_avg,dt_avg,tmin,tmax,startTime);
    u_tavg = do_avg(dirpath,OUTN_U_AVG,Nx,Ny,Nlay,n0_avg,N_avg,dt_avg,tmin,tmax,startTime);
    
    %%% Compute transports
    Ttot(n_am,n_tm) = mean(sum(sum(hu_tavg,3),2),1)*dy;
    Tbt(n_am,n_tm) = mean(sum(u_tavg(:,:,end).*(-hhb).*dy,2),1);
    Tbc(n_am,n_tm) = Ttot(n_am,n_tm) - Tbt(n_am,n_tm);
    
  end
end

%%% Make figure
figure(1);
semilogx(tau_mean,Tbc(3,:)/1e6,'o-');
hold on;
semilogx(tau_mean,Tbc(1,:)/1e6,'o-');
semilogx(tau_mean,Tbc(2,:)/1e6,'o-');
semilogx(tau_mean,Tbc(4,:)/1e6,'o-');
% semilogx(tau_mean,Tbc(5,:)/1e6,'o-');
hold off;
xlabel('Wind Stress (N/m^2)');
ylabel('Baroclinic transport (Sv)');
set(gca,'XTick',tau_mean);
set(gca,'XLim',[0 0.3]);
title(['C_d = ',num2str(quad_drag),', r_b = ',num2str(lin_drag),' m/s']);
grid on;
legend(['Taabw=',num2str(AABW_mean(2)),' Sv'],['Taabw=',num2str(AABW_mean(1)),' Sv'],['Taabw=',num2str(AABW_mean(3)),' Sv'],'Location','SouthEast');
print('-dpng',fullfile('figures',['TransportSensitivity_Cd=',num2str(quad_drag),'_rb=',num2str(lin_drag),'_',datestr(datenum(datetime('now')),'ddmmmyyyy'),'.png']));

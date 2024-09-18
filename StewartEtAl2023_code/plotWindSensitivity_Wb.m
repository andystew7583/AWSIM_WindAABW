%%%
%%% plotWindSensitivity_Wb.m
%%%
%%% Plots sensitivity of ACC transport to wind for various topographic widths.
%%%

%%% Load static definitions
constants;

%%% Directory to store runs
local_home_dir = '/Volumes/Kilchoman/UCLA/Projects/AWSIM_WindAABW/runs';

%%% Spinup simulations are long and produce no diagnostic, diagnostic
%%% simulations output high-frequency diagnostics to resolve the forcing
%%% period
is_spinup = true;

%%% Grid resolution 
Ny = 128;
Nlay = 2;

%%% Averaging period
tmin = 180.5*t1year;
tmax = 200.5*t1year;

%%% Parameters defining the batch of runs to plot
tau_mean = [0.01 0.017 0.03 0.05 0.1 0.17 0.3];
tau_pert = 0;
tau_freq = 0;
% AABW_mean = [-1.5 -.75 0 .75 1.5];
% AABW_mean = [-1.5 0 1.5];
AABW_mean = 0;
AABW_pert = 0;
AABW_freq = 0;
quad_drag = 2e-3;
lin_drag = 0e-4;
topog_width = [40 80 150 300 600];
topog_height = 1500;
% quad_drag = 0e-3;
% lin_drag = 2e-4;
N_Wb = length(topog_width);
N_tm = length(tau_mean);

%%% Loop over runs and compute transport
Ttot = zeros(N_Wb,N_tm);
Tbt = zeros(N_Wb,N_tm);
Tbc = zeros(N_Wb,N_tm);
for n_Wb = 1:N_Wb
  for n_tm = 1:N_tm
    
    [n_Wb n_tm]
    %%% Simulation name
    run_name = constructRunName (is_spinup,Ny,Nlay, ...
                            tau_mean(n_tm),tau_pert,tau_freq, ...
                            AABW_mean,AABW_pert,AABW_freq, ...
                            quad_drag,lin_drag,topog_width(n_Wb),topog_height);
    loadParams;
    
    %%% Read time-mean zonal flux   
    Nframes = Nframes + n0;
    n0 = 0;
    startTime = 0;
    h_tavg = do_avg(dirpath,OUTN_H,Nx,Ny,Nlay,n0,Nframes,dt_s,tmin,tmax,startTime);
    u_tavg = do_avg(dirpath,OUTN_U,Nx,Ny,Nlay,n0,Nframes,dt_s,tmin,tmax,startTime);
    hu_tavg = h_tavg.*u_tavg;
    
    %%% Compute transports
    Ttot(n_Wb,n_tm) = mean(sum(sum(hu_tavg,3),2),1)*dy;
    Tbt(n_Wb,n_tm) = mean(sum(u_tavg(:,:,end).*(-hhb).*dy,2),1);
    Tbc(n_Wb,n_tm) = Ttot(n_Wb,n_tm) - Tbt(n_Wb,n_tm);
    
  end
end

%%% Make figure
figure(1);
for n_Wb=1:N_Wb
  semilogx(tau_mean,Ttot(n_Wb,:)/1e6,'o-');
  if (n_Wb == 1)
  hold on;
  end
end
hold off;
xlabel('Wind Stress (N/m^2)');
ylabel('Total transport (Sv)');
set(gca,'XTick',tau_mean);
set(gca,'XLim',[0 0.3]);
title(['H_b = ',num2str(topog_height),' m']);
grid on;
legstr = {};
for n_Wb=1:N_Wb
  legstr = {legstr{:},['Wb=',num2str(topog_width(n_Wb)),' km']}; 
end
legend(legstr,'Location','SouthEast');
print('-dpng',fullfile('figures',['TotalTransportSensitivity_Hb=',num2str(topog_height),'_',datestr(datenum(datetime('now')),'ddmmmyyyy'),'.png']));



%%% Make figure
figure(2);
for n_Wb=1:N_Wb
  semilogx(tau_mean,Tbt(n_Wb,:)/1e6,'o-');
  if (n_Wb == 1)
  hold on;
  end
end
hold off;
xlabel('Wind Stress (N/m^2)');
ylabel('Barotropic transport (Sv)');
set(gca,'XTick',tau_mean);
set(gca,'XLim',[0 0.3]);
title(['H_b = ',num2str(topog_height),' m']);
grid on;
legstr = {};
for n_Wb=1:N_Wb
  legstr = {legstr{:},['Wb=',num2str(topog_width(n_Wb)),' km']}; 
end
legend(legstr,'Location','SouthEast');
print('-dpng',fullfile('figures',['BTTransportSensitivity_Hb=',num2str(topog_height),'_',datestr(datenum(datetime('now')),'ddmmmyyyy'),'.png']));



%%% Make figure
figure(3);
for n_Wb=1:N_Wb
  semilogx(tau_mean,Tbc(n_Wb,:)/1e6,'o-');
  if (n_Wb == 1)
  hold on;
  end
end
hold off;
xlabel('Wind Stress (N/m^2)');
ylabel('Baroclinic transport (Sv)');
set(gca,'XTick',tau_mean);
set(gca,'XLim',[0 0.3]);
title(['H_b = ',num2str(topog_height),' m']);
grid on;
legstr = {};
for n_Wb=1:N_Wb
  legstr = {legstr{:},['Wb=',num2str(topog_width(n_Wb)),' km']}; 
end
legend(legstr,'Location','SouthEast');
print('-dpng',fullfile('figures',['BCTransportSensitivity_Hb=',num2str(topog_height),'_',datestr(datenum(datetime('now')),'ddmmmyyyy'),'.png']));

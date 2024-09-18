%%% Select run
run_name = 'ACC_AABW_Ny256_Nlay2_tauM0.5_tauP0_tauF0_wDiaM0_wDiaP0_wDiaF0_Cd0.000e+00_rb2.000e-04_diags'; %%% Run to analyze
local_home_dir = '/Volumes/Kilchoman/UCLA/Projects/AWSIM_WindAABW/runs';
dirpath = fullfile(local_home_dir,run_name);

%%% Load parameters   
loadParams;
gtild = reshape(cumsum(gg),[Nlay 1 1]);

%%% Averaging period
tmin = 0.5*t1year;
tmax = 30.5*t1year;

%%% Time-average required model output variables
MM_tavg = do_avg(dirpath,OUTN_M_AVG,Nx,Ny,Nlay,n0_avg,N_avg,dt_avg,tmin,tmax,startTime);
hu_tavg = do_avg(dirpath,OUTN_HU_AVG,Nx,Ny,Nlay,n0_avg,N_avg,dt_avg,tmin,tmax,startTime);
hv_tavg = do_avg(dirpath,OUTN_HV_AVG,Nx,Ny,Nlay,n0_avg,N_avg,dt_avg,tmin,tmax,startTime);

Psi = -cumsum(hu_tavg,2)*dy;
div = (hu_tavg([2:Nx 1],:,:) - hu_tavg(:,:,:))/dx + (hv_tavg(:,[2:Ny 1],:) - hv_tavg(:,:,:))/dy;

%%% Plot lower layer circulation
figure(101);
pcolor(XX_h,YY_h,MM_tavg(:,:,2));
shading interp;
colorbar;
colormap(cmocean('balance'));


%%% Plot lower layer circulation
figure(102);
clf;
set(gcf,'Position',[400 400 850 400]);
pcolor(XX_h/1000,YY_h/1000,Psi(:,:,2)/1e6);
shading interp;
cbhandle = colorbar;
set(get(cbhandle,'title'),'String','Sv');
colormap(cmocean('balance'));
caxis([-150 150]);
set(gca,'FontSize',14);
set(gca,'Position',[0.08 0.1 0.82 0.82]);
xlabel('x (km)');
ylabel('y (km)');
title('\tau_m_a_x = 0.5 N/m^2');


%%% Plot lower layer circulation
figure(103);
clf;
set(gcf,'Position',[400 400 850 400]);
pcolor(XX_h/1000,YY_h/1000,-(MM_tavg(:,:,2)-mean(mean(MM_tavg(:,1,2)))));
shading interp;
cbhandle = colorbar;
set(get(cbhandle,'title'),'String','m^2/s^2');
colormap(cmocean('balance'));
% caxis([-150 150]);
set(gca,'FontSize',14);
set(gca,'Position',[0.08 0.1 0.82 0.82]);
xlabel('x (km)');
ylabel('y (km)');
title('\tau_m_a_x = 0.5 N/m^2');

%%% Plot lower layer circulation
figure(104);
clf;
set(gcf,'Position',[400 400 850 400]);
pcolor(XX_h/1000,YY_h/1000,div(:,:,2));
shading interp;
cbhandle = colorbar;
set(get(cbhandle,'title'),'String','m/s');
colormap(cmocean('balance'));
% caxis([-150 150]);
set(gca,'FontSize',14);
set(gca,'Position',[0.08 0.1 0.82 0.82]);
xlabel('x (km)');
ylabel('y (km)');
title('\tau_m_a_x = 0.5 N/m^2');

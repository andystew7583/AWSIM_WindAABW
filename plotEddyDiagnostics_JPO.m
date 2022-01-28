%%%
%%% plotEddyDiagnostics_JPO.m
%%%
%%% Plots diagnostics of the eddy field for our JPO paper.
%%%

%%% Run to load
run_name = 'ACC_AABW_Ny256_Nlay2_tauM0.1_tauP0_tauF0_wDiaM0_wDiaP0_wDiaF0_Cd2.000e-03_rb0.000e+00_diags';

%%% Load parameters   
local_home_dir = '/Volumes/Kilchoman/UCLA/Projects/AWSIM_WindAABW/runs';
prod_dir = fullfile('./products');
loadParams;
dirpath = fullfile(local_home_dir,run_name);
gtild = reshape(cumsum(gg),[Nlay 1 1]);


%%% Load pre-computed eddy diagnostics
load(fullfile(prod_dir,['kap_nu_',run_name,'.mat']));



















%%% Plotting options
fontsize = 14;
axpos = zeros(4,4);
axpos(1,:) = [0.11 0.71 .85 .26];
axpos(2,:) = [0.11 0.38 .85 .26];
axpos(3,:) = [0.11 0.05 .85 .26];
axlabels = {'(a)','(b)','(c)','(d)'};
tau_ticks = [0.01 0.017 0.03 0.05 0.1 0.17 0.3];
rho0 = 1000;

figure(103);
clf;
set(gcf,'Position',[382   306   700   800]);



EKE_thresh = quantile(EKE_zavg(:),0.75);

subplot('Position',axpos(1,:));
pcolor(XX_h/m1km,YY_h/m1km,EKE_zavg);
shading interp;
hold on;
[C,h]=contour(XX_h/m1km,YY_h/m1km,-hhb,[500 1500 2500 3250 3500 3750],'EdgeColor',[.5 .5 .5],'LineWidth',1);
contour(XX_h/m1km,YY_h/m1km,EKE_zavg,[EKE_thresh EKE_thresh],'EdgeColor','k');
hold off;
% xlabel('Longitude $x$ (km)','interpreter','latex');
ylabel('Latitude $y$ (km)','interpreter','latex');
colorbar;
caxis([0 0.025]);
colormap(gca,cmocean('amp',25));
set(gca,'FontSize',fontsize);
title('Depth-averaged EKE (m^2/s^2)');

subplot('Position',axpos(2,:));
pcolor(XX_h/m1km,YY_h/m1km,kap_map);
shading interp;
hold on;
[C,h]=contour(XX_h/m1km,YY_h/m1km,-hhb,[500 1500 2500 3250 3500 3750],'EdgeColor',[.5 .5 .5],'LineWidth',1);
contour(XX_h/m1km,YY_h/m1km,EKE_zavg,[EKE_thresh EKE_thresh],'EdgeColor','k');
hold off;
% xlabel('Longitude $x$ (km)','interpreter','latex');
ylabel('Latitude $y$ (km)','interpreter','latex');
colorbar;
caxis([-1000 1000]);
colormap(gca,cmocean('balance'));
set(gca,'FontSize',fontsize);
title('Transient eddy diffusivity, \kappa (m^2/s)');

subplot('Position',axpos(3,:));
pcolor(XX_h/m1km,YY_h/m1km,100*r_kap_map.^2);
shading interp;
hold on;
[C,h]=contour(XX_h/m1km,YY_h/m1km,-hhb,[500 1500 2500 3250 3500 3750],'EdgeColor',[.5 .5 .5],'LineWidth',1);
contour(XX_h/m1km,YY_h/m1km,EKE_zavg,[EKE_thresh EKE_thresh],'EdgeColor','k');
hold off;
xlabel('Longitude $x$ (km)','interpreter','latex');
ylabel('Latitude $y$ (km)','interpreter','latex');
colorbar;
caxis([0 100]);
colormap(gca,cmocean('amp'));
set(gca,'FontSize',fontsize);
title('Variance explained by \kappa (%)');


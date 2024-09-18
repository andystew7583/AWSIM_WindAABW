%%%
%%% plotMeanStreamfunction_JPO.m
%%%
%%% Plots wind forcing, instantaneous model state and bathymetry for our JPO paper.
%%%

%%% Run to load
run_name = 'ACC_AABW_Ny256_Nlay2_tauM0.1_tauP0_tauF0_wDiaM0_wDiaP0_wDiaF0_Cd2.000e-03_rb0.000e+00_diags';

%%% Load parameters   
local_home_dir = '/Volumes/Kilchoman/UCLA/Projects/AWSIM_WindAABW/runs';
prod_dir = fullfile(local_home_dir,'AWSIM_WindAABW_products');
loadParams;
dirpath = fullfile(local_home_dir,run_name);
gtild = reshape(cumsum(gg),[Nlay 1 1]);

%%% Averaging period
tmin = 0.5*t1year;
tmax = 30.5*t1year;

%%% Load hu
hu_tavg = do_avg(dirpath,OUTN_HU_AVG,Nx,Ny,Nlay,n0_avg,N_avg,dt_avg,tmin,tmax,startTime);

%%% Calculate streamfunction
Psi = -cumsum(hu_tavg,2)*dy;





%%% Plotting options
scrsz = get(0,'ScreenSize');
framepos = [209   632   900   800];
fontsize = 14;
ax_pos = zeros(2,4);
ax_pos(1,:) = [0.07 0.55 0.9 0.4];
ax_pos(2,:) = [0.26 0.45 0.65 0.5];
cb_pos = [0.94 0.45 0.015 0.5];
lab_size = [0.05 0.03];


%%% Set up the frame
figure(107);
clf; 
set(gcf,'Position',framepos);
set(gcf,'Color','w');
cntr = 0;



%%% Plot wind stress
cntr = cntr + 1;
ax = axes('position',ax_pos(cntr,:));  
pcolor(XX_h/1000,YY_h/1000,Psi(:,:,1)/1e6);
shading interp;
colorbar
caxis([-60 60]);
colormap(cmocean('balance',60));
ylabel('Latitude y (km)');
xlabel('Longitude x (km)');
set(gca,'FontSize',fontsize);
annotation('textbox',[ax_pos(cntr,1)-0.06 ax_pos(cntr,2)-0.06 lab_size],'String','(a)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

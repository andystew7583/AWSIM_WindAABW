%%%
%%% plotMeanEddyIFS.m
%%%
%%% Plots fluxes across latitude lines across a range of experiments.
%%%

%%% Load static definitions
constants;

%%% Directory to store runs
local_home_dir = '/Volumes/Kilchoman/UCLA/Projects/AWSIM_WindAABW/runs';
prod_dir = fullfile('./products');

%%% Averaging period
tmin = 0.5*t1year;
tmax = 30.5*t1year;




%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameter selection %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ny = 256;
Nx = Ny*2;
Nlay = 2;
is_spinup = false;
tau_mean = [0.01 0.013 0.017 0.022 0.03 0.039 0.05 0.07 0.1 0.13 0.17 0.22 0.3 0.39 0.5];
tau_pert = 0;
tau_freq = 0;
AABW_mean = 0;
AABW_pert = 0;
AABW_freq = 0;
quad_drag = [.5e-3 1e-3 1.5e-3 2e-3 2.5e-3 3e-3 3.5e-3 4e-3];
% quad_drag = 2e-3;
% quad_drag = 0e-3;
lin_drag = 0e-4;
topog_width = 150;
topog_height = 1000;

%%% Index of reference simulation
idx_ref_tm = 9;
idx_ref_Cd = 4;
idx_lat = 1:Ny;

%%% Loop over runs
N_tm = length(tau_mean);
N_Cd = length(quad_drag);
kap_map_batch = zeros(Nx,Ny,N_Cd,N_tm);
diff_curl_twa_batch = zeros(Nx,Ny,N_Cd,N_tm);
MM_surf_batch = zeros(Nx,Ny,N_Cd,N_tm);
tau_cntr_batch = zeros(N_Cd,N_tm);
hgradM_eddy_cntr_batch = zeros(N_Cd,N_tm);
hgradM_mean_cntr_batch = zeros(N_Cd,N_tm);
eddyforce_cntr_batch = zeros(N_Cd,N_tm);
Klat_cntr_batch = zeros(N_Cd,N_tm);
for n_tm=1:N_tm
  for n_Cd=1:N_Cd
%   for n_Cd=idx_ref_Cd
    
    %%% Generate simulation name
    run_name = constructRunName (is_spinup,Ny,Nlay, ...
                tau_mean(n_tm),tau_pert,tau_freq, ...
                AABW_mean,AABW_pert,AABW_freq, ...
                quad_drag(n_Cd),lin_drag,topog_width,topog_height);

    %%% Load parameters   
    loadParams;
    dirpath = fullfile(local_home_dir,run_name);
    gtild = reshape(cumsum(gg),[Nlay 1 1]);
    f0 = mean(mean(2*Omega_z));

    %%% Load eddy diffusivity and viscosity maps
    load(fullfile(prod_dir,['kap_nu_',run_name,'.mat']));
    kap_map_batch(:,:,n_Cd,n_tm) = kap_map;

    %%% Time-average required model output variables
    hh_tavg = do_avg(dirpath,OUTN_H_AVG,Nx,Ny,Nlay,n0_avg,N_avg,dt_avg,tmin,tmax,startTime);
    MM_tavg = do_avg(dirpath,OUTN_M_AVG,Nx,Ny,Nlay,n0_avg,N_avg,dt_avg,tmin,tmax,startTime);
    hu_tavg = do_avg(dirpath,OUTN_HU_AVG,Nx,Ny,Nlay,n0_avg,N_avg,dt_avg,tmin,tmax,startTime);
    hv_tavg = do_avg(dirpath,OUTN_HV_AVG,Nx,Ny,Nlay,n0_avg,N_avg,dt_avg,tmin,tmax,startTime);
    husq_tavg = do_avg(dirpath,OUTN_HUU_AVG,Nx,Ny,Nlay,n0_avg,N_avg,dt_avg,tmin,tmax,startTime);
    hvsq_tavg = do_avg(dirpath,OUTN_HVV_AVG,Nx,Ny,Nlay,n0_avg,N_avg,dt_avg,tmin,tmax,startTime);
    huv_tavg = do_avg(dirpath,OUTN_HUV_AVG,Nx,Ny,Nlay,n0_avg,N_avg,dt_avg,tmin,tmax,startTime);
    hdMdx_tavg = -do_avg(dirpath,OUTN_UMOM_GRADM,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
    hdMdy_tavg = -do_avg(dirpath,OUTN_VMOM_GRADM,Nx,Ny,Nlay,n0_avg_hv,N_avg_hv,dt_avg_hv,tmin,tmax,startTime);

    %%% Compute mean layer thicknesses on cell edges
    hh_w_tavg = 0.5*(hh_tavg(1:Nx,:,:)+hh_tavg([Nx 1:Nx-1],:,:)); %%% N.B. assumes double periodicity
    hh_s_tavg = 0.5*(hh_tavg(:,1:Ny,:)+hh_tavg(:,[Ny 1:Ny-1],:)); 
    hh_q_tavg = 0.25.*(hh_tavg(1:Nx,1:Ny,:)+hh_tavg(1:Nx,[Ny 1:Ny-1],:)+hh_tavg([Nx 1:Nx-1],1:Ny,:)+hh_tavg([Nx 1:Nx-1],[Ny 1:Ny-1],:));

    %%% Thickness-weighted average velocities
    uu_twa = hu_tavg ./ hh_w_tavg;
    vv_twa = hv_tavg ./ hh_s_tavg;

    %%% Compute eddy form stress terms
    hdMdx_mean = hh_w_tavg.*(MM_tavg(1:Nx,:,:)-MM_tavg([Nx 1:Nx-1],:,:))/dx;
    hdMdx_eddy = hdMdx_tavg - hdMdx_mean;
    hdMdy_mean = hh_s_tavg.*(MM_tavg(:,1:Ny,:)-MM_tavg(:,[Ny 1:Ny-1],:))/dy;
    hdMdy_eddy = hdMdy_tavg - hdMdy_mean;

    %%% Compute eddy momentum flux terms
    husq_mean = hh_w_tavg.*0.5.*(uu_twa(1:Nx,:,:).^2+uu_twa([2:Nx 1],:,:).^2);
    hvsq_mean = hh_s_tavg.*0.5.*(vv_twa(:,1:Ny,:).^2+vv_twa(:,[2:Ny 1],:).^2);
    huv_mean = 0.5.*(uu_twa(1:Nx,:,:)+uu_twa([Nx 1:Nx-1],:,:)) ...
                    .* 0.5.*(vv_twa(:,1:Ny,:)+vv_twa(:,[Ny 1:Ny-1],:)) ...
                    .* hh_q_tavg;                  
    husq_eddy = husq_tavg - husq_mean;
    hvsq_eddy = hvsq_tavg - hvsq_mean;
    huv_eddy = huv_tavg - huv_mean;

    %%% Compute advective eddy forcing of the mean flow
    dx_husq_eddy = (husq_eddy([2:Nx 1],:,:)-husq_eddy([Nx 1:Nx-1],:,:))/(2*dx);
    dy_hvsq_eddy = (hvsq_eddy(:,[2:Ny 1],:)-hvsq_eddy(:,[Ny 1:Ny-1],:))/(2*dy);
    dx_huv_eddy = (huv_eddy([2:Nx 1],:,:)-huv_eddy(1:Nx,:,:))/dx;
    dy_huv_eddy = (huv_eddy(:,[2:Ny 1],:)-huv_eddy(:,1:Ny,:))/dy;
    eddyforce_x = dx_husq_eddy + dy_huv_eddy;
    eddyforce_y = dx_huv_eddy + dy_hvsq_eddy;

    %%% Depth-integrated EKE
    EKE = 0.25*(husq_eddy(1:Nx,:,:)+husq_eddy([2:Nx 1],:,:)+hvsq_eddy(:,1:Ny,:)+hvsq_eddy(:,[2:Ny 1],:));
    EKE_zint = sum(EKE,3);
    EKE_zavg = EKE_zint ./ sum(hh_tavg,3);
    idx_EKE = find(EKE_zavg>quantile(EKE_zavg(:),0.75));

    %%% Remove redundant time dimension in wind stress vector
    taux = squeeze(mean(taux,1));
    tauy = squeeze(mean(tauy,1));

    %%% Along-streamline averages
    %%% TODO do averages over latitudes
    tau_cntr_batch(n_Cd,n_tm) = sum(mean(taux(:,idx_lat),2)*dx);
    hgradM_eddy_cntr_batch(n_Cd,n_tm) = sum(mean(hdMdx_eddy(:,idx_lat),2)*dx);
    hgradM_mean_cntr_batch(n_Cd,n_tm) = sum(mean(hdMdx_mean(:,idx_lat),2)*dx);
    eddyforce_cntr_batch(n_Cd,n_tm) = sum(mean(eddyforce_x(:,idx_lat),2)*dx);
    Klat_cntr_batch(n_Cd,n_tm) = hgradM_eddy_cntr_batch(n_Cd,n_tm)/Lx / (f0^2/gg(2)*mean(mean(uu_twa(:,idx_lat,1)-uu_twa(:,idx_lat,2),2),1));

  end
end

%%%%%%%%%%%%%%%%%%
%%% MAKE PLOTS %%%
%%%%%%%%%%%%%%%%%%

%%% Plotting options

fontsize = 14;
axpos = zeros(2,4);
axpos(1,:) = [0.075 0.11 .41 .85];
axpos(2,:) = [0.575 0.11 .41 .85];
lab_size = [0.05 0.03];
tau_ticks = [5e-3 0.01 0.017 0.03 0.05 0.1 0.17 0.3 0.5];
rho0 = 1000;

defaultcolororder = get(gca,'ColorOrder');
colororder = zeros(8,3);
colororder(1:idx_ref_Cd-1,:) = defaultcolororder(1:idx_ref_Cd-1,:);
colororder(idx_ref_Cd,:) = [0 0 0];
colororder(idx_ref_Cd+1:N_Cd,:) = defaultcolororder(idx_ref_Cd:N_Cd-1,:);

markersize = 14;
markershapes = {'>','o','*','<','v','d','^','s','x','+'};

figure(110);
clf;
set(gcf,'Position',[382   306   900   400]);

subplot('Position',axpos(1,:));
loglog(tau_mean,rho0*tau_cntr_batch(idx_ref_Cd,:)/Lx,'o-','LineWidth',1.5);
hold on;
loglog(tau_mean,rho0*hgradM_mean_cntr_batch(idx_ref_Cd,:)/Lx,'o-','LineWidth',1.5);
loglog(tau_mean,rho0*hgradM_eddy_cntr_batch(idx_ref_Cd,:)/Lx,'o-','LineWidth',1.5);
% loglog(tau_mean,rho0*eddyforce_cntr_batch(idx_ref_Cd,:)/Lx,'o-','LineWidth',1.5);
% loglog(tau_mean,rho0*(hgradM_eddy_cntr_batch(idx_ref_Cd,:)+hgradM_mean_cntr_batch(idx_ref_Cd,:)+eddyforce_cntr_batch(idx_ref_Cd,:))/Lx,'o-','LineWidth',1.5);
% loglog([.1 .1],[5e-3 .5],'k:','LineWidth',2);
hold off;
ylabel('Total transport (Sv)');
set(gca,'XTick',tau_ticks);
set(gca,'XLim',[0 0.5]);
% set(gca,'YLim',[10 90]);
set(gca,'FontSize',fontsize);
grid on;
xlabel('Wind stress maximum (N/m$^2$)','interpreter','latex');
ylabel('Area-averaged zonal stress (N/m$^2$)','interpreter','latex');
leghandle = legend('Wind stress','Mean IFS','Eddy IFS');
set(leghandle,'interpreter','latex','Location','NorthWest');
grid on;
annotation('textbox',[axpos(1,1)-0.07 axpos(1,2)-0.04 lab_size],'String','(a)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

subplot('Position',axpos(2,:));
for n_Cd = 1:N_Cd
  semilogx(tau_mean,Klat_cntr_batch(n_Cd,:),'o-','Color',colororder(n_Cd,:),'Marker',markershapes{n_Cd},'MarkerFaceColor',colororder(n_Cd,:));
  if (n_Cd == 1)
    hold on;
  end
end
hold off;
axis tight
grid on;
set(gca,'XLim',[0 0.5]);
set(gca,'XTick',tau_ticks);
xlabel('Wind stress maximum (N/m$^2$)','interpreter','latex');
ylabel('Meridional eddy diffusivity (m$^2$/s)','interpreter','latex');
set(gca,'FontSize',fontsize);

legstr = {};
for n_Cd=1:N_Cd
  legstr = {legstr{:},['Cd=',num2str(quad_drag(n_Cd)*1000),'$\times$ 10$^{-3}$']}; 
end
leghandle = legend(legstr,'Location','SouthWest');
set(leghandle,'interpreter','latex');

annotation('textbox',[axpos(2,1)-0.07 axpos(2,2)-0.04 lab_size],'String','(b)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
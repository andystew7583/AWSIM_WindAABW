%%%
%%% plotCrossStreamlineTransport_JPO.m
%%%
%%% Plots fluxes across streamlines across a range of experiments for our
%%% JPO paper.
%%%

%%% Load static definitions
constants;

%%% Directory to store runs
local_home_dir = '/Volumes/Kilchoman/UCLA/Projects/AWSIM_WindAABW/runs';
prod_dir = fullfile('./products');

%%% Averaging period
tmin = 0.5*t1year;
tmax = 30.5*t1year;

%%% Index of reference simulation
idx_ref_tm = 9;
idx_ref_Cd = 4;




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
for n_tm=1:N_tm
  for n_Cd=1:N_Cd
    
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
    husq_mean = hh_w_tavg.*uu_twa.^2;
    hvsq_mean = hh_s_tavg.*vv_twa.^2;
    huv_mean = 0.5.*(vv_twa(1:Nx,:,:)+vv_twa([Nx 1:Nx-1],:,:)) ...
                    .* 0.5.*(uu_twa(:,1:Ny,:)+uu_twa(:,[Ny 1:Ny-1],:)) ...
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

    %%% Momentum curl budget quantities
    curl_tau = (tauy(1:Nx,:,1)-tauy([Nx 1:Nx-1],:,1))/dx - (taux(:,1:Ny,1)-taux(:,[Ny 1:Ny-1],1))/dy;
    curl_tau(:,1) = 0;
    curl_hgradM_eddy = (hdMdy_eddy(1:Nx,:,1)-hdMdy_eddy([Nx 1:Nx-1],:,1))/dx - (hdMdx_eddy(:,1:Ny,1)-hdMdx_eddy(:,[Ny 1:Ny-1],1))/dy;
    curl_hgradM_eddy(:,1) = 0;
    curl_hgradM_eddy(:,2,:) = (hdMdy_eddy(1:Nx,2,1)-hdMdy_eddy([Nx 1:Nx-1],2,1))/dx - (hdMdx_eddy(:,2,1))/dy;
    curl_hgradM_mean = (hdMdy_mean(1:Nx,:,1)-hdMdy_mean([Nx 1:Nx-1],:,1))/dx - (hdMdx_mean(:,1:Ny,1)-hdMdx_mean(:,[Ny 1:Ny-1],1))/dy;
    curl_hgradM_mean(:,1) = 0;
    curl_hgradM_mean(:,2,:) = (hdMdy_mean(1:Nx,2,1)-hdMdy_mean([Nx 1:Nx-1],2,1))/dx - (hdMdx_mean(:,2,1))/dy;
    curl_eddyforce = (eddyforce_y(1:Nx,:,1)-eddyforce_y([Nx 1:Nx-1],:,1))/dx - (eddyforce_x(:,1:Ny,1)-eddyforce_x(:,[Ny 1:Ny-1],1))/dy;
    curl_eddyforce(:,1) = 0;
    curl_eddyforce(:,2,:) = (eddyforce_y(1:Nx,2,1)-eddyforce_y([Nx 1:Nx-1],2,1))/dx - (eddyforce_x(:,2,1))/dy;
    curl_twa = (vv_twa(1:Nx,:,:)-vv_twa([Nx 1:Nx-1],:,:))/dx - (uu_twa(:,1:Ny,:)-uu_twa(:,[Ny 1:Ny-1],:))/dy;
    curl_twa(:,1,:) = 0;
    curl_twa(:,2,:) = (vv_twa(1:Nx,2,:)-vv_twa([Nx 1:Nx-1],2,:))/dx - (uu_twa(:,2,:))/dy;
    diff_curl_twa = diff(curl_twa,1,3);
    diff_curl_twa_batch(:,:,n_Cd,n_tm) = diff_curl_twa;

    %%% Select streamline for calculations
    MM_surf = MM_tavg(:,:,1);  
    MM_surf(:,2:Ny,:) = 0.25*(MM_surf(1:Nx,1:Ny-1,:)+MM_surf([Nx 1:Nx-1],1:Ny-1,:)+MM_surf(1:Nx,2:Ny,:)+MM_surf([Nx 1:Nx-1],2:Ny,:));
    MM_surf(:,1,:) = 0.5*(MM_surf(1:Nx,1,:)+MM_surf([Nx 1:Nx-1],1,:));
    MM_cntr = median(MM_surf(:)); 
    MM_surf_batch(:,:,n_Cd,n_tm) = MM_surf;

    %%% Along-streamline averages
    tau_cntr_batch(n_Cd,n_tm) = sum(sum(curl_tau(MM_surf<MM_cntr)*dx*dy));
    hgradM_eddy_cntr_batch(n_Cd,n_tm) = sum(sum(curl_hgradM_eddy(MM_surf<MM_cntr)*dx*dy));
    hgradM_mean_cntr_batch(n_Cd,n_tm) = sum(sum(curl_hgradM_mean(MM_surf<MM_cntr)*dx*dy));
    eddyforce_cntr_batch(n_Cd,n_tm) = sum(sum(curl_eddyforce(MM_surf<MM_cntr)*dx*dy));  

  end
end



%%% Compute estimated eddy contributions
hgradM_kap_cntr_batch = zeros(N_Cd,N_tm);
hgradM_refkap_cntr_batch = zeros(N_Cd,N_tm);
hgradM_refflow_cntr_batch = zeros(N_Cd,N_tm);
for n_tm=1:N_tm
  for n_Cd=1:N_Cd
    
    

    %%% Estimated eddy flux using kappa and mean flow from each simulation  
    MM_surf = MM_surf_batch(:,:,n_Cd,n_tm);
    MM_cntr = median(MM_surf(:));
    IPT_kap = kap_map_batch(:,:,n_Cd,n_tm).*(f0^2/gg(2)).*diff_curl_twa_batch(:,:,n_Cd,n_tm);  
    hgradM_kap_cntr_batch(n_Cd,n_tm) = sum(sum(IPT_kap(MM_surf<MM_cntr)*dx*dy));

    %%% Estimated eddy flux using reference kappa  
    MM_surf = MM_surf_batch(:,:,n_Cd,n_tm);
    MM_cntr = median(MM_surf(:));
    IPT_kap = kap_map_batch(:,:,idx_ref_Cd,idx_ref_tm).*(f0^2/gg(2)).*diff_curl_twa_batch(:,:,n_Cd,n_tm);  
    hgradM_refkap_cntr_batch(n_Cd,n_tm) = sum(sum(IPT_kap(MM_surf<MM_cntr)*dx*dy));

    %%% Estimated eddy flux using reference mean flow  
    MM_surf = MM_surf_batch(:,:,idx_ref_Cd,idx_ref_tm); %%% Reference Montgomery potential
    MM_cntr = median(MM_surf(:));
    IPT_kap = kap_map_batch(:,:,n_Cd,n_tm).*(f0^2/gg(2)).*diff_curl_twa_batch(:,:,idx_ref_Cd,idx_ref_tm);  
    hgradM_refflow_cntr_batch(n_Cd,n_tm) = sum(sum(IPT_kap(MM_surf<MM_cntr)*dx*dy));
  
  end
end






%%%%%%%%%%%%%%%%%%
%%% MAKE PLOTS %%%
%%%%%%%%%%%%%%%%%%

%%% Plotting options

fontsize = 14;
axpos = zeros(4,4);
axpos(1,:) = [0.075 0.56 .41 .42];
axpos(2,:) = [0.575 0.56 .41 .42];
axpos(3,:) = [0.075 0.06 .41 .42];
axpos(4,:) = [0.575 0.06 .41 .42];
axlabels = {'(a)','(b)','(c)','(d)'};
lab_size = [0.05 0.03];
tau_ticks = [5e-3 0.01 0.017 0.03 0.05 0.1 0.17 0.3 0.5];
rho0 = 1000;
linewidth = 1.5;

figure(106);
clf;
set(gcf,'Position',[382   306   888   679]);

subplot('Position',axpos(1,:));
loglog(tau_mean,-rho0*tau_cntr_batch(idx_ref_Cd,:)/Lx,'o-','LineWidth',1.5);
hold on;
loglog(tau_mean,-rho0*hgradM_eddy_cntr_batch(idx_ref_Cd,:)/Lx,'o-','LineWidth',1.5);
loglog(tau_mean,-rho0*eddyforce_cntr_batch(idx_ref_Cd,:)/Lx,'o-','LineWidth',1.5);
loglog(tau_mean,-rho0*(hgradM_eddy_cntr_batch(idx_ref_Cd,:)+eddyforce_cntr_batch(idx_ref_Cd,:))/Lx,'o-','LineWidth',1.5);
loglog([.1 .1],[5e-3 .5],'k:','LineWidth',2);
hold off;
axis([0.01 .5 5e-3 .5]);
set(gca,'XTick',tau_ticks);
set(gca,'YTick',tau_ticks);
set(gca,'FontSize',fontsize);
xlabel('Wind stress maximum (N/m$^2$)','interpreter','latex');
ylabel('Along-streamline forcing (N/m$^2$)','interpreter','latex');
leghandle = legend('Wind stress','Eddy IFS','Eddy advection','Total eddy force');
set(leghandle,'interpreter','latex','Location','NorthWest');
grid on;

colororder = get(gca,'ColorOrder');

subplot('Position',axpos(2,:));
loglog(tau_mean,-rho0*hgradM_eddy_cntr_batch(idx_ref_Cd,:)/Lx,'o-','Color',colororder(5,:),'LineWidth',1.5);
hold on;
% semilogx(tau_mean,rho0*hgradM_kap_cntr_batch(idx_ref_Cd,:)/Lx,'o-');
loglog(tau_mean,rho0*hgradM_refkap_cntr_batch(idx_ref_Cd,:)/Lx,'o-','Color',colororder(6,:),'LineWidth',1.5);
loglog(tau_mean,rho0*hgradM_refflow_cntr_batch(idx_ref_Cd,:)/Lx,'o-','Color',colororder(7,:),'LineWidth',1.5);
loglog([.1 .1],[5e-3 .5],'k:','LineWidth',2);
hold off;
axis([0.01 .5 5e-3 .5]);
set(gca,'XTick',tau_ticks);
set(gca,'YTick',tau_ticks);
set(gca,'FontSize',fontsize);
xlabel('Wind stress maximum (N/m$^2$)','interpreter','latex');
%   'Eddy IFS, reconstructed from $\kappa$ and $\overline{\mathbf{u}}$', ...
leghandle = legend('Eddy IFS, diagnosed', ...
  'Reconstructed eddy IFS, $\kappa = \kappa_{\mathrm{ref}}$', ...
  'Reconstructed eddy IFS, $\overline{\mathbf{u}} = \overline{\mathbf{u}}_{\mathrm{ref}}$');
set(leghandle,'interpreter','latex','Location','NorthWest');
grid on;

%%% Scatter plot options
markersize = 14;
markershapes = {'>','*','<','o','v','d','^','s'};
defaultcolororder = get(gca,'ColorOrder');
colororder = zeros(8,3);
colororder(1:idx_ref_Cd-1,:) = defaultcolororder(1:idx_ref_Cd-1,:);
colororder(idx_ref_Cd,:) = [0 0 0];
colororder(idx_ref_Cd+1:N_Cd,:) = defaultcolororder(idx_ref_Cd:N_Cd-1,:);

lrmse = sqrt(mean((log10(-rho0*hgradM_eddy_cntr_batch(:)/Lx)-log10(rho0*hgradM_refkap_cntr_batch(:)/Lx)).^2));
subplot('Position',axpos(3,:));
for n_Cd = 1:N_Cd
  scatter(-rho0*hgradM_eddy_cntr_batch(n_Cd,:)/Lx,rho0*hgradM_refkap_cntr_batch(n_Cd,:)/Lx,markersize,'filled',markershapes{n_Cd},'MarkerEdgeColor',colororder(n_Cd,:),'MarkerFaceColor',colororder(n_Cd,:));
  if (n_Cd == 1)
    hold on;
  end
end
plot([1e-3 0.3],[1e-3 0.3],'k--');
hold off;
set(gca,'FontSize',fontsize);
xlabel('Eddy IFS, diagnosed (N/m$^2$)','interpreter','latex');
ylabel('Reconstructed eddy IFS, $\kappa = \kappa_{\mathrm{ref}}$ (N/m$^2$)','interpreter','latex');
box on;
grid on;
set(gca,'XTick',tau_ticks);
set(gca,'YTick',tau_ticks);
set(gca,'XScale','log');
set(gca,'YScale','log');
axis([4e-3 0.3 4e-3 0.3]);
text(0.07,0.005,['LRMSE = ',num2str(lrmse,'%.3f')],'FontSize',fontsize);

legstr = {};
for n_Cd=1:N_Cd
  legstr = {legstr{:},['Cd=',num2str(quad_drag(n_Cd)*1000),'$\times$ 10$^{-3}$']}; 
end
leghandle = legend(legstr,'Location','NorthWest');
set(leghandle,'interpreter','latex');

lrmse = sqrt(mean((log10(-rho0*hgradM_eddy_cntr_batch(:)/Lx)-log10(rho0*hgradM_refflow_cntr_batch(:)/Lx)).^2));
subplot('Position',axpos(4,:));
for n_Cd = 1:N_Cd
  scatter(-rho0*hgradM_eddy_cntr_batch(n_Cd,:)/Lx,rho0*hgradM_refflow_cntr_batch(n_Cd,:)/Lx,markersize,'filled',markershapes{n_Cd},'MarkerEdgeColor',colororder(n_Cd,:),'MarkerFaceColor',colororder(n_Cd,:));
  if (n_Cd == 1)
    hold on;
  end
end
plot([0.001 0.3],[0.001 0.3],'k--');
hold off;
set(gca,'FontSize',fontsize);
xlabel('Eddy IFS, diagnosed (N/m$^2$)','interpreter','latex');
ylabel('Reconstructed eddy IFS, $\overline{\mathbf{u}} = \overline{\mathbf{u}}_{\mathrm{ref}}$ (N/m$^2$)','interpreter','latex');
box on;
grid on;
set(gca,'XTick',tau_ticks);
set(gca,'YTick',tau_ticks);
set(gca,'XScale','log');
set(gca,'YScale','log');
axis([4e-3 0.3 4e-3 0.3]);
text(0.07,0.005,['LRMSE = ',num2str(lrmse,'%.3f')],'FontSize',fontsize);


%%% Add axis labels
for cntr = 1:size(axpos,1)
  annotation('textbox',[axpos(cntr,1)-0.07 axpos(cntr,2)-0.05 lab_size],'String',axlabels{cntr},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
end

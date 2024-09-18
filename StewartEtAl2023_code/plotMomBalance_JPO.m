%%%
%%% plotMomBalance_JPO.m
%%%
%%% Plots the model momentum balance for our JPO paper.
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

%%% Averaging period
tmin = 0.5*t1year;
tmax = 30.5*t1year;

%%% Remove redundant time dimension in wind stress vector
taux = squeeze(mean(taux,1));
tauy = squeeze(mean(tauy,1));

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

%%% Average u-momentum diagnostics  
UMom_PVadvection = do_avg(dirpath,OUTN_UMOM_Q,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
UMom_KEgradient = do_avg(dirpath,OUTN_UMOM_GRADKE,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
UMom_dhdt = do_avg(dirpath,OUTN_UMOM_DHDT,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
UMom_hypervisc = do_avg(dirpath,OUTN_UMOM_A4,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
UMom_quadBotDrag = do_avg(dirpath,OUTN_UMOM_CDBOT,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
UMom_windStress = do_avg(dirpath,OUTN_UMOM_WIND,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
UMom_advection = UMom_PVadvection + UMom_KEgradient + UMom_dhdt;

%%% Average v-momentum diagnostics  
VMom_PVadvection = do_avg(dirpath,OUTN_VMOM_Q,Nx,Ny,Nlay,n0_avg_hv,N_avg_hv,dt_avg_hv,tmin,tmax,startTime);
VMom_KEgradient = do_avg(dirpath,OUTN_VMOM_GRADKE,Nx,Ny,Nlay,n0_avg_hv,N_avg_hv,dt_avg_hv,tmin,tmax,startTime);
VMom_dhdt = do_avg(dirpath,OUTN_VMOM_DHDT,Nx,Ny,Nlay,n0_avg_hv,N_avg_hv,dt_avg_hv,tmin,tmax,startTime);
VMom_hypervisc = do_avg(dirpath,OUTN_VMOM_A4,Nx,Ny,Nlay,n0_avg_hv,N_avg_hv,dt_avg_hv,tmin,tmax,startTime);
VMom_quadBotDrag = do_avg(dirpath,OUTN_VMOM_CDBOT,Nx,Ny,Nlay,n0_avg_hv,N_avg_hv,dt_avg_hv,tmin,tmax,startTime);
VMom_windStress = do_avg(dirpath,OUTN_VMOM_WIND,Nx,Ny,Nlay,n0_avg_hv,N_avg_hv,dt_avg_hv,tmin,tmax,startTime);
VMom_advection = VMom_PVadvection + VMom_KEgradient + VMom_dhdt;

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
  
%%% Surface Montgomery potential
MM_surf = MM_tavg(:,:,1);
MM_thresh = median(MM_surf(:));
MM_grid = gg(1)*(-.25:0.005:.3);
MM_surf_q = MM_surf;
MM_surf_q(:,2:Ny,:) = 0.25*(MM_surf(1:Nx,1:Ny-1,:)+MM_surf([Nx 1:Nx-1],1:Ny-1,:)+MM_surf(1:Nx,2:Ny,:)+MM_surf([Nx 1:Nx-1],2:Ny,:));
MM_surf_q(:,1,:) = 0.5*(MM_surf(1:Nx,1,:)+MM_surf([Nx 1:Nx-1],1,:));

%%% Momentum curl budget quantities  
curl_quadBotDrag = (VMom_quadBotDrag(1:Nx,:,:)-VMom_quadBotDrag([Nx 1:Nx-1],:,:))/dx - (UMom_quadBotDrag(:,1:Ny,:)-UMom_quadBotDrag(:,[Ny 1:Ny-1],:))/dy;
curl_quadBotDrag(:,1,:) = 0;
curl_quadBotDrag(:,2,:) = (VMom_quadBotDrag(1:Nx,2,:)-VMom_quadBotDrag([Nx 1:Nx-1],2,:))/dx - (UMom_quadBotDrag(:,2,:))/dy;
curl_quadBotDrag = curl_quadBotDrag(:,:,Nlay);
curl_hypervisc = (VMom_hypervisc(1:Nx,:,:)-VMom_hypervisc([Nx 1:Nx-1],:,:))/dx - (UMom_hypervisc(:,1:Ny,:)-UMom_hypervisc(:,[Ny 1:Ny-1],:))/dy;
curl_hypervisc(:,1,:) = 0;
curl_hypervisc(:,2,:) = (VMom_hypervisc(1:Nx,2,:)-VMom_hypervisc([Nx 1:Nx-1],2,:))/dx - (UMom_hypervisc(:,2,:))/dy;
curl_advection = (VMom_advection(1:Nx,:,:)-VMom_advection([Nx 1:Nx-1],:,:))/dx - (UMom_advection(:,1:Ny,:)-UMom_advection(:,[Ny 1:Ny-1],:))/dy;
curl_advection(:,1,:) = 0;
curl_advection(:,2,:) = (VMom_advection(1:Nx,2,:)-VMom_advection([Nx 1:Nx-1],2,:))/dx - (UMom_advection(:,2,:))/dy;
curl_tau = (tauy(1:Nx,:)-tauy([Nx 1:Nx-1],:))/dx - (taux(:,1:Ny)-taux(:,[Ny 1:Ny-1]))/dy;
curl_tau(:,1) = 0;
curl_hgradM_eddy = (hdMdy_eddy(1:Nx,:,:)-hdMdy_eddy([Nx 1:Nx-1],:,:))/dx - (hdMdx_eddy(:,1:Ny,:)-hdMdx_eddy(:,[Ny 1:Ny-1],:))/dy;
curl_hgradM_eddy(:,1,:) = 0;
curl_hgradM_eddy(:,2,:) = (hdMdy_eddy(1:Nx,2,:)-hdMdy_eddy([Nx 1:Nx-1],2,:))/dx - (hdMdx_eddy(:,2,:))/dy;
curl_hgradM_mean = (hdMdy_mean(1:Nx,:,:)-hdMdy_mean([Nx 1:Nx-1],:,:))/dx - (hdMdx_mean(:,1:Ny,:)-hdMdx_mean(:,[Ny 1:Ny-1],:))/dy;
curl_hgradM_mean(:,1,:) = 0;
curl_eddyforce = (eddyforce_y(1:Nx,:,:)-eddyforce_y([Nx 1:Nx-1],:,:))/dx - (eddyforce_x(:,1:Ny,:)-eddyforce_x(:,[Ny 1:Ny-1],:))/dy;
curl_eddyforce(:,1,:) = 0;
curl_eddyforce(:,2,:) = (eddyforce_y(1:Nx,2,:)-eddyforce_y([Nx 1:Nx-1],2,:))/dx - (eddyforce_x(:,2,:))/dy;
curl_hgradM_tavg = (hdMdy_tavg(1:Nx,:,:)-hdMdy_tavg([Nx 1:Nx-1],:,:))/dx - (hdMdx_tavg(:,1:Ny,:)-hdMdx_tavg(:,[Ny 1:Ny-1],:))/dy;
curl_hgradM_tavg(:,1,:) = 0;
curl_hgradM_tavg(:,2,:) = (hdMdy_tavg(1:Nx,2,:)-hdMdy_tavg([Nx 1:Nx-1],2,:))/dx - (hdMdx_tavg(:,2,:))/dy;
BPT = sum(curl_hgradM_tavg,3);

%%% For reconstruction of eddy fluxes
f0 = mean(mean(2*Omega_z));
curl_twa = (vv_twa(1:Nx,:,:)-vv_twa([Nx 1:Nx-1],:,:))/dx - (uu_twa(:,1:Ny,:)-uu_twa(:,[Ny 1:Ny-1],:))/dy;
curl_twa(:,1,:) = 0;
curl_twa(:,2,:) = (vv_twa(1:Nx,2,:)-vv_twa([Nx 1:Nx-1],2,:))/dx - (uu_twa(:,2,:))/dy;
diff_curl_twa = diff(curl_twa,1,3);
IPT_kap = kap_map.*(f0^2/gg(2)).*diff_curl_twa;  

%%% Along-streamline averages
tau_cntr = 0*MM_grid;
tfs_cntr = 0*MM_grid;
drag_cntr = 0*MM_grid;
visc_cntr_bt = 0*MM_grid;
visc_cntr_bc = 0*MM_grid;
adv_cntr_bt = 0*MM_grid;
adv_cntr_bc = 0*MM_grid;
hgradM_tavg_cntr_bc = 0*MM_grid;
hgradM_eddy_cntr_bc = 0*MM_grid;
hgradM_kap_cntr = 0*MM_grid;
eddyforce_cntr_bt = 0*MM_grid;
eddyforce_cntr_bc = 0*MM_grid;
for m=1:length(MM_grid)
  
  MM_cntr = MM_grid(m);
   
  tau_cntr(m) = -sum(sum(curl_tau(MM_surf_q<MM_cntr)*dx*dy));
  tfs_cntr(m) = -sum(sum(BPT(MM_surf_q<MM_cntr)*dx*dy));    
  drag_cntr(m) = -sum(sum(curl_quadBotDrag(MM_surf_q<MM_cntr)*dx*dy));
  
  curl_advection_bt = sum(curl_advection,3);
  adv_cntr_bt(m) = -sum(sum(curl_advection_bt(MM_surf_q<MM_cntr)*dx*dy));  
  curl_advection_bc = curl_advection(:,:,1);
  adv_cntr_bc(m) = -sum(sum(curl_advection_bc(MM_surf_q<MM_cntr)*dx*dy));  
  
  curl_hgradM_tavg_bc = curl_hgradM_tavg(:,:,1);
  hgradM_tavg_cntr_bc(m) = -sum(sum(curl_hgradM_tavg_bc(MM_surf_q<MM_cntr)*dx*dy));  
  curl_hgradM_eddy_bc = curl_hgradM_eddy(:,:,1);
  hgradM_eddy_cntr_bc(m) = -sum(sum(curl_hgradM_eddy_bc(MM_surf_q<MM_cntr)*dx*dy));  
     
  curl_eddyforce_bt = sum(curl_eddyforce,3);
  eddyforce_cntr_bt(m) = -sum(sum(curl_eddyforce_bt(MM_surf_q<MM_cntr)*dx*dy));  
  curl_eddyforce_bc = curl_eddyforce(:,:,1);
  eddyforce_cntr_bc(m) = -sum(sum(curl_eddyforce_bc(MM_surf_q<MM_cntr)*dx*dy));  
  
  curl_visc_bt = sum(curl_hypervisc,3);
  visc_cntr_bt(m) = -sum(sum(curl_visc_bt(MM_surf_q<MM_cntr)*dx*dy));  
  curl_visc_bc = curl_hypervisc(:,:,1);
  visc_cntr_bc(m) = -sum(sum(curl_visc_bc(MM_surf_q<MM_cntr)*dx*dy));  
    
  hgradM_kap_cntr(m) = sum(sum(IPT_kap(MM_surf_q<MM_cntr)*dx*dy));
    
end



%%% Plotting options
fontsize = 14;
axpos = zeros(5,4);
axpos(1,:) = [0.21 0.72 .65 .25];
axpos(2,:) = [0.07 0.44 .4 .2];
axpos(4,:) = [0.07 0.18 .4 .2];
axpos(3,:) = [0.57 0.44 .4 .2];
axpos(5,:) = [0.57 0.18 .4 .2];
legpos1 = [0.17 0.01 0.2 0.1];
legpos2 = [0.67 0.02 0.2 0.1];
axlabels = {'(a)','(b)','(c)','(d)','(e)'};
lab_size = [0.05 0.03];
tau_ticks = [0.01 0.017 0.03 0.05 0.1 0.17 0.3];
rho0 = 1000;
linewidth = 2;

figure(102);
clf;
set(gcf,'Position',[382   306   900   900]);





subplot('Position',axpos(1,:));
pcolor(XX_h/m1km,YY_h/m1km,MM_surf/gg(1));
shading interp;
hold on;
[C,h]=contour(XX_h/m1km,YY_h/m1km,-hhb,[500 1500 2500 3250 3500 3750],'EdgeColor',[.5 .5 .5],'LineWidth',1);
contour(XX_h/m1km,YY_h/m1km,MM_surf/gg(1),[MM_thresh MM_thresh]/gg(1),'EdgeColor','k');
hold off;
xlabel('Longitude $x$ (km)','interpreter','latex');
ylabel('Latitude $y$ (km)','interpreter','latex');
colorbar;
caxis([-.25 .3]);
colormap(gca,haxby(22));
set(gca,'FontSize',fontsize);
title('Mean sea surface height (m)');

colororder = get(gca,'ColorOrder');

subplot('Position',axpos(2,:));
plot(yy_h/m1km,mean(taux,1)*rho0,'Color',colororder(1,:),'LineWidth',2);
hold on;
plot(yy_h/m1km,mean(sum(-hdMdx_tavg,3),1)*rho0,'Color',colororder(2,:),'LineWidth',2);
plot(yy_h/m1km,mean(sum(UMom_advection,3),1)*rho0,'Color',colororder(3,:),'LineWidth',2);
plot(yy_h/m1km,mean(sum(-eddyforce_x,3),1)*rho0,'Color',colororder(3,:),'LineStyle','--','LineWidth',2);
plot(yy_h/m1km,mean(sum(UMom_quadBotDrag,3),1)*rho0,'Color',colororder(4,:),'LineWidth',2);
% plot(yy_h/m1km,mean(sum(UMom_hypervisc,3),1)*rho0,'Color',colororder(5,:),'LineWidth',2);
plot(yy_h/m1km,mean(taux+sum(-hdMdx_tavg+UMom_quadBotDrag+UMom_advection+UMom_hypervisc,3),1)*rho0,'Color',[.7 .7 .7],'LineWidth',2);
hold off;
set(gca,'FontSize',fontsize);
ylabel('N/m$^2$','interpreter','latex');
xlabel('Latitude (km)','interpreter','latex');
set(gca,'YLim',[-.13 .13]);
title('Depth-integrated momentum balance');
grid on;

subplot('Position',axpos(3,:));
plot(yy_h/m1km,mean(taux,1)*rho0,'Color',colororder(1,:),'LineWidth',2);
hold on;
plot(yy_h/m1km,mean(-(hdMdx_tavg(:,:,1)-hdMdx_eddy(:,:,1)),1)*rho0,'Color',colororder(2,:),'LineWidth',2);
plot(yy_h/m1km,mean(-hdMdx_eddy(:,:,1),1)*rho0,'Color',colororder(2,:),'LineStyle','--','LineWidth',2);
plot(yy_h/m1km,mean(UMom_advection(:,:,1)+eddyforce_x(:,:,1),1)*rho0,'Color',colororder(3,:),'LineWidth',2);
plot(yy_h/m1km,mean(-eddyforce_x(:,:,1),1)*rho0,'Color',colororder(3,:),'LineStyle','--','LineWidth',2);
% plot(yy_h/m1km,mean(UMom_hypervisc(:,:,1),1)*rho0,'Color',colororder(5,:),'LineWidth',2);
plot(yy_h/m1km,mean(taux-hdMdx_tavg(:,:,1)+UMom_advection(:,:,1)+UMom_hypervisc(:,:,1),1)*rho0,'Color',[.7 .7 .7],'LineWidth',2);
hold off;
set(gca,'FontSize',fontsize);
xlabel('Latitude (km)','interpreter','latex');
set(gca,'YLim',[-.13 .13]);
title('Upper layer momentum balance');
grid on;

subplot('Position',axpos(4,:));
plot(MM_grid/gg(1),tau_cntr*rho0/Lx,'Color',colororder(1,:),'LineWidth',2);
hold on;
plot(MM_grid/gg(1),-tfs_cntr*rho0/Lx,'Color',colororder(2,:),'LineWidth',2);
plot(MM_grid/gg(1),adv_cntr_bt*rho0/Lx,'Color',colororder(3,:),'LineWidth',2);
plot(MM_grid/gg(1),-eddyforce_cntr_bt*rho0/Lx,'Color',colororder(3,:),'LineStyle','--','LineWidth',2);
plot(MM_grid/gg(1),drag_cntr*rho0/Lx,'Color',colororder(4,:),'LineWidth',2);
% plot(MM_grid/gg(1),visc_cntr_bt*rho0/Lx,'Color',colororder(5,:),'LineWidth',2);
plot(MM_grid/gg(1),(tau_cntr-tfs_cntr+adv_cntr_bt+drag_cntr+visc_cntr_bt)*rho0/Lx,'Color',[.7 .7 .7],'LineWidth',2);
plot([MM_thresh, MM_thresh]/gg(1),[-.13 .13],'k--');
hold off;
set(gca,'FontSize',fontsize);
set(gca,'YLim',[-.13 .13])
ylabel('N/m$^2$','interpreter','latex');
set(gca,'XLim',[-.25 .3]);
xlabel('Mean sea surface height (m)','interpreter','latex');
grid on;

leghandle = legend('Wind stress','Topographic form stress','Advection (mean)','Advection (eddy)','Bottom friction','Sum');
set(leghandle,'Position',legpos1);

subplot('Position',axpos(5,:));
plot(MM_grid/gg(1),tau_cntr*rho0/Lx,'Color',colororder(1,:),'LineWidth',2);
hold on;
plot(MM_grid/gg(1),-(hgradM_tavg_cntr_bc-hgradM_eddy_cntr_bc)*rho0/Lx,'Color',colororder(2,:),'LineWidth',2);
plot(MM_grid/gg(1),-hgradM_eddy_cntr_bc*rho0/Lx,'Color',colororder(2,:),'LineStyle','--','LineWidth',2);
plot(MM_grid/gg(1),-hgradM_kap_cntr*rho0/Lx,'Color',colororder(2,:),'LineStyle',':','LineWidth',2);
plot(MM_grid/gg(1),(adv_cntr_bc+eddyforce_cntr_bc)*rho0/Lx,'Color',colororder(3,:),'LineWidth',2);
plot(MM_grid/gg(1),-eddyforce_cntr_bc*rho0/Lx,'Color',colororder(3,:),'LineStyle','--','LineWidth',2);
% plot(MM_grid/gg(1),visc_cntr_bc*rho0/Lx,'Color',colororder(5,:),'LineWidth',2);
plot(MM_grid/gg(1),(tau_cntr-hgradM_tavg_cntr_bc+adv_cntr_bc+visc_cntr_bc)*rho0/Lx,'Color',[.7 .7 .7],'LineWidth',1);
plot([MM_thresh, MM_thresh]/gg(1),[-.13 .13],'k--');
hold off;
set(gca,'FontSize',fontsize);
xlabel('Mean sea surface height (m)','interpreter','latex');
set(gca,'YLim',[-.13 .13]);
set(gca,'XLim',[-.25 .3]);
grid on;

leghandle = legend('Wind stress','Interfacial form stress (mean)','Interfacial form stress (eddy)','Interfacial form stress (reconstructed)','Advection (mean)','Advection (eddy)','Sum');
set(leghandle,'Position',legpos2);

%%% Add axis labels
for cntr = 1:size(axpos,1)
  annotation('textbox',[axpos(cntr,1)-0.05 axpos(cntr,2)-0.04 lab_size],'String',axlabels{cntr},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
end
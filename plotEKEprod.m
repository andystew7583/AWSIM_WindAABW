%%%
%%% plotEKEprod.m
%%%
%%% Calculates and plots EKE production rates.
%%%

%%% Run to load
% run_name = 'ACC_AABW_Ny128_Nlay2_tauM0.1_tauP0_tauF0_wDiaM0_wDiaP0_wDiaF0_Cd2.000e-03_rb0.000e+00_diags';
run_name = 'ACC_AABW_ML_doubleMOC_hires';
% run_name = 'ACC_AABW_ML_doubleMOC';

%%% Load parameters   
local_home_dir = '/Volumes/Kilchoman/UCLA/Projects/AWSIM/runs';
% local_home_dir = '/Volumes/Kilchoman/UCLA/Projects/AWSIM_WindAABW/runs';
prod_dir = fullfile(local_home_dir,'AWSIM_WindAABW_products');
loadParams;
dirpath = fullfile(local_home_dir,run_name);
gtild = reshape(cumsum(gg),[Nlay 1 1]);

%%% Load time-mean products
load([run_name,'_tavg.mat']);

%%% Thickness-weighted average velocities
uu_twa = hu_tavg ./ hh_w_tavg;
vv_twa = hv_tavg ./ hh_s_tavg;
usq_twa = husq_tavg ./ hh_w_tavg;
vsq_twa = hvsq_tavg ./ hh_s_tavg;
uu_bol = uu_twa - uu_tavg;
vv_bol = vv_twa - vv_tavg;

%%% Compute eddy pressure flux terms
hphi_mean = hh_tavg.*phi_tavg;
hphi_eddy = hphi_tavg - hphi_mean;
dx_hphi_eddy = (hphi_eddy(1:Nx,:,:)-hphi_eddy([Nx 1:Nx-1],:,:))/dx;
dy_hphi_eddy = (hphi_eddy(:,1:Ny,:)-hphi_eddy(:,[Ny 1:Ny-1],:))/dy;

%%% Compute eddy form stress terms
pdedx_mean = 0.5*(phi_int_tavg(1:Nx,:,:)+phi_int_tavg([Nx 1:Nx-1],:,:)).*(eta_tavg(1:Nx,:,:)-eta_tavg([Nx 1:Nx-1],:,:))/dx;
pdedx_eddy = pdedx_tavg - pdedx_mean;
pdedy_mean = 0.5*(phi_int_tavg(:,1:Ny,:)+phi_int_tavg(:,[Ny 1:Ny-1],:)).*(eta_tavg(:,1:Ny,:)-eta_tavg(:,[Ny 1:Ny-1],:))/dy;
pdedy_eddy = pdedy_tavg - pdedy_mean;

%%% Compute eddy momentum flux terms
husq_mean = hh_w_tavg.*uu_twa.^2; %hh_tavg.*0.5.*(uu_twa(1:Nx,:,:).^2+uu_twa([2:Nx 1],:,:).^2);
hvsq_mean = hh_s_tavg.*vv_twa.^2; %hh_tavg.*0.5.*(vv_twa(:,1:Ny,:).^2+vv_twa(:,[2:Ny 1],:).^2);
huv_mean = 0.5.*(vv_twa(1:Nx,:,:)+vv_twa([Nx 1:Nx-1],:,:)) ...
                .* 0.5.*(uu_twa(:,1:Ny,:)+uu_twa(:,[Ny 1:Ny-1],:)) ...
                .* 0.25.*(hh_tavg(1:Nx,1:Ny,:)+hh_tavg(1:Nx,[Ny 1:Ny-1],:)+hh_tavg([Nx 1:Nx-1],1:Ny,:)+hh_tavg([Nx 1:Nx-1],[Ny 1:Ny-1],:));                  
husq_eddy = husq_tavg - husq_mean;
hvsq_eddy = hvsq_tavg - hvsq_mean;
huv_eddy = huv_tavg - huv_mean;
dx_husq_eddy = (husq_eddy([2:Nx 1],:,:)-husq_eddy([Nx 1:Nx-1],:,:))/(2*dx);
dy_hvsq_eddy = (hvsq_eddy(:,[2:Ny 1],:)-hvsq_eddy(:,[Ny 1:Ny-1],:))/(2*dy);
dx_huv_eddy = (huv_eddy([2:Nx 1],:,:)-huv_eddy(1:Nx,:,:))/dx;
dy_huv_eddy = (huv_eddy(:,[2:Ny 1],:)-huv_eddy(:,1:Ny,:))/dy;
% dx_husq_eddy = (husq_eddy(1:Nx,:,:)-husq_eddy([Nx 1:Nx-1],:,:))/dx;
% dy_hvsq_eddy = (hvsq_eddy(:,1:Ny,:)-hvsq_eddy(:,[Ny 1:Ny-1],:))/dy;
% dx_huv_eddy = (huv_eddy([2:Nx 1],:,:)-huv_eddy(1:Nx,:,:))/dx;
% dy_huv_eddy = (huv_eddy(:,[2:Ny 1],:)-huv_eddy(:,1:Ny,:))/dy;

%%% Total eddy pressure forcing
hdMdx_mean = hh_w_tavg.*(MM_tavg(1:Nx,:,:)-MM_tavg([Nx 1:Nx-1],:,:))/dx;
hdMdy_mean = hh_s_tavg.*(MM_tavg(:,1:Ny,:)-MM_tavg(:,[Ny 1:Ny-1],:))/dy;  
hdMdx_eddy = hdMdx_tavg - hdMdx_mean;
hdMdy_eddy = hdMdy_tavg - hdMdy_mean;

%%% Eddy-Montgomery potential correlations
edMdx_p_mean = eta_w_tavg(:,:,1:Nlay).*(MM_tavg(1:Nx,:,:)-MM_tavg([Nx 1:Nx-1],:,:))/dx;
edMdy_p_mean = eta_s_tavg(:,:,1:Nlay).*(MM_tavg(:,1:Ny,:)-MM_tavg(:,[Ny 1:Ny-1],:))/dy;     
edMdx_m_mean = eta_w_tavg(:,:,2:Nlay+1).*(MM_tavg(1:Nx,:,:)-MM_tavg([Nx 1:Nx-1],:,:))/dx;
edMdy_m_mean = eta_s_tavg(:,:,2:Nlay+1).*(MM_tavg(:,1:Ny,:)-MM_tavg(:,[Ny 1:Ny-1],:))/dy; 
edMdx_p_eddy = edMdx_p_tavg - edMdx_p_mean;
edMdy_p_eddy = edMdy_p_tavg - edMdy_p_mean;
edMdx_m_eddy = edMdx_m_tavg - edMdx_m_mean;
edMdy_m_eddy = edMdy_m_tavg - edMdy_m_mean;

PEtoEKE_v1 = uu_twa.*hdMdx_eddy+vv_twa.*hdMdy_eddy;
PEtoEKE_v2 = uu_twa .* ((hphi_eddy(1:Nx,:,:)-hphi_eddy([Nx 1:Nx-1],:,:))/dx + diff(pdedx_eddy,1,3)) ...
           + vv_twa .* ((hphi_eddy(:,1:Ny,:)-hphi_eddy(:,[Ny 1:Ny-1],:))/dy + diff(pdedy_eddy,1,3));
PEtoEKE_v3 = - hphi_eddy .* (uu_twa([2:Nx 1],:,:)-uu_twa([1:Nx],:,:))/dx + uu_twa.*diff(pdedx_eddy,1,3) ...
             - hphi_eddy .* (vv_twa(:,[2:Ny 1],:)-vv_twa(:,[1:Ny],:))/dy + vv_twa.*diff(pdedy_eddy,1,3);  
PEtoEKE_v4_iso = uu_twa .* (hphi_eddy(1:Nx,:,:)-hphi_eddy([Nx 1:Nx-1],:,:))/dx + vv_twa .* (hphi_eddy(:,1:Ny,:)-hphi_eddy(:,[Ny 1:Ny-1],:))/dy;
PEtoEKE_v4_dia = - pdedx_eddy(:,:,2:Nlay) .* diff(uu_twa,1,3) - pdedy_eddy(:,:,2:Nlay) .* diff(vv_twa,1,3);
PEtoEKE_v5_iso = - hphi_eddy .* (uu_twa([2:Nx 1],:,:)-uu_twa([1:Nx],:,:))/dx ...
                 - hphi_eddy .* (vv_twa(:,[2:Ny 1],:)-vv_twa(:,[1:Ny],:))/dy;
PEtoEKE_v5_dia = - pdedx_eddy(:,:,2:Nlay) .* diff(uu_twa,1,3) - pdedy_eddy(:,:,2:Nlay) .* diff(vv_twa,1,3);
PEtoEKE_v6 = uu_twa(:,:,2:Nlay).*edMdx_p_eddy(:,:,2:Nlay) + vv_twa(:,:,2:Nlay).*edMdy_p_eddy(:,:,2:Nlay) ...
           - uu_twa(:,:,1:Nlay-1).*edMdx_m_eddy(:,:,1:Nlay-1) - vv_twa(:,:,1:Nlay-1).*edMdy_m_eddy(:,:,1:Nlay-1);
MKEtoEKE = uu_twa.*dx_husq_eddy + uu_twa.*dy_huv_eddy + vv_twa.*dy_hvsq_eddy + vv_twa.*dx_huv_eddy;

%%% Plot eddy forcing
fignum = 1;
figpos = [272   805   889   372];
axpos = [0.08 0.11 0.83 0.79];
fontsize = 14;
yaxlim = [-1 1]*1e-4;
outfdir = '~/Desktop/IsopycnalEnergetics_23Mar2022';

figure(fignum);
fignum = fignum+1;
pcolor(XX_h/1000,YY_h/1000,sum(PEtoEKE_v1,3));
shading interp
colorbar;
colormap redblue;
caxis(yaxlim);
title(['Baroclinic production, Aiki, depth-integrated (m^3/s^3), total = ',num2str(sum(1000*PEtoEKE_v1(:)*dx*dy),'%.2e'),' J/s']);
set(gcf,'Position',figpos);
set(gca,'Position',axpos);
set(gca,'FontSize',fontsize);
xlabel('x (km)');
ylabel('y (km)');
% print('-dpng',fullfile(outfdir,'Aiki_DepthIntegrated.png'));


figure(fignum);
fignum = fignum+1;
pcolor(XX_h/1000,YY_h/1000,sum(PEtoEKE_v2,3));
shading interp
colorbar;
colormap redblue;
caxis(yaxlim);
title('Baroclinic production, v2');


figure(fignum);
fignum = fignum+1;
pcolor(XX_h/1000,YY_h/1000,sum(PEtoEKE_v3,3));
shading interp
colorbar;
colormap redblue;
caxis(yaxlim);
title('Baroclinic production, v3');

figure(fignum);
fignum = fignum+1;
pcolor(XX_h/1000,YY_h/1000,sum(PEtoEKE_v4_iso,3)+sum(PEtoEKE_v4_dia,3));
shading interp
colorbar;
colormap redblue;
caxis(yaxlim);
title('Baroclinic production, v4');

figure(fignum);
fignum = fignum+1;
pcolor(XX_h/1000,YY_h/1000,sum(PEtoEKE_v5_iso,3)+sum(PEtoEKE_v5_dia,3));
shading interp
colorbar;
colormap redblue;
caxis(yaxlim);
title('Baroclinic production, v5');

figure(fignum);
fignum = fignum+1;
pcolor(XX_h/1000,YY_h/1000,sum(PEtoEKE_v6,3));
shading interp
colorbar;
colormap redblue;
caxis(yaxlim);
title(['Baroclinic production, Stewart, depth-integrated (m^3/s^3), total = ',num2str(sum(1000*PEtoEKE_v6(:)*dx*dy),'%.2e'),' J/s']);
set(gcf,'Position',figpos);
set(gca,'Position',axpos);
set(gca,'FontSize',fontsize);
xlabel('x (km)');
ylabel('y (km)');
% print('-dpng',fullfile(outfdir,'Stewart_DepthIntegrated.png'));

figure(fignum);
fignum = fignum+1;
pcolor(XX_h/1000,YY_h/1000,PEtoEKE_v1(:,:,1));
shading interp
colorbar;
colormap redblue;
caxis(yaxlim);
title(['Baroclinic production, Aiki, layer 1 (m^3/s^3), total = ',num2str(sum(sum(1000*PEtoEKE_v1(:,:,1)*dx*dy)),'%.2e'),' J/s']);
set(gcf,'Position',figpos);
set(gca,'Position',axpos);
set(gca,'FontSize',fontsize);
xlabel('x (km)');
ylabel('y (km)');
% print('-dpng',fullfile(outfdir,'Aiki_Layer1.png'));

figure(fignum);
fignum = fignum+1;
pcolor(XX_h/1000,YY_h/1000,PEtoEKE_v1(:,:,2));
shading interp
colorbar;
colormap redblue;
caxis(yaxlim);
title(['Baroclinic production, Aiki, layer 2 (m^3/s^3), total = ',num2str(sum(sum(1000*PEtoEKE_v1(:,:,2)*dx*dy)),'%.2e'),' J/s']);
set(gcf,'Position',figpos);
set(gca,'Position',axpos);
set(gca,'FontSize',fontsize);
xlabel('x (km)');
ylabel('y (km)');
% print('-dpng',fullfile(outfdir,'Aiki_Layer2.png'));

figure(fignum);
fignum = fignum+1;
pcolor(XX_h/1000,YY_h/1000,PEtoEKE_v1(:,:,3));
shading interp
colorbar;
colormap redblue;
caxis(yaxlim);
title(['Baroclinic production, Aiki, layer 3 (m^3/s^3), total = ',num2str(sum(sum(1000*PEtoEKE_v1(:,:,3)*dx*dy)),'%.2e'),' J/s']);
set(gcf,'Position',figpos);
set(gca,'Position',axpos);
set(gca,'FontSize',fontsize);
xlabel('x (km)');
ylabel('y (km)');
% print('-dpng',fullfile(outfdir,'Aiki_Layer3.png'));

figure(fignum);
fignum = fignum+1;
pcolor(XX_h/1000,YY_h/1000,PEtoEKE_v6(:,:,1));
shading interp
colorbar;
colormap redblue;
caxis(yaxlim);
title(['Baroclinic production, Stewart, isopycnal 1 (m^3/s^3), total = ',num2str(sum(sum(1000*PEtoEKE_v6(:,:,1)*dx*dy)),'%.2e'),' J/s']);
set(gcf,'Position',figpos);
set(gca,'Position',axpos);
set(gca,'FontSize',fontsize);
xlabel('x (km)');
ylabel('y (km)');
% print('-dpng',fullfile(outfdir,'Stewart_Isopycnal1.png'));

figure(fignum);
fignum = fignum+1;
pcolor(XX_h/1000,YY_h/1000,PEtoEKE_v6(:,:,2));
shading interp
colorbar;
colormap redblue;
caxis(yaxlim);
title(['Baroclinic production, Stewart, isopycnal 2 (m^3/s^3), total = ',num2str(sum(sum(1000*PEtoEKE_v6(:,:,2)*dx*dy)),'%.2e'),' J/s']);
set(gcf,'Position',figpos);
set(gca,'Position',axpos);
set(gca,'FontSize',fontsize);
xlabel('x (km)');
ylabel('y (km)');
% print('-dpng',fullfile(outfdir,'Stewart_Isopycnal2.png'));

figure(fignum);
fignum = fignum+1;
pcolor(XX_h/1000,YY_h/1000,sum(MKEtoEKE,3));
shading interp
colorbar;
colormap redblue;
caxis(yaxlim);
title('MKE to EKE');

figure(fignum);
fignum = fignum+1;
pcolor(XX_h/1000,YY_h/1000,MKEtoEKE(:,:,1));
shading interp
colorbar;
colormap redblue;
caxis(yaxlim);
title('MKE to EKE, layer 1');

figure(fignum);
fignum = fignum+1;
pcolor(XX_h/1000,YY_h/1000,MKEtoEKE(:,:,2));
shading interp
colorbar;
colormap redblue;
caxis(yaxlim);
title('MKE to EKE, layer 2');

figure(fignum);
fignum = fignum+1;
pcolor(XX_h/1000,YY_h/1000,MKEtoEKE(:,:,3));
shading interp
colorbar;
colormap redblue;
caxis(yaxlim);
title('MKE to EKE, layer 3');

figure(fignum);
fignum = fignum+1;
pcolor(XX_h/1000,YY_h/1000,-pdedt_tavg(:,:,2));
shading interp
colorbar;
colormap redblue;
caxis(yaxlim);
title(['Downward eddy energy flux, isopycnal 1 (m^3/s^3), total = ',num2str(sum(sum(1000*-pdedt_tavg(:,:,2)*dx*dy)),'%.2e'),' J/s']);
set(gcf,'Position',figpos);
set(gca,'Position',axpos);
set(gca,'FontSize',fontsize);
xlabel('x (km)');
ylabel('y (km)');
% print('-dpng',fullfile(outfdir,'EnergyFlux_Isopycnal1.png'));

figure(fignum);
fignum = fignum+1;
pcolor(XX_h/1000,YY_h/1000,-pdedt_tavg(:,:,3));
shading interp
colorbar;
colormap redblue;
caxis(yaxlim);
title(['Downward eddy energy flux, isopycnal 2 (m^3/s^3), total = ',num2str(sum(sum(1000*-pdedt_tavg(:,:,2)*dx*dy)),'%.2e'),' J/s']);
set(gcf,'Position',figpos);
set(gca,'Position',axpos);
set(gca,'FontSize',fontsize);
xlabel('x (km)');
ylabel('y (km)');
% print('-dpng',fullfile(outfdir,'EnergyFlux_Isopycnal2.png'));

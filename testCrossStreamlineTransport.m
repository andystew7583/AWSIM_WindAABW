%%%
%%% testCrossStreamlineTransport.m
%%%
%%% Test script for diagnostics of fluxes across streamlines
%%%

%%% Run to load
run_name = 'ACC_AABW_Ny128_Nlay2_tauM0.1_tauP0_tauF0_wDiaM0_wDiaP0_wDiaF0_Cd4.000e-03_rb0.000e+00_diags';
% run_name = 'ACC_AABW_ML_randWdia_randTau_white';

%%% Load parameters   
local_home_dir = '/Volumes/Kilchoman/UCLA/Projects/AWSIM_WindAABW/runs';
prod_dir = fullfile(local_home_dir,'AWSIM_WindAABW_products');
loadParams;
dirpath = fullfile(local_home_dir,run_name);
gtild = reshape(cumsum(gg),[Nlay 1 1]);

%%% Averaging period
tmin = 0.5*t1year;
tmax = 30.5*t1year;

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

%%% Momentum curl budget quantities
curl_tau = (tauy(1:Nx,:,:)-tauy([Nx 1:Nx-1],:,:))/dx - (taux(:,1:Ny,:)-taux(:,[Ny 1:Ny-1],:))/dy;
curl_tau(:,1) = 0;
curl_hgradM_eddy = (hdMdy_eddy(1:Nx,:,:)-hdMdy_eddy([Nx 1:Nx-1],:,:))/dx - (hdMdx_eddy(:,1:Ny,:)-hdMdx_eddy(:,[Ny 1:Ny-1],:))/dy;
curl_hgradM_eddy(:,1) = 0;
curl_hgradM_mean = (hdMdy_mean(1:Nx,:,:)-hdMdy_mean([Nx 1:Nx-1],:,:))/dx - (hdMdx_mean(:,1:Ny,:)-hdMdx_mean(:,[Ny 1:Ny-1],:))/dy;
curl_hgradM_mean(:,1) = 0;
curl_eddyforce = (eddyforce_y(1:Nx,:,:)-eddyforce_y([Nx 1:Nx-1],:,:))/dx - (eddyforce_x(:,1:Ny,:)-eddyforce_x(:,[Ny 1:Ny-1],:))/dy;
curl_eddyforce(:,1) = 0;
curl_twa = (vv_twa(1:Nx,:,:)-vv_twa([Nx 1:Nx-1],:,:))/dx - (uu_twa(:,1:Ny,:)-uu_twa(:,[Ny 1:Ny-1],:))/dy;
curl_twa(:,1,:) = 0;
diff_curl_twa = diff(curl_twa,1,3);



fignum = 1;

figure(fignum);
fignum = fignum + 1;
contourf(XX_h,YY_h,MM_tavg(:,:,1),30);
colorbar;

figure(fignum)
fignum = fignum + 1;
pcolor(XX_h,YY_h,sum(hdMdx_eddy.*uu_twa+hdMdy_eddy.*vv_twa,3));
shading interp;
colorbar;
colormap redblue;
caxis([-3 3]*1e-3);

figure(fignum)
fignum = fignum + 1;
pcolor(XX_h,YY_h,sum(eddyforce_x.*uu_twa+eddyforce_y.*vv_twa,3));
shading interp;
colorbar;
colormap redblue;
caxis([-3 3]*1e-3);




MM_surf = MM_tavg(:,:,1);
MM_cntr = median(MM_surf(:));

figure(fignum);
fignum = fignum + 1;
pcolor(XX_h,YY_h,MM_surf);
shading interp;
colorbar;
colormap redblue;
hold on;
contour(XX_h,YY_h,MM_surf,[MM_cntr MM_cntr],'EdgeColor','k');
hold off

tau_cntr = sum(sum(curl_tau(MM_surf<MM_cntr)*dx*dy));
hgradM_eddy_cntr = sum(sum(curl_hgradM_eddy(MM_surf<MM_cntr)*dx*dy));
hgradM_mean_cntr = sum(sum(curl_hgradM_mean(MM_surf<MM_cntr)*dx*dy));
eddyforce_cntr = sum(sum(curl_eddyforce(MM_surf<MM_cntr)*dx*dy));
diff_twa_cntr = sum(sum(diff_curl_twa(MM_surf<MM_cntr)*dx*dy));

load(fullfile('products',['kap_nu_',run_name,'.mat']));
IPT_kap = kap_map.*(f0^2/gg(2)).*diff_curl_twa;
IPT_kap(:,1) = 0;
hgradM_kap_cntr = sum(sum(IPT_kap(MM_surf<MM_cntr)*dx*dy));

%%%
%%% Estimates eddy viscosity and diffusivity from time-mean model output.
%%%

%%% Run to load
run_name = 'ACC_AABW_Ny256_Nlay2_tauM0.1_tauP0_tauF0_wDiaM0_wDiaP0_wDiaF0_Cd2.000e-03_rb0.000e+00_diags'
% run_name = 'ACC_AABW_ML_randWdia_randTau_white';

%%% Load parameters   
local_home_dir = '/Volumes/Kilchoman/UCLA/Projects/AWSIM_WindAABW/runs';
% local_home_dir = '/Volumes/Kilchoman/UCLA/Projects/AWSIM/runs';
prod_dir = fullfile(local_home_dir,'AWSIM_WindAABW_products');
loadParams;
dirpath = fullfile(local_home_dir,run_name);
gtild = reshape(cumsum(gg),[Nlay 1 1]);

%%% Averaging period
tmin = 0.5*t1year;
tmax = 30.5*t1year;

%%% Time-average required model output variables
hh_tavg = do_avg(dirpath,OUTN_H_AVG,Nx,Ny,Nlay,n0_avg,N_avg,dt_avg,tmin,tmax,startTime);
uu_tavg = do_avg(dirpath,OUTN_U_AVG,Nx,Ny,Nlay,n0_avg,N_avg,dt_avg,tmin,tmax,startTime);
vv_tavg = do_avg(dirpath,OUTN_V_AVG,Nx,Ny,Nlay,n0_avg,N_avg,dt_avg,tmin,tmax,startTime);
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
                .* 0.25.*(hh_tavg(1:Nx,1:Ny,:)+hh_tavg(1:Nx,[Ny 1:Ny-1],:)+hh_tavg([Nx 1:Nx-1],1:Ny,:)+hh_tavg([Nx 1:Nx-1],[Ny 1:Ny-1],:));                  
husq_eddy = husq_tavg - husq_mean;
hvsq_eddy = hvsq_tavg - hvsq_mean;
huv_eddy = huv_tavg - huv_mean;
%%% TODO think about grid locations
dx_husq_eddy = (husq_eddy(1:Nx,:,:)-husq_eddy([Nx 1:Nx-1],:,:))/dx;
dy_hvsq_eddy = (hvsq_eddy(:,1:Ny,:)-hvsq_eddy(:,[Ny 1:Ny-1],:))/dy;
dx_huv_eddy = (huv_eddy([2:Nx 1],:,:)-huv_eddy(1:Nx,:,:))/dx;
dy_huv_eddy = (huv_eddy(:,[2:Ny 1],:)-huv_eddy(:,1:Ny,:))/dy;

%%% Depth-integrated EKE
EKE = 0.25*(husq_eddy(1:Nx,:,:)+husq_eddy([2:Nx 1],:,:)+hvsq_eddy(:,1:Ny,:)+hvsq_eddy(:,[2:Ny 1],:));
EKE_zint = sum(EKE,3);
EKE_zavg = EKE_zint ./ sum(hh_tavg,3);
  
  
%%% Eddy diffusivity estimate
f0 = mean(mean(2*Omega_z));
IPT_eddy = (hdMdy_eddy(1:Nx,:,1)-hdMdy_eddy([Nx 1:Nx-1],:,1))/dx - (hdMdx_eddy(:,1:Ny,1)-hdMdx_eddy(:,[Ny 1:Ny-1],1))/dy;
IPT_eddy(:,1) = 0;
curl_twa = (vv_twa(1:Nx,:,:)-vv_twa([Nx 1:Nx-1],:,:))/dx - (uu_twa(:,1:Ny,:)-uu_twa(:,[Ny 1:Ny-1],:))/dy;
curl_twa(:,1,:) = 0;
diff_curl_twa = diff(curl_twa,1,3);

%%% Plot eddy forcing
fignum = 1;

figure(fignum);
fignum = fignum + 1;
pcolor(XX_h,YY_h,IPT_eddy);
shading interp;
colorbar;

figure(fignum);
fignum = fignum + 1;
pcolor(XX_h,YY_h,hdMdx_tavg(:,:,1));
shading interp;
colorbar;

figure(fignum);
fignum = fignum + 1;
pcolor(XX_h,YY_h,hdMdx_mean(:,:,1));
shading interp;
colorbar;

figure(fignum);
fignum = fignum + 1;
pcolor(XX_h,YY_h,hdMdx_eddy(:,:,1));
shading interp;
colorbar;

%%% Another attempt, removing small EKE values
denom = diff_curl_twa(:);
numer = -(gg(2)/f0^2) * IPT_eddy(:);
cond = EKE_zavg;
idx = find(cond<median(cond)+1*iqr(cond));
denom(idx) = [];
numer(idx) = [];
cond(idx) = [];
kap_vort_1 = denom(:) \ numer(:)
kap_vort_2 = 1 / (numer(:) \ denom(:))
length(numer)
length(diff_curl_twa(:))

X = [denom,numer];
stdX = std(X);
[coeff,score,latent,tsquared,explained,mu] = pca(X./std(X));
t = [min(score(:,1)), max(score(:,1))];
dirVect = coeff(:,1);
kap_pca = std(numer)/std(denom);%/coeff(1,1)

figure(fignum);
fignum = fignum+1;
% scatter(numer(:),denom(:),10,EKE_zint(:))
scatter(denom(:),numer(:),10,cond(:));
hold on;
plot([min(diff_curl_twa(:)) max(diff_curl_twa(:))],kap_vort_1*[min(diff_curl_twa(:)) max(diff_curl_twa(:))],'k--');
plot([min(diff_curl_twa(:)) max(diff_curl_twa(:))],kap_vort_2*[min(diff_curl_twa(:)) max(diff_curl_twa(:))],'r--');
plot([min(diff_curl_twa(:)) max(diff_curl_twa(:))],kap_pca*[min(diff_curl_twa(:)) max(diff_curl_twa(:))],'k-');
hold off
colormap(cmocean('amp'));
colorbar;
title('Modified IPT/TWA curl relationship');



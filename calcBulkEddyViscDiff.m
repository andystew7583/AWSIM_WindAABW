%%%
%%% calcBulkEddyViscDiff.m
%%%
%%% Estimates bulk eddy viscosity and diffusivity from time-mean model output.
%%%
function [kap_bulk,nu_bulk,r_kap_bulk,r_nu_bulk,EKE_zavg] = calcEddyViscDiff (local_home_dir,run_name)

  %%% Load parameters    
  loadParams;
  
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
%   idx_EKE = find(EKE_zavg>0);
 
  %%% For eddy diffusivity estimate  
  f0 = mean(mean(2*Omega_z));
  IPT_eddy = (hdMdy_eddy(1:Nx,:,1)-hdMdy_eddy([Nx 1:Nx-1],:,1))/dx - (hdMdx_eddy(:,1:Ny,1)-hdMdx_eddy(:,[Ny 1:Ny-1],1))/dy;
  IPT_eddy(:,1) = 0;
  curl_twa = (vv_twa(1:Nx,:,:)-vv_twa([Nx 1:Nx-1],:,:))/dx - (uu_twa(:,1:Ny,:)-uu_twa(:,[Ny 1:Ny-1],:))/dy;
  curl_twa(:,1,:) = 0;
  diff_curl_twa = diff(curl_twa,1,3);
  
  %%% For eddy viscosity estimate
  layidx = 2;
  div_Uz_eddy = (eddyforce_y(1:Nx,:,layidx)-eddyforce_y([Nx 1:Nx-1],:,layidx))/dx - (eddyforce_x(:,1:Ny,layidx)-eddyforce_x(:,[Ny 1:Ny-1],layidx))/dy;
  div_Uz_eddy = div_Uz_eddy ./ hh_q_tavg(:,:,layidx);
  div_Uz_eddy(:,1:2) = 0;
  div_Uz_eddy(:,Ny-1:Ny) = 0;
  zeta_twa = (vv_twa(1:Nx,:,layidx)-vv_twa([Nx 1:Nx-1],:,layidx))/dx - (uu_twa(:,1:Ny,layidx)-uu_twa(:,[Ny 1:Ny-1],layidx))/dy;    
  del2_zeta_twa = zeros(Nx,Ny);
  del2_zeta_twa(:,2:Ny-1) = (zeta_twa([2:Nx 1],2:Ny-1) - 2*zeta_twa(1:Nx,2:Ny-1) + zeta_twa([Nx 1:Nx-1],2:Ny-1)) / dx^2 ...
                    + (zeta_twa(:,3:Ny) - 2*zeta_twa(:,2:Ny-1) + zeta_twa(:,1:Ny-2)) / dy^2;
  del2_zeta_twa(:,1:2) = 0;
  del2_zeta_twa(:,Ny-1:Ny) = 0;

  %%% Compute eddy diffusivity  
  numer = -(gg(2)/f0^2)*IPT_eddy(idx_EKE);
  denom = diff_curl_twa(idx_EKE);
  X = [denom,numer];  
  [coeff,score,roots] = pca(X./std(X));  
  kap_bulk = std(numer)/std(denom)*coeff(2,1)/coeff(1,1);
%   kap_bulk = denom \ numer;
  r_kap_bulk = corr(denom,numer);
  
  %%% Compute eddy viscosity
  numer = -div_Uz_eddy(idx_EKE);
  denom = del2_zeta_twa(idx_EKE);  
  X = [denom,numer];  
  [coeff,score,roots] = pca(X./std(X));  
  nu_bulk = std(numer)/std(denom)*coeff(2,1)/coeff(1,1);
%   nu_bulk = denom \ numer;
  r_nu_bulk = corr(denom,numer);

end


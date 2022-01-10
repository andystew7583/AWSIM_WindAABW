%%%
%%% Estimates eddy viscosity and diffusivity from time-mean model output.
%%%
function [kap0,nu0] = calcEddyViscDiff (local_home_dir,run_name)

  %%% Load parameters    
  loadParams;

  %%% Load time-mean products
  load([run_name,'_tavg.mat']);
  load([run_name,'_QGtavg.mat']);

  %%% Thickness-weighted average velocities
  uu_twa = hu_tavg ./ hh_w_tavg;
  vv_twa = hv_tavg ./ hh_s_tavg;
  usq_twa = husq_tavg ./ hh_w_tavg;
  vsq_twa = hvsq_tavg ./ hh_s_tavg;
  uu_bol = uu_twa - uu_tavg;
  vv_bol = vv_twa - vv_tavg;
  dx_h_mean = (hh_tavg(1:Nx,:,:)-hh_tavg([Nx 1:Nx-1],:,:))/dx;
  dy_h_mean = (hh_tavg(:,1:Ny,:)-hh_tavg(:,[Ny 1:Ny-1],:))/dy;
  dx_hu_mean = (hh_w_tavg([2:Nx 1],:,:).*uu_tavg([2:Nx 1],:,:)-hh_w_tavg(1:Nx,:,:).*uu_tavg(1:Nx,:,:))/dx;
  dy_hv_mean = (hh_s_tavg(:,[2:Ny 1],:).*vv_tavg(:,[2:Ny 1],:)-hh_s_tavg(:,1:Ny,:).*vv_tavg(:,1:Ny,:))/dy;

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
  husq_mean = hh_tavg.*0.5.*(uu_twa(1:Nx,:,:).^2+uu_twa([2:Nx 1],:,:).^2);
  hvsq_mean = hh_tavg.*0.5.*(vv_twa(:,1:Ny,:).^2+vv_twa(:,[2:Ny 1],:).^2);
  huv_mean = 0.5.*(uu_twa(1:Nx,:,:)+uu_twa([Nx 1:Nx-1],:,:)) ...
                  .* 0.5.*(vv_twa(:,1:Ny,:)+vv_twa(:,[Ny 1:Ny-1],:)) ...
                  .* 0.25.*(hh_tavg(1:Nx,1:Ny,:)+hh_tavg(1:Nx,[Ny 1:Ny-1],:)+hh_tavg([Nx 1:Nx-1],1:Ny,:)+hh_tavg([Nx 1:Nx-1],[Ny 1:Ny-1],:));                  
  husq_eddy = husq_tavg - husq_mean;
  hvsq_eddy = hvsq_tavg - hvsq_mean;
  huv_eddy = huv_tavg - huv_mean;
  dx_husq_eddy = (husq_eddy(1:Nx,:,:)-husq_eddy([Nx 1:Nx-1],:,:))/dx;
  dy_hvsq_eddy = (hvsq_eddy(:,1:Ny,:)-hvsq_eddy(:,[Ny 1:Ny-1],:))/dy;
  dx_huv_eddy = (huv_eddy([2:Nx 1],:,:)-huv_eddy(1:Nx,:,:))/dx;
  dy_huv_eddy = (huv_eddy(:,[2:Ny 1],:)-huv_eddy(:,1:Ny,:))/dy;

  %%% Total eddy pressure forcing
  hdMdx_mean = hh_w_tavg.*(MM_tavg(1:Nx,:,:)-MM_tavg([Nx 1:Nx-1],:,:))/dx;
  hdMdy_mean = hh_s_tavg.*(MM_tavg(:,1:Ny,:)-MM_tavg(:,[Ny 1:Ny-1],:))/dy;  
  hdMdx_eddy = hdMdx_tavg - hdMdx_mean;
  hdMdy_eddy = hdMdy_tavg - hdMdy_mean;

  %%% Eddy diffusivity estimate  
  f0 = mean(mean(2*Omega_z));
  IPT_eddy = (pdedy_eddy(1:Nx,:,2:Nlay)-pdedy_eddy([Nx 1:Nx-1],:,2:Nlay))/dx - (pdedx_eddy(:,1:Ny,2:Nlay)-pdedx_eddy(:,[Ny 1:Ny-1],2:Nlay))/dy;
  curl_twa = (vv_twa(1:Nx,:,:)-vv_twa([Nx 1:Nx-1],:,:))/dx - (uu_twa(:,1:Ny,:)-uu_twa(:,[Ny 1:Ny-1],:))/dy;
  diff_curl_twa = diff(curl_twa,1,3);
%   curl_min = 1e-8;
%   diff_curl_twa(abs(diff_curl_twa)<curl_min) = sign(diff_curl_twa(abs(diff_curl_twa)<curl_min))*curl_min;
%   kap_vort_1 = diff_curl_twa(:) \ (-(gg(2)/f0^2) * IPT_eddy(:));
%   kap_vort_2 = 1 / ((-(gg(2)/f0^2) * IPT_eddy(:)) \ diff_curl_twa(:));
%   kap0 = mean([kap_vort_1 kap_vort_2]);
  numer = -(gg(2)/f0^2)*IPT_eddy(:);
  denom = diff_curl_twa(:);
  X = [denom,numer];  
  [coeff,score,roots] = pca(X./std(X));  
  kap0 = std(numer)/std(denom)/coeff(1,1);

  %%% Eddy viscosity estimate
  uzg_eddy = uzg_tavg - ug_tavg.*zetag_tavg;
  vzg_eddy = vzg_tavg - vg_tavg.*zetag_tavg;
  div_Uzg_eddy = (uzg_eddy([2:Nx 1],:,:) - uzg_eddy([Nx 1:Nx-1],:,:)) / (2*dx);
  div_Uzg_eddy(:,2:Ny-1,:) = div_Uzg_eddy(:,2:Ny-1,:) + (vzg_eddy(:,3:Ny,:) - vzg_eddy(:,1:Ny-2,:)) / (2*dy);
  div_Uzg_eddy(:,1,:) = div_Uzg_eddy(:,1,:) + (vzg_eddy(:,2,:) - vzg_eddy(:,1,:)) / dy;
  div_Uzg_eddy(:,Ny,:) = div_Uzg_eddy(:,Ny,:) + (vzg_eddy(:,Ny,:) - vzg_eddy(:,Ny-1,:)) / dy;
  del2_zetag_mean = zeros(Nx,Ny,Nlay);
  del2_zetag_mean(:,2:Ny-1,:) = (zetag_tavg([2:Nx 1],2:Ny-1,:) - 2*zetag_tavg(1:Nx,2:Ny-1,:) + zetag_tavg([Nx 1:Nx-1],2:Ny-1,:)) / dx^2 ...
                    + (zetag_tavg(:,3:Ny,:) - 2*zetag_tavg(:,2:Ny-1,:) + zetag_tavg(:,1:Ny-2,:)) / dy^2;
  del2_zetag_mean(:,1,:) = (zetag_tavg([2:Nx 1],1,:) - 2*zetag_tavg(1:Nx,1,:) + zetag_tavg([Nx 1:Nx-1],1,:)) / dx^2;
  del2_zetag_mean(:,Ny,:) = (zetag_tavg([2:Nx 1],Ny,:) - 2*zetag_tavg(1:Nx,Ny,:) + zetag_tavg([Nx 1:Nx-1],Ny,:)) / dx^2;

  kidx = 1;
%   iidx = find(xx_h>1000*m1km);
  % iidx = find(xx_h>1000*m1km & xx_h<2000*m1km);
  iidx = 1:Nx;
  jidx = find(yy_h>200*m1km & yy_h < Ly-200*m1km);
  numer = -div_Uzg_eddy(iidx,jidx,kidx);
  denom = del2_zetag_mean(iidx,jidx,kidx);
  numer = numer(:);
  denom = denom(:);
%   nu_1 = 1/ (numer(:) \ denom(:));
%   nu_2 = (denom(:) \ numer(:));
%   nu0 = mean([nu_1 nu_2]);
  X = [denom,numer];  
  [coeff,score,roots] = pca(X./std(X));  
  nu0 = std(numer)/std(denom)/coeff(1,1);

end


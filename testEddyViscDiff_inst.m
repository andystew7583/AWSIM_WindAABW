%%%
%%% Estimates eddy viscosity and diffusivity from time-mean model output.
%%%

%%% Run to load
run_name = 'ACC_AABW_Ny128_Nlay2_tauM0.1_tauP0_tauF0_wDiaM0_wDiaP0_wDiaF0_Cd2.000e-03_rb0.000e+00_diags'
% run_name = 'ACC_AABW_ML_randWdia_randTau_white';

%%% Load parameters   
local_home_dir = '/Volumes/Kilchoman/UCLA/Projects/AWSIM_WindAABW/runs';
% local_home_dir = '/Volumes/Kilchoman/UCLA/Projects/AWSIM/runs';
prod_dir = fullfile(local_home_dir,'AWSIM_WindAABW_products');
loadParams;
dirpath = fullfile(local_home_dir,run_name);
gtild = reshape(cumsum(gg),[Nlay 1 1]);

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

%%% Depth-integrated EKE
EKE_zint = sum(husq_eddy+hvsq_eddy,3);
EKE_zavg = EKE_zint ./ sum(hh_tavg,3);

%%% Eddy diffusivity estimate
f0 = mean(mean(2*Omega_z));
% kappa = (1e-2/1e-8*pdedx_eddy(:,:,2:Nlay) .* diff(uu_twa,1,3) + 1e-2/1e-8*pdedy_eddy(:,:,2:Nlay) .* diff(vv_twa,1,3)) ./ (diff(uu_twa,1,3).^2 + diff(vv_twa,1,3).^2);
grad_h_squared = (dx_h_mean(:,:,1).^2 + dy_h_mean(:,:,1).^2);
grad_h_thresh = sqrt(1e-7);
grad_h_squared(grad_h_squared<grad_h_thresh^2) = sign(grad_h_squared(grad_h_squared<grad_h_thresh^2))*grad_h_thresh^2;
kappa = - ( hh_w_tavg(:,:,1).*uu_bol(:,:,1).*dx_h_mean(:,:,1) + hh_s_tavg(:,:,1).*vv_bol(:,:,1).*dy_h_mean(:,:,1) ) ./ grad_h_squared;

% kap0 = - sum(sum( hh_w_tavg(:,:,1).*uu_bol(:,:,1).*dx_h_mean(:,:,1) + hh_s_tavg(:,:,1).*vv_bol(:,:,1).*dy_h_mean(:,:,1) )) ./ sum(sum(dx_h_mean(:,:,1).^2 + dy_h_mean(:,:,1).^2))
% kap0 = sum(sum( hh_tavg(:,:,1).*(dx_hu_mean(:,:,1)+dy_hv_mean(:,:,1)) )) ./ sum(sum(dx_h_mean(:,:,1).^2+dy_h_mean(:,:,1).^2));
PEtoEKE1 =-  pdedx_eddy(:,:,2:Nlay) .* diff(uu_twa,1,3) - pdedy_eddy(:,:,2:Nlay) .* diff(vv_twa,1,3) + sum(dx_hphi_eddy.*uu_twa+dy_hphi_eddy.*vv_twa,3);
% kap0 = (gg(2)/f0^2)*sum(sum(pdedx_eddy(:,:,2:Nlay) .* diff(uu_twa,1,3) + pdedy_eddy(:,:,2:Nlay) .* diff(vv_twa,1,3))) / sum(sum(diff(uu_twa,1,3).^2 + diff(vv_twa,1,3).^2))
PEtoEKE = hdMdx_eddy.*uu_twa + hdMdy_eddy.*vv_twa;

IPT_eddy = (pdedy_eddy(1:Nx,:,2:Nlay)-pdedy_eddy([Nx 1:Nx-1],:,2:Nlay))/dx - (pdedx_eddy(:,1:Ny,2:Nlay)-pdedx_eddy(:,[Ny 1:Ny-1],2:Nlay))/dy;
curl_twa = (vv_twa(1:Nx,:,:)-vv_twa([Nx 1:Nx-1],:,:))/dx - (uu_twa(:,1:Ny,:)-uu_twa(:,[Ny 1:Ny-1],:))/dy;
diff_curl_twa = diff(curl_twa,1,3);
curl_min = 1e-8;
diff_curl_twa(abs(diff_curl_twa)<curl_min) = sign(diff_curl_twa(abs(diff_curl_twa)<curl_min))*curl_min;
kappa_vort = -(gg(2)/f0^2) * IPT_eddy ./ diff_curl_twa;

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
dx_zetag_mean = (zetag_tavg([2:Nx 1],:,:) - zetag_tavg([Nx 1:Nx-1],:,:)) / (2*dx);
dy_zetag_mean = zeros(Nx,Ny,Nlay);
dy_zetag_mean(:,2:Ny-1,:) = (zetag_tavg(:,3:Ny,:) - zetag_tavg(:,1:Ny-2,:)) / (2*dy);
dy_zetag_mean(:,1,:) =  (zetag_tavg(:,2,:) - zetag_tavg(:,1,:)) / dy;
dy_zetag_mean(:,Ny,:) = (zetag_tavg(:,Ny,:) - zetag_tavg(:,Ny-1,:)) / dy;


%%% Plot eddy forcing
fignum = 1;

figure(fignum);
fignum = fignum+1;
pcolor(XX_h,YY_h,kappa);
shading interp
colorbar;
colormap redblue;
title('Eddy diffusivity');

figure(fignum);
fignum = fignum+1;
pcolor(XX_h,YY_h,sum(PEtoEKE1,3));
shading interp
colorbar;
colormap redblue;
title('PE to EKE');

figure(fignum);
fignum = fignum+1;
pcolor(XX_h,YY_h,sum(PEtoEKE,3));
shading interp
colorbar;
colormap redblue;
title('PE to EKE');


figure(fignum);
fignum = fignum+1;
pcolor(XX_h,YY_h,PEtoEKE(:,:,1));
shading interp
colorbar;
colormap redblue;
title('PE to EKE, layer 1');


figure(fignum);
fignum = fignum+1;
pcolor(XX_h,YY_h,PEtoEKE(:,:,2));
shading interp
colorbar;
colormap redblue;
title('PE to EKE, layer 2');

figure(fignum);
fignum = fignum+1;
plot(yy_h,sum(hdMdx_eddy(:,:,1)*dx/Lx*1000,1));
hold on;
plot(yy_h,sum(pdedx_eddy(:,:,2)*dx/Lx*1000,1));
hold off;
title('Eddy interfactial form stress');

figure(fignum);
fignum = fignum+1;
plot(yy_h,sum(hdMdx_mean(:,:,1)*dx/Lx*1000,1));
hold on;
plot(yy_h,sum(pdedx_mean(:,:,2)*dx/Lx*1000,1));
hold off;
title('Mean interfactial form stress');

figure(fignum);
fignum = fignum+1;
pcolor(XX_h,YY_h,IPT_eddy);
shading interp
colorbar;
colormap redblue;
title('Eddy IPT');

figure(fignum);
fignum = fignum+1;
pcolor(XX_h,YY_h,curl_twa(:,:,1));
shading interp
colorbar;
colormap redblue;
title('Curl of TWA velocity');

figure(fignum);
fignum = fignum+1;
pcolor(XX_h,YY_h,kappa_vort);
shading interp
colorbar;
colormap redblue;
caxis([-2000 2000]);
title('Kappa from curl');


figure(fignum);
fignum = fignum+1;
pcolor(XX_h,YY_h,hh_tavg(:,:,1));
shading interp
colorbar;
colormap redblue;
title('Mean upper layer thickness');



%%% Test estimate of kappa via linear regression
figure(fignum);
fignum = fignum+1;
scatter(diff_curl_twa(:),-(gg(2)/f0^2) * IPT_eddy(:));
kap_vort_1 = diff_curl_twa(:) \ (-(gg(2)/f0^2) * IPT_eddy(:))
kap_vort_2 = 1 / ((-(gg(2)/f0^2) * IPT_eddy(:)) \ diff_curl_twa(:))
hold on;
plot([min(diff_curl_twa(:)) max(diff_curl_twa(:))],kap_vort_1*[min(diff_curl_twa(:)) max(diff_curl_twa(:))],'k--');
plot([min(diff_curl_twa(:)) max(diff_curl_twa(:))],kap_vort_2*[min(diff_curl_twa(:)) max(diff_curl_twa(:))],'r--');
plot([min(diff_curl_twa(:)) max(diff_curl_twa(:))],mean([kap_vort_1 kap_vort_2])*[min(diff_curl_twa(:)) max(diff_curl_twa(:))],'r--');
hold off;
shading interp
colorbar;
colormap redblue;
title('IPT/TWA curl relationship');




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



% denom = diff_curl_twa(:);
% numer = -(gg(2)/f0^2) * IPT_eddy(:);
% corr_radius = 200*m1km;
% 
% kap_map = zeros(Nx,Ny);
% for i=1:Nx
%   i
%   for j=1:Ny
%    
%     %%% Create modified x-matrix to account for periodicity
%     Xdiff = XX_h-XX_h(i,j);
%     idx = abs(Xdiff)>Lx/2;
%     XX_mod = XX_h;
%     XX_mod(idx) = XX_mod(idx)-Lx*sign(Xdiff(idx));
%     
%     %%% Find points within specified radius
%     R = sqrt((XX_mod-XX_h(i,j)).^2 + (YY_h-YY_h(i,j)).^2);
%     idx = find(R<corr_radius);   
% 
%     %%% Compute kappa;
%     X = [denom(idx),numer(idx)];
%     std(X)
%     [coeff,score,latent,tsquared,explained,mu] = pca(X./std(X));
%     kap_map(i,j) = std(numer(idx))/std(denom(idx))*coeff(2,1)/coeff(1,1);
%     coeff(1,1)
% %     kap_map(i,j) = denom(idx) \ numer(idx);    
%     kap_map(i,j)
%     
%     figure(fignum);
%     fignum = fignum+1;
%     scatter(denom(idx),numer(idx));
%     hold on;
%     plot([min(denom(idx)) max(denom(idx))],kap_map(i,j)*[min(denom(idx)) max(denom(idx))],'k-');     
%     hold off;
%     
%   end  
% end
% 
% figure(fignum);
% fignum = fignum+1;
% pcolor(XX_h,YY_h,kap_map);
% shading interp
% colormap(cmocean('amp'));






EKE_thresh = [1e-4:1e-4:0.1]; 
corr_loEKE = 0*EKE_thresh;
corr_hiEKE = 0*EKE_thresh;
pval_loEKE = 0*EKE_thresh;
pval_hiEKE = 0*EKE_thresh;
kap1_hiEKE = 0*EKE_thresh;
kap2_hiEKE = 0*EKE_thresh;
for n=1:length(EKE_thresh)
  denom = diff_curl_twa(:);
  numer = -(gg(2)/f0^2) * IPT_eddy(:);
  idx_lo = find(EKE_zavg<EKE_thresh(n));
  idx_hi = find(EKE_zavg>=EKE_thresh(n));
  if (isempty(idx_lo) || isempty(idx_hi))
    continue;
  end
  [r,p] = corr(numer(idx_lo),denom(idx_lo));
  corr_loEKE(n) = r;
  pval_loEKE(n) = p;
  [r,p] = corr(numer(idx_hi),denom(idx_hi));
  corr_hiEKE(n) = r;
  pval_hiEKE(n) = p;
  kap1_hiEKE(n) = 1 / (numer(idx_hi) \ denom(idx_hi));
  kap2_hiEKE(n) = denom(idx_hi) \ numer(idx_hi);
end



figure(fignum);
fignum = fignum+1;
plot(EKE_thresh,corr_loEKE)
hold on
plot(EKE_thresh,corr_hiEKE)
hold off;

figure(fignum);
fignum = fignum+1;
plot(EKE_thresh,kap1_hiEKE)
hold on
plot(EKE_thresh,kap2_hiEKE)
hold off;

figure(fignum);
fignum = fignum+1;
pcolor(XX_h,YY_h,EKE_zavg);
shading interp
colormap(cmocean('amp'));
colorbar
hold on
contour(XX_h,YY_h,EKE_zavg,[median(EKE_zavg(:)) median(EKE_zavg(:))],'EdgeColor','k')
hold off
colormap(cmocean('amp',20));
% caxis([0 0.1])


logEKE_min = -3.5;
logEKE_max = -1;
dlogEKE = 0.05;
logEKE_band = logEKE_min:dlogEKE:logEKE_max;
kap1_band = NaN*logEKE_band;
kap2_band = NaN*logEKE_band;
for n=1:length(logEKE_band)
  denom = diff_curl_twa(:);
  numer = -(gg(2)/f0^2) * IPT_eddy(:);
  idx_band = find((log10(EKE_zavg)>=logEKE_band(n)-dlogEKE/2) & (log10(EKE_zavg)<logEKE_band(n)+dlogEKE/2));
  if (isempty(idx_band))
    continue;
  end
  [r,p] = corr(numer(idx_band),denom(idx_band));
  kap1_band(n) = 1 / (numer(idx_band) \ denom(idx_band));
  kap2_band(n) = denom(idx_band) \ numer(idx_band);  
end

figure(fignum);
fignum = fignum+1;
plot(10.^logEKE_band,kap1_band);
hold on;
plot(10.^logEKE_band,kap2_band);
hold off;


iidx = find(xx_h>1000*m1km);
% iidx = find(xx_h>1000*m1km & xx_h<2000*m1km);
% iidx = 1:Nx;
jidx = find(yy_h>400*m1km & yy_h < Ly-400*m1km);
denom = diff_curl_twa(iidx,jidx,:);
numer = -(gg(2)/f0^2) * IPT_eddy(iidx,jidx,:);
% denom(abs(numer)<1e-12) = sign(denom(abs(numer)<1e-12))*1e-12;
kap_geog_1 = 1/ (numer(:) \ denom(:))
kap_geog_2 = (denom(:) \ numer(:))
figure(fignum);
fignum = fignum+1;
scatter(denom(:),numer(:));
hold on;
plot([min(denom(:)) max(denom(:))],nu_1*[min(denom(:)) max(denom(:))],'k--');
plot([min(denom(:)) max(denom(:))],nu_2*[min(denom(:)) max(denom(:))],'k--');
hold off;
corr(denom(:),numer(:))
mean([kap_geog_1 kap_geog_2])


tau_vals = [0.01 0.017 0.03 0.05 0.1 0.17 0.3];
kap_vals = [195 326 367 410 468 509 627];
figure(fignum);
fignum = fignum+1;
plot(tau_vals.^.5,kap_vals);
title('Transient eddy diffusivity vs. wind');




figure(fignum);
fignum = fignum+1;
pcolor(XX_h,YY_h,del2_zetag_mean(:,:,1));
shading interp
colorbar;
colormap redblue;
caxis([-2 2]*1e-14)
title('Laplacian of mean geostrophic vorticity');

figure(fignum);
fignum = fignum+1;
pcolor(XX_h,YY_h,div_Uzg_eddy(:,:,1));
shading interp
colorbar;
colormap redblue;
caxis([-5 5]*1e-12)
title('Divergence of geostrophic eddy vorticity flux');

kidx = 1;
iidx = find(xx_h>1000*m1km);
% iidx = find(xx_h>1000*m1km & xx_h<2000*m1km);
% iidx = 1:Nx;
jidx = find(yy_h>400*m1km & yy_h < Ly-400*m1km);
numer = -div_Uzg_eddy(iidx,jidx,kidx);
denom = del2_zetag_mean(iidx,jidx,kidx);
% denom(abs(numer)<1e-12) = sign(denom(abs(numer)<1e-12))*1e-12;
nu_1 = 1/ (numer(:) \ denom(:))
nu_2 = (denom(:) \ numer(:))
figure(fignum);
fignum = fignum+1;
scatter(denom(:),numer(:));
hold on;
plot([min(denom(:)) max(denom(:))],nu_1*[min(denom(:)) max(denom(:))],'k--');
plot([min(denom(:)) max(denom(:))],nu_2*[min(denom(:)) max(denom(:))],'k--');
hold off;
corr(denom(:),numer(:))
mean([nu_1 nu_2])

numer = -div_Uzg_eddy(:,:,kidx);
denom = del2_zetag_mean(:,:,kidx);
figure(fignum);
fignum = fignum+1;
pcolor(XX_h,YY_h,numer./denom);
shading interp
colorbar;
colormap redblue;
caxis([-3000 3000])
title('nu estimated from divergence');

% iidx = find(xx_h>1000*m1km);
% jidx = find(yy_h>200*m1km & yy_h < Ly-200*m1km);
% numer = - (uzg_eddy(iidx,jidx,1).*dx_zetag_mean(iidx,jidx,1) + vzg_eddy(iidx,jidx,1).*dy_zetag_mean(iidx,jidx,1));
% denom = dx_zetag_mean(iidx,jidx,1).^2 + dy_zetag_mean(iidx,jidx,1).^2;
% nu_1 = 1/ (numer(:) \ denom(:))
% nu_2 = (denom(:) \ numer(:))
% figure(fignum);
% fignum = fignum+1;
% scatter(denom(:),numer(:));
% hold on;
% plot([min(denom(:)) max(denom(:))],nu_1*[min(denom(:)) max(denom(:))],'k--');
% plot([min(denom(:)) max(denom(:))],nu_2*[min(denom(:)) max(denom(:))],'k--');
% hold off;


%%% Another attempt, removing small EKE values
kidx = 1;
numer = -div_Uzg_eddy(:,:,kidx);
denom = del2_zetag_mean(:,:,kidx);
numer = numer(:);
denom = denom(:);
cond = EKE_zavg(:);
idx = find(cond<median(cond)+0*iqr(cond));
denom(idx) = [];
numer(idx) = [];
cond(idx) = [];
nu_vort_1 = denom(:) \ numer(:)
nu_vort_2 = 1 / (numer(:) \ denom(:))
X = [denom,numer];  
[coeff,score,roots] = pca(X./std(X));  
nu_pca = std(numer)/std(denom)/coeff(1,1)

figure(fignum);
fignum = fignum+1;
% scatter(numer(:),denom(:),10,EKE_zint(:))
scatter(denom(:),numer(:),10,cond(:));
hold on;
plot([min(denom(:)) max(denom(:))],nu_vort_1*[min(denom(:)) max(denom(:))],'k--');
plot([min(denom(:)) max(denom(:))],nu_vort_2*[min(denom(:)) max(denom(:))],'r--');
plot([min(denom(:)) max(denom(:))],mean([nu_vort_1 nu_vort_2])*[min(denom(:)) max(denom(:))],'r--');
plot([min(denom(:)) max(denom(:))],nu_pca*[min(denom(:)) max(denom(:))],'k-');
hold off
colormap(cmocean('amp'));
colorbar;
title('nu estimated from divergence for high EKE');


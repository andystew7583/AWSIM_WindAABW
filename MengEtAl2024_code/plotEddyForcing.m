%%%
%%% plotEddyForcing.m
%%%
%%% Calculates and plots eddy forcing terms from time-mean model output.
%%%

%%% Run to load
run_name = 'ACC_AABW_ML_randWdia_randTau_white';

%%% Load parameters   
local_home_dir = '/Volumes/Kilchoman/UCLA/Projects/AWSIM/runs';
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

%%% Plot eddy forcing
fignum = 1;

figure(fignum);
fignum = fignum+1;
pcolor(XX_h,YY_h,-dx_hphi_eddy(:,:,1));
shading interp
colorbar;
colormap redblue;
title('Zonal eddy pressure flux convergence');

figure(fignum);
fignum = fignum+1;
pcolor(XX_h,YY_h,-(dx_husq_eddy(:,:,1)+dy_huv_eddy(:,:,1)));
shading interp
colorbar;
colormap redblue;
title('Zonal momentum flux convergence');

figure(fignum);
fignum = fignum+1;
pcolor(XX_h,YY_h,pdedx_eddy(:,:,1)-pdedx_eddy(:,:,2));
shading interp
colorbar;
colormap redblue;
title('Zonal isopycnal form stress');

figure(fignum);
fignum = fignum+1;
pcolor(XX_h,YY_h,-hdMdx_eddy(:,:,1));
shading interp
colorbar;
colormap redblue;
title('Zonal eddy pressure forcing');

figure(fignum);
fignum = fignum+1;
pcolor(XX_h,YY_h,-dy_hphi_eddy(:,:,1));
shading interp
colorbar;
colormap redblue;
title('Meridional eddy pressure flux convergence');

figure(fignum);
fignum = fignum+1;
pcolor(XX_h,YY_h,-(dx_huv_eddy(:,:,1)+dy_hvsq_eddy(:,:,1)));
shading interp
colorbar;
colormap redblue;
title('Meridional momentum flux convergence');

figure(fignum);
fignum = fignum+1;
pcolor(XX_h,YY_h,pdedy_eddy(:,:,1)-pdedy_eddy(:,:,2));
shading interp
colorbar;
colormap redblue;
title('Meridional isopycnal form stress');

figure(fignum);
fignum = fignum+1;
pcolor(XX_h,YY_h,-hdMdy_eddy(:,:,1));
shading interp
colorbar;
colormap redblue;
title('Meridional eddy pressure forcing');


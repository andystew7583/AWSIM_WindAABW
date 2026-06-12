%%%
%%% calcVortBalance.m
%%%
%%% Computes time series of wind stress and form stresses.
%%%

%%% Load constant parameters
constants;
Ny = 128;
Nlay = 3;
tau_mean = [0.15];
tau_pert = 0;
tau_freq = 0;
AABW_mean = 0;
AABW_pert = 0;
AABW_freq = 0;
quad_drag = 2e-3;
lin_drag = 0e-4;  
topog_width = 150;
topog_height = 1000;
rough_topog = true;
n_E = 1;
run_name = constructRunName (false,Ny,Nlay, ...
                                  tau_mean,tau_pert,tau_freq, ...
                                  AABW_mean,AABW_pert,AABW_freq, ...
                                  quad_drag,lin_drag,topog_width,topog_height,rough_topog,n_E);

%%% Load parameters   
local_home_dir = '/Volumes/Stewart-RAID1-A/UCLA/Projects/AWSIM_WindAABW/runs_varywind';
prod_dir = fullfile(local_home_dir,'products');
loadParams;
dirpath = fullfile(local_home_dir,run_name);
gtild = reshape(cumsum(gg),[Nlay 1 1]);
rho0 = 1000;

%%% Max time at which to load transports
tmin = 0;
tmax = 15*t1year;

%%% Load time-averaged momentum budget terms
Umom_hgradM = do_avg (dirpath,OUTN_UMOM_GRADM,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,0);
Vmom_hgradM = do_avg (dirpath,OUTN_VMOM_GRADM,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,0);
Umom_tau = do_avg (dirpath,OUTN_UMOM_WIND,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,0);
Vmom_tau = do_avg (dirpath,OUTN_VMOM_WIND,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,0);
Umom_CdBot = do_avg (dirpath,OUTN_UMOM_CDBOT,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,0); 
Vmom_CdBot = do_avg (dirpath,OUTN_VMOM_CDBOT,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,0); 

Umom_dhdt = do_avg (dirpath,OUTN_UMOM_DHDT,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,0); 
Umom_q = do_avg (dirpath,OUTN_UMOM_Q,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,0); 
Umom_gradKE = do_avg (dirpath,OUTN_UMOM_GRADKE,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,0); 
Umom_adv = Umom_q + Umom_gradKE + Umom_dhdt;

Vmom_dhdt = do_avg (dirpath,OUTN_VMOM_DHDT,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,0); 
Vmom_q = do_avg (dirpath,OUTN_VMOM_Q,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,0); 
Vmom_gradKE = do_avg (dirpath,OUTN_VMOM_GRADKE,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,0); 
Vmom_adv = Vmom_q + Vmom_gradKE + Vmom_dhdt;

%%% Load time-averaged products
huv_tavg = do_avg (dirpath,OUTN_VMOM_DHDT,Nx,Ny,Nlay,n0_avg,N_avg,dt_avg,tmin,tmax,0);
hv_tavg = do_avg (dirpath,OUTN_HV_AVG,Nx,Ny,Nlay,n0_avg,N_avg,dt_avg,tmin,tmax,0);
hu_tavg = do_avg (dirpath,OUTN_HU_AVG,Nx,Ny,Nlay,n0_avg,N_avg,dt_avg,tmin,tmax,0);

%%% Coriolis term
ff = 2*Omega_z;
beta = 2*(Omega_z(1,end)-Omega_z(1,1)) / (yy_q(1,end)-yy_q(1,1));
Umom_cori = repmat(ff,[1 1 Nlay]);
Umom_cori = 0.5*(Umom_cori(1:Nx,:,:)+Umom_cori(2:Nx+1,:,:));
Umom_cori(:,1:Ny,:) = Umom_cori(:,1:Ny,:).*hv_tavg;
Umom_cori(:,Ny+1,:) = 0;
Umom_cori = 0.5*(Umom_cori(:,1:Ny,:)+Umom_cori(:,2:Ny+1,:));
Vmom_cori = repmat(ff(1:Nx,:),[1 1 Nlay]);
Vmom_cori = 0.5*(Vmom_cori(:,1:Ny,:)+Vmom_cori(:,2:Ny+1,:));
Vmom_cori = - Vmom_cori.*hu_tavg;
Vmom_cori = 0.5*(Vmom_cori(1:Nx,:,:)+Vmom_cori([Nx 1:Nx-1],:,:));


curl_hgradM = ( (Vmom_hgradM(1:Nx,1:Ny,:) - Vmom_hgradM([Nx 1:Nx-1],1:Ny,:)) / dx ) ...
            - ( (Umom_hgradM(1:Nx,1:Ny,:) - Umom_hgradM(1:Nx,[Ny 1:Ny-1],:)) / dy );
curl_hgradM(:,1,:) = 0;
          
curl_tau = ( (Vmom_tau(1:Nx,1:Ny,:) - Vmom_tau([Nx 1:Nx-1],1:Ny,:)) / dx ) ...
         - ( (Umom_tau(1:Nx,1:Ny,:) - Umom_tau(1:Nx,[Ny 1:Ny-1],:)) / dy );
curl_tau(:,1,:) = 0;
       
curl_adv = ( (Vmom_adv(1:Nx,1:Ny,:) - Vmom_adv([Nx 1:Nx-1],1:Ny,:)) / dx ) ...
         - ( (Umom_adv(1:Nx,1:Ny,:) - Umom_adv(1:Nx,[Ny 1:Ny-1],:)) / dy );
curl_adv(:,1,:) = 0;

curl_CdBot = ( (Vmom_CdBot(1:Nx,1:Ny,:) - Vmom_CdBot([Nx 1:Nx-1],1:Ny,:)) / dx ) ...
         - ( (Umom_CdBot(1:Nx,1:Ny,:) - Umom_CdBot(1:Nx,[Ny 1:Ny-1],:)) / dy );
curl_CdBot(:,1,:) = 0;

curl_cori = ( (Vmom_cori(1:Nx,1:Ny,:) - Vmom_cori([Nx 1:Nx-1],1:Ny,:)) / dx ) ...
          - ( (Umom_cori(1:Nx,1:Ny,:) - Umom_cori(1:Nx,[Ny 1:Ny-1],:)) / dy );      
curl_cori(:,1,:) = 0;
curl_cori = beta.*0.5.*(hv_tavg(:,:,:)+hv_tavg([Nx 1:Nx-1],:,:));

fignum = 0;
       
fignum = fignum + 1;
figure(fignum);
pcolor(XX_q(1:Nx,1:Ny)/1000,YY_q(1:Nx,1:Ny)/1000,sum(curl_tau,3));
shading interp;
colorbar
colormap redblue;
caxis([-1 1]*1e-8);

fignum = fignum + 1;
figure(fignum);
pcolor(XX_q(1:Nx,1:Ny)/1000,YY_q(1:Nx,1:Ny)/1000,sum(curl_hgradM,3));
shading interp;
colorbar
colormap redblue;
caxis([-1 1]*1e-8);


fignum = fignum + 1;
figure(fignum);
pcolor(XX_q(1:Nx,1:Ny)/1000,YY_q(1:Nx,1:Ny)/1000,sum(curl_adv,3));
shading interp;
colorbar
colormap redblue;
caxis([-1 1]*1e-8);


fignum = fignum + 1;
figure(fignum);
pcolor(XX_q(1:Nx,1:Ny)/1000,YY_q(1:Nx,1:Ny)/1000,sum(curl_CdBot,3));
shading interp;
colorbar
colormap redblue;
caxis([-1 1]*1e-8);

fignum = fignum + 1;
figure(fignum);
pcolor(XX_q(1:Nx,1:Ny)/1000,YY_q(1:Nx,1:Ny)/1000,sum(curl_cori,3));
shading interp;
colorbar
colormap redblue;
caxis([-1 1]*1e-8);

fignum = fignum + 1;
figure(fignum);
pcolor(XX_q(1:Nx,1:Ny)/1000,YY_q(1:Nx,1:Ny)/1000,sum(curl_adv+curl_hgradM,3));
shading interp;
colorbar
colormap redblue;
caxis([-1 1]*1e-8);



fignum = fignum + 1;
figure(fignum);
pcolor(XX_q(1:Nx,1:Ny)/1000,YY_q(1:Nx,1:Ny)/1000,curl_hgradM(:,:,1));
shading interp;
colorbar
colormap redblue;
caxis([-1 1]*1e-8);


fignum = fignum + 1;
figure(fignum);
pcolor(XX_q(1:Nx,1:Ny)/1000,YY_q(1:Nx,1:Ny)/1000,sum(curl_hgradM(:,:,1:2),3));
shading interp;
colorbar
colormap redblue;
caxis([-1 1]*1e-8);

fignum = fignum + 1;
figure(fignum);
pcolor(XX_q(1:Nx,1:Ny)/1000,YY_q(1:Nx,1:Ny)/1000,curl_cori(:,:,1));
shading interp;
colorbar
colormap redblue;
caxis([-1 1]*1e-8);

fignum = fignum + 1;
figure(fignum);
plot(yy_q(1:Ny),cumsum(sum(sum(curl_hgradM,3)*dx,1)*dy,2));
hold on;
plot(yy_q(1:Ny),cumsum(sum(sum(curl_tau,3)*dx,1)*dy,2));
plot(yy_q(1:Ny),sum(sum(Umom_hgradM,3)*dx,1));
hold off
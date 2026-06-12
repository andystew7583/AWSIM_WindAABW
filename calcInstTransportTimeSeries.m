%%%
%%% calcInstTransportTimeSeries.m
%%%
%%% Computes MOC and SSH/OBP time series.
%%%

%%% Run to load
run_name = 'ACC_AABW_Ny128_Nlay3_tauM0.15_tauP0_tauF0_wDiaM0_wDiaP0_wDiaF0_Cd2.000e-03_rb0.000e+00_2xwind_E10_diags'; %%% Run to analyze

%%% Load parameters   
local_home_dir = '/Volumes/Stewart-RAID1-A/UCLA/Projects/AWSIM_WindAABW/runs_varywind';
prod_dir = fullfile(local_home_dir,'products');
loadParams;
dirpath = fullfile(local_home_dir,run_name);
gtild = reshape(cumsum(gg),[Nlay 1 1]);
rho0 = 1000;

%%% Define range of latitudes over which to compute AABW transport
ymin = 100*m1km;
ymax = Ly-100*m1km;






%%% Time series of wind forcing
[tauNrecs tauNrecs_found] = readparam(params_file,'tauNrecs','%lf'); 
[tauPeriod tauPeriod_found] = readparam(params_file,'tauPeriod','%lf'); 
tauFile = 'taux.dat';
tauMax = zeros(1,tauNrecs);
tauTimes = (0:1:tauNrecs-1)/tauNrecs*tauPeriod;
fid = fopen(fullfile(dirpath,tauFile),'r','b');
if (fid == -1)
  error(['Could not open ',tauFile]);
end
taux = zeros(tauNrecs,Nx,Ny);
for j=1:Ny
  for i=1:Nx
    taux(:,i,j) = fread(fid,[tauNrecs 1],'real*8','ieee-le');
  end
end
for n=1:tauNrecs
  tauMax(n) = max(max(squeeze(taux(n,:,:))*rho0));
end
fclose(fid);

%%% Time series of "true" AABW export
[wDiaNrecs wDiaNrecs_found] = readparam(params_file,'wDiaNrecs','%lf'); 
[wDiaPeriod wDiaPeriod_found] = readparam(params_file,'wDiaPeriod','%lf'); 
wDiaFile = 'wDiaFile.dat';
wDiaInt = zeros(1,wDiaNrecs);
wDiaTimes = (0:1:wDiaNrecs-1)/wDiaNrecs*wDiaPeriod;
fid = fopen(fullfile(dirpath,wDiaFile),'r','b');
if (fid == -1)
  error(['Could not open ',wDiaFile]);
end
wDia = zeros(wDiaNrecs,Nlay+1,Nx,Ny);
for j=1:Ny
  for i=1:Nx
    wDia(:,:,i,j) = fread(fid,[wDiaNrecs Nlay+1],'real*8','ieee-le');
  end
end
for n=1:wDiaNrecs
  wDiaInt(n) = sum(sum(squeeze(wDia(n,Nlay,:,:))*dx*dy));
end
fclose(fid);



%%% Precompute iteration numbers and output times
iters = n0:1:n0+Nframes-1;
tt = startTime + (iters-n0)*dt_s;

%%% At each time iteration...
cntr = 0;
hu_int = NaN*ones(Nlay,Nx,Nframes);
hv_int = NaN*ones(Nlay,Ny,Nframes);
hv_mean_int = NaN*ones(Nlay,Ny,Nframes);
hv_stand_int = NaN*ones(Nlay,Ny,Nframes);
Tacc = zeros(1,Nframes);
Taabw = zeros(1,Nframes);
Taaiw = zeros(1,Nframes);
h = zeros(Nlay,Nx,Ny);
v = zeros(Nlay,Nx,Ny);
u = zeros(Nlay,Nx,Ny);
for n=iters

  cntr = cntr + 1; 
      
  %%% Load instantaneous model state
  for k=1:Nlay
    data_file = fullfile(dirpath,[OUTN_H,num2str(k-1),'_n=',num2str(n),'.dat']);
    h(k,:,:) = readOutputFile(data_file,Nx,Ny);
  end  
  for k=1:Nlay
    data_file = fullfile(dirpath,[OUTN_U,num2str(k-1),'_n=',num2str(n),'.dat']);
    u(k,:,:) = readOutputFile(data_file,Nx,Ny);
  end
  for k=1:Nlay
    data_file = fullfile(dirpath,[OUTN_V,num2str(k-1),'_n=',num2str(n),'.dat']);
    v(k,:,:) = readOutputFile(data_file,Nx,Ny);
  end
  h_v = 0.5*(h(:,:,1:Ny)+h(:,:,[2:Ny 1])); %%% Interpolate to v-points
  h_u = 0.5*(h(:,1:Nx,:)+h(:,[2:Nx 1],:)); %%% Interpolate to u-points 

  %%% Meridional transport decomposition
  hv_int(:,:,cntr) = squeeze(sum(h_v.*v.*dx,2));
  hv_mean_int(:,:,cntr) = squeeze(mean(h_v,2) .* mean(v,2) * Lx);
  hv_stand_int(:,:,cntr) = hv_int(:,:,cntr) - hv_mean_int(:,:,cntr);  

  %%% Zonal transport
  hu_int(:,:,cntr) = squeeze(sum(h_u.*u.*dy,3));  

  %%% AABW export
  yminidx = find(yy_v>ymin,1,'first');
  ymaxidx = find(yy_v>ymax,1,'first');
  Taabw(cntr) = squeeze(mean(hv_int(Nlay,yminidx:ymaxidx,cntr),2));

  %%% AAIW export
  if (Nlay > 2)
    Taaiw(cntr) = squeeze(mean(hv_int(1,yminidx:ymaxidx,cntr),2));
  end

  %%% ACC transport
  Tacc(cntr) = sum(hu_int(:,1,cntr),1);

end

%%% Write chunk to .mat file
save(fullfile(prod_dir,[run_name,'_InstTrans.mat']), ...
  'wDiaTimes','wDiaInt', ...
  'tauTimes','tauMax', ...
  'XX_h','YY_h','hhb','gg', ...
  'tt', ...
  'hu_int','hv_int','hv_mean_int','hv_stand_int', ...
  'Tacc','Taabw','Taaiw', ...
  '-v7.3');

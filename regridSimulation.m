%%%
%%% regridSimulation.m
%%%
%%% Takes output and simulations setup and interpolates to a different grid
%%% size, then creates new output files that can be used to initialize a
%%% new simulation. Ghost points are used to impose walls or periodic
%%% boundary conditions.

local_home_dir = '/Volumes/Stewart-RAID1-A/UCLA/Projects/AWSIM/runs';
run_name = 'ACC_AABW_ML_doubleMOC_hires';
srcIter = 200*73;
% Nx_i = 1024;
% Ny_i = 512;
% destDir = fullfile(local_home_dir,'ACC_AABW_ML_doubleMOC_veryhires');
Nx_i = 256;
Ny_i = 128;
destDir = fullfile(local_home_dir,'ACC_AABW_ML_doubleMOC_lores');
destIter = 0;

%%% Load source simulation parameters
loadParams;

%%% Length of time over which to regrid forcing files, starting at
%%% srcIter
tlen_forc = 100*t1year;
tmin_forc = srcIter * dt;
tmax_forc = tmin_forc + tlen_forc;

%%% Check grid spacings match
dx_i = Lx/Nx_i;
dy_i = Ly/Ny_i;    

%%% Extended source grids that accommodate ghost points
[xx_h yy_h XX_h YY_h] = createmesh(-0.5*dx,Lx+0.5*dx,Nx+2,-0.5*dy,Ly+0.5*dy,Ny+2);      
[xx_u yy_u XX_u YY_u] = createmesh(-dx,Lx,Nx+2,-0.5*dy,Ly+0.5*dy,Ny+2);  
[xx_v yy_v XX_v YY_v] = createmesh(-0.5*dx,Lx+0.5*dx,Nx+2,-dy,Ly,Ny+2);
[xx_q yy_q XX_q YY_q] = createmesh(0,Lx,Nx+1,0,Ly,Ny+1);

%%% Destination grids  
[xx_h_i yy_h_i XX_h_i YY_h_i] = createmesh(0.5*dx_i,Lx-0.5*dx_i,Nx_i,0.5*dy_i,Ly-0.5*dy_i,Ny_i);      
[xx_u_i yy_u_i XX_u_i YY_u_i] = createmesh(0,Lx-dx_i,Nx_i,0.5*dy_i,Ly-0.5*dy_i,Ny_i);  
[xx_v_i yy_v_i XX_v_i YY_v_i] = createmesh(0.5*dx_i,Lx-0.5*dx_i,Nx_i,0,Ly-dy_i,Ny_i);    
[xx_q_i yy_q_i XX_q_i YY_q_i] = createmesh(0,Lx,Nx_i+1,0,Ly,Ny_i+1);

%%% Coriolis parameter
Omega_z = readDataFile(params_file,dirpath,'OmegazFile',Nx+1,Ny+1,zeros(Nx+1,Ny+1)); 
Omega_z_i = reshape(interp2(XX_q',YY_q',Omega_z',XX_q_i',YY_q_i','linear')',[1 Nx_i+1,Ny_i+1]);
fname = fullfile(destDir,'Omegaz.dat');
writeDataFile(fname,Omega_z_i);

%%% Bathymetry  
hhb = zeros(Nx+2,Ny+2);  
hhb(2:Nx+1,2:Ny+1)= readDataFile(params_file,dirpath,'hbFile',Nx,Ny,-1000*ones(Nx,Ny)); 
hhb(:,1) = hhb(:,2); %%% Reflection at N/S boundaries
hhb(:,Ny+2) = hhb(:,Ny+1);
hhb(1,:) = hhb(Nx+1,:); %%% Periodic E/W boundaries
hhb(Ny+2,:) = hhb(2,:);        
hhb_i = interp2(XX_h',YY_h',hhb',XX_h_i',YY_h_i','linear')';
fname = fullfile(destDir,'etab.dat');
writeDataFile(fname,hhb_i);

%%% Relaxation time scale (layer thickness)
hTime = zeros(Nx+2,Ny+2);  
hTime(2:Nx+1,2:Ny+1)= readDataFile(params_file,dirpath,'hTimeFile',Nx,Ny,-1*ones(Nx,Ny)); 
hTime(:,1) = hTime(:,2); %%% Reflection at N/S boundaries
hTime(:,Ny+2) = hTime(:,Ny+1);
hTime(1,:) = hTime(Nx+1,:); %%% Periodic E/W boundaries
hTime(Ny+2,:) = hTime(2,:);        
hTime_i = interp2(XX_h',YY_h',hTime',XX_h_i',YY_h_i','linear')';
fname = fullfile(destDir,'hTime.dat');
writeDataFile(fname,hTime_i);

%%% Relxation target (layer thickness)  
hRelax = readDataFile3D(params_file,dirpath,'hRelaxFile',Nlay,Nx,Ny,zeros(Nlay,Nx,Ny)); 
hRelax_i = zeros(Nlay,Nx_i,Ny_i);
for k=1:Nlay
  hRelax_rec = zeros(Nx+2,Ny+2);  
  hRelax_rec(2:Nx+1,2:Ny+1) = squeeze(hRelax(k,:,:));
  hRelax_rec(:,1) = hRelax_rec(:,2); %%% Reflection at N/S boundaries
  hRelax_rec(:,Ny+2) = hRelax_rec(:,Ny+1);
  hRelax_rec(1,:) = hRelax_rec(Nx+1,:); %%% Periodic E/W boundaries
  hRelax_rec(Ny+2,:) = hRelax_rec(2,:);        
  hRelax_i(k,:,:) = reshape(interp2(XX_h',YY_h',hRelax_rec',XX_h_i',YY_h_i','linear')',[1 Nx_i Ny_i]);
end
fname = fullfile(destDir,'hRelax.dat');
writeDataFile(fname,hRelax_i);

%%% Wind stress forcing
[tauPeriod tauPeriod_found] = readparam(params_file,'tauPeriod','%lf');
if (~tauPeriod_found)
  tauPeriod = 0;
end
[tauNrecs tauNrecs_found] = readparam(params_file,'tauNrecs','%lf');
if (~tauNrecs_found)
  tauNrecs = 1;
end
taux = readDataFile3D(params_file,dirpath,'tauxFile',tauNrecs,Nx,Ny,zeros(tauNrecs,Nx,Ny));   
tauTimes = (0:1:tauNrecs-1)/tauNrecs*tauPeriod;
tau_idx = find((tauTimes >= tmin_forc) & (tauTimes < tmax_forc));
taux = taux(tau_idx,:,:);
tauNrecs_i = length(tau_idx);
taux_i = zeros(tauNrecs_i,Nx_i,Ny_i);
for n=1:tauNrecs_i  
  taux_rec = zeros(Nx+2,Ny+2);
  taux_rec(2:Nx+1,2:Ny+1)= squeeze(taux(n,:,:));
  taux_rec(:,1) = - taux_rec(:,2); %%% 0 at N/S boundaries
  taux_rec(:,Ny+2) = - taux_rec(:,Ny+1);
  taux_rec(1,:) = taux_rec(Nx+1,:); %%% Periodic E/W boundaries
  taux_rec(Ny+2,:) = taux_rec(2,:);        
  taux_i(n,:,:) = reshape(interp2(XX_u',YY_u',taux_rec',XX_u_i',YY_u_i','linear')',[1 Nx_i,Ny_i]);
end  
fname = fullfile(destDir,'taux.dat');
writeDataFile(fname,taux_i);
clear('taux');
clear('taux_i');

%%% Diapycnal velocity
[wDiaPeriod wDiaPeriod_found] = readparam(params_file,'wDiaPeriod','%lf');  
[wDiaNrecs wDiaNrecs_found] = readparam(params_file,'wDiaNrecs','%lf');  
wDiaFile = 'wDiaFile.dat';
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
fclose(fid);
wDiaTimes = (0:1:wDiaNrecs-1)/wDiaNrecs*wDiaPeriod;
wDia_idx = find((wDiaTimes >= tmin_forc) & (wDiaTimes < tmax_forc));
wDia = wDia(wDia_idx,:,:,:);
wDiaNrecs_i = length(wDia_idx);
wDia_i = zeros(wDiaNrecs_i,Nlay+1,Nx_i,Ny_i);
m = 0;
for n=1:wDiaNrecs_i  
  for k=1:Nlay+1
    wDia_rec = zeros(Nx+2,Ny+2);
    wDia_rec(2:Nx+1,2:Ny+1)= squeeze(wDia(n,k,:,:));
    wDia_rec(:,1) = wDia_rec(:,2); %%% Reflection at N/S boundaries
    wDia_rec(:,Ny+2) = wDia_rec(:,Ny+1);
    wDia_rec(1,:) = wDia_rec(Nx+1,:); %%% Periodic E/W boundaries
    wDia_rec(Ny+2,:) = wDia_rec(2,:);        
    wDia_i(n,k,:,:) = reshape(interp2(XX_h',YY_h',wDia_rec',XX_h_i',YY_h_i','linear')',[1 1,Nx_i,Ny_i]);
  end
end  
fname = fullfile(destDir,'wDiaFile.dat');
writeDataFile(fname,wDia_i);
clear('wDia');
clear('wDia_i');

%%% Loop over layers to regrid model setate
uu = zeros(Nx+2,Ny+2);
vv = zeros(Nx+2,Ny+2);
hh = zeros(Nx+2,Ny+2);
for k=1:Nlay

  %%% Load model state in kth layer, allowing space for ghost points
  data_file = fullfile(dirpath,createOutputFilename(OUTN_U,srcIter,k));
  uu(2:Nx+1,2:Ny+1) = readOutputFile(data_file,Nx,Ny);
  data_file = fullfile(dirpath,createOutputFilename(OUTN_V,srcIter,k));
  vv(2:Nx+1,2:Ny+1) = readOutputFile(data_file,Nx,Ny);    
  data_file = fullfile(dirpath,createOutputFilename(OUTN_H,srcIter,k));
  hh(2:Nx+1,2:Ny+1) = readOutputFile(data_file,Nx,Ny);    

  %%% Fill in ghost points depending on boundary conditions
  if (useWallNS)
    uu(:,1) = uu(:,2);
    uu(:,Ny+2) = uu(:,Ny+1);
    vv(:,1) = vv(:,2);
    vv(:,Ny+2) = 0;
    hh(:,1) = hh(:,2);
    hh(:,Ny+2) = hh(:,Ny+1);
  else
    uu(:,1) = uu(:,Ny+1);
    uu(:,Ny+2) = uu(:,2);
    vv(:,1) = vv(:,Ny+1);
    vv(:,Ny+2) = vv(:,2);
    hh(:,1) = hh(:,Ny+1);
    hh(:,Ny+2) = hh(:,2);
  end    
  if (useWallEW)
    uu(1,:) = uu(2,:);
    uu(Nx+2,:) = 0;
    vv(1,:) = vv(2,:);
    vv(Nx+2,:) = vv(Nx+1,:);
    hh(1,:) = hh(2,:);
    hh(Nx+2,:) = hh(Nx+1,:);
  else
    uu(1,:) = uu(Nx+1,:);
    uu(Nx+2,:) = uu(2,:);
    vv(1,:) = vv(Nx+1,:);
    vv(Nx+2,:) = vv(2,:);
    hh(1,:) = hh(Nx+1,:);
    hh(Nx+2,:) = hh(2,:);
  end    

  %%% Do interpolation
  uu_i = interp2(XX_u',YY_u',uu',XX_u_i',YY_u_i','linear')';
  vv_i = interp2(XX_v',YY_v',vv',XX_v_i',YY_v_i','linear')';
  hh_i = interp2(XX_h',YY_h',hh',XX_h_i',YY_h_i','linear')';   

  %%% Enforce boundary conditions if required
  if (useWallNS)
    vv_i(:,1) = 0;
  end
  if (useWallEW)
    uu_i(1,:) = 0;
  end

  %%% Write output files
  fname = fullfile(destDir,createOutputFilename(OUTN_U,destIter,k));
  writeDataFile(fname,uu_i);
  fname = fullfile(destDir,createOutputFilename(OUTN_V,destIter,k));
  writeDataFile(fname,vv_i);
  fname = fullfile(destDir,createOutputFilename(OUTN_H,destIter,k));
  writeDataFile(fname,hh_i);

end
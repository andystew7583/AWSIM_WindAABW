%%%
%%% plotForcingMOCRelationships.m
%%%
%%% Examines relationships between.
%%%

%%% Make sure we have access to AWSIM tools
setPaths
constants;

%%% Options
run_name = 'ACC_AABW_ML_doubleMOC_hires'; %%% Run to analyze
tmin = 50*t1year; %%% Actual analysis period
tmax = 500*t1year;
Dsmooth = 73; %%% Interval between smoothed times. N.B. 73 snapshots/year
Nsmooth = 5*73; %%% Smoothing window width. N.B. 73 snapshots/year
Nchunks = 10;

%%% Load parameters   
local_home_dir = '/Volumes/Kilchoman/UCLA/Projects/AWSIM/runs';
prod_dir = fullfile(local_home_dir,'AWSIM_WindAABW_products');
loadParams;
dirpath = fullfile(local_home_dir,run_name);
gtild = reshape(cumsum(gg),[Nlay 1 1]);

%%% Load pre-computed transport time series
tt = [];
Tacc = [];
Taabw = [];
Taaiw = [];
for n=1:Nchunks
  
  %%% Load data for this chunk
  load(fullfile(prod_dir,[run_name,'_InstTrans_chunk',num2str(n),'.mat']), ...    
    'tauTimes','tauMax', ...
    'tt_chunk','Tacc_chunk','Taabw_chunk','Taaiw_chunk');
  tt = [tt tt_chunk];
  Tacc = [Tacc Tacc_chunk];
  Taabw = [Taabw Taabw_chunk];
  Taaiw = [Taaiw Taaiw_chunk];
    
end
Nt = size(Taabw,2);

%%% Time series of imposed diapycnal fluxes
[wDiaNrecs wDiaNrecs_found] = readparam(params_file,'wDiaNrecs','%lf'); 
[wDiaPeriod wDiaPeriod_found] = readparam(params_file,'wDiaPeriod','%lf'); 
wDiaFile = 'wDiaFile.dat';
wDiaInt = zeros(Nlay+1,wDiaNrecs);
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
  wDiaInt(:,n) = sum(sum(wDia(n,:,:,:)*dx*dy,3),4)';
end
fclose(fid);

%%% Forcing is periodic on a 500-year time scale - add additional points at
%%% the end to ensure linear interpolation proceeds correctly
tauTimes = [tauTimes 500*t1year];
tauMax = [tauMax tauMax(1)];
wDiaTimes = [wDiaTimes 500*t1year];
wDiaInt = [wDiaInt wDiaInt(:,1)];

%%% Compute forcing contemporaneous with MOC, as done in the model
tau = interp1(tauTimes,tauMax,tt,'linear');
Faaiw = interp1(wDiaTimes,wDiaInt(2,:),tt,'linear');
Faabw = interp1(wDiaTimes,wDiaInt(3,:),tt,'linear');

%%% Restrict to desired time range
range_idx = find(tt>=tmin & tt<=tmax);
tt = tt(range_idx);
Tacc = Tacc(range_idx);
Taaiw = Taaiw(range_idx);
Taabw = Taabw(range_idx);
tau = tau(range_idx);
Faaiw = Faaiw(range_idx);
Faabw = Faabw(range_idx);

%%% Compute coherence

[cxy,f] = mscohere(tau,Tacc,[],[],[],73)
figure(4);semilogx(f,cxy);
figure(5);plotyy(tt,smooth(tau,3000),tt,smooth(Tacc,3000));

[cxy,f] = mscohere(tau,Taabw,[],[],[],73)
figure(6);semilogx(f,cxy);
figure(7);plotyy(tt,smooth(tau,73),tt,smooth(Taabw,73));

[cxy,f] = mscohere(tau,Taaiw,[],[],[],73)
figure(8);semilogx(f,cxy);
figure(9);plotyy(tt,smooth(tau,73),tt,smooth(Taaiw,73));

[cxy,f] = mscohere(Faaiw,Taaiw,[],[],[],73)
figure(10);semilogx(f,cxy);
figure(11);plotyy(tt,smooth(Faaiw,3000),tt,smooth(Taaiw,3000));

[cxy,f] = mscohere(Faaiw,Taabw,[],[],[],73)
figure(12);semilogx(f,cxy);
figure(13);plotyy(tt,smooth(Faaiw,3000),tt,smooth(Taabw,3000));

[cxy,f] = mscohere(Faaiw,Tacc,[],[],[],73)
figure(14);semilogx(f,cxy);
figure(15);plotyy(tt,smooth(Faaiw,3000),tt,smooth(Tacc,3000));

[cxy,f] = mscohere(Faabw,Taaiw,[],[],[],73)
figure(16);semilogx(f,cxy);
figure(17);plotyy(tt,smoothdata(Faabw,'gauss',3000),tt,smoothdata(Taaiw,'gauss',3000));

[cxy,f] = mscohere(Faabw,Taabw,[],[],[],73)
figure(18);semilogx(f,cxy);
figure(19);plotyy(tt,smoothdata(Faabw,'gauss',3000),tt,smoothdata(Taabw,'gauss',3000));

[cxy,f] = mscohere(Faabw,Tacc,[],[],[],73)
figure(20);semilogx(f,cxy);
figure(21);plotyy(tt,smoothdata(Faabw,'gauss',3000),tt,smoothdata(Tacc,'gauss',3000));

%%%
%%% plotMOChovmoeller.m
%%%
%%% Plots hovmoeller diagram of MOC strength.
%%%

%%% Options
run_name = 'ACC_AABW_ML_doubleMOC'; %%% Run to analyze
% run_name = 'ACC_AABW_ML_randWdia_randTau_white'; %%% Run to analyze
% run_name = 'ACC_AABW_ML_randWdia_randTau_white_Nlay2'; %%% Run to analyze
tmin = 100*t1year; %%% Actual analysis period
tmax = 500*t1year;
Dsmooth = 73; %%% Interval between smoothed times. N.B. 73 snapshots/year
Nsmooth = 20*73; %%% Smoothing window width. N.B. 73 snapshots/year
Ntimes = length(tt);

%%% Load parameters   
% % local_home_dir = '/Volumes/Kilchoman/UCLA/Projects/AWSIM/runs';
local_home_dir = '/Volumes/Stewart-RAID1-A/UCLA/Projects/AWSIM/runs';
% local_home_dir = '/data2/astewart/AWSIM/runs';
prod_dir = fullfile(local_home_dir,'AWSIM_WindAABW_products');
loadParams;
dirpath = fullfile(local_home_dir,run_name);
gtild = reshape(cumsum(gg),[Nlay 1 1]);

%%% Load pre-computed transport time series
load(fullfile(prod_dir,[run_name,'_InstTrans.mat']));

%%% Extract transports
Taaiw = squeeze(hv_int(1,:,:));
Taabw = squeeze(hv_int(3,:,:));


%%% Smoothing
Nsmooth2 = round(Nsmooth/2);
Dsmooth2 = round(Dsmooth/2);
idx_smooth = round(Dsmooth/2:Dsmooth:Ntimes-Dsmooth/2);
tt_smooth = NaN*ones(1,length(idx_smooth));
Taabw_smooth = NaN*ones(Ny,length(idx_smooth));
Taaiw_smooth = NaN*ones(Ny,length(idx_smooth));
for m=1:length(idx_smooth)
  m
  if ((idx_smooth(m) < Nsmooth2+1) || (idx_smooth(m)>Ntimes-Nsmooth2))
    continue;
  end    
  tt_smooth(m) = mean(tt(idx_smooth(m)-Nsmooth2:idx_smooth(m)+Nsmooth2));
  Taabw_smooth(:,m) = mean(Taabw(:,idx_smooth(m)-Nsmooth2:idx_smooth(m)+Nsmooth2),2);  
  Taaiw_smooth(:,m) = mean(Taaiw(:,idx_smooth(m)-Nsmooth2:idx_smooth(m)+Nsmooth2),2);  
end

tidx = (tt_smooth > tmin) & (tt_smooth < tmax);
Taabw_anom = Taabw_smooth-repmat(nanmean(Taabw_smooth(:,tidx),2),[1 length(tt_smooth)]);
Taaiw_anom = Taaiw_smooth-repmat(nanmean(Taaiw_smooth(:,tidx),2),[1 length(tt_smooth)]);

[TT,YY] = meshgrid(tt_smooth,yy_h);
figure(1);
pcolor(YY,TT/t1year,Taaiw_anom/1e6);
shading interp;
colorbar;
colormap redblue;
caxis([-1.5 1.5]);


[TT,YY] = meshgrid(tt_smooth,yy_h);
figure(2);
pcolor(YY,TT/t1year,Taaiw_anom/1e6);
shading interp;
colorbar;
colormap redblue;
caxis([-1.5 1.5]);


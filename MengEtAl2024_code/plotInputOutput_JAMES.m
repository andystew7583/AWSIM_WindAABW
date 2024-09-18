%%%
%%% plotInputOutput_JAMES.m
%%%
%%% Plots wind forcing, instantaneous model state and bathymetry for our JAMES paper.
%%%

%%% Run to load
run_name = 'ACC_AABW_ML_doubleMOC_hires';
tmin = 0*t1year; %%% Actual analysis period
tmax = 500*t1year;
Nchunks = 10;

%%% Load parameters   
local_home_dir = '/Volumes/Kilchoman/UCLA/Projects/AWSIM/runs';
% local_home_dir = '/data2/astewart/AWSIM/runs';
prod_dir = fullfile(local_home_dir,'AWSIM_WindAABW_products');
loadParams;
dirpath = fullfile(local_home_dir,run_name);
gtild = reshape(cumsum(gg),[Nlay 1 1]);

%%% To store data
ssh_map = [];
obp_map = [];
tt_full = [];
Tacc_full = [];
Taabw_full = [];
Taaiw_full = [];

%%% Load pre-computed transport time series
for n=1:Nchunks
  n
  load(fullfile(prod_dir,[run_name,'_InstTrans_chunk',num2str(n),'.mat']),'Taabw_chunk','Taaiw_chunk','Tacc_chunk','tt_chunk');
  tt_full = [tt_full tt_chunk];
  Tacc_full = [Tacc_full Tacc_chunk];
  Taabw_full = [Taabw_full Taabw_chunk];
  Taaiw_full = [Taaiw_full Taaiw_chunk];
  if n == 2
    load(fullfile(prod_dir,[run_name,'_InstTrans_chunk',num2str(n),'.mat']),'obp_chunk','ssh_chunk');
    ssh_map = ssh_chunk(:,:,100);
    obp_map = obp_chunk(:,:,100);
  end
end




%%% Plotting options
scrsz = get(0,'ScreenSize');
framepos = [209   632  900   900];
fontsize = 14;
ax_pos = zeros(3,4);
ax_pos(1,:) = [0.06 0.54 0.9 0.46];
ax_pos(2,:) = [0.07 0.36 0.39 0.1];
ax_pos(3,:) = [0.07 0.06 0.39 0.26];
ax_pos(4,:) = [0.56 0.36 0.4 0.1];
ax_pos(5,:) = [0.56 0.21 0.4 0.1];
ax_pos(6,:) = [0.56 0.06 0.4 0.1];
cb_pos = [0.92 0.55 0.015 0.4];
lab_size = [0.05 0.03];
colororder = get(gca,'ColorOrder');
omega_color = colororder(7,:);
ACC_color = colororder(5,:);
MOC_color = colororder(1,:);

%%% Set up the frame
figure(202);
clf; 
set(gcf,'Position',framepos);
set(gcf,'Color','w');
cntr = 0;


plot(tt,Taabw);
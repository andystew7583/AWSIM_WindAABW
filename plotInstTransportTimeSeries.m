%%%
%%% plotInstTransportTimeSeries.m
%%%
%%% Plots time series of AABW export and looks at correlations/regressions.
%%%

%%% Options
run_name = 'ACC_AABW_ML_randWdia_randTau_white_Nlay2'; %%% Run to analyze
tmin = 50*t1year; %%% Actual analysis period
tmax = 300*t1year;
Dsmooth = 73; %%% Interval between smoothed times. N.B. 73 snapshots/year
Nsmooth = 5*73; %%% Smoothing window width. N.B. 73 snapshots/year

%%% Load parameters   
local_home_dir = '/Volumes/Kilchoman/UCLA/Projects/AWSIM/runs';
prod_dir = fullfile(local_home_dir,'AWSIM_WindAABW_products');
loadParams;
dirpath = fullfile(local_home_dir,run_name);
gtild = reshape(cumsum(gg),[Nlay 1 1]);

%%% Load pre-computed transport time series
load(fullfile(prod_dir,[run_name,'_InstTrans.mat']));

%%% Compute meridional average
yminidx = find(yy_v>100*m1km,1,'first');
ymaxidx = find(yy_v>Ly-100*m1km,1,'first');
hv_yint = squeeze(mean(hv_int(:,yminidx:ymaxidx,:),2));
hv_mean_yint = squeeze(mean(hv_mean_int(:,yminidx:ymaxidx,:),2));
hv_stand_yint = squeeze(mean(hv_stand_int(:,yminidx:ymaxidx,:),2));

%%% Just plot the transports
for m=1:Nlay
  figure(1+m)
  plot(tt(:)/t1year,hv_mean_yint(m,:));
  hold on
  plot(tt(:)/t1year,hv_stand_yint(m,:));
  plot(tt(:)/t1year,hv_yint(m,:));
  hold off
  legend('Mean','Standing','Total');
  title(['Layer ',num2str(m)]);
end

%%% Define AABW time series
Taabw = hv_yint(Nlay,:);
Nt = size(Taabw,2);

%%% Geostrophic reconstruction. Needs to be computed for each snapshot to
%%% capture the eddy contribution. Should only work in 2 layers.
Tdyn = zeros(1,Nt);
ff = 2*Omega_z(1,1:Ny);  
for n=1:Nt
  n  
  h1_r = (g*ssh(:,:,n) - obp(:,:,n) - (gg(1)+gg(2))*hhb) / gg(2);
%   h2_r = (obp(:,:,n) + g*hhb - g*ssh(:,:,n)) / gg(2);
  M1_r = ssh(:,:,n)*g;
%   M2_r = obp(:,:,n) + (gg(1)+gg(2))*hhb; %%% Strictly 2-layer version
  Mbot_r = obp(:,:,n) + gtild(Nlay)*hhb; %%% Exact for bottom isopycnal M in any number of layers
  if (Nlay == 2)
    hbot_r = (obp(:,:,n) + g*hhb - g*ssh(:,:,n)) / gg(2); %%% Estimate for bottom isopycnal h in 2 layers
  end
  if (Nlay == 3)
    H1 = 750+500*YY_h/Ly;
    hbot_r = (obp(:,:,n) + (gg(1)+gg(2))*hhb - g*ssh(:,:,n) + gg(2)*H1) / gg(3); %%% Layer thickness estimates not bad but transport way too high
%     hbot_r = (obp(:,:,n) + gg(1)*hhb - g*ssh(:,:,n)) / gg(3); %%%% Layer thickness estimates way off but transport spot on.
  end
  hv_xint_r = squeeze( sum( ...
     0.5*(hbot_r(1:Nx,:)+hbot_r([Nx 1:Nx-1],:)) ...
  .* (Mbot_r(1:Nx,:)-Mbot_r([Nx 1:Nx-1],:))./ff ...
  ,1));
  Tdyn(n) = mean(hv_xint_r(yminidx:ymaxidx));
end


%%% Smoothing
Nsmooth2 = round(Nsmooth/2);
Dsmooth2 = round(Dsmooth/2);
idx_smooth = round(Dsmooth/2:Dsmooth:Nt-Dsmooth/2);
tt_smooth = NaN*ones(size(idx_smooth));
Taabw_smooth = NaN*ones(size(idx_smooth));
Tdyn_smooth = 0*idx_smooth;
ssh_smooth = NaN*ones(Nx,Ny,length(idx_smooth));
obp_smooth = NaN*ones(Nx,Ny,length(idx_smooth));
for m=1:length(idx_smooth)
  m
  if ((idx_smooth(m) < Nsmooth2+1) || (idx_smooth(m)>Nt-Nsmooth2))
    continue;
  end    
  tt_smooth(m) = mean(tt(idx_smooth(m)-Nsmooth2:idx_smooth(m)+Nsmooth2));
  Taabw_smooth(m) = mean(Taabw(idx_smooth(m)-Nsmooth2:idx_smooth(m)+Nsmooth2));
  ssh_smooth(:,:,m) = mean(ssh(:,:,idx_smooth(m)-Nsmooth2:idx_smooth(m)+Nsmooth2),3);
  obp_smooth(:,:,m) = mean(obp(:,:,idx_smooth(m)-Nsmooth2:idx_smooth(m)+Nsmooth2),3);
  Tdyn_smooth(m) = mean(Tdyn(idx_smooth(m)-Nsmooth2:idx_smooth(m)+Nsmooth2));  
end


%%% Limits on data range

tminidx = find(tt_smooth>tmin,1,'first');
tmaxidx = find(tt_smooth>tmax,1,'first');
tt_smooth = tt_smooth(tminidx:tmaxidx);
ssh_smooth = ssh_smooth(:,:,tminidx:tmaxidx);
obp_smooth = obp_smooth(:,:,tminidx:tmaxidx);
Taabw_smooth = Taabw_smooth(tminidx:tmaxidx);
Tdyn_smooth = Tdyn_smooth(tminidx:tmaxidx);
Nt_smooth = length(tt_smooth);

tminidx = find(tt>tmin,1,'first');
tmaxidx = find(tt>tmax,1,'first');
tt = tt(tminidx:tmaxidx);
ssh = ssh(:,:,tminidx:tmaxidx);
obp = obp(:,:,tminidx:tmaxidx);
Taabw = Taabw(tminidx:tmaxidx);
Tdyn = Tdyn(tminidx:tmaxidx);
Nt = length(tt);


%%% Pointwise correlations
sshcorr = zeros(Nx,Ny);
obpcorr = zeros(Nx,Ny);
sshobpcorr = zeros(Nx,Ny);
for i=1:Nx
  i
  for j=1:Ny
    sshcorr(i,j) = corr(squeeze(ssh_smooth(i,j,:)),Taabw_smooth'); 
    obpcorr(i,j) = corr(squeeze(obp_smooth(i,j,:)),Taabw_smooth');  
    sshobpcorr(i,j) = corr(squeeze(ssh_smooth(i,j,:)),squeeze(obp_smooth(i,j,:)));
  end
end

figure(7);
pcolor(XX_h,YY_h,sshcorr)
shading interp
colorbar

figure(8);
pcolor(XX_h,YY_h,obpcorr)
shading interp
colorbar

figure(9);
pcolor(XX_h,YY_h,sshobpcorr)
shading interp
colorbar

figure(10);
pcolor(XX_h,YY_h,std(ssh_smooth,0,3))
shading interp
colorbar

figure(11);
pcolor(XX_h,YY_h,std(obp_smooth,0,3))
shading interp
colorbar





%%%%%%%%%%%%%%%%%%
%%% Regression %%%
%%%%%%%%%%%%%%%%%%

%%% First compute means and anomalies
sshmean = mean(ssh_smooth,3);
obpmean = mean(obp_smooth,3);
sshanom = ssh_smooth - repmat(sshmean,[1 1 Nt_smooth]);
obpanom = obp_smooth - repmat(obpmean,[1 1 Nt_smooth]);
Tmean = mean(Taabw_smooth);
Tanom = Taabw_smooth - Tmean;

%%% Compute regressions of SSH and OBP on MOC
sshreg = zeros(Nx,Ny);
for i=1:Nx
  i
  for j=1:Ny
    sshreg(i,j) = reshape(Tanom,[Nt_smooth 1])\squeeze(sshanom(i,j,:));
  end
end
obpreg = zeros(Nx,Ny);
for i=1:Nx
  i
  for j=1:Ny
    obpreg(i,j) = reshape(Tanom,[Nt_smooth 1])\squeeze(obpanom(i,j,:));
  end
end

figure(12);
pcolor(XX_h,YY_h,sshreg*1e6)
shading interp
colorbar

figure(13);
pcolor(XX_h,YY_h,obpreg./mean(obp,3)*1e6)
shading interp
colorbar





%%% Reconstruct transport time series via linear regression
Trec_ssh = zeros(1,Nt_smooth);
for n=1:Nt_smooth
  Trec_ssh(n) = reshape(sshreg,1,[])' \ reshape(sshanom(:,:,n),1,[])';
end
Trec_obp = zeros(1,Nt_smooth);
for n=1:Nt_smooth
  Trec_obp(n) = reshape(obpreg./obpmean,1,[])' \ reshape(obpanom(:,:,n)./obpmean,1,[])';
end

[r,p] = corr(Trec_ssh',Tanom')

figure(14);
plot(tt_smooth/t1year,Trec_ssh);
hold on;
plot(tt_smooth/t1year,Tanom);
hold off;
legend('T_a_n_o_m','T_s_s_h');
title(['r = ',num2str(r)]);

[r,p] = corr(Trec_obp',Tanom')

figure(15);
plot(tt_smooth/t1year,Tanom);
hold on;
plot(tt_smooth/t1year,Trec_obp);
hold off;
legend('T_a_n_o_m','T_o_b_p');
title(['r = ',num2str(r)]);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Dynamical reconstruction %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SIFS = zeros(Ny,Nt_smooth);
for n=1:Nt_smooth
  SIFS(:,n) = squeeze( sum( ...
     0.5*(ssh_smooth(1:Nx,:,n)+ssh_smooth([Nx 1:Nx-1],:,n)) ...
  .* (hhb(1:Nx,:)-hhb([Nx 1:Nx-1],:)) ...
  ,1));
end
TFS = zeros(Ny,Nt_smooth);
for n=1:Nt_smooth
  TFS(:,n) = squeeze( sum( ...
     0.5*(obp_smooth(1:Nx,:,n)+obp_smooth([Nx 1:Nx-1],:,n)) ...
  .* (hhb(1:Nx,:)-hhb([Nx 1:Nx-1],:)) ...
  ,1));
end

SIFS_xyavg = mean(SIFS(yminidx:ymaxidx,:))/Lx*1000;

[r,p] = corr(SIFS_xyavg',Taabw_smooth');

figure(16);
scatter(SIFS_xyavg,Taabw_smooth/1e6);
xlabel('SIFS (N/m^2)');
ylabel('T_A_A_B_W (smoothed)');
title(['r = ',num2str(r)]);



[r,p] = corr(Tdyn_smooth',Taabw_smooth');

figure(17);
plot(tt_smooth/t1year,Taabw_smooth/1e6)
hold on;
plot(tt_smooth/t1year,Tdyn_smooth/1e6)
hold off;
legend('T_A_A_B_W','T_d_y_n');
title(['r = ',num2str(r)]);


[r,p] = corr(Tdyn',Taabw');

figure(18);
plot(tt/t1year,Taabw/1e6)
hold on;
plot(tt/t1year,Tdyn/1e6)
hold off;
legend('T_A_A_B_W','T_d_y_n');
title(['r = ',num2str(r)]);






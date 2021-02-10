
%%% Load parameters   
% local_home_dir = '/Volumes/LaCie/UCLA/Projects/AWSIM_WindAABW/runs';
% run_name = 'test_run';
local_home_dir = '/Volumes/Kilchoman/UCLA/Projects/AWSIM/runs';
run_name = 'ACC_AABW_ML_randWdia_Nlay2';
loadParams;
dirpath = fullfile(local_home_dir,run_name);
gtild = reshape(cumsum(gg),[Nlay 1 1]);


%%% Time series of "true" AABW export
[wDiaNrecs h0_found] = readparam(params_file,'wDiaNrecs','%lf'); 
[wDiaPeriod h0_found] = readparam(params_file,'wDiaPeriod','%lf'); 
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


% minidx = find(yy_v>900*m1km,1,'first');
% maxidx = find(yy_v>Ly-900*m1km,1,'last');
tend = 45*t1year;

%%% At each time iteration...
cntr = 0;
tt = zeros(1,Nframes);
hv_int = NaN*ones(Nlay,Ny,Nframes);
hv_mean_int = NaN*ones(Nlay,Ny,Nframes);
hv_stand_int = NaN*ones(Nlay,Ny,Nframes);
h = zeros(Nlay,Nx,Ny);
v = zeros(Nlay,Nx,Ny);
ssh = zeros(Nx,Ny,Nframes);
obp = zeros(Nx,Ny,Nframes);
for n=n0+1:1:n0+Nframes-1   

  %%% Current simulation time    
  t = startTime + (n-n0)*dt_s;
  
  if (t>tend)
    break;
  end
  
  cntr = cntr + 1
  
  tt(cntr) = t;
 
  %%% Load instantaneous model state
  for k=1:Nlay
    data_file = fullfile(dirpath,[OUTN_H,num2str(k-1),'_n=',num2str(n),'.dat']);
    h(k,:,:) = readOutputFile(data_file,Nx,Ny);
  end  
  for k=1:Nlay
    data_file = fullfile(dirpath,[OUTN_V,num2str(k-1),'_n=',num2str(n),'.dat']);
    v(k,:,:) = readOutputFile(data_file,Nx,Ny);
  end
  h(:,:,2:Ny) = 0.5*(h(:,:,1:Ny-1)+h(:,:,2:Ny)); %%% Interpolate to v-points
  data_file = fullfile(dirpath,[OUTN_PI,'_n=',num2str(n),'.dat']);
  pi = readOutputFile(data_file,Nx,Ny);
  
  %%% Calculate sea surface height and ocean bottom pressure
  ssh(:,:,cntr) = pi/g;
  obp(:,:,cntr) = pi + squeeze(sum(gtild.*h,1));
  
  %%% Meridional transport decomposition
  hv_int(:,:,cntr) = squeeze(sum(h.*v.*dx,2));
  hv_mean_int(:,:,cntr) = squeeze(mean(h,2) .* mean(v,2) * Lx);
  hv_stand_int(:,:,cntr) = hv_int(:,:,cntr) - hv_mean_int(:,:,cntr);  
  
end

%%% Remove missing data
tt = tt(1:cntr);
hv_int = hv_int(:,:,1:cntr);
hv_mean_int = hv_mean_int(:,:,1:cntr);
hv_stand_int = hv_stand_int(:,:,1:cntr);
ssh = ssh(:,:,1:cntr);
obp = obp(:,:,1:cntr);

%%% Limits on data range
tminidx = find(tt>50*t1year,1,'first');
tmaxidx = find(tt>200*t1year,1,'first');

%%% Compute meridional average
yminidx = find(yy_v>100*m1km,1,'first');
ymaxidx = find(yy_v>Ly-100*m1km,1,'first');
hv_yint = squeeze(mean(hv_int(:,yminidx:ymaxidx,:),2));
hv_mean_yint = squeeze(mean(hv_mean_int(:,yminidx:ymaxidx,:),2));
hv_stand_yint = squeeze(mean(hv_stand_int(:,yminidx:ymaxidx,:),2));


figure(2);
plot(tt(:)/t1year,hv_mean_yint(1,:));
hold on
plot(tt(:)/t1year,hv_stand_yint(1,:));
plot(tt(:)/t1year,hv_yint(1,:));
hold off
legend('Mean','Standing','Total');
title('Layer 1');

figure(3);
plot(tt(:)/t1year,hv_mean_yint(2,:));
hold on
plot(tt(:)/t1year,hv_stand_yint(2,:));
plot(tt(:)/t1year,hv_yint(2,:));
hold off
legend('Mean','Standing','Total');
title('Layer 2');

figure(4);
plot(tt(:)/t1year,hv_mean_yint(3,:));
hold on
plot(tt(:)/t1year,hv_stand_yint(3,:));
plot(tt(:)/t1year,hv_yint(3,:));
hold off
legend('Mean','Standing','Total');
title('Layer 3');

%%% FFT of layer three mean transport
Tfft = abs(fft(hv_yint(3,tminidx:tmaxidx)));
Tfft(1) = 0;
figure(5);
plot(Tfft(1:length(Tfft)/2))

%%% Correlation between transports at different latitudes
figure(6);
plot(yy_h,(corr(squeeze(hv_int(3,:,tminidx:tmaxidx))',squeeze(hv_int(3,Ny/4,tminidx:tmaxidx)))))

%%% Pointwise correlations
sshcorr = zeros(Nx,Ny);
for i=1:Nx
  i
  for j=1:Ny
    sshcorr(i,j) = corr(squeeze(ssh(i,j,tminidx:tmaxidx)),hv_yint(3,tminidx:tmaxidx)');
  end
end
obpcorr = zeros(Nx,Ny);
for i=1:Nx
  i
  for j=1:Ny
    obpcorr(i,j) = corr(squeeze(obp(i,j,tminidx:tmaxidx)),hv_yint(3,tminidx:tmaxidx)');
  end
end
sshobpcorr = zeros(Nx,Ny);
for i=1:Nx
  i
  for j=1:Ny
    sshobpcorr(i,j) = corr(squeeze(ssh(i,j,tminidx:tmaxidx)),squeeze(obp(i,j,tminidx:tmaxidx)));
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
pcolor(XX_h,YY_h,std(ssh(:,:,tminidx:tmaxidx),0,3))
shading interp
colorbar

figure(11);
pcolor(XX_h,YY_h,std(obp(:,:,tminidx:tmaxidx),0,3))
shading interp
colorbar





%%%%%%%%%%%%%%%%%%
%%% Regression %%%
%%%%%%%%%%%%%%%%%%

%%% First compute means and anomalies
sshmean = mean(ssh(:,:,tminidx:tmaxidx),3);
obpmean = mean(obp(:,:,tminidx:tmaxidx),3);
sshanom = ssh - repmat(sshmean,[1 1 cntr]);
obpanom = obp - repmat(obpmean,[1 1 cntr]);
Tmean = mean(hv_yint(3,tminidx:tmaxidx));
Tanom = hv_yint(3,:) - Tmean;

%%% Take running averages
W_win = 1; %%% N.B. 73 snapshots/year
Ntot = length(tminidx:tmaxidx);
Nwin = floor(Ntot/W_win);
ssh_win = zeros(Nx,Ny,Nwin);
obp_win = zeros(Nx,Ny,Nwin);
T_win = zeros(1,Nwin);
t_win = zeros(1,Nwin);
for m=1:Nwin
  ssh_win(:,:,m) = mean(sshanom(:,:,tminidx+(m-1)*W_win:tminidx+m*W_win-1),3);
  obp_win(:,:,m) = mean(obpanom(:,:,tminidx+(m-1)*W_win:tminidx+m*W_win-1),3);
  T_win(m) = mean(Tanom(tminidx+(m-1)*W_win:tminidx+m*W_win-1));
  t_win(m) = mean(tt(tminidx+(m-1)*W_win:tminidx+m*W_win-1));
end

%%% Compute regressions of SSH and OBP on MOC
sshreg = zeros(Nx,Ny);
for i=1:Nx
  i
  for j=1:Ny
    sshreg(i,j) = squeeze(T_win')\squeeze(ssh_win(i,j,:));
  end
end
obpreg = zeros(Nx,Ny);
for i=1:Nx
  i
  for j=1:Ny
    obpreg(i,j) = squeeze(T_win')\squeeze(obp_win(i,j,:));
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


%%%%%%%%%%%%%%%%%%
%%% Regression %%%
%%%%%%%%%%%%%%%%%%

%%% First compute means and anomalies
sshmean = mean(ssh(:,:,tminidx:tmaxidx),3);
obpmean = mean(obp(:,:,tminidx:tmaxidx),3);
sshanom = ssh - repmat(sshmean,[1 1 cntr]);
obpanom = obp - repmat(obpmean,[1 1 cntr]);
Tmean = mean(hv_yint(3,tminidx:tmaxidx));
Tanom = hv_yint(3,:) - Tmean;




%%% Take running averages
W_win = 1; %%% N.B. 73 snapshots/year
Ntot = length(tminidx:tmaxidx);
Nwin = floor(Ntot/W_win);
ssh_win = zeros(Nx,Ny,Nwin);
obp_win = zeros(Nx,Ny,Nwin);
T_win = zeros(1,Nwin);
t_win = zeros(1,Nwin);
for m=1:Nwin
  ssh_win(:,:,m) = mean(sshanom(:,:,tminidx+(m-1)*W_win:tminidx+m*W_win-1),3);
  obp_win(:,:,m) = mean(obpanom(:,:,tminidx+(m-1)*W_win:tminidx+m*W_win-1),3);
  T_win(m) = mean(Tanom(tminidx+(m-1)*W_win:tminidx+m*W_win-1));
  t_win(m) = mean(tt(tminidx+(m-1)*W_win:tminidx+m*W_win-1));
end




%%% Compute regressions of SSH and OBP on MOC
sshreg = zeros(Nx,Ny);
for i=1:Nx
  i
  for j=1:Ny
    sshreg(i,j) = squeeze(T_win')\squeeze(ssh_win(i,j,:));
  end
end
obpreg = zeros(Nx,Ny);
for i=1:Nx
  i
  for j=1:Ny
    obpreg(i,j) = squeeze(T_win')\squeeze(obp_win(i,j,:));
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
Trec_ssh = zeros(1,Nwin);
for n=1:Nwin
  Trec_ssh(n) = reshape(sshreg,1,[])' \ reshape(sshanom(:,:,n),1,[])';
end
Trec_obp = zeros(1,Nwin);
for n=1:Nwin
  Trec_obp(n) = reshape(obpreg./obpmean,1,[])' \ reshape(obpanom(:,:,n)./obpmean,1,[])';
end

figure(14);
plot(t_win/t1year,Trec_ssh);

figure(15);
plot(t_win/t1year,Trec_obp);




%%% TODO multilinear regression?

%%% TODO can we dynamically reconstruct the deep geostrophic transport? in
%%% 2 layers yes, but 3 layers might be trickier

%%% TODO how about correlation with surface-induced topographic form stress

%%% TODO try regression on SSP-OBP difference? Approximates thickness of
%%% lowest layer?

%%% TODO run a case with actual MOC variability
%%% TODO run a 2-layer case and predict MOC dynamically




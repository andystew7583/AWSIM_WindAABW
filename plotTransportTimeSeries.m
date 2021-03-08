
%%% Load parameters   
% local_home_dir = '/Volumes/LaCie/UCLA/Projects/AWSIM_WindAABW/runs';
% run_name = 'test_run';
local_home_dir = '/Volumes/Kilchoman/UCLA/Projects/AWSIM/runs';
run_name = 'test_wdia_long_N256';
loadParams;
dirpath = fullfile(local_home_dir,run_name);

layer = 1;
% minidx = find(yy_v>100*m1km,1,'first');
% maxidx = find(yy_v>Ly-100*m1km,1,'last');
minidx = find(yy_v>900*m1km,1,'first');
maxidx = find(yy_v>Ly-900*m1km,1,'last');

%%% At each time iteration...
cntr = 1;
tt = zeros(1,N_avg_hu);
taux_int = NaN*ones(1,N_avg_hu);
hv_int = NaN*ones(3,N_avg_hu);
hv_mean_int = NaN*ones(3,N_avg_hu);
hv_stand_int = NaN*ones(3,N_avg_hu);
hv_eddy_int = NaN*ones(3,N_avg_hu);
hv = zeros(Nlay,Nx,Ny);
h = zeros(Nlay,Nx,Ny);
v = zeros(Nlay,Nx,Ny);
for n=n0_avg_hu+1:1:n0_avg_hu+N_avg_hu-1   

  %%% Current simulation time    
  t = startTime + (n-n0_avg_hu)*dt_avg_hu;
  tt(cntr) = t;

  data_file = fullfile(dirpath,[OUTN_UMOM_WIND,num2str(layer-1),'_n=',num2str(n),'.dat']);
  taux = readOutputFile(data_file,Nx,Ny);
  if (isempty(taux))
    break;
  end
  taux_int(cntr) = sum(sum(taux*dx*dy));
  
  for k=1:Nlay
    data_file = fullfile(dirpath,[OUTN_HV_AVG,num2str(k-1),'_n=',num2str(n),'.dat']);
    hv(k,:,:) = readOutputFile(data_file,Nx,Ny);
  end
  hv_int(:,cntr) = mean(sum(hv(:,:,minidx:maxidx)*dx,2),3);
%   hv_int(:,cntr) = min(sum(hv(:,:,minidx:maxidx)*dx,2),[],3);
  
  for k=1:Nlay
    data_file = fullfile(dirpath,[OUTN_H_AVG,num2str(k-1),'_n=',num2str(n),'.dat']);
    h(k,:,:) = readOutputFile(data_file,Nx,Ny);
  end  
  for k=1:Nlay
    data_file = fullfile(dirpath,[OUTN_V_AVG,num2str(k-1),'_n=',num2str(n),'.dat']);
    v(k,:,:) = readOutputFile(data_file,Nx,Ny);
  end

  hv_mean_int(:,cntr) = mean(mean(h(:,:,minidx:maxidx),2).*mean(v(:,:,minidx:maxidx),2)*Lx,3);
  hv_stand_int(:,cntr) = mean(sum(h(:,:,minidx:maxidx).*v(:,:,minidx:maxidx)*dx,2),3) - hv_mean_int(:,cntr);
  hv_eddy_int(:,cntr) = hv_int(:,cntr) - hv_mean_int(:,cntr) - hv_stand_int(:,cntr);
  
  cntr = cntr + 1;
  
end

figure(1);
minidx = find(tt>2*t1year,1);
maxidx = find(~isnan(taux_int),1,'last'); %%% N.B> maxidx redefined here
plot(tt(1:maxidx)/t1year,taux_int(1:maxidx));

figure(2);
plot(tt(1:maxidx)/t1year,hv_mean_int(1,1:maxidx));
hold on
plot(tt(1:maxidx)/t1year,hv_stand_int(1,1:maxidx));
plot(tt(1:maxidx)/t1year,hv_eddy_int(1,1:maxidx));
plot(tt(1:maxidx)/t1year,hv_int(1,1:maxidx));
hold off
legend('Mean','Standing','Eddy','Total');

figure(3);
plot(tt(1:maxidx)/t1year,hv_mean_int(2,1:maxidx));
hold on
plot(tt(1:maxidx)/t1year,hv_stand_int(2,1:maxidx));
plot(tt(1:maxidx)/t1year,hv_eddy_int(2,1:maxidx));
plot(tt(1:maxidx)/t1year,hv_int(2,1:maxidx));
hold off
legend('Mean','Standing','Eddy','Total');

figure(4);
plot(tt(1:maxidx)/t1year,hv_mean_int(3,1:maxidx));
hold on
plot(tt(1:maxidx)/t1year,hv_stand_int(3,1:maxidx));
plot(tt(1:maxidx)/t1year,hv_eddy_int(3,1:maxidx));
plot(tt(1:maxidx)/t1year,hv_int(3,1:maxidx));
hold off
legend('Mean','Standing','Eddy','Total');

figure(5);
lag = 2;
scatter(taux_int(minidx+lag:maxidx),hv_int(3,minidx:maxidx-lag))
[r,p] = corr(taux_int(minidx+lag:maxidx)',hv_int(3,minidx:maxidx-lag)')

figure(6)
lag = 0;
scatter(hv_mean_int(1,minidx+lag:maxidx),hv_int(2,minidx:maxidx-lag))
[r,p] = corr(hv_mean_int(1,minidx+lag:maxidx)',hv_int(2,minidx:maxidx-lag)')

figure(7)
lag = 2;
scatter(hv_mean_int(1,minidx+lag:maxidx),hv_int(3,minidx:maxidx-lag))
[r,p] = corr(hv_mean_int(1,minidx+lag:maxidx)',hv_int(3,minidx:maxidx-lag)')
  
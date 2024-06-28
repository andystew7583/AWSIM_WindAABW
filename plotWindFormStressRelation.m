%%%
%%% plotWindFormStressRelation.m
%%%
%%% Compares form stress variability with imposed wind forcing.
%%%

%%% Options
run_name = 'ACC_AABW_ML_doubleMOC'; %%% Run to analyze
% run_name = 'ACC_AABW_ML_randWdia_randTau_white'; %%% Run to analyze
% run_name = 'ACC_AABW_ML_randWdia_randTau_white_Nlay2'; %%% Run to analyze
t1year = 86400*365;
rho0=1000;
tmin = 50*t1year; %%% Actual analysis period
tmax = 500*t1year;

%%% Load parameters   
local_home_dir = '/Volumes/Kilchoman/UCLA/Projects/AWSIM/runs';
% local_home_dir = '/data2/astewart/AWSIM/runs';
prod_dir = fullfile(local_home_dir,'AWSIM_WindAABW_products');
loadParams;
dirpath = fullfile(local_home_dir,run_name);
gtild = reshape(cumsum(gg),[Nlay 1 1]);

%%% Load pre-computed transport time series
load(fullfile(prod_dir,[run_name,'_MomBalance.mat']));

Nt = size(formStress,2);

%%% Interfacial form stress
ifs1 = formStress(1,:)/Lx/Ly*rho0;
ifs2 = formStress(2,:)/Lx/Ly*rho0;

%%% Topographic form stress
tfs = formStress(end,:)/Lx/Ly*rho0;


tauMean_5day = interp1([tauTimes tmax],[tauMax tauMax(1)],tt,'linear') / 2; %%% Divide max by 2 to get mean for sin^2 profile
wDiaInt_5day = interp1([wDiaTimes tmax],[wDiaInt wDiaInt(1)],tt,'linear');




figure(1);
plot(tt,tauMean_5day);
hold on;
plot(tt,formStress(end,:)/Lx/Ly*rho0);
hold off

figure(2);
scatter(tauMean_5day,tfs);
corrcoef(tauMean_5day,tfs)


figure(3);
plot(tt/t1year,tauMean_5day);
hold on;
plot(tt/t1year,formStress(1,:)/Lx/Ly*rho0);
plot(tt/t1year,formStress(2,:)/Lx/Ly*rho0);
plot(tt/t1year,formStress(3,:)/Lx/Ly*rho0);
hold off

figure(4);
scatter(tauMean_5day/2,formStress(1,:)/Lx/Ly*rho0);


modes = [0:1:Nt/2-1 -Nt/2:1:-1];
TT = tt(end) ./ modes;
omega = 2*pi./TT;

tau_fft = fft(tauMean_5day) / length(tauMean_5day);
tfs_fft = fft(tfs) / length(tfs);
ifs_fft = fft(ifs1) / length(ifs1);

C = real(tau_fft.*conj(tfs_fft));
Q = imag(tau_fft.*conj(tfs_fft));
Stau = tau_fft.*conj(tau_fft);
Stfs = tfs_fft.*conj(tfs_fft);
gamma = (C.^2 + Q.^2) ./ (Stau.*Stfs);
figure(5);
semilogx(TT/t1year,gamma,'.');

covar = 0.5*(tau_fft.*conj(tfs_fft)+conj(tau_fft).*tfs_fft);
covar = 2*covar(1:end/2); %%% Cut off half of wavenumbers, multiply by 2 to compensate
covar(1) = 0;
figure(6)
semilogx(TT(1:end/2)/t1year,cumsum(covar)/(std(tauMean_5day)*std(tfs)));

covar = 0.5*(tau_fft.*conj(ifs_fft)+conj(tau_fft).*ifs_fft) ./ (abs(tau_fft).*abs(ifs_fft));
covar(1) = 0;
figure(7);
semilogx(TT/t1year,covar,'.');

covar = 0.5*(tau_fft.*conj(ifs_fft)+conj(tau_fft).*ifs_fft);
covar = 2*covar(1:end/2); %%% Cut off half of wavenumbers, multiply by 2 to compensate
covar(1) = 0;
figure(8)
semilogx(TT(1:end/2)/t1year,cumsum(covar)/(std(tauMean_5day)*std(ifs1)));

smoothlen = [6:6:6000];
rr = 0*smoothlen;
for n=1:length(smoothlen)
  tau_smooth = smooth(tauMean_5day,smoothlen(n))';
  tfs_smooth = smooth(tfs,smoothlen(n))';
  rr(n) = corr(tau_smooth',tfs_smooth');
end
smoothlen = [0 smoothlen];
rr = [corr(tauMean_5day',tfs') rr];
figure(9);
semilogx(smoothlen/73,rr.^2);
xlabel('Smoothing window (years)');
ylabel('r^2');
title('Wind stress versus eddy topographic form stress');
set(gca,'FontSize',14);

smoothlen = [6:6:6000];
rr = 0*smoothlen;
for n=1:length(smoothlen)
  tau_smooth = smooth(tauMean_5day,smoothlen(n))';
  ifs1_smooth = smooth(ifs1,smoothlen(n))';
  rr(n) = corr(tau_smooth',ifs1_smooth');
end
smoothlen = [0 smoothlen];
rr = [corr(tauMean_5day',ifs1') rr];
figure(10);
semilogx(smoothlen/73,rr.^2);
xlabel('Smoothing window (years)');
ylabel('r^2');
title('Wind stress versus eddy interfacial form stress (upper interface)');
set(gca,'FontSize',14);

smoothlen = [6:6:6000];
rr = 0*smoothlen;
for n=1:length(smoothlen)
  tau_smooth = smooth(tauMean_5day,smoothlen(n))';
  ifs2_smooth = smooth(ifs2,smoothlen(n))';
  rr(n) = corr(tau_smooth',ifs2_smooth');
end
smoothlen = [0 smoothlen];
rr = [corr(tauMean_5day',ifs2') rr];
figure(11);
semilogx(smoothlen/73,rr.^2);
xlabel('Smoothing window (years)');
ylabel('r^2');
title('Wind stress versus eddy interfacial form stress (lower interface)');
set(gca,'FontSize',14);

tau_smooth = smooth(tauMean_5day,5*73)';
ifs2_smooth = smooth(ifs2,5*73)';
figure(12);
plot(tt,tau_smooth);
hold on;
plot(tt,ifs2_smooth);
hold off;
corr(tau_smooth',ifs2_smooth')^2
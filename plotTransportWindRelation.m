%%%
%%% plotTransportWindRelation.m
%%%
%%% Compares MOC variability with imposed wind forcing.
%%%

%%% Options
% run_name = 'ACC_AABW_ML_doubleMOC'; %%% Run to analyze
run_name = 'ACC_AABW_ML_randWdia_randTau_white'; %%% Run to analyze
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
load(fullfile(prod_dir,[run_name,'_InstTrans.mat']));

%%% Define AABW time series
Nt = size(Taabw,2);








%%% Power spectra of wind
tauMax_fft = fft(tauMax);
tauMax_fft = tauMax_fft(1:length(tauMax_fft)/2);
tauMax_fft(1) = 0;
Nfft = length(tauMax_fft);
TT = tmax./(0:1:Nfft-1);

figure(1);
loglog(TT/t1year,abs(tauMax_fft.^2));

figure(2);
semilogx(TT/t1year,cumsum(abs(tauMax_fft.^2),'reverse'));




%%% Power spectra of AABW transport
Taabw_fft = fft(Taabw);
Taabw_fft = Taabw_fft(1:length(Taabw_fft)/2);
Taabw_fft(1) = 0;
Nfft = length(Taabw_fft);
TT = tmax./(0:1:Nfft-1);

figure(3);
loglog(TT/t1year,abs(Taabw_fft.^2));

figure(4);
semilogx(TT/t1year,cumsum(abs(Taabw_fft.^2),'reverse'));



tauMax_5day = interp1([tauTimes tmax],[tauMax tauMax(1)],tt,'linear');
wDiaInt_5day = interp1([wDiaTimes tmax],[wDiaInt wDiaInt(1)],tt,'linear');

tauMax_smooth = smooth(tauMax_5day,730)';
tauMax_pert = tauMax_5day - tauMax_smooth;
Taabw_smooth = smooth(Taabw,730)';
Taabw_pert = Taabw - Taabw_smooth;




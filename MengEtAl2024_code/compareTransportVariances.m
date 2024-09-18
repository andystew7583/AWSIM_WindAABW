%%%
%%% compareTransportVariances.m
%%%
%%% Compares spectral content of diagnosed transports at high vs. low
%%% resolution.
%%%

%%% Paths
local_home_dir = '/Volumes/Kilchoman/UCLA/Projects/AWSIM/runs';
prod_dir = fullfile(local_home_dir,'AWSIM_WindAABW_products');

%%% Select chunk
m = 2;

%%% Load high-resolution data
run_name = 'ACC_AABW_ML_doubleMOC_hires';
load(fullfile(prod_dir,[run_name,'_InstTrans_chunk',num2str(2),'.mat']),'tt_chunk','Tacc_chunk','Taabw_chunk','Taaiw_chunk');
Tacc_hires = Tacc_chunk;
Taaiw_hires = Taaiw_chunk;
Taabw_hires = Taabw_chunk;

%%% Load low-resolution data
run_name = 'ACC_AABW_ML_doubleMOC';
load(fullfile(prod_dir,[run_name,'_InstTrans_chunk',num2str(2),'.mat']),'tt_chunk','Tacc_chunk','Taabw_chunk','Taaiw_chunk');

%%% Compute spectra
Taaiw_chunk_fft = fft(Taaiw_chunk)/length(tt_chunk); %%% Ensures that sum of squared Fourier components is equal to real-space variance
Taaiw_chunk_fft(1) = 0;
Taaiw_chunk_fft(end/2+1:end) = [];
Taaiw_hires_fft = fft(Taaiw_hires)/length(tt_chunk);
Taaiw_hires_fft(1) = 0;
Taaiw_hires_fft(end/2+1:end) = [];
Tacc_chunk_fft = fft(Tacc_chunk)/length(tt_chunk);
Tacc_chunk_fft(1) = 0;
Tacc_chunk_fft(end/2+1:end) = [];
Tacc_hires_fft = fft(Tacc_hires)/length(tt_chunk);
Tacc_hires_fft(1) = 0;
Tacc_hires_fft(end/2+1:end) = [];
Taabw_chunk_fft = fft(Taabw_chunk)/length(tt_chunk);
Taabw_chunk_fft(1) = 0;
Taabw_chunk_fft(end/2+1:end) = [];
Taabw_hires_fft = fft(Taabw_hires)/length(tt_chunk);
Taabw_hires_fft(1) = 0;
Taabw_hires_fft(end/2+1:end) = [];

nn = 0:length(tt_chunk)/2-1;
periods = (tt_chunk(end)-tt_chunk(1)) ./ nn / t1year;

figure(6);
plot(nn,flip(cumsum(abs(flip(Taaiw_hires_fft)).^2)./cumsum(abs(flip(Taaiw_chunk_fft)).^2)))
xlabel('Mode number');
ylabel('Ratio of cumulative variances');
title('AAIW transport');

figure(7);
plot(nn,flip(cumsum(abs(flip(Tacc_hires_fft)).^2)./cumsum(abs(flip(Tacc_chunk_fft)).^2)))
xlabel('Mode number');
ylabel('Ratio of cumulative variances');
title('ACC transport');

figure(8);
plot(nn,flip(cumsum(abs(flip(Taabw_hires_fft)).^2)./cumsum(abs(flip(Taabw_chunk_fft)).^2)))
xlabel('Mode number');
ylabel('Ratio of cumulative variances');
title('AAIW transport');
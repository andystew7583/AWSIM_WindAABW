%%%
%%% calcMomBalanceTimeSeries.m
%%%
%%% Computes time series of wind stress and form stresses.
%%%

%%% Load constant parameters
constants;
Ny = 128;
Nlay = 3;
tau_mean = [0.15];
% tau_pert = 0.075;
tau_pert = 0;
% tau_freq = t1year * 2^0;
tau_freq = 0;
% AABW_mean = 1.5;
AABW_mean = 0;
AABW_pert = 0;
AABW_freq = 0;
quad_drag = 2e-3;
lin_drag = 0e-4;  
topog_width = 150;
topog_height = 1000;
rough_topog = true;
rough_topog_height = 150;
double_wind = true;
double_ridge = false;
n_E_batch = 1:10;

run_name = constructRunName (false,Ny,Nlay, ...
                                  tau_mean,tau_pert,tau_freq, ...
                                  AABW_mean,AABW_pert,AABW_freq, ...
                                  quad_drag,lin_drag,...
                                  topog_width,topog_height,rough_topog,rough_topog_height,...
                                  double_ridge,double_wind,n_E_batch(end));
  
%%% Write to .mat file
rho0 = 1000;
f0 = -1e-4;
local_home_dir = '/Volumes/Stewart-RAID1-A/UCLA/Projects/AWSIM_WindAABW/runs_varywind';
prod_dir = fullfile(local_home_dir,'products');
load(fullfile(prod_dir,[run_name,'_MomBalance.mat']))

%%% Sample plot
figure(1);
plot(tt/t1year,mean(surfStress,3)/Lx/Ly);
hold on;
plot(tt/t1year,mean(formStress(1,:,:),3)/Lx/Ly);
plot(tt/t1year,mean(formStress(2,:,:),3)/Lx/Ly);
plot(tt/t1year,mean(formStress(3,:,:),3)/Lx/Ly);
hold off;
legend('Surface stress','IFS_u_p_p_e_r','IFS_l_o_w_e_r','TFS');
title('Momentum balance, channel-averaged');
xlabel('Time (years)')
ylabel('N/m^2');

f0 = 2*mean(Omega_z(:));
%%% Sample plot
figure(2);
plot(tt/t1year,mean(surfStress,3)/Ly/rho0/abs(f0)/1e6);
hold on;
plot(tt/t1year,mean(MOC(1,:,:),3)/1e6);
plot(tt/t1year,mean(MOC(2,:,:),3)/1e6);
plot(tt/t1year,mean(MOC(3,:,:),3)/1e6);
hold off;
legend('T_E_k_m_a_n','T_1','T_2','T_3');
title('Overturning, channel-averaged');
xlabel('Time (years)')
ylabel('Sv');

%%% Sample plot
figure(3);
plot(tt/t1year,mean(surfStress_mid,3)/Lx/(Ly/2));
hold on;
plot(tt/t1year,mean(formStress_mid(1,:,:),3)/Lx/(Ly/2));
plot(tt/t1year,mean(formStress_mid(2,:,:),3)/Lx/(Ly/2));
plot(tt/t1year,mean(formStress_mid(3,:,:),3)/Lx/(Ly/2));
hold off;
legend('Surface stress','IFS_u_p_p_e_r','IFS_l_o_w_e_r','TFS');
title('Momentum balance in channel center');
xlabel('Time (years)')
ylabel('N/m^2');

%%% Sample plot
figure(4);
plot(tt/t1year,mean(surfStress_mid,3)/(Ly/2)/rho0/abs(f0)/1e6);
hold on;
plot(tt/t1year,mean(MOC_mid(1,:,:),3)/1e6);
plot(tt/t1year,mean(MOC_mid(2,:,:),3)/1e6);
plot(tt/t1year,mean(MOC_mid(3,:,:),3)/1e6);
hold off;
legend('T_E_k_m_a_n','T_1','T_2','T_3');
title('Overturning, channel center');
xlabel('Time (years)')
ylabel('Sv');

%%% Sample plot
figure(5);
plot(tt/t1year,mean(Tacc,3)/1e6);
hold on;
plot(tt/t1year,mean(Tacc_bc,3)/1e6);
plot(tt/t1year,mean(Tacc_bt,3)/1e6);
hold off;
legend('Total transport','Baroclinic transport','Barotropic transport');
title('ACC transport');
xlabel('Time (years)')
ylabel('Sv');


%%% Sample plot
figure(6);
plot(tt/t1year,mean(formStress(3,:,:),3)/Lx/Ly);
hold on;
% plot(tt/t1year,mean(formStress_ridge(1,:,:)+formStress_bumps(1,:,:),3)/Lx/Ly);
plot(tt/t1year,mean(formStress_ridge(1,:,:),3)/Lx/Ly);
plot(tt/t1year,mean(formStress_bumps(1,:,:),3)/Lx/Ly);
hold off;
legend('TFS','TFS (ridge)','TFS (bumps)');
title('TFS decomposition, channel-averaged');
xlabel('Time (years)')
ylabel('N/m^2');
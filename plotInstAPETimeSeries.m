%%%
%%% plotInstAPETimeSeries.m
%%%
%%% Plots time series of domain integrated APE. 
%%% N.B. ONLY WORKS FOR 2 ISOPYCNAL LAYERS WITH A RIGID LID.
%%%

%%% Options
% run_name = 'ACC_AABW_Ny256_Nlay2_tauM0.5_tauP0_tauF0_wDiaM0_wDiaP0_wDiaF0_Cd5.000e-04_rb0.000e+00_diags'; %%% Run to analyze
run_name = 'ACC_AABW_Ny256_Nlay2_tauM0.5_tauP0_tauF0_wDiaM0_wDiaP0_wDiaF0_Cd0.000e+00_rb2.000e-04_diags'; %%% Run to analyze
local_home_dir = '/Volumes/Kilchoman/UCLA/Projects/AWSIM_WindAABW/runs';
dirpath = fullfile(local_home_dir,run_name);

%%% Load parameters   
loadParams;


%%% At each time iteration...
tt = zeros(1,Nframes);
APE = zeros(1,Nframes);
Tacc = zeros(1,Nframes);
cntr = 0;
for n=n0:1:n0+Nframes-1   

  %%% Current simulation time    
  t = startTime + (n-n0)*dt_s;
  cntr = cntr + 1;
  tt(cntr) = t;
  
  %%% Load zonal velocity and layer thicknesses
  uu = zeros(Nlay,Nx,Ny);
  hh = zeros(Nlay,Nx,Ny);
  for k=1:Nlay 

    %%% Load kth layer thickness
    data_file = fullfile(dirpath,[OUTN_H,num2str(k-1),'_n=',num2str(n),'.dat']);
    hh_tmp = readOutputFile(data_file,Nx,Ny);      
    if (isempty(hh_tmp))
      hh_tmp = NaN(Nx,Ny);
    end
    hh(k,:,:) = hh_tmp;
    
    %%% Load kth zonal velocity
    data_file = fullfile(dirpath,[OUTN_U,num2str(k-1),'_n=',num2str(n),'.dat']);
    uu_tmp = readOutputFile(data_file,Nx,Ny);      
    if (isempty(uu_tmp))
      uu_tmp = NaN(Nx,Ny);
    end
    uu(k,:,:) = uu_tmp;

  end
  
  %%% Calculate layer surface height
  layer = 2;
  eta = hhb;
  for k=Nlay:-1:layer
  
    %%% Add layer thickness to eta
    eta = eta + squeeze(hh(k,:,:));

  end
   
  %%% Compute APE
  APE(cntr) = sum(sum(0.5*gg(2)*(eta-mean(mean(eta))).^2*dx*dy));
  Tacc(cntr) = mean(sum(sum(uu.*hh*dy,1),3),2);
           
end

figure(10);
plot(tt/t1year,APE);
xlabel('Time');
ylabel('APE');

figure(11);
plot(tt/t1year,Tacc/1e6);
xlabel('Time');
ylabel('ACC transport (Sv)');

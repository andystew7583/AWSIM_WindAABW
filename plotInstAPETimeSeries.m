%%%
%%% plotInstAPETimeSeries.m
%%%
%%% Plots time series of domain integrated APE. 
%%% N.B. ONLY WORKS FOR 2 ISOPYCNAL LAYERS WITH A RIGID LID.
%%%

%%% Options
run_name = 'ACC_AABW_Ny256_Nlay2_tauM0.01_tauP0_tauF0_wDiaM0_wDiaP0_wDiaF0_Cd5.000e-04_rb0.000e+00_spinup'; %%% Run to analyze
local_home_dir = '/Volumes/Kilchoman/UCLA/Projects/AWSIM_WindAABW/runs';
dirpath = fullfile(local_home_dir,run_name);

%%% Load parameters   
loadParams;


%%% At each time iteration...
tt = zeros(1,Nframes);
APE = zeros(1,Nframes);
cntr = 0;
for n=n0:1:n0+Nframes-1   

  %%% Current simulation time    
  t = startTime + (n-n0)*dt_s;
  cntr = cntr + 1;
  tt(cntr) = t;
  
  %%% Calculate layer surface height
  layer = 2;
  eta = hhb;
  for k=Nlay:-1:layer

    %%% Load kth layer thickness
    data_file = fullfile(dirpath,[OUTN_H,num2str(k-1),'_n=',num2str(n),'.dat']);
    hh = readOutputFile(data_file,Nx,Ny);      
    if (isempty(hh))
      hh = NaN(Nx,Ny);
    end

    %%% Add layer thickness to eta
    eta = eta + hh;

  end
   
  %%% Compute APE
  APE(cntr) = sum(sum(0.5*gg(2)*(eta-mean(mean(eta))).^2*dx*dy));
           
end

figure(10);
plot(tt/t1year,APE);
xlabel('Time');
ylabel('APE');



%%% 
%%% create_batch.m
%%%
%%% Creates a batch of simulations with varying model input parameters.
%%%

%%% Load static definitions
constants;

%%% Directory to store runs
local_home_dir = '/Volumes/Stewart-RAID1-A/UCLA/Projects/AWSIM_WindAABW/runs_varywind';

%%% Remote cluster directory
uname = 'astewart';
% uname = 'andrewst';
cluster_addr = 'caolila.atmos.ucla.edu';
% cluster_addr = 'hoffman2.idre.ucla.edu';
cluster_home_dir = '/jbod/astewart/AWSIM_WindAABW/runs_varywind';
% cluster_home_dir = '/u/scratch/a/andrewst/AWSIM_WindAABW/runs';

%%% Spinup simulations are long and produce no diagnostics; diagnostic
%%% simulations output high-frequency diagnostics to resolve the forcing
%%% period
is_spinup = false;

%%% For non-spinup runs (i.e. production runs), start from the end of a
%%% corresponding spinup run with no time-variation in the forcing
start_from_steady_forcing = true;

%%% Set true to extend a previous run
extend_run = false;

%%% Grid resolution 
Ny = 128;
Nlay = 3;

%%% N.B. Batches run so far:
%%%
%%% N=128, steady wind with tau=0.15, no rough topog, 10 ensemble members
%%% spinup runs
%%% N=128, steady wind with tau=0.15, Taabw=1.5, no rough topog, 10 ensemble members
%%% spinup runs
%%% N=128, oscillating wind with tau=0.15, dtau=0.075, no rough topog, 10 ensemble members
%%% production runs

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameter selection %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Reference values
% tau_mean = 0.15;
% tau_pert = 0.05;
% tau_freq = 1*t1year;
% AABW_mean = 1.5;
% AABW_pert = 0;
% AABW_freq = 0;
% topog_width = 150;
% topog_height = 1000;

%%% Perturbation ensemble with rough topography
% Nensemble = 10;
% tau_mean = [0.15];
% tau_pert = 0.075;
% tau_freq = t1year .* 2.^[-3:1:4];
% AABW_mean = 0;
% AABW_pert = 0;
% AABW_freq = 0;
% quad_drag = 2e-3;
% lin_drag = 0e-4;  
% topog_width = 150;
% topog_height = 1000;
% rough_topog = true;

%%% Steady forcing ensemble with rough topography
% Nensemble = 10;
% tau_mean = [0.15];
% tau_pert = 0;
% tau_freq = 0;
% AABW_mean = 0;
% AABW_pert = 0;
% AABW_freq = 0;
% quad_drag = 2e-3;
% lin_drag = 0e-4;  
% topog_width = 150;
% topog_height = 1000;
% rough_topog = true;

%%% Perturbation ensemble with smooth topography
% Nensemble = 10;
% tau_mean = [0.15];
% tau_pert = 0.075;
% tau_freq = t1year .* 2.^[-3:1:4];
% AABW_mean = 0;
% AABW_pert = 0;
% AABW_freq = 0;
% quad_drag = 2e-3;
% lin_drag = 0e-4;  
% topog_width = 150;
% topog_height = 1000;
% rough_topog = false;

%%% Steady forcing ensemble with smooth topography
% Nensemble = 10;
% tau_mean = [0.15];
% tau_pert = 0;
% tau_freq = 0;
% AABW_mean = 0;
% AABW_pert = 0;
% AABW_freq = 0;
% quad_drag = 2e-3;
% lin_drag = 0e-4;  
% topog_width = 150;
% topog_height = 1000;
% rough_topog = false;

%%% Perturbation ensemble with smooth topography plus AABW
% Nensemble = 10;
% tau_mean = [0.15];
% tau_pert = 0.075;
% tau_freq = t1year .* 2.^[-3:1:4];
% AABW_mean = 1.5;
% AABW_pert = 0;
% AABW_freq = 0;
% quad_drag = 2e-3;
% lin_drag = 0e-4;  
% topog_width = 150;
% topog_height = 1000;
% rough_topog = false;

%%% Steady forcing ensemble with smooth topography plus AABW
% Nensemble = 10;
% tau_mean = [0.15];
% tau_pert = 0;
% tau_freq = 0;
% AABW_mean = 1.5;
% AABW_pert = 0;
% AABW_freq = 0;
% quad_drag = 2e-3;
% lin_drag = 0e-4;  
% topog_width = 150;
% topog_height = 1000;
% rough_topog = false;

%%% Perturbation ensemble with rough topography plus AABW
Nensemble = 10;
tau_mean = [0.15];
tau_pert = 0.075;
tau_freq = t1year .* 2.^[-3:1:4];
AABW_mean = 1.5;
AABW_pert = 0;
AABW_freq = 0;
quad_drag = 2e-3;
lin_drag = 0e-4;  
topog_width = 150;
topog_height = 1000;
rough_topog = true;

%%% Steady forcing ensemble with rough topography plus AABW
% Nensemble = 10;
% tau_mean = [0.15];
% tau_pert = 0;
% tau_freq = 0;
% AABW_mean = 1.5;
% AABW_pert = 0;
% AABW_freq = 0;
% quad_drag = 2e-3;
% lin_drag = 0e-4;  
% topog_width = 150;
% topog_height = 1000;
% rough_topog = true;


%%% Script files
run_batch_fname = 'run_batch.sh';
run_batch_file = fopen(fullfile(local_home_dir,run_batch_fname),'w');
upload_batch_fname = 'upload_batch.sh';
upload_batch_file = fopen(fullfile(local_home_dir,upload_batch_fname),'w');

%%% Start upload file
fprintf(upload_batch_file,'rsync -av --update %s ',run_batch_fname);

%%% Create runs
for n_tm=1:length(tau_mean)
  for n_tp = 1:length(tau_pert)
    for n_tf = 1:length(tau_freq)
      for n_am = 1:length(AABW_mean)
        for n_ap = 1:length(AABW_pert)
          for n_af = 1:length(AABW_freq)
            for n_Cd = 1:length(quad_drag)
              for n_rb = 1:length(lin_drag)
                for n_Wb = 1:length(topog_width)
                  for n_Hb = 1:length(topog_height)
                    for n_E = 1:Nensemble

                      %%% Generate simulation name
                      run_name = constructRunName (is_spinup,Ny,Nlay, ...
                                  tau_mean(n_tm),tau_pert(n_tp),tau_freq(n_tf), ...
                                  AABW_mean(n_am),AABW_pert(n_ap),AABW_freq(n_af), ...
                                  quad_drag(n_Cd),lin_drag(n_rb),topog_width(n_Wb),topog_height(n_Hb),rough_topog,n_E);

                      %%% Identify previous simulation from which to copy the restart file
                      if (is_spinup)

                        %%% Start low-res run from a previous long integration in the
                        %%% same geometry
                        if (Ny == 128)

                          if (extend_run)
                            dir_pickup = local_home_dir;
                            run_name_pickup = run_name;
                            pickup_iter = findLastOutput(dir_pickup,run_name_pickup);
                            restart_idx = pickup_iter;
                            end_time = 500*t1year;
                          else
                            dir_pickup = '/Volumes/Kilchoman/UCLA/Projects/AWSIM/runs';
                            run_name_pickup = 'ACC_AABW_ML_doubleMOC_hires';
                            pickup_iter = 7300;  
                            restart_idx = -1;
                            if (tau_mean(n_tm)<0.05)
                              end_time = 400*t1year;
                            else
                              end_time = 200*t1year;
                            end
                          end

                        %%% Start hi-res run from the end of the low-res run
                        else %%% N=256               

                          if (extend_run)
                            dir_pickup = local_home_dir;
                            run_name_pickup = run_name;
                            pickup_iter = findLastOutput(dir_pickup,run_name_pickup);
                            restart_idx = pickup_iter;
                            end_time = 100*t1year;
                          else
                            dir_pickup = local_home_dir;
                            run_name_pickup = constructRunName (true,Ny/2,Nlay, ...
                                    tau_mean(n_tm),tau_pert(n_tp),tau_freq(n_tf), ...
                                    AABW_mean(n_am),AABW_pert(n_ap),AABW_freq(n_af), ...
                                    quad_drag(n_Cd),lin_drag(n_rb),topog_width(n_Wb),topog_height(n_Hb),rough_topog,n_E);
                            pickup_iter = findLastOutput(dir_pickup,run_name_pickup);        
                            restart_idx = 0;
                            end_time = 100*t1year;
                          end

                        end

                      %%% Diagnostic runs start from the end of the corresponding
                      %%% spinup runs
                      else
                        
                        dir_pickup = local_home_dir;
                        
                        if (start_from_steady_forcing)
                          run_name_pickup = constructRunName (true,Ny,Nlay, ...
                                  tau_mean(n_tm),0,0, ...
                                  AABW_mean(n_am),0,0, ...
                                  quad_drag(n_Cd),lin_drag(n_rb),topog_width(n_Wb),topog_height(n_Hb),rough_topog,n_E);
                        else                        
                          run_name_pickup = constructRunName (true,Ny,Nlay, ...
                                    tau_mean(n_tm),tau_pert(n_tp),tau_freq(n_tf), ...
                                    AABW_mean(n_am),AABW_pert(n_ap),AABW_freq(n_af), ...
                                    quad_drag(n_Cd),lin_drag(n_rb),topog_width(n_Wb),topog_height(n_Hb),rough_topog,n_E);
                        end
                        pickup_iter = findLastOutput(dir_pickup,run_name_pickup);    
                        restart_idx = 0;
                        end_time = 30*t1year; %%% Default - only used for steady simulations

                      end

                      %%% Create simulation directory and input files
                      rng(n_E)
                      setparams (local_home_dir,run_name,is_spinup,Ny,Nlay, ...
                                  tau_mean(n_tm),tau_pert(n_tp),tau_freq(n_tf), ...
                                  AABW_mean(n_am),AABW_pert(n_ap),AABW_freq(n_af), ...
                                  quad_drag(n_Cd),lin_drag(n_rb),topog_width(n_Wb),topog_height(n_Hb),rough_topog, ...
                                  restart_idx, end_time);

                      %%% Copy/regrid pickup files to the simulation directory
                      if (~extend_run && (restart_idx >= 0))
                        regridOutput(dir_pickup,run_name_pickup,pickup_iter,2*Ny,Ny,fullfile(local_home_dir,run_name),0);  
                      end

                      %%% Add lines to run_batch file to execute this simulation  
                      fprintf(run_batch_file,'cd %s\n',run_name);
                      fprintf(run_batch_file,'sh Make_fftw.sh\n');
                      fprintf(run_batch_file,'sh Run.sh\n');
                      fprintf(run_batch_file,'cd ..\n');

                      %%% Add text to upload_batch file to upload this
                      %%% simulation
                      fprintf(upload_batch_file,'%s ',run_name);

                    end
                  end
                end
              end
            end            
          end
        end
      end
    end
  end
end

%%% Finish upload command
fprintf(upload_batch_file,'%s@%s:%s\n',uname,cluster_addr,cluster_home_dir);

%%% Close open file handles
fclose(run_batch_file);
fclose(upload_batch_file);
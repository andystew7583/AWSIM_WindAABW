%%% 
%%% create_batch.m
%%%
%%% Creates a batch of simulations with varying model input parameters.
%%%

%%% Load static definitions
constants;

%%% Directory to store runs
local_home_dir = '/Volumes/Kilchoman/UCLA/Projects/AWSIM_WindAABW/runs';

%%% Remote cluster directory
uname = 'andrewst';
cluster_addr = 'hoffman2.idre.ucla.edu';
cluster_home_dir = '/u/scratch/a/andrewst/AWSIM_WindAABW/runs';

%%% Spinup simulations are long and produce no diagnostic, diagnostic
%%% simulations output high-frequency diagnostics to resolve the forcing
%%% period
is_spinup = false;

%%% Set true to extend a previous run
extend_run = false;

%%% Grid resolution 
Ny = 256;
Nlay = 2;

%%% N.B. Batches run so far:
%%% N=128
%%% - tau_mean = [0.01 0.017 0.03 0.05 0.1 0.17 0.3], AABW_mean = 0,
%%% quad_drag = 2e-3, lin_drag = 0e-4
%%% spinup run, diagnostics run
%%%
%%% - tau_mean = [0.01 0.017 0.03 0.05 0.1 0.17 0.3], AABW_mean = 0,
%%% quad_drag = 0e-3, lin_drag = 2e-4
%%% spinup run, diagnostics run
%%%
%%% - tau_mean = [0.01 0.017 0.03 0.05 0.1 0.17 0.3], AABW_mean = [-1.5 -.75
%%% .75 1.5], quad_drag = 2e-3, lin_drag = 0e-4
%%% ???
%%%
%%% - tau_mean = [0.01 0.017 0.03 0.05 0.1 0.17 0.3], AABW_mean = [-1.5 
%%% 1.5], quad_drag = 0e-3, lin_drag = 2e-4
%%% ???
%%%
%%% - tau_mean = [0.01 0.017 0.03 0.05 0.1 0.17 0.3], AABW_mean = 0,
%%% quad_drag = 2e-3, lin_drag = 0e-4, topog_width
%%% spinup and diagnostics run
%%%
%%% - tau_mean = [0.01 0.017 0.03 0.05 0.1 0.17 0.3], AABW_mean = 0,
%%% quad_drag = 2e-3, lin_drag = 0e-4, topog_width = [40 80 300 600];
%%% ???
%%%
%%% - tau_mean = [0.01 0.017 0.03 0.05 0.1 0.17 0.3], AABW_mean = 0,
%%% quad_drag = 2e-3, lin_drag = 0e-4, topog_width = [40 80 150 300 600], topog_height = [500 1500];
%%% ???
%%%
%%% - tau_mean = 0.1, AABW_mean = 0, quad_drag = [.5e-3 1e-3 1.5e-3 2.5e-3
%%% 3e-3 3.5e-3 4e-3], lin_drag = 0, topog_width = 150, topog_height = 1000
%%% spinup run, diagnostics run
%%%
%%% - tau_mean = 0.1, AABW_mean = 0, quad_drag = 0 
%%% lin_drag = [1e-4 3e-4 4e-4 5e-4 6e-4 7e-4 8e-4 9e-4 10e-4], topog_width = 150, topog_height = 1000
%%% spinup run, diagnostics run
%%%
%%% - tau_mean =  [0.013 0.022 0.039 0.07 0.13 0.22], AABW_mean = 0,
%%% quad_drag = 2e-3, lin_drag = 0e-4, topog_width = 150, topog_height =
%%% 1000
%%% spinup run, diagnostics run
%%%
%%% - tau_mean =  [0.013 0.022 0.039 0.07 0.13 0.22], AABW_mean = 0,
%%% quad_drag = 0e-3, lin_drag = 2e-4, topog_width = 150, topog_height =
%%% 1000
%%% spinup run, diagnostics run
%%%
%%% - tau_mean =  [0.01 0.013 0.017 0.022 0.03 0.039 0.05 0.07 0.13 0.17 0.22 0.3], AABW_mean = 0,
%%% quad_drag = [.5e-3 1e-3 1.5e-3 2.5e-3 3e-3 3.5e-3 4e-3], lin_drag = 0e-4, topog_width = 150, topog_height =
%%% 1000
%%% spinup run, diagnostics run
%%%
%%% - tau_mean =  [0.01 0.013 0.017 0.022 0.03 0.039 0.05 0.07 0.13 0.17 0.22 0.3], AABW_mean = 0,
%%% quad_drag = 0, lin_drag = [1e-4 3e-4 4e-4 5e-4 6e-4 7e-4 8e-4 9e-4 10e-4], topog_width = 150, topog_height =
%%% 1000
%%% spinup run, diagnostics run
%%%
%%% - tau_mean =  [0.39 0.5], AABW_mean = 0,
%%% quad_drag = [.5e-3 1e-3 1.5e-3 2e-3 2.5e-3 3e-3 3.5e-3 4e-3], lin_drag = 0e-4, topog_width = 150, topog_height =
%%% 1000
%%% spinup running on hoffman
%%%
%%% - tau_mean =  [0.39 0.5], AABW_mean = 0,
%%% quad_drag = 0, lin_drag = [1e-4 2e-4 3e-4 4e-4 5e-4 6e-4 7e-4 8e-4 9e-4 10e-4], topog_width = 150, topog_height =
%%% 1000
%%% spinup running on hoffman
%%%
%%% N = 256
%%% - tau_mean = [0.01 0.017 0.03 0.05 0.1 0.17 0.3], AABW_mean = 0,
%%% quad_drag = 2e-3, lin_drag = 0e-4, topog_width = 150, topog_height = 1000
%%% spinup run, diagnostics run
%%%
%%% - tau_mean =  [0.013 0.022 0.039 0.07 0.13 0.22], AABW_mean = 0,
%%% quad_drag = 2e-3, lin_drag = 0e-4, topog_width = 150, topog_height =
%%% 1000
%%% spinup run, diagnostics running on hoffman
%%%
%%% - tau_mean =  [0.01 0.013 0.017 0.022 0.03 0.039 0.05 0.07 0.1 0.13 0.17 0.22 0.3], AABW_mean = 0,
%%% quad_drag = [.5e-3 1e-3 1.5e-3 2.5e-3 3e-3 3.5e-3 4e-3], lin_drag = 0e-4, topog_width = 150, topog_height =
%%% 1000
%%% spinup running on hoffman
%%%
%%% - tau_mean =  [0.01 0.013 0.017 0.022 0.03 0.039 0.05 0.07 0.1 0.13 0.17 0.22 0.3], AABW_mean = 0,
%%% lin_drag = [1e-4 2e-4 3e-4 4e-4 5e-4 6e-4 7e-4 8e-4 9e-4 10e-4], lin_drag = 0e-4, topog_width = 150, topog_height =
%%% 1000
%%% spinup running on hoffman

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

%%% Perturbation values
% tau_mean = [0.03 0.05 0.1 0.17 0.3];
% tau_mean = [0.01 0.017 0.03 0.05 0.1 0.17 0.3];
% tau_mean = [0.01 0.017 0.03 0.05 0.1];
% tau_mean = [0.01 0.017];
% tau_mean = 0.1;
% tau_mean = [0.013 0.022 0.039 0.07 0.13 0.22];
% tau_mean = [0.01 0.013 0.017 0.022 0.03 0.039 0.05 0.07 0.13 0.17 0.22 0.3];
tau_mean = [0.01 0.013 0.017 0.022 0.03 0.039 0.05 0.07 0.1 0.13 0.17 0.22 0.3];
% tau_mean = [0.39 0.5];
tau_pert = 0;
tau_freq = 0;
% AABW_mean = [-1.5 1.5];
% AABW_mean = [-1.5 0 1.5];
% AABW_mean = [-.75 .75];
% AABW_mean = [-1.5 -.75 0 .75 1.5];
AABW_mean = 0;
AABW_pert = 0;
AABW_freq = 0;
% quad_drag = 2e-3;
% quad_drag = 0e-3;
quad_drag = [.5e-3 1e-3 1.5e-3 2.5e-3 3e-3 3.5e-3 4e-3];
% quad_drag = [.5e-3 1e-3 1.5e-3 2e-3 2.5e-3 3e-3 3.5e-3 4e-3];
lin_drag = 0e-4;  
% lin_drag = [1e-4 3e-4# 4e-4 5e-4 6e-4 7e-4 8e-4 9e-4 10e-4];
% lin_drag = [1e-4 2e-4 3e-4 4e-4 5e-4 6e-4 7e-4 8e-4 9e-4 10e-4];
% quad_drag = 0e-3;
% lin_drag = 2e-4;
topog_width = 150;
% topog_width = [40 80 150 300 600];
topog_height = 1000;
% topog_height = [500 1500];

%%% Wind perturbation batch
% tau_mean = [0 0.05 0.1 0.15 0.2];
% tau_pert = [0.05 0.1];
% tau_freq = t1day*[1 3 10 30 100 t1year 3*t1year 10*t1year 30*t1year];
% AABW_mean = 1.5;
% AABW_pert = 0;
% AABW_freq = 0;

%%% AABW perturbation batch
% tau_mean = 0.15;
% tau_pert = 0;
% tau_freq = 0;
% AABW_mean = [0 0.5 1 1.5 2];
% AABW_pert = [0.5 1];
% AABW_freq = t1day*[1 3 10 30 100 t1year 3*t1year 10*t1year 30*t1year];



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

                    %%% Generate simulation name
                    run_name = constructRunName (is_spinup,Ny,Nlay, ...
                                tau_mean(n_tm),tau_pert(n_tp),tau_freq(n_tf), ...
                                AABW_mean(n_am),AABW_pert(n_ap),AABW_freq(n_af), ...
                                quad_drag(n_Cd),lin_drag(n_rb),topog_width(n_Wb),topog_height(n_Hb));

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

                        dir_pickup = local_home_dir;
                        run_name_pickup = constructRunName (true,Ny/2,Nlay, ...
                                tau_mean(n_tm),tau_pert(n_tp),tau_freq(n_tf), ...
                                AABW_mean(n_am),AABW_pert(n_ap),AABW_freq(n_af), ...
                                quad_drag(n_Cd),lin_drag(n_rb),topog_width(n_Wb),topog_height(n_Hb));
                        pickup_iter = findLastOutput(dir_pickup,run_name_pickup);        
                        restart_idx = 0;
                        end_time = 100*t1year;

                      end

                    %%% Diagnostic runs start from the end of the corresponding
                    %%% spinup runs
                    else

                      dir_pickup = local_home_dir;
                      run_name_pickup = constructRunName (true,Ny,Nlay, ...
                                tau_mean(n_tm),tau_pert(n_tp),tau_freq(n_tf), ...
                                AABW_mean(n_am),AABW_pert(n_ap),AABW_freq(n_af), ...
                                quad_drag(n_Cd),lin_drag(n_rb),topog_width(n_Wb),topog_height(n_Hb));
                      pickup_iter = findLastOutput(dir_pickup,run_name_pickup);    
                      restart_idx = 0;
                      end_time = 30*t1year;

                    end

                    %%% Create simulation directory and input files
                    setparams (local_home_dir,run_name,is_spinup,Ny,Nlay, ...
                                tau_mean(n_tm),tau_pert(n_tp),tau_freq(n_tf), ...
                                AABW_mean(n_am),AABW_pert(n_ap),AABW_freq(n_af), ...
                                quad_drag(n_Cd),lin_drag(n_rb),topog_width(n_Wb),topog_height(n_Hb), ...
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

%%% Finish upload command
fprintf(upload_batch_file,'%s@%s:%s\n',uname,cluster_addr,cluster_home_dir);

%%% Close open file handles
fclose(run_batch_file);
fclose(upload_batch_file);
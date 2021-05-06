%%% 
%%% create_batch.m
%%%
%%% Creates a batch of simulations with varying model input parameters.
%%%

%%% Load static definitions
constants;

%%% Directory to store runs
local_home_dir = './runs';

%%% Spinup simulations are long and produce no diagnostic, diagnostic
%%% simulations output high-frequency diagnostics to resolve the forcing
%%% period
is_spinup = true;

%%% Grid resolution 
Ny = 128;
Nlay = 2;



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

%%% Reference values
tau_mean = [0.01 0.017 0.03 0.05 0.1 0.17 0.3];
tau_pert = 0;
tau_freq = 0;
AABW_mean = [-1.5 0 1.5];
AABW_pert = 0;
AABW_freq = 0;
% quad_drag = 2e-3;
% lin_drag = 0e-4;
quad_drag = 0e-3;
lin_drag = 2e-4;

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

%%% Create runs
for n_tm=1:length(tau_mean)
  for n_tp = 1:length(tau_pert)
    for n_tf = 1:length(tau_freq)
      for n_am = 1:length(AABW_mean)
        for n_ap = 1:length(AABW_pert)
          for n_af = 1:length(AABW_freq)
            for n_Cd = 1:length(quad_drag)
              for n_rb = 1:length(lin_drag)

                %%% Generate simulation name
                run_name = constructRunName (is_spinup,Ny,Nlay, ...
                            tau_mean(n_tm),tau_pert(n_tp),tau_freq(n_tf), ...
                            AABW_mean(n_am),AABW_pert(n_ap),AABW_freq(n_af), ...
                            quad_drag(n_Cd),lin_drag(n_rb));

                %%% Create simulation directory and input files
                setparams (local_home_dir,run_name,is_spinup,Ny,Nlay, ...
                            tau_mean(n_tm),tau_pert(n_tp),tau_freq(n_tf), ...
                            AABW_mean(n_am),AABW_pert(n_ap),AABW_freq(n_af), ...
                            quad_drag(n_Cd),lin_drag(n_rb));

                %%% Identify previous simulation from which to copy the restart file
                if (is_spinup)

                  %%% Start low-res run from a previous long integration in the
                  %%% same geometry
                  if (Ny == 128)

                    dir_pickup = '../AWSIM/runs';
                    run_name_pickup = 'ACC_AABW_ML_doubleMOC_hires';
                    pickup_iter = 7300;                

                  %%% Start hi-res run from the end of the low-res run
                  else %%% N=256               

                    dir_pickup = local_home_dir;
                    run_name_pickup = constructRunName (true,Ny/2,Nlay, ...
                            tau_mean(n_tm),tau_pert(n_tp),tau_freq(n_tf), ...
                            AABW_mean(n_am),AABW_pert(n_ap),AABW_freq(n_af), ...
                            quad_drag(n_Cd),lin_drag(n_rb));
                    pickup_iter = 100;               

                  end

                %%% Diagnostic runs start from the end of the corresponding
                %%% spinup runs
                else

                  dir_pickup = local_home_dir;
                  run_name_pickup = constructRunName (false,Ny,Nlay, ...
                            tau_mean(n_tm),tau_pert(n_tp),tau_freq(n_tf), ...
                            AABW_mean(n_am),AABW_pert(n_ap),AABW_freq(n_af), ...
                            quad_drag(n_Cd),lin_drag(n_rb));
                  pickup_iter = 100;              

                end

                %%% Copy/regrid pickup files to the simulation directory
                regridOutput(dir_pickup,run_name_pickup,pickup_iter,2*Ny,Ny,fullfile(local_home_dir,run_name),0);  

                %%% Add lines to run_batch file to execute this simulation  
                fprintf(run_batch_file,'cd %s\n',run_name);
                fprintf(run_batch_file,'sh Make_fftw.sh\n');
                fprintf(run_batch_file,'sh Run.sh\n');
                fprintf(run_batch_file,'cd ..\n');

              end
            end            
          end
        end
      end
    end
  end
end

%%% Close open file handles
fclose(run_batch_file);
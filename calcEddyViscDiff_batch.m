%%%
%%% calcEddyViscDiff_batch.m
%%%
%%% Calculate eddy viscosity and diffusivity maps over a batch of
%%% simulations.
%%%

%%% Load static definitions
constants;

%%% Directory to store runs
local_home_dir = '/Volumes/Kilchoman/UCLA/Projects/AWSIM_WindAABW/runs';

%%% Spinup simulations are long and produce no diagnostic, diagnostic
%%% simulations output high-frequency diagnostics to resolve the forcing
%%% period
is_spinup = false;

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
% topog_width = 150;
% topog_height = 1000;

%%% Perturbation values
% tau_mean = [0.01 0.013 0.017 0.022 0.03 0.039 0.05 0.07 0.13 0.17 0.22 0.3];
% tau_mean = [0.17 0.22 0.3];
tau_mean = [0.39 0.5];
% tau_mean = 0.1;
tau_pert = 0;
tau_freq = 0;
AABW_mean = 0;
AABW_pert = 0;
AABW_freq = 0;
% quad_drag = 2e-3;
% quad_drag = 0e-3;
% quad_drag = [.5e-3 1e-3 1.5e-3 2.5e-3 3e-3 3.5e-3 4e-3];
quad_drag = [.5e-3 1e-3 1.5e-3 2e-3 2.5e-3 3e-3 3.5e-3 4e-3];
lin_drag = 0e-4;
% lin_drag = [1e-4 2e-4 3e-4 4e-4 5e-4 6e-4 7e-4 8e-4 9e-4 10e-4];
topog_width = 150;
topog_height = 1000;

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
                              
                    %%% Compute eddy diffusivity and viscosity and save to
                    %%% file
                    [kap_map,nu_map,r_kap_map,r_nu_map,p_kap_map,p_nu_map] = calcEddyViscDiff(local_home_dir,run_name); 
                    [kap_bulk,nu_bulk,r_kap_bulk,r_nu_bulk,EKE_zavg] = calcBulkEddyViscDiff(local_home_dir,run_name);
                    save(fullfile('products',['kap_nu_',run_name,'.mat']), ...
                      'kap_map','nu_map','r_kap_map','r_nu_map','p_kap_map','p_nu_map', ...
                      'kap_bulk','nu_bulk','r_kap_bulk','r_nu_bulk','EKE_zavg');

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

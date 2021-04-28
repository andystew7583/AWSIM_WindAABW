%%%
%%% constructRunName.m
%%%
%%% Convenience function to construct simulation names consistently.
%%%
function run_name = constructRunName (is_spinup,N,...
                        tau_mean,tau_pert,tau_freq, ...
                        AABW_mean,AABW_pert,AABW_freq)

  %%% Load definitions
  constants;
                      
  %%% Base simulation name
  run_name = 'ACC_AABW';
  
  %%% Resolution
  run_name = [run_name,'_N',num2str(N)];
  
  %%% Parameter values
  run_name = [run_name,'_tauM',num2str(tau_mean)];
  run_name = [run_name,'_tauP',num2str(tau_pert)];
  run_name = [run_name,'_tauF',num2str(tau_freq/t1day)];
  run_name = [run_name,'_wDiaM',num2str(AABW_mean)];
  run_name = [run_name,'_wDiaP',num2str(AABW_pert)];
  run_name = [run_name,'_wDiaF',num2str(AABW_freq/t1day)];
  
  %%% Distinguish spinup runs from runs with online averaging
  if (is_spinup)  
    run_name = [run_name,'_spinup']; 
  else
    run_name = [run_name,'_diags']; 
  end

end
%%%
%%% constructRunName.m
%%%
%%% Convenience function to construct simulation names consistently.
%%%
function run_name = constructRunName ( ...
    is_spinup,grid_size,num_layers, ...
    tau_mean,tau_pert,tau_freq, ...
    AABW_mean,AABW_pert,AABW_freq, ...
    quad_drag,lin_drag,topog_width,topog_height,ensemble_id)

  %%% Load definitions
  constants;
                      
  %%% Base simulation name
  run_name = 'ACC_AABW';
  
  %%% Resolution
  run_name = [run_name,'_Ny',num2str(grid_size)];
  run_name = [run_name,'_Nlay',num2str(num_layers)];
  
  %%% Parameter values
  run_name = [run_name,'_tauM',num2str(tau_mean)];
  run_name = [run_name,'_tauP',num2str(tau_pert)];
  run_name = [run_name,'_tauF',num2str(tau_freq/t1day)];
  run_name = [run_name,'_wDiaM',num2str(AABW_mean)];
  run_name = [run_name,'_wDiaP',num2str(AABW_pert)];
  run_name = [run_name,'_wDiaF',num2str(AABW_freq/t1day)];
  run_name = [run_name,'_Cd',num2str(quad_drag,'%.3e')];
  run_name = [run_name,'_rb',num2str(lin_drag,'%.3e')];  
  if (topog_width ~= 150)
    run_name = [run_name,'_Wb',num2str(topog_width)];
  end  
  if (topog_height ~= 1000)
    run_name = [run_name,'_Hb',num2str(topog_height)];
  end  
  run_name = [run_name,'_E',num2str(ensemble_id)];
  
  %%% Distinguish spinup runs from runs with online averaging
  if (is_spinup)  
    run_name = [run_name,'_spinup']; 
  else
    run_name = [run_name,'_diags']; 
  end

end
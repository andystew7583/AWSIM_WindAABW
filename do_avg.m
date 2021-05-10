%%%
%%% do_avg
%%%
%%% Convenience function to do time averaging
%%%
function var_avg = do_avg (dirpath,var_name,Nx,Ny,Nlay,n0,N,dt,tmin,tmax,startTime)
    
  constants;

  %%% Loop through output files to calculate average
  var_avg = zeros(Nx,Ny,Nlay);
  n_avg = 0;
  for n=n0+1:1:n0+N

    %%% Simulation time at the end of the averaging period
    t = startTime + (n-n0)*dt;

    %%% Check that the output falls within the desired time frame
    if ((tmin >= 0 && t < tmin) || (tmax >= 0 && t > tmax))
      continue;
    end

    %%% Increment averaging counter
    n_avg = n_avg + 1;

    %%% Add diagnostics to average
    for k=1:Nlay      
      if (strcmp(var_name,OUTN_PI_AVG))
        data_file = fullfile(dirpath,[var_name,'_n=',num2str(n),'.dat']);      
      else
        data_file = fullfile(dirpath,[var_name,num2str(k-1),'_n=',num2str(n),'.dat']);      
      end
      var_avg(:,:,k) = var_avg(:,:,k) + reshape(readOutputFile(data_file,Nx,Ny),[Nx Ny 1]);      
    end        


  end

  %%% Divide by number of files read to calculate time-mean
  var_avg = var_avg/n_avg;

end
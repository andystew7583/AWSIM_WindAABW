%%%
%%% timeavg_QG.m
%%%
%%% Reads in data from the output of 'AWSIM' and computes the time
%%% averages of various QG quantities between tmin and tmax.
%%%
function timeavg_QG (local_home_dir,run_name,tmin,tmax) 

  %%% Load parameters  
  loadParams;    
  
  %%% To store instantaneous output  
  hh = zeros(Nx,Ny,Nlay);
  eta = zeros(Nx,Ny,Nlay+1);  
  MM = zeros(Nx,Ny,Nlay);  
  
  %%% To store averages
  psig_tavg = zeros(Nx,Ny,Nlay);
  ug_tavg = zeros(Nx,Ny,Nlay);
  vg_tavg = zeros(Nx,Ny,Nlay);
  zetag_tavg = zeros(Nx,Ny,Nlay);
  uzg_tavg = zeros(Nx,Ny,Nlay);
  vzg_tavg = zeros(Nx,Ny,Nlay);
  
  %%% Calculate gsum, the sum of reduced gravities down to the layer in
  %%% question. This definition gives gsum_k*rho0 = g*rho_k, so this is
  %%% effectively a measure of the density in each layer.
  gsum = cumsum(gg);
  f0 = mean(mean(2*Omega_z));
  
  %%% At each time iteration...
  navg = 0;
  t_total = 0;
  for n=n0:1:n0+Nframes-1 
    n
    %%% Current simulaxtion time
    t = startTime + (n-n0)*dt_s;
    
    if (t < tmin || t > tmax)
      continue;
    end
    
    %%% Load data for this layer
    for k = 1:Nlay                     
      data_file = fullfile(dirpath,createOutputFilename(OUTN_H,n,k));
      hh(:,:,k) = readOutputFile(data_file,Nx,Ny);                         
    end    
    
    %%% Calculate surface pressure
    if (useRL)            
      data_file = fullfile(dirpath,createOutputFilename(OUTN_PI,n,-1));
      pi = readOutputFile(data_file,Nx,Ny);                     
    else
      pi = zeros(Nx,Ny);
    end
    
    %%% Calculate Montgomery potential       
    MM(:,:,1) = pi;    
    for k=2:Nlay
      MM(:,:,k) = MM(:,:,k-1) + gg(k)*eta(:,:,k);          
    end     
    for k=1:Nlay
      MM(:,:,k) = MM(:,:,k) - gsum(k) .* h0.^4 ./ hh(:,:,k).^3 ./ 3;
    end
    
    %%% Calculate layer surface heights and layer mid-depths
    eta(:,:,Nlay+1) = hhb;      
    for k=Nlay:-1:1               
      eta(:,:,k) = eta(:,:,k+1) + hh(:,:,k);
    end
       
    %%% Geostrophic velocity and vorticity estimates
    psig = MM / f0;
    vg = (psig([2:Nx 1],:,:) - psig([Nx 1:Nx-1],:,:)) / (2*dx);
    ug = zeros(Nx,Ny,Nlay);
    ug(:,2:Ny-1,:) = - (psig(:,3:Ny,:) - psig(:,1:Ny-2,:)) / (2*dy);
    ug(:,1,:) = - (psig(:,2,:) - psig(:,1,:)) / dy;
    ug(:,Ny,:) = - (psig(:,Ny,:) - psig(:,Ny-1,:)) / dy;
    zetag = zeros(Nx,Ny,Nlay);
    zetag(:,2:Ny-1,:) = (psig([2:Nx 1],2:Ny-1,:) - 2*psig(1:Nx,2:Ny-1,:) + psig([Nx 1:Nx-1],2:Ny-1,:)) / dx^2 ...
                      + (psig(:,3:Ny,:) - 2*psig(:,2:Ny-1,:) + psig(:,1:Ny-2,:)) / dy^2;
    zetag(:,1,:) = (psig([2:Nx 1],1,:) - 2*psig(1:Nx,1,:) + psig([Nx 1:Nx-1],1,:)) / dx^2;
    zetag(:,Ny,:) = (psig([2:Nx 1],Ny,:) - 2*psig(1:Nx,Ny,:) + psig([Nx 1:Nx-1],Ny,:)) / dx^2;

    %%% Add to averages
    psig_tavg = psig_tavg + psig;
    ug_tavg = ug_tavg + ug;
    vg_tavg = vg_tavg + vg;
    zetag_tavg = zetag_tavg + zetag;
    uzg_tavg = uzg_tavg + ug.*zetag;
    vzg_tavg = vzg_tavg + vg.*zetag;
    
    %%% Increment counter
    navg = navg + 1;    
    
  end
  

  %%% Divide by number of iterations summed to obtain average 
  %%% N.B. Assumes uniform temporal spacing)
  psig_tavg = psig_tavg/navg;
  ug_tavg = ug_tavg/navg;
  vg_tavg = vg_tavg/navg;
  zetag_tavg = zetag_tavg/navg;
  uzg_tavg = uzg_tavg/navg;      
  vzg_tavg = vzg_tavg/navg;      
  
  %%% Save to a matlab file
  save([run_name,'_QGtavg.mat'], ...
    'tmin','tmax','dt_s', ...
    'ug_tavg','vg_tavg','psig_tavg', ...
    'zetag_tavg','uzg_tavg','vzg_tavg');
  
end

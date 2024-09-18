%%%
%%% calcEddyForcingTimeSeries.m
%%%
%%% Computes time series of spatially-filtered eddy forcing.
%%%

%%% Run to load
run_name = 'ACC_AABW_ML_doubleMOC';

%%% Load parameters   
local_home_dir = '/Volumes/Kilchoman/UCLA/Projects/AWSIM/runs';
prod_dir = fullfile(local_home_dir,'AWSIM_WindAABW_products');
loadParams;
dirpath = fullfile(local_home_dir,run_name);
gtild = reshape(cumsum(gg),[Nlay 1 1]);
rho0 = 1000;
gsum = cumsum(gg);

%%% Max time at which to load transports
tend = 200*t1year;

%%% Define coarse grid
cfac = 4;
Nx_coarse = Nx/cfac;
Ny_coarse = Ny/cfac;

%%% At each time iteration...
cntr = 0;
tt = zeros(1,Nframes);
hh = zeros(Nx,Ny,Nlay);
eta = zeros(Nx,Ny,Nlay+1);
MM = zeros(Nx,Ny,Nlay);
pi = zeros(Nx,Ny);
hh_coarse = zeros(Nx_coarse,Ny_coarse,Nlay,Nframes);
pi_coarse = zeros(Nx_coarse,Ny_coarse,Nframes);
MM_coarse = zeros(Nx_coarse,Ny_coarse,Nlay,Nframes);
eta_coarse = zeros(Nx_coarse,Ny_coarse,Nlay+1,Nframes);
Fx_eddy = zeros(Nx_coarse,Ny_coarse,Nlay,Nframes);
Fy_eddy = zeros(Nx_coarse,Ny_coarse,Nlay,Nframes);
Fx_mean = zeros(Nx_coarse,Ny_coarse,Nlay,Nframes);
Fy_mean = zeros(Nx_coarse,Ny_coarse,Nlay,Nframes);
for n=n0+1:1:n0+Nframes-1   

  %%% Current simulation time    
  t = startTime + (n-n0)*dt_s;
  
  if (t>tend)
    break;
  end
  
  cntr = cntr + 1
  
  tt(n) = t;
 
  %%% Load instantaneous model state
  for k=1:Nlay
    data_file = fullfile(dirpath,[OUTN_H,num2str(k-1),'_n=',num2str(n),'.dat']);
    hh(:,:,k) = readOutputFile(data_file,Nx,Ny);
  end  
%   for k=1:Nlay
%     data_file = fullfile(dirpath,[OUTN_U,num2str(k-1),'_n=',num2str(n),'.dat']);
%     uu(:,:,k) = readOutputFile(data_file,Nx,Ny);
%   end
%   for k=1:Nlay
%     data_file = fullfile(dirpath,[OUTN_V,num2str(k-1),'_n=',num2str(n),'.dat']);
%     vv(:,:,k) = readOutputFile(data_file,Nx,Ny);
%   end
  data_file = fullfile(dirpath,[OUTN_PI,'_n=',num2str(n),'.dat']);
  pi = readOutputFile(data_file,Nx,Ny);
  
  
  %%% Calculate layer surface heights and layer mid-depths
  eta(:,:,1) = 0;
  for k=2:Nlay+1
    eta(:,:,k) = eta(:,:,k-1) - hh(:,:,k-1);
  end
  eta(:,:,Nlay+1) = hhb;      
  for k=Nlay:-1:1               
    eta(:,:,k) = eta(:,:,k+1) + hh(:,:,k);
  end
    
  %%% Calculate Montgomery potential       
  MM(:,:,1) = pi;    
  for k=2:Nlay
    MM(:,:,k) = MM(:,:,k-1) + gg(k)*eta(:,:,k);          
  end     
  for k=1:Nlay
    MM(:,:,k) = MM(:,:,k) - gsum(k) .* h0.^4 ./ hh(:,:,k).^3 ./ 3;
  end
  
  %%% Compute variables on coarse grid
  for i=1:Nx_coarse
    for j=1:Ny_coarse 
      hh_coarse(i,j,:,n) = mean(mean(hh((i-1)*cfac+1:i*cfac,(j-1)*cfac+1:j*cfac,:),1),2);
      pi_coarse(i,j,n) = mean(mean(pi((i-1)*cfac+1:i*cfac,(j-1)*cfac+1:j*cfac),1),2);
      MM_coarse(i,j,:,n) = mean(mean(MM((i-1)*cfac+1:i*cfac,(j-1)*cfac+1:j*cfac,:),1),2);
      eta_coarse(i,j,:,n) = mean(mean(eta((i-1)*cfac+1:i*cfac,(j-1)*cfac+1:j*cfac,:),1),2);
    end
  end
  
  %%% Layer thickness on cell edges
  hh_w = 0.5*(hh(1:Nx,:,:)+hh([Nx 1:Nx-1],:,:)); 
  hh_s = 0.5*(hh(:,1:Ny,:)+hh(:,[Ny 1:Ny-1],:)); 
  hh_w_coarse = 0.5*(hh_coarse(1:Nx_coarse,:,:,n)+hh_coarse([Nx_coarse 1:Nx_coarse-1],:,:,n)); 
  hh_s_coarse = 0.5*(hh_coarse(:,1:Ny_coarse,:,n)+hh_coarse(:,[Ny_coarse 1:Ny_coarse-1],:,n)); 
  
  %%% Products  
  hdMdx = hh_w.*(MM(1:Nx,:,:)-MM([Nx 1:Nx-1],:,:))/dx;
  hdMdy = hh_s.*(MM(:,1:Ny,:)-MM(:,[Ny 1:Ny-1],:))/dy;      
  hdMdx_coarse = hh_w_coarse.*(MM_coarse(1:Nx_coarse,:,:,n)-MM_coarse([Nx_coarse 1:Nx_coarse-1],:,:,n))/dx;
  hdMdy_coarse = hh_s_coarse.*(MM_coarse(:,1:Ny_coarse,:,n)-MM_coarse(:,[Ny_coarse 1:Ny_coarse-1],:,n))/dy;
  
  %%% Eliminate forcing on walls (gradients are spurious there)
  hdMdy(:,1,:) = 0;
  hdMdy_coarse(:,1,:) = 0;
    
  %%% Interpolate to cell centers
  hdMdx = 0.5*(hdMdx(1:Nx,:,:)+hdMdx([2:Nx 1],:,:));
  hdMdy = 0.5*(hdMdy(:,1:Ny,:)+hdMdy(:,[2:Ny 1],:));
  hdMdx_coarse = 0.5*(hdMdx_coarse(1:Nx_coarse,:,:)+hdMdx_coarse([2:Nx_coarse 1],:,:));
  hdMdy_coarse = 0.5*(hdMdy_coarse(:,1:Ny_coarse,:)+hdMdy_coarse(:,[2:Ny_coarse 1],:));
  
  %%% Compute eddy forcing terms
  Fx_mean(:,:,:,n) = hdMdx_coarse;
  Fy_mean(:,:,:,n) = hdMdy_coarse;
  for i=1:Nx_coarse
    for j=1:Ny_coarse       
      Fx_eddy(i,j,:,n) = mean(mean(hdMdx((i-1)*cfac+1:i*cfac,(j-1)*cfac+1:j*cfac,:),1),2);
      Fy_eddy(i,j,:,n) = mean(mean(hdMdy((i-1)*cfac+1:i*cfac,(j-1)*cfac+1:j*cfac,:),1),2);
    end
  end
  Fx_eddy(:,:,:,n) = Fx_eddy(:,:,:,n) - Fx_mean(:,:,:,n);
  Fy_eddy(:,:,:,n) = Fy_eddy(:,:,:,n) - Fy_mean(:,:,:,n);
       
end

%%% Remove missing data
tt = tt(1:cntr);
hh_coarse = hh_coarse(:,:,:,1:cntr);
pi_coarse = pi_coarse(:,:,1:cntr);
MM_coarse = MM_coarse(:,:,:,1:cntr);
eta_coarse = eta_coarse(:,:,:,1:cntr);
Fx_mean = Fx_mean(:,:,:,1:cntr);
Fy_mean = Fy_mean(:,:,:,1:cntr);
Fx_eddy = Fx_eddy(:,:,:,1:cntr);
Fy_eddy = Fy_eddy(:,:,:,1:cntr);

%%% Generate coarse grid
dx_coarse = cfac*dx;
dy_coarse = cfac*dy;
xx_h_coarse = dx_coarse/2:dx_coarse:Lx-dx_coarse/2;
yy_h_coarse = dy_coarse/2:dy_coarse:Ly-dy_coarse/2;
[YY_h_coarse,XX_h_coarse] = meshgrid(yy_h_coarse,xx_h_coarse);

%%% Save to .mat file
save(fullfile(prod_dir,[run_name,'_EddyForcing.mat']), ...
  'cfac','tt','xx_h_coarse','yy_h_coarse','XX_h_coarse','YY_h_coarse', ...
  'hh_coarse','pi_coarse','MM_coarse','eta_coarse', ...
  'Fx_mean','Fy_mean','Fx_eddy','Fy_eddy', ...
  '-v7.3');

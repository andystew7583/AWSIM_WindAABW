%%%
%%% genBathy.m
%%%
%%% Generates a random bathymetry. 
%%%
%%% lambda0 - fundamental wavelength
%%% H0 - mean depth
%%% Hrms - root-mean-square deviation from H0
%%% pow - exponent for high-wavenumber portion of bathymetric spectrum
%%%       (reference value is -3)
%%% Nx,Ny - grid dimensions
%%% Lx,Ly - domain sizes
%%% 
function H = genBathy (lambda0,H0,Hrms,pow,Nx,Ny,Lx,Ly) 

  %%% Load constant parameters
  constants;
  
  %%% Spectral grids
  k = [0:1:Nx/2-1,-Nx/2:1:-1]; %%% Zonal wavenumber
  K_xk = 2*pi.*(k)./Lx;
  l = [0:1:Ny/2-1,-Ny/2:1:-1]; %%% Meridional wavenumber
  K_yl = 2*pi.*(l)./Ly;
  [K_ykl,K_xkl] = meshgrid(K_yl, K_xk);
  
  %%% Physical grids
  x_grid = (0:1:Nx-1);
  y_grid = (0:1:Ny-1);
  x_grid = Lx/Nx .* x_grid;
  y_grid = Ly/Ny .* y_grid;
  
  %%% Fundamental wavenumber  
  K_0 = 2*pi/lambda0; 

  %%% Amplitude is exponential about K0, and phases are random. N.B. here
  %%% we only define the amplitude up to a constant - below we constrain it
  %%% based on the RMS KE.
  K = sqrt(K_xkl.^2 + K_ykl.^2);
  theta = 2 .* pi .* rand(Nx,Ny);
  H_fft = K.^(-1/2).*(K + 4*K_0).^(pow) .* exp(1i*theta);

  %%% Avoids infinite mode-0 amplitude 
  H_fft(1,1) = 0;

  %%% Transform to real space
  H = Nx*Ny*real(ifft2(H_fft));
  
  %%% Normalize and add reference depth
  H = H * Hrms/rms(H(:)) + H0;   

end
%%%
%%% calcQGStreamfunction.m
%%%
%%% Computes baroclinic modes for a 3-layer shallow water fluid.
%%%
%%% Inputs:
%%%   pi      Nx x Ny matrix of surface dynamic pressures.
%%%   eta     Nx x Ny x (Nlay+1) matrix of layer surface elevations.
%%%   HH      Vector (length: Nlay) of reference layer thicknesses.
%%%   gg      Vector (length: Nlay) of reduced gravities. The first element is
%%%           assumed to correspond to the upper surface of the ocean, so
%%%           it is ignored under the rigid lid approximation.
%%%   f0      Reference Coriolis parameter.
%%%   
%%% Outputs:
%%%   psig    Nx x Ny x Nlay matrix of QG streamfunctions.
%%%   etag    Nx x Ny x (Nlay+1) matrix of QG layer interfaces.
%%%
function [psig,etag] = calcQGStreamfunction (pi,eta,HH,gg,f0)

  %%% Grids
  Nx = size(eta,1);
  Ny = size(eta,2);
  Nlay = size(eta,3)-1;

  %%% Error checking
  if ((size(pi,1) ~= Nx) || (size(pi,2)~=Ny) || (size(pi,3)~=1) || (length(HH)~=Nlay) || (length(gg)~=Nlay))
    error('Incompatible input dimensions in calcQGStreamfunction');
  end

  %%% Geostrophic layer surface elevations
  etag = eta;
  for k=2:Nlay+1
    etag(:,:,k) = etag(:,:,k)+sum(HH(1:k-1));
  end

  %%% Calculate QG streamfunction in each layer
  psig(:,:,1) = pi./f0;
  for k=2:Nlay
    psig(:,:,k) = psig(:,:,k-1) + (gg(k)./f0) .* etag(:,:,k);
  end

end


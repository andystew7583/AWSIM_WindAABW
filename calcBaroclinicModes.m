%%%
%%% calcBaroclinicModes.m
%%%
%%% Computes baroclinic modes for a 3-layer shallow water fluid.
%%%
%%% Inputs:
%%%   psig    Nx x Ny x 3 matrix of geostrophic streamfunction (or any other
%%%           layerwise quantity).
%%%   HH      Vector (length: 3) of reference layer thicknesses.
%%%   gg      Vector (length: 3) of reduced gravities. The first element is
%%%           assumed to correspond to the upper surface of the ocean, so
%%%           it is ignored under the rigid lid approximation.
%%%   
%%% Outputs:
%%%   Mpsig   Baroclinic mode amplitudes
%%%   EE      Matrix of baroclinic mode structure functions
%%%   eps     3 x 3 x 3 interaction matrix for 3-mode interactions.
%%%
function [Mpsig,EE,eps] = calcBaroclinicModes (psig,HH,gg)  

  %%% Error checking
  if ((size(psig,3) ~= 3) || (length(HH)~=3) || (length(gg)~=3))
    error('calcBaroclinicModes requires exactly 3 isopycnal layers');
  end

  %%% Grid dimensions
  Nx = size(psig,1);
  Ny = size(psig,2);
  Nlay = 3;
  HH = reshape(HH,[3 1]);  

  %%% Solve for eigenvalues of baroclinic mode decomposition (inverse wave
  %%% speeds)
  Gamma = 0*HH;
  Gamma_vert = 1/(gg(2))*(1/HH(1)+1/HH(2)) + 1/gg(3)*(1/HH(2)+1/HH(3));
  Gamma_disc = ( 1/gg(3)*(1/HH(2)+1/HH(3)) - 1/(gg(2))*(1/HH(1)+1/HH(2)) )^2 + 4/(gg(2)*gg(3)*HH(2)^2);
  Gamma(1) = 0;
  Gamma(2) = 0.5* (Gamma_vert - sqrt(Gamma_disc));
  Gamma(3) = 0.5* (Gamma_vert + sqrt(Gamma_disc));

  %%% Construct vector of baroclinic modal structures
  EE = zeros(3,3);
  for n=1:3
    vv = [ 1/(1-Gamma(n)*gg(2)*HH(1)) ; ...
           1 ; ...
           1/(1-Gamma(n)*gg(3)*HH(3)) ];
    EE(:,n) = vv / sqrt(sum(vv.^2.*HH)); %%% Normalize so inner product of modes mn and n is delta_mn
  end
  
  %%% Construct the mode interaction matrix
  eps = zeros(Nlay,Nlay);
  for n = 1:Nlay
    for m = 1:Nlay
      for l = 1:Nlay
        eps(l,m,n) = sum(EE(:,l).*EE(:,m).*EE(:,n).*HH);
      end
    end
  end

  %%% Decompose streamfunction and vorticity into baroclinic modes
  Mpsig = zeros(Nx,Ny,Nlay);
  for i=1:Nx
    for j=1:Ny            
      Mpsig(i,j,:) = reshape(EE\squeeze(psig(i,j,:)),[1 1 Nlay]);  
    end
  end

end




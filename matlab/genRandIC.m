%%%
%%% genRandIC.m
%%%
%%% Generates a random initial condition for the geostrophic streamfunction
%%% with characteristic spectral wavelength 'lambda' and kinetic energy 
%%% density E0. Nx and Ny define the grid size, and Ly the meridional 
%%% domain length. geff is the effective gravity at the surface.
%%% 
function psi = genRandIC (lambda,E0,Nx,Ny,Lx,Ly) 

  %%% Load constant parameters
  constants;
  
  %%% Spectral grids 
  k = [0:1:ceil(Nx/2)-1,-floor(Nx/2):1:-1]; %%% Zonal wavenumber
  K_xk = 2*pi.*(k)./Lx;
  l = [0:1:ceil(Ny/2)-1,-floor(Ny/2):1:-1]; %%% Meridional wavenumber
  K_yl = 2*pi.*(l)./Ly;
  [K_xkl, K_ykl] = meshgrid(K_yl, K_xk);
  
  
  %%% Physical grids
  x_grid = (0:1:Nx-1);
  y_grid = (0:1:Ny-1);
  x_grid = Lx/Nx .* x_grid;
  y_grid = Ly/Ny .* y_grid;
  
  %%% Parameters defining shape of energy spectrum
  K_0 = 2*pi/lambda; %%% Most energetic wavenumber
  W = K_0/8; %%% Exponential width of energy band in wavenumber space

  %%% Amplitude is exponential about K0, and phases are random. N.B. here
  %%% we only define the amplitude up to a constant - below we constrain it
  %%% based on the RMS KE.
  K = sqrt(K_xkl.^2 + K_ykl.^2);
  theta = 2 .* pi .* rand(Nx,Ny);
  psi_fft = K.^(-1).*exp(-((K-K_0)/W).^2) .* exp(1i*theta);
  
  %%% Avoids infinite mode-0 amplitude 
  psi_fft(1,1) = 0;

  %%% Spectral energy
  E_fft = 0.25.*K.^2.*abs(psi_fft).^2;
  E_mean = sum(sum(E_fft));
  
  %%% This ratio that when multiplied onto initial_E will give an E matrix with E_rms = E0
  alphaE = E0/E_mean;                     

  %%% Renormalize eta in spectral space to ensure E_rms=E0
  psi_fft = sqrt(alphaE).*psi_fft;
  
  %%% Transform back to real space
  psi = Nx*Ny*real(ifft2(psi_fft));


  %plot(K_xk);
  %p = surf(etafft);
  %p = surf(K);
  %p = surf(real(ifft2(etafft)));
  %p = surf(X,Y, real(ifft2(etafft)));
  %shading interp;
  %p = pcolor(X,Y,real(ifft2(etafft)));
  %set(p,'LineStyle','none');


end
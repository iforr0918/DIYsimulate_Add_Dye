%%%
%%% setparams.m
%%%
%%% Sets parameters for DIYsimulate experiments. 
%%% 
function setparams (local_home_dir,run_name)
  
  %%% Load common matlab scripts/functions
  addpath ../matlab_common;
  
  %%% Load constant parameters
  constants;
  
  %%% Run directory
  run_name = strtrim(run_name); 
  local_home_dir = strtrim(local_home_dir); 
  local_run_dir = fullfile(local_home_dir,run_name);
  mkdir(local_run_dir);
  pfname = fullfile(local_run_dir,[run_name,'_in']); 
  
  %%% Cluster config
  use_cluster = false;
  use_intel = use_cluster;
  uname = 'astewart';
  cluster_addr = 'tip.atmos.ucla.edu';
  cluster_home_dir = 'Desktop/QGSWchannel/runs';
   

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%% PARAMETERS %%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  %%% To store parameters
  paramTypes;
  PARAMS = {};
            
  amult = 1; %%% Number of subdivisions of the annulus
  rmin = 0.17/10; %%% Inner wall radius
  rmax = 0.172; %%% Outer wall radius
  dr = 0.0025; %%% Radial grid spacing in m
%   dr = 0.001; %%% Radial grid spacing in m
  H = 0.05; %%% Water depth in m
  Nr = ceil((rmax-rmin)/dr) + 1 %%% r-gridpoints
%   Na = ceil(2*pi*(rmax+rmin)/2/amult/dr) %%% theta-gridpoints
  Na = ceil(2*pi*sqrt(rmax*rmin)/amult/dr) %%% theta-gridpoints
  tmax = 180; %%% Integration time    
  nu = 1e-6; %%% Actual fluid viscosity  
  f = 1; %%% Background rotation rate (rad/s)      
  
  %%% Tracer parameters
%   Ntracs_r = 30;
%   Ntracs_a = 180;
  Ntracs_r = 40;
  Ntracs_a = 200;
  
  %%% Spin-down parameters
  lambdaK = 0.08; %%% Initial eddy wavelength (m)
%   lambdaK = 0.04; %%% Initial eddy wavelength (m)
%   E0 = 0.00005; %%% Initial eddy energy (m^s/s^2)
   E0 = 0; %%% Initial eddy energy (m^s/s^2)
  
  %%% Azimuthal flow parameters
  zeta0 = -0.2; %%% Initial relative vorticity
  %zeta0 = 0; %%% Initial relative vorticity
  psi0Init = 0.25*zeta0*(rmax^2-rmin^2); %%% Initial along-channel transport (c.f. Stewart et al. 2014)
  
  %%% Numerical viscosity is chosen to be order 1 at the grid scale over 
  %%% one dynamical time scale (1/f), but cannot be smaller than the
  %%% actual fluid viscosity
  Ah = max(abs(f)/5*dr^2,nu)
  
  %%% Bottom drag based on inverse Ekman spindown time using actual
  %%% viscosity
  kappa = sqrt(nu*f)/H
%   kappa = 0;
  
  %%% Frequency of model output in s
  savefreq = 0.5;
  savefreqP = 0.5;
        
  %%% Grids  
  da = 2*pi/amult/Na;
  rr = rmin:dr:rmax;
  aa = 0:da:2*pi/amult-da;    
  [AA,RR] = meshgrid(aa,rr);  
  XX = RR.*cos(AA);
  YY = RR.*sin(AA);
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% OTHER PARAMETERS %%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
     
  %%% Estimate maximum flow speed
  Umax = max(sqrt(2*E0),2*abs(psi0Init/(rmax-rmin)));

  %%% Set timestep   
  dt = 0.25*min(dr,rmin*da)/(Umax);    
    
  %%% Correct dt so that tmax is an integer number of time steps.    
  Nt = ceil(tmax/dt)+1;      
         
  %%% Define parameters
  PARAMS = addParameter(PARAMS,'Nr',Nr,PARM_INT);
  PARAMS = addParameter(PARAMS,'Na',Na,PARM_INT);  
  PARAMS = addParameter(PARAMS,'Nt',Nt,PARM_INT);
  PARAMS = addParameter(PARAMS,'Np',Ntracs_r*Ntracs_a,PARM_INT);
  PARAMS = addParameter(PARAMS,'tmax',tmax,PARM_REALF);
  PARAMS = addParameter(PARAMS,'savefreq',savefreq,PARM_INT);
  PARAMS = addParameter(PARAMS,'savefreqP',savefreqP,PARM_INT);  
  PARAMS = addParameter(PARAMS,'rmin',rmin,PARM_REALF);
  PARAMS = addParameter(PARAMS,'rmax',rmax,PARM_REALF);
  PARAMS = addParameter(PARAMS,'f',f,PARM_REALF); 
  PARAMS = addParameter(PARAMS,'H',H,PARM_REALF);  
  PARAMS = addParameter(PARAMS,'amult',amult,PARM_INT);
  PARAMS = addParameter(PARAMS,'kappa',kappa,PARM_REALF);      
  PARAMS = addParameter(PARAMS,'nu',Ah,PARM_REALF);      
  PARAMS = addParameter(PARAMS,'transportInit',psi0Init,PARM_REALF);      
    
  
  %%%%%%%%%%%%%%%%%%%%%%
  %%%%% BATHYMETRY %%%%%
  %%%%%%%%%%%%%%%%%%%%%%
  
  %%% Nothing
  hh = zeros(Nr,Na);
  
  %%% Ridge
%   Hb = 0.1;
%   Wb = pi/16;
%   Ab = pi;
%   hh = Hb*exp(-((AA-Ab)/Wb).^2);

  %%% Isolated bump
%   X_bump = 0;
%   Y_bump = 0.3;
%   W_bump = 0.05;
%   H_bump = 0.05;
%   hh = H_bump*exp(-((XX-X_bump)/W_bump).^2-((YY-Y_bump)/W_bump).^2);
%   

  %%% TANH PUCK
  X_puck = 0;
  Y_puck = (rmax -rmin)/2;
  H_puck = H/8;
  Rad_puck = (rmax -rmin)/8;
  puck_slope_coef = 250;
  d_to_puck = sqrt((XX-X_puck).^2+(YY-Y_puck).^2);
  hh = (H_puck/2)*(tanh(puck_slope_coef.*(Rad_puck - d_to_puck))+1);
  mesh(XX, YY, hh)
  zlim([0 H])
  %view(90,5)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% DYE INITIAL COND %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%% Nothing
  redred = zeros(Nr,Na);
  
  %%% Isolated exp red splat
   X_splat = 0;
   Y_splat = -(rmax-rmin)/2;;
   W_splat = 0.02;
   H_splat = 0.5;
   redred = H_splat*exp(-((XX-X_splat)/W_splat).^2-((YY-Y_splat)/W_splat).^2);
   %enforce BC
   redred(1,:) = 0;
   redred(end,:) = 0;
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% INITIAL CONDITIONS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%% Generate random streamfunction on a square grid
  xx_ic = -rmax:dr/4:rmax;
  yy_ic = -rmax:dr/4:rmax;
  [YY_ic,XX_ic]=meshgrid(yy_ic,xx_ic);
  N_ic = length(xx_ic);
  psi = genRandIC(lambdaK,E0,N_ic,N_ic,2*rmax,2*rmax);  
  
  %%% Interpolate to annulus
  psi = interp2(XX_ic',YY_ic',psi,XX,YY,'linear');
  
  %%% Enforce boundary conditions, and taper close to boundaries
  psi =  (1 - exp(-((RR-rmin)./(lambdaK/2)).^2)) .* (1 - exp(-((RR-rmax)./(lambdaK/2)).^2)) .* psi;
  psi(1,:) = 0;
  psi(end,:) = 0;
    
  %%% Plot sea surface
  figure(2);
  contourf(XX,YY,psi);
  colorbar;
  
  %%% Set relative vorticity
  %%% NOTE: Stays zero on the boundaries to enforce no-stress BC
  vort = 0*psi;
  vort(2:Nr-1,:) = ( psi(3:Nr,:) - 2*psi(2:Nr-1,:) + psi(1:Nr-2,:) ) / dr^2 ...
                   + ( psi(3:Nr,:) - psi(1:Nr-2,:) ) ./ RR(2:Nr-1,1:Na) ./ (2*dr) ...
                   + ( psi(2:Nr-1,[2:Na 1]) - 2*psi(2:Nr-1,1:Na) + psi(2:Nr-1,[Na 1:Na-1]) ) ./ RR(2:Nr-1,1:Na).^2 ./ da^2;
  
  %%% Add initial flow component
  vort = vort + zeta0;
  
  %%% Boundary conditions
  vort([1 Nr],:) = 0;
         
  %%% Plot vorticity
  figure(3);  
  contourf(XX,YY,vort);   
  colorbar;
  colormap redblue;
  
  %%% Calculate kinetic energy of initial state
  meanE = 0.5*mean(mean(psi.*vort))
  
  %%% Initial tracer positions
  tracPos = zeros(Ntracs_r*Ntracs_a,2);
  dr_trac = (rmax-rmin)/Ntracs_r;
  da_trac = 2*pi/amult/Ntracs_a;
  for i=1:Ntracs_r
    for j=1:Ntracs_a      
      tracPos((i-1)*Ntracs_a+j,1) = rmin + (i-0.5)*dr_trac;
      tracPos((i-1)*Ntracs_a+j,2) = j*da_trac;
    end
  end
  
  figure(4)
  tracX = tracPos(:,1).*cos(tracPos(:,2));
  tracY = tracPos(:,1).*sin(tracPos(:,2));
  scatter(tracX,tracY,'.');
  
  
  
  
  
  
  
  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% CREATE INPUT FILES %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
  
%   %%% Initial streamfunction 
%   psiInitFile = 'psiInit.dat';
%   writeDataFile(fullfile(local_run_dir,psiInitFile),psi);
%   PARAMS = addParameter(PARAMS,'psiInitFile',psiInitFile,PARM_STR);  
  
  %%% Initial vorticity 
  vortInitFile = 'vortInit.dat';
  writeDataFile(fullfile(local_run_dir,vortInitFile),vort);
  PARAMS = addParameter(PARAMS,'vortInitFile',vortInitFile,PARM_STR);  

  %%% Bathymetry  
  hhFile = 'hh.dat';          
  writeDataFile(fullfile(local_run_dir,hhFile),hh); 
  PARAMS = addParameter(PARAMS,'bathyFile',hhFile,PARM_STR); 

  %%% Tracer  
  tracFile = 'tracInitFile.dat';          
  writeDataFile(fullfile(local_run_dir,tracFile),tracPos'); 
  PARAMS = addParameter(PARAMS,'tracInitFile',tracFile,PARM_STR); 
  
  %%% Red Tracer  
  redTracFile = 'redTracInitFile.dat';          
  writeDataFile(fullfile(local_run_dir,redTracFile),redred'); 
  PARAMS = addParameter(PARAMS,'redTracInitFile',redTracFile,PARM_STR); 

  %%% Create the input parameter file
  writeParamFile(pfname,PARAMS);    
  

  %%% Create a run script
  createRunScript (  local_home_dir, ...
                     run_name, ...
                     model_code_dir, ...
                     exec_name, ...
                     use_intel, ...
                     false, ... 
                     use_cluster, ...
                     uname, ...
                     cluster_addr, ...
                     cluster_home_dir);

end

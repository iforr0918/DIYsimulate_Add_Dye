%%%
%%% loadParams.m
%%%
%%% Loads a set of typically-used parameters for QGSWchannel model runs. The
%%% variables 'run_name' and 'local_home_dir' must be specified before calling
%%% this script.
%%%

%%% Load constant parameter definitions
constants;

%%% Parameter and data file names
run_name = strtrim(run_name);
dirpath = fullfile(local_home_dir,run_name);
params_file = fullfile(dirpath,[run_name,'_in']);
t_file = fullfile(dirpath,OUTN_TFILE);

%%% Grids
[Nr Nr_found] = readparam(params_file,'Nr','%u');
[Na Na_found] = readparam(params_file,'Na','%u');  
[rmin rmin_found] = readparam(params_file,'rmin','%lf');  
[rmax rmax_found] = readparam(params_file,'rmax','%lf');  
[amult amult_found] = readparam(params_file,'amult','%u');  
if ((~Nr_found) || (~Na_found) || (~rmin_found) || (~rmax_found))
  error('Could not read spatial grid parameters');
end
if (~amult_found)
  amult = 1;
end
dr = (rmax-rmin)/(Nr-1);
da = 2*pi/amult/Na;
rr = rmin:dr:rmax;
aa = 0:da:2*pi/amult-da;    
[AA,RR] = meshgrid(aa,rr);  
XX = RR.*cos(AA);
YY = RR.*sin(AA);

%%% Parameters related to number of iterations
[Nt Nt_found] = readparam(params_file,'Nt','%u');
[tmax tmax_found] = readparam(params_file,'tmax','%lf');
[dt_s dt_s_found] = readparam(params_file,'savefreq','%lf');
[dt_p dt_p_found] = readparam(params_file,'savefreqP','%lf');
if ((~Nt_found) || (~tmax_found))
  error('Could not read temporal grid parameters');
end
if (~dt_s_found)
  dt_s = tmax/Nt; %%% Defaults to the time step, dt
end
if (~dt_p_found)
  dt_p = dt_s; %%% Defaults to dt_s
end
Nframes = round(tmax / dt_s) + 1;

%%% Tracers
[Np Np_found] = readparam(params_file,'Np','%u');  
if (~Np_found)
  Np = 0;
end

%%% Tank depth
[H H_found] = readparam(params_file,'H','%lf'); 
if (~H_found)
  H = 1000;
end 

%%% Coriolis parameter
[f f_found] = readparam(params_file,'f','%lf'); 
if (~f_found)
  f = 0;
end

%%% Drag coefficient
[kappa kappa_found] = readparam(params_file,'kappa','%lf'); 
if (~kappa_found)
  kappa = 0;
end      

%%% Load topography
hh = readDataFile(params_file,dirpath,'bathyFile',Nr,Na,zeros(Nr,Na)); 
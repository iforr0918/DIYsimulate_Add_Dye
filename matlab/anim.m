function M = anim (local_home_dir,run_name,t_start,t_end, ...
            field,contours,show_tracers,show_velocity,show_topog)
%%%
%%% USAGE: anim (local_home_dir, run_name, t_start, t_end, field, contours, show_tracers, show_velocity, show_topog)
%%%
%%% Creates an animation of the 'QGChannel' simulation output.
%%%
%%% Arguments:
%%% local_home_dir - Path to directory that contains simulation directory,
%%%                  e.g. if the simulation files are in /Desktop/runs/test_run 
%%%                  then this parameter should be '/Desktop/runs/'.
%%% run_name - Name of simulation directory, e.g. in previous example, this
%%%            would be 'test_run'.
%%% t_start - Simulation time (in seconds) at which to start the movie.
%%% t_end - Simulation time (in seconds) at which to end the movie.
%%% field - Selects what to plot. Options are 
%%%                'zeta' (relative vorticity)
%%%                'pv' (potential vorticity)
%%%                'psi' (streamfunction)
%%%                'speed' (absolute flow speed)
%%% contours - Vector of contour levels to plot. If just an integer number
%%%            N is specified then N arbitrary contours will be plotted.
%%% show_tracers - If set 'true' then the passive tracer particles in the
%%%                simulation will be displayed. If set to 'false' then
%%%                they will not.
%%% show_velocity - If set 'true' then the velocity vectors will be
%%%                 indicated using arrows.
%%% show_topog - If set 'true' then white lines will be plotted to indicate
%%%              the top and bottom of the continental slope
%%%

  %%% Load model parameters
  loadParams;
 
  %%% Plotting options
  scrsz = get(0,'ScreenSize');
  fontsize = 22;
  plotloc = [0.07 0.07 0.88 0.87];
  framepos = [scrsz(3)/4 0 scrsz(4) scrsz(4)];

  %%% Make a movie of the data
  handle = figure(1);
  set(handle,'Position',framepos);
  clf;
  axes('FontSize',fontsize);
  M = moviein(Nframes);
  
  %%% At each time iteration...
  for n=1:1:Nframes
    
    %%% Current simulation time
    t = n*dt_s;    
    if ((t_start > 0 && t < t_start) || ((t_end>0) && (t > t_end)))
      continue;
    end        
    
    %%% Storage
    psi = [];
    vort = [];
    pv = [];
    
    %%% Read tracer positions
    if (Np > 0)
      data_file = fullfile(dirpath,[OUTN_TRACER,'_n=',num2str(n),'.dat']);      
      tracers = readOutputFile(data_file,2,Np);                 
    end
    
    %%% Load streamfunction
    data_file = fullfile(dirpath,[OUTN_PSI,'_n=',num2str(n),'.dat']);
    psi = readOutputFile(data_file,Nr,Na);
    
    %%% Stop if we've run out of data
    if (isempty(psi))
      break;
    end
    
    %%% Calculate the flow velocity
    ur = - (psi(:,[2:Na 1])-psi(:,[Na 1:Na-1])) ./ RR ./ da;
    ua = zeros(Nr,Na);
    ua(2:Nr-1,:) = (psi(3:Nr,1:Na)-psi(1:Nr-2,1:Na)) / dr;    
    uu = ur.*cos(AA) - ua.*sin(AA);
    vv = ua.*cos(AA) + ur.*sin(AA);
    

    %%% Modify tracer zonal positions so that they lie in [0 2pi/amult)
    for p=1:Np      
      tracers(2,p) = tracers(2,p) - 2*pi/amult * floor(tracers(1,p)/(2*pi/amult));                                        
    end    
        
    %%% Used to check that the specified field exists
    str_match = false;
    
    %%% Plot relative vorticity
    if (strcmp(field,'zeta'))
      
      %%% Load potential vorticity
      data_file = fullfile(dirpath,[OUTN_PV,'_n=',num2str(n),'.dat']);
      pv = readOutputFile(data_file,Nr,Na);  
      
      %%% Compute the relative vorticity 
      vort =  pv - f*hh/H;                     

      str_match = true;
%       contourf(XX/1000,YY/1000,vort/abs(f0),contours,'EdgeColor','None');              
      pcolor(XX,YY,vort/abs(f));
      shading interp;
      vortmax = max(max(abs(vort/abs(f))));
      caxis([-vortmax vortmax]);
      colormap(cmocean('balance'));
      title(strcat(['Vortex Rossby number at t=',num2str(t,'%3.1f'),' seconds']),'FontSize',fontsize,'Interpreter','latex');
      
    end
    
    %%% Plot potential vorticity
    if (strcmp(field,'pv'))
      
      %%% Load potential vorticity
      data_file = fullfile(dirpath,[OUTN_PV,'_n=',num2str(n),'.dat']);
      pv = readOutputFile(data_file,Nr,Na);
      
      str_match = true;
%       contourf(XX/1000,YY/1000,pv,contours,'EdgeColor','None');              
      pcolor(XX,YY,pv);
      shading interp;
      colormap(cmocean('balance'));
      title(strcat(['Potential vorticity (s$^{-1}$) at t=',num2str(t,'%3.1f'),' seconds']),'FontSize',fontsize,'Interpreter','latex');
      
    end
    
    %%% Plot streamfunction
    if (strcmp(field,'psi'))
               
      str_match = true;
%       contourf(XX,YY,psi,contours,'EdgeColor','None');              
      pcolor(XX,YY,psi);
      shading interp;
      colormap jet;
      title(strcat(['Streamfunction (m$^2$/s) at t=',num2str(t,'%3.1f'),' seconds']),'FontSize',fontsize,'Interpreter','latex');
%       clim 
    end
    
    %%% Plot zonal velocity
    if (strcmp(field,'uvel'))

      str_match = true;
      contourf(XX_u,YY_u,ur,contours,'EdgeColor','None');              
      colormap redblue;
      title(strcat(['Zonal velocity (m/s) at t=',num2str(t,'%3.1f'),' seconds']),'FontSize',fontsize,'Interpreter','latex');
      
    end
    
    %%% Plot absolute speed
    if (strcmp(field,'speed'))
      
      %%% Calculate flow speed
      speed = sqrt(uu.^2+vv.^2);
      
      str_match = true;
%       contourf(XX,YY,speed,contours,'EdgeColor','None');              
      pcolor(XX,YY,speed);
      shading interp;
      colormap hot;
      title(strcat(['Flow speed (m/s) at t=',num2str(t,'%3.1f'),' seconds']),'FontSize',fontsize,'Interpreter','latex');
      
    end
    
    %%% Restrict colorbar limits
    if (length(contours)>1)
      caxis([min(contours),max(contours)]);
    end
      
    
    %%% Not a recognized field name
    if (~str_match)
      error(['Field not recognized: ',field]);
    end
               
    %%% Contour the velocity vectors
    if (show_velocity)
      hold on;
      substep = 5;
%       substep = 10;
      quiver(XX(1:substep:Nr,1:substep:Na),YY(1:substep:Nr,1:substep:Na),uu(1:substep:Nr,1:substep:Na),vv(1:substep:Nr,1:substep:Na),1.5,'Color','k');
      hold off;
    end
    
    %%% Plot the top and bottom of the continental slope
    if (show_topog)
      hold on;    
      [C,h] = contour(XX,YY,hh,[0.02:0.02:0.1],'EdgeColor','k');
      clabel(C,h);
      hold off;
    end
    
    %%% Plot walls      
    hold on;    
    plot(XX(1,:),YY(1,:),'k-','LineWidth',3);
    plot(XX(Nr,:),YY(Nr,:),'k-','LineWidth',3);
    hold off;
    
    %%% Plot tracer particles
    if (show_tracers && (Np > 0))
      hold on;      
%       handle = plot(tracers(1,:)/1000,tracers(2,:)/1000,'w.');
%       tracX = tracers(:,1).*cos(tracers(:,2));
%       tracY = tracers(:,1).*sin(tracers(:,2));
      tracX = tracers(1,:).*cos(tracers(2,:));
      tracY = tracers(1,:).*sin(tracers(2,:));
%       scatter(tracX,tracY,'.');
      handle = plot(tracX,tracY,'k.');
%       set(handle,'MarkerSize',10);    
      hold off;
    end
    
    %%% Add other labels to the plot    
%     axis([0 Lx/1000 0 max(yy)/1000]);
%     axis([150 250 250 350]);
    xlabel('x (m)','Interpreter','latex');
    ylabel('y (m)','Interpreter','latex');        
    set(gca,'Position',plotloc);
    handle = colorbar;    
    set(handle,'FontSize',fontsize);    
    set(gca,'FontSize',fontsize);    
    nextframe = getframe(gcf);
    M(n) = nextframe;                
    
  end
    
end
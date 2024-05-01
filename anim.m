%%%
%%% anim.m
%%%
%%% Reads in data from the output of 'AWSIM' and makes a movie of the
%%% layer output.
%%%
function M = anim (local_home_dir,run_name,var,layer,tmin,tmax) 

  %%% Load parameters   
  loadParams;
  dirpath = fullfile(local_home_dir,run_name);
  
  %%% Make a movie of the data
  clf;
  axes('FontSize',16);
  M = moviein(1);      
  Mcnt = 1;
  
  %%% At each time iteration...
  for n=n0:1:n0+Nframes-1   
    
    %%% Current simulation time    
    t = startTime + (n-n0)*dt_s;
    tt(Mcnt) = t;
    
    if ((tmin >= 0 && t < tmin) || (tmax >= 0 && t > tmax))
      continue;
    end
    n
    
    %%% Set up figure
    figure(3); 
        
    switch (var)
      
      %%% Contour plot of relative vorticity
      case 'z'                
        
        %%% Load u
        data_file = fullfile(dirpath,[OUTN_U,num2str(layer-1),'_n=',num2str(n),'.dat']);
        uu = readOutputFile(data_file,Nx,Ny);
        
        %%% Load v
        data_file = fullfile(dirpath,[OUTN_V,num2str(layer-1),'_n=',num2str(n),'.dat']);
        vv = readOutputFile(data_file,Nx,Ny);            
        
        %%% Load h
        data_file = fullfile(dirpath,[OUTN_H,num2str(layer-1),'_n=',num2str(n),'.dat']);
        hh = readOutputFile(data_file,Nx,Ny); 
                     
        %%% Calculate the relative vorticity
        zeta = zeros(Nx+1,Ny+1);            
        zeta(2:Nx,2:Ny) = (uu(1:Nx-1,1:Ny-1)-uu(1:Nx-1,2:Ny)) / dy + (vv(2:Nx,2:Ny)-vv(1:Nx-1,2:Ny)) / dx;                      
        if (~useWallNS)
          zeta(:,Ny+1) = zeta(:,1);
        end
        zeta(Nx+1,:) = zeta(1,:);
        Ro = zeta./(2*Omega_z);
        max(max(abs(Ro)))        

        %%% Make the plot
        pcolor(XX_q/1000,YY_q/1000,Ro);
        shading interp;
        hold on
        contour(XX_h/1000,YY_h/1000,hhb,[-4000:500:-500],'EdgeColor','k');
        hold off;
        colorbar;
        colormap(cmocean('balance'));
        caxis([-0.4 0.4]);        
        title(strcat(['Ro=\zeta/|f0| at t=',num2str(t/t1year,'%.2f'),' years']));        
        xlabel('x (km)');
        ylabel('y (km)');                   
        set(gca,'FontSize',14);
        
      %%% Contour plot of potential vorticity
      case 'pv'
        
        %%% Load u
        data_file = fullfile(dirpath,[OUTN_U,num2str(layer-1),'_n=',num2str(n),'.dat']);
        uu = readOutputFile(data_file,Nx,Ny);
        
        %%% Load v
        data_file = fullfile(dirpath,[OUTN_V,num2str(layer-1),'_n=',num2str(n),'.dat']);
        vv = readOutputFile(data_file,Nx,Ny);    
        
        %%% Load h
        data_file = fullfile(dirpath,[OUTN_H,num2str(layer-1),'_n=',num2str(n),'.dat']);
        hh = readOutputFile(data_file,Nx,Ny);      
        
        %%% Calculate the relative vorticity
        zeta = zeros(Nx+1,Ny+1);            
        zeta(2:Nx,2:Ny) = (uu(1:Nx-1,1:Ny-1)-uu(1:Nx-1,2:Ny)) / dy + (vv(2:Nx,2:Ny)-vv(1:Nx-1,2:Ny)) / dx;                                         
        
        %%% Calculate PV        
        ff = 2*Omega_z;        
        hh_q = NaN*ones(Nx+1,Ny+1);
        hh_q(2:Nx,2:Ny) = 0.25 * (hh(2:Nx,1:Ny-1)+hh(2:Nx,2:Ny)+hh(1:Nx-1,1:Ny-1)+hh(1:Nx-1,2:Ny));        
        pv = (ff+zeta)./hh_q;                   
                
        %%% Make the plot
        pcolor(XX_q/1000,YY_q/1000,pv);
        pcolor(XX_q/1000,YY_q/1000,log10(abs(pv)));
        shading interp;
        colorbar;
%         colormap(pmkmp(100,'Swtth'));
        colormap(cmocean('balance'));
        title(strcat(['PV at t=',num2str(t/t1year,'%.2f'),' years']));        
        xlabel('x (km)');
        ylabel('y (km)');         
        
      %%% Plot layer surface height
      case 'e'
        
        %%% Calculate layer surface height
        eta = hhb;
        for k=Nlay:-1:layer

          %%% Load kth layer thickness
          data_file = fullfile(dirpath,[OUTN_H,num2str(k-1),'_n=',num2str(n),'.dat']);
          hh = readOutputFile(data_file,Nx,Ny);      
          
          %%% Add layer thickness to eta
          eta = eta + hh;
          
        end
        
        %%% Make the plot
        pcolor(XX_h/1000,YY_h/1000,eta);
        shading interp;
        colorbar;
        title(strcat(['t=',num2str(t/t1day,'%.2f'),' days']));
        colormap jet;      

      %%% Plot planetary pv
      case 'w'
        
        %%% Load h
        data_file = fullfile(dirpath,[OUTN_W,num2str(layer-1),'_n=',num2str(n),'.dat']);
        ww = readOutputFile(data_file,Nx,Ny);    
    
        %%% Make the plot
        pcolor(XX_q/1000,YY_q/1000,ww);
        shading interp;
        colorbar;
        title(strcat(['t=',num2str(t/t1day,'%.2f'),' days']));
        colormap jet;       


      %%% Plot planetary pv
      case 'qp'
        
        %%% Load h
        data_file = fullfile(dirpath,[OUTN_H,num2str(layer-1),'_n=',num2str(n),'.dat']);
        hh = readOutputFile(data_file,Nx,Ny);    
        
        %%% Calculate planetary PV        
        ff = 2 * 2*Omega_z;
        hh_q = 0*ff;
        hh_q(1:Nx,2:Ny) = 0.25 * (hh(1:Nx,1:Ny-1)+hh(1:Ny,2:Ny)+hh([Nx 1:Nx-1],1:Ny-1)+hh([Nx 1:Nx-1],2:Ny));
        hh_q(1:Nx,Ny+1) = 0.5 * (hh(1:Ny,Ny)+hh([Nx 1:Nx-1],Ny));
        hh_q(1:Nx,1) = 0.5 * (hh(1:Nx,1)+hh([Nx 1:Nx-1],1));
        qp = ff./hh_q;
    
        %%% Make the plot
        pcolor(XX_q/1000,YY_q/1000,qp);
        shading interp;
        colorbar;
        title(strcat(['t=',num2str(t/t1day,'%.2f'),' days']));
        colormap jet;       

      %%% Plot layer thickness
      case 'h'
        
        %%% Load h
        data_file = fullfile(dirpath,[OUTN_H,num2str(layer-1),'_n=',num2str(n),'.dat']);
        hh = readOutputFile(data_file,Nx,Ny);     
    
        %%% Make the plot
        pcolor(XX_h/1000,YY_h/1000,hh);
        shading interp;
        colorbar;
        title(strcat(['t=',num2str(t/t1day,'%.2f'),' days']));
        colormap jet; 
%         caxis([0 1000])
        
      %%% Plot tracer
      case 'b'
        
        %%% Load h
        data_file = fullfile(dirpath,[OUTN_B,num2str(layer-1),'_n=',num2str(n),'.dat']);
        bb = readOutputFile(data_file,Nx,Ny);   
    
        %%% Make the plot
        pcolor(XX_h/1000,YY_h/1000,bb);
        shading interp;
        colorbar;
        title(strcat(['t=',num2str(t/t1day,'%.2f'),' days']));
        colormap jet;
       
      %%% Plot surface pressure
      case 'p'
        
        %%% Calculate layer surface height
        eta = zeros(Nx,Ny,Nlay+1);
        eta(:,:,Nlay+1) = hhb;        
        for k=Nlay:-1:1

          %%% Load kth layer thickness
          data_file = fullfile(dirpath,[OUTN_H,num2str(k-1),'_n=',num2str(n),'.dat']);
          hh = readOutputFile(data_file,Nx,Ny);    
          
          %%% Add layer thickness to eta
          eta(:,:,k) = eta(:,:,k+1) + hh;
          
        end
        
        %%% Calculate uppermost layer presure
        if (useRL)        
          %%% Load surface pressure
          data_file = fullfile(dirpath,[OUTN_PI,'_n=',num2str(n),'.dat']);
          pi = readOutputFile(data_file,Nx,Ny); 
        else
          pi = gg(1)*eta(:,:,1);
        end
        
        %%% Calculate hydrostatic pressure
        for k=2:layer
          pi = pi + gg(k)*eta(:,:,k);          
        end                
    
        %%% Make the plot
        pcolor(XX_h/1000,YY_h/1000,pi);
        hold on
        [C,h] = contour(XX_h/1000,YY_h/1000,hhb,[-4000:500:-1000 -750],'EdgeColor','k');        
        caxis([min(min(pi)) max(max(pi))]);
        hold off;
        shading interp;
        colorbar;
        title(strcat(['t=',num2str(t/t1day,'%.2f'),' days']));
        colormap jet;        
%         caxis([-200 200]);
       

      %%% Plot zonal velocity
      case 'u'
        
        %%% Load u
        data_file = fullfile(dirpath,[OUTN_U,num2str(layer-1),'_n=',num2str(n),'.dat']);
        uu = readOutputFile(data_file,Nx,Ny);       
    
        %%% Make the plot
        pcolor(XX_u/1000,YY_u/1000,uu);        
        shading interp;
        hold on
        [C,h] = contour(XX_h/1000,YY_h/1000,hhb,[-4000:500:-1000 -750],'EdgeColor','k');
        clabel(C,h);
        hold off;
        colorbar;
        title(strcat(['t=',num2str(t/t1day,'%.2f'),' days']));
        colormap jet;        
        caxis([-max(max(abs(uu))) max(max(abs(uu)))]);
        
      %%% Plot meridional velocity
      case 'v'
        
        %%% Load v
        data_file = fullfile(dirpath,[OUTN_V,num2str(layer-1),'_n=',num2str(n),'.dat']);
        vv = readOutputFile(data_file,Nx,Ny);      
    
        %%% Make the plot
        pcolor(XX_v/1000,YY_v/1000,vv);
        shading interp;
        colorbar;
        title(strcat(['t=',num2str(t/t1day,'%.2f'),' days']));
        colormap jet;        
        caxis([-max(max(abs(vv))) max(max(abs(vv)))]);
        
        
      case 'vcorr'
        
        %%% Load v
        uu = zeros(Nx,Ny,Nlay);
        vv = zeros(Nx,Ny,Nlay);
        for k=1:Nlay
          data_file = fullfile(dirpath,[OUTN_U,num2str(k-1),'_n=',num2str(n),'.dat']);
          uu(:,:,k) = readOutputFile(data_file,Nx,Ny);      
          data_file = fullfile(dirpath,[OUTN_V,num2str(k-1),'_n=',num2str(n),'.dat']);
          vv(:,:,k) = readOutputFile(data_file,Nx,Ny);      
        end
        uabs = sqrt(uu.^2+vv.^2);
        usurf = vv(:,:,1);
        usub = vv(:,:,3);
        scatter(usurf(:),usub(:))
        xlabel('Surface velocity');
        ylabel('Subsurface velocity');
        
        r = corr(usurf(:),usub(:));
        disp(r^2);

      %%% Plot zonally-averaged zonal velocity
      case 'ua'
        
        %%% Load u
        uu = zeros(Nx,Ny,Nlay);
        for k=1:Nlay
          data_file = fullfile(dirpath,[OUTN_U,num2str(k-1),'_n=',num2str(n),'.dat']);
          uu(:,:,k) = readOutputFile(data_file,Nx,Ny);    
        end        
        
        %%% Make the plot
        plot(yy_u/1000,mean(uu(:,:,1),1));        
        hold on
        for k=2:Nlay
          plot(yy_u/1000,mean(uu(:,:,k),1));        
        end
        plot(yy_u/1000,0*yy_u,'k--');
        hold off
        title(strcat(['t=',num2str(t/t1day,'%.2f'),' days']));        
        axis([0 Ly/1000 -.3 .3]);
        
        
      %%% Plot zonally-averaged meridional velocity
      case 'va'
        
        %%% Load v
        data_file = fullfile(dirpath,[OUTN_V,num2str(layer-1),'_n=',num2str(n),'.dat']);
        vv = readOutputFile(data_file,Nx,Ny);    
            
        %%% Make the plot
        plot(yy_v/1000,mean(vv,1));        
        hold on
        plot(yy_v/1000,0*yy_v,'k--');
        hold off
        title(strcat(['t=',num2str(t/t1day,'%.2f'),' days']));        
        axis([0 Ly/1000 -.02 .02]);
    
      %%% Zonally-averaged isopycnals
      case 'i'

        %%% Calculate layer surface height
        eta = zeros(Nx,Ny,Nlay+1);
        eta(:,:,Nlay+1) = hhb;
        for k=Nlay:-1:1

          %%% Load kth layer thickness
          data_file = fullfile(dirpath,[OUTN_H,num2str(k-1),'_n=',num2str(n),'.dat']);
          hh = readOutputFile(data_file,Nx,Ny);    
          
          %%% Add layer thickness to eta
          eta(:,:,k) = eta(:,:,k+1) + hh;
          
        end                
    
        %%% Make the plot
        plot(yy_h, squeeze(mean(eta,1)));                
        title(strcat(['t=',num2str(t/t1day,'%.2f'),' days']));        
        xlabel('x');
        ylabel('z');        
        
      %%% 3D isopycnals
      case 'i3d'

        %%% Calculate layer surface height
        eta = zeros(Nx,Ny,Nlay+1);
        eta(:,:,Nlay+1) = hhb;
        hh = zeros(Nx,Ny,Nlay);
        zeta = zeros(Nx+1,Ny+1,Nlay);            
        for k=Nlay:-1:1

          %%% Load kth layer thickness
          data_file = fullfile(dirpath,[OUTN_H,num2str(k-1),'_n=',num2str(n),'.dat']);
          hh(:,:,k) = readOutputFile(data_file,Nx,Ny);    
          
          %%% Add layer thickness to eta
          eta(:,:,k) = eta(:,:,k+1) + hh(:,:,k);
          
          %%% Calculate the relative vorticity
          data_file = fullfile(dirpath,[OUTN_U,num2str(layer-1),'_n=',num2str(n),'.dat']);
          uu = readOutputFile(data_file,Nx,Ny);
          data_file = fullfile(dirpath,[OUTN_V,num2str(layer-1),'_n=',num2str(n),'.dat']);
          vv = readOutputFile(data_file,Nx,Ny);           
          zeta(2:Nx,2:Ny,k) = (uu(1:Nx-1,1:Ny-1)-uu(1:Nx-1,2:Ny)) / dy + (vv(2:Nx,2:Ny)-vv(1:Nx-1,2:Ny)) / dx;                      
          
        end
        
        data_file = fullfile(dirpath,[OUTN_PI,'_n=',num2str(n),'.dat']);
        pi = readOutputFile(data_file,Nx,Ny); 

        %%% Correct boundaries
        if (~useWallNS)
          zeta(:,Ny+1,:) = zeta(:,1,:);
        end
        zeta(Nx+1,:,:) = zeta(1,:,:);        
        Ro = zeta(:,:,1)./(2*Omega_z);

        %%% Calculate PV        
        ff = 2*Omega_z;        
        hh_q = NaN*ones(Nx+1,Ny+1,Nlay);
        hh_q(2:Nx,2:Ny,:) = 0.25 * (hh(2:Nx,1:Ny-1,:)+hh(2:Nx,2:Ny,:)+hh(1:Nx-1,1:Ny-1,:)+hh(1:Nx-1,2:Ny,:));        
        pv = (ff+zeta)./hh_q; 
        min(min(pv))
        max(max(pv))
        
        %%% Make the plot
        clf;
        pbaspect([2,1,1]);
        eta_plot = eta(:,:,1);        
%         p = surface(XX_h/1000,YY_h/1000,eta_plot,Ro);
%         p = surface(XX_h/1000,YY_h/1000,eta_plot,sqrt(uu.^2+vv.^2));
        p = surface(XX_h/1000,YY_h/1000,eta_plot,pv(:,:,1));
%         p = surface(XX_h/1000,YY_h/1000,eta_plot,pi/gg(1));
        p.FaceColor = 'texturemap';
%         contourf(XX_h/1000,YY_h/1000,eta_plot,pi/gg(1),[-.5:0.05:.5]);
%         colormap(cmocean('curl'));
%         colormap(cmocean('balance'));
%         colormap(cmocean('amp'));
        colormap(pmkmp(100,'Swtth'));
%         colormap haxby;
%         caxis([-.3 .3]);
%         caxis([0 .5]);
        caxis([-1e-7 -4e-8]);
%         caxis([-.3 .3]);
        p.EdgeColor = 'none';         
        alpha(p,1);
        hold on;
        eta_plot = eta(:,:,2);
        p = surface(XX_h/1000,YY_h/1000,eta_plot);
        p.FaceColor = [48 129 238]/256;
        p.EdgeColor = 'none';        
        alpha(p,0.7);
%         eta_plot = eta(:,:,3);
%         p = surface(XX_h/1000,YY_h/1000,eta_plot);
%         p.FaceColor = [24 60 139]/256;
%         p.EdgeColor = 'none';        
%         alpha(p,0.7);
        p = surface(XX_h/1000,YY_h/1000,eta(:,:,Nlay+1));
        p.FaceColor = [139,69,19]/256;
        p.EdgeColor = 'none';          
        hold off;        
        view(60,32);
%         view(32,44);
        lighting gouraud;
        camlight('headlight');        
        title(strcat(['t=',num2str(t/t1year,'%.2f'),' years']));        
        xlabel('x (km)');
        ylabel('y (km)');
        zlabel('z (m)');
        set(gca,'FontSize',16);
        set(gcf,'Color','w');        
        
    end
    
             
    %%% Store the image in the movie buffer
    nextframe = getframe(gcf);        
    M(Mcnt) = nextframe;       
    Mcnt = Mcnt + 1;
    
  end
  
end

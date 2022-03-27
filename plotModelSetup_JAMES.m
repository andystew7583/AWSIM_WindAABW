%%%
%%% plotModelSetup_JPO.m
%%%
%%% Plots wind forcing, instantaneous model state and bathymetry for our JAMES paper.
%%%

%%% Run to load
run_name = 'ACC_AABW_ML_doubleMOC_hires';

%%% Load parameters   
local_home_dir = '/Volumes/Kilchoman/UCLA/Projects/AWSIM/runs';
loadParams;
dirpath = fullfile(local_home_dir,run_name);
gtild = reshape(cumsum(gg),[Nlay 1 1]);

%%% Output file number to load
layer = 1;
n = 3750;

%%% Physical parameters
rho0 = 1000;
Ln = 100*m1km;

%%% Load forcing data
dt_tau = tauPeriod/tauNrecs;
tt_tau = 0:dt_tau:(tauNrecs-1)*dt_tau;
[wDiaPeriod wDiaPeriod_found] = readparam(params_file,'wDiaPeriod','%lf');
if (~wDiaPeriod_found)
  wDiaPeriod = 0;
end
[wDiaNrecs wDiaNrecs_found] = readparam(params_file,'wDiaNrecs','%lf');
if (~wDiaNrecs_found)
  wDiaNrecs = 1;
end
fid = fopen(fullfile(dirpath,'wDiaFile.dat'),'r','b'); %%% Read diapycnal velocity data
if (fid == -1)
  error(['Could not open wDiaFile']);
end
wDia = fread(fid,wDiaNrecs*(Nlay+1)*Nx*Ny,'real*8','ieee-le');
fclose(fid);
if (length(wDia) ~= Nx*Ny*wDiaNrecs*(Nlay+1)) %% Check all data was read
  error(['Insufficient data found in ',paramFile]);
end
wDia = reshape(wDia,[wDiaNrecs Nlay+1 Nx Ny]); %%% Reshape to required 3D shape
dt_wDia = wDiaPeriod/wDiaNrecs;
tt_wDia = 0:dt_wDia:(wDiaNrecs-1)*dt_wDia;

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



%%% Plotting options
scrsz = get(0,'ScreenSize');
framepos = [209   632  900   900];
fontsize = 14;
ax_pos = zeros(3,4);
ax_pos(1,:) = [0.06 0.54 0.9 0.46];
ax_pos(2,:) = [0.07 0.36 0.39 0.1];
ax_pos(3,:) = [0.07 0.06 0.39 0.26];
ax_pos(4,:) = [0.56 0.36 0.4 0.1];
ax_pos(5,:) = [0.56 0.21 0.4 0.1];
ax_pos(6,:) = [0.56 0.06 0.4 0.1];
cb_pos = [0.92 0.55 0.015 0.4];
lab_size = [0.05 0.03];


%%% Set up the frame
figure(201);
clf; 
set(gcf,'Position',framepos);
set(gcf,'Color','w');
cntr = 0;



%%% Plot wind stress
cntr = cntr + 1;
ax = axes('position',ax_pos(cntr,:));  
set(gca,'ZDir','reverse');
pbaspect([2,1,1]);
eta_plot = eta(:,:,1);        
%         p = surface(XX_h/1000,YY_h/1000,eta_plot,Ro);
%         p = surface(XX_h/1000,YY_h/1000,eta_plot,sqrt(uu.^2+vv.^2));
p = surface(XX_h/1000,YY_h/1000,-eta_plot,pv(:,:,1));
%         p = surface(XX_h/1000,YY_h/1000,eta_plot,pi/gg(1));
p.FaceColor = 'texturemap';
%         contourf(XX_h/1000,YY_h/1000,eta_plot,pi/gg(1),[-.5:0.05:.5]);
colormap(cmocean('curl'));
%         colormap(cmocean('amp'));
%         colormap(pmkmp(100,'Swtth'));
%         colormap haxby;
%         caxis([-.3 .3]);
%         caxis([0 .5]);
caxis([-2e-7 -7e-8]);
%         caxis([-.5 .5]);
p.EdgeColor = 'none';         
alpha(p,1);
hold on;
eta_plot = eta(:,:,2);
p = surface(XX_h/1000,YY_h/1000,-eta_plot);
p.FaceColor = [48 129 238]/256;
p.EdgeColor = 'none';        
alpha(p,0.7);
eta_plot = eta(:,:,3);
p = surface(XX_h/1000,YY_h/1000,-eta_plot);
p.FaceColor = [24 60 139]/256;
p.EdgeColor = 'none';        
alpha(p,0.7);
p = surface(XX_h/1000,YY_h/1000,-eta(:,:,Nlay+1));
p.FaceColor = [139,69,19]/256;
p.EdgeColor = 'none';          
hold off;        
view(30,32);
%         view(32,44);
lighting gouraud;
camlight('headlight');        
title(strcat(['t=',num2str(t/t1year,'%.2f'),' years']));        
xlabel('Longitude x (km)');
ylabel('Latitude y (km)');
zlabel('Depth z (m)');
set(gca,'FontSize',fontsize);
set(gcf,'Color','w');    
cbhandle = colorbar;
set(cbhandle,'Position',cb_pos);
grid on;
axis([0 3200 0 1600 0 4000]);
annotation(figure1,'textbox',...
  [0.906349206349207 0.954678719620628 0.0407142857142857 0.0328571428571428],...
  'String',{'1/ms'},...
  'LineStyle','none',...
  'FontSize',14,...
  'FitBoxToText','off');




%%% Plot wind profile
cntr = cntr + 1;
ax = axes('position',ax_pos(cntr,:));  
plot(yy_h/m1km,squeeze(mean(taux(:,1,:),1))*rho0);
ylabel('Wind stress (N/m^2)');
set(gca,'FontSize',fontsize);

%%% Plot bathymetry and pycnocline depth
cntr = cntr + 1;
ax = axes('position',ax_pos(cntr,:));  
plot(yy_h/m1km,-hhb(1,:),'Color',[0.41,0.19,0.07],'LineWidth',1.5);
hold on;
plot(yy_h/m1km,min(-hhb,[],1),'Color',[0.41,0.19,0.07],'LineWidth',1.5,'LineStyle','--');
plot(yy_h/m1km,-mean(eta(:,:,2),1),'-','LineWidth',1,'Color',[48 129 238]/256);
plot(yy_h/m1km,-mean(eta(:,:,3),1),'-','LineWidth',1,'Color',[24 60 139]/256);
plot([Ln Ln]/m1km,[0 4000],'k--');
fill([Ly-Ln Ly Ly Ly-Ln Ly-Ln]/m1km,[0 0 4000 4000 0],ones(1,5),'FaceColor',[.5 .5 .5],'EdgeColor','None','FaceAlpha',0.5);
hold off;
axis([0 Ly/m1km 0 4000]);
set(gca,'YDir','reverse');
xlabel('Latitude y (km)');
ylabel('Depth (m)');
set(gca,'FontSize',fontsize);
text(1550,1500,'Sponge','Rotation',270,'FontSize',14);

%%% Plot wind stress time series
cntr = cntr + 1;
ax = axes('position',ax_pos(cntr,:));  
plot(tt_tau/t1year,max(taux(:,1,:),[],3)*rho0);
ylabel('\tau_m_a_x (N/m^2)');
set(gca,'FontSize',fontsize);

%%% Plot AAIW formation time series
cntr = cntr + 1;
ax = axes('position',ax_pos(cntr,:));  
plot(tt_wDia/t1year,sum(sum(wDia(:,2,:,:),3),4)*dx*dy/1e6);
ylabel('\varpi_A_A_I_W (Sv)');
set(gca,'FontSize',fontsize);

%%% Plot AABW formation time series
cntr = cntr + 1;
ax = axes('position',ax_pos(cntr,:));  
plot(tt_wDia/t1year,-sum(sum(wDia(:,3,:,:),3),4)*dx*dy/1e6);
ylabel('\varpi_A_A_B_W (Sv)');
set(gca,'FontSize',fontsize);
xlabel('Time (years)');

%%% Add annotations
figure1 = gcf;


% Create arrow
annotation(figure1,'arrow',[0.0822222222222222 0.0822222222222222],...
  [0.195666666666667 0.162222222222222],'Color','r','LineWidth',4,'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.0833333333333333 0.0833333333333333],...
  [0.250111111111111 0.29],'Color','r','LineWidth',4,'HeadStyle','plain');

% Create doublearrow
annotation(figure1,'doublearrow',[0.276666666666667 0.276666666666667],...
  [0.126666666666667 0.0581005071461503]);

% Create textbox
annotation(figure1,'textbox',...
  [0.28031746031746 0.0721454916683131 0.0471428571428573 0.0375714285714286],...
  'String',{'H_r_i_d_g_e'},...
  'LineStyle','none',...
  'FontSize',14,...
  'FitBoxToText','off');

% Create doublearrow
annotation(figure1,'doublearrow',[0.262222222222222 0.262222222222222],...
  [0.358888888888889 0.458888888888889]);

% Create textbox
annotation(figure1,'textbox',...
  [0.0729206349206347 0.123137192913127 0.0592857142857143 0.0385714285714287],...
  'String',{'\varpi_A_A_B_W'},...
  'LineStyle','none',...
  'FontSize',14,...
  'FitBoxToText','off');

% Create textbox
annotation(figure1,'textbox',...
  [0.271809523809524 0.38647052624646 0.0592857142857143 0.0385714285714286],...
  'String',{'\tau_m_a_x'},...
  'LineStyle','none',...
  'FontSize',14,...
  'FitBoxToText','off');

% Create textbox
annotation(figure1,'textbox',...
  [0.0729206349206347 0.285359415135349 0.0592857142857142 0.0385714285714285],...
  'String',{'\varpi_A_A_I_W'},...
  'LineStyle','none',...
  'FontSize',14,...
  'FitBoxToText','off');



% Create textbox
annotation(figure1,'textbox',...
  [0.0144444444444445 0.281111111111111 0.05 0.03],'String','(c)',...
  'LineStyle','none',...
  'Interpreter','latex',...
  'FontSize',16,...
  'FitBoxToText','off');

% Create textbox
annotation(figure1,'textbox',...
  [0.0122222222222222 0.488888888888889 0.05 0.03],'String','(b)',...
  'LineStyle','none',...
  'Interpreter','latex',...
  'FontSize',16,...
  'FitBoxToText','off');

% Create textbox
annotation(figure1,'textbox',...
  [0.503333333333333 0.466666666666667 0.0500000000000002 0.03],...
  'String','(d)',...
  'LineStyle','none',...
  'Interpreter','latex',...
  'FontSize',16,...
  'FitBoxToText','off');

% Create textbox
annotation(figure1,'textbox',...
  [0.503333333333333 0.307777777777778 0.0500000000000002 0.03],...
  'String','(e)',...
  'LineStyle','none',...
  'Interpreter','latex',...
  'FontSize',16,...
  'FitBoxToText','off');

% Create textbox
annotation(figure1,'textbox',...
  [0.501111111111111 0.161111111111111 0.05 0.03],'String','(f)',...
  'LineStyle','none',...
  'Interpreter','latex',...
  'FontSize',16,...
  'FitBoxToText','off');

% Create textbox
annotation(figure1,'textbox',...
  [0.131111111111111 0.956666666666667 0.0500000000000001 0.03],...
  'String','(a)',...
  'LineStyle','none',...
  'Interpreter','latex',...
  'FontSize',16,...
  'FitBoxToText','off');





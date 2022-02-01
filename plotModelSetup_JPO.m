%%%
%%% plotModelSetup_JPO.m
%%%
%%% Plots wind forcing, instantaneous model state and bathymetry for our JPO paper.
%%%

%%% Run to load
run_name = 'ACC_AABW_Ny256_Nlay2_tauM0.1_tauP0_tauF0_wDiaM0_wDiaP0_wDiaF0_Cd2.000e-03_rb0.000e+00_diags';

%%% Load parameters   
local_home_dir = '/Volumes/Kilchoman/UCLA/Projects/AWSIM_WindAABW/runs';
prod_dir = fullfile(local_home_dir,'AWSIM_WindAABW_products');
loadParams;
dirpath = fullfile(local_home_dir,run_name);
gtild = reshape(cumsum(gg),[Nlay 1 1]);

%%% Averaging period
tmin = 0.5*t1year;
tmax = 30.5*t1year;

%%% Output file number to load
layer = 1;
n = 22;

%%% Physical parameters
rho0 = 1000;
Ln = 100*m1km;

%%% Remove redundant time dimension in wind stress vector
taux = squeeze(mean(taux,1));
tauy = squeeze(mean(tauy,1));

%%% Load u
data_file = fullfile(dirpath,[OUTN_U,num2str(layer-1),'_n=',num2str(n),'.dat']);
uu = readOutputFile(data_file,Nx,Ny);

%%% Load v
data_file = fullfile(dirpath,[OUTN_V,num2str(layer-1),'_n=',num2str(n),'.dat']);
vv = readOutputFile(data_file,Nx,Ny);    

%%% Load h
data_file = fullfile(dirpath,[OUTN_H,num2str(layer-1),'_n=',num2str(n),'.dat']);
hh = readOutputFile(data_file,Nx,Ny);      

%%% Load h
data_file = fullfile(dirpath,[OUTN_PI,'_n=',num2str(n),'.dat']);
pp = readOutputFile(data_file,Nx,Ny);      

%%% Calculate the relative vorticity
zeta = zeros(Nx+1,Ny+1);            
zeta(2:Nx,2:Ny) = (uu(1:Nx-1,1:Ny-1)-uu(1:Nx-1,2:Ny)) / dy + (vv(2:Nx,2:Ny)-vv(1:Nx-1,2:Ny)) / dx;                                         

%%% Calculate PV        
ff = 2*Omega_z;        
hh_q = NaN*ones(Nx+1,Ny+1);
hh_q(2:Nx,2:Ny) = 0.25 * (hh(2:Nx,1:Ny-1)+hh(2:Nx,2:Ny)+hh(1:Nx-1,1:Ny-1)+hh(1:Nx-1,2:Ny));        
pv = (ff+zeta)./hh_q;                   
  



%%% Plotting options
scrsz = get(0,'ScreenSize');
framepos = [209   632   900   482];
fontsize = 14;
ax_pos = zeros(3,4);
ax_pos(1,:) = [0.07 0.45 0.1 0.5];
ax_pos(2,:) = [0.26 0.45 0.65 0.5];
ax_pos(3,:) = [0.26 0.09 0.65 0.28];
cb_pos = [0.94 0.45 0.015 0.5];
lab_size = [0.05 0.03];


%%% Set up the frame
figure(101);
clf; 
set(gcf,'Position',framepos);
set(gcf,'Color','w');
cntr = 0;



%%% Plot wind stress
cntr = cntr + 1;
ax = axes('position',ax_pos(cntr,:));  
plot(rho0*taux(1,:),yy_h/m1km,'LineWidth',1.5);
hold on;
plot([0 0],[0 Ly/m1km],'k--','LineWidth',0.5);
hold off;
axis([0 0.1 0 Ly/m1km]);
ylabel('Latitude y (km)');
xlabel('\tau^(^x^) (N/m^2)');
% title('Wind stress');
set(gca,'FontSize',fontsize);
annotation('textbox',[ax_pos(cntr,1)-0.06 ax_pos(cntr,2)-0.06 lab_size],'String','(a)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

%%% Plot PV snapshot
cntr = cntr + 1;
ax = axes('position',ax_pos(cntr,:));  
pcolor(XX_q/1000,YY_q/1000,pv);
shading interp;
hold on;
[C,h]=contour(XX_h/1000,YY_h/1000,-hhb,[500 1500 2500 3250 3500 3750],'EdgeColor',[.3 .3 .3],'LineWidth',1);
fill([0 Lx Lx 0 0]/m1km,[Ly-Ln Ly-Ln Ly Ly Ly-Ln]/m1km,ones(1,5),'FaceColor',[.5 .5 .5],'EdgeColor','None','FaceAlpha',0.5);
% plot([0 Lx Lx 0 0]/m1km,[Ly-Ln Ly-Ln Ly Ly Ly-Ln]/m1km,'k-');
% clabel(C,h);
% plot([-Lx/2/m1km Lx/2/m1km Lx/2/m1km -Lx/2/m1km -Lx/2/m1km],[-Ly/2/m1km -Ly/2/m1km Ly/2/m1km Ly/2/m1km -Ly/2/m1km],'k-');
hold off;
% axis([-Lx/2/m1km Lx/2/m1km -Ly/2/m1km Ly/2/m1km]);
set(gca,'FontSize',fontsize);
colormap(pmkmp(100,'Swtth'));
cbhandle = colorbar;
set(cbhandle,'Position',cb_pos);
caxis([min(pv(:)) max(pv(:))]);
% title('Model bathymetry');
annotation('textbox',[ax_pos(cntr,1)-0.06 ax_pos(cntr,2)-0.06 lab_size],'String','(b)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
text(Lx/2/m1km,(Ly-Ln/2)/m1km,'Sponge','FontSize',10);

%%% Plot bathymetry and pycnocline depth
cntr = cntr + 1;
ax = axes('position',ax_pos(cntr,:));  
plot(xx_h/1000,-hhb(:,Ny/2),'Color',[0.41,0.19,0.07],'LineWidth',1.5);
hold on;
plot(xx_h/1000,mean(hh(:,:,1),2),'b-','LineWidth',1);
% plot(xx_h/1000,mean(pp(:,:),2)/gg(1)*1000,'b-','LineWidth',1);
hold off;
% set(gca,'XLim',[-Lx/2/m1km Lx/2/m1km]);
axis([0 Lx/m1km 0 4000]);
set(gca,'YDir','reverse');
xlabel('Longitude x (km)');
ylabel('Depth (m)');
set(gca,'FontSize',fontsize);
annotation('textbox',[ax_pos(cntr,1)-0.06 ax_pos(cntr,2)-0.06 lab_size],'String','(c)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');


%%% Add annotations
figure1 = gcf;

annotation(figure1,'line',[0.86 0.867714285714286],...
  [0.351857142857143 0.351857142857143]);

% Create doublearrow
annotation(figure1,'doublearrow',[0.0719047619047623 0.171111111111111],...
  [0.695068168346176 0.695020746887967]);

% Create textbox
annotation(figure1,'textbox',...
  [0.0995873015873014 0.712026081802016 0.0592857142857143 0.0385714285714286],...
  'String',{'\tau_m_a_x'},...
  'LineStyle','none',...
  'FontSize',14,...
  'FitBoxToText','off');

% Create textbox
annotation(figure1,'textbox',...
  [0.535873015873016 0.121034380557202 0.0471428571428573 0.0375714285714286],...
  'String',{'H_r_i_d_g_e'},...
  'LineStyle','none',...
  'FontSize',14,...
  'FitBoxToText','off');

% Create doublearrow
annotation(figure1,'doublearrow',[0.528888888888889 0.528888888888889],...
  [0.170124481327801 0.0892116182572614]);

% Create doublearrow
annotation(figure1,'doublearrow',[0.422380952380953 0.505555555555555],...
  [0.209136336692353 0.20954356846473]);

% Create textbox
annotation(figure1,'textbox',...
  [0.440158730158731 0.233997036158862 0.0721428571428571 0.0385714285714285],...
  'String',{'W_r_i_d_g_e'},...
  'LineStyle','none',...
  'FontSize',14,...
  'FitBoxToText','off');

% Create textbox
annotation(figure1,'textbox',...
  [0.926349206349207 0.964678719620628 0.0407142857142857 0.0328571428571428],...
  'String',{'1/ms'},...
  'LineStyle','none',...
  'FontSize',14,...
  'FitBoxToText','off');



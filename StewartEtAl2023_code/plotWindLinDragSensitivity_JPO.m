%%%
%%% plotWindLinDragSensitivity_JPO.m
%%%
%%% Plots sensitivity of ACC transport and eddy properties to wind and drag
%%% for our JPO paper. This version of the code uses our linear drag
%%% experiments.
%%%

%%% Load static definitions
constants;

%%% Directory to store runs
local_home_dir = '/Volumes/Kilchoman/UCLA/Projects/AWSIM_WindAABW/runs';
prod_dir = fullfile('./products');

%%% Spinup simulations are long and produce no diagnostic, diagnostic
%%% simulations output high-frequency diagnostics to resolve the forcing
%%% period
is_spinup = false;

%%% Grid resolution 
Ny = 256;
Nlay = 2;

%%% Averaging period
tmin = 0.5*t1year;
tmax = 30.5*t1year;

%%% Index of reference simulation
idx_ref_tm = 9;
idx_ref_rb = 3;

%%% Parameters defining the batch of runs to plot
% tau_mean = [0.01 0.017 0.03 0.05 0.1 0.17 0.3];
tau_mean = [0.01 0.013 0.017 0.022 0.03 0.039 0.05 0.07 0.1 0.13 0.17 0.22 0.3 0.39 0.5];
tau_pert = 0;
tau_freq = 0;
% AABW_mean = [-1.5 -.75 0 .75 1.5];
AABW_mean = 0;
% AABW_mean = [-1.5 0 1.5];
AABW_pert = 0;
AABW_freq = 0;
% quad_drag = 2e-3;
% quad_drag = [.5e-3 1e-3 1.5e-3 2e-3 2.5e-3 3e-3 3.5e-3 4e-3];
% lin_drag = 0e-4;
lin_drag = [2e-4 3e-4 4e-4 5e-4 6e-4 7e-4 8e-4 9e-4 10e-4];
quad_drag = 0e-3;
% lin_drag = 2e-4;
topog_width = 150;
topog_height = 1000;
N_tm = length(tau_mean);
N_rb = length(lin_drag);

%%% Loop over runs and compute transport
Ttot = zeros(N_rb,N_tm);
Tbt = zeros(N_rb,N_tm);
Tbc = zeros(N_rb,N_tm);
kap = zeros(N_rb,N_tm);
r_kap = zeros(N_rb,N_tm);
SKE_tot = zeros(N_rb,N_tm);
EKE_tot = zeros(N_rb,N_tm);
for n_tm = 1:N_tm
  for n_rb = 1:N_rb
    
    [n_tm,n_rb]
    %%% Simulation name
    run_name = constructRunName (is_spinup,Ny,Nlay, ...
                            tau_mean(n_tm),tau_pert,tau_freq, ...
                            AABW_mean,AABW_pert,AABW_freq, ...
                            quad_drag,lin_drag(n_rb),topog_width,topog_height);
    loadParams;
    
    %%% Read time-mean zonal flux    
    hu_tavg = do_avg(dirpath,OUTN_HU_AVG,Nx,Ny,Nlay,n0_avg,N_avg,dt_avg,tmin,tmax,startTime);
    hv_tavg = do_avg(dirpath,OUTN_HV_AVG,Nx,Ny,Nlay,n0_avg,N_avg,dt_avg,tmin,tmax,startTime);
    hh_tavg = do_avg(dirpath,OUTN_H_AVG,Nx,Ny,Nlay,n0_avg,N_avg,dt_avg,tmin,tmax,startTime);
    uu_tavg = do_avg(dirpath,OUTN_U_AVG,Nx,Ny,Nlay,n0_avg,N_avg,dt_avg,tmin,tmax,startTime);
    vv_tavg = do_avg(dirpath,OUTN_V_AVG,Nx,Ny,Nlay,n0_avg,N_avg,dt_avg,tmin,tmax,startTime);
    husq_tavg = do_avg(dirpath,OUTN_HUU_AVG,Nx,Ny,Nlay,n0_avg,N_avg,dt_avg,tmin,tmax,startTime);
    hvsq_tavg = do_avg(dirpath,OUTN_HVV_AVG,Nx,Ny,Nlay,n0_avg,N_avg,dt_avg,tmin,tmax,startTime);
  
    %%% Compute mean layer thicknesses on cell edges
    hh_w_tavg = 0.5*(hh_tavg(1:Nx,:,:)+hh_tavg([Nx 1:Nx-1],:,:)); %%% N.B. assumes double periodicity
    hh_s_tavg = 0.5*(hh_tavg(:,1:Ny,:)+hh_tavg(:,[Ny 1:Ny-1],:)); 
  
    %%% Thickness-weighted average velocities
    uu_twa_t = hu_tavg ./ hh_w_tavg;
    vv_twa_t = hv_tavg ./ hh_s_tavg;
    
    uu_twa_x = mean(hu_tavg,1) ./ mean(hh_w_tavg,1);
    vv_twa_x = mean(hv_tavg,1) ./ mean(hh_s_tavg,1);
    uu_sw = uu_twa_t - repmat(uu_twa_x,[Nx 1]);
    vv_sw = vv_twa_t - repmat(vv_twa_x,[Nx 1]);
    SKE_u = 0.5*sum(hh_w_tavg.*uu_sw.^2,1);
    SKE_v = 0.5*sum(hh_s_tavg.*vv_sw.^2,1);
    SKE_tot(n_rb,n_tm) = (sum(sum(SKE_u*dx*dy)) + sum(sum(SKE_v*dx*dy))) / sum(sum((-hhb)*dx*dy));
  
    %%% Compute transports
    Ttot(n_rb,n_tm) = mean(sum(sum(hu_tavg,3),2),1)*dy;
    Tbt(n_rb,n_tm) = mean(sum(uu_tavg(:,:,end).*(-hhb).*dy,2),1);
    Tbc(n_rb,n_tm) = Ttot(n_rb,n_tm) - Tbt(n_rb,n_tm);
    
    %%% Calculate at the ridge - yields almost identical results
%     hhb_u = 0.5*(hhb(1:Nx,:)+hhb([Nx 1:Nx-1],:)); 
%     idx_topog = find(hhb_u(:,Ny/2)==max(hhb_u(:,Ny/2)));
%     Ttot(n_rb,n_tm) = sum(sum(squeeze(hu_tavg(idx_topog,:,:)),2),1)*dy;
%     Tbt(n_rb,n_tm) = sum(uu_tavg(idx_topog,:,end).*(-hhb(idx_topog,:)).*dy,2);
%     Tbc(n_rb,n_tm) = Ttot(n_rb,n_tm) - Tbt(n_rb,n_tm);
    
%     %%% Compute diffusivity
%     [kap_bulk,nu_bulk,r_kap_bulk,r_nu_bulk,EKE_zavg] = calcBulkEddyViscDiff(local_home_dir,run_name);
%     kap(n_rb,n_tm) = kap_bulk;
%     r_kap(n_rb,n_tm) = r_kap_bulk;

    %%% Load eddy diffusivity and viscosity maps
    load(fullfile(prod_dir,['kap_nu_',run_name,'.mat']));
    EKE_thresh = quantile(EKE_zavg(:),0.75);
    kap(n_rb,n_tm) = mean(kap_map(EKE_zavg>EKE_thresh));  
    EKE_tot(n_rb,n_tm) = sum(sum(EKE_zavg.*(-hhb)*dx*dy)) / sum(sum((-hhb)*dx*dy));
    
  end
end

 

%%% Plotting options
fontsize = 14;
axpos = zeros(4,4);
axpos(1,:) = [0.12 0.71 .84 .28];
axpos(2,:) = [0.12 0.38 .84 .28];
axpos(3,:) = [0.12 0.05 .84 .28];
axlabels = {'(a)','(b)','(c)'};
lab_size = [0.05 0.03];
tau_ticks = [0.01 0.017 0.03 0.05 0.1 0.17 0.3 0.5];
% tau_ticks = [0:0.05:0.5];
rho0 = 1000;

defaultcolororder = zeros(8,3);
defaultcolororder(1:7,:) = get(gca,'ColorOrder');
defaultcolororder(8,:) = [.3 .3 .3];
colororder = zeros(9,3);
colororder(1:idx_ref_rb-1,:) = defaultcolororder(1:idx_ref_rb-1,:);
colororder(idx_ref_rb,:) = [0 0 0];
colororder(idx_ref_rb+1:N_rb,:) = defaultcolororder(idx_ref_rb:N_rb-1,:);

markersize = 10;
markershapes = {'>','o','*','<','v','d','^','s','x','+'};

figure(107);
clf;
set(gcf,'Position',[382   306   500   1000]);



%%% Make figure
subplot('Position',axpos(1,:));
for n_rb=1:N_rb
  semilogx(tau_mean,Ttot(n_rb,:)/1e6,'o-','Color',colororder(n_rb,:),'Marker',markershapes{n_rb},'MarkerFaceColor',colororder(n_rb,:));
  if (n_rb == 1)
  hold on;
  end
end
hold off;
% xlabel('Wind Stress (N/m^2)');
ylabel('Total transport (Sv)');
set(gca,'XTick',tau_ticks);
set(gca,'XLim',[0 0.5]);
set(gca,'YLim',[-10 100]);
set(gca,'FontSize',fontsize);
grid on;
annotation('textbox',[axpos(1,1)-0.105 axpos(1,2)-0.04 lab_size],'String','(a)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');


%%% Make figure
subplot('Position',axpos(2,:));
for n_rb=1:N_rb
%   semilogx(tau_mean,Tbt(n_rb,:)/1e6,'-','Color',colororder(n_rb,:),'Marker',markershapes{n_rb},'MarkerFaceColor',colororder(n_rb,:));
  plot(tau_mean,Tbt(n_rb,:)/1e6,'-','Color',colororder(n_rb,:),'Marker',markershapes{n_rb},'MarkerFaceColor',colororder(n_rb,:));
  if (n_rb == 1)
  hold on;
  end
end
hold off;
% xlabel('Wind Stress (N/m^2)');
ylabel('Barotropic transport (Sv)');
set(gca,'XTick',tau_ticks);
set(gca,'XLim',[0 0.5]);
set(gca,'YLim',[-10 100]);
set(gca,'FontSize',fontsize);
grid on;
legstr = {};
for n_rb=1:N_rb
  legstr = {legstr{:},['r$_{\mathrm{b}}$=',num2str(lin_drag(n_rb)*10000,'%.1f'),'$\times$10$^{-4}$']}; 
end
leghandle = legend(legstr,'Location','NorthWest');
set(leghandle,'interpreter','latex');
annotation('textbox',[axpos(2,1)-0.105 axpos(2,2)-0.04 lab_size],'String','(b)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

%%% Make figure
subplot('Position',axpos(3,:));
for n_rb=1:N_rb
  semilogx(tau_mean,Tbc(n_rb,:)/1e6,'o-','Color',colororder(n_rb,:),'Marker',markershapes{n_rb},'MarkerFaceColor',colororder(n_rb,:));
  if (n_rb == 1)
  hold on;
  end
end
hold off;
xlabel('Wind Stress (N/m^2)');
ylabel('Baroclinic transport (Sv)');
set(gca,'XTick',tau_ticks);
set(gca,'XLim',[0 0.5]);
set(gca,'YLim',[-10 100]);
set(gca,'FontSize',fontsize);
grid on;
annotation('textbox',[axpos(3,1)-0.105 axpos(3,2)-0.04 lab_size],'String','(c)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');








% 
% %%% Make figure
% figure(2);
% for n_tm=1:N_tm
%   plot(quad_drag,Ttot(:,n_tm)/1e6,'o-');
%   if (n_tm == 1)
%   hold on;
%   end
% end
% hold off;
% xlabel('Quadratic drag coefficient');
% ylabel('Total transport (Sv)');
% grid on;
% legstr = {};
% for n_tm=1:N_tm
%   legstr = {legstr{:},['\tau=',num2str(tau_mean(n_tm)),' N/m^2']}; 
% end
% legend(legstr,'Location','SouthEast');
% 
% 
% 
% 
% %%% Make figure
% figure(4);
% for n_tm=1:N_tm
%   plot(quad_drag,Tbt(:,n_tm)/1e6,'o-');
%   if (n_tm == 1)
%   hold on;
%   end
% end
% hold off;
% xlabel('Quadratic drag coefficient');
% ylabel('Barotropic transport (Sv)');
% grid on;
% legstr = {};
% for n_tm=1:N_tm
%   legstr = {legstr{:},['\tau=',num2str(tau_mean(n_tm)),' N/m^2']}; 
% end
% legend(legstr,'Location','SouthEast');
% 
% 
% 
% %%% Make figure
% figure(6);
% for n_tm=1:N_tm
%   plot(quad_drag,Tbc(:,n_tm)/1e6,'o-');
%   if (n_tm == 1)
%   hold on;
%   end
% end
% hold off;
% xlabel('Quadratic drag coefficient');
% ylabel('Baroclinic transport (Sv)');
% grid on;
% legstr = {};
% for n_tm=1:N_tm
%   legstr = {legstr{:},['\tau=',num2str(tau_mean(n_tm)),' N/m^2']}; 
% end
% legend(legstr,'Location','SouthEast');










figure(108);
clf;
set(gcf,'Position',[1082   306   500   1000]);






%%% Make figure
subplot('Position',axpos(1,:));
for n_rb=1:N_rb
  loglog(tau_mean,kap(n_rb,:),'o-','Color',colororder(n_rb,:),'Marker',markershapes{n_rb},'MarkerFaceColor',colororder(n_rb,:));
%   plot(tau_mean,kap(n_Cd,:),'o-');
  if (n_rb == 1)
  hold on;
  end
end
loglog(tau_mean(7:10),400*(tau_mean(7:10)/0.05).^.5,'k--');
loglog(tau_mean(7:10),100*(tau_mean(7:10)/0.05).^.8,'k--');
% loglog(tau_mean(7:10),1000*(tau_mean(7:10)/0.1).^.4,'k--');
% loglog(tau_mean(2:4),200*(tau_mean(2:4)/0.017).^1,'k--');
hold off;
% xlabel('Wind Stress (N/m^2)');
ylabel('Transient eddy diffusivity (m^2/s)');
set(gca,'XTick',tau_ticks);
set(gca,'XLim',[0 0.5]);
set(gca,'FontSize',fontsize);
grid on;
annotation('textbox',[axpos(1,1)-0.105 axpos(1,2)-0.04 lab_size],'String','(a)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
text(0.075,650,'$\tau^{0.5}$','interpreter','latex','FontSize',fontsize);
text(0.075,100,'$\tau^{0.8}$','interpreter','latex','FontSize',fontsize);
% text(0.07,1000,'$\tau^{0.4}$','interpreter','latex','FontSize',fontsize);
% text(0.017,150,'$\tau^{1}$','interpreter','latex','FontSize',fontsize);



%%% Make figure
subplot('Position',axpos(2,:));
for n_rb=1:N_rb
  loglog(tau_mean,EKE_tot(n_rb,:),'o-','Color',colororder(n_rb,:),'Marker',markershapes{n_rb},'MarkerFaceColor',colororder(n_rb,:));
  if (n_rb == 1)
  hold on;
  end
end
loglog(tau_mean(9:12),1.2e-3*(tau_mean(9:12)/0.1).^.6,'k--');
loglog(tau_mean(7:10),8e-3*(tau_mean(7:10)/0.1).^1,'k--');
hold off;
% xlabel('Wind Stress (N/m^2)');
ylabel('Transient eddy kinetic energy (m^2/s^2)');
set(gca,'XTick',tau_ticks);
set(gca,'XLim',[0 0.5]);
set(gca,'FontSize',fontsize);
grid on;
legstr = {};
for n_rb=1:N_rb
  legstr = {legstr{:},['r$_{\mathrm{b}}$=',num2str(lin_drag(n_rb)*10000,'%.1f'),'$\times$10$^{-4}$']};  
end
leghandle = legend(legstr,'Location','NorthWest');
set(leghandle,'interpreter','latex');
annotation('textbox',[axpos(2,1)-0.105 axpos(2,2)-0.04 lab_size],'String','(b)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
text(0.14,1.2e-3,'$\tau^{0.6}$','interpreter','latex','FontSize',fontsize);
text(0.07,8e-3,'$\tau^{1}$','interpreter','latex','FontSize',fontsize);

%%% Make figure
subplot('Position',axpos(3,:));
for n_rb=1:N_rb
%   semilogx(tau_mean,SKE_tot(n_Cd,:),'o-');
  loglog(tau_mean,SKE_tot(n_rb,:),'o-','Color',colororder(n_rb,:),'Marker',markershapes{n_rb},'MarkerFaceColor',colororder(n_rb,:));
  if (n_rb == 1)
  hold on;
  end
end
loglog(tau_mean(7:10),1.5e-3*(tau_mean(7:10)/0.1).^1,'k--');
loglog(tau_mean(7:10),3.5e-4*(tau_mean(7:10)/0.1).^1.4,'k--');
hold off;
xlabel('Wind Stress (N/m^2)');
ylabel('Standing wave kinetic energy (m^2/s^2)');
set(gca,'XTick',tau_ticks);
set(gca,'XLim',[0 0.5]);
set(gca,'FontSize',fontsize);
grid on;
annotation('textbox',[axpos(3,1)-0.105 axpos(3,2)-0.04 lab_size],'String','(c)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
text(0.07,1.5e-3,'$\tau^1$','interpreter','latex','FontSize',fontsize);
text(0.07,1.3e-4,'$\tau^{1.4}$','interpreter','latex','FontSize',fontsize);


% 
% %%% Make figure
% figure(8);
% for n_tm=1:N_tm
%   plot(quad_drag,kap(:,n_tm),'o-');
%   if (n_tm == 1)
%   hold on;
%   end
% end
% hold off;
% xlabel('Quadratic drag coefficient');
% ylabel('Transient eddy diffusivity (m^2/s)');
% title(['C_d = ',num2str(quad_drag),', r_b = ',num2str(lin_drag),' m/s']);
% grid on;
% legstr = {};
% for n_tm=1:N_tm
%   legstr = {legstr{:},['\tau=',num2str(tau_mean(n_tm)),' N/m^2']}; 
% end
% legend(legstr,'Location','SouthEast');
% 
% 
% %%% Make figure
% figure(10);
% for n_tm=1:N_tm
%   plot(quad_drag,EKE_tot(:,n_tm),'o-');
%   if (n_tm == 1)
%   hold on;
%   end
% end
% hold off;
% xlabel('Quadratic drag coefficient');
% ylabel('EKE (m^2/s^2)');
% title(['C_d = ',num2str(quad_drag),', r_b = ',num2str(lin_drag),' m/s']);
% grid on;
% legstr = {};
% for n_tm=1:N_tm
%   legstr = {legstr{:},['\tau=',num2str(tau_mean(n_tm)),' N/m^2']}; 
% end
% legend(legstr,'Location','SouthEast');
% 
% 
% %%% Make figure
% figure(12);
% for n_tm=1:N_tm
%   plot(quad_drag,SKE_tot(:,n_tm),'o-');
%   if (n_tm == 1)
%   hold on;
%   end
% end
% hold off;
% xlabel('Quadratic drag coefficient');
% ylabel('SKE (m^2/s^2)');
% title(['C_d = ',num2str(quad_drag),', r_b = ',num2str(lin_drag),' m/s']);
% grid on;
% legstr = {};
% for n_tm=1:N_tm
%   legstr = {legstr{:},['\tau=',num2str(tau_mean(n_tm)),' N/m^2']}; 
% end
% legend(legstr,'Location','SouthEast');
% 
% 

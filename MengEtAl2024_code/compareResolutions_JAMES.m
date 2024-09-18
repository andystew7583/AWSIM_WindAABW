%%%
%%% compareResolutions_JAMES.m
%%%
%%% Compares transport/eddy properties between resolutions with identical
%%% forcing.
%%%


% run_name = 'ACC_AABW_ML_doubleMOC_veryhires';
% [KE_veryhires,PE_veryhires,E_veryhires,Z_veryhires,t_veryhires] = readEZfile(local_home_dir,run_name);
% 
% tlen = t_veryhires(end);
% tlen = 20*t1year;
% 
% [Tacc_veryhires,Taabw_veryhires,Taaiw_veryhires,EKE_veryhires,tt_veryhires, ...
%   Taaiw_mean_veryhires,Taaiw_eddy_veryhires,Taabw_mean_veryhires,Taabw_eddy_veryhires] ...
%     = calcMOCTimeSeries(local_home_dir,run_name,0,t_veryhires(end));
% 
% run_name = 'ACC_AABW_ML_doubleMOC_hires';
% [KE_hires,PE_hires,E_hires,Z_hires,t_hires] = readEZfile(local_home_dir,run_name);
% [Tacc_hires,Taabw_hires,Taaiw_hires,EKE_hires,tt_hires, ...
%   Taaiw_mean_hires,Taaiw_eddy_hires,Taabw_mean_hires,Taabw_eddy_hires] ...
%      = calcMOCTimeSeries(local_home_dir,run_name,200*t1year,200*t1year+tlen);
% 
% run_name = 'ACC_AABW_ML_doubleMOC_lores';
% [KE_lores,PE_lores,E_lores,Z_lores,t_lores] = readEZfile(local_home_dir,run_name);
% [Tacc_lores,Taabw_lores,Taaiw_lores,EKE_lores,tt_lores, ...
%   Taaiw_mean_lores,Taaiw_eddy_lores,Taabw_mean_lores,Taabw_eddy_lores] ...
%      = calcMOCTimeSeries(local_home_dir,run_name,0,tlen);
% 
% 
% save('CompareResolutions.mat', ...
%   'KE_lores','PE_lores','E_lores','Z_lores','t_lores', ...
%   'Tacc_lores','Taabw_lores','Taaiw_lores','EKE_lores','tt_lores', ...
%   'KE_hires','PE_hires','E_hires','Z_hires','t_hires', ...
%   'Tacc_hires','Taabw_hires','Taaiw_hires','EKE_hires','tt_hires', ...
%   'KE_veryhires','PE_veryhires','E_veryhires','Z_veryhires','t_veryhires', ...
%   'Tacc_veryhires','Taabw_veryhires','Taaiw_veryhires','EKE_veryhires','tt_veryhires', ...
%   'Taaiw_mean_veryhires','Taaiw_eddy_veryhires','Taabw_mean_veryhires','Taabw_eddy_veryhires', ...
%   'Taaiw_mean_hires','Taaiw_eddy_hires','Taabw_mean_hires','Taabw_eddy_hires', ...
%   'Taaiw_mean_lores','Taaiw_eddy_lores','Taabw_mean_lores','Taabw_eddy_lores' ...
%   );

load('CompareResolutions.mat');


%%% Load horizontal grids
run_name = 'ACC_AABW_ML_doubleMOC_lores';
loadParams;
yy_lores = yy_v;
run_name = 'ACC_AABW_ML_doubleMOC_hires';
loadParams;
yy_hires = yy_v;
run_name = 'ACC_AABW_ML_doubleMOC_veryhires';
loadParams;
yy_veryhires = yy_v;


figure(30);
plot(yy_lores,Taabw_eddy_lores);
hold on
plot(yy_hires,Taabw_eddy_hires);
plot(yy_veryhires,Taabw_eddy_veryhires);
hold off;

figure(31);
plot(yy_lores,Taaiw_eddy_lores);
hold on
plot(yy_hires,Taaiw_eddy_hires);
plot(yy_veryhires,Taaiw_eddy_veryhires);
hold off;

figure(12);
clf;
idx = find(diff(t_lores)<0,1,'last');
plot(t_lores(idx+1:end)/t1year,PE_lores(idx+1:end)-PE_lores(idx+1));
hold on;
idxrange = find(t_hires>=200*t1year & t_hires<300*t1year);
plot(t_hires(idxrange)/t1year-200,PE_hires(idxrange)-PE_hires(idxrange(1)))
idx = find(diff(t_veryhires)<0,1,'last');
plot(t_veryhires(idx+1:end)/t1year,PE_veryhires(idx+1:end)-PE_veryhires(idx+1))
hold off;
set(gca,'XLim',[0 tlen/t1year]);
legend('d=12km','d=6km','d=3km');

figure(13);
clf;
idx = find(diff(t_lores)<0,1,'last');
plot(t_lores(idx+1:end)/t1year,KE_lores(idx+1:end))
hold on;
idxrange = find(t_hires>=200*t1year & t_hires<300*t1year);
plot(t_hires(idxrange)/t1year-200,KE_hires(idxrange))
idx = find(diff(t_veryhires)<0,1,'last');
plot(t_veryhires(idx+1:end)/t1year,KE_veryhires(idx+1:end));
hold off;
set(gca,'XLim',[0 tlen/t1year]);
legend('d=12km','d=6km','d=3km');

figure(14);
plot(tt_lores/t1year,Taabw_lores/1e6);
hold on;
plot(tt_hires/t1year-200,Taabw_hires/1e6);
plot(tt_veryhires/t1year,Taabw_veryhires/1e6);
hold off;
set(gca,'XLim',[0 tlen/t1year]);
legend('d=12km','d=6km','d=3km');


figure(19);
plot(tt_lores/t1year,Tacc_lores/1e6);
hold on;
plot(tt_hires/t1year-200,Tacc_hires/1e6);
plot(tt_veryhires/t1year,Tacc_veryhires/1e6);
hold off;
set(gca,'XLim',[0 tlen/t1year]);
legend('d=12km','d=6km','d=3km');

figure(15);
plot(tt_lores/t1year,EKE_lores);
hold on;
plot(tt_hires/t1year-200,EKE_hires);
plot(tt_veryhires/t1year,EKE_veryhires);
hold off;
set(gca,'XLim',[0 tlen/t1year]);
legend('d=12km','d=6km','d=3km');

res_plot = [1 2 3];
EKE_plot = [mean(EKE_lores) mean(EKE_hires) mean(EKE_veryhires)] ./ sum(sum(-hhb*dx*dy));
Taaiw_plot = [std(Taaiw_lores) std(Taaiw_hires) std(Taaiw_veryhires)];
Taabw_plot = [std(Taabw_lores) std(Taabw_hires) std(Taabw_veryhires)];
Tacc_plot = [std(Tacc_lores) std(Tacc_hires) std(Tacc_veryhires)];

figure(16);
plot(res_plot,EKE_plot,'o-');

figure(17);
plot(res_plot,Taaiw_plot,'o-');



figure(20);
plot(res_plot,Taabw_plot,'o-');


figure(18);
plot(res_plot,Tacc_plot/1e6,'o-');



rho0 = 1000;
axpos = zeros(4,4);
framepos = [ 1001         797         810         542 ];
axpos(1,:) = [0.08 0.57 0.4 0.38];
axpos(2,:) = [0.58 0.57 0.4 0.38];
axpos(3,:) = [0.08 0.09 0.4 0.38];
axpos(4,:) = [0.58 0.09 0.4 0.38];
axlabels = {'(a)','(b)','(c)','(d)','(e)','(f)'};
lab_size = [0.05 0.05];
fontsize = 14;
xticklabels = {'256x128','512x256','1024x512'};
figure(21);
set(gcf,'Position',framepos);
subplot('Position',axpos(1,:));
plot(res_plot,EKE_plot,'o-','LineWidth',1.5);
title('Volume-averaged eddy kinetic energy');
set(gca,'XTick',res_plot);
set(gca,'XTickLabel',xticklabels);
set(gca,'FontSize',fontsize);
ylabel('m^2/s^2');
axis([0.5 3.5 0 0.01]);
grid on;
subplot('Position',axpos(2,:));
plot(res_plot,Tacc_plot/1e6,'o-','LineWidth',1.5);
title('ACC transport std. dev.');
set(gca,'XTick',res_plot);
set(gca,'XTickLabel',xticklabels);
set(gca,'FontSize',fontsize);
axis([0.5 3.5 0 10]);
ylabel('Sv');
grid on;
subplot('Position',axpos(3,:));
plot(res_plot,Taaiw_plot/1e6,'o-','LineWidth',1.5);
title('AAIW transport std. dev.');
set(gca,'XTick',res_plot);
set(gca,'XTickLabel',xticklabels);
set(gca,'FontSize',fontsize);
xlabel('Grid size');
axis([0.5 3.5 0 1.2]);
ylabel('Sv');
grid on;
subplot('Position',axpos(4,:));
plot(res_plot,Taabw_plot/1e6,'o-','LineWidth',1.5);
title('AABW transport std. dev.');
set(gca,'XTick',res_plot);
set(gca,'XTickLabel',xticklabels);
set(gca,'FontSize',fontsize);
xlabel('Grid size');
axis([0.5 3.5 0 1.2]);
ylabel('Sv');
grid on;

rho0 = 1000;
axpos = zeros(4,4);
framepos = [ 1001         797         810         800 ];
axpos(1,:) = [0.08 0.71 0.4 0.26];
axpos(2,:) = [0.58 0.71 0.4 0.26];
axpos(3,:) = [0.08 0.385 0.4 0.26];
axpos(4,:) = [0.58 0.385 0.4 0.26];
axpos(5,:) = [0.08 0.06 0.4 0.26];
axpos(6,:) = [0.58 0.06 0.4 0.26];
fontsize = 14;
xticklabels = {'256x128','512x256','1024x512'};
figure(22);
clf;
set(gcf,'Position',framepos);
subplot('Position',axpos(1,:));
plot(res_plot,EKE_plot,'o-','LineWidth',1.5);
title('Volume-averaged eddy kinetic energy');
set(gca,'XTick',res_plot);
set(gca,'XTickLabel',xticklabels);
set(gca,'FontSize',fontsize);
ylabel('m^2/s^2');
axis([0.5 3.5 0 0.01]);
grid on;
subplot('Position',axpos(2,:));
plot(res_plot,Tacc_plot/1e6,'o-','LineWidth',1.5);
title('ACC transport std. dev.');
set(gca,'XTick',res_plot);
set(gca,'XTickLabel',xticklabels);
set(gca,'FontSize',fontsize);
axis([0.5 3.5 0 10]);
ylabel('Sv');
grid on;
subplot('Position',axpos(3,:));
plot(res_plot,Taaiw_plot/1e6,'o-','LineWidth',1.5);
title('AAIW transport std. dev.');
set(gca,'XTick',res_plot);
set(gca,'XTickLabel',xticklabels);
set(gca,'FontSize',fontsize);
xlabel('Grid size');
axis([0.5 3.5 0 1.2]);
ylabel('Sv');
grid on;
subplot('Position',axpos(4,:));
plot(res_plot,Taabw_plot/1e6,'o-','LineWidth',1.5);
title('AABW transport std. dev.');
set(gca,'XTick',res_plot);
set(gca,'XTickLabel',xticklabels);
set(gca,'FontSize',fontsize);
xlabel('Grid size');
axis([0.5 3.5 0 1.2]);
ylabel('Sv');
grid on;
subplot('Position',axpos(5,:));
plot(yy_lores/m1km,Taabw_eddy_lores/1e6,'LineWidth',1.5);
hold on
plot(yy_hires/m1km,Taabw_eddy_hires/1e6,'LineWidth',1.5);
plot(yy_veryhires/m1km,Taabw_eddy_veryhires/1e6,'LineWidth',1.5);
hold off;
set(gca,'FontSize',fontsize);
xlabel('y (km)');
ylabel('Eddy transpor of AABW (Sv)'); 
legend(xticklabels);
subplot('Position',axpos(6,:));
plot(yy_lores/m1km,Taaiw_eddy_lores/1e6,'LineWidth',1.5);
hold on
plot(yy_hires/m1km,Taaiw_eddy_hires/1e6,'LineWidth',1.5);
plot(yy_veryhires/m1km,Taaiw_eddy_veryhires/1e6,'LineWidth',1.5);
hold off;
set(gca,'FontSize',fontsize);
xlabel('y (km)');
ylabel('Eddy transport of AAIW (Sv)');

%%% Add axis labels
for cntr = 1:size(axpos,1)
  annotation('textbox',[axpos(cntr,1)-0.07 axpos(cntr,2)-0.06 lab_size],'String',axlabels{cntr},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
end


function [Tacc,Taabw,Taaiw,EKE,tt_sub,...
  Taaiw_mean,Taaiw_eddy,Taabw_mean,Taabw_eddy] ...
    = calcMOCTimeSeries (local_home_dir,run_name,tmin,tmax)

  %%% Load parameters for this experiment
  loadParams;
  
  %%% Define range of latitudes over which to compute AABW transport
  ymin = 100*m1km;
  ymax = Ly-100*m1km;
  yminidx = find(yy_v>ymin,1,'first');
  ymaxidx = find(yy_v>ymax,1,'first');
  
  %%% Precompute iteration numbers and output times
  iters = n0+1:1:n0+Nframes-1;
  tt = startTime + (iters-n0)*dt_s;  
  subidx = find((tt>=tmin) & (tt<tmax));
  tt_sub = tt(subidx);
  sublen = length(subidx);
  
  %%% At each time iteration...
  cntr = 0;  
  Tacc = zeros(1,sublen);
  Taabw = zeros(1,sublen);
  Taaiw = zeros(1,sublen);
  TKE = zeros(1,sublen);
  MKE = 0; 
  h = zeros(Nlay,Nx,Ny);
  v = zeros(Nlay,Nx,Ny);
  u = zeros(Nlay,Nx,Ny);  
  hu_tavg = zeros(Nlay,Nx,Ny);  
  hv_tavg = zeros(Nlay,Nx,Ny);  
  hw_tavg = zeros(Nlay,Nx,Ny);  
  hs_tavg = zeros(Nlay,Nx,Ny);  
  u_tavg = zeros(Nlay,Nx,Ny);  
  v_tavg = zeros(Nlay,Nx,Ny);  
  enddata = false;
  for n=subidx
    
    disp(['Iteration ',num2str(n)]);
    
    %%% No need to proceed if we've run out of data
    if (enddata)
      continue;
    end
        
    %%% Load instantaneous model state
    for k=1:Nlay
      data_file = fullfile(dirpath,[OUTN_H,num2str(k-1),'_n=',num2str(n),'.dat']);
      tmp = readOutputFile(data_file,Nx,Ny);
      if (isempty(tmp))
        tmp = NaN*ones(Nx,Ny);
      end
      h(k,:,:) = tmp;
    end  
    for k=1:Nlay
      data_file = fullfile(dirpath,[OUTN_U,num2str(k-1),'_n=',num2str(n),'.dat']);
      tmp = readOutputFile(data_file,Nx,Ny);
      if (isempty(tmp))
        tmp = NaN*ones(Nx,Ny);
      end      
      u(k,:,:) = tmp;
    end
    for k=1:Nlay
      data_file = fullfile(dirpath,[OUTN_V,num2str(k-1),'_n=',num2str(n),'.dat']);
      tmp = readOutputFile(data_file,Nx,Ny);
      if (isempty(tmp))
        tmp = NaN*ones(Nx,Ny);
      end      
      v(k,:,:) = tmp;
    end
    h_v = 0.5*(h(:,:,1:Ny)+h(:,:,[2:Ny 1])); %%% Interpolate to v-points
    h_u = 0.5*(h(:,1:Nx,:)+h(:,[2:Nx 1],:)); %%% Interpolate to u-points    
    
    %%% No need to proceed if we've run out of data
    if (enddata)
      continue;
    end
    
    %%% Increment iteration counter
    cntr = cntr + 1;

    %%% Meridional transport decomposition
    hv_int = squeeze(sum(h_v.*v.*dx,2));

    %%% Zonal transport
    hu_int = squeeze(sum(h_u.*u.*dy,3));  

    %%% Compute transports     
    Taabw(cntr) = squeeze(mean(hv_int(Nlay,yminidx:ymaxidx),2));    
    Taaiw(cntr) = squeeze(mean(hv_int(1,yminidx:ymaxidx),2));
    Tacc(cntr) = sum(hu_int(:,1),1);
    
    %%% Compute total kinetic energy
    TKE(cntr) = sum(sum(sum(0.5*(h_v.*v.^2 + h_u.*u.^2)*dx*dy,3),2),1);
    
    %%% Increment time averages
    hu_tavg = hu_tavg + h_u.*u;
    hv_tavg = hv_tavg + h_v.*v;
    hw_tavg = hw_tavg + h_u;
    hs_tavg = hs_tavg + h_v;
    u_tavg = u_tavg + u;
    v_tavg = v_tavg + v;

  end
  
  %%% Compute time averages
  hu_tavg = hu_tavg / cntr;
  hv_tavg = hv_tavg / cntr;
  hw_tavg = hw_tavg / cntr;
  hs_tavg = hs_tavg / cntr;
  u_tavg = u_tavg / cntr;
  v_tavg = v_tavg / cntr;
  
  %%% Compute mean kinetic energy
  MKE = 0.5*(hu_tavg.^2./hw_tavg + hv_tavg.^2./hs_tavg);
  MKE = sum(sum(sum(MKE*dx*dy,3),2),1);
  
  hv_mean = hs_tavg .* v_tavg;
  hv_mean_int = squeeze(sum(hv_mean.*dx,2));  
  hv_tavg_int = squeeze(sum(hv_tavg.*dx,2));  
  Taaiw_mean = hv_mean_int(1,:);
  Taabw_mean = hv_mean_int(Nlay,:);
  Taaiw_eddy = hv_tavg_int(1,:) - Taaiw_mean;
  Taabw_eddy = hv_tavg_int(Nlay,:) - Taabw_mean;
  
  %%% Compute EKE
  EKE = TKE - MKE;
  
end
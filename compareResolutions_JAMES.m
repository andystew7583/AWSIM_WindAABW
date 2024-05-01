

figure(12);
clf;
run_name = 'ACC_AABW_ML_doubleMOC_lores';
[KE,PE,E,Z,t]=readEZfile(local_home_dir,run_name);
idx = find(diff(t)<0,1,'last');
plot(t(idx+1:end)/t1year,PE(idx+1:end)-PE(idx+1));
hold on;
run_name = 'ACC_AABW_ML_doubleMOC_hires';
[KE,PE,E,Z,t]=readEZfile(local_home_dir,run_name);
idxrange = find(t>=200*t1year & t<300*t1year);
plot(t(idxrange)/t1year-200,PE(idxrange)-PE(idxrange(1)))
run_name = 'ACC_AABW_ML_doubleMOC_veryhires';
[KE,PE,E,Z,t]=readEZfile(local_home_dir,run_name);
idx = find(diff(t)<0,1,'last');
plot(t(idx+1:end)/t1year,PE(idx+1:end)-PE(idx+1))
hold off;
set(gca,'XLim',[0 100]);
legend('d=12km','d=6km','d=3km');

figure(13);
clf;
run_name = 'ACC_AABW_ML_doubleMOC_lores';
[KE,PE,E,Z,t]=readEZfile(local_home_dir,run_name);
idx = find(diff(t)<0,1,'last');
plot(t(idx+1:end)/t1year,KE(idx+1:end))
hold on;
run_name = 'ACC_AABW_ML_doubleMOC_hires';
[KE,PE,E,Z,t]=readEZfile(local_home_dir,run_name);
idxrange = find(t>=200*t1year & t<300*t1year);
plot(t(idxrange)/t1year-200,KE(idxrange))
run_name = 'ACC_AABW_ML_doubleMOC_veryhires';
[KE,PE,E,Z,t]=readEZfile(local_home_dir,run_name);
idx = find(diff(t)<0,1,'last');
plot(t(idx+1:end)/t1year,KE(idx+1:end));
hold off;
set(gca,'XLim',[0 100]);
legend('d=12km','d=6km','d=3km');
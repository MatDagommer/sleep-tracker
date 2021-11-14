% TP_Dreem_Sleep_complete
% Novembre 2018
%
%
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A - Initialisation and visualisation

%info to fullfill
record_ref = 2112262;
userAge = 20;
folderDreem = pwd;

%folders and files
folderRecords = fullfile(folderDreem,'Records');
folderNight = fullfile(folderRecords, num2str(record_ref));

filerecord = fullfile(folderNight, [num2str(record_ref) '_signals.mat']);
filehypno = fullfile(folderNight, [num2str(record_ref) '_hypno.mat']);
filespectro = fullfile(folderNight, [num2str(record_ref) '_spectrogram.mat']);

%load data
load(filerecord);
load(filehypno);
load(filespectro);


%% Visualisation
% - use the example below to visualise signals of your night on Matlab
% - you can also display signals you will compute or tag specif moments of the night

%signals to display
display(labels_eeg);
signals{1} = eeg{1};            titles{1} = labels_eeg{1};      y_borders{1} = [-200 200]*1e4;
signals{2} = eeg{2};            titles{2} = labels_eeg{2};      y_borders{2} = [-200 200]*1e4;
signals{3} = eeg{5};            titles{3} = labels_eeg{5};      y_borders{3} = [-200 200];
signals{4} = pulse_oximeter;    titles{4} = 'pulse oxymeter';   y_borders{4} = [-200 200]*1e3;
signals{5} = breathing;         titles{5} = 'respiration';      y_borders{5} = [-3 3]*1e2;


%Viewer initialization
vdr = ViewerDreemRecord(signals,'DisplayWindow',[20 20]);
for ch=1:vdr.nb_channel
    vdr.set_title(titles{ch}, ch);
    vdr.set_ymin(ch, y_borders{ch}(1)); 
    vdr.set_ymax(ch, y_borders{ch}(2));
end
% show slow waves on virtual channel
vdr.add_rectangles(slowwaves_start, slowwaves_duration, 3); 

%plot
vdr.run_window


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% B - Description of your night

%% Sleep metrics
% Compute those sleep metrics on your night
% Eventually, they can be retrieved manually with csv files

%SOL = Sleep Onset Latency
% two first consecutives epochs of sleep, ie first period of sleep which is at least 1 min long
idx = find(SleepStage(:,3)~=5 & SleepStage(:,2)-SleepStage(:,1)>=60,1);
sol = SleepStage(idx,1);

%WASO
% total duration of wake period after sleep onset and before waking
idx = find(SleepStage(:,3)==5 & SleepStage(:,1)>=sol);
idx(idx==size(SleepStage,1))=[];
waso = sum(SleepStage(idx,2) - SleepStage(idx,1));

%TST - total sleep time
idx = find(SleepStage(:,3)~=5);
tst = sum(SleepStage(idx,2) - SleepStage(idx,1));

%TIB - time in bed
% total duration of record
tib = SleepStage(end,2); % total time of record

%Sleep efficiency
% percentage of total sleep on duration of record
sleep_efficiency = 100 * tst / tib;

%last sleep - end of sleep period
idx = find(SleepStage(:,3)~=5,1,'last');
night_duration = SleepStage(idx,2);


%% Sleep stage duration and ratio (ratio of each substage on TST)
% Compute those ratios and durations
% can be retrieved manually with csv files

%N1
idx = find(SleepStage(:,3)==1);
n1_duration = sum(SleepStage(idx,2) - SleepStage(idx,1));
n1_ratio = n1_duration / tst;

%N2
idx = find(SleepStage(:,3)==2);
n2_duration = sum(SleepStage(idx,2) - SleepStage(idx,1));
n2_ratio = n2_duration / tst;

%N3
idx = find(SleepStage(:,3)==3);
n3_duration = sum(SleepStage(idx,2) - SleepStage(idx,1));
n3_ratio = n3_duration / tst;

%REM
idx = find(SleepStage(:,3)==4);
rem_duration = sum(SleepStage(idx,2) - SleepStage(idx,1));
rem_ratio = rem_duration / tst;


%% Slow waves : quantification and density across the night

%number and rate
% compute the total number of slow waves
nb_slowwaves = length(slowwaves_start); % total number of slow waves

% compute the slow wave rate (per min) = number of slow waves / duration of sleep
slowwaves_rate = nb_slowwaves / (tst/60); % slow waves per min of sleep


% write an algorythm, using a 'for loop' or a 'while loop' to compute the
% evolution of slow-waves density along the night 

%parameters to compute slow waves density
record_duration = max(eeg{1}(:,1));
windowsize = 60;
nb_minute = floor(record_duration/windowsize);

timestamps = 0:nb_minute;
density_slowwaves = [];
for t=1:length(timestamps)
    density_slowwaves(t) = sum(slowwaves_start>=timestamps(t)*60 & slowwaves_start<(timestamps(t)+1)*60);
end


%% Plot the temporal description of your night
% spectrogram
% heart rate
% hypnogram
% slow waves density


figure, hold on

% Spectrogram of virtual channel (time in hours)
subplot(4,1,1), hold on
imagesc(times_spectro/3600, freq_spectro, log(spectrogram)'), hold on
axis xy, ylabel('frequency'), hold on
set(gca, 'xlim', [0 max(times_spectro/3600)], 'ylim', [0.5 15]);
caxis([-5 10]),
title(['Spectrogram' label])

% heart rate evolution (time in hours)
subplot(4,1,2), hold on
plot(heart_rate(:,1)/3600, heart_rate(:,2), 'k', 'linewidth',2), hold on,
xlim([0 max(heart_rate(:,1)/3600)]), ylim([40 90]),
line([0.18 0.18], ylim, 'color','r', 'LineWidth', 2);
line([1.78 1.78], ylim, 'color','r', 'LineWidth', 2);
line([3.43 3.43], ylim, 'color','r', 'LineWidth', 2);
line([5.29 5.29], ylim, 'color','r', 'LineWidth', 2);
line([6.85 6.85], ylim, 'color','r', 'LineWidth', 2);
ylabel('pulse/min')
title('Heart rate'); 

% Hypnogram (time in hours)
subplot(4,1,3), hold on
ylabel_substage = {'N3','N2','N1','REM','WAKE'};
ytick_substage = [1 1.5 2 3 4]; %ordinate in graph
plot(t_hypno/3600, y_hypno,'k', 'linewidth',2), hold on,
xlim([0 max(t_hypno/3600)]), ylim([0.5 5]), set(gca,'Ytick',ytick_substage,'YTickLabel',ylabel_substage), hold on,
line([0.18 0.18], ylim, 'color','r', 'LineWidth', 2);
line([1.78 1.78], ylim, 'color','r', 'LineWidth', 2);
line([3.43 3.43], ylim, 'color','r', 'LineWidth', 2);
line([5.29 5.29], ylim, 'color','r', 'LineWidth', 2);
line([6.85 6.85], ylim, 'color','r', 'LineWidth', 2);
title('Hypnogram');

%density of slow waves (time in hours)
subplot(4,1,4), hold on
plot(timestamps/60, density_slowwaves, 'k', 'linewidth',2), hold on,
xlim([0 max(timestamps/60)]), xlabel('Time (h)'), ylabel('density of slow waves per min')

line([0.18 0.18], ylim, 'color','r', 'LineWidth', 2);
line([1.78 1.78], ylim, 'color','r', 'LineWidth', 2);
line([3.43 3.43], ylim, 'color','r', 'LineWidth', 2);
line([5.29 5.29], ylim, 'color','r', 'LineWidth', 2);
line([6.85 6.85], ylim, 'color','r', 'LineWidth', 2);

title('Density of slow-waves'); 


% % how to add lines between sleep cycles on a graph
% % example for a cycle border at 1.3h
% subplot(4,1,3), hold on
% line([1.3 1.3], ylim, 'color','r');

%% Random Plots



%% Homeostasis


% quantify slow waves in each cycle: 
% - numbers of slow-waves
% - sleep cycle duration
% - slow waves rate per cycle (number/duration)

cycles_border = [sol 1.78*3600 3.43*3600 5.29*3600 6.85*3600 night_duration];
for i=1:length(cycles_border)-1
    nb_slowwaves_cycle(i) = sum(slowwaves_start>=cycles_border(i) & slowwaves_start<cycles_border(i+1)); 
    duration_cycle(i) = (cycles_border(i+1) - cycles_border(i))/3600;
end
slowwaves_density_cycle = nb_slowwaves_cycle ./ duration_cycle;

% quantify slow waves at the beginning and at the end of the night 
% - in the 2 hours after sleep onset
% - in the 2 last hours before waking up

nb_slowwaves_beginning = sum(slowwaves_start>=sol & slowwaves_start<sol+2*3600); 
nb_slowwaves_end = sum(slowwaves_start>=night_duration-2*3600 & slowwaves_start<night_duration); 


% You have already computed and plot the slow waves density evolution with time 
% - make a linear regression on this curve
% - find the slope
% - plot the line of this regression

%linear regression
x_density = timestamps/60;
y_density = density_slowwaves;

p_density = polyfit(x_density, y_density, 1);
y_slope = polyval(p_density,x_density);

%slope
mySlope = p_density(1);

%plot the line
subplot(4,1,4), hold on
plot(x_density, y_slope, 'b', 'linewidth',1), hold on,
ylim([0 50])


%% respiration

% - via the hypnogram and the variable SleepStage, choose periods of REM, N3 and Wake
% - compute and plot the fast-fourier transform (fft) of the breathing
% during these period (use the function fft_signaux_headband.m)

idx = find(SleepStage(:,3)==5 & SleepStage(:,2) - SleepStage(:,1)>60, 1)
wake_period = [SleepStage(idx,1) SleepStage(idx,2)]

idx = find(SleepStage(:,3)==3 & SleepStage(:,2) - SleepStage(:,1)>60, 1)
N3_period = [SleepStage(idx,1) SleepStage(idx,2)]

idx = find(SleepStage(:,3)==4 & SleepStage(:,2) - SleepStage(:,1)>60, 1)
rem_period = [SleepStage(idx,1) SleepStage(idx,2)]

figure, hold on

%[power_fft, frequencies_fft] = fft_signaux_headband(breathing, 1440, 4110);
[power_fft, frequencies_fft] = fft_signaux_headband(breathing, wake_period(1), wake_period(2));
plot(frequencies_fft(1:500),power_fft(1:500),'r'); hold on

%[power_fft, frequencies_fft] = fft_signaux_headband(breathing, 17010, 17700);
[power_fft, frequencies_fft] = fft_signaux_headband(breathing, N3_period(1), N3_period(2));
plot(frequencies_fft(1:500),power_fft(1:500),'b'); hold on

%[power_fft, frequencies_fft] = fft_signaux_headband(breathing, 20000, 21000);
[power_fft, frequencies_fft] = fft_signaux_headband(breathing, rem_period(1), rem_period(2));
plot(frequencies_fft(1:200),power_fft(1:200),'g'); hold on

ylabel('densité')
xlabel('fréquence (Hz)')
legend('Eveil', 'N3', 'REM')
title('Transformée de Fourier du rythme de respiration (accéléromètre)'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% C - Compare to a large database 
% In this section, you will compare the data above with the data provided
% by Dreem database of users
% 

%% Sleep metrics

%load sleep metrics data from DREEM
load(fullfile(folderDreem, 'population_data','PopulationData.mat'));

%plot options
scattersize = 25; % size of each point of the scatter plot
scatter_color = [0 0.5 0]; % color of each point of the scatter plot


figure, hold on

% TST Total sleep time
% - plot the data from the dreem database of TST (var: 'subject_sleep')
% - plot the linear regression for these data
% - is their a correlation between age and total sleep time ?
% - insert your data in this graph  

%plot init
subplot(2,2,1), hold on
xlabel('âge'), ylabel('durée (h)'), title('Total sleep time (TST)')
%data from DREEM
scatter(subject_age, subject_sleep, scattersize, scatter_color), hold on
% linear regression
pfit = polyfit(subject_age, subject_sleep,1);%fit
yfit = pfit(1)*subject_age + pfit(2);
plot(subject_age, yfit, 'r', 'LineWidth', 2), hold on
% correlation ?
[rc, pv] = corrcoef(subject_age,subject_sleep); % correlations
rc = rc(1,2); pv=pv(1,2);
disp(['correlation r = ' num2str(round(rc,2)) ' - p = ' num2str(pv)])
if pv>0.05
    disp('no significant correlation')
end

idx = find(SleepStage(:,3)==5 & SleepStage(:,1)>=sol);
idx(idx==size(SleepStage,1))=[];
waso = sum(SleepStage(idx,2) - SleepStage(idx,1));

%TST - total sleep time
idx = find(SleepStage(:,3)~=5);
tst = sum(SleepStage(idx,2) - SleepStage(idx,1));

%TIB - time in bed
% total duration of record
tib = SleepStage(end,2); % total time of record

%Sleep efficiency
% percentage of total sleep on duration of record
sleep_efficiency = 100 * tst / tib;

%last sleep - end of sleep period
idx = find(SleepStage(:,3)~=5,1,'last');
night_duration = SleepStage(idx,2);

plot(20, tst/3600, 'bx', 'linewidth', 3, 'MarkerSize', 20)

% Sleep Efficiency 
%plot init
subplot(2,2,2), hold on
xlabel('âge'), ylabel('durée (h)'), title('Sleep Efficiency (SE)')
%data from DREEM
scatter(subject_age, subject_sleepeff, scattersize, scatter_color), hold on
% linear regression
pfit = polyfit(subject_age, subject_sleepeff,1);%fit
yfit = pfit(1)*subject_age + pfit(2);
plot(subject_age, yfit, 'r', 'LineWidth', 2), hold on
% correlation ?
[rc, pv] = corrcoef(subject_age,subject_sleepeff); % correlations
rc = rc(1,2); pv=pv(1,2);
disp(['correlation r = ' num2str(round(rc,2)) ' - p = ' num2str(pv)])
if pv>0.05
    disp('no significant correlation')
end

plot(20, sleep_efficiency, 'bx', 'linewidth', 3, 'MarkerSize', 20)

% Sleep Onset Latency
%plot init
subplot(2,2,3), hold on
xlabel('âge'), ylabel('durée (h)'), title('Sleep Onset Latency (SOL)')
%data from DREEM
scatter(subject_age, subject_sol, scattersize, scatter_color), hold on
% linear regression
pfit = polyfit(subject_age, subject_sol,1);%fit
yfit = pfit(1)*subject_age + pfit(2);
plot(subject_age, yfit, 'r', 'LineWidth', 2), hold on
% correlation ?
[rc, pv] = corrcoef(subject_age,subject_sol); % correlations
rc = rc(1,2); pv=pv(1,2);
disp(['correlation r = ' num2str(round(rc,2)) ' - p = ' num2str(pv)])
if pv>0.05
    disp('no significant correlation')
end

plot(20, sol, 'bx', 'linewidth', 3, 'MarkerSize', 20)

% WASO
%plot init
subplot(2,2,4), hold on
xlabel('âge'), ylabel('durée (h)'), title('Wake After Sleep Onset (WASO)')
%data from DREEM
scatter(subject_age, subject_waso, scattersize, scatter_color), hold on
% linear regression
pfit = polyfit(subject_age, subject_waso,1);%fit
yfit = pfit(1)*subject_age + pfit(2);
plot(subject_age, yfit, 'r', 'LineWidth', 2), hold on
% correlation ?
[rc, pv] = corrcoef(subject_age,subject_waso); % correlations
rc = rc(1,2); pv=pv(1,2);
disp(['correlation r = ' num2str(round(rc,2)) ' - p = ' num2str(pv)])
if pv>0.05
    disp('no significant correlation')
end

plot(20, waso, 'bx', 'linewidth', 3, 'MarkerSize', 20)

%% Ondes lentes
%load sleep metrics data from DREEM
load(fullfile(folderDreem, 'population_data','PopulationDataSlowWaves.mat'));

%plot options
scattersize = 25; % size of each point of the scatter plot
scatter_color = [0 0.5 0]; % color of each point of the scatter plot


figure, hold on

% nb slow waves
%plot init
subplot(3,1,1), hold on
xlabel('âge'), ylabel('nombre ondes'), title('Nombre total ondes lentes')
%data from DREEM
scatter(subject_age, subject_nbslowwaves, scattersize, scatter_color), hold on
% linear regression
pfit = polyfit(subject_age, subject_nbslowwaves,1);%fit
yfit = pfit(1)*subject_age + pfit(2);
plot(subject_age, yfit, 'r', 'LineWidth', 2), hold on
% correlation ?
[rc, pv] = corrcoef(subject_age,subject_nbslowwaves); % correlations
rc = rc(1,2); pv=pv(1,2);
disp(['correlation r = ' num2str(round(rc,2)) ' - p = ' num2str(pv)])
if pv>0.05
    disp('no significant correlation')
end

plot(20, nb_slowwaves, 'bx', 'linewidth', 3, 'MarkerSize', 20)



% densité d'ondes lentess 
%plot init
subplot(3,1,2), hold on
xlabel('âge'), ylabel('densité (/min)'), title('Densité Ondes Lentes')
%data from DREEM
scatter(subject_age, subject_swdensity, scattersize, scatter_color), hold on
% linear regression
pfit = polyfit(subject_age, subject_swdensity,1);%fit
yfit = pfit(1)*subject_age + pfit(2);
plot(subject_age, yfit, 'r', 'LineWidth', 2), hold on
% correlation ?
[rc, pv] = corrcoef(subject_age,subject_swdensity); % correlations
rc = rc(1,2); pv=pv(1,2);
disp(['correlation r = ' num2str(round(rc,2)) ' - p = ' num2str(pv)])
if pv>0.05
    disp('no significant correlation')
end

plot(20, nb_slowwaves/night_duration*60, 'bx', 'linewidth', 3, 'MarkerSize', 20)


% pente de l'homeostasie
%plot init
subplot(3,1,3), hold on
xlabel('âge'), ylabel('durée (h)'), title('Pente Homéostasie')
%data from DREEM
scatter(subject_age, subject_densityslope, scattersize, scatter_color), hold on
% linear regression
pfit = polyfit(subject_age, subject_densityslope,1);%fit
yfit = pfit(1)*subject_age + pfit(2);
plot(subject_age, yfit, 'r', 'LineWidth', 2), hold on
% correlation ?
[rc, pv] = corrcoef(subject_age,subject_densityslope); % correlations
rc = rc(1,2); pv=pv(1,2);
disp(['correlation r = ' num2str(round(rc,2)) ' - p = ' num2str(pv)])
if pv>0.05
    disp('no significant correlation')
end

plot(20, mySlope, 'bx', 'linewidth', 3, 'MarkerSize', 20)


%% Sleep stage duration and ratio (ratio of each substage on TST)
% Compute those ratios and durations
% can be retrieved manually with csv files

%N1
idx = find(SleepStage(:,3)==1);
n1_duration = sum(SleepStage(idx,2) - SleepStage(idx,1));
n1_ratio = n1_duration / tst;

%N2
idx = find(SleepStage(:,3)==2);
n2_duration = sum(SleepStage(idx,2) - SleepStage(idx,1));
n2_ratio = n2_duration / tst;

%N3
idx = find(SleepStage(:,3)==3);
n3_duration = sum(SleepStage(idx,2) - SleepStage(idx,1));
n3_ratio = n3_duration / tst;

%REM
idx = find(SleepStage(:,3)==4);
rem_duration = sum(SleepStage(idx,2) - SleepStage(idx,1));
rem_ratio = rem_duration / tst;


%% Slow waves : quantification and density across the night

%number and rate
% compute the total number of slow waves
nb_slowwaves = length(slowwaves_start); % total number of slow waves

% compute the slow wave rate (per min) = number of slow waves / duration of sleep
slowwaves_rate = nb_slowwaves / (tst/60); % slow waves per min of sleep


% write an algorythm, using a 'for loop' or a 'while loop' to compute the
% evolution of slow-waves density along the night 

%parameters to compute slow waves density
record_duration = max(eeg{1}(:,1));
windowsize = 60;
nb_minute = floor(record_duration/windowsize);

timestamps = 0:nb_minute;
density_slowwaves = [];
for t=1:length(timestamps)
    density_slowwaves(t) = sum(slowwaves_start>=timestamps(t)*60 & slowwaves_start<(timestamps(t)+1)*60);
end


%% Plot the temporal description of your night
% spectrogram
% heart rate
% hypnogram
% slow waves density


figure, hold on

% Spectrogram of virtual channel (time in hours)
subplot(4,1,1), hold on
imagesc(times_spectro/3600, freq_spectro, log(spectrogram)'), hold on
axis xy, ylabel('frequency'), hold on
set(gca, 'xlim', [0 max(times_spectro/3600)], 'ylim', [0.5 15]);
caxis([-5 10]),
title(['Spectrogram' label])

% heart rate evolution (time in hours)
subplot(4,1,2), hold on
plot(heart_rate(:,1)/3600, heart_rate(:,2), 'k', 'linewidth',2), hold on,
xlim([0 max(heart_rate(:,1)/3600)]), ylim([40 90]),
ylabel('pulse/min')
title('Heart rate'); 

% Hypnogram (time in hours)
subplot(4,1,3), hold on
ylabel_substage = {'N3','N2','N1','REM','WAKE'};
ytick_substage = [1 1.5 2 3 4]; %ordinate in graph
plot(t_hypno/3600, y_hypno,'k', 'linewidth',2), hold on,
xlim([0 max(t_hypno/3600)]), ylim([0.5 5]), set(gca,'Ytick',ytick_substage,'YTickLabel',ylabel_substage), hold on,
title('Hypnogram');

%density of slow waves (time in hours)
subplot(4,1,4), hold on
plot(timestamps/60, density_slowwaves, 'k', 'linewidth',2), hold on,
xlim([0 max(timestamps/60)]), xlabel('Time (h)'), ylabel('density of slow waves per min')
title('Density of slow-waves'); 


% % how to add lines between sleep cycles on a graph
% % example for a cycle border at 1.3h
% subplot(4,1,3), hold on
% line([1.3 1.3], ylim, 'color','r');



%% Homeostasis


% quantify slow waves in each cycle: 
% - numbers of slow-waves
% - sleep cycle duration
% - slow waves rate per cycle (number/duration)

cycles_border = [sol 1.3*3600 2.5*3600 5.3*3600 7.6*3600 night_duration];
for i=1:length(cycles_border)-1
    nb_slowwaves_cycle(i) = sum(slowwaves_start>=cycles_border(i) & slowwaves_start<cycles_border(i+1)); 
    duration_cycle(i) = (cycles_border(i+1) - cycles_border(i))/3600;
end
slowwaves_density_cycle = nb_slowwaves_cycle ./ duration_cycle;

% quantify slow waves at the beginning and at the end of the night 
% - in the 2 hours after sleep onset
% - in the 2 last hours before waking up

nb_slowwaves_beginning = sum(slowwaves_start>=sol & slowwaves_start<sol+2*3600); 
nb_slowwaves_end = sum(slowwaves_start>=night_duration-2*3600 & slowwaves_start<night_duration); 


% You have already computed and plot the slow waves density evolution with time 
% - make a linear regression on this curve
% - find the slope
% - plot the line of this regression

%linear regression
x_density = timestamps/60;
y_density = density_slowwaves;

p_density = polyfit(x_density, y_density, 1);
y_slope = polyval(p_density,x_density);

%slope
mySlope = p_density(1);

%plot the line
subplot(4,1,4), hold on
plot(x_density, y_slope, 'b', 'linewidth',1), hold on,
ylim([0 30])


%% respiration

% - via the hypnogram and the variable SleepStage, choose periods of REM, N3 and Wake
% - compute and plot the fast-fourier transform (fft) of the breathing
% during these period (use the function fft_signaux_headband.m)


figure, hold on

[power_fft, frequencies_fft] = fft_signaux_headband(breathing, 1440, 4110);
plot(frequencies_fft,power_fft,'r'); hold on

[power_fft, frequencies_fft] = fft_signaux_headband(breathing, 17010, 17700);
plot(frequencies_fft,power_fft,'b'); hold on

[power_fft, frequencies_fft] = fft_signaux_headband(breathing, 20000, 21000);
plot(frequencies_fft,power_fft,'g'); hold on



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% C - Compare to a large database 
% In this section, you will compare the data above with the data provided
% by Dreem database of users
% 

%% Sleep metrics

%load sleep metrics data from DREEM
load(fullfile(folderDreem, 'population_data','PopulationData.mat'));

%plot options
scattersize = 25; % size of each point of the scatter plot
scatter_color = [0.3 0.3 0.3]; % color of each point of the scatter plot

figure, hold on

% TST Total sleep time
% - plot the data from the dreem database of TST (var: 'subject_sleep')
% - plot the linear regression for these data
% - is their a correlation between age and total sleep time ?
% - insert your data in this graph  

%plot init
subplot(2,2,1), hold on
xlabel('age'), ylabel('duration (h)'), title('Total sleep time (TST)')
%data from DREEM
scatter(subject_age, subject_sleep, scattersize, scatter_color, 'filled'), hold on
% linear regression
pfit = polyfit(subject_age, subject_sleep,1);%fit
yfit = pfit(1)*subject_age + pfit(2);
plot(subject_age, yfit), hold on
% correlation ?
[rc, pv] = corrcoef(subject_age,subject_sleep); % correlations
rc = rc(1,2); pv=pv(1,2);
disp(['correlation r = ' num2str(round(rc,2)) ' - p = ' num2str(pv)])
if pv>0.05
    disp('no significant correlation')
end

% insert your data

% insert your data
plot(myAge,tst/3600,'bx','linewidth',3,'MarkerSize',20), %my point




%% Slow waves metrics







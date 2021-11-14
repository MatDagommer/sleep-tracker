% TP_Dreem_Sleep
% Novembre 2018
%
%
%


%% INIT

%info to fullfill
record_ref = 2112262;
userAge = 20;

%folders and files
folderDreem = pwd;
folderRecords = fullfile(folderDreem,'Records');
folderNight = fullfile(folderRecords, num2str(record_ref));

filerecord = fullfile(folderNight, [num2str(record_ref) '_signals.mat']);
filehypno = fullfile(folderNight, [num2str(record_ref) '_hypno.mat']);
filespectro = fullfile(folderNight, [num2str(record_ref) '_spectrogram.mat']);

%load data
load(filerecord);
load(filehypno);
load(filespectro);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A - Visualisation
% - use https://dreem-viewer.rythm.co
%  OR
% - go to TP_Dreem_Sleep_visualisation.m



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% B - Description of your night


%% Sleep metrics
% Compute those sleep metrics on your night
% Eventually, they can be retrieved manually with csv files


% SOL = Sleep Onset Latency
% two first consecutives epochs of sleep, ie first period of sleep which is at least 1 min long

idx = find(SleepStage(:,3)~=5 & SleepStage(:,2)-SleepStage(:,1)>=60,1);
sol = SleepStage(idx,1); 

% waso - wake after sleep onset
% total duration of wake period after sleep onset and before waking
idx = find(SleepStage(:,3)==5 & SleepStage(:,1)>=sol);
idx(idx==size(SleepStage,1))=[];
waso = sum(SleepStage(idx,2) - SleepStage(idx,1));

% tst - total sleep time

idx = find(SleepStage(:,3)~=5);
tst = sum(SleepStage(idx,2) - SleepStage(idx,1));

% tib - time in bed
% total duration of record

tib = SleepStage(end,2);

% sleep_efficiency
% percentage of total sleep on duration of record

sleep_efficiency = 100 * tst / tib;

% night_duration
% last sleep moment - end of sleep period

idx = find(SleepStage(:,3)~=5,1,'last');
night_duration = SleepStage(idx,2);

%% Sleep stage duration and ratio (ratio of each substage on TST)
% Compute those ratios and durations
% can be retrieved manually with csv files


%N1 (n1_duration n1_ratio)


%N2 (n2_duration n2_ratio)


%N3 (n3_duration n3_ratio)


%REM (rem_duration rem_ratio)



%% Slow waves : quantification and density across the night

% number and rate
% - compute the total number of slow waves

% - compute the slow wave rate (per min) = number of slow waves / duration of sleep

nbTotalSlowWaves = size(slowwaves_start);
density = nbTotalSlowWaves/night_duration;


nbIntervals1min = floor(night_duration/60);
slowwaves_density = zeros(nbIntervals1min,1);



% write an algorythm, using a 'for loop' or a 'while loop' to compute the
% evolution of slow-waves density (per min) along the night :
% slowwaves_density = f(t)


for i=1:nbIntervals1min
    for j=1:nbTotalSlowWaves
        if i*nbIntervals1min <= slowwaves_start(j) < (i+1)*nbIntervals1min
            slowwaves_density(i) = slowwaves_density(i) + 1;
        end
    end
end


%% Plot the temporal description of your night
% Uncomment the following section to plot your :
% - spectrogram
% - heart rate
% - hypnogram
% Add the data you have just computed on slow waves density in the last graph :
% use the plot function:  plot(t, slowwaves_density)
% NB : time in hours


% % % UNCOMMENT the following section to plot hypnogram 
% % BEGINNING OF SECTION 
% 
% figure, hold on
% 
% % Spectrogram of virtual channel (time in hours)
% subplot(4,1,1), hold on
% imagesc(times_spectro/3600, freq_spectro, log(spectrogram)'), hold on
% axis xy, ylabel('frequency'), hold on
% set(gca, 'xlim', [0 max(times_spectro/3600)], 'ylim', [0.5 15]);
% caxis([-5 10]),
% title(['Spectrogram - ' label])
% 
% % heart rate evolution (time in hours)
% subplot(4,1,2), hold on
% plot(heart_rate(:,1)/3600, heart_rate(:,2), 'k', 'linewidth',2), hold on,
% xlim([0 max(heart_rate(:,1)/3600)]), ylim([40 90]),
% ylabel('pulse/min')
% title('Heart rate'); 
% 
% % Hypnogram (time in hours)
% subplot(4,1,3), hold on
% ylabel_substage = {'N3','N2','N1','REM','WAKE'};
% ytick_substage = [1 1.5 2 3 4]; %ordinate in graph
% plot(t_hypno/3600, y_hypno,'k', 'linewidth',2), hold on,
% xlim([0 max(t_hypno/3600)]), ylim([0.5 5]), set(gca,'Ytick',ytick_substage,'YTickLabel',ylabel_substage), hold on,
% title('Hypnogram');
% 
% %density of slow waves (time in hours)
% % - add your data on slow-waves density
% subplot(4,1,4), hold on
% 
% xlabel('Time (h)'), ylabel('density of slow waves per min'), ylim([0 30]),
% title('Density of slow-waves'); 
% 
% 
% %  END OF SECTION


% Using your graph, determine the borders between your sleep cycles
% plot them on the graph

%HELP:
% % how to add lines between sleep cycles on the hypnogram
% % example for a cycle border at 1.3h
% subplot(4,1,3), hold on
% line([1.3 1.3], ylim, 'color','r');



%% Homeostasis

% Slow-waves and SWA naturally decrease along the night: we will quantify this phenomenon

% Quantify slow waves in each cycle: 
% - numbers of slow-waves
% - sleep cycle duration
% - slow waves rate per cycle (number/duration)


% Quantify slow waves at the beginning and at the end of the night 
% - in the 2 hours after sleep onset
% - in the 2 last hours before waking up



% You have already computed and plot the slow waves density evolution with time 
% - make a linear regression on this curve
% - find the slope (mySlope)
% - plot the line of this regression

%linear regression


%mySlope


%plot the line
subplot(4,1,4), hold on



%% respiration

% - via the hypnogram and the variable SleepStage, choose periods of REM, N3 and Wake
% - compute and plot the fast-fourier transform (fft) of the breathing
% during these period (use the function fft_signaux_headband.m)




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
% - plot the data from the Dreem database of TST (var: 'subject_sleep')
% - plot the linear regression for these data
% - is their a correlation between age and total sleep time ?
% - insert your data in this graph  

%plot init
subplot(2,2,1), hold on

xlabel('age'), ylabel('duration (h)'), title('Total sleep time (TST)')

%data from DREEM (use scatter.m)


% linear regression




% correlation ?


% insert your data



%DO THE SAME WITH OTHER Sleep metrics
% subplot these graphs together

% Sleep Efficiency (var: 'subject_sleepeff')
% SOL Sleep Onset Latency (var: 'subject_sol')
% WASO wake after sleep onset (var: 'subject_waso')

% N1 percentage (var: 'subject_N1percentage')
% N2 percentage (var: 'subject_N2percentage')
% N3 percentage (var: 'subject_N3percentage')
% REM percentage (var: 'subject_REMpercentage')



%% Slow waves metrics

%load slow waves metrics data from DREEM
load(fullfile(folderDreem, 'population_data','PopulationDataSlowWaves.mat'));


%DO THE SAME WITH Slow waves metrics
% subplot these graphs together

% Number of Slow Waves (var: 'subject_nbslowwaves')
% Slow Waves Rate (var: 'subject_swdensity')
% Slope of slow waves density (var: 'subject_densityslope')











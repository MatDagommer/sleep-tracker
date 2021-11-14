% TP_Dreem_Sleep_visualisation
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

%load data
load(filerecord);


%% A - Visualisation
% - use the example below to visualise signals of your night on Matlab
% - you can also display signals you will compute or tag specif moments of the night

%signals to display
display(labels_eeg);
signals{1} = eeg{1};            titles{1} = labels_eeg{1};      y_borders{1} = [-200 200];
signals{2} = eeg{2};            titles{2} = labels_eeg{2};      y_borders{2} = [-200 200];
signals{3} = eeg{5};            titles{3} = labels_eeg{5};      y_borders{3} = [-200 200];
signals{4} = pulse_oximeter;    titles{4} = 'pulse oxymeter';   y_borders{4} = [-50 50]*4;
signals{5} = breathing;         titles{5} = 'respiration';      y_borders{5} = [-3 3];


%Viewer initialization
vdr = ViewerDreemRecord(signals,'DisplayWindow',[20 20]);
for ch=1:vdr.nb_channel
    vdr.set_title(titles{ch}, ch);
    vdr.set_ymin(ch, y_borders{ch}(1)); 
    vdr.set_ymax(ch, y_borders{ch}(2));
    
    
end
% show slow waves on virtual channel
vdr.add_rectangles(slowwaves_start, slowwaves_duration, 3); 

%launch and plot the viewer
vdr.run_window



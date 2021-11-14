%LoadPathTPDreem

res=pwd;

cd(res)
addpath(genpath(res))

eval(['cd(''',res,''')'])

clear res

%colors
cmap = colormap('jet');
set(groot,'DefaultFigureColorMap',cmap);
close all
clear cmap

%window for figures
set(0,'DefaultFigureWindowStyle','docked')
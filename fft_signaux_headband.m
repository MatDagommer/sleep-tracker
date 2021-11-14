function [power_fft, frequencies_fft] = fft_signaux_headband(signal, start_time, stop_time)
%
%    [power_fft, frequencies_fft] = fft_signaux_headband(signal, start_time, stop_time)
%
% INPUT:
% - signal                  2-dimensions vector ([timestamps values]): signals on which the fft is applied. e.g: EEG, respiration... 
% - start_time              start timestamp of the sample to analyze
% - end_time                end timestamp of the sample to analyze 
%
% OUTPUT:
% - power_fft               power values of the fft              
% - frequencies_fft         corresponding frequencies
%


% init
x_sig = signal(:,1);
y_sig = signal(:,2);
fs = round(1/median(diff(x_sig)));

%restrict to period
idx = x_sig>=start_time & x_sig<=stop_time;
y_sig = y_sig(idx);
if mod(length(y_sig),2)
    y_sig = y_sig(1:end-1);
end

%% FFT
L = length(y_sig);
Y = y_sig - mean(y_sig);
Y = fft(Y);

P2 = abs(Y/L);
P1 = P2(1:(L/2)+1);
P1(2:end) = 2*P1(2:end);

%result
frequencies_fft = fs * (0:(L/2)) / L;
power_fft = smooth(P1,5);



end
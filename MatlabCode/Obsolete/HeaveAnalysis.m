%%
clc;clearvars;close all;





addpath ./Mausund200701_221241
load Heave.mat;
load GpsFix.mat;
rmpath ./Mausund200701_221241
%%


heave_signal = Heave.value(Heave.src_ent == 39);
heave_timestamp = Heave.timestamp(Heave.src_ent == 39) - Heave.timestamp(1);
m = 6*120;
c = 1;

lpFilt = designfilt('lowpassfir','PassbandFrequency',0.3, ...
    'StopbandFrequency',0.5,'PassbandRipple',0.5, ...
    'StopbandAttenuation',65,'DesignMethod','kaiserwin','SampleRate',2);
filtSig = filter(lpFilt,heave_signal);

%filtSig = heave_signal;
freqz = zeros(length(10000 : m : length(heave_timestamp) - 10000),1);
speeds = zeros(length(10000 : m : length(heave_timestamp) - 10000),1);
%smoothsignal = filter(lpFilt,heave_small_signal);


for i = 10000 : m : length(heave_timestamp) - 10000
    signal_timest_vec = heave_timestamp(i:i+m);
    heave_small_signal = filtSig(i:i+m);

    [pks,locs] = findpeaks(heave_small_signal,signal_timest_vec,'MinPeakProminence',0.1,'MinPeakHeight',0.1 ,'MinPeakDistance',3);
    % where you feed the heave signal and its timestamp vector.
    avg_periods_from_peaks = mean(diff(locs));

    avg_freq_hz = 1./avg_periods_from_peaks;
    avg_freq_radians_per_second = 2*pi*avg_freq_hz;
    freqz(c) = avg_freq_radians_per_second;
    speeds(c) = mean(GpsFix.sog(i:i+m));
    c = c +1 ;
    
    
end
figpeaks = figure;
plot(signal_timest_vec, heave_small_signal);
hold on
plot(locs, pks, 'pg', 'MarkerFaceColor','red','MarkerSize',15)
xlabel('Time [s]')
ylabel('Heave Displacement [m]', 'Interpreter','latex', 'FontSize', 11)
title(['Heave Displacement and Peaks'])
grid on

avg_periods_from_peaks = mean(diff(locs));

avg_freq_hz = 1./avg_periods_from_peaks;
avg_freq_radians_per_second = 2*pi*avg_freq_hz;

figure;
scatter(freqz, speeds)

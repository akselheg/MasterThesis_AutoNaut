%%
clc;clearvars;close all;
addpath ./Mausund200701_221241
load Heave.mat;
rmpath ./Mausund200701_221241
%%
heave_signal = Heave.value(Heave.src_ent == 39);
heave_timestamp = Heave.timestamp(Heave.src_ent == 39) - Heave.timestamp(1);
signal_timest_vec = heave_timestamp(3000:3500);
heave_small_signal = heave_signal(3000:3500);
[pks,locs] = findpeaks(heave_small_signal,signal_timest_vec,'MinPeakProminence',0.01,'MinPeakHeight',0.01 ,'MinPeakDistance',1);
% where you feed the heave signal and its timestamp vector.

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
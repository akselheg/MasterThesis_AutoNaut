clearvars; close all; clc;
avrager = 60*6;
Fs = 2;
T = 1/Fs;
L = 2*avrager;
hzz1= [];
hzz2= [];
rms1 = [];
rms2 = [];
ampl = [];
sog_data = [];
%% load data
for i = 1:7
    disp('Loading new data')
    %% load data
    if i == 1
        path = './Mausund200701_181204/';
        addpath(path);
        gpsFix = load('GpsFix.mat');
        RelativeWind = load('RelativeWind.mat');
        EulerAngles = load('EulerAngles.mat');
        Heave = load('./Mausund200701_181204/Heave.mat');
        
        disp('Done loading data')
        rmpath(path)
    end
    if i == 2
        path = './Mausund200701_221241/';
        addpath(path);
        gpsFix = load('GpsFix.mat');
        RelativeWind = load('RelativeWind.mat');
        EulerAngles = load('EulerAngles.mat');
        Heave = load('Heave.mat');
        rmpath(path)
        disp('Done loading data')
    end
    if i == 3
        path = './Mausund200703_080820/';
        addpath(path);
        gpsFix = load('GpsFix.mat');
        RelativeWind = load('RelativeWind.mat');
        EulerAngles = load('EulerAngles.mat');
        Heave = load('./Mausund200703_080820/Heave.mat');
        rmpath(path)
        disp('Done loading data')
    end
    if i == 4
        path = './Mausund200703_132548/';
        addpath(path);
        gpsFix = load('GpsFix.mat');
        RelativeWind = load('RelativeWind.mat');
        EulerAngles = load('EulerAngles.mat');
        Heave = load('./Mausund200703_132548/Heave.mat');
        rmpath(path)
        disp('Done loading data')
    end
    if i == 5
        path = './Mausund200705_120030/';
        addpath(path);
        gpsFix = load('GpsFix.mat');
        RelativeWind = load('RelativeWind.mat');
        EulerAngles = load('EulerAngles.mat');
        Heave = load('./Mausund200705_120030/Heave.mat');
        rmpath(path)
        disp('Done loading data')
    end
    if i == 6
        path = './Mausund200706_154608/';
        addpath(path);
        gpsFix = load('GpsFix.mat');
        RelativeWind = load('RelativeWind.mat');
        EulerAngles = load('EulerAngles.mat');
        Heave = load('./Mausund200706_154608/Heave.mat');
        rmpath(path)
        disp('Done loading data')
    end
     if i == 7
        path = './Mausund200709_53748/';
        addpath(path);
        gpsFix = load('GpsFix.mat');
        RelativeWind = load('RelativeWind.mat');
        EulerAngles = load('EulerAngles.mat');
        Heave = load('./Mausund200709_53748/Heave.mat');
        rmpath(path)
        disp('Done loading data')
     end
gpsFix = gpsFix.GpsFix;
EulerAngles = EulerAngles.EulerAngles;
Heave = Heave.Heave;
heave = Heave.value;
pitch = EulerAngles.phi;
time = gpsFix.timestamp;
heave1 = heave(1:2:end);
heave2 = heave(2:2:end);


%%
% testsample = pitch(1000-avrager:1000+avrager+L-1);
% Y = fft(testsample,L)/L;
% Y(1) = 0;
% ampl= 2*abs(Y(1:ceil(L/2)));
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% f = Fs * (1:(L/2))/L;
% [mp, i] = max(abs(Y).^2);
% hz = f(i);
% hzz1 = cat(1,hzz1,hz);
% freq = meanfreq(testsample, Fs);
% stop;
    for m = (10*120) : length(time) - (10*120)
        if ~mod(round((m),1),avrager)
        sog = mean(gpsFix.sog(m-avrager:m+avrager));
        sog_data = cat(1, sog_data,sog);
    %%
        testsample1 = pitch(m-avrager:m+avrager+L-1);
        Y = fft(testsample1,L)/L;
        ampl = cat(1,ampl, mean(Y));
        Y(1) = 0;
        %ampl= 2*abs(Y(1:ceil(L/2)));
        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        f = Fs * (1:(L/2))/L;
        [mp, i] = max(abs(Y(1:ceil(L/2))).^2);
        hz = meanfreq(testsample1, Fs);
        hzz1 = cat(1,hzz1, hz);
        rms1 = cat(1, rms1, rms(testsample1));
        %
        testsample2 = heave2(m-avrager:m+avrager+L-1);
        Y = fft(testsample2,L)/L;
        Y(1) = 0;
        ampl= 2*abs(Y(1:ceil(L/2)));
        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        f = Fs * (1:(L/2))/L;

        [mp, i] = max(abs(Y(1:ceil(L/2))).^2);
        hz2 = meanfreq(testsample2, Fs);
        hzz2 = cat(1,hzz2, hz2);
        rms2 = cat(1, rms2, rms(testsample2));
        end
    end
end
figure;
scatter(hzz1,sog_data)
figure;
scatter(hzz2,sog_data)
figure;
scatter(rms1,sog_data)
figure;
scatter(rms2,sog_data)
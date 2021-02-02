clearvars; close all; clc;
avrager = 2*120;
Fs = 2;
T = 1/Fs;
L = 2*avrager;
hzz1= [];
hzz2= [];
rms1 = [];
rms2 = [];
ampl = [];
wind_data = [];
sog_data = [];
%% load data
for i = 1:4
    disp('Loading new data')
    %% load data
    if i == 1
        path = './Mausund200701_181204/';
        addpath(path);
        gpsFix = load('GpsFix.mat');
        RelativeWind = load('RelativeWind.mat');
        EulerAngles = load('EulerAngles.mat');
        Heave = load('Heave.mat');
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
        disp('Done loading data')
        rmpath(path)
    end
    if i == 3
        path = './Mausund200703_080820/';
        addpath(path);
        gpsFix = load('GpsFix.mat');
        RelativeWind = load('RelativeWind.mat');
        EulerAngles = load('EulerAngles.mat');
        Heave = load('Heave.mat');
        disp('Done loading data')
        rmpath(path)
    end
    if i == 4
        path = './Mausund200703_132548/';
        addpath(path);
        gpsFix = load('GpsFix.mat');
        RelativeWind = load('RelativeWind.mat');
        EulerAngles = load('EulerAngles.mat');
        Heave = load('Heave.mat');
        disp('Done loading data')
        rmpath(path)
    end
    if i == 5
        path = './Mausund200705_120030/';
        addpath(path);
        gpsFix = load('GpsFix.mat');
        RelativeWind = load('RelativeWind.mat');
        EulerAngles = load('EulerAngles.mat');
        Heave = load('Heave.mat');
        disp('Done loading data')
        rmpath(path)
    end
    if i == 6
        path = './Mausund200706_154608/';
        addpath(path);
        gpsFix = load('GpsFix.mat');
        RelativeWind = load('RelativeWind.mat');
        EulerAngles = load('EulerAngles.mat');
        Heave = load('Heave.mat');
        disp('Done loading data')
        rmpath(path)
    end
    if i == 7
        path = './Mausund200709_53748/';
        addpath(path);
        gpsFix = load('GpsFix.mat');
        RelativeWind = load('RelativeWind.mat');
        EulerAngles = load('EulerAngles.mat');
        Heave = load('Heave.mat');
        disp('Done loading data')
        rmpath(path)
    end
    if i == 8
        path = './Mausund200710_092034/';
        addpath(path);
        addpath(path);
        gpsFix = load('GpsFix.mat');
        RelativeWind = load('RelativeWind.mat');
        EulerAngles = load('EulerAngles.mat');
        Heave = load('Heave.mat');
        disp('Done loading data')
        rmpath(path)
    end
    if i == 9
        path = './Mausund200703_215938/';
        addpath(path);
        gpsFix = load('GpsFix.mat');
        RelativeWind = load('RelativeWind.mat');
        EulerAngles = load('EulerAngles.mat');
        Heave = load('Heave.mat');
        disp('Done loading data')
        rmpath(path)
    end
    if i == 10
        path = './Mausund200712_220202/';
        addpath(path);
        gpsFix = load('GpsFix.mat');
        RelativeWind = load('RelativeWind.mat');
        EulerAngles = load('EulerAngles.mat');
        Heave = load('Heave.mat');
        disp('Done loading data')
        rmpath(path)
    end
gpsFix = gpsFix.GpsFix;
EulerAngles = EulerAngles.EulerAngles;
RelativeWind = RelativeWind.RelativeWind;
Heave = Heave.Heave;
heave = Heave.value;
pitch = EulerAngles.phi;
time = gpsFix.timestamp;
heave1 = heave(1:2:end);
heave2 = heave(2:2:end);
messuredRelWindDir = interp1(RelativeWind.timestamp, ssa(RelativeWind.angle,'deg' ),gpsFix.timestamp);
messuredRelWindSpeed = interp1(RelativeWind.timestamp, RelativeWind.speed,gpsFix.timestamp);


%%
% testsample = pitch(1000-avrager:1000+avrager+L-1);
% Y = fft(testsample,L)/L;
% Y(1) = 0;
% ampl= 2*abs(Y(1:ceil(L/2)));
% P2 = abs(Y/L); rms(testsample1)
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% f = Fs * (1:(L/2))/L;
% [mp, i] = max(abs(Y).^2);
% hz = f(i);
% hzz1 = cat(1,hzz1,hz);
% freq = meanfreq(testsample, Fs);
% stop;
    for m = (10*120) :avrager: length(gpsFix.sog) - (10*120)
        sog = mean(gpsFix.sog(m-avrager:m+avrager));
        sog_data = cat(1, sog_data,sog);
        wind = mean(messuredRelWindSpeed(m-avrager:m+avrager));
        wind_data= cat(1,wind_data,wind);
    %%  
        testsample1 = pitch(m-avrager:m+avrager);
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
        hzz1 = cat(1,hzz1, f(i));
        rms1 = cat(1, rms1,rms(testsample1));
        %
        testsample2 = heave2(m-avrager:m+avrager);
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
%%
Plot2Dim(sog_data,hzz1, ' ',' ')
Plot2Dim(sog_data,rms1, ' ',' ')
%%
% Plot3Dim(sog_data, hzz1, rms1, 0.05, 0.09, true, 'Ampl', 'Pitch Frequency', 'Vg')
% Plot3Dim(sog_data, hzz2, rms2, 0.4, 0.6, true, 'Ampl', 'Heave Frequency', 'Vg')
% Plot3Dim(sog_data, rms1, hzz1, 0.25, 0.4, true, 'Freq', 'Pitch Amplitude', 'Vg')
% Plot3Dim(sog_data, rms2, hzz2, 0.1, 0.15, true, 'Freq', 'Heave Amplitude', 'Vg')
% %%
% 
% X = [hzz1 rms1 rms2 wind_data ones(length(sog_data),1)];
% w1 = (X'*X)\(X'*sog_data);
% Mdl1 = fitrgp(X(:,1:end-1), sog_data, 'KernelFunction', 'matern52');
% PlotLinear(sog_data, w1, X, 'Vg');
% PlotGaus(sog_data, Mdl1, X(:,1:end-1), 'Vg')

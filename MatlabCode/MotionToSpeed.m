%% Clear Workspace
clc; clearvars;close all;



pitchRatesFreq = [];
pitchRateAmplitudes = [];
pitchFreq = [];
pitchAmplitudes = [];
rollRateFreq = [];
rollRateAmplitudes = [];
rollFreq = [];
rollAmplitudes = [];
heavefreq = [];
heaveamplitude = [];
Sogs = [];

for inputData = 1 : 3
%% Load Data
    if inputData == 1 
        path = './Mausund200701_181204';
    end
    if inputData == 2 
        path = './Mausund200703_062402';
    end
    if inputData == 3 
        path = './Mausund200703_080820';
    end

    
    addpath(path)
    load AngularVelocity.mat
    load Acceleration.mat
    load GpsFix.mat
    load Heave.mat
    load EulerAngles.mat
    rmpath(path)

    %% 
    %plot(Acceleration.x)

    % 123 Hz 

    fs = 123;                   % sampling frequency
    minutes = 1;                % 30 minutes data
    h = 1/fs;                   % sampling time
    N = minutes*60*fs;                 
    t = (0:N-1)*h;   

    timestamp = AngularVelocity.timestamp - AngularVelocity.timestamp(1);
    GpsFix.timestamp = GpsFix.timestamp - GpsFix.timestamp(1);
    Heave.timestamp = Heave.timestamp(Heave.src_ent==39);
    Heave.timestamp = Heave.timestamp - Heave.timestamp(1);
    %Acceleration.z = 

    SpeedOverGround = interp1(GpsFix.timestamp, GpsFix.sog, timestamp);
    HeaveEsttimator = Heave.value(Heave.src_ent==39);
    HeaveEst = interp1(Heave.timestamp, HeaveEsttimator, timestamp, 'pchip');
    EulerAngles.timestamp = EulerAngles.timestamp - EulerAngles.timestamp(1);
    roll = interp1(EulerAngles.timestamp, EulerAngles.phi, timestamp, 'pchip');
    pitch = interp1(EulerAngles.timestamp, EulerAngles.theta, timestamp, 'pchip');

    omega_0 = 0.7823;
    sigma = sqrt(2.5997);      



    M = 200;
    for i = 10*N:N:length(AngularVelocity.y) - 10*N

        X = AngularVelocity.y(i:i+N)';
        n = 2^nextpow2(N);      % pad the input with trailing zeros
        Y = fft(X,n);         % compute the FFT
        P2 = abs(Y/N);          % double-sided spectrum of each signal
        P1 = P2(:,1:n/2+1);     % single-sided spectrum of each signal
        P1(:,2:end-1) = 2*P1(:,2:end-1);
        P1(1) = 0;
        f = 0:(fs/n):(fs/2-fs/n); w = 2*pi*f;  % frequency vectors      
        fitfun = fit(w(1:M)',(P1(1:M)/max(P1(1:M)))','gauss1');% no of samples to plot
        [~,idx] = max(fitfun(w(1:M)));
        pitchRatesFreq = cat(1,pitchRatesFreq,w(idx));
        pitchRateAmplitudes = cat(1,pitchRateAmplitudes, sqrt(2)*rms(X - mean(X)));

        X = AngularVelocity.x(i:i+N)';
        n = 2^nextpow2(N);      % pad the input with trailing zeros
        Y = fft(X,n);         % compute the FFT
        P2 = abs(Y/N);          % double-sided spectrum of each signal
        P1 = P2(:,1:n/2+1);     % single-sided spectrum of each signal
        P1(:,2:end-1) = 2*P1(:,2:end-1);
        P1(1) = 0;
        f = 0:(fs/n):(fs/2-fs/n); w = 2*pi*f;  % frequency vectors       
        fitfun = fit(w(1:M)',(P1(1:M)/max(P1(1:M)))','gauss1');% no of samples to plot
        [~,idx] = max(fitfun(w(1:M)));
        rollRateFreq = cat(1,rollRateFreq,w(idx));
        rollRateAmplitudes = cat(1,rollRateAmplitudes, sqrt(2)*rms(X - mean(X)));
        
        X = roll(i:i+N)';
        n = 2^nextpow2(N);      % pad the input with trailing zeros
        Y = fft(X,n);         % compute the FFT
        P2 = abs(Y/N);          % double-sided spectrum of each signal
        P1 = P2(:,1:n/2+1);     % single-sided spectrum of each signal
        P1(:,2:end-1) = 2*P1(:,2:end-1);
        P1(1:7) = zeros(1,7);
        f = 0:(fs/n):(fs/2-fs/n); w = 2*pi*f;  % frequency vectors       
        fitfun = fit(w(1:M)',(P1(1:M)/max(P1(1:M)))','gauss1');% no of samples to plot
        [~,idx] = max(fitfun(w(1:M)));
        rollFreq = cat(1,rollFreq,w(idx));
        rollAmplitudes = cat(1,rollAmplitudes, sqrt(2)*rms(X - mean(X)));
        
        X = pitch(i:i+N)';
        n = 2^nextpow2(N);      % pad the input with trailing zeros
        Y = fft(X,n);         % compute the FFT
        P2 = abs(Y/N);          % double-sided spectrum of each signal
        P1 = P2(:,1:n/2+1);     % single-sided spectrum of each signal
        P1(:,2:end-1) = 2*P1(:,2:end-1);
        P1(1:7) = zeros(1,7);
        %plot(w(1:M)',(P1(1:M)/max(P1(1:M)))')
        f = 0:(fs/n):(fs/2-fs/n); w = 2*pi*f;  % frequency vectors       
        fitfun = fit(w(1:M)',(P1(1:M)/max(P1(1:M)))','gauss1');% no of samples to plot
        [~,idx] = max(fitfun(w(1:M)));
        pitchFreq = cat(1,pitchFreq,w(idx));
        pitchAmplitudes = cat(1,pitchAmplitudes, sqrt(2)*rms(X - mean(X)));
        
        
        X = HeaveEst(i:i+N)';
        n = 2^nextpow2(N);      % pad the input with trailing zeros
        Y = fft(X,n);         % compute the FFT
        P2 = abs(Y/N);          % double-sided spectrum of each signal
        P1 = P2(:,1:n/2+1);     % single-sided spectrum of each signal
        P1(:,2:end-1) = 2*P1(:,2:end-1);
        P1(1:7) = zeros(1,7);
        f = 0:(fs/n):(fs/2-fs/n); w = 2*pi*f;  % frequency vectors  
        %plot(w(1:M)',(P1(1:M)/max(P1(1:M)))')
        fitfun = fit(w(1:M)',(P1(1:M)/max(P1(1:M)))','gauss1');% no of samples to plot
        [~,idx] = max(fitfun(w(1:M)));
        heavefreq = cat(1,heavefreq,w(idx));
        heaveamplitude = cat(1,heaveamplitude, sqrt(2)*rms(X - mean(X)));
        
        
        Sogs = cat(1,Sogs, mean(SpeedOverGround(i:i+N)));
    end
end

%%
figure;scatter(pitchRatesFreq,Sogs);
figure;scatter(rollRateFreq,Sogs);
figure;scatter(pitchRateAmplitudes,Sogs);
figure;scatter(rollRateAmplitudes,Sogs);
figure;scatter(heavefreq,Sogs);
figure;scatter(heaveamplitude,Sogs);
figure;scatter(rollFreq,Sogs);
figure;scatter(rollAmplitudes,Sogs);
figure;scatter(pitchFreq,Sogs);
figure;scatter(pitchAmplitudes,Sogs);

%%
MdlInput = [pitchRatesFreq rollRateFreq pitchRateAmplitudes ...
    rollRateAmplitudes heavefreq heaveamplitude rollFreq ...
    rollAmplitudes pitchFreq pitchAmplitudes];
%MdlInput = [pitchamplitudes rollamplitudes  heaveamplitude];
%%
Mdl1 = fitrgp(MdlInput, (Sogs), 'KernelFunction', 'matern52');
%%
PlotGaus(Sogs, Mdl1,MdlInput,'Vg')

%%
CorrData = [Sogs MdlInput];
corrCoefs = corrcoef(CorrData);
figure;
% yvalues = {'Vg','pitchamplitudes',...
%     'rollamplitudes',  'heaveamplitude'};
% xvalues = {'Vg','pitchamplitudes',...
%     'rollamplitudes',  'heaveamplitude'};
yvalues = {'Vg','pitchRatesFreq','rollRateFreq','pitchRateAmplitudes',...
    'rollRateAmplitudes', 'heavefreq', 'heaveamplitude','rollFreq', ...
    'rollAmplitudes', 'pitchFreq', 'pitchAmplitudes'};
xvalues = {'Vg','pitchRatesFreq','rollRateFreq','pitchRateAmplitudes',...
    'rollRateAmplitudes', 'heavefreq', 'heaveamplitude','rollFreq', ...
    'rollAmplitudes', 'pitchFreq', 'pitchAmplitudes'};
h = heatmap(xvalues,yvalues,corrCoefs);
h.Title = 'Correlation Matrix';
%% TODo
% Heave
% Gaussian and test

%%
% 
% 
% 
% 
% a = [];b= [];speed = [];
% for i = 31: 60 : length(GpsFix.height) - 31
%     b = [b meanfreq(GpsFix.height(i-30:i+30)-mean(GpsFix.height(i-30:i+30)))];
% end
% for i = 31: 60 : length(GpsFix.height) - 31
%     a = [a rms(GpsFix.height(i-30:i+30)-mean(GpsFix.height(i-30:i+30)))];
% end
% for i = 31: 60 : length(GpsFix.height) - 31
%     speed = [speed mean(GpsFix.sog(i-30:i+30))];
% enduntitled fit 1
% 
% hpFilt = designfilt('highpassfir','StopbandFrequency',0.2, ...
%     'PassbandFrequency',0.3,'PassbandRipple',0.5, ...
%     'StopbandAttenuation',65,'DesignMethod','kaiserwin', 'SampleRate',126);
% 
% lpFilt = designfilt('lowpassfir','PassbandFrequency',3, ...
%     'StopbandFrequency',4,'PassbandRipple',0.5, ...
%     'StopbandAttenuation',65,'DesignMethod','kaiserwin','SampleRate',126);
% fvtool(lpFilt)
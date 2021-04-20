%% Clear Workspace
clc; clearvars;close all;
addpath '..'
addpath '../AnalysisFiles'



pitchRatesFreq = [];
pitchRateAmplitudes = [];
pitchFreq = [];
pitchAmplitudes = [];
rollRateFreq = [];
rollRateAmplitudes = [];
rollFreq = [];
rollAmplitudes = [];
heavefreq = [];
heaveAmplitude = [];
heaveDotAmpl = [];
heaveDotfreq = [];
heaveeAmpl = [];
heaveeFreq = [];
Sogs = [];

for inputData = 1 : 8
%% Load Data
    if inputData == 1 
        path = '../Mausund200701_181204';
    elseif inputData == 2 
        path = '../Mausund200703_062402';
    elseif inputData == 3 
        path = '../Mausund200703_080820';
    elseif inputData == 4
        path = '../Mausund200701_221241';
    elseif inputData == 5
        path = '../Mausund200703_132548'; 
    elseif inputData == 6
        path = '../Mausund200703_215938';
    elseif inputData == 7
        path = '../Mausund200706_154608';
    elseif inputData == 8
        path = '../Mausund200709_53748';
    end 
    addpath(path)
    load AngularVelocity.mat;
    load Acceleration.mat;
    load GpsFix.mat;
    load Heave.mat;
    load EulerAngles.mat;
    rmpath(path)
    
    %% 
    %plot(Acceleration.x)

    % 123 Hz 

    fs = 123;                   % sampling frequency
    minutes = 3;                % 30 minutes data
    h = 1/fs;                   % sampling time
    N = minutes*60*fs;                 
    t = (0:N-1)*h;   

    timestamp = AngularVelocity.timestamp - AngularVelocity.timestamp(1);
    
    GpsFix.timestamp = GpsFix.timestamp - GpsFix.timestamp(1);
    Heave.timestamp = Heave.timestamp(Heave.src_ent==39);
    Heave.timestamp = Heave.timestamp - Heave.timestamp(1);
    
    SpeedOverGround = interp1(GpsFix.timestamp, GpsFix.sog, timestamp);
    HeaveEsttimator = Heave.value(Heave.src_ent==39);
    HeaveEsttimator = smooth(HeaveEsttimator);
    HeaveEst = interp1(Heave.timestamp-Heave.timestamp(1), HeaveEsttimator, timestamp);
    roll = smooth(EulerAngles.phi);
    roll = interp1(EulerAngles.timestamp -EulerAngles.timestamp(1) , roll, timestamp);
    pitch = smooth(EulerAngles.theta);
    pitch = interp1(EulerAngles.timestamp - EulerAngles.timestamp(1), pitch, timestamp);
  
    
    M = 4000;
    AngularVelocity.y = smooth(AngularVelocity.y);
    AngularVelocity.x = smooth(AngularVelocity.x);
    for i = 10*N:N:length(AngularVelocity.y) - 10*N

        X = (double(AngularVelocity.y(i:3:i+N)));
        time = timestamp(i:3:i+N);
        [pks,locs] = findpeaks(X,time,'MinPeakProminence',0.04 ,'MinPeakDistance',0.1);
        avg_periods_from_peaks = mean(diff(locs));

        avg_freq_hz = 1./avg_periods_from_peaks;
        avg_freq_radians_per_second = 2*pi*avg_freq_hz;
        
        
        pitchRatesFreq = cat(1,pitchRatesFreq,avg_freq_radians_per_second);
        pitchRateAmplitudes = cat(1,pitchRateAmplitudes, 2*sqrt(2)*rms(X - mean(X)));
        
        
        
        
%         figpeaks = figure;
%         plot(time, X);
%         hold on
%         plot(locs, pks, 'pg', 'MarkerFaceColor','red','MarkerSize',15)
%         xlabel('Time [s]')
%         ylabel('Heave Displacement [m]', 'Interpreter','latex', 'FontSize', 11)
%         title(['Heave Displacement and Peaks'])
%         grid on
        
%         n = 2^nextpow2(N);      % pad the input with trailing zeros
%         Y = fft(X,n);           % compute the FFT
%         P2 = abs(Y/N);          % double-sided spectrum of each signal
%         P1 = P2(:,1:n/2+1);     % single-sided spectrum of each signal
%         P1(:,2:end-1) = 2*P1(:,2:end-1);
%         %P1(1) = 0;
%         f = 0:(fs/n):(fs/2-fs/n); w = 2*pi*f;  % frequency vectors      
%         fitfun = fit(w(1:M)',(P1(1:M)/max(P1(1:M)))','gauss1');
%         [~,idx] = max(fitfun(w(1:M)));

        
        X = (double(AngularVelocity.x(i:3:i+N)));
        time = timestamp(i:3:i+N);
        [pks,locs] = findpeaks(X,time,'MinPeakProminence',0.04 ,'MinPeakDistance',1);
        avg_periods_from_peaks = mean(diff(locs));

        avg_freq_hz = 1./avg_periods_from_peaks;
        avg_freq_radians_per_second = 2*pi*avg_freq_hz;
        rollRateFreq = cat(1,rollRateFreq,avg_freq_radians_per_second);
        rollRateAmplitudes = cat(1,rollRateAmplitudes, 2*sqrt(2)*rms(X - mean(X)));
        

        
        X = smooth(double(roll(i:3:i+N)));
        time = timestamp(i:3:i+N);
        [pks,locs] = findpeaks(X,time,'MinPeakProminence',0.01 ,'MinPeakDistance',1);
        avg_periods_from_peaks = mean(diff(locs));

        avg_freq_hz = 1./avg_periods_from_peaks;
        avg_freq_radians_per_second = 2*pi*avg_freq_hz;
        rollFreq = cat(1,rollFreq,avg_freq_radians_per_second);
        rollAmplitudes = cat(1,rollAmplitudes, sqrt(2)*rms(X - mean(X)));
        
%         figure;
%         plot(time, X);
%         hold on
%        
%         plot(locs, pks, 'pg', 'MarkerFaceColor','red','MarkerSize',15)
%         xlabel('Time [s]')
%         ylabel('Heave Displacement [m]', 'Interpreter','latex', 'FontSize', 11)
%         grid
%         pause;
%         close;

        
        X = smooth(pitch(i:3:i+N));
        time = timestamp(i:3:i+N);
        [pks,locs] = findpeaks(X,time,'MinPeakProminence',0.01 ,'MinPeakDistance',1);
        avg_periods_from_peaks = mean(diff(locs));

        avg_freq_hz = 1./avg_periods_from_peaks;
        avg_freq_radians_per_second = 2*pi*avg_freq_hz;
        pitchFreq = cat(1,pitchFreq,avg_freq_radians_per_second);
        pitchAmplitudes = cat(1,pitchAmplitudes, sqrt(2)*rms(X - mean(X)));
%                 figure;
%         plot(time, X);
%         hold on
%        
%         plot(locs, pks, 'pg', 'MarkerFaceColor','red','MarkerSize',15)
%         xlabel('Time [s]')
%         ylabel('Heave Displacement [m]', 'Interpreter','latex', 'FontSize', 11)
%         grid
%         pause;
%         close;
%         figure(1)
%         plot(w(1:M)',(P1(1:M)/max(P1(1:M)))')
%         hold on
%         plot(fitfun, w(1:M)',(P1(1:M)/max(P1(1:M)))')
%         pause(0.1)
%         hold off
%         %stop;
        
        X = (double(HeaveEst(i:3:i+N)));
        
        time = timestamp(i:3:i+N);
        [pks,locs] = findpeaks(X,time,'MinPeakProminence',0.1,'MinPeakDistance',3);
        avg_periods_from_peaks = mean(diff(locs));
%         plot(time,X)
%         pause(0.1)
        avg_freq_hz = 1./avg_periods_from_peaks;
        avg_freq_radians_per_second = 2*pi*avg_freq_hz;
        heavefreq = cat(1,heavefreq,avg_freq_radians_per_second);
        heaveAmplitude = cat(1,heaveAmplitude, 2*sqrt(2)*rms(X - mean(X)));
%         figure;
%         plot(time, X);
%         hold on
%        
%         plot(locs, pks, 'pg', 'MarkerFaceColor','red','MarkerSize',15)
%         xlabel('Time [s]')
%         ylabel('Heave Displacement [m]', 'Interpreter','latex', 'FontSize', 11)
%         grid
%         pause;
%         close;
        
%         figure;
%         plot(time, X);
%         hold on
%        
%         plot(locs, pks, 'pg', 'MarkerFaceColor','red','MarkerSize',15)
%         xlabel('Time [s]')
%         ylabel('Heave Displacement [m]', 'Interpreter','latex', 'FontSize', 11)
%         grid
%         pause;
%         close;

%         figure(1)
%         plot(w(1:M)',(P1(1:M)/max(P1(1:M)))')
%         hold on
%         plot(fitfun, w(1:M)',(P1(1:M)/max(P1(1:M)))')
%         pause(0.1)+N)');
        
%         hold off
%         %stop;
        Sogs = cat(1,Sogs, mean(SpeedOverGround(i:i+N)));
        
%         X = smooth(double(heaveDot(i:3:i+N)));
%         time = timestamp(i:3:i+N);
%         [pks,locs] = findpeaks(X,time,'MinPeakProminence',0.01,'MinPeakHeight',0.01 ,'MinPeakDistance',1);
%         avg_periods_from_peaks = mean(diff(locs));
% 
%         avg_freq_hz = 1./avg_periods_from_peaks;
%         avg_freq_radians_per_second = 2*pi*avg_freq_hz;
%         heaveDotfreq = cat(1,heaveDotfreq,avg_freq_radians_per_second);
%         heaveDotAmpl = cat(1,heaveDotAmpl, sqrt(2)*rms(X - mean(X)));
%         %figure(12);plot(w(1:M)',(P1(1:M)/max(P1(1:M)))')
%         
%         X = smooth(double(heavee(i:3:i+N)));
%         time = timestamp(i:3:i+N);
%         [pks,locs] = findpeaks(X,time,'MinPeakProminence',0.01,'MinPeakHeight',0.01 ,'MinPeakDistance',1);
%         avg_periods_from_peaks = mean(diff(locs));
% 
%         avg_freq_hz = 1./avg_periods_from_peaks;
%         avg_freq_radians_per_second = 2*pi*avg_freq_hz;
%         heaveeFreq = cat(1,heaveeFreq,avg_freq_radians_per_second);
%         heaveeAmpl = cat(1,heaveeAmpl, sqrt(2)*rms(X - mean(X)));
    end
end

%%
figure;scatter(pitchRatesFreq,Sogs);
figure;scatter(rollRateFreq,Sogs);
figure;scatter(pitchRateAmplitudes,Sogs);
figure;scatter(rollRateAmplitudes,Sogs);
figure;scatter(heavefreq,Sogs);
figure;scatter(heaveAmplitude,Sogs);
figure;scatter(rollFreq,Sogs);
figure;scatter(rollAmplitudes,Sogs);
figure;scatter(pitchFreq,Sogs);
figure;scatter(pitchAmplitudes,Sogs);
% figure;scatter(heaveDotfreq,Sogs);
% figure;scatter(heaveDotAmpl,Sogs);
% figure;scatter(heaveeFreq,Sogs);
% figure;scatter(heaveeAmpl,Sogs);



%%
%MdlInput = [  pitchRateAmplitudes ...
%     rollRateAmplitudes  ...
%     rollAmplitudes  pitchAmplitudes  ...
%     heavefreq  heaveamplitude];
%MdlInput = [heavefreq heaveamplitude];
% MdlInput = [pitchRatesFreq rollRateFreq pitchRateAmplitudes rollRateAmplitudes heavefreq heaveAmplitude ...
%     rollFreq  rollAmplitudes pitchFreq pitchAmplitudes];
MdlInput = [pitchRateAmplitudes rollRateAmplitudes heavefreq heaveAmplitude ...
      rollAmplitudes  pitchAmplitudes];
%%
Mdl1 = fitrgp(MdlInput, Sogs, 'KernelFunction', 'matern52','BasisFunction', 'linear','Sigma', 0.072985);%, ...
%     'OptimizeHyperparameters' ,'auto', 'HyperparameterOptimizationOptions',...
%     struct('AcquisitionFunctionName','expected-improvement-plus'));
%%
PlotGaus(Sogs, Mdl1,MdlInput,'Vg')
%%
X_ML = MdlInput';
[MyNet, performance, e, tr] = neuralNet(X_ML,Sogs', 12);
figure; plotregression(Sogs', MyNet(X_ML));

% %%
% CorrData = [Sogs MdlInput];
% corrCoefs = corrcoef(CorrData);
% figure;
% % yvalues = {'Vg','pitchamplitudes',...
% %     'rollamplitudes',  'heaveamplitude'};
% % xvalues = {'Vg','pitchamplitudes',...
% %     'rollamplitudes',  'heaveamplitude'};
% yvalues = {'Vg','pitchRatesFreq','rollRateFreq','pitchRateAmplitudes',...
%     'rollRateAmplitudes', 'heavefreq', 'heaveamplitude', 'rollFreq','rollAmplitudes',...
%     'pitchFreq', 'pitchAmplitudes'};
% xvalues = {'Vg','pitchRatesFreq','rollRateFreq','pitchRateAmplitudes',...
%     'rollRateAmplitudes', ' heavefreq', 'heaveamplitude', 'rollFreq','rollAmplitudes',...
%     'pitchFreq', 'pitchAmplitudes'};
% h = heatmap(xvalues,yvalues,corrCoefs);
% h.Title = 'Correlation Matrix';
%%
MotionToSpeedTest;



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
%% Clear Workspace
clc; clearvars; close all;
addpath '..'
addpath '../AnalysisFiles'


%% Constants
fs = 2;
n = 6*60*fs;
pad = 15*60*fs;

%% Data to save
lon_data = [];
lat_data = [];
sog_data = [];
CurrentSpeed_data = [];
Currentdir_data = [];
WaveDir_data = [];
WindDir_data = [];
WindSpeed_data = [];
WaveSize_data = [];
waveHz_data = [];

test_lon_data = [];
test_lat_data = [];
test_sog_data = [];
test_CurrentSpeed_data = [];
test_Currentdir_data = [];
test_WaveDir_data = [];
test_WindDir_data = [];
test_WindSpeed_data = [];
test_WaveSize_data = [];
test_waveHz_data = [];

for dataRun = 1:13
    %% Data Loading
    disp('Loading Weather data data')
    if dataRun == 1
        path = '../Mausund200701_181204';
        load '../Weather/weatherData_2020-7-1_2020-7-2';
        load '../Weather/currentweatherData_2020-7-1_2020-7-3';
        
    elseif dataRun == 2
        path = '../Mausund200701_221241';
    elseif dataRun == 3
        path = '../Mausund200703_062402';
        load '../Weather/weatherData_2020-7-3_2020-7-4';
        load '../Weather/currentweatherData_2020-7-3_2020-7-4';       
    elseif dataRun == 4
        path = '../Mausund200703_080820';   
    elseif dataRun == 5
        path = '../Mausund200703_215938';
    elseif dataRun == 6
        path = '../Mausund200705_120030';
        load '../Weather/weatherData_2020-7-5_2020-7-5';
        load '../Weather/currentweatherData_2020-7-5_2020-7-5' 
    elseif dataRun == 7
        path = '../Mausund200706_154608';
        load '../Weather/weatherData_2020-7-6_2020-7-6';
        load '../Weather/currentweatherData_2020-7-6_2020-7-6' 
    elseif dataRun == 8
        path = '../Mausund200709_53748';
        load '../Weather/weatherData_2020-7-9_2020-7-9';
        load '../Weather/currentweatherData_2020-7-9_2020-7-9' 
    elseif dataRun == 9
        path = '../Frohavet210312_100504';
        load '../Weather/weatherData_2021-3-12_2021-3-13';
        load '../Weather/currentweatherData_2021-3-12_2021-3-13'; 
    elseif dataRun == 10
        path = '../Frohavet210312_172709';
    elseif dataRun == 11
        path = '../Frohavet210313_083148';
        load '../Weather/weatherData_2021-3-13_2021-3-13';
        load '../Weather/currentweatherData_2021-3-13_2021-3-13'; 
        
    elseif dataRun == 12
        path = '../Frohavet210314_125306';
        load '../Weather/weatherData_2021-3-14_2021-3-15';
        load '../Weather/currentweatherData_2021-3-14_2021-3-15';
    elseif dataRun == 13
        path = '../Frohavet210314_204158';  
    else
        path = '../Mausund200703_132548';
        load '../Weather/weatherData_2020-7-3_2020-7-4';
        load '../Weather/currentweatherData_2020-7-3_2020-7-4';   
    end
    disp('done')
    disp('Loading AutoNaut data')
    addpath(path);
    load 'RelativeWind.mat';
    load 'AbsoluteWind.mat';
    load 'GpsFix.mat';
    load 'Heave.mat';
    load 'EulerAngles';
    rmpath(path);
    disp('done')
    % RelativeWind,GpsFix,Heave,EulerAngles,AbsoluteWind


    %% Interpolation and gathering correct data
    RelativeWind.angle = interp1(RelativeWind.timestamp, ssa(RelativeWind.angle,'deg' ),GpsFix.timestamp);
    RelativeWind.speed = interp1(RelativeWind.timestamp, RelativeWind.speed,GpsFix.timestamp);
    RelativeWind.timestamp = interp1(RelativeWind.timestamp, RelativeWind.timestamp,GpsFix.timestamp);
    
    AbsoluteWind.dir = interp1(AbsoluteWind.timestamp, ssa(AbsoluteWind.dir,'deg' ),GpsFix.timestamp);
    AbsoluteWind.speed = interp1(AbsoluteWind.timestamp, AbsoluteWind.speed,GpsFix.timestamp);
    AbsoluteWind.timestamp = interp1(AbsoluteWind.timestamp, AbsoluteWind.timestamp,GpsFix.timestamp);
    
    if Heave.src_ent==45
        HeaveValue = smooth(Heave.value(Heave.src_ent==45));
        HeaveTime = Heave.timestamp(Heave.src_ent==45);
    else
        HeaveValue = smooth(Heave.value(Heave.src_ent==46));
        HeaveTime = Heave.timestamp(Heave.src_ent==46);
    end
    
    %% Set off more space 
    idx_offset = length(lon_data);
    test_idx_offset = length(test_lon_data);
    newPadding = length(find(mod(1:length(pad: n : length(GpsFix.sog) - pad),5) ~= 0));
    test_newPadding = length(find(mod(1:length(pad: n : length(GpsFix.sog) - pad),5) == 0));
    
    lon_data = cat(1, lon_data, zeros(newPadding, 1));
    lat_data = cat(1, lat_data, zeros(newPadding, 1));
    sog_data = cat(1, sog_data, zeros(newPadding, 1));
    CurrentSpeed_data = cat(1, CurrentSpeed_data, zeros(newPadding, 1));
    Currentdir_data = cat(1, Currentdir_data, zeros(newPadding, 1));
    WaveDir_data = cat(1, WaveDir_data, zeros(newPadding, 1));
    WindDir_data = cat(1, WindDir_data, zeros(newPadding, 1));
    WindSpeed_data = cat(1, WindSpeed_data, zeros(newPadding, 1));
    WaveSize_data = cat(1, WaveSize_data, zeros(newPadding, 1));
    waveHz_data = cat(1, waveHz_data, zeros(newPadding, 1));
    
    test_lon_data = cat(1, test_lon_data, zeros(test_newPadding, 1));
    test_lat_data = cat(1, test_lat_data, zeros(test_newPadding, 1));
    test_sog_data = cat(1, test_sog_data, zeros(test_newPadding, 1));
    test_CurrentSpeed_data = cat(1, test_CurrentSpeed_data, zeros(test_newPadding, 1));
    test_Currentdir_data = cat(1, test_Currentdir_data, zeros(test_newPadding, 1));
    test_WaveDir_data = cat(1, test_WaveDir_data, zeros(test_newPadding, 1));
    test_WindDir_data = cat(1, test_WindDir_data, zeros(test_newPadding, 1));
    test_WindSpeed_data = cat(1, test_WindSpeed_data, zeros(test_newPadding, 1));
    test_WaveSize_data = cat(1, test_WaveSize_data, zeros(test_newPadding, 1));
    test_waveHz_data = cat(1, test_waveHz_data, zeros(test_newPadding, 1));
    
    if  strcmp(path,'../Mausund200701_221241/') % Remove end of this mission
        offff = 100;
    else
        offff = 15;
    end

    %% Main loop
    idx = 1;
    test_idx = 1;
    testingdata = 1;
    disp('Start Going through data')
    for i= pad: n : length(GpsFix.sog) - (120*offff)
        curr_hour = floor(double(GpsFix.utc_time(i))/3600) ...
            + 24*(double(GpsFix.utc_day(i)-GpsFix.utc_day(1)));

        % Latidtude and longitude position of the vessel
        lat = mean(rad2deg(GpsFix.lat(i:i+n)));
        lon = mean(rad2deg(GpsFix.lon(i:i+n)));

        % Heading, Cog and Sog
        cog = rad2deg(ssa(mean(GpsFix.cog(i:i+n))));
        psi = rad2deg(ssa(mean(EulerAngles.psi(i:i+n))));
        sog = mean(GpsFix.sog(i:i+n));

        % Find wave direction in Forecast data
        error_map = sqrt((latitudeMapWave - lat).^2 + (longitudeMapWave - lon).^2);
        [x,y] = find(error_map == min(error_map, [], 'all'));

        % Find Current info in Current Forecast data
        error_map = sqrt((latitudeCurrentMap - lat).^2 + (longitudeCurrentMap - lon).^2);
        [xcurrent,ycurrent] = find(error_map == min(error_map, [], 'all'));
        currentNorthCur = currentNorth(xcurrent,ycurrent,curr_hour + 1);
        currentEastCur = currentEast(xcurrent,ycurrent,curr_hour + 1);
        
        Vc = [currentNorthCur; currentEastCur];
        VcDir = atan2d( Vc(2), Vc(1));

        % Wind data
        WindDir = ssa(mean(RelativeWind.angle(i:i+n)-180),'deg');
        WindSpeed = mean(AbsoluteWind.speed(i:i+n));

        % Wave approx 
        curHeave = HeaveValue(i:i+n);
        time = HeaveTime(i:i+n);
        [pks,locs] = findpeaks(curHeave,time,'MinPeakProminence',0.1 ,'MinPeakDistance',1);
        avg_periods_from_peaks = mean(diff(locs));
        avg_freq_hz = 1./avg_periods_from_peaks;
        avg_freq_radians_per_second = 2*pi*avg_freq_hz;
        
        % Save data
        if mod(testingdata,5) 
            sog_data(idx+idx_offset) = sog;
            lon_data(idx + idx_offset) = lon;
            lat_data(idx + idx_offset) = lat;
            WindDir_data(idx + idx_offset) = WindDir;
            WindSpeed_data(idx + idx_offset) = WindSpeed;
            CurrentSpeed_data(idx + idx_offset) = norm(Vc);
            Currentdir_data(idx + idx_offset) = ssa(psi - VcDir,'deg');
            WaveDir_data(idx + idx_offset) = ssa(psi - waveDir(x,y,curr_hour+1),'deg');
            waveHz_data(idx + idx_offset) = avg_freq_radians_per_second;
            WaveSize_data(idx + idx_offset) = 2*sqrt(2)*rms(curHeave-mean(curHeave));
            idx = idx + 1;
        else
            test_sog_data(test_idx+test_idx_offset) = sog;
            test_lon_data(test_idx + test_idx_offset) = lon;
            test_lat_data(test_idx + test_idx_offset) = lat;
            test_WindDir_data(test_idx + test_idx_offset) = WindDir;
            test_WindSpeed_data(test_idx + test_idx_offset) = WindSpeed;
            test_CurrentSpeed_data(test_idx + test_idx_offset) = norm(Vc);
            test_Currentdir_data(test_idx + test_idx_offset) = ssa(psi - VcDir,'deg');
            test_WaveDir_data(test_idx + test_idx_offset) = ssa(psi - waveDir(x,y,curr_hour+1),'deg');
            test_waveHz_data(test_idx + test_idx_offset) = avg_freq_radians_per_second;
            test_WaveSize_data(test_idx + test_idx_offset) = 2*sqrt(2)*rms(curHeave-mean(curHeave));
            test_idx = test_idx + 1;
        end
        testingdata = testingdata + 1;
    end
    disp('done')
end


%% Plots for analysis
figure;scatter(CurrentSpeed_data, sog_data);ylabel 'Sog';xlabel 'current speed';
figure;scatter(Currentdir_data, sog_data);ylabel 'Sog';xlabel 'current dir';
figure;scatter(CurrentSpeed_data.*cos(deg2rad(Currentdir_data)), sog_data);ylabel 'Sog';xlabel 'Current speed surge';
figure;scatter(WaveDir_data, sog_data); ylabel 'Sog';xlabel 'wave dir';
figure;scatter(WindDir_data, sog_data); ylabel 'Sog';xlabel 'wind dir';
figure;scatter(WindSpeed_data, sog_data); ylabel 'Sog';xlabel 'wind speed';
figure;scatter(WindSpeed_data.*cos(deg2rad(WindDir_data)), sog_data);ylabel 'Sog';xlabel 'Wind speed surge';
figure;scatter(WaveSize_data, sog_data); ylabel 'Sog';xlabel 'wave size';
figure;scatter(waveHz_data, sog_data); ylabel 'Sog';xlabel 'wave Hz';

%% Correlation matrix
X = [WaveSize_data waveHz_data  abs(cos(deg2rad(WaveDir_data)))  ...
   WindSpeed_data.*cos(deg2rad(WindDir_data)) CurrentSpeed_data.*cos(deg2rad(Currentdir_data)) ...
   ones(length(sog_data),1)];
CorrData = [sog_data X(:, 1:end-1)];
corrCoefs = corrcoef(CorrData);
PlotHeat(corrCoefs,'Vg');

%% Fit linear Regression model
w1 = (X'*X)\(X'*sog_data);

%% Fit Gaussian Regression model
X_gauss = [CurrentSpeed_data Currentdir_data WaveDir_data WindDir_data ...
    WindSpeed_data WaveSize_data waveHz_data];
Mdl1 = fitrgp(X_gauss, sog_data, 'KernelFunction', 'ardmatern52','BasisFunction', 'linear','Sigma', 0.037317, ...% );%, ...
     'OptimizeHyperparameters' ,'auto', 'HyperparameterOptimizationOptions',...
     struct('AcquisitionFunctionName','expected-improvement-plus'));
%%
X_test = [test_WaveSize_data test_waveHz_data  abs(cos(deg2rad(test_WaveDir_data)))  ...
   test_WindSpeed_data.*cos(deg2rad(test_WindDir_data)) test_CurrentSpeed_data.*cos(deg2rad(test_Currentdir_data)) ...
   ones(length(test_sog_data),1)];


X_gauss_test = [test_CurrentSpeed_data test_Currentdir_data ...
    test_WaveDir_data test_WindDir_data test_WindSpeed_data ...
    test_WaveSize_data test_waveHz_data];


%% Bootstrap
My_indices = randperm(length(X_gauss));
somen = floor(length(X_gauss)/5);
X_gauss2 = X_gauss(My_indices(1:somen), :);
X_gauss3 = X_gauss(My_indices(somen + 1: 2*somen), :);
X_gauss4 = X_gauss(My_indices(2*somen + 1: 3*somen), :);
X_gauss5 = X_gauss(My_indices(3*somen + 1: 4*somen), :);
X_gauss6 = X_gauss(My_indices(4*somen + 1: end), :);
sog_data2 = sog_data(My_indices(1:somen), :);
sog_data3 = sog_data(My_indices(somen + 1: 2*somen), :);
sog_data4 = sog_data(My_indices(2*somen + 1: 3*somen), :);
sog_data5 = sog_data(My_indices(3*somen + 1: 4*somen), :);
sog_data6 = sog_data(My_indices(4*somen + 1: end), :);
% Mdl2 = fitrgp(X_gauss2, sog_data2, 'KernelFunction', 'ardmatern52','BasisFunction', 'linear','Sigma',  0.039647, ...% );%, ...
%      'OptimizeHyperparameters' ,'auto', 'HyperparameterOptimizationOptions',...
%      struct('AcquisitionFunctionName','expected-improvement-plus'));
% Mdl3 = fitrgp(X_gauss3, sog_data3, 'KernelFunction', 'ardmatern52','BasisFunction', 'linear','Sigma', 0.062457, ...% );%, ...
%      'OptimizeHyperparameters' ,'auto', 'HyperparameterOptimizationOptions',...
%      struct('AcquisitionFunctionName','expected-improvement-plus'));
% Mdl4 = fitrgp(X_gauss4, sog_data4, 'KernelFunction', 'ardmatern52','BasisFunction', 'linear','Sigma', 0.062545, ...% );%, ...
%      'OptimizeHyperparameters' ,'auto', 'HyperparameterOptimizationOptions',...
%      struct('AcquisitionFunctionName','expected-improvement-plus'));
% Mdl5 = fitrgp(X_gauss5, sog_data5, 'KernelFunction', 'ardmatern52','BasisFunction', 'linear','Sigma', 0.04399, ...% );%, ...
%      'OptimizeHyperparameters' ,'auto', 'HyperparameterOptimizationOptions',...
%      struct('AcquisitionFunctionName','expected-improvement-plus'));
% Mdl6 = fitrgp(X_gauss6, sog_data6, 'KernelFunction', 'ardmatern52','BasisFunction', 'linear','Sigma', 0.051829, ...% );%, ...
%      'OptimizeHyperparameters' ,'auto', 'HyperparameterOptimizationOptions',...
%      struct('AcquisitionFunctionName','expected-improvement-plus'));
% 
% Mdl7 = fitrgp(X_gauss, sog_data, 'KernelFunction', 'ardmatern32','BasisFunction', 'linear','Sigma', 0.051829, ...% );%, ...
%      'OptimizeHyperparameters' ,'auto', 'HyperparameterOptimizationOptions',...
%      struct('AcquisitionFunctionName','expected-improvement-plus'));
% Mdl8 = fitrgp(X_gauss, sog_data, 'KernelFunction', 'ardsquaredexponential','BasisFunction', 'linear','Sigma', 0.051829, ...% );%, ...
%      'OptimizeHyperparameters' ,'auto', 'HyperparameterOptimizationOptions',...
%      struct('AcquisitionFunctionName','expected-improvement-plus'));

%% Plot Model Data
PlotLinear(test_sog_data,w1,X_test,'Vg')
PlotGaus(test_sog_data, Mdl1, X_gauss_test,'Vg')

[pred,~, yci] = predict(Mdl1, X_gauss_test);
figure;
plot(1:length(test_sog_data),test_sog_data,'r.');
hold on
plot(1:length(test_sog_data),pred);
plot(1:length(test_sog_data),(yci(:,1)),'k:');
plot(1:length(test_sog_data),(yci(:,2)),'k:');
xlabel('x');
ylabel('y');

out2 = predict(Mdl2, X_gauss_test);
out3 = predict(Mdl3, X_gauss_test);
out4 = predict(Mdl4, X_gauss_test);
out5 = predict(Mdl5, X_gauss_test);
out6 = predict(Mdl6, X_gauss_test);
bootOut = (out2 + out3 + out4 + out5 + out6)/5;
figure;
plot(1:length(test_sog_data),test_sog_data,'r.');
hold on
plot(1:length(test_sog_data),bootOut);

%% Run Testscript
testForecastBased;
testSensorBased;

%% Clearing and ending script
rmpath '../AnalysisFiles'
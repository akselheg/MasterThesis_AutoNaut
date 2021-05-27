%% Clear Workspace
clc; clearvars; close all;
addpath '..'
addpath './AnalysisFiles'


%% Constants
fs = 2;
n = 3*60*fs;
pad = 5*60*fs;
first = true;

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

for dataRun = 1:4
    %% Data Loading
    disp('Loading Weather data data')
    if dataRun == 1
        path = './Trondheim082546';
        load './Weather/currentweatherData_2020-2-20_2020-2-20.mat';
        load './Weather/weatherData_2020-2-20_2020-2-20.mat';
        Heavemask = 42;
    elseif dataRun == 2
        path = 'Trondheim094058';
    elseif dataRun == 3
        path = 'Trondheim101916';
    elseif dataRun == 4
        path = 'Trondheim101916';
    elseif dataRun == 5
        path = 'Trondheim101916';
        load './Weather/currentweatherData_2020-5-28_2020-5-29.mat';
        load './Weather/weatherData_2020-5-28_2020-5-28.mat';
        Heavemask = 39;
    end
    disp('done')
        
    disp('Loading AutoNaut data')
    addpath(path);
    load 'RelativeWind.mat';
    %load 'AbsoluteWind.mat';
    load 'GpsFix.mat';
    load 'Heave.mat';
    load 'EulerAngles';
    rmpath(path);
    disp('done')
    % RelativeWind,GpsFix,Heave,EulerAngles,AbsoluteWind


    %% Interpolation and correcting data
    RelativeWind.angle = interp1(RelativeWind.timestamp, ssa(RelativeWind.angle,'deg' ),GpsFix.timestamp);
    RelativeWind.speed = interp1(RelativeWind.timestamp, RelativeWind.speed,GpsFix.timestamp);
    RelativeWind.timestamp = interp1(RelativeWind.timestamp, RelativeWind.timestamp,GpsFix.timestamp);
    
%     AbsoluteWind.dir = interp1(AbsoluteWind.timestamp, ssa(AbsoluteWind.dir,'deg' ),GpsFix.timestamp);
%     AbsoluteWind.speed = interp1(AbsoluteWind.timestamp, AbsoluteWind.speed,GpsFix.timestamp);
%     AbsoluteWind.timestamp = interp1(AbsoluteWind.timestamp, AbsoluteWind.timestamp,GpsFix.timestamp);

    HeaveValue = smooth(Heave.value(Heave.src_ent==Heavemask));
    HeaveTime = Heave.timestamp(Heave.src_ent==Heavemask);

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


    %% Main loop
    idx = 1;
    test_idx = 1;
    testingdata = 1;
    disp('Start Going through data')
    for i= pad: n : length(GpsFix.sog) - pad
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
        WindSpeed = mean(RelativeWind.speed(i:i+n));

        % Wave approx 
        X = HeaveValue(i:i+n);
        time = HeaveTime(i:i+n);
        [pks,locs] = findpeaks(X,time,'MinPeakProminence',0.01,'MinPeakHeight',0.01 ,'MinPeakDistance',1);
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
            WaveSize_data(idx + idx_offset) = 2*sqrt(2)*rms(X-mean(X));
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
            test_WaveSize_data(test_idx + test_idx_offset) = 2*sqrt(2)*rms(X-mean(X));
            test_idx = test_idx + 1;
        end
        testingdata = testingdata + 1;
    end
    disp('done')
end


%% Plotting for analysis
figure;scatter(CurrentSpeed_data, sog_data);ylabel 'Sog';xlabel 'current speed';
figure;scatter(Currentdir_data, sog_data);ylabel 'Sog';xlabel 'current dir';
figure;scatter(CurrentSpeed_data.*cos(deg2rad(Currentdir_data)), sog_data);ylabel 'Sog';xlabel 'Current speed surge';
figure;scatter(WaveDir_data, sog_data); ylabel 'Sog';xlabel 'wave dir';
figure;scatter(WindDir_data, sog_data); ylabel 'Sog';xlabel 'wind dir';
figure;scatter(WindSpeed_data, sog_data); ylabel 'Sog';xlabel 'wind speed';
figure;scatter(WindSpeed_data.*cos(deg2rad(WindDir_data)), sog_data);ylabel 'Sog';xlabel 'Wind speed surge';
figure;scatter(WaveSize_data, sog_data); ylabel 'Sog';xlabel 'wave size';
figure;scatter(waveHz_data, sog_data); ylabel 'Sog';xlabel 'wave Hz';

% Heatmap
X = [CurrentSpeed_data.*cos(deg2rad(Currentdir_data)) WaveDir_data  ...
   WindSpeed_data.*cos(deg2rad(WindDir_data))  WaveSize_data  waveHz_data ...
   ones(length(sog_data),1)];
CorrData = [sog_data X(:, 1:end-1)];
corrCoefs = corrcoef(CorrData);

figure;
labels = {'Vg','CurrentSpeedsurge','WaveDir',...
    'WindSpeedSurge', 'WaveSize', 'waveHz'};
h = heatmap(labels,labels,corrCoefs);
h.Title = 'Correlation Matrix';


w1 = (X'*X)\(X'*sog_data);

X_test = [test_CurrentSpeed_data.*cos(deg2rad(test_Currentdir_data)) test_WaveDir_data  ...
   test_WindSpeed_data.*cos(deg2rad(test_WindDir_data)) test_WaveSize_data  test_waveHz_data ...
   ones(length(test_sog_data),1)];
PlotLinear(test_sog_data,w1,X_test,'Vg')

%% Gaussian Regression model
X_gauss = [CurrentSpeed_data Currentdir_data WaveDir_data WindDir_data ...
    WindSpeed_data WaveSize_data waveHz_data];
Mdl1 = fitrgp(X_gauss, sog_data, 'KernelFunction', 'matern52');
X_gauss_test = [test_CurrentSpeed_data test_Currentdir_data ...
    test_WaveDir_data test_WindDir_data test_WindSpeed_data ...
    test_WaveSize_data test_waveHz_data];
PlotGaus(test_sog_data, Mdl1, X_gauss_test,'Vg');

%% Machine Learning model
X_ML = X_gauss';
[MyNet, performance, e, tr] = neuralNet(X_ML,sog_data', 10);
X_ML_test = X_gauss_test';
figure; plotregression(sog_data', MyNet(X_ML));
figure; plotregression(test_sog_data', MyNet(X_ML_test));

%% Run Testscript
%testNewAnalysis;

%% Clearing and ending script
rmpath './AnalysisFiles'
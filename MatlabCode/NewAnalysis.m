%% Clear Workspace
clc; clearvars; close all;

%% Constants
fs = 2;
n = 3*60*fs;
pad = 15*60*fs;
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

%% Data Loading

disp('Loading AutoNaut data')
path = 'Mausund200701_181204';
addpath(path);
    load 'AbsoluteWind.mat';
    load 'GpsFix.mat';
    load 'Heave.mat';
    load 'EulerAngles';
rmpath(path);
disp('done')
disp('Loading Weather data data')
load './Weather/weatherData_2020-7-1_2020-7-2';
load './Weather/currentweatherData_2020-7-1_2020-7-3';
disp('done')

%% Interpolation and correcting data
AbsoluteWind.dir = interp1(AbsoluteWind.timestamp, ssa(AbsoluteWind.dir,'deg' ),GpsFix.timestamp);
AbsoluteWind.speed = interp1(AbsoluteWind.timestamp, AbsoluteWind.speed,GpsFix.timestamp);
AbsoluteWind.timestamp = interp1(AbsoluteWind.timestamp, AbsoluteWind.timestamp,GpsFix.timestamp);

HeaveValue = Heave.value(Heave.src_ent==39);
HeaveTime = Heave.timestamp(Heave.src_ent==39);



%% Set off more space 
idx_offset = length(lon_data);
lon_data = cat(1, lon_data, zeros(length(pad: n : length(GpsFix.sog) - pad), 1));
lat_data = cat(1, lat_data, zeros(length(pad: n : length(GpsFix.sog) - pad), 1));
sog_data = cat(1, sog_data, zeros(length(pad: n : length(GpsFix.sog) - pad), 1));
CurrentSpeed_data = cat(1, CurrentSpeed_data, zeros(length(pad: n : length(GpsFix.sog) - pad), 1));
Currentdir_data = cat(1, Currentdir_data, zeros(length(pad: n : length(GpsFix.sog) - pad), 1));
WaveDir_data = cat(1, WaveDir_data, zeros(length(pad: n : length(GpsFix.sog) - pad), 1));
WindDir_data = cat(1, WindDir_data, zeros(length(pad: n : length(GpsFix.sog) - pad), 1));
WindSpeed_data = cat(1, WindSpeed_data, zeros(length(pad: n : length(GpsFix.sog) - pad), 1));
WaveSize_data = cat(1, WaveSize_data, zeros(length(pad: n : length(GpsFix.sog) - pad), 1));
waveHz_data = cat(1, waveHz_data, zeros(length(pad: n : length(GpsFix.sog) - pad), 1));


%% Main loop
idx = 1;
for i= pad: n : length(GpsFix.sog) - pad
    curr_hour = floor(double(GpsFix.utc_time(i))/3600) ...
        + 24*(double(GpsFix.utc_day(i)-GpsFix.utc_day(1)));
    
    % Latidtude and longitude position of the vessel
    lat = mean(rad2deg(GpsFix.lat(i:i+n)));
    lon = mean(rad2deg(GpsFix.lon(i:i+n)));
    lon_data(idx + idx_offset) = lon;
    lat_data(idx + idx_offset) = lat;
    

    % Heading, Cog and Sog
    cog = rad2deg(ssa(mean(GpsFix.cog(i:i+n))));
    psi = rad2deg(ssa(mean(EulerAngles.psi(i:i+n))));
    sog = mean(GpsFix.sog(i:i+n));
    sog_data(idx+idx_offset) = sog;

    % Find wave direction in Forecast data
    error_map = sqrt((latitudeMapWave - lat).^2 + (longitudeMapWave - lon).^2);
    [x,y] = find(error_map == min(error_map, [], 'all'));
    WaveDir_data(idx + idx_offset) = ssa(psi - waveDir(x,y,curr_hour+1),'deg');
    
    
    % Find Current info in Current Forecast data
    error_map = sqrt((latitudeCurrentMap - lat).^2 + (longitudeCurrentMap - lon).^2);
    [xcurrent,ycurrent] = find(error_map == min(error_map, [], 'all'));
    currentNorthCur = currentNorth(xcurrent,ycurrent,curr_hour + 1);
    currentEastCur = currentEast(xcurrent,ycurrent,curr_hour + 1);
    Vc = [currentNorthCur; currentEastCur];
    CurrentSpeed_data(idx + idx_offset) = norm(Vc);
    
    VcDir = atan2d( Vc(2), Vc(1));
    Currentdir_data(idx + idx_offset) = ssa(psi - VcDir,'deg');
    
    % Wind data
    WindDir_data(idx + idx_offset) = ssa(psi - mean(AbsoluteWind.dir(i:i+n)),'deg');
    WindSpeed_data(idx + idx_offset) = mean(AbsoluteWind.speed(i:i+n));
    
    % Wave approx 
    
    X = HeaveValue(i:i+n);
    time = HeaveTime(i:i+n);
    [pks,locs] = findpeaks(X,time,'MinPeakProminence',0.1,'MinPeakHeight',0.1 ,'MinPeakDistance',2);
    avg_periods_from_peaks = mean(diff(locs));

    avg_freq_hz = 1./avg_periods_from_peaks;
    avg_freq_radians_per_second = 2*pi*avg_freq_hz;
    
    waveHz_data(idx + idx_offset) = avg_freq_radians_per_second;
    WaveSize_data(idx + idx_offset) = sqrt(2)*rms(X - mean(X));
    
    
    % increase idx
    idx = idx + 1;
end

%% Plotting and Analysis
figure;scatter(CurrentSpeed_data, sog_data)
figure;scatter(Currentdir_data, sog_data)
figure;scatter(WaveDir_data, sog_data)
figure;scatter(WindDir_data, sog_data)
figure;scatter(WindSpeed_data, sog_data)
figure;scatter(WaveSize_data, sog_data)
figure;scatter(waveHz_data, sog_data)


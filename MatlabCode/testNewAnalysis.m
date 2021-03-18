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
disp('Loading Weather data data')
path = './Mausund200703_132548';
load './Weather/weatherData_2020-7-3_2020-7-4';
load './Weather/currentweatherData_2020-7-3_2020-7-4';   
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

%% Interpolation and correcting data
RelativeWind.angle = interp1(RelativeWind.timestamp, ssa(RelativeWind.angle,'deg' ),GpsFix.timestamp);
RelativeWind.speed = interp1(RelativeWind.timestamp, RelativeWind.speed,GpsFix.timestamp);
RelativeWind.timestamp = interp1(RelativeWind.timestamp, RelativeWind.timestamp,GpsFix.timestamp);
AbsoluteWind.dir = interp1(AbsoluteWind.timestamp, ssa(AbsoluteWind.dir,'deg' ),GpsFix.timestamp);
AbsoluteWind.speed = interp1(AbsoluteWind.timestamp, AbsoluteWind.speed,GpsFix.timestamp);
AbsoluteWind.timestamp = interp1(AbsoluteWind.timestamp, AbsoluteWind.timestamp,GpsFix.timestamp);


HeaveValue = smooth(Heave.value(Heave.src_ent==39));
HeaveTime = Heave.timestamp(Heave.src_ent==39);

%% Set off more space 
idx_offset = length(lon_data);
newPadding = length(1:length(pad: n : length(GpsFix.sog) - pad));

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
%% Main loop
idx = 1;
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
    WindSpeed = mean(AbsoluteWind.speed(i:i+n));

    % Wave approx 
    X = HeaveValue(i:i+n);
    time = HeaveTime(i:i+n);
    [pks,locs] = findpeaks(X,time,'MinPeakProminence',0.1,'MinPeakHeight',0.1 ,'MinPeakDistance',1);
    avg_periods_from_peaks = mean(diff(locs));
    avg_freq_hz = 1./avg_periods_from_peaks;
    avg_freq_radians_per_second = 2*pi*avg_freq_hz;

    % Save Data
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

end
disp('done')
%% Linear Regression model
X = [CurrentSpeed_data.*cos(deg2rad(Currentdir_data)) WaveDir_data  ...
   WindSpeed_data.*cos(deg2rad(WindDir_data))  WaveSize_data  waveHz_data ...
   ones(length(sog_data),1)];
PlotLinear(sog_data,w1,X,'Vg')

%% Gaussian Regression model
X_gauss = [CurrentSpeed_data Currentdir_data WaveDir_data WindDir_data ...
    WindSpeed_data WaveSize_data waveHz_data];
PlotGaus(sog_data, Mdl1, X_gauss,'Vg')

%% Machine Learning model
X_ML = X_gauss';
figure; plotregression(sog_data', MyNet(X_ML));
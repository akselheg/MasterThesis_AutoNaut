%% Constants
fs = 2;
n = 3*60*fs;
pad = 15*60*fs;
first = true;

%% Data to save
new_lon_data = [];
new_lat_data = [];
new_sog_data = [];
new_CurrentSpeed_data = [];
new_Currentdir_data = [];
new_WaveDir_data = [];
new_WindDir_data = [];
new_WindSpeed_data = [];
new_WaveSize_data = [];
new_waveHz_data = [];

%% Data Loading
disp('Loading Weather data data')
path = '../Mausund200703_132548';
load '../Weather/weatherData_2020-7-3_2020-7-4';
load '../Weather/currentweatherData_2020-7-3_2020-7-4';   
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


HeaveValue = smooth(Heave.value(Heave.src_ent==46));
HeaveTime = Heave.timestamp(Heave.src_ent==46);

%% Set off more space 
idx_offset = length(new_lon_data);
newPadding = length(1:length(pad: n : length(GpsFix.sog) - pad));

new_lon_data = cat(1, new_lon_data, zeros(newPadding, 1));
new_lat_data = cat(1, new_lat_data, zeros(newPadding, 1));
new_sog_data = cat(1, new_sog_data, zeros(newPadding, 1));
new_CurrentSpeed_data = cat(1, new_CurrentSpeed_data, zeros(newPadding, 1));
new_Currentdir_data = cat(1, new_Currentdir_data, zeros(newPadding, 1));
new_WaveDir_data = cat(1, new_WaveDir_data, zeros(newPadding, 1));
new_WindDir_data = cat(1, new_WindDir_data, zeros(newPadding, 1));
new_WindSpeed_data = cat(1, new_WindSpeed_data, zeros(newPadding, 1));
new_WaveSize_data = cat(1, new_WaveSize_data, zeros(newPadding, 1));
new_waveHz_data = cat(1, new_waveHz_data, zeros(newPadding, 1));
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
    curHeave = HeaveValue(i:i+n);
    time = HeaveTime(i:i+n);
    [pks,locs] = findpeaks(curHeave,time,'MinPeakProminence',0.1,'MinPeakDistance',1);
    avg_periods_from_peaks = mean(diff(locs));
    avg_freq_hz = 1./avg_periods_from_peaks;
    avg_freq_radians_per_second = 2*pi*avg_freq_hz;

    % Save Data
    new_sog_data(idx+idx_offset) = sog;
    new_lon_data(idx + idx_offset) = lon;
    new_lat_data(idx + idx_offset) = lat;
    new_WindDir_data(idx + idx_offset) = WindDir;
    new_WindSpeed_data(idx + idx_offset) = WindSpeed;
    new_CurrentSpeed_data(idx + idx_offset) = norm(Vc);
    new_Currentdir_data(idx + idx_offset) = ssa(psi - VcDir,'deg');
    new_WaveDir_data(idx + idx_offset) = ssa(psi - waveDir(x,y,curr_hour+1),'deg');
    new_waveHz_data(idx + idx_offset) = avg_freq_radians_per_second;
    new_WaveSize_data(idx + idx_offset) = 2*sqrt(2)*rms(curHeave-mean(curHeave));
    idx = idx + 1;

end
disp('done')
% addpath '../AnalysisFiles'
%% Linear Regression model
new_X = [new_CurrentSpeed_data.*cos(deg2rad(new_Currentdir_data)) new_WaveDir_data  ...
   new_WindSpeed_data.*cos(deg2rad(new_WindDir_data))  new_WaveSize_data  new_waveHz_data ...
   ones(length(new_sog_data),1)];
PlotLinear(new_sog_data,w1,new_X,'Vg')

%% Gaussian Regression model
new_X_gauss = [new_CurrentSpeed_data new_Currentdir_data new_WaveDir_data new_WindDir_data ...
    new_WindSpeed_data new_WaveSize_data new_waveHz_data];
PlotGaus(new_sog_data, Mdl1, new_X_gauss,'Vg')
%%
[pred,~, yci] = predict(Mdl1, new_X_gauss);
figure;
plot(1:length(new_sog_data),new_sog_data,'r.');
hold on
plot(1:length(new_sog_data),pred);
plot(1:length(new_sog_data),(yci(:,1)),'k:');
plot(1:length(new_sog_data),(yci(:,2)),'k:');
xlabel('x');
ylabel('y');
%%
out2 = predict(Mdl2, new_X_gauss);
out3 = predict(Mdl3, new_X_gauss);
out4 = predict(Mdl4, new_X_gauss);
out5 = predict(Mdl5, new_X_gauss);
out6 = predict(Mdl6, new_X_gauss);
bootOut = (out2 + out3 + out4 + out5 + out6)/5;
figure;
plot(1:length(new_sog_data),new_sog_data,'r.');
hold on
plot(1:length(new_sog_data),bootOut);


%% Machine Learning model
new_X_ML = new_X_gauss';
figure; plotregression(new_sog_data', MyNet(new_X_ML));
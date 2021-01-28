%% Clear workspace
clc; clearvars; close all;

%% Data That can be downladed from neptus that are relevant
% AbsoluteWind,Depth,DesiredHeading,DesiredPath,DesiredSpeed,DesiredZ,GpsFix,RelativeWind,RemoteSensorInfo,EstimatedState,EulerAngles
% GpsFix,RelativeWind,EulerAngles,Heave
% Data to be saved for plots
addpath './AnalysisFiles/';
lon_data = [];
lat_data = [];
sog_data = [];
ForecastWaveSize_data = [];
messuredRelWindDir_data = [];
messuredRelWindSpeed_data = [];
relWaveDir_data = [];
ForecastWaveFreq_data = [];
ForcastWindSpeed_data = [];
CurrentSpeed_data = [];
Vr_data = [];
test_sog_data = [];
test_ForecastWaveSize_data = [];
test_ForecastWaveFreq_data = [];
test_relWaveDir_data = [];
test_ForcastWindSpeed_data = [];
test_CurrentSpeed_data = [];
test_Vr_data = [];

xmax = 0; ymax = 0; ymin = inf; xmin = inf;
avrager = 6*60; % average over x min
count = 1;
for i = 2:10
    disp('Loading new data')
    %% load data
    if i == 1
        path = './Mausund200701_181204/';
        addpath(path);
        gpsFix = load('GpsFix.mat');
        RelativeWind = load('RelativeWind.mat');
        EulerAngles = load('EulerAngles.mat');
        rmpath(path)
        load('./Weather/weatherData_2020-7-1_2020-7-2.mat') % Must be downloaded locally
        load('./Weather/currentweatherData_2020-7-1_2020-7-3.mat') % Must be downloaded locally
        disp('Done loading data')
    end
    if i == 2
        path = './Mausund200701_221241/';
        addpath(path);
        gpsFix = load('GpsFix.mat');
        RelativeWind = load('RelativeWind.mat');
        EulerAngles = load('EulerAngles.mat');
        load('./Weather/weatherData_2020-7-1_2020-7-2.mat')
        load('./Weather/currentweatherData_2020-7-1_2020-7-3.mat')
        rmpath(path)
        disp('Done loading data')
    end
    if i == 3
        path = './Mausund200703_080820/';
        addpath(path);
        gpsFix = load('GpsFix.mat');
        RelativeWind = load('RelativeWind.mat');
        EulerAngles = load('EulerAngles.mat');
        load('./Weather/weatherData_2020-7-3_2020-7-4.mat')
        load('./Weather/currentweatherData_2020-7-3_2020-7-4.mat')
        rmpath(path)
        disp('Done loading data')
    end
    if i == 4
        path = './Mausund200703_132548/';
        addpath(path);
        gpsFix = load('GpsFix.mat');
        RelativeWind = load('RelativeWind.mat');
        EulerAngles = load('EulerAngles.mat');
        load('./Weather/weatherData_2020-7-3_2020-7-4.mat')
        load('./Weather/currentweatherData_2020-7-3_2020-7-4.mat')
        rmpath(path)
        disp('Done loading data')
    end
    if i == 5
        path = './Mausund200705_120030/';
        addpath(path);
        gpsFix = load('GpsFix.mat');
        RelativeWind = load('RelativeWind.mat');
        EulerAngles = load('EulerAngles.mat');
        rmpath(path)
        load('./Weather/weatherData_2020-7-5_2020-7-5.mat')
        load('./Weather/currentweatherData_2020-7-5_2020-7-5.mat')
        disp('Done loading data')
    end
    if i == 6
        path = './Mausund200706_154608/';
        addpath(path);
        gpsFix = load('GpsFix.mat');
        RelativeWind = load('RelativeWind.mat');
        EulerAngles = load('EulerAngles.mat');
        rmpath(path)
        load('./Weather/weatherData_2020-7-6_2020-7-6.mat')
        load('./Weather/currentweatherData_2020-7-6_2020-7-6.mat')
        disp('Done loading data')
    end
    if i == 7
        path = './Mausund200709_53748/';
        addpath(path);
        gpsFix = load('GpsFix.mat');
        RelativeWind = load('RelativeWind.mat');
        EulerAngles = load('EulerAngles.mat');
        rmpath(path)
        load('./Weather/weatherData_2020-7-9_2020-7-9.mat')
        load('./Weather/currentweatherData_2020-7-9_2020-7-9.mat')
        disp('Done loading data')
    end
    if i == 8
        path = './Mausund200710_092034/';
        addpath(path);
        gpsFix = load('GpsFix.mat');
        RelativeWind = load('RelativeWind.mat');
        EulerAngles = load('EulerAngles.mat');
        rmpath(path)
        load('./Weather/weatherData_2020-7-10_2020-7-10.mat')
        load('./Weather/currentweatherData_2020-7-10_2020-7-10.mat')
        disp('Done loading data')
    end
    if i == 9
        path = './Mausund200703_215938/';
        addpath(path);
        gpsFix = load('GpsFix.mat');
        gpsFix.GpsFix.sog = gpsFix.GpsFix.sog(1:54441);
        RelativeWind = load('RelativeWind.mat');
        EulerAngles = load('EulerAngles.mat');
        load('./Weather/weatherData_2020-7-3_2020-7-4.mat')
        load('./Weather/currentweatherData_2020-7-3_2020-7-4.mat')
        rmpath(path)
        disp('Done loading data')
    end
    if i == 10
        path = './Mausund200712_220202/';
        addpath(path);
        gpsFix = load('GpsFix.mat');
        gpsFix.GpsFix.sog = gpsFix.GpsFix.sog(1:70080);
        RelativeWind = load('RelativeWind.mat');
        EulerAngles = load('EulerAngles.mat');
        rmpath(path)
        load('./Weather/weatherData_2020-7-12_2020-7-12.mat')
        load('./Weather/currentweatherData_2020-7-12_2020-7-13.mat')
        disp('Done loading data')
    end
    %% Format and interpolations
    gps_data = gpsFix.GpsFix;
    windData = RelativeWind.RelativeWind;
    EulerAngles = EulerAngles.EulerAngles;
    EulerAngles.psi = ssa(EulerAngles.psi,'deg');
    messuredRelWindDir = interp1(windData.timestamp, ssa(windData.angle,'deg' ),gps_data.timestamp);
    messuredRelWindSpeed = interp1(windData.timestamp, windData.speed,gps_data.timestamp);
    old_hour = 100000;
    disp('Done formating')
    disp('Start running through data')
    %% run
    for m = (10*120) : 2*avrager: length(gps_data.sog) - (10*120)
        curr_hour = floor(double(gps_data.utc_time(m))/3600) ...
            + 24*(double(gps_data.utc_day(m)-gps_data.utc_day(1)));

        % Latidtude and longitude position of the vessel
        lat = mean(rad2deg(gps_data.lat(m-avrager:m+avrager)));
        lon = mean(rad2deg(gps_data.lon(m-avrager:m+avrager)));

        % Heading, Cog and Sog
        cog = rad2deg(mean(gps_data.cog(m-avrager:m+avrager)));
        psi = rad2deg(mean(EulerAngles.psi(m-avrager:m+avrager)));
        sog = mean(gps_data.sog(m-avrager:m+avrager));

        % Find position in wave data
        error_map = sqrt((latitudeMapWave - lat).^2 + (longitudeMapWave - lon).^2);
        [x,y] = find(error_map == min(error_map, [], 'all'));
        if x > xmax; xmax = x; end
        if y > ymax; ymax = y; end
        if x < xmin; xmin = x; end
        if y < ymin; ymin = y; end
        lon_data = cat(1,lon_data, lon);
        lat_data = cat(1,lat_data, lat);

        % Find position in Current data
        error_map = sqrt((latitudeCurrentMap - lat).^2 + (longitudeCurrentMap - lon).^2);
        [xcurrent,ycurrent] = find(error_map == min(error_map, [], 'all'));

        % Wind and wave directions and size at given time and postion
        curWaveDir = ssa(waveDir(x,y,curr_hour+1),'deg');
        curWindDir = ssa(windDir(x,y,curr_hour+1),'deg');
        ForcastWindSpeed = windSpeed(x,y,curr_hour + 1);

        % Wave frequency at given time and position           
        if waveHZ(x,y,curr_hour+1) < 0.1 
            disp([num2str(waveHZ(x,y,curr_hour+1)) num2str(lat) num2str(lon)])
        end

        % Current vector at given time and position
        currentNorthCur = currentNorth(xcurrent,ycurrent,curr_hour+1);
        currentEastCur = currentEast(xcurrent,ycurrent,curr_hour+1);
        Vc = [currentNorthCur; currentEastCur];

        % Velocity vector of the vessel
        Vg = [sog*cos(deg2rad(cog)); sog*sin(deg2rad(cog))];
        Vr = Vg - Vc;

        % Angle between velocity and current direction
        VcDir = atan2d( Vc(2), Vc(1));
        VrDir = atan2d( Vr(2), Vr(1));
        CurVsVelAnglre = VrDir- VcDir;

        % magnitude of the current
        currentSpeed = norm(Vc);
        VrSpeed = norm(Vr);

        % Messured wind speed and direction relative to the vessel
        curMessuredRelWindDir = mean(messuredRelWindDir(m-avrager:m+avrager));
        curMessuredRelWindSpeed = mean(messuredRelWindSpeed(m-avrager:m+avrager));
        relWaveDir = ssa(psi- curWaveDir - 180, 'deg');
        ForecastWaveFreq =  1/waveHZ(x,y,curr_hour+1);
        ForecastWaveSize = waveSize(x, y, curr_hour + 1);
        currentSurge = currentSpeed*cos(ssa(deg2rad(VcDir - psi)));

        % Save current data
        if mod(count,5)
            ForecastWaveFreq_data = cat(1,ForecastWaveFreq_data, ForecastWaveFreq);            
            if waveSize(x, y, curr_hour + 1) < 0.001
                disp('Some error in forecast data size. Using previous')
                ForecastWaveSize_data = cat(1, ForecastWaveSize_data, ForecastWaveSize_data(end));
            else
                ForecastWaveSize_data = cat(1, ForecastWaveSize_data, ForecastWaveSize);
            end
            messuredRelWindDir_data = cat(1, messuredRelWindDir_data, curMessuredRelWindDir);
            messuredRelWindSpeed_data = cat(1, messuredRelWindSpeed_data, curMessuredRelWindSpeed);
            Vr_data = cat(1, Vr_data,VrSpeed);
            sog_data = cat(1, sog_data,sog);
            ForcastWindSpeed_data = cat(1, ForcastWindSpeed_data, ForcastWindSpeed*cos(ssa(deg2rad(curWindDir - psi))));
            CurrentSpeed_data = cat(1, CurrentSpeed_data, currentSurge);
            relWaveDir_data = cat(1, relWaveDir_data, relWaveDir);
        else
            if waveSize(x, y, curr_hour + 1) < 0.001
               disp('Some error in forecast data size. Using previous')
               test_ForecastWaveSize_data = cat(1,test_ForecastWaveSize_data,test_ForecastWaveSize_data(end));
            else
               test_ForecastWaveSize_data = cat(1,test_ForecastWaveSize_data, ForecastWaveSize);
            end
            test_Vr_data = cat(1, test_Vr_data,VrSpeed);
            test_sog_data = cat(1, test_sog_data,sog);
            test_ForecastWaveFreq_data = cat(1,test_ForecastWaveFreq_data, ForecastWaveFreq);
            test_ForcastWindSpeed_data = cat(1,test_ForcastWindSpeed_data, ForcastWindSpeed*cos(ssa(deg2rad(curWindDir - psi))));
            test_CurrentSpeed_data = cat(1,test_CurrentSpeed_data, currentSurge);
            test_relWaveDir_data = cat(1,test_relWaveDir_data, relWaveDir);
        end
        if old_hour ~= curr_hour 
            str = sprintf('| Day: %d  | Hour: %d \t|', ...
                (floor(curr_hour/24)+1) + gps_data.utc_day(1)-1, (mod(curr_hour,24)));
            disp(str)
            old_hour = curr_hour;
        end
        count = count + 1;
    end
    disp('Run Success')
    figure(1)
    plot(lon_data,lat_data, 'b')
    hold on
    scatter(lon_data(1), lat_data(1), 'g')
    scatter(lon_data(end), lat_data(end), 'r')
    lon_data = [];
    lat_data = [];
    pause(0.01)
end
%%
quiver(longitudeMapWave,latitudeMapWave, windDir(:,:,2),windSpeed(:,:,2))
meanSog = mean(sog_data);

%% Fit Linear model
X = [ForecastWaveSize_data ForecastWaveFreq_data ((deg2rad(relWaveDir_data)))  ...
   ForcastWindSpeed_data CurrentSpeed_data ones(length(sog_data),1)];
w1 = (X'*X)\(X'*sog_data);
X_test = [test_ForecastWaveSize_data test_ForecastWaveFreq_data ((deg2rad(test_relWaveDir_data)))  ...
    test_ForcastWindSpeed_data test_CurrentSpeed_data ones(length(test_sog_data),1)]; 
CorrData1 = [[sog_data;test_sog_data] [X(:, 1:end-1);X_test(:, 1:end-1)]];
corrCoefs1 = corrcoef(CorrData1);
w2 = (X'*X)\(X'*Vr_data);
CorrData2 = [[Vr_data;test_Vr_data] [X(:, 1:end-1);X_test(:, 1:end-1)]];
corrCoefs2 = corrcoef(CorrData2);
%% Fit Gaussian model
rng('default')
Mdl1 = fitrgp(X(:,1:end-1), sog_data,'Basis','linear' ,...
       'FitMethod','exact','PredictMethod','exact', 'KernelFunction', 'matern52');
Mdl2 = fitrgp(X(:,1:end-1), Vr_data,'Basis','linear',...
       'FitMethod','exact','PredictMethod','exact', 'KernelFunction', 'matern52');


%% Test of models on test dataset
PlotLinear(test_sog_data, w1,X_test,'Vg')
PlotLinear(test_Vr_data, w2,X_test,'Vr')
PlotGaus(test_sog_data, Mdl1,X_test(:,1:end-1),'Vg')
PlotGaus(test_Vr_data, Mdl2,X_test(:,1:end-1),'Vr')

%% Test of models on training dataset
PlotLinear(sog_data,w1,X,'Vg')
PlotLinear(Vr_data,w2,X,'Vr')
PlotGaus(sog_data,Mdl1,X(:,1:end-1),'Vg')
PlotGaus(Vr_data,Mdl2,X(:,1:end-1),'Vr')

%% Plot correlation matrices
PlotHeat(corrCoefs1,'Vg')
PlotHeat(corrCoefs2,'Vr')

%% Relevant plots
disp('Plotting Data')
% speed as a func of relative wave angle of attack with wave size coloured 
Plot3Dim([sog_data; test_sog_data], [relWaveDir_data; test_relWaveDir_data],...
    [ForecastWaveSize_data; test_ForecastWaveSize_data], 1.2, 1.6, ...
    false, 'Wave Size', 'Wave angle of attack [deg]', 'Vg [m/s]')
% speed as a func of relative wave angle of attack with wave frequency coloured 
Plot3Dim([sog_data; test_sog_data], [relWaveDir_data; test_relWaveDir_data],...
    [ForecastWaveFreq_data; test_ForecastWaveFreq_data], 0.13, 0.16, ...
    false, 'Wave Hz', 'Wave angle of attack [deg]', 'Vg [m/s]')
% speed as a func of wave frequency with wave size coloured 
Plot3Dim([sog_data; test_sog_data], [ForecastWaveFreq_data; test_ForecastWaveFreq_data],...
    [ForecastWaveSize_data; test_ForecastWaveSize_data], 1.2, 1.6, ...
    true, 'Wave Size', 'Wave Frequency [Hz]', 'Vg [m/s]')
% speed as a func of wave size with wave frequency coloured 
Plot3Dim([sog_data; test_sog_data], [ForecastWaveSize_data; test_ForecastWaveSize_data], ...
    [ForecastWaveFreq_data; test_ForecastWaveFreq_data], 0.13, 0.16, ...
    true, 'Wave Hz', 'Wave Size [m]', 'Vg [m/s]')
% speed as a func of wave frequency
Plot2Dim([sog_data; test_sog_data],  [ForecastWaveFreq_data; test_ForecastWaveFreq_data], ...
    'Wave Frequency [Hz]', 'Vg [m/s]')
% speed as a func of wave size 
Plot2Dim([sog_data; test_sog_data], [ForecastWaveSize_data; test_ForecastWaveSize_data], ...
    'Wave Size [m]', 'Vg [m/s]')
% speed as a func of rellative wind direction with relative wind speed coloured 
Plot3Dim(sog_data, messuredRelWindDir_data, messuredRelWindSpeed_data, 3, 6, ...
    false, 'Wind Speed', 'Measured Relative wind direction [deg]', 'Vg [m/s]')
% speed as a func of current speed in surge direction of body 
Plot2Dim([sog_data; test_sog_data], [CurrentSpeed_data; test_CurrentSpeed_data], ...
    'Current Speed in Surge Direction [m/s]', 'Vg [m/s]')
% speed as a func of current speed in surge direction of body 
Plot2Dim([sog_data; test_sog_data],[ForcastWindSpeed_data; test_ForcastWindSpeed_data], ...
    'Wind Speed in Surge Direction [m/s]',  'Vg [m/s]')

%% Redo for testing on different dataset
% Data to be saved for plots
lon_data = [];
lat_data = [];
sog_data = [];
ForecastWaveSize_data = [];
messuredRelWindDir_data = [];
messuredRelWindSpeed_data = [];
relWaveDir_data = [];
ForecastWaveFreq_data = [];
ForcastWindSpeed_data = [];
CurrentDir_data = [];
CurrentSpeed_data = [];
Vr_data = [];
VcDir_data = [];

xmax = 0; ymax = 0; ymin = inf; xmin = inf;
avrager = 6*60; % average over x min
count = 1;
disp('Loading new data')
%% load data
path = './Mausund200701_181204/';
addpath(path);
gpsFix = load('GpsFix.mat');
RelativeWind = load('RelativeWind.mat');
EulerAngles = load('EulerAngles.mat');
rmpath(path)
load('./Weather/weatherData_2020-7-1_2020-7-2.mat') % Must be downloaded locally
load('./Weather/currentweatherData_2020-7-1_2020-7-3.mat') % Must be downloaded locally
disp('Done loading data')
%% Format and interpolations
gps_data = gpsFix.GpsFix;
windData = RelativeWind.RelativeWind;
EulerAngles = EulerAngles.EulerAngles;
EulerAngles.psi = ssa(EulerAngles.psi,'deg');
messuredRelWindDir = interp1(windData.timestamp, ssa(windData.angle,'deg' ),gps_data.timestamp);
messuredRelWindSpeed = interp1(windData.timestamp, windData.speed,gps_data.timestamp);
old_hour = 10000;
disp('Done formating')
disp('Start running through data')
%% run
for m = (10*120) : 2*avrager:length(gps_data.sog) - (10*120)
    curr_hour = floor(double(gps_data.utc_time(m))/3600) ...
        + 24*(double(gps_data.utc_day(m)-gps_data.utc_day(1)));

    % Latidtude and longitude position of the vessel
    lat = mean(rad2deg(gps_data.lat(m-avrager:m+avrager)));
    lon = mean(rad2deg(gps_data.lon(m-avrager:m+avrager)));

    % Heading, Cog and Sog
    cog = rad2deg(mean(gps_data.cog(m-avrager:m+avrager)));
    psi = rad2deg(mean(EulerAngles.psi(m-avrager:m+avrager)));
    sog = mean(gps_data.sog(m-avrager:m+avrager));

    % Find position in wave data
    error_map = sqrt((latitudeMapWave - lat).^2 + (longitudeMapWave - lon).^2);
    [x,y] = find(error_map == min(error_map, [], 'all'));
    if x > xmax; xmax = x; end
    if y > ymax; ymax = y; end
    if x < xmin; xmin = x; end
    if y < ymin; ymin = y; end
    lon_data = cat(1,lon_data, lon);
    lat_data = cat(1,lat_data, lat);

    % Find position in Current data
    error_map = sqrt((latitudeCurrentMap - lat).^2 + (longitudeCurrentMap - lon).^2);
    [xcurrent,ycurrent] = find(error_map == min(error_map, [], 'all'));

    % Wind and wave directions and size at given time and postion
    curWaveDir = ssa(waveDir(x,y,curr_hour+1),'deg');
    curWindDir = ssa(windDir(x,y,curr_hour+1),'deg');
    ForcastWindSpeed = windSpeed(x,y,curr_hour + 1);


    % Wave frequency at given time and position

    if waveHZ(x,y,curr_hour+1) < 0.1 
        disp([num2str(waveHZ(x,y,curr_hour+1)) num2str(lat) num2str(lon)])
    end
    % Current vector at given time and position
    currentNorthCur = currentNorth(xcurrent,ycurrent,curr_hour+1);
    currentEastCur = currentEast(xcurrent,ycurrent,curr_hour+1);
    Vc = [currentNorthCur; currentEastCur];

    % Velocity vector of the vessel
    Vg = [sog*cos(deg2rad(cog)); sog*sin(deg2rad(cog))];
    Vr = Vg - Vc;
    % Angle between velocity and current direction
    VcDir = atan2d( Vc(2), Vc(1));
    VrDir = atan2d( Vr(2), Vr(1));
    VcDir_data = cat(1,VcDir_data, VcDir);

    % magnitude of the current
    currentSpeed = norm(Vc);
    VrSpeed = norm(Vr);

    % Messured wind speed and direction relative to the vessel
    curMessuredRelWindDir = mean(messuredRelWindDir(m-avrager:m+avrager));
    curMessuredRelWindSpeed = mean(messuredRelWindSpeed(m-avrager:m+avrager));

    % Save current data
    ForecastWaveFreq_data = cat(1,ForecastWaveFreq_data, 1/waveHZ(x,y,curr_hour+1));            
    if waveSize(x, y, curr_hour + 1) < 0.001
        disp('Some error in forecast data size. Using previous')
        ForecastWaveSize_data = cat(1, ForecastWaveSize_data, ForecastWaveSize_data(end));
    else
        ForecastWaveSize_data = cat(1, ForecastWaveSize_data, waveSize(x, y, curr_hour + 1));
    end
    Vr_data = cat(1, Vr_data,VrSpeed);
    sog_data = cat(1, sog_data,sog);
    relWaveDir_data = cat(1, relWaveDir_data, ssa(psi- curWaveDir - 180, 'deg'));
    messuredRelWindDir_data = cat(1, messuredRelWindDir_data, curMessuredRelWindDir);
    messuredRelWindSpeed_data = cat(1, messuredRelWindSpeed_data, curMessuredRelWindSpeed);
    ForcastWindSpeed_data = cat(1, ForcastWindSpeed_data, ForcastWindSpeed*cos(ssa(deg2rad(curWindDir - psi))));
    CurrentDir_data = cat(1, CurrentDir_data, ssa(VcDir - psi,'deg'));
    CurrentSpeed_data = cat(1, CurrentSpeed_data, currentSpeed*cos(ssa(deg2rad(VcDir - psi))));

    if old_hour ~= curr_hour
        str = sprintf('| Day: %d  | Hour: %d \t|', ...
            (floor(curr_hour/24)+1) + gps_data.utc_day(1)-1, (mod(curr_hour,24)));
        disp(str)
        old_hour = curr_hour;
    end
    count = count + 1;

end
disp('Run Success')
%%
X_testt = [ForecastWaveSize_data ForecastWaveFreq_data ((deg2rad(relWaveDir_data)))...
   ForcastWindSpeed_data CurrentSpeed_data, ones(length(sog_data),1)];
%%
PlotLinear(sog_data, w1,X_testt,'Vg')
%%
PlotLinear(Vr_data, w2, X_testt,'Vr')
%%
PlotGaus(sog_data, Mdl1,X_testt(:,1:end-1),'Vg')
%%
PlotGaus(Vr_data, Mdl2,X_testt(:,1:end-1),'Vr')
%%
%%
% disp('Doing Neural')
% nninputs =  double([X; X_test])';
% test_nninputs =  double(X_testt)';
% nntargets = double([Y1; Y1_test])';
% performance = inf;
% for i = 1:10
%     [temp_net, temp_perform, temp_netError, temp_netTrainState] =...
%         neuralNet(nninputs,nntargets , 4);
%     if temp_perform < performance
%         net = temp_net; performance = t%%

%         netError = temp_netError; netTrainState = temp_netTrainState;
%     end
% end
% disp(['Performance: ' num2str(performance)])
% figure, ploterrhist(netError)
% out = net(test_nninputs);
% e = gsubtract(double(sog_data)',out);
% figure, ploterrhist(e)
% figure, plotregression(double(sog_data),out)
%%
disp('Script done')



%% Clear workspace
clc; clearvars; close all;

%rehash; rehash path;rehash toolbox; rehash toolboxcache;

%% Data That can be downladed from neptus that are relevant
% AbsoluteWind,Depth,DesiredHeading,DesiredPath,DesiredSpeed,DesiredZ,GpsFix,RelativeWind,RemoteSensorInfo,EstimatedState,EulerAngles
% GpsFix,RelativeWind,EulerAngles,Heave
% Data to be saved for plots
addpath '..';
addpath '../AnalysisFiles/';

%% Constant
fs = 2;
n = 6*60*fs;
pad = 15*60*fs;

%% Data to be saved
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
ForcastWindDir_data = [];
CurrentDir_data =[];
Vr_data = [];

test_sog_data = [];
test_ForecastWaveSize_data = [];
test_ForecastWaveFreq_data = [];
test_relWaveDir_data = [];
test_ForcastWindSpeed_data = [];
test_CurrentSpeed_data = [];
test_ForcastWindDir_data = [];
test_CurrentDir_data =[];
test_Vr_data = [];
actuation = [];
test_actuation = [];

xmax = 0; ymax = 0; ymin = inf; xmin = inf; % Probably not neccesary, but if vessel outside of data-range
count = 1; % For every fifth 
for i = 1:13
    disp('Loading new data')
    %% load data
    if i == 1
        path = '../Mausund200701_181204/';
        load('../Weather/weatherData_2020-7-1_2020-7-2.mat') % Must be downloaded locally
        load('../Weather/currentweatherData_2020-7-1_2020-7-3.mat') % Must be downloaded locally
    elseif i == 2
        path = '../Mausund200701_221241/';
    elseif i == 3
        path = '../Mausund200703_062402';
        load '../Weather/weatherData_2020-7-3_2020-7-4';
        load '../Weather/currentweatherData_2020-7-3_2020-7-4';      
    elseif i == 4 
        path = '../Mausund200703_080820/';
    elseif i == 5 
        path = '../Mausund200703_215938/';
    elseif i == 6 
        path = '../Mausund200705_120030/';
        load('../Weather/weatherData_2020-7-5_2020-7-5.mat')
        load('../Weather/currentweatherData_2020-7-5_2020-7-5.mat')
    elseif i == 7 
        path = '../Mausund200706_154608/';
        load('../Weather/weatherData_2020-7-6_2020-7-6.mat')
        load('../Weather/currentweatherData_2020-7-6_2020-7-6.mat')
    elseif i == 8 
        path = '../Mausund200709_53748/';
        load('../Weather/weatherData_2020-7-9_2020-7-9.mat')
        load('../Weather/currentweatherData_2020-7-9_2020-7-9.mat')
    elseif i == 9
        path = '../Frohavet210312_100504';
        load '../Weather/weatherData_2021-3-12_2021-3-13';
        load '../Weather/currentweatherData_2021-3-12_2021-3-13'; 
    elseif i == 10
        path = '../Frohavet210312_172709';
    elseif i == 11
        path = '../Frohavet210313_083148';
        load '../Weather/weatherData_2021-3-13_2021-3-13';
        load '../Weather/currentweatherData_2021-3-13_2021-3-13'; 
        
    elseif i == 12
        path = '../Frohavet210314_125306';
        load '../Weather/weatherData_2021-3-14_2021-3-15';
        load '../Weather/currentweatherData_2021-3-14_2021-3-15';
    elseif i == 13
        path = '../Frohavet210314_204158';  
    end
    addpath(path);
    gpsFix = load('GpsFix.mat');
    RelativeWind = load('RelativeWind.mat');
    EulerAngles = load('EulerAngles.mat');
    rmpath(path)
    %% Format and interpolations
    disp('Formating')
    gps_data = gpsFix.GpsFix;
    windData = RelativeWind.RelativeWind;
    EulerAngles = EulerAngles.EulerAngles;
    EulerAngles.psi = ssa(EulerAngles.psi,'deg');
    messuredRelWindDir = interp1(windData.timestamp, ssa(windData.angle,'deg' ),gps_data.timestamp);
    messuredRelWindSpeed = interp1(windData.timestamp, windData.speed,gps_data.timestamp);
    old_hour = 100000;
    disp('Done formating')
    
    if  strcmp(path,'../Mausund200701_221241/') % Remove end of this mission
        offff = 100;
    else
        offff = 15;
    end
    %% run
    disp('Start running through data')
    for m = pad : n: length(gps_data.sog) - (offff*120)
        curr_hour = floor(double(gps_data.utc_time(m))/3600) ...
            + 24*(double(gps_data.utc_day(m)-gps_data.utc_day(1)));

        % Latidtude and longitude position of the vessel
        lat = mean(rad2deg(gps_data.lat(m : m + n)));
        lon = mean(rad2deg(gps_data.lon(m : m + n)));

        % Heading, Cog and Sog
        cog = rad2deg(mean(gps_data.cog(m : m + n)));
        psi = rad2deg(mean(EulerAngles.psi(m : m + n)));
        sog = mean(gps_data.sog(m : m + n));

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
        CurVsVelAnglre = ssa(VrDir- VcDir, 'deg');

        % magnitude of the current
        currentSpeed = norm(Vc);
        VrSpeed = norm(Vr);

        % Messured wind speed and direction relative to the vessel
        curMessuredRelWindDir = mean(messuredRelWindDir(m : m + n));
        curMessuredRelWindSpeed = mean(messuredRelWindSpeed(m : m + n));
        relWaveDir = ssa(psi - curWaveDir - 180, 'deg');
        ForecastWaveFreq =  2*pi/waveHZ(x,y,curr_hour+1);
        ForecastEncWaveFreq = abs(ForecastWaveFreq -((ForecastWaveFreq^2*sog*cos(deg2rad(relWaveDir))/9.81)));
        ForecastWaveSize = waveSize(x, y, curr_hour + 1);
        currentSurge = currentSpeed*cos(ssa(deg2rad(VcDir - psi)));
        windSurge = ForcastWindSpeed*cos(ssa(deg2rad(curWindDir - psi)));
        
        % Save current data
        if mod(count,5)
            ForecastWaveFreq_data = cat(1,ForecastWaveFreq_data, ForecastEncWaveFreq);            
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
            ForcastWindSpeed_data = cat(1, ForcastWindSpeed_data, ForcastWindSpeed);
            ForcastWindDir_data = cat(1, ForcastWindDir_data, ssa(curWindDir - psi,'deg'));
            CurrentSpeed_data = cat(1, CurrentSpeed_data, currentSpeed);
            CurrentDir_data = cat(1, CurrentDir_data, ssa(VcDir - psi,'deg'));
            relWaveDir_data = cat(1, relWaveDir_data, relWaveDir);

        else
            test_ForecastWaveSize_data = cat(1,test_ForecastWaveSize_data, ForecastWaveSize);
            test_Vr_data = cat(1, test_Vr_data,VrSpeed);
            test_sog_data = cat(1, test_sog_data,sog);
            test_ForecastWaveFreq_data = cat(1,test_ForecastWaveFreq_data, ForecastEncWaveFreq);
            
            test_ForcastWindSpeed_data = cat(1, test_ForcastWindSpeed_data, ForcastWindSpeed);
            test_ForcastWindDir_data = cat(1, test_ForcastWindDir_data, ssa(curWindDir - psi,'deg'));
            test_CurrentSpeed_data = cat(1, test_CurrentSpeed_data, currentSpeed);
            test_CurrentDir_data = cat(1, test_CurrentDir_data, ssa(VcDir - psi,'deg'));
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
    
    % Plot for visualising mission
    figure(1)
    geoplot(lat_data,lon_data, 'b')
    hold on
    geoscatter(lat_data(1), lon_data(1), 'g')
    geoscatter(lat_data(end), lon_data(end), 'r')
    lon_data = [];
    lat_data = [];
    pause(0.01)
end

%% Fit Linear model
X_linear = [ForecastWaveSize_data ForecastWaveFreq_data abs(cos(deg2rad(relWaveDir_data))) ...
    ForcastWindSpeed_data.*cos(deg2rad(ForcastWindDir_data))  CurrentSpeed_data.*cos(deg2rad(CurrentDir_data)) ...
    ones(length(sog_data),1)];
w1 = (X_linear'*X_linear)\(X_linear'*sog_data);
X_linear_test = [test_ForecastWaveSize_data test_ForecastWaveFreq_data abs(cos(deg2rad(test_relWaveDir_data))) ...
    test_ForcastWindSpeed_data.*cos(deg2rad(test_ForcastWindDir_data))  test_CurrentSpeed_data.*cos(deg2rad(test_CurrentDir_data)) ...
    ones(length(test_sog_data),1)]; 

%% Correlation matrix
CorrData1 = [[sog_data;test_sog_data] [X_linear(:, 1:end-1);X_linear_test(:, 1:end-1)]];
corrCoefs1 = corrcoef(CorrData1);
PlotHeat(corrCoefs1,'Vg')

%% Gaussian Process Regression model Inputs
X_Gauss = [CurrentSpeed_data CurrentDir_data relWaveDir_data ForcastWindDir_data ...
    ForcastWindSpeed_data ForecastWaveSize_data ForecastWaveFreq_data];
X_Gauss_test = [test_CurrentSpeed_data test_CurrentDir_data test_relWaveDir_data test_ForcastWindDir_data ...
    test_ForcastWindSpeed_data test_ForecastWaveSize_data test_ForecastWaveFreq_data];
%% Alteration for Bootstrap aggregation
My_indices = randperm(length(X_Gauss));
somen = floor(length(X_Gauss)/5);
X_gauss2 = X_Gauss(My_indices(1:somen), :);
X_gauss3 = X_Gauss(My_indices(somen + 1: 2*somen), :);
X_gauss4 = X_Gauss(My_indices(2*somen + 1: 3*somen), :);
X_gauss5 = X_Gauss(My_indices(3*somen + 1: 4*somen), :);
X_gauss6 = X_Gauss(My_indices(4*somen + 1: end), :);
sog_data2 = sog_data(My_indices(1:somen), :);
sog_data3 = sog_data(My_indices(somen + 1: 2*somen), :);
sog_data4 = sog_data(My_indices(2*somen + 1: 3*somen), :);
sog_data5 = sog_data(My_indices(3*somen + 1: 4*somen), :);
sog_data6 = sog_data(My_indices(4*somen + 1: end), :);

%% Learn all models
Mdl1 = fitrgp(X_Gauss, sog_data, 'KernelFunction', 'ardmatern52','BasisFunction', 'linear','Sigma', 0.042154, ...
    'OptimizeHyperparameters' ,'auto', 'HyperparameterOptimizationOptions',...
    struct('AcquisitionFunctionName','expected-improvement-plus'));
%%
Mdl2 = fitrgp(X_gauss2, sog_data2, 'KernelFunction', 'ardmatern52','BasisFunction', 'linear','Sigma',  0.039647, ...% );%, ...
     'OptimizeHyperparameters' ,'auto', 'HyperparameterOptimizationOptions',...
     struct('AcquisitionFunctionName','expected-improvement-plus'));
Mdl3 = fitrgp(X_gauss3, sog_data3, 'KernelFunction', 'ardmatern52','BasisFunction', 'linear','Sigma', 0.062457, ...% );%, ...
     'OptimizeHyperparameters' ,'auto', 'HyperparameterOptimizationOptions',...
     struct('AcquisitionFunctionName','expected-improvement-plus'));
Mdl4 = fitrgp(X_gauss4, sog_data4, 'KernelFunction', 'ardmatern52','BasisFunction', 'linear','Sigma', 0.062545, ...% );%, ...
     'OptimizeHyperparameters' ,'auto', 'HyperparameterOptimizationOptions',...
     struct('AcquisitionFunctionName','expected-improvement-plus'));
Mdl5 = fitrgp(X_gauss5, sog_data5, 'KernelFunction', 'ardmatern52','BasisFunction', 'linear','Sigma', 0.04399, ...% );%, ...
     'OptimizeHyperparameters' ,'auto', 'HyperparameterOptimizationOptions',...
     struct('AcquisitionFunctionName','expected-improvement-plus'));
Mdl6 = fitrgp(X_gauss6, sog_data6, 'KernelFunction', 'ardmatern52','BasisFunction', 'linear','Sigma', 0.051829, ...% );%, ...
     'OptimizeHyperparameters' ,'auto', 'HyperparameterOptimizationOptions',...
     struct('AcquisitionFunctionName','expected-improvement-plus'));
 %%
 Mdl7 = fitrgp(X_Gauss, sog_data, 'KernelFunction', 'ardmatern32','BasisFunction', 'linear','Sigma', 0.051829, ...% );%, ...
     'OptimizeHyperparameters' ,'auto', 'HyperparameterOptimizationOptions',...
     struct('AcquisitionFunctionName','expected-improvement-plus'));
 Mdl8 = fitrgp(X_Gauss, sog_data, 'KernelFunction', 'ardsquaredexponential','BasisFunction', 'linear','Sigma', 0.051829, ...% );%, ...
     'OptimizeHyperparameters' ,'auto', 'HyperparameterOptimizationOptions',...
     struct('AcquisitionFunctionName','expected-improvement-plus'));

%% Test of models on test dataset
PlotLinear(test_sog_data, w1,X_linear_test,'Vg')
PlotGaus(test_sog_data, Mdl1,X_Gauss_test,'Vg')

%% Test of models on training dataset
PlotLinear(sog_data,w1,X_linear,'Vg')
PlotGaus(sog_data, Mdl1,X_Gauss,'Vg')

%% Relevant plots
disp('Plotting Data')
% speed as a func of relative wave angle of attack with wave size coloured 
% Plot3Dim([sog_data; test_sog_data], [relWaveDir_data; test_relWaveDir_data],...
%     [ForecastWaveSize_data; test_ForecastWaveSize_data], 1.2, 1.6, ...
%     false, 'Vg [m/s]','Wave angle of attack [deg]',  'Wave Size')
% % speed as a func of relative wave angle of attack with wave frequency coloured 
% Plot3Dim([sog_data; test_sog_data], [relWaveDir_data; test_relWaveDir_data],...
%     [ForecastWaveFreq_data; test_ForecastWaveFreq_data], 6.5, 8, ...
%     false, 'Vg [m/s]', 'Wave angle of attack [deg]', 'Wave peak period')
% % speed as a func of wave frequency with wave size coloured 
% Plot3Dim([sog_data; test_sog_data], [ForecastWaveFreq_data; test_ForecastWaveFreq_data],...
%     [ForecastWaveSize_data; test_ForecastWaveSize_data], 1.2, 1.6, ...
%     true, 'Vg [m/s]', 'Wave Frequency [Hz]','Wave Size' )
% % speed as a func of wave size with wave frequency coloured 
% Plot3Dim([sog_data; test_sog_data], [ForecastWaveSize_data; test_ForecastWaveSize_data], ...
%     [ForecastWaveFreq_data; test_ForecastWaveFreq_data], 6.5, 8, ...
%     true, 'Vg [m/s]', 'Wave Size [m]', 'Wave peak period')
% % speed as a func of wave frequency
% Plot2Dim([sog_data; test_sog_data],  [ForecastWaveFreq_data; test_ForecastWaveFreq_data], ...
%     'Wave Frequency [Hz]', 'Vg [m/s]' )
% % speed as a func of wave size 
% Plot2Dim([sog_data; test_sog_data], [ForecastWaveSize_data; test_ForecastWaveSize_data], ...
%      'Wave Size [m]','Vg [m/s]')
% % speed as a func of rellative wind direction with relative wind speed coloured 
% Plot3Dim(sog_data, messuredRelWindDir_data, messuredRelWindSpeed_data, 3, 6, ...
%     false, 'Vg [m/s]', 'Measured Relative wind direction [deg]','Wind Speed' )
% % speed as a func of current speed in surge direction of body 
% Plot2Dim([sog_data; test_sog_data], [CurrentSpeed_data; test_CurrentSpeed_data], ...
%     'Current Speed in Surge Direction [m/s]','Vg [m/s]')
% % speed as a func of current speed in surge direction of body 
% Plot2Dim([sog_data; test_sog_data],[ForcastWindSpeed_data; test_ForcastWindSpeed_data], ...
%       'Wind Speed in Surge Direction [m/s]','Vg [m/s]')

%% Redo for testing on different dataset
%% New Data to be saved for plots
new_lon_data = [];
new_lat_data = [];
new_sog_data = [];
new_ForecastWaveSize_data = [];
new_messuredRelWindDir_data = [];
new_messuredRelWindSpeed_data = [];
new_relWaveDir_data = [];
new_ForecastWaveFreq_data = [];
new_ForcastWindSpeed_data = [];
new_CurrentSpeed_data = [];
new_ForcastWindDir_data = [];
new_CurrentDir_data =[];
new_Vr_data = [];


%% load data
disp('Loading new data')
path = '../Mausund200703_132548/';
addpath(path);
gpsFix = load('GpsFix.mat');
RelativeWind = load('RelativeWind.mat');
EulerAngles = load('EulerAngles.mat');
rmpath(path)
load '../Weather/weatherData_2020-7-3_2020-7-4';
load '../Weather/currentweatherData_2020-7-3_2020-7-4';  
disp('Done loading new data')

%% Format and interpolations
disp('Formating')
gps_data = gpsFix.GpsFix;
windData = RelativeWind.RelativeWind;
EulerAngles = EulerAngles.EulerAngles;
EulerAngles.psi = ssa(EulerAngles.psi,'deg');
messuredRelWindDir = interp1(windData.timestamp, ssa(windData.angle,'deg'),gps_data.timestamp);
messuredRelWindSpeed = interp1(windData.timestamp, windData.speed,gps_data.timestamp);
old_hour = 100000;
disp('Done formating')

%% run
disp('Start running through data')
for m = pad : n:length(gps_data.sog) - pad
    curr_hour = floor(double(gps_data.utc_time(m))/3600) ...
        + 24*(double(gps_data.utc_day(m)-gps_data.utc_day(1)));

    % Latidtude and longitude position of the vessel
    lat = mean(rad2deg(gps_data.lat(m : m + n)));
    lon = mean(rad2deg(gps_data.lon(m : m + n)));

    % Heading, Cog and Sog
    cog = rad2deg(mean(gps_data.cog(m : m + n)));
    psi = rad2deg(mean(EulerAngles.psi(m : m + n)));
    sog = mean(gps_data.sog(m : m + n));

    % Find position in wave data
    error_map = sqrt((latitudeMapWave - lat).^2 + (longitudeMapWave - lon).^2);
    [x,y] = find(error_map == min(error_map, [], 'all'));
    if x > xmax; xmax = x; end
    if y > ymax; ymax = y; end
    if x < xmin; xmin = x; end
    if y < ymin; ymin = y; end
    new_lon_data = cat(1,new_lon_data, lon);
    new_lat_data = cat(1,new_lat_data, lat);

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
    CurVsVelAnglre = ssa(VrDir- VcDir,'deg');

    % magnitude of the current
    currentSpeed = norm(Vc);
    VrSpeed = norm(Vr);

    % Messured wind speed and direction relative to the vessel
    curMessuredRelWindDir = mean(messuredRelWindDir(m : m + n));
    curMessuredRelWindSpeed = mean(messuredRelWindSpeed(m : m + n));
    relWaveDir = ssa(psi - curWaveDir - 180, 'deg');
    ForecastWaveFreq =  2*pi/waveHZ(x,y,curr_hour+1);
    ForecastEncWaveFreq = abs(ForecastWaveFreq -((ForecastWaveFreq^2*sog*cos(deg2rad(relWaveDir))/9.81)));
    ForecastWaveSize = waveSize(x, y, curr_hour + 1);
    currentSurge = currentSpeed*cos(ssa(deg2rad(VcDir - psi)));
    windSurge = ForcastWindSpeed*cos(ssa(deg2rad(curWindDir - psi)));


    % Save data
    new_ForecastWaveFreq_data = cat(1,new_ForecastWaveFreq_data, ForecastEncWaveFreq);            
    new_ForecastWaveSize_data = cat(1, new_ForecastWaveSize_data, ForecastWaveSize);
    new_messuredRelWindDir_data = cat(1, new_messuredRelWindDir_data, curMessuredRelWindDir);
    new_messuredRelWindSpeed_data = cat(1, new_messuredRelWindSpeed_data, curMessuredRelWindSpeed);
    new_Vr_data = cat(1, new_Vr_data,VrSpeed);
    new_sog_data = cat(1, new_sog_data,sog);
    new_ForcastWindSpeed_data = cat(1, new_ForcastWindSpeed_data, ForcastWindSpeed);
    new_ForcastWindDir_data = cat(1, new_ForcastWindDir_data, ssa(curWindDir - psi,'deg'));
    new_CurrentSpeed_data = cat(1, new_CurrentSpeed_data, currentSpeed);
    new_CurrentDir_data = cat(1, new_CurrentDir_data, ssa(VcDir - psi,'deg'));
    new_relWaveDir_data = cat(1, new_relWaveDir_data, relWaveDir);

    if old_hour ~= curr_hour
        str = sprintf('| Day: %d  | Hour: %d \t|', ...
            (floor(curr_hour/24)+1) + gps_data.utc_day(1)-1, (mod(curr_hour,24)));
        disp(str)
        old_hour = curr_hour;
    end
end
disp('Run Success')
%%
X_linear_new = [new_ForecastWaveSize_data new_ForecastWaveFreq_data abs(cos(deg2rad(new_relWaveDir_data))) ...
    new_ForcastWindSpeed_data.*cos(deg2rad(new_ForcastWindDir_data)) new_CurrentSpeed_data.*cos(deg2rad(new_CurrentDir_data))  ...
          ones(length(new_sog_data),1)];
X_Gauss_new = [new_CurrentSpeed_data new_CurrentDir_data new_relWaveDir_data new_ForcastWindDir_data ...
    new_ForcastWindSpeed_data new_ForecastWaveSize_data new_ForecastWaveFreq_data];
%%
PlotLinear(new_sog_data, w1,X_linear_new,'Vg')
%%
PlotGaus(new_sog_data, Mdl1,X_Gauss_new,'Vg')
%%
PlotGaus(new_sog_data, Mdl7,X_Gauss_new,'Vg')
PlotGaus(new_sog_data, Mdl8,X_Gauss_new,'Vg')

%% Bootstrap aggregating
[out2,~,yc2]= predict(Mdl2, X_Gauss_new);
[out3,~,yc3 ]= predict(Mdl3, X_Gauss_new);
[out4 ,~,yc4 ]= predict(Mdl4, X_Gauss_new);
[out5 ,~,yc5 ]= predict(Mdl5, X_Gauss_new);
[out6,~,yc6 ]= predict(Mdl6, X_Gauss_new);
bootOut = (out2 + out3 + out4 + out5 + out6)/5;

%% Stolen from Analysis file
diff1 = [];
sog_MSE = 0;

mean_output = mean(bootOut);
sog_summ1 = 0;
sog_summ2 = 0;
sog_summ3 = 0;
for i = 1:length(new_sog_data)
    out2 = predict(Mdl2, X_Gauss_new(i,:));
    out3 = predict(Mdl3, X_Gauss_new(i,:));
    out4 = predict(Mdl4, X_Gauss_new(i,:));
    out5 = predict(Mdl5, X_Gauss_new(i,:));
    out6 = predict(Mdl6, X_Gauss_new(i,:));
    out = (out2 + out3 + out4 + out5 + out6)/5;
    sog_summ1 = sog_summ1 + (out-mean_output)*(new_sog_data(i) - mean(new_sog_data));
    sog_summ2 = sog_summ2 + (out-mean_output)^2;
    sog_summ3 = sog_summ3 + (new_sog_data(i) - mean(new_sog_data))^2;
    diff1 = cat(1,diff1, out - new_sog_data(i));
    sog_MSE = sog_MSE + (new_sog_data(i) - out)^2;
end
sog_r = sog_summ1/sqrt(sog_summ2*sog_summ3);
sog_RMSE = sqrt(sog_MSE/length(new_sog_data));
figure;
scatter(linspace(1,1,length(diff1)),diff1)
hold on
boxplot(diff1)
tittel  = join(['Error between guessed ', 'Vg',  ' and actual ', 'Vg']);
title(tittel)
hold off
%%
figure;
scatter(new_sog_data, bootOut) 
hold on
[p] = polyfit(new_sog_data,bootOut,1);
x1 = linspace(min(new_sog_data),max(new_sog_data), length(new_sog_data));
y1 = polyval(p,x1);
plot(x1,y1)
plot(x1,x1, 'k--')
legend('Data',' Fit', 'Y = T', 'Location', 'NorthWest')
tittel = join(['RMSE = ', num2str(sog_RMSE,4),  ', R = ', num2str(sog_r,4)]);
xlabel(join(['Actual ', string, ' [m/s]']))
ylabel(join(['Predicted ', string, ' [m/s]']))
title(tittel)

%%
disp('Script done')



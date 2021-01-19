%% Clear workspace
clc; clearvars; close all;

%% Data That can be downladed from neptus that are relevant
% AbsoluteWind,Depth,DesiredHeading,DesiredPath,DesiredSpeed,DesiredZ,GpsFix,RelativeWind,RemoteSensorInfo,EstimatedState,EulerAngles
% GpsFix,RelativeWind,EulerAngles

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
for i = 2:11
    disp('Loading new data')
    %% load data
    if i == 1
        path = './Mausund200701_181204/';
        addpath(path);
        gpsFix = load('GpsFix.mat');
        RelativeWind = load('RelativeWind.mat');
        EulerAngles = load('EulerAngles.mat');
        rmpath(path)
        load('weatherData_2020-7-1_2020-7-2.mat') % Must be downloaded locally
        load('currentweatherData_2020-7-1_2020-7-3.mat') % Must be downloaded locally
        disp('Done loading data')
    end
    if i == 2
        path = './Mausund200701_221241/';
        addpath(path);
        gpsFix = load('GpsFix.mat');
        RelativeWind = load('RelativeWind.mat');
        EulerAngles = load('EulerAngles.mat');
        load('weatherData_2020-7-1_2020-7-2.mat')
        load('currentweatherData_2020-7-1_2020-7-3.mat')
        rmpath(path)
        disp('Done loading data')
    end
    if i == 3
        path = './Mausund200703_080820/';
        addpath(path);
        gpsFix = load('GpsFix.mat');
        RelativeWind = load('RelativeWind.mat');
        EulerAngles = load('EulerAngles.mat');
        load('weatherData_2020-7-3_2020-7-4.mat')
        load('currentweatherData_2020-7-3_2020-7-4.mat')
        rmpath(path)
        disp('Done loading data')
    end
    if i == 4
        path = './Mausund200703_132548/';
        addpath(path);
        gpsFix = load('GpsFix.mat');
        RelativeWind = load('RelativeWind.mat');
        EulerAngles = load('EulerAngles.mat');
        load('weatherData_2020-7-3_2020-7-4.mat')
        load('currentweatherData_2020-7-3_2020-7-4.mat')
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
        load('weatherData_2020-7-5_2020-7-5.mat')
        load('currentweatherData_2020-7-5_2020-7-5.mat')
        disp('Done loading data')
    end
    if i == 6
        path = './Mausund200706_154608/';
        addpath(path);
        gpsFix = load('GpsFix.mat');
        RelativeWind = load('RelativeWind.mat');
        EulerAngles = load('EulerAngles.mat');
        rmpath(path)
        load('weatherData_2020-7-6_2020-7-6.mat')
        load('currentweatherData_2020-7-6_2020-7-6.mat')
        disp('Done loading data')
    end
     if i == 7
        path = './Mausund200709_53748/';
        addpath(path);
        gpsFix = load('GpsFix.mat');
        RelativeWind = load('RelativeWind.mat');
        EulerAngles = load('EulerAngles.mat');
        rmpath(path)
        load('weatherData_2020-7-9_2020-7-9.mat')
        load('currentweatherData_2020-7-9_2020-7-9.mat')
        disp('Done loading data')
     end
    if i == 8
        path = './Trondheim082546/';
        addpath(path);
        gpsFix = load('GpsFix.mat');
        RelativeWind = load('RelativeWind.mat');
        EulerAngles = load('EulerAngles.mat');
        rmpath(path)
        load('weatherData_2020-2-20_2020-2-20.mat')
        load('currentweatherData_2020-2-20_2020-2-20.mat')
        disp('Done loading data')
    end
    if i == 9
        path = './Trondheim111728/';
        addpath(path);
        gpsFix = load('GpsFix.mat');
        RelativeWind = load('RelativeWind.mat');
        EulerAngles = load('EulerAngles.mat');
        rmpath(path)
        disp('Done loading data')
    end
    if i == 10
        path = './Trondheum094058/';
        addpath(path);
        gpsFix = load('GpsFix.mat');
        RelativeWind = load('RelativeWind.mat');
        EulerAngles = load('EulerAngles.mat');
        rmpath(path)
        disp('Done loading data')
    end
    if i == 11
        path = './Trondheum101916/';
        addpath(path);
        gpsFix = load('GpsFix.mat');
        RelativeWind = load('RelativeWind.mat');
        EulerAngles = load('EulerAngles.mat');
        rmpath(path)
        disp('Done loading data')
    end
    %% Format and interpolations
    gps_data = gpsFix.GpsFix;
    windData = RelativeWind.RelativeWind;
    EulerAngles = EulerAngles.EulerAngles;
    EulerAngles.psi = ssa(EulerAngles.psi,'deg');
    messuredRelWindDir = interp1(windData.timestamp, ssa(windData.angle,'deg' ),gps_data.timestamp);
    messuredRelWindSpeed = interp1(windData.timestamp, windData.speed,gps_data.timestamp);
    first = true;
    disp('Done formating')
    disp('Start running through data')
    %% run
    for m = (10*120) : length(gps_data.sog) - (10*120)
        if ~mod(gps_data.utc_time(m),avrager)
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
            CurVsVelAnglre = VrDir- VcDir;
            
            % magnitude of the current
            currentSpeed = norm(Vc);
            VrSpeed = norm(Vr);
           
            % Messured wind speed and direction relative to the vessel
            curMessuredRelWindDir = mean(messuredRelWindDir(m-avrager:m+avrager));
            curMessuredRelWindSpeed = mean(messuredRelWindSpeed(m-avrager:m+avrager));
            
            %currentDir = rad2deg(acos(dot(Sog_vector,wave_current_vector)/(norm(Sog_vector)*norm(wave_current_vector))));
            
            % Save current data
            if mod(count,5)
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
            else
                if waveSize(x, y, curr_hour + 1) < 0.001
                   test_ForecastWaveSize_data = cat(1,test_ForecastWaveSize_data,test_ForecastWaveSize_data(end));
                else
                   test_ForecastWaveSize_data = cat(1,test_ForecastWaveSize_data, waveSize(x, y, curr_hour + 1));
                end
                test_Vr_data = cat(1, test_Vr_data,VrSpeed);
                test_sog_data = cat(1, test_sog_data,sog);
                test_ForecastWaveFreq_data = cat(1,test_ForecastWaveFreq_data, 1/waveHZ(x,y,curr_hour+1));
                test_relWaveDir_data = cat(1,test_relWaveDir_data, ssa(ssa(curWaveDir - psi, 'deg') , 'deg'));
                test_ForcastWindSpeed_data = cat(1,test_ForcastWindSpeed_data, ForcastWindSpeed*cos(ssa(deg2rad(curWindDir - psi))));
                test_CurrentSpeed_data = cat(1,test_CurrentSpeed_data, currentSpeed*cos(ssa(deg2rad(VcDir - psi))));
            end
            if ~mod(gps_data.utc_time(m),3600) || first
                str = sprintf('| Day: %d  | Hour: %d \t|', ...
                    (floor(curr_hour/24)+1) + gps_data.utc_day(1)-1, (mod(curr_hour,24)));
                disp(str)
                first = false;
            end
            count = count + 1;
        end
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
meanSog = mean(sog_data);
%%
X = [ForecastWaveSize_data ForecastWaveFreq_data abs(cos(deg2rad(relWaveDir_data))) ...
    ForcastWindSpeed_data CurrentSpeed_data ones(length(sog_data),1)];
Y1 = sog_data;
Y1_test = test_sog_data;
w1 = (X'*X)\(X'*Y1);
X_test = [test_ForecastWaveSize_data test_ForecastWaveFreq_data abs(cos(deg2rad(test_relWaveDir_data))) ...
    test_ForcastWindSpeed_data test_CurrentSpeed_data ones(length(test_sog_data),1)]; 
CorrData1 = [Y1 X(:, 1:end-1)];
corrCoefs1 = corrcoef(CorrData1);
Y2 = Vr_data;
Y2_test = test_Vr_data;
w2 = (X'*X)\(X'*Y2);
CorrData2 = [Y2 X(:, 1:end-1)];
corrCoefs2 = corrcoef(CorrData2);
%%
diff1 = [];
sog_MSE = 0;
output = [];

for i = 1:length(test_sog_data)
    out = w1'*X_test(i, :)';
    output = cat(1, output, out);
end
mean_output = mean(output);
sog_summ1 = 0;
sog_summ2 = 0;
sog_summ3 = 0;
for i = 1:length(test_sog_data)
    out = w1'*X_test(i, :)';
    sog_summ1 = sog_summ1 + (out-mean_output)*(test_sog_data(i) - mean(test_sog_data));
    sog_summ2 = sog_summ2 + (out-mean_output)^2;
    sog_summ3 = sog_summ3 + (test_sog_data(i) - mean(test_sog_data))^2;
    diff1 = cat(1,diff1, out - test_sog_data(i));
    sog_MSE = sog_MSE + (test_sog_data(i) - out)^2;
end
sog_r = sog_summ1/sqrt(sog_summ2*sog_summ3);
sog_RMSE = sqrt(sog_MSE/length(test_sog_data));
figure;
scatter(linspace(1,1,length(diff1)),diff1)
hold on
boxplot(diff1)
title 'Error between guessed Vg and actual Vg'

hold off
%%
figure;
scatter(output, test_sog_data) 
hold on
[p] = polyfit(output,test_sog_data,1);
x1 = linspace(min(output),max(output), length(output));
y1 = polyval(p,x1);
plot(x1,y1)
plot(x1,x1, 'k--')
legend('Data',' Fit', 'Y = T', 'Location', 'SouthEast')
string = join(['RMSE = ', num2str(sog_RMSE,4),  ', R = ', num2str(sog_r,4)]);
ylabel 'Actual Vg [m/s]'
xlabel 'Predicted Vg [m/s]'
title(string)
%%
diff2 = [];
vr_MSE = 0;
output = [];
for i = 1:length(test_Vr_data)
    out = w2'*X_test(i, :)';
    output = cat(1, output, out);
end
mean_output = mean(output);
vr_summ1 = 0;
vr_summ2 = 0;
vr_summ3 = 0;
for i = 1:length(test_Vr_data)
    out = w2'*X_test(i, :)';
    diff2 = cat(1,diff2,out - test_Vr_data(i));
    vr_summ1 = vr_summ1 + (out-mean_output)*(test_Vr_data(i) - mean(test_Vr_data));
    vr_summ2 = vr_summ2 + (out-mean_output)^2;
    vr_summ3 = vr_summ3 + (test_Vr_data(i) - mean(test_Vr_data))^2;
    vr_MSE = vr_MSE + (out - test_Vr_data(i))^2;
end
vr_r = vr_summ1/sqrt(vr_summ2*vr_summ3);
vr_RMSE = sqrt(vr_MSE/length(test_Vr_data));
figure;
scatter(linspace(1,1,length(diff2)),diff2)
hold on
boxplot(diff2)
title 'Error between guessed Vr and actual Vr'
hold off
%%
figure;
scatter(output, test_Vr_data) 
hold on
p = polyfit(output,test_Vr_data,1);
x1 = linspace(min(output),max(output), length(output));
y1 = polyval(p,x1);
plot(x1,y1)
plot(x1,x1, 'k--')
legend('Data',' Fit', 'Y = T', 'Location', 'SouthEast')
string = join(['RMSE = ', num2str(vr_RMSE,4),  ', R = ', num2str(vr_r,4)]);
ylabel 'Actual Vr [m/s]'
xlabel 'Predicted Vr [m/s]'
title(string)
%%
figure
yvalues = {'Vg','ForecastWaveSize','ForecastWaveFreq','relWaveDir',...
    'WindSurgeSpeed', 'CurrentSurgeSpeed'};
xvalues = {'Vg','ForecastWaveSize','ForecastWaveFreq','relWaveDir',...
    'WindSurgeSpeed', 'CurrentSurgeSpeed'};
h = heatmap(xvalues,yvalues,corrCoefs1);
h.Title = 'Correlation Matrix';

%%
figure
yvalues = {'Vr','ForecastWaveSize','ForecastWaveFreq','relWaveDir',...
    'WindSurgeSpeed', 'CurrentSurgeSpeed'};
xvalues = {'Vr','ForecastWaveSize','ForecastWaveFreq','relWaveDir',...
    'WindSurgeSpeed', 'CurrentSurgeSpeed'};
h = heatmap(xvalues,yvalues,corrCoefs2);
h.Title = 'Correlation Matrix';
%%
disp('Plotting Data')
table1 = [];table2 = [];table3 = [];
for i = 1: length(ForecastWaveSize_data)
    if ForecastWaveSize_data(i) < 1.2
        table1 = cat(1,table1,[relWaveDir_data(i) Vr_data(i)]);
    elseif ForecastWaveSize_data(i) < 1.6
        table2 = cat(1,table2,[relWaveDir_data(i) Vr_data(i)]);
    else
        table3 = cat(1,table3,[relWaveDir_data(i) Vr_data(i)]);
    end
end
figure;
scatter(table1(:,1), table1(:,2))
hold on 
scatter(table2(:,1), table2(:,2))
scatter(table3(:,1), table3(:,2))
legend('Wave Size < 1.2', '1.2 < Wave Size < 1.6', 'Wave Size > 1.6')
xlabel 'Wave angle of attack [deg]',ylabel 'Vr [m/s]';
hold off

%% 

table1 = [];table2 = [];table3 = [];
for i = 1: length(ForecastWaveFreq_data)
    if ForecastWaveFreq_data(i) < 0.13
        table1 = cat(1,table1,[relWaveDir_data(i) Vr_data(i)]);
    elseif ForecastWaveFreq_data(i) < 0.16
        table2 = cat(1,table2,[relWaveDir_data(i) Vr_data(i)]);
    else
        table3 = cat(1,table3,[relWaveDir_data(i) Vr_data(i)]);
    end
end
figure;
scatter(table1(:,1), table1(:,2))
hold on 
scatter(table2(:,1), table2(:,2))
scatter(table3(:,1), table3(:,2))
legend('Wave Hz < 0.13', '0.13 < Wave Hz < 0.16', 'Wave Hz > 1.16')
xlabel 'Wave angle of attack [deg]',ylabel 'Vr [m/s]';
hold off
%%
table1 = [];table2 = [];table3 = [];
for i = 1: length(ForecastWaveSize_data)
    if ForecastWaveSize_data(i) < 1.2
        table1 = cat(1,table1,[ForecastWaveFreq_data(i) Vr_data(i)]);
    elseif ForecastWaveSize_data(i) < 1.6
        table2 = cat(1,table2,[ForecastWaveFreq_data(i) Vr_data(i)]);
    else
        table3 = cat(1,table3,[ForecastWaveFreq_data(i) Vr_data(i)]);
    end
end
figure;
scatter(table1(:,1), table1(:,2))
hold on 
scatter(table2(:,1), table2(:,2))
scatter(table3(:,1), table3(:,2))
legend('Wave Size < 1.2', '1.2 < Wave Size < 1.6', 'Wave Size > 1.6')
xlabel 'Wave Frequency [Hz]',ylabel 'Vr [m/s]';
p = polyfit(ForecastWaveFreq_data,Vr_data,1);
x1 = linspace(min(ForecastWaveFreq_data),max(ForecastWaveFreq_data), length(ForecastWaveFreq_data));
y1 = polyval(p,x1);
plot(x1,y1)
legend('Wave Size < 1.2', '1.2 < Wave Size < 1.6', 'Wave Size > 1.6', 'Linear Regression')
hold off
%%
table1 = [];table2 = [];table3 = [];
for i = 1: length(ForecastWaveSize_data)
    if ForecastWaveFreq_data(i) < 0.13
        table1 = cat(1,table1,[ForecastWaveSize_data(i) Vr_data(i)]);
    elseif ForecastWaveFreq_data(i) < 0.16
        table2 = cat(1,table2,[ForecastWaveSize_data(i) Vr_data(i)]);
    else
        table3 = cat(1,table3,[ForecastWaveSize_data(i) Vr_data(i)]);
    end
end
figure;
scatter(table1(:,1), table1(:,2))
hold on 
scatter(table2(:,1), table2(:,2))
scatter(table3(:,1), table3(:,2))

xlabel 'Wave Size [m]',ylabel 'Vr [m/s]';
p = polyfit(ForecastWaveSize_data,Vr_data,1);
x1 = linspace(min(ForecastWaveSize_data),max(ForecastWaveSize_data), length(ForecastWaveSize_data));
y1 = polyval(p,x1);
plot(x1,y1)
legend('Wave Hz < 0.13', '0.13 < Wave Hz < 0.16', 'Wave Hz > 0.16', 'Linear Regression')
hold off
%%
table1 = [];table2 = [];table3 = [];
for i = 1: length(messuredRelWindSpeed_data)
    if messuredRelWindSpeed_data(i) < 3
        table1 = cat(1,table1,[messuredRelWindDir_data(i) Vr_data(i)]);
    elseif messuredRelWindSpeed_data(i) < 6
        table2 = cat(1,table2,[messuredRelWindDir_data(i) Vr_data(i)]);
    else
        table3 = cat(1,table3,[messuredRelWindDir_data(i) Vr_data(i)]);
    end
end
figure;
scatter(table1(:,1), table1(:,2))
hold on 
scatter(table2(:,1), table2(:,2))
scatter(table3(:,1), table3(:,2))
legend('Wind Speed < 3', '3 < Wind Speed < 6', 'Wind Speed > 6')
xlabel 'Measured Relative wind direction',ylabel 'Vr';
hold off
%%
figure;
scatter(CurrentSpeed_data, Vr_data)
hold on 
p = polyfit(CurrentSpeed_data, Vr_data, 1);
x1 = linspace(min(CurrentSpeed_data), max(CurrentSpeed_data), length(CurrentSpeed_data));
y1 = polyval(p,x1);
plot(x1,y1)
%legend('Curent speed < 0.15', '0.15< Curent speed <0.3', 'Curent speed > 0.3')
xlabel 'Current Speed in Surge Direction [m/s]',ylabel 'Vr [m/s]';
hold off
pause(0.01)
%%
figure;
scatter(ForcastWindSpeed_data, Vr_data)
hold on 
p = polyfit(ForcastWindSpeed_data, Vr_data, 1);
x1 = linspace(min(ForcastWindSpeed_data), max(ForcastWindSpeed_data), length(ForcastWindSpeed_data));
y1 = polyval(p,x1);
plot(x1,y1)
xlabel 'Wind Speed in Surge Direction [m/s]',ylabel 'Vr [m/s]';
hold off
pause(0.01)
%% Redo for testing

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
test_sog_data = [];
test_ForecastWaveSize_data = [];
test_ForecastWaveFreq_data = [];
test_relWaveDir_data = [];
test_ForcastWindSpeed_data = [];
test_CurrentSpeed_data = [];
test_Vr_data = [];

xmax = 0; ymax = 0; ymin = inf; xmin = inf;
avrager = 6*60; % average over x min
count = 2;
disp('Loading new data')
%% load data
path = './Mausund200701_181204/';
addpath(path);
gpsFix = load('GpsFix.mat');
RelativeWind = load('RelativeWind.mat');
EulerAngles = load('EulerAngles.mat');
rmpath(path)
load('weatherData_2020-7-1_2020-7-2.mat')
load('currentweatherData_2020-7-1_2020-7-3.mat')
disp('Done loading data')
%% Format and interpolations
gps_data = gpsFix.GpsFix;
windData = RelativeWind.RelativeWind;
EulerAngles = EulerAngles.EulerAngles;
EulerAngles.psi = ssa(EulerAngles.psi,'deg');
messuredRelWindDir = interp1(windData.timestamp, ssa(windData.angle,'deg' ),gps_data.timestamp);
messuredRelWindSpeed = interp1(windData.timestamp, windData.speed,gps_data.timestamp);
first = true;
disp('Done formating')
disp('Start running through data')
%% run
for m = (10*120) : length(gps_data.sog) - (10*120)
    if ~mod(gps_data.utc_time(m),avrager)
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
        CurVsVelAnglre = VrDir- VcDir;

        % magnitude of the current
        currentSpeed = norm(Vc);
        VrSpeed = norm(Vr);

        % Messured wind speed and direction relative to the vessel
        curMessuredRelWindDir = mean(messuredRelWindDir(m-avrager:m+avrager));
        curMessuredRelWindSpeed = mean(messuredRelWindSpeed(m-avrager:m+avrager));

        %currentDir = rad2deg(acos(dot(Sog_vector,wave_current_vector)/(norm(Sog_vector)*norm(wave_current_vector))));

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
        if ~mod(gps_data.utc_time(m),3600) || first
            str = sprintf('| Day: %d  | Hour: %d \t|', ...
                (floor(curr_hour/24)+1) + gps_data.utc_day(1)-1, (mod(curr_hour,24)));
            disp(str)
            first = false;
        end
        count = count + 1;
    end
end
disp('Run Success')
%%
X_testt = [ForecastWaveSize_data ForecastWaveFreq_data abs(cos(deg2rad(relWaveDir_data))) ...
    ForcastWindSpeed_data CurrentSpeed_data, ones(length(sog_data),1)];
%%
diff1 = [];
sog_MSE = 0;
output = [];

for i = 1:length(sog_data)
    out = w1'*X_testt(i, :)';
    output = cat(1, output, out);
end
mean_output = mean(output);
sog_summ1 = 0;
sog_summ2 = 0;
sog_summ3 = 0;
for i = 1:length(sog_data)
    out = w1'*X_testt(i, :)';
    sog_summ1 = sog_summ1 + (out-mean_output)*(sog_data(i) - mean(sog_data));
    sog_summ2 = sog_summ2 + (out-mean_output)^2;
    sog_summ3 = sog_summ3 + (sog_data(i) - mean(sog_data))^2;
    diff1 = cat(1,diff1, out - sog_data(i));
    sog_MSE = sog_MSE + (sog_data(i) - out)^2;
end
sog_r = sog_summ1/sqrt(sog_summ2*sog_summ3);
sog_RMSE = sqrt(sog_MSE/length(sog_data));
figure;
scatter(linspace(1,1,length(diff1)),diff1)
hold on
boxplot(diff1)
title 'Error between guessed Vg and actual Vg'
hold off
%%
figure;
scatter(output, sog_data) 
hold on
[p] = polyfit(output,sog_data,1);
x1 = linspace(min(output),max(output), length(output));
y1 = polyval(p,x1);
plot(x1,y1)
plot(x1,x1, 'k--')
legend('Data',' Fit', 'Y = T', 'Location', 'SouthEast')
string = join(['RMSE = ', num2str(sog_RMSE,4),  ', R = ', num2str(sog_r,4)]);
ylabel 'Actual Vg [m/s]'
xlabel 'Predicted Vg [m/s]'
title(string)
%%
diff2 = [];
vr_MSE = 0;
output = [];
for i = 1:length(Vr_data)
    out = w2'*X_testt(i, :)';
    output = cat(1, output, out);
end
mean_output = mean(output);
vr_summ1 = 0;
vr_summ2 = 0;
vr_summ3 = 0;
for i = 1:length(Vr_data)
    out = w2'*X_testt(i, :)';
    diff2 = cat(1,diff2,out - Vr_data(i));
    vr_summ1 = vr_summ1 + (out-mean_output)*(Vr_data(i) - mean(Vr_data));
    vr_summ2 = vr_summ2 + (out-mean_output)^2;
    vr_summ3 = vr_summ3 + (Vr_data(i) - mean(Vr_data))^2;
    vr_MSE = vr_MSE + (out - Vr_data(i))^2;
end
vr_r = vr_summ1/sqrt(vr_summ2*vr_summ3);
vr_RMSE = sqrt(vr_MSE/length(Vr_data));
figure;
scatter(linspace(1,1,length(diff2)),diff2)
hold on
boxplot(diff2)
title 'Error between guessed Vr and actual Vr'
hold off
%%
figure;
scatter(output, Vr_data) 
hold on
[p,S] = polyfit(output,Vr_data,1);
x1 = linspace(min(output),max(output), length(output));
y1 = polyval(p,x1);
plot(x1,y1)
plot(x1,x1, 'k--')
legend('Data',' Fit', 'Y = T', 'Location', 'SouthEast')
string = join(['RMSE = ', num2str(vr_RMSE,4),  ', R = ', num2str(vr_r,4)]);
ylabel 'Actual Vr [m/s]'
xlabel 'Predicted Vr [m/s]'
title(string)

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
%         net = temp_net; performance = temp_perform;
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



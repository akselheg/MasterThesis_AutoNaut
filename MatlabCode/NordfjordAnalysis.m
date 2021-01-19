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
avrager = 3*60; % average over x min
count = 1;
disp('Loading new data')
path = './Nordfjord082813/';
addpath(path);
gpsFix = load('GpsFix.mat');
RelativeWind = load('RelativeWind.mat');
EulerAngles = load('EulerAngles.mat');
rmpath(path)
load('weatherData_2020-5-28_2020-5-28.mat')
load('currentweatherData_2020-5-28_2020-5-29.mat')
disp('Done loading data')


%% Format and interpolations
gps_data = gpsFix.GpsFix;
windData = RelativeWind.RelativeWind;
EulerAngles = EulerAngles.EulerAngles;
EulerAngles.psi = ssa(EulerAngles.psi,'deg');
messuredRelWindDir = interp1(windData.timestamp, ssa(windData.angle,'deg' ),gps_data.timestamp);
messuredRelWindSpeed = interp1(windData.timestamp, windData.speed,gps_data.timestamp);
[latSize,longSize] = size(latitudeMapWave);
[curlatSize,curlongSize] = size(latitudeCurrentMap);
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
                disp('eeeeee')
                ForecastWaveSize_data = cat(1, ForecastWaveSize_data, ForecastWaveSize_data(end));
            else
                ForecastWaveSize_data = cat(1, ForecastWaveSize_data, waveSize(x, y, curr_hour + 1));
            end
            Vr_data = cat(1, Vr_data,VrSpeed);
            sog_data = cat(1, sog_data,sog);
            relWaveDir_data = cat(1, relWaveDir_data, ssa(psi- curWaveDir - 180, 'deg'));


            messuredRelWindDir_data = cat(1, messuredRelWindDir_data, curMessuredRelWindDir);
            messuredRelWindSpeed_data = cat(1, messuredRelWindSpeed_data, curMessuredRelWindSpeed);

            ForcastWindSpeed_data = cat(1, ForcastWindSpeed_data, ForcastWindSpeed*cos(ssa(deg2rad(psi-curWindDir))));
            CurrentDir_data = cat(1, CurrentDir_data, ssa(psi - VcDir,'deg'));
            CurrentSpeed_data = cat(1, CurrentSpeed_data, currentSpeed*cos(ssa(deg2rad(psi-VcDir))));
        else
            if waveSize(x, y, curr_hour + 1) < 0.001
               test_ForecastWaveSize_data = cat(1,test_ForecastWaveSize_data,test_ForecastWaveSize_data(end));
            else
               test_ForecastWaveSize_data = cat(1,test_ForecastWaveSize_data, waveSize(x, y, curr_hour + 1));
            end
            test_Vr_data = cat(1, test_Vr_data,VrSpeed);
            test_sog_data = cat(1, test_sog_data,sog);
            test_ForecastWaveFreq_data = cat(1,test_ForecastWaveFreq_data, 1/waveHZ(x,y,curr_hour+1));
            test_relWaveDir_data = cat(1,test_relWaveDir_data, ssa(ssa(psi- curWaveDir, 'deg') , 'deg'));
            test_ForcastWindSpeed_data = cat(1,test_ForcastWindSpeed_data, ForcastWindSpeed*cos(ssa(deg2rad(psi-curWindDir))));
            test_CurrentSpeed_data = cat(1,test_CurrentSpeed_data, currentSpeed*cos(ssa(deg2rad(psi-VcDir))));
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
%%
i = 1; 
while i < length(ForecastWaveSize_data)
    if isnan(ForecastWaveSize_data(i)) ||isnan(ForecastWaveFreq_data(i)) ||...
            isnan(relWaveDir_data(i)) || isnan(ForcastWindSpeed_data(i)) ||...
            isnan(CurrentSpeed_data(i))
        ForecastWaveSize_data = ForecastWaveSize_data([1:i-1 i+1:end]);
        ForecastWaveFreq_data = ForecastWaveFreq_data([1:i-1 i+1:end]);
        relWaveDir_data = relWaveDir_data([1:i-1 i+1:end]);
        ForcastWindSpeed_data = ForcastWindSpeed_data([1:i-1 i+1:end]);
        CurrentSpeed_data = CurrentSpeed_data([1:i-1 i+1:end]);
        sog_data = sog_data([1:i-1 i+1:end]);
        Vr_data = Vr_data([1:i-1 i+1:end]);
    else
        i = i + 1;
    end
end
i = 1; 

while i < length(test_ForecastWaveSize_data)
    if isnan(test_ForecastWaveSize_data(i)) ||isnan(test_ForecastWaveFreq_data(i)) ||...
            isnan(test_relWaveDir_data(i)) || isnan(test_ForcastWindSpeed_data(i)) ||...
            isnan(test_CurrentSpeed_data(i))
        test_ForecastWaveSize_data = test_ForecastWaveSize_data([1:i-1 i+1:end]);
        test_ForecastWaveFreq_data = test_ForecastWaveFreq_data([1:i-1 i+1:end]);
        test_relWaveDir_data = test_relWaveDir_data([1:i-1 i+1:end]);
        test_ForcastWindSpeed_data = test_ForcastWindSpeed_data([1:i-1 i+1:end]);
        test_CurrentSpeed_data = test_CurrentSpeed_data([1:i-1 i+1:end]);
        test_sog_data = test_sog_data([1:i-1 i+1:end]);
        test_Vr_data = test_Vr_data([1:i-1 i+1:end]);
    else
        i = i + 1;
    end
end
X = [ForecastWaveSize_data ForecastWaveFreq_data abs(relWaveDir_data) ...
    ForcastWindSpeed_data CurrentSpeed_data];
Y1 = sog_data;
Y1_test = test_sog_data;
w1 = (X'*X)\(X'*Y1);
X_test = [test_ForecastWaveSize_data test_ForecastWaveFreq_data abs(test_relWaveDir_data) ...
    test_ForcastWindSpeed_data test_CurrentSpeed_data];
CorrData1 = [Y1 X];
corrCoefs1 = corrcoef(CorrData1);
Y2 = Vr_data;
Y2_test = test_Vr_data;
w2 = (X'*X)\(X'*Y2);
CorrData2 = [Y2 X];
corrCoefs2 = corrcoef(CorrData2);
%%
diff1 = [];
for i = 1:length(test_sog_data)
    out = w1'*X_test(i, :)';
    diff1 = cat(1,diff1, out - test_sog_data(i));
end
figure;
scatter(linspace(1,1,length(diff1)),diff1)
hold on
boxplot(diff1)
title 'Error between guessed SOG and actual SOG'
hold off
%%
diff2 = [];
for i = 1:length(test_Vr_data)
    out = w2'*X_test(i, :)';
    diff2 = cat(1,diff2,out - test_Vr_data(i));
end
figure;
scatter(linspace(1,1,length(diff2)),diff2)
hold on
boxplot(diff2)
title 'Error between guessed |Vr| and actual |Vr|'
hold off
%%
figure
yvalues = {'Sog','ForecastWaveSize','ForecastWaveFreq','relWaveDir',...
    'WindSurgeSpeed', 'CurrentSurgeSpeed'};
xvalues = {'Sog','ForecastWaveSize','ForecastWaveFreq','relWaveDir',...
    'WindSurgeSpeed', 'CurrentSurgeSpeed'};
h = heatmap(xvalues,yvalues,corrCoefs1);
h.Title = 'Covariance Matrix';

%%
figure
yvalues = {'|Vr|','ForecastWaveSize','ForecastWaveFreq','relWaveDir',...
    'WindSurgeSpeed', 'CurrentSurgeSpeed'};
xvalues = {'|Vr|','ForecastWaveSize','ForecastWaveFreq','relWaveDir',...
    'WindSurgeSpeed', 'CurrentSurgeSpeed'};
h = heatmap(xvalues,yvalues,corrCoefs2);
h.Title = 'Covariance Matrix';
%%

%% 
table1 = [];table2 = [];table3 = [];
for i = 1: length(ForecastWaveFreq_data)
    if ForecastWaveFreq_data(i) < 0.7
        table1 = cat(1,table1,[relWaveDir_data(i) Vr_data(i)]);
    elseif ForecastWaveFreq_data(i) < 0.75
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
legend('Hz < 0.7', '0.7 < Hz < 0.75', 'Hz > 0.75')
xlabel 'Relative wave direction',ylabel 'SOG';
hold off
%%
figure;
scatter(ForecastWaveSize_data, Vr_data)
hold on
xlabel 'Wave Size [m]',ylabel '|V_r| [m/s]';
p = polyfit(ForecastWaveSize_data,Vr_data,1);
x1 = linspace(min(ForecastWaveSize_data),max(ForecastWaveSize_data), length(ForecastWaveSize_data));
y1 = polyval(p,x1);
plot(x1,y1)
legend('Wave Hz < 0.7', '0.7 < Wave Hz < 0.75', 'Wave Hz > 0.75', 'Linear Regression')
hold off
%%
% table1 = [];table2 = [];table3 = [];
% for i = 1: length(messuredRelWindSpeed_data)
%     if messuredRelWindSpeed_data(i) < 20
%         table1 = cat(1,table1,[messuredRelWindDir_data(i) Vr_data(i)]);
%     elseif messuredRelWindSpeed_data(i) < 6
%         table2 = cat(1,table2,[messuredRelWindDir_data(i) Vr_data(i)]);
%     else
%         table3 = cat(1,table3,[messuredRelWindDir_data(i) Vr_data(i)]);
%     end
% end
% figure;
% scatter(table1(:,1), table1(:,2))
% hold on 
% %scatter(table2(:,1), table2(:,2))
% %scatter(table3(:,1), table3(:,2))
% %legend('Wind Speed < 3', '3 < Wind Speed < 6', 'Wind Speed > 6')
% xlabel 'Measured Relative wind direction',ylabel 'SOG';
% hold off
%%
figure;
scatter(CurrentSpeed_data, Vr_data)
hold on 
p = polyfit(CurrentSpeed_data, Vr_data, 1);
x1 = linspace(min(CurrentSpeed_data), max(CurrentSpeed_data), length(CurrentSpeed_data));
y1 = polyval(p,x1);
plot(x1,y1)
%legend('Curent speed < 0.15', '0.15< Curent speed <0.3', 'Curent speed > 0.3')
xlabel 'Current Speed in Surge Direction',ylabel '|V_r|';
hold off
pause(0.01)
%%
figure;
n = find(min(ForcastWindSpeed_data) ==  ForcastWindSpeed_data);
tempdata = ForcastWindSpeed_data([1:n-1 n+1:end]);
tempdata2 = Vr_data([1:n-1 n+1:end]);
scatter(ForcastWindSpeed_data, Vr_data)
hold on 
p = polyfit(tempdata, tempdata2, 1);
x1 = linspace(min(ForcastWindSpeed_data), max(ForcastWindSpeed_data), length(ForcastWindSpeed_data));
y1 = polyval(p,x1);
plot(x1,y1)
xlabel 'Wind Speed in Surge Direction [m/s]',ylabel 'SOG [m/s]';
hold off
%%
disp('Script done')

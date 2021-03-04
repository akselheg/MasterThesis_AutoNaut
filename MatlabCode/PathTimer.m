
clc;clearvars;close all;
r = 6371000;
% startLat = 64.17; startLon = 9;
% endLat = 64.4; endLon = 10;
startLat = 64.2298855766667; startLon = 8.89639759166667;
endLat = 64.3847535783333; endLon = 8.88621290666667;
distance = Haversine_deg(startLat, startLon, endLat, endLon,r);
angle = bearing(startLat,startLon,endLat,endLon);
n = floor(distance/100);
waypoints = zeros(floor(distance/20),2);
disp(num2str(rad2deg(angle)))
[lat_points, long_points] = pathCreator([startLat, endLat]', [startLon, endLon]', n);

load('SpeedEstimator.mat');
load('./Weather/weatherData_2020-7-3_2020-7-4.mat')
load('./Weather/currentweatherData_2020-7-3_2020-7-4.mat')
load('./Mausund200703_132548/GpsFix.mat')
startTime = GpsFix.utc_time(1);
time = startTime;
hours = floor(double(GpsFix.utc_time(1))/3600) ...
+ 24*(double(GpsFix.utc_day(1)-GpsFix.utc_day(1)));
%%
avrspeed = 0;
curhour = 0;
for i=1:n - 1
    error_map = sqrt((latitudeMapWave - lat_points(i)).^2 + (longitudeMapWave - long_points(i)).^2);
    [x,y] = find(error_map == min(error_map, [], 'all'));
    error_map = sqrt((latitudeCurrentMap - lat_points(i)).^2 + (longitudeCurrentMap - long_points(i)).^2);
    [xcurrent,ycurrent] = find(error_map == min(error_map, [], 'all'));
    curWaveSize = waveSize(x,y,hours+1);
    curWaveDir = waveDir(x,y,hours+1);
    curWaveHz = waveHZ(x,y,hours+1);
    curWindDir = windDir(x,y,hours+1);
    curWindSpeed = windSpeed(x,y,hours+1);
    curCurrentEast = currentEast(xcurrent,ycurrent,hours+1);
    curCurrentNorth = currentNorth(xcurrent,ycurrent,hours+1);
    if isnan(curWaveSize) || isnan(curCurrentEast)
        disp('oh nononono')
    end
    
    Vc = [curCurrentNorth; curCurrentEast];
    VcDir = atan2d( Vc(2), Vc(1));
    
    angle = bearing(lat_points(i),long_points(i),lat_points(i+1),long_points(i+1));
    disp(num2str(rad2deg(angle)))
    
    windSurge = curWindSpeed*cos(ssa(deg2rad(curWindDir) - angle));
    currentSurge = norm(Vc)*cos(ssa(deg2rad(VcDir) - angle));
    relWaveDir = ssa(rad2deg(angle)- curWaveDir - 180, 'deg');
    
    X = [curWaveSize curWaveHz (cos(deg2rad(relWaveDir)))  ...
        windSurge currentSurge];
    dist = Haversine_deg(lat_points(i),long_points(i),lat_points(i+1),long_points(i+1),r);
    %disp(num2str((dist)));
    speed = predict(Mdl1, X);
    avrspeed = avrspeed + speed;
    time = time + dist/speed;
    hours = floor(time/3600);
end
%%
elapseTime = time-startTime;
hours = floor(elapseTime/3600);
min = floor(mod(elapseTime,3600)/60);
sec = round(mod(elapseTime,60));
disp(['Total time used: ', num2str(hours), ' hours, ' num2str(min),' min, ', num2str(sec), ' sec.' ])
elapseTime = GpsFix.utc_time(end)-startTime;
hours = floor(elapseTime/3600);
min = floor(mod(elapseTime,3600)/60);
sec = round(mod(elapseTime,60));
disp(['Total time used: ', num2str(hours), ' hours, ' num2str(min),' min, ', num2str(sec), ' sec.' ])
avrspeed = avrspeed/(n-1);
figure;
geoplot(lat_points,long_points)
hold on
geoplot(rad2deg(GpsFix.lat),rad2deg(GpsFix.lon))
geoscatter(rad2deg(GpsFix.lat(1)),rad2deg(GpsFix.lon(1)),'g')
geoscatter(rad2deg(GpsFix.lat(end)),rad2deg(GpsFix.lon(end)),'r')
a = latitudeMapWave(1:100,100:250);
b = longitudeMapWave(1:100,100:250);
%geoscatter(a(:),b(:), 'w')
%quiver(longitudeMapWave,latitudeMapWave, windDir(:,:,2),windSpeed(:,:,2))
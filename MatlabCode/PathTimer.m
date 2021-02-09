
clc;clearvars;close all;
r = 6371000;

startx = 64.2298855766667; starty = 8.89639759166667;
endx = 64.3847535783333; endy = 8.88621290666667;
distance = Haversine_deg(startx, starty, endx, endy,r);
angle = bearing(startx,starty,endx,endy);
n = floor(distance/100);
waypoints = zeros(floor(distance/100),2);

[lat_points, long_points] = pathCreator([startx, endx]', [starty, endy]', n);
time = 0;
load('SpeedEstimator.mat');
load('./Weather/weatherData_2020-7-1_2020-7-2.mat')
load('./Weather/currentweatherData_2020-7-1_2020-7-3.mat')
avrspeed = 0;
curhour = 0;
hours = 0;
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
    %disp(num2str(rad2deg(angle)))
    
    windSurge = curWindSpeed*cos(ssa(deg2rad(curWindDir) - angle));
    currentSurge = norm(Vc)*cos(ssa(deg2rad(VcDir) - angle));
    relWaveDir = ssa(rad2deg(angle)- curWaveDir - 180, 'deg');
    
    X = [curWaveSize curWaveHz (cos(deg2rad(relWaveDir)))  ...
        windSurge currentSurge];
    dist = Haversine_deg(lat_points(i),long_points(i),lat_points(i+1),long_points(i+1),r);
    disp(num2str((dist)))
    speed = predict(Mdl1, X);
    avrspeed = avrspeed + speed;
    time = time + dist/speed;
    hours = floor(time/3600);
end
%%
min = floor(mod(time,3600)/60);
sec = round(mod(time,60));
disp(['Total time used: ', num2str(hours), ' hours, ' num2str(min),' min, ', num2str(sec), ' sec.' ])
avrspeed = avrspeed/(n-1);
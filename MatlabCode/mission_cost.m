function [time,v,sun] = mission_cost(input)
%Function calculating the mission cost of a mission
% 
%     INPUT: input - An array describing Longitude and Latitude for each
%            waypoint and the control input for each waypint
% 
%     OUTPUT: time - Variable describing the duration of the mission
%                    TODO: Soon to describe mission costqfgwervwrqwrfr\!!!!!!!
%             v    - Array including velocity of the vessel between each
%                    waypoint. 


% persistent current_dir current_speed wind_dir wind_speed wave_size wave_dir latitude_map longitude_map wave_Hz
% persistent lat_start long_start lat_end long_end
% persistent goals_lat goals_long goals_rad goals_val goals_num

load './Weather/weatherData_2020-7-1_2020-7-2.mat';
load './Weather/currentweatherData_2020-7-1_2020-7-3.mat';

load('constraints.mat')
lat_start = start_lat;
long_start = start_long;
lat_end = end_lat;
long_end = end_long;
load('measurement_goals.mat')
goals_lat = [];
goals_long = [];
goals_rad = rad_goals;
goals_val = [];
goals_num = [];
[inputs,~] = size(input);
latitude = [lat_start; input(1:inputs/3);lat_end];
longitude = [long_start; input(inputs/3 + 1: inputs*(2/3)); long_end];
sensor = [0;input(inputs*(2/3) + 1:inputs);0];

% implement solar radiation map
% implement energy states
% implement data transmission map
% tune all cost parameters
% ???
% sucsess
[steps,~] = size(longitude);
deg2rad = pi/180;
theGoals = goals_val;
r = 6371000;
t = zeros(1,steps);
v = zeros(steps-1,1);
sun =zeros(1,steps-1);

[lat_size_map,long_size_map] = size(latitude_map);
val = 0;



for i = 1:steps - 1
    seconds = t(i);
    hour = round(seconds/3600)*0 + 10;
    hour_real = seconds/3600 + 10;
    lat = latitude(i);
    long = longitude(i);
    
    % Find position on weather map
    long_lat = repmat(lat,lat_size_map,long_size_map);
    long_long = repmat(long,lat_size_map,long_size_map);
    pos_diff = Haversine_deg(long_lat,long_long,latitude_map,longitude_map,r);
    [min_dist, index] = minmat(pos_diff);
    
    % Find Mission goal Values. 
    goal_dist = Haversine_deg(goals_long, goals_lat, repmat(long,goals_num,1),repmat(lat,goals_num,1),r);
    inside  = (goal_dist<goals_rad)*sensor(i);
    val = val + sum(inside.*theGoals);
    theGoals = theGoals.*(~inside);
    
    % - Get current
%     curr_dir = current_dir(index(1),index(2),hour)*deg2rad;
%     curr_spd = current_speed(index(1),index(2),hour);
    % - Get wind
    wnd_dir = wind_dir(index(1),index(2),hour)*deg2rad;
    wnd_spd = wind_speed(index(1),index(2),hour);
    % - Get waves
    ww_dir = wave_dir(index(1),index(2),hour)*deg2rad;
    ww_spd = wave_size(index(1),index(2),hour);
% % %     waves = [ww_spd*cos(ww_dir); ww_spd*sin(ww_dir)];
    % - Describe full function (possibly newton method??
    course = bearing(latitude(i)*deg2rad,longitude(i)*deg2rad,latitude(i+1)*deg2rad,longitude(i+1)*deg2rad);
    testSum = ww_spd + wnd_spd;
    if isnan(testSum)
        t(steps) = inf;
        break
    end
    X = [ww_spd curWaveHz (cos(ssa(course- deg2rad(curWaveDir - 180))))  ...
        windSurge currentSurge];
    v(i) = predict(Mdl1, X);
    distance = Haversine_deg(latitude(i+1),longitude(i+1),latitude(i),longitude(i),r);
    if (v(i) < 0)||(~isreal(v(i)))
        t(steps) = inf;
        break
    end
    t(i+1) = t(i) + distance/v(i);
    
    % Calculating Battery energy
    sun(i) = cos(solar_angle(long,lat,hour_real/24));
end

time = t(steps)/3600 - val;
end


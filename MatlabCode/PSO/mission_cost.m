function [time,v] = mission_cost(input)
%Function calculating the mission cost of a mission
% 
%     INPUT: input - An array describing Longitude and Latitude for each
%            waypoint and the control input for each waypint
% 
%     OUTPUT: time - Variable describing the duration of the mission
%                    TODO: Soon to describe mission costqfgwervwrqwrfr\!!!!!!!
%             v    - Array including velocity of the vessel between each
%                    waypoint. 


persistent currentEast currentNorth windDir windSpeed waveSize waveDir latitudeMapWave ...
    latitudeCurrentMap longitudeCurrentMap longitudeMapWave waveHZ
persistent lat_start long_start lat_end long_end
persistent goals_lat goals_long goals_rad goals_val goals_num
persistent Mdl1

load 'speedEst.mat'

if isempty(latitudeMapWave)
    %path = '../Mausund200703_132548';
    
    load '../Weather/weatherData_2020-7-3_2020-7-4';
    load '../Weather/currentweatherData_2020-7-3_2020-7-4';
end

% error_map = sqrt((latitudeMapWave - lat).^2 + (longitudeMapWave - lon).^2);
% [x,y] = find(error_map == min(error_map, [], 'all'));
% 
% error_map = sqrt((latitudeCurrentMap - lat).^2 + (longitudeCurrentMap - lon).^2);
% [xcurrent,ycurrent] = find(error_map == min(error_map, [], 'all'));
%if isempty()

if isempty(lat_start)
    load('constraints.mat')
    lat_start = start_lat;
    long_start = start_long;
    lat_end = end_lat;
    long_end = end_long;
    load('measurement_goals.mat')
    goals_lat = lat_goals;
    goals_long = long_goals;
    goals_rad = rad_goals;
    goals_val = val_goals;
    goals_num = num_goals;
end
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

[lat_size_map,long_size_map] = size(latitudeMapWave);
val = 0;



    for step = 1:steps - 1
        seconds = t(step);
        hour = round(seconds/3600)*0 + 10;
        hour_real = seconds/3600 + 10;
        lat = latitude(step);
        long = longitude(step);
        error_map = sqrt((latitudeMapWave - lat).^2 + (longitudeMapWave - long).^2);
        [xx,yy] = find(error_map == min(error_map, [], 'all'));
        error_map = sqrt((latitudeCurrentMap - lat).^2 + (longitudeCurrentMap - long).^2);
        [xcurrent,ycurrent] = find(error_map == min(error_map, [], 'all'));

        % Find position on weather map
%         long_lat = repmat(lat,lat_size_map,long_size_map);
%         long_long = repmat(long,lat_size_map,long_size_map);

        % Find Mission goal Values. 
        goal_dist = Haversine_deg(goals_long, goals_lat, long_end,lat_end,r);
        inside  = (goal_dist<goals_rad)*sensor(step);
        val = val + sum(inside.*theGoals);
        theGoals = theGoals.*(~inside);

        course = bearing(latitude(step), longitude(step), latitude(step+1), longitude(step+1));
        if length(course) > 1
            course
        end
        
        % - Get current
    %     curr_dir = current_dir(index(1),index(2),hour)*deg2rad;
    %     curr_spd = current_speed(index(1),index(2),hour);
        % - Get wind
        wnd_dir = windDir(xx(1),yy(1),hour);
        wnd_spd = windSpeed(xx(1),yy(1),hour);
        % - Get waves
        ww_dir = waveDir(xx(1),yy(1),hour);
        ww_spd = waveSize(xx(1),yy(1),hour);
        curWaveHz = waveHZ(xx(1),yy(1),hour);
        currentDir = atan2(currentEast(xcurrent(1),ycurrent(1),hour), currentNorth(xcurrent(1),ycurrent(1),hour));
        currentSpeed = norm(currentEast(xcurrent(1),ycurrent(1),hour), currentNorth(xcurrent(1),ycurrent(1),hour));
        % % %     waves = [ww_spd*cos(ww_dir); ww_spd*sin(ww_dir)];
        % - Describe full function (possibly newton method??
        testSum = ww_spd + wnd_spd;
        relWaveDir = ssa(rad2deg(course) - ww_dir - 180, 'deg');
%         if isnan(testSum)
%             t(steps) = inf;
%             
%             break
%         end
        if length(currentSpeed) >1
            currentSpeed
        elseif length(rad2deg(ssa(currentDir-course)))> 1
            rad2deg(ssa(currentDir-course))
        elseif length(relWaveDir)> 1
            xx
            yy
            relWaveDir
            course
            rad2deg(course)
            ww_dir
        elseif length(ssa(wnd_dir - rad2deg(course),'deg'))> 1
            ssa(wnd_dir - rad2deg(course),'deg')
        elseif length(wnd_spd)> 1
            wnd_spd
            hour
        end
            
        X = [currentSpeed rad2deg(ssa(currentDir-course)) relWaveDir ssa(wnd_dir - rad2deg(course),'deg') ...
            wnd_spd(1) ww_spd curWaveHz];
        v(step) = predict(Mdl1, X);
        distance = Haversine_deg(latitude(step+1),longitude(step+1),latitude(step),longitude(step),r);
        if (v(step) < 0)||(~isreal(v(step)))
            t(steps) = inf;
            disp('halla')
            break
        end
        t(step+1) = t(step) + distance/v(step);

        % Calculating Battery energy
    end

    time = t(steps)/3600 - val;
end


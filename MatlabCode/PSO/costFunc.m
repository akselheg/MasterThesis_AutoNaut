
function elapseTime = costFunc(x)
    persistent curcurrentEast curcurrentNorth curwindDir curwindSpeed ...
        curwaveSize curwaveDir curlatitudeMapWave curlatitudeCurrentMap ...
        curlongitudeCurrentMap curlongitudeMapWave curwaveHZ
    persistent startLat endLat startLon endLon
    persistent GpsFix
    persistent Mdl1 Mdl2 Mdl3 Mdl4 Mdl5 Mdl6 w1
    if isempty(Mdl1)
        load 'speedEst2.mat' 
        load 'w1.mat'
    end
    if isempty(GpsFix)
        load('../Mausund200703_132548/GpsFix.mat')
    end
    r = 6371000;
    if isempty(curlatitudeMapWave)
        load('../Weather/weatherData_2020-7-3_2020-7-4.mat')
        load('../Weather/currentweatherData_2020-7-3_2020-7-4.mat')
        load('../Mausund200703_132548/GpsFix.mat')
        curcurrentEast =  currentEast;
        curcurrentNorth = currentNorth;
        curwindDir = windDir;
        curwindSpeed = windSpeed;
        curwaveSize = waveSize; 
        curwaveDir = waveDir;
        curlatitudeMapWave = latitudeMapWave;
        curlatitudeCurrentMap = latitudeCurrentMap;
        curlongitudeCurrentMap = longitudeCurrentMap;
        curlongitudeMapWave = longitudeMapWave;
        curwaveHZ = waveHZ;
    end
    if isempty(startLat)
        load 'Waypoints.mat';
    end

    %%
    inlat_points = x(1:length(x)/2);
    inlong_points = x(length(x)/2 + 1:end);
    lat_points = [startLat inlat_points endLat];
    long_points = [startLon  inlong_points endLon];
    
    time = GpsFix.utc_time(1);
    hours = floor(double(GpsFix.utc_time(1))/3600);
    n = length(lat_points);
    if hours > 47, hours = 47; end
    for i=1:n - 1
        angle = bearing(lat_points(i),long_points(i),lat_points(i+1),long_points(i+1));
        dist = Haversine_deg(lat_points(i),long_points(i),lat_points(i+1),long_points(i+1),r);
        error_map = sqrt((curlatitudeMapWave - lat_points(i)).^2 + (curlongitudeMapWave - long_points(i)).^2);
        [x,y] = find(error_map == min(error_map, [], 'all'));
        error_map = sqrt((curlatitudeCurrentMap - lat_points(i)).^2 + (curlongitudeCurrentMap - long_points(i)).^2);
        [xcurrent,ycurrent] = find(error_map == min(error_map, [], 'all'));
        cWaveSize = curwaveSize(x(1),y(1),hours+1);
        cWaveDir = ssa(rad2deg(angle) - curwaveDir(x(1),y(1),hours+1), 'deg');
        cWaveHz = 2*pi/curwaveHZ(x(1),y(1),hours+1);
        cWindDir = ssa(curwindDir(x(1),y(1),hours+1) - rad2deg(angle),'deg') ;
        cWindSpeed = curwindSpeed(x(1),y(1),hours+1);
        cCurrentEast = curcurrentEast(xcurrent(1),ycurrent(1),hours+1);
        cCurrentNorth = curcurrentNorth(xcurrent(1),ycurrent(1),hours+1);
        
        if ~ (isnan(cWaveSize(1)) || isnan(cCurrentEast(1)))
            Vc = [cCurrentNorth; cCurrentEast];
            VcDir = atan2d( Vc(2), Vc(1));
            
            X = [norm(Vc) ssa(rad2deg(angle) - VcDir,'deg') cWaveDir cWindDir ...
                     cWindSpeed cWaveSize cWaveHz];
            out2 = predict(Mdl2, X);
            out3 = predict(Mdl3, X);
            out4 = predict(Mdl4, X);
            out5 = predict(Mdl5, X);
            out6 = predict(Mdl6, X);
            speed = (out2 + out3 + out4 + out5 + out6)/5;
%             X = [cWaveSize cWaveHz abs(cos(deg2rad(cWaveDir)))  cWindSpeed*cos(cWindDir)  ...
%                 norm(Vc)* cos(ssa(angle - deg2rad(VcDir))) 1];
%             speed = w1'*X';
        else
            disp('Not good')
        end
        
        time = time + dist/speed;
        hours = floor(time/3600);
        if hours > 47, hours = 47; end
    end
    %%
    elapseTime = time-GpsFix.utc_time(1);
end
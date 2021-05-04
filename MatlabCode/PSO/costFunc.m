
function elapseTime = costFunc(x)
    persistent curcurrentEast curcurrentNorth curwindDir curwindSpeed ...
        curwaveSize curwaveDir curlatitudeMapWave curlatitudeCurrentMap ...
        curlongitudeCurrentMap curlongitudeMapWave curwaveHZ
    persistent startLat endLat startLon endLon
    persistent GpsFix
    persistent Mdl1
    if isempty(Mdl1)
        load 'speedEst.mat' 
    end
    if isempty(GpsFix)
        load('../Mausund200703_132548/GpsFix.mat')
    end
    r = 6371000;
    if isempty(curlatitudeMapWave)
        'hei'
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
    avrspeed = 0;
    
    %curhour = 0;
    for i=1:n - 1

        angle = bearing(lat_points(i),long_points(i),lat_points(i+1),long_points(i+1));
        dist = Haversine_deg(lat_points(i),long_points(i),lat_points(i+1),long_points(i+1),r);
        %%
        error_map = sqrt((curlatitudeMapWave - lat_points(i)).^2 + (curlongitudeMapWave - long_points(i)).^2);
        [x,y] = find(error_map == min(error_map, [], 'all'));
        error_map = sqrt((curlatitudeCurrentMap - lat_points(i)).^2 + (curlongitudeCurrentMap - long_points(i)).^2);
        [xcurrent,ycurrent] = find(error_map == min(error_map, [], 'all'));
        cWaveSize = curwaveSize(x(1),y(1),hours+1);
        cWaveDir = ssa(rad2deg(angle) - curwaveDir(x(1),y(1),hours+1) - 180, 'deg');
        cWaveHz = curwaveHZ(x(1),y(1),hours+1);
        cWindDir = ssa(curwindDir(x(1),y(1),hours+1) - rad2deg(angle),'deg') ;
        cWindSpeed = curwindSpeed(x(1),y(1),hours+1);
        cCurrentEast = curcurrentEast(xcurrent(1),ycurrent(1),hours+1);
        cCurrentNorth = curcurrentNorth(xcurrent(1),ycurrent(1),hours+1);
        if ~ (isnan(cWaveSize(1)) || isnan(cCurrentEast(1)))

            Vc = [cCurrentNorth; cCurrentEast];
            VcDir = atan2( Vc(2), Vc(1));

            X = [norm(Vc) rad2deg(ssa(VcDir-angle)) cWaveDir cWindDir  ...
                    cWindSpeed cWaveSize cWaveHz];

            %disp(num2str((dist)));
            speed = predict(Mdl1, X);
        end
        avrspeed = avrspeed + speed;
        time = time + dist/speed;
        %hours = floor(time/3600)
    end
    %%
    elapseTime = time-GpsFix.utc_time(1);
end
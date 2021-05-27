%% Clear workspace
clearvars; close all; clc;

%% Load Validation mission
load('../Mausund200703_132548/GpsFix.mat')

%% Set start and goal positions and number of waypoints between
startLat = rad2deg(GpsFix.lat(1));
endLat = rad2deg(GpsFix.lat(end));
startLon = rad2deg(GpsFix.lon(1));
endLon = rad2deg(GpsFix.lon(end));
nWaypoints = 15;
save('Waypoints.mat','startLat','endLat','startLon','endLon');

%% Create the Guess (Straight line
n = nWaypoints + 2;
[lat_points, long_points] = pathCreator([startLat, endLat]', [startLon, endLon]', n);
lat_points(2:end-1) = lat_points(2:end-1);
long_points(2:end-1) = long_points (2:end-1);
straight = costFunc([lat_points(2:end-1); long_points(2:end-1)]');

%% constrainting the  problem
latmin = min(lat_points);
lonmin = min(long_points);
latmax = max(lat_points);
lonmax = max(long_points);
l = max(latmax-latmin,lonmax-lonmin);

%% PSO size
nParticle = 500; % 500
MaxIter = 200; % 200

%% seed and generate random paths and option for PSO
seed = [lat_points(2:end-1); long_points(2:end-1)];
x0 = (repmat(seed,1,nParticle) + l/2*[randn(nWaypoints*2,nParticle)])';
option = optimoptions('particleswarm','SwarmSize',nParticle, 'Display', 'iter', ...
    'MaxIterations', MaxIter, 'InitialSwarmMatrix', x0,...
    'PlotFcn','pswplotbestf', 'InitialSwarmSpan', seed');

%% Run PSO
best = inf;
for i =  1:4 % Running 4 itterations to obtain best solution.
    
    x0 = (repmat(seed,1,nParticle) + l/2*[randn(nWaypoints*2,nParticle)])';
    option = optimoptions('particleswarm','SwarmSize',nParticle, 'Display', 'iter', ...
        'MaxIterations', MaxIter, 'InitialSwarmMatrix', x0,...
        'PlotFcn','pswplotbestf', 'InitialSwarmSpan', seed');
    x = particleswarm(@costFunc,2*(nWaypoints),[latmin-l*ones(1,(nWaypoints)) lonmin-l*ones(1,(nWaypoints))],...
        [latmax+l*ones(1,(nWaypoints)) lonmax+l*ones(1,(nWaypoints))], option);
    
    act = GpsFix.utc_time(end) - GpsFix.utc_time(1);
    
    Optimal = costFunc(x);
    if Optimal < best
        best = Optimal;
        Xopt = x;
    end

    %% Plot this itteration
    result_lat = [startLat x(1:(nWaypoints)) endLat];
    result_long = [startLon x((nWaypoints) +1:end) endLon];
    figure;
    geoplot(result_lat,result_long)
    geolimits([latmin-l/2 latmax+l/2],[lonmin-l/2 lonmax+(l/2)])
end

%% Plot straight path
figure;
geoplot(lat_points,long_points)
geolimits([latmin-(l/2) latmax+(l/2)],[lonmin-(l/2) lonmax+(l/2)])

%% Plot AutoNaut Path
figure;
geoplot(rad2deg(GpsFix.lat(2:50:end-1)),rad2deg(GpsFix.lon(2:50:end-1)))
geolimits([latmin-l/2 latmax+l/2],[lonmin-l/2 lonmax+(l/2)])

%% Plot best results
result_lat = [startLat Xopt(1:(nWaypoints)) endLat];
result_long = [startLon Xopt((nWaypoints) +1:end) endLon];
figure;
geoplot(result_lat,result_long)
geolimits([latmin-l/2 latmax+l/2],[lonmin-l/2 lonmax+(l/2)])
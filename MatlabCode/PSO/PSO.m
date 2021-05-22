clearvars; close all; clc;

load('../Mausund200703_132548/GpsFix.mat')
startLat = rad2deg(GpsFix.lat(1));
endLat = rad2deg(GpsFix.lat(end));
startLon = rad2deg(GpsFix.lon(1));
endLon = rad2deg(GpsFix.lon(end));
dist = Haversine_deg(startLat,startLon,endLat,endLon,6371000);
nWaypoints = 15;
n = nWaypoints + 2;
save('Waypoints.mat','startLat','endLat','startLon','endLon');
[lat_points, long_points] = pathCreator([startLat, endLat]', [startLon, endLon]', n);
lat_points(2:end-1) = lat_points(2:end-1);% + 0.05*rand(n-2,1);
long_points(2:end-1) = long_points (2:end-1);% + 0.05*rand(n-2,1);
straight = costFunc([lat_points(2:end-1); long_points(2:end-1)]');
latmin = min(lat_points);
lonmin = min(long_points);
latmax = max(lat_points);
lonmax = max(long_points);
l = max(latmax-latmin,lonmax-lonmin);
%tic;
costFunc([lat_points(2:end-1); long_points(2:end-1)]');
%toc;
nParticle = 500; % 500
MaxIter = 200; % 200
seed = [lat_points(2:end-1); long_points(2:end-1)];
x0 = (repmat(seed,1,nParticle) + l/2*[randn(nWaypoints*2,nParticle)])';
%seed = [lat_points(2:end-1)',long_points(2:end-1)', ones(1,n-2)]';
option = optimoptions('particleswarm','SwarmSize',nParticle, 'Display', 'iter', ...
    'MaxIterations', MaxIter, 'InitialSwarmMatrix', x0,...
    'PlotFcn','pswplotbestf', 'InitialSwarmSpan', seed');
best = inf;

for i =  1:4
    x0 = (repmat(seed,1,nParticle) + l/2*[randn(nWaypoints*2,nParticle)])';
    %seed = [lat_points(2:end-1)',long_points(2:end-1)', ones(1,n-2)]';
    option = optimoptions('particleswarm','SwarmSize',nParticle, 'Display', 'iter', ...
        'MaxIterations', MaxIter, 'InitialSwarmMatrix', x0,...
        'PlotFcn','pswplotbestf', 'InitialSwarmSpan', seed');
    tic;
    
    x = particleswarm(@costFunc,2*(nWaypoints),[latmin-l*ones(1,(nWaypoints)) lonmin-l*ones(1,(nWaypoints))],...
        [latmax+l*ones(1,(nWaypoints)) lonmax+l*ones(1,(nWaypoints))], option);
    toc;
    %pred = costFunc([rad2deg(GpsFix.lat(2:50:end-1)); rad2deg(GpsFix.lon(2:50:end-1))]');
    act = GpsFix.utc_time(end) - GpsFix.utc_time(1);
    Optimal = costFunc(x);
    if Optimal < best
        best = Optimal;
        Xopt = x;
    end
    straight = costFunc([lat_points(2:end-1); long_points(2:end-1)]');

    %%
    [inputs, ~] = size(x);
    result_lat = [startLat x(1:(nWaypoints)) endLat];
    result_long = [startLon x((nWaypoints) +1:end) endLon];
    %%
    %%
    figure;
    geoplot(result_lat,result_long)
    geolimits([latmin-l/2 latmax+l/2],[lonmin-l/2 lonmax+(l/2)])
end
figure;
geoplot(lat_points,long_points)
geolimits([latmin-(l/2) latmax+(l/2)],[lonmin-(l/2) lonmax+(l/2)])
%%
figure;
PSOgeoplot(rad2deg(GpsFix.lat(2:50:end-1)),rad2deg(GpsFix.lon(2:50:end-1)))
 geolimits([latmin-l/2 latmax+l/2],[lonmin-l/2 lonmax+(l/2)])
%%
[inputs, ~] = size(Xopt);
result_lat = [startLat Xopt(1:(nWaypoints)) endLat];
result_long = [startLon Xopt((nWaypoints) +1:end) endLon];

figure;
geoplot(result_lat,result_long)
 geolimits([latmin-l/2 latmax+l/2],[lonmin-l/2 lonmax+(l/2)])
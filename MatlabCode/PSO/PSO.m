clear all; close all; clc;

load('../Mausund200703_132548/GpsFix.mat')
startLat = 64;
endLat = 66.5;
startLon = 8;
endLon = 10;
dist = Haversine_deg(startLat,startLon,endLat,endLon,6371000);
n = 15;
save('Waypoints.mat','startLat','endLat','startLon','endLon');
[lat_points, long_points] = pathCreator([startLat, endLat]', [startLon, endLon]', n);
lat_points(2:end-1) = lat_points(2:end-1);% + 0.05*rand(n-2,1);
long_points(2:end-1) = long_points (2:end-1);% + 0.05*rand(n-2,1);
costFunc([lat_points(2:end-1); long_points(2:end-1)]');
tic;
default = costFunc([lat_points(2:end-1); long_points(2:end-1)]');
toc;
latmin = min(lat_points);
lonmin = min(long_points);
latmax = max(lat_points);
lonmax = max(long_points);
nParticles = 10;
maxIter = 10;
contInput = 2*(n-2);
boolInput = n-2;
seed = [lat_points(2:end-1)',long_points(2:end-1)', ones(1,n-2)]';
option = optimoptions('particleswarm','SwarmSize',500, 'Display', 'iter', 'MaxIterations', 200);
tic;
x = particleswarm(@costFunc,2*(n-2),[lat_points(2:end-1)' lonmin-((lonmax -lonmin)/2)*ones(1,(n-2))],...
    [lat_points(2:end-1)' lonmax+((lonmax-lonmin))/2*ones(1,(n-2))], option);
toc;


%[x_2,x_val_2] = hybrid_PSO(@costFunc, contInput,boolInput,seed,nParticles, maxIter, GpsFix.utc_time(1) );
%%
[inputs, ~] = size(x);
result_lat = [startLat x(1:(n-2)) endLat];
result_long = [startLon x((n-2) +1:end) endLon];
%%
figure(1)
plot([startLat result_lat endLat],[startLon result_long endLon])
% %%
% figure(2)
% X = repmat((1:maxIter)',1,nParticles+1);
% defaultMtx = repmat(default,maxIter,1);
% plot(X,[x_val_2,defaultMtx]);
%%
figure(3)
geoplot(result_lat,result_long)
geolimits([55 70],[-5 10])
figure(4)
geoplot(lat_points,long_points)
geolimits([55 70],[-5 10])

clc;clearvars;close all;


startx = 60; starty = 6;
endx = 60.2; endy = 6.2;
radius=6371;
lat1=startx*pi/180;
lat2=endx*pi/180;
lon1=starty*pi/180;
lon2=endy*pi/180;
deltaLat=lat2-lat1;
deltaLon=lon2-lon1;
a=sin((deltaLat)/2)^2 + cos(lat1)*cos(lat2) * sin(deltaLon/2)^2;
c=2*atan2(sqrt(a),sqrt(1-a));
d1km=radius*c;    %Haversine distance
x=deltaLon*cos((lat1+lat2)/2);
y=deltaLat;
d2km=radius*sqrt(x*x + y*y); %Pythagoran distance
    
plot([startx endx], [starty endy])

psi = rad2deg(atan2(endy -starty, endx-startx));
speed = 0.5;

function d2km = latlon2km(lat1,lon1, lat2,lon2)
    lat1=lat1*pi/180;
    lat2=lat2*pi/180;
    lon1=lon1*pi/180;
    lon2=lon2*pi/180;
    deltaLat=lat2-lat1;
    deltaLon=lon2-lon1;
    a=sin((deltaLat)/2)^2 + cos(lat1)*cos(lat2) * sin(deltaLon/2)^2;
    c=2*atan2(sqrt(a),sqrt(1-a));
    d1km=radius*c;    %Haversine distance
    x=deltaLon*cos((lat1+lat2)/2);
    y=deltaLat;
    d2km=radius*sqrt(x*x + y*y); %Pythagoran distance
end
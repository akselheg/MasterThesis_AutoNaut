
clc;clearvars;
[lat lon] = TryDownloadWeather(1,1,1,1,1,1);

function [lat lon] = TryDownloadWeather(x0,y0,x1,y1, day,time)
    waveData = "https://thredds.met.no/thredds/dodsC/fou-hi/mywavewam800vhf/mywavewam800_vestlandet.an.2019052818.nc";
    %waveInfo = ncinfo(waveData);    
    lat = ncread(waveData,'latitude');
    lon = ncread(waveData,'longitude');
end
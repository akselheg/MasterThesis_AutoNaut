
clc;clearvars;close all;
r = 6371000;

startx = 60; starty = 6;
endx = 60.2; endy = 6.2;
distance = Haversine_deg(startx, starty, endx, endy,r);

for i startz
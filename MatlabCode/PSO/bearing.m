function [bearing] = bearing(startLat,startLong,endLat,endLong)
startLat = (pi/180)*startLat;
endLat = (pi/180)*endLat;
startLong = (pi/180)*startLong;
endLong = (pi/180)*endLong;
deltaLong = endLong - startLong;
bearing = atan2(sin(deltaLong)*cos(endLat) , cos(startLat)*sin(endLat) - sin(startLat)*cos(endLat)*cos(deltaLong));
end
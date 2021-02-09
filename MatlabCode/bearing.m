function [bearing] = bearing(startLat,startLong,endLat,endLong)
deltaLong = endLong - startLong;
bearing = atan2(sin(deltaLong)*cos(endLat) , cos(startLat)*sin(endLat) - sin(startLat)*cos(endLat)*cos(deltaLong));
end
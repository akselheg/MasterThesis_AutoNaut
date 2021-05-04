function [distance] = Haversine_deg(lat_start,long_start,lat_end,long_end,radius)
    lat_start = (pi/180)*lat_start;
    lat_end = (pi/180)*lat_end;
    long_start = (pi/180)*long_start;
    long_end = (pi/180)*long_end;
    
    distance = 2*radius*asin(sqrt(sin((lat_end - lat_start)/2).^2 + cos(lat_start).*cos(lat_end).*sin((long_end - long_start)/2).^2));
    
end


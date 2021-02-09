function [lat_points, long_points] = pathCreator(goals_lat, goals_long, n_points)
%PATHCREATOR - Adds a set of equally spaced points between the goals to
%describe a path that intersects all goals including by a first order
%splined path
% UNFINISHED!!!! TODO: Fix lost points

n_goals = size(goals_lat,1);
lines = n_goals -1;

new_points = n_points - n_goals;
distances = Haversine_deg(goals_lat(1:end-1),goals_long(1:end-1),goals_lat(2:end),goals_long(2:end),6371*10^3); 
tot_distance = sum(distances,'all');

distance_sets = tot_distance/(new_points +1);

distance_set = tril(ones(new_points,new_points))*repmat(distance_sets,new_points,1);
goals_distance= tril(ones(lines,lines))*distances;

[X,Y] = meshgrid(distance_set,goals_distance);

point_partition_sum= sum((X < Y),2);

difference_mtx = diag(-ones(lines-1,1),-1) + diag(ones(lines,1));
point_partition = difference_mtx*point_partition_sum;

lat_points  = [];
long_points = [];

    for i = 1:(n_goals-1)

        lat_points = [lat_points;goals_lat(i)];
        long_points = [long_points;goals_long(i)];


        intermediate_lat = (1:point_partition(i))*(goals_lat(i+1) - goals_lat(i))/(point_partition(i) +1) + goals_lat(i);
        intermediate_long = (1:point_partition(i))*(goals_long(i+1) - goals_long(i))/(point_partition(i) +1) + goals_long(i);

        lat_points = [lat_points; intermediate_lat'];
        long_points = [long_points; intermediate_long'];



    end
    
lat_points = [lat_points; goals_lat(n_goals)];
long_points = [long_points; goals_long(n_goals)];
end


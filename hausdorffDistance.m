function dH = hausdorffDistance(points, index_set_one, index_set_two)
    % points is a n*3 matrix, n is the number of points
    % index_set_one and index_set_two are the index of the points in the
    % two sets
    % dH is the hausdorff distance between the two sets
    
    % get the points in the two sets
    set_one = points(index_set_one, :);
    set_two = points(index_set_two, :);
    
    % get the distance matrix
    distance_matrix = pdist2(set_one, set_two);
    
    % get the hausdorff distance
    dH = max(min(distance_matrix, [], 2));
end
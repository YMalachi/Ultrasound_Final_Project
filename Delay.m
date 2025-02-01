function [time_delay_vec] = Delay(element_number, element_distances, c, focus_point)
    % focus point is a vector that looks like that [x,0,z]
    time_delay_vec = zeros(1,element_number);
    x_co = focus_point(1);
    y_co = focus_point(2);
    z_co = focus_point(3);
    dist_to_focus_from_center = sqrt(x_co^2 + z_co^2);
    for curr_ele = 1:element_number % assuming there are always an even number of elements, for symmetry
        dist_to_focus_from_curr_ele = sqrt((element_distances(curr_ele)-x_co)^2+z_co^2);
        delay_ele = (dist_to_focus_from_center - dist_to_focus_from_curr_ele)/c; % c is [m/sec]
        time_delay_vec(curr_ele) = delay_ele;
    end
    min_time = min(time_delay_vec);
    time_delay_vec = time_delay_vec ;
end
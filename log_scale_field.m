function [log_scaled_field] = log_scale_field(field)
    normalized_pressure_field = [];
    for idx = 1:length(field(1,:))
        normalized_pressure_field(idx) = norm(field(:,idx));
    end
    scaled_pressure_field = reshape(normalized_pressure_field,[],100); % scaling
    min_scaled_pressure_field = min(scaled_pressure_field,[], "all");
    max_scaled_pressure_field = max(scaled_pressure_field,[], "all");
    delta_max_min = max_scaled_pressure_field - min_scaled_pressure_field;
    linear_norm_scaled_pressure_field = (scaled_pressure_field - min_scaled_pressure_field)/(delta_max_min);
    
    % log scale
    a = (10^-1.5);
    b = 1-a;
    log_scaled_field = transpose(20*log10(b*linear_norm_scaled_pressure_field+a));
end
% phantom size
x_size = transpose(linspace(-15/1000, 15/1000, 100)); 
z_size = transpose(linspace(0/1000, 50/1000, 100));
y_size = zeros(100,1);
[x, y, z] = meshgrid(x_size, 0, z_size);
phantom_us = zeros(100);
phantom_us(80,24) = 256; % this is [-8 0 40]
phantom_us(80,77) = 256;

x_lim = (x_size(end)-x_size(1))+x_size(1);
z_lim = (z_size(end)-z_size(1))+z_size(1);
[~,x_rand_pos1] = min(abs(rand(1,1)*x_lim-x_size));
[~,z_rand_pos1] = min(abs(rand(1,1)*z_lim-z_size));
[~,x_rand_pos2] = min(abs(rand(1,1)*x_lim-x_size));
[~,z_rand_pos2] = min(abs(rand(1,1)*z_lim-z_size));

phantom_us(x_rand_pos1, z_rand_pos1) = 100;
phantom_us(x_rand_pos2, z_rand_pos2) = 100;

if plot_flag
    % show the phantom. note that the random scaterers will not be the same
    % (obviously)
    figure;
    imagesc(x_size*1000, z_size*1000, phantom_us);
    colormap("gray");
end

scaterers_mat = [
    80, 24;
    80, 77;
x_rand_pos1, z_rand_pos1;
x_rand_pos2, z_rand_pos2
                ];
amp_vec = [zeros(100,1);256;256;10;10];
points_x = [x_size ; scaterers_mat(:,1)];
points_y = [y_size ; 0 ; 0 ; 0 ; 0];
points_z = [z_size ; scaterers_mat(:,2)];
points = [points_x, points_y, points_z];

scan_echo = zeros(1,100);
% simulating 100 line scan
for line = 1:100
    if line <= 50
        xdc_focus(linear_128_trans, 0, [x_field_vector(line), 0, 40/1000]);
        xdc_apodization(linear_128_trans,0,left_half_apod);
        fprintf("Processing:%.1f", line);
        line_echo = calc_scat(linear_128_trans, linear_128_trans_r,points,amp_vec);
        scan_echo(1:length(line_echo),line) = line_echo;
    end
    if line > 50
        xdc_focus(linear_128_trans, 0, [x_field_vector(line), 0, 40/1000]);
        xdc_apodization(linear_128_trans,0,right_half_apod);
        fprintf("Processing:%.d\n", line);
        line_echo = calc_scat(linear_128_trans, linear_128_trans_r,scaterers_mat,amp_vec);
        scan_echo(1:length(line_echo),line) = line_echo;
    end
end


% hilbert transform
rf_data = zeros(size(scan_echo)); 
for line = 1:size(scan_echo, 2) 
    rf_data(:, line) = hilbert(scan_echo(:, line));
    fprintf("Hilberting... %.d\n",line)
end
envelope_rf = abs(rf_data);
log_scaled_envelope_rf = log_scale_field(rf_data);

% plot
if plot_flag
    figure;
    imagesc(x_field_vector*1000, z_field_vector*1000, log_scaled_envelope_rf);
    title('Phantom Pressure Field Image, Log Scaled, 100 Scan Lines', FontSize=13);
    xlabel('X [mm]');
    ylabel('Z [mm]');
    colormap("gray");
    colorbar;
end

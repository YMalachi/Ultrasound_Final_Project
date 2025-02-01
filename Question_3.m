clear; clc; close all;
%% flags
plot_flag = 1; 

field_init(0);

%% params
% initiating the transducer and params
linear_128_trans = xdc_linear_array(128,0.1/1000,1/1000,0.05/1000,1,1,[0 0 80/1000]);
linear_128_trans_r = xdc_linear_array(128,0.1/1000,1/1000,0.05/1000,1,1,[0 0 80/1000]);
f0 = 3e6; % main frequency
t0 = 1/f0; % period time
fs = 100e6; % sampling rate
ts = 1/fs; % sampling time
set_sampling(fs);
c = 1540; % m/sec, speed of sound in soft tissue

% initiating impulse
time_vec = 0:ts:t0; % for one cycle
excitation = sin(2*pi*time_vec*f0);
impulse = sin(2*pi*time_vec*f0);
% setting it to the transducers
xdc_excitation(linear_128_trans, excitation);
xdc_impulse(linear_128_trans, impulse);
xdc_excitation(linear_128_trans_r, excitation);
xdc_impulse(linear_128_trans_r, impulse);

%% 3a

right_half_apod = [zeros(1,64) ones(1,64)];
left_half_apod = [ones(1,64) zeros(1,64)];
odd_ele_apod = repmat([1, 0], 1, 64);
even_ele_apod = repmat([0, 1], 1, 64);

% calculating distances & delays
pitch = 0.15/1000;

% half and half distances
% the calculations are done this way because we need half of the
% transducer's delay to be calculated as if the middle is the point between
% the 32nd and 33rd element but it needs to have the right spatial
% dimensions for delay calculations as well... 
half_and_half_ele_vec = -64:1:64;
half_and_half_ele_vec(65) = [];
left_half_dist(1:32) = half_and_half_ele_vec(1:32) * pitch + 0.5*pitch;
left_half_dist(33:64) = half_and_half_ele_vec(33:64) * pitch - 0.5*pitch;
right_half_dist(1:32) = half_and_half_ele_vec(65:96) * pitch + 0.5*pitch;
right_half_dist(33:64) = half_and_half_ele_vec(97:128) * pitch - 0.5*pitch;
% calculating delays
delay_left = Delay(64, left_half_dist, c, [-8 0 40]/1000);
delay_right = Delay(64, right_half_dist, c, [8 0 40]/1000);
delay_vec_h_a_h = [delay_left delay_right];

% odd/even distances
odd_even_ele_vec_idx = -64:1:64;
odd_even_ele_vec_idx(65) = [];
odd_even_distances(1:64) = odd_even_ele_vec_idx(1:64) * pitch + 0.5*pitch;
odd_even_distances(65:128) = odd_even_ele_vec_idx(65:128) * pitch - 0.5*pitch;
odd_dist = odd_even_distances(1:2:128);
even_dist = odd_even_distances(2:2:128);
% calculating delays
odd_delay_vec = Delay(64, odd_dist, c, [-8 0 40]/1000);
even_delay_vec = Delay(64, even_dist, c, [8 0 40]/1000);
odd_even_delay = reshape([odd_delay_vec;even_delay_vec], 1, []);

% calculating pressure fields
final_pressure_field_h_and_h = set_and_combine_fields(linear_128_trans, linear_128_trans_r, left_half_apod, right_half_apod, delay_vec_h_a_h);
final_pressure_field_odd_even = set_and_combine_fields(linear_128_trans, linear_128_trans_r, odd_ele_apod, even_ele_apod, odd_even_delay);

% plotting
x_field_vector = linspace(-15/1000, 15/1000, 100);
z_field_vector = linspace(5/1000,50/1000, 100);
log_pressure_field_h_and_h = log_scale_field(final_pressure_field_h_and_h);
log_pressure_field_odd_even = log_scale_field(final_pressure_field_odd_even);
% lateral cut
z_lateral = 40/1000; % m
z_lateral_idx = find(z_field_vector == z_lateral); 
z_lateral_vector_h_a_h = log_pressure_field_h_and_h(z_lateral_idx,:) + 30;
z_lateral_vector_odd_even = log_pressure_field_odd_even(z_lateral_idx,:) + 30;

if plot_flag
    figure;
    subplot(1,2,1);
    imagesc(x_field_vector*1000, z_field_vector*1000, log_pressure_field_h_and_h);
    colormap('hot');
    colorbar;
    title('Pressure Field "Half & Half" Apodization', FontSize=13);
    xlabel('X [mm]');
    ylabel('Z [mm]');

    subplot(1,2,2);
    imagesc(x_field_vector*1000, z_field_vector*1000, log_pressure_field_odd_even);
    colormap('hot');
    colorbar;
    title('Pressure Field "Odd-Even" Apodization', FontSize=13);
    xlabel('X [mm]');
    ylabel('Z [mm]');
    
    % lateral cut plot
    figure;
    subplot(1,2,1);
    plot(x_field_vector, z_lateral_vector_h_a_h);
    title('Lateral Cut of Pressure Field at 40mm Depth, "Half & Half" Apodization', FontSize=13);
    xlabel('X [mm]');
    ylabel('Amplitude [dB]');

    subplot(1,2,2);
    plot(x_field_vector, z_lateral_vector_odd_even);
    title('Lateral Cut of Pressure Field at 40mm Depth, "Odd-Even" Apodization', FontSize=13);
    xlabel('X [mm]');
    ylabel('Amplitude [dB]');
end


%% 3b

% creating phantom 
% Define the parameters
num_points = 98;  % Number of random points
x_lim = [-15/1000 15/1000]; % Range of x coordinates
z_lim = [0/1000 50/1000]; % Range of z coordinates
x_axis = linspace(x_lim(1), x_lim(2), 100);
z_axis = linspace(z_lim(1), z_lim(2), 100);
scatter_idx = {[8/1000 0 40/1000; -8/1000 0 40/1000], [10; 10]}; % Scatter points and amplitudes
num_lines = 100; % Number of scanning lines
x_scan = transpose(linspace(x_lim(1), x_lim(2), num_lines)); % X coordinates of scanning lines
y_scan = zeros(num_lines,1);
z_scan = ones(num_lines,1) * 40/1000;
focus_point = [x_scan, y_scan, z_scan]; % Focal points for each scanning line

% Generate random points within the specified range
random_x = (x_lim(2) - x_lim(1)) * rand(num_points, 1) + x_lim(1);
random_z = (z_lim(2) - z_lim(1)) * rand(num_points, 1) + z_lim(1);
random_y = zeros(num_points, 1); % Set y coordinate to 0 for all points
random_amplitudes = zeros(num_points, 1); % Set amplitudes to 0 for all points

% Concatenate scatter points with random points
points_x = [random_x; scatter_idx{1}(:, 1)];
points_y = [random_y; scatter_idx{1}(:, 2)];
points_z = [random_z; scatter_idx{1}(:, 3)];
points = [points_x, points_y, points_z];
amplitudes = [random_amplitudes; scatter_idx{2}];

image_data=zeros(1,num_lines);
% Loop through each scanning line
for line = 1:num_lines
    if line <= 50
        xdc_apodization(linear_128_trans,0,left_half_apod);
        xdc_apodization(linear_128_trans_r,0,left_half_apod);
        xdc_focus(linear_128_trans, 0, focus_point(line,:));
        xdc_focus(linear_128_trans_r,0,focus_point(line,:));
        v = calc_scat(linear_128_trans, linear_128_trans_r, points, amplitudes);
        fprintf("processing %.d\n",line);
        % Store the result
        image_data(1:max(size(v)),line)=v;
    end
    if line > 50
        xdc_apodization(linear_128_trans,0,right_half_apod);
        xdc_apodization(linear_128_trans_r,0,right_half_apod);
        xdc_focus(linear_128_trans, 0, focus_point(line,:));
        xdc_focus(linear_128_trans_r,0,focus_point(line,:));
        v = calc_scat(linear_128_trans, linear_128_trans_r, points, amplitudes);
        fprintf("processing %.d\n",line);
        % Store the result
        image_data(1:max(size(v)),line)=v;
    end
end

% Apply Hilbert transform to calculate envelope
hilb_data = hilbert(image_data);
envelope_data = abs(hilb_data);
if image_data == hilb_data
    disp("hilbert didnt work");
end

% Normalize the envelope data
envelope_data_norm=(envelope_data-min(envelope_data,[],'all'))/(max(envelope_data,[],'all')-min(envelope_data,[],'all'));

% Convert to logarithmic scale of 30dB
a = 10^(-1.5);
b = 1-a;
envelope_data_log = 20*log10(a.*envelope_data_norm+b);

% Display the image
figure;
imagesc(x_axis*1000, z_axis*1000, envelope_data_log);
colormap(gray);
xlabel('X [mm]');
ylabel('Z [mm]');
title('Phantom Image - 100 Scanning Lines', FontSize=13);
colorbar;


%% local functions
function [final_pressure_field] = set_and_combine_fields(trans, trans_r, apod_1, apod_2, delay_vec)
    apodizations = {apod_1, apod_2};
    pressure_fields = {};
    for i = 1:2
        % initiating vectors for field matrix
        x_field_vector = linspace(-15/1000, 15/1000, 100);
        z_field_vector = linspace(5/1000 ,50/1000, 100);
        [x, y, z] = meshgrid(x_field_vector, 0, z_field_vector);
        points = [x(:) y(:) z(:)];

        % setting focus point
        xdc_focus_times(trans, 0, delay_vec);
        xdc_focus_times(trans_r, 0, delay_vec);
        xdc_apodization(trans, 0, apodizations{i});
        xdc_apodization(trans_r, 0, apodizations{i});

        % calculating the pressure field
        pressure_fields{i} = calc_hhp(trans, trans_r, points);
    end
    final_pressure_field = pressure_fields{1} + pressure_fields{2};
end


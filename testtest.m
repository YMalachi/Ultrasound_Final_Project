% Initialize Field
close all;
clear all;
clc;

addpath('C:\Users\noash\Downloads\Field_II_ver_3_22_windows');
field_init(0);

%% Parameters
f0 = 3e6;                      % Transducer center frequency [Hz]
fs = 100e6;                    % Sampling frequency [Hz]
c = 1482;                      % Speed of sound [m/s]
lambda = c / f0;               % Wavelength [m]
width = 0.1 / 1000;            % Width of element [m]
element_height = 1 / 1000;     % Height of element [m]
kerf = 0.05 / 1000;             % Kerf [m]
focus = [0, 0, 40] / 1000;     % Fixed focal point [m]
N_tx_elements = 128;           % Number of transmit elements
N_rx_elements = 128;           % Number of receive elements

% Define the impulse response
impulse_response = sin(2 * pi * f0 * (0:1/fs:2/f0)) .* hann(length(0:1/fs:2/f0))';

tx=xdc_linear_array(N_tx_elements,width,element_height,kerf,1,1,focus);
xdc_impulse(tx,impulse_response);
rx=xdc_linear_array(N_rx_elements,width,element_height,kerf,1,1,focus);
xdc_impulse(rx,impulse_response);


%% creating phantom

num_random_points = 100; % Number of random points (water)
scatter_amplitude = 10; % Amplitude of point scatterers
% Define the 2D phantom boundaries
x_min = -15 / 1000; % Minimum x-coordinate [m]
x_max = 15 / 1000;  % Maximum x-coordinate [m]
z_min = 0 / 1000;   % Minimum z-coordinate [m]
z_max = 50 / 1000;  % Maximum z-coordinate [m]

% Generate random water points
random_x = (x_max - x_min) * rand(num_random_points, 1) + x_min; % Random x-coordinates
random_z = (z_max - z_min) * rand(num_random_points, 1) + z_min; % Random z-coordinates
random_y = zeros(num_random_points, 1); % y-coordinates (always 0 for 2D phantom)
random_amplitudes = zeros(num_random_points, 1); % Amplitude for water points

% Define the scatterers
scatterers_x = [-8, 8] / 1000; % x-coordinates of scatterers [m]
scatterers_z = [40, 40] / 1000; % z-coordinates of scatterers [m]
scatterers_y = [0, 0]; % y-coordinates of scatterers [m]
scatterer_amplitudes = [scatter_amplitude; scatter_amplitude]; % Amplitude for scatterers

% Combine random points and scatterers into a single phantom
phantom_x = [random_x; scatterers_x'];
phantom_z = [random_z; scatterers_z'];
phantom_y = [random_y; scatterers_y'];
phantom_positions = [...
    [random_x; scatterers_x'], ... % x-coordinates
    [random_y; scatterers_y'], ... % y-coordinates
    [random_z; scatterers_z']];    % z-coordinates

phantom_amplitudes = [random_amplitudes; scatterer_amplitudes];% Display the phantom

% figure;
% scatter(phantom_x * 1000, phantom_z * 1000, 10, phantom_amplitudes, 'filled');
% colormap('jet');
% colorbar;
% xlabel('x [mm]');
% ylabel('z [mm]');
% title('2D Phantom');


%% Define Scanning Lines
num_scan_lines = 100; % Total number of scanning lines
num_elements = 128; % Total number of transducer elements
% Define Element Positions
element_positions = zeros(N_tx_elements, 1);
for i = 1:N_tx_elements
    element_positions(i) = (i - N_tx_elements / 2) * (kerf + width);
end
% Element positions (relative to array center)
element_positions = ((1:N_tx_elements) - (N_tx_elements + 1) / 2) * (width + kerf);


% Define Focus Points
left_focus = [-8, 0, 40] / 1000;   % Left focus (in meters)
right_focus = [8, 0, 40] / 1000;   % Right focus (in meters)

left_elements = 1:(num_elements / 2); % Left half of the transducer
%delays_left = Delay(length(left_elements), c, left_focus, element_positions(left_elements));

right_elements = (num_elements / 2 + 1):num_elements; % Right half of the transducer
%delays_right = Delay(length(right_elements), c, right_focus, element_positions(right_elements));

% Combine delays into full aperture delays
%delays_all = zeros(1, N_tx_elements);  % Single row, one column per element
%delays_all(left_elements) = delays_left;   % Apply delays for left group
%delays_all(right_elements) = delays_right; % Apply delays for right group

% Define focus times (static focus, valid from time 0)
times_all = 0;  % Single scalar for static focus

% Apply focus times and delays
%xdc_focus_times(tx, times_all, delays_all);

x_range = linspace(-15, 15, num_scan_lines) / 1000; % X coordinates of scanning lines in meters
z_range = linspace(-15, 15, num_scan_lines) / 1000; % X coordinates of scanning lines in meters
[X, Z] = meshgrid(x_range, z_range);
points = [X(:), zeros(numel(X), 1), Z(:)];

% Left group
weights_left = zeros(1, N_tx_elements);
weights_left(left_elements) = 1;  % Activate left elements only
xdc_apodization(tx, 0, weights_left);
%xdc_focus_times(tx, times_all, delays_all);  % Apply focus times
pressure_field_left = calc_hhp(tx,rx, points);


% Right group
weights_right = zeros(1, N_tx_elements);
weights_right(right_elements) = 1;  % Activate right elements only
xdc_apodization(tx, 0, weights_right);
%xdc_focus_times(tx, times_all, delays_all);  % Apply focus times
pressure_field_right = calc_hhp(tx,rx, points);


xdc_apodization(rx, 0, weights_left);
%xdc_focus_times(rx, times_all, delays_all);  % Apply focus times
xdc_apodization(rx, 0, weights_right);
%xdc_focus_times(rx, times_all, delays_all);  % Apply focus times

% Combine fields
pressure_field_combined = pressure_field_left + pressure_field_right;


%% Normalize and Convert to dB Scale
normalized_field = vecnorm(pressure_field_combined, 2, 1);
image_field = reshape(normalized_field, [length(z_range), length(x_range)]);
image_field_lin = image_field / max(image_field(:));

a = (10^(-30/20));
b = 1-a;
hp_image_dB = 20*log10(b*image_field_lin+a);






focus_depth = 40 / 1000; % Focus depth [m]

% Define focal points for scanning lines
focus_points = [x_range; zeros(1, num_scan_lines); repmat(focus_depth, 1, num_scan_lines)]; % [x, y, z]

%% Simulate Scanning Lines
rf_data = zeros(5000, num_scan_lines); % Preallocate RF data matrix (5000 samples per line, adjust as needed)
times = zeros(1, num_scan_lines); % Preallocate times array

for i = 1:num_scan_lines
    if i <= num_scan_lines / 2
        % Use the left half of the transducer for the first 50 lines
        active_elements = left_elements;
        focus_point = left_focus;

    else
        % Use the right half of the transducer for the last 50 lines
        active_elements = right_elements;
        focus_point = right_focus;
    end

  % Define apodization for the current scanning line
    apodization = zeros(1, num_elements);
    apodization(active_elements) = 1; % Activate only the selected half of the transducer
    xdc_apodization(tx, 0, apodization);
    xdc_apodization(rx, 0, apodization);
% Set focus for the current line
    xdc_focus(tx, 0, focus_points(:, i)');
    xdc_focus(rx, 0, focus_points(:, i)');

    % Calculate the received response for the current scanning line
    [v, t1] = calc_scat(tx, rx, phantom_positions, phantom_amplitudes);
   % Store the received RF data and time
    rf_data(1:numel(v), i) = v; % Store the RF data (pad with zeros if v is shorter)
    times(i) = t1;
end



%% Apply Hilbert Transform
% Envelope detection using Hilbert transform
rf_envelope = abs(hilbert(rf_data));

%% Normalize the Envelope Data
% Normalize between 0 and 1
rf_envelope_normalized = rf_envelope / max(rf_envelope(:));


rf_envelope_dB = 20 * log10(rf_envelope_normalized + eps); % Convert to dB scale
rf_envelope_dB(rf_envelope_dB < -30) = -30; % Clip below -30 dB



z_axis = linspace(z_min, z_max, size(rf_data, 1)); % Define the axial (z) axis

%% Display the Result
% Use x_scan for lateral axis and z_axis for axial positions
figure;
imagesc(x_range * 1000, z_axis * 1000, rf_envelope_dB); % Convert axes to mm
colormap('gray'); % Use grayscale for the display
colorbar;
xlabel('Lateral Position [mm]');
ylabel('Axial Position [mm]');
title('B-mode Image of the Phantom');


function delays = Delay(num_elements, c, focus_point, element_positions)
    % INPUT:
    % num_elements - Number of elements in the transducer
    % c - Speed of sound [m/s]
    % focus_point - [x, y, z] coordinates of the focal point [m]
    % element_positions - [x, y, z] positions of the elements [m]

    % Extract Z coordinates of the elements
    delays=zeros(num_elements,1);
    for i = 1:num_elements
        delays(i) = (sqrt(focus_point(1)^2 + focus_point(3)^2)-sqrt((element_positions(i)-focus_point(1))^2 + focus_point(3)^2))/c
    end
end
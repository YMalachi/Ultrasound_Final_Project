%% Initialize Environment

clear all;
clc;
addpath('C:\Users\idans\OneDrive\שולחן העבודה\לימודים\שנה ג\אולטרסאונד\Field_II_ver_3_30_windows');

% Add Field II library
field_init(0);

%% Simulation Parameters
f0 = 3e6;                    % Transducer center frequency [Hz]
fs = 100e6;                  % Sampling frequency [Hz]
c = 1482;                    % Speed of sound [m/s]
width = 0.1 / 1000;          % Width of element [m]
kerf = 0.05 / 1000;           % Kerf [m]
element_height = 1 / 1000;   % Height of element [m]
focus_depth = 40 / 1000;     % Fixed focal depth [m]
N_tx_elements = 128;         % Number of elements in the transmit aperture

% Set simulation parameters
set_sampling(fs);
set_field('c', c);
set_field('use_rectangles', 1);

%% Transducer Setup
tx = xdc_linear_array(N_tx_elements, width, element_height, kerf, 1, 1, [0, 0, focus_depth]);
impulse_response = sin(2 * pi * f0 * (0:1/fs:2/f0)) .* hann(2 * fs / f0)';
xdc_impulse(tx, impulse_response);

rx = xdc_linear_array(N_tx_elements, width, element_height, kerf, 1, 1, [0, 0, focus_depth]);

excitation = sin(2 * pi * f0 * (0:1/fs:1/f0));
xdc_excitation(tx, excitation);

%% Define Focus Points
left_focus = [-8, 0, 40] / 1000;   % Left focus (in meters)
right_focus = [8, 0, 40] / 1000;   % Right focus (in meters)

%% Define Element Positions
element_positions = zeros(N_tx_elements, 1);
for i = 1:N_tx_elements
    element_positions(i) = (i - N_tx_elements / 2) * (kerf + width);
end

%% Define Delays for Left and Right Groups
% Element positions (relative to array center)
element_positions = ((1:N_tx_elements) - (N_tx_elements + 1) / 2) * (width + kerf);

% Left group
left_elements = 1:(N_tx_elements / 2);  % Indices for left group
delays_left = Delay(length(left_elements), c, left_focus, element_positions(left_elements));

% Right group
right_elements = (N_tx_elements / 2 + 1):N_tx_elements;  % Indices for right group
delays_right = Delay(length(right_elements), c, right_focus, element_positions(right_elements));

% Combine delays into full aperture delays
delays_all = zeros(1, N_tx_elements);  % Single row, one column per element
delays_all(left_elements) = delays_left;   % Apply delays for left group
delays_all(right_elements) = delays_right; % Apply delays for right group

% Define focus times (static focus, valid from time 0)
times_all = 0;  % Single scalar for static focus

% Apply focus times and delays
xdc_focus_times(tx, times_all, delays_all);

%% Apply Focus Times and Compute Fields
% Define grid for the acoustic field
x_range = linspace(-15, 15, 100) / 1000;  % Lateral range in meters
z_range = linspace(5, 50, 100) / 1000;    % Axial range in meters
[X, Z] = meshgrid(x_range, z_range);
points = [X(:), zeros(numel(X), 1), Z(:)];

% Left group
weights_left = zeros(1, N_tx_elements);
weights_left(left_elements) = 1;  % Activate left elements only
xdc_apodization(tx, 0, weights_left);
xdc_focus_times(tx, times_all, delays_all);  % Apply focus times

xdc_apodization(rx, 0, weights_left);
xdc_focus_times(rx, times_all, delays_all);  % Apply focus times

pressure_field_left = calc_hhp(tx,rx, points);

% Right group
weights_right = zeros(1, N_tx_elements);
weights_right(right_elements) = 1;  % Activate right elements only
xdc_apodization(tx, 0, weights_right);
xdc_focus_times(tx, times_all, delays_all);  % Apply focus times

xdc_apodization(rx, 0, weights_right);
xdc_focus_times(rx, times_all, delays_all);  % Apply focus times

pressure_field_right = calc_hhp(tx,rx, points);

% Combine fields
pressure_field_combined = pressure_field_left + pressure_field_right;


%% Normalize and Convert to dB Scale
normalized_field = vecnorm(pressure_field_combined, 2, 1);
image_field = reshape(normalized_field, [length(z_range), length(x_range)]);
image_field_lin = image_field / max(image_field(:));

a = (10^(-30/20));
b = 1-a;
hp_image_dB = 20*log10(b*image_field_lin+a);

%% Plot the Combined Acoustic Field
figure;
imagesc(x_range * 1000, z_range * 1000, hp_image_dB);  % Convert to mm
colormap('hot');
colorbar;
xlabel('Lateral Position [mm]');
ylabel('Axial Position [mm]');
title('Acoustic Field (Left and Right Focus with Delays)');

%% Extract Lateral Section at Depth 40 mm
% Find the index for z = 40 mm
z_target = 40 / 1000;  % Convert mm to meters
[~, z_index] = min(abs(z_range - z_target));  % Find the nearest index

% Extract the lateral section
lateral_section = image_field(z_index, :);  % Extract the row at depth 40 mm

% Normalize the lateral section
lateral_section_norm = lateral_section / max(lateral_section);

% Plot the Lateral Section
figure;
plot(x_range * 1000, lateral_section_norm, 'LineWidth', 1.5);  % Convert x_range to mm
xlabel('Lateral Position [mm]');
ylabel('Normalized Amplitude');
title('Lateral Section at Depth of 40 mm');
grid on;

%% %%%%%%%%%%%%%%%%%%%%%%%% new way new me %%%%%%%%%%%%%%%%%%%%%%%%%%

%% Division by Even and Odd Elements
% Define even and odd elements
even_elements = 2:2:N_tx_elements;  % Indices of even elements
odd_elements = 1:2:N_tx_elements;   % Indices of odd elements

% Define focus points for even and odd groups
focus_odd = [-8, 0, 40] / 1000;   % Odd group focus
focus_even = [8, 0, 40] / 1000;   % Even group focus

% Calculate delays for odd and even groups
delays_odd = Delay(length(odd_elements), c, focus_odd, element_positions(odd_elements));
delays_even = Delay(length(even_elements), c, focus_even, element_positions(even_elements));

% Create combined delays and apodization
delays_all = zeros(1, N_tx_elements);
delays_all(odd_elements) = delays_odd;  % Assign delays to odd elements
delays_all(even_elements) = delays_even;  % Assign delays to even elements

apodization = zeros(1, N_tx_elements);
apodization(odd_elements) = 1;  % Activate odd elements
apodization(even_elements) = 1;  % Activate even elements

%% Apply Focus and Compute the Field
% Set focus times and apodization
xdc_apodization(tx, 0, apodization);
xdc_focus_times(tx, 0, delays_all);

xdc_apodization(rx, 0, apodization);
xdc_focus_times(rx, 0, delays_all);

% Define grid for the acoustic field
x_range = linspace(-15, 15, 100) / 1000; % Lateral range in meters
z_range = linspace(5, 50, 100) / 1000;   % Axial range in meters
[X, Z] = meshgrid(x_range, z_range);
points = [X(:), zeros(numel(X), 1), Z(:)];

% Calculate the combined field
pressure_field = calc_hhp(tx,rx, points);

% Normalize and reshape field
normalized_field = vecnorm(pressure_field, 2, 1);
image_field = reshape(normalized_field, [length(z_range), length(x_range)]);

%% Normalize and Convert to dB Scale
image_field_lin = image_field / max(image_field(:));
a = (10^(-30/20));
b = 1-a;
hp_image_dB = 20*log10(b*image_field_lin+a);
%% Plot the Acoustic Field
figure;
imagesc(x_range * 1000, z_range * 1000, hp_image_dB);
colormap('hot');
colorbar;
xlabel('Lateral Position [mm]');
ylabel('Axial Position [mm]');
title('Acoustic Field (Odd and Even Elements with Focus)');

%% Extract and Plot Lateral Section at 40 mm
% Find the index for z = 40 mm
z_target = 40 / 1000;  % Convert mm to meters
[~, z_index] = min(abs(z_range - z_target));  % Find the nearest index

% Extract the lateral section
lateral_section = image_field(z_index, :);  % Extract the row at depth 40 mm

% Normalize the lateral section
lateral_section_norm = lateral_section / max(lateral_section);

% Plot the lateral section
figure;
plot(x_range * 1000, lateral_section_norm, 'LineWidth', 1.5);  % Convert x_range to mm
xlabel('Lateral Position [mm]');
ylabel('Normalized Amplitude');
title('Lateral Section at Depth of 40 mm (Odd and Even Elements)');
grid on;


%%
function delays = Delay(num_elements, c, focus_point, element_positions)
    % INPUT:
    % num_elements - Number of elements in the group
    % c - Speed of sound [m/s]
    % focus_point - [x, y, z] coordinates of the focal point [m]
    % element_positions - Relative positions of the elements [m]

    delays = zeros(num_elements, 1);
    for i = 1:num_elements
        delays(i) = (sqrt(focus_point(1)^2 + focus_point(3)^2) - ...
                     sqrt((element_positions(i) - focus_point(1))^2 + focus_point(3)^2)) / c;
    end
end

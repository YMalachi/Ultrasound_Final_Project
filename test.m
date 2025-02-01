
%% 3.b

% parameters
f0=3.5e6;                   %  Transducer center frequency [Hz]
fs=100e6;                   %  Sampling frequency [Hz]
c=1540;                     %  Speed of sound [m/s]
lambda=c/f0;                %  Wavelength [m]
width=0.1/1000;             %  Width of element
element_height=1/1000;      %  Height of element [m]
kerf=0.05/1000;             %  Kerf [m]
focus=[0 0 60]/1000;        %  Fixed focal point [m]
N_tx_elements=128;          %  Number of physical elements in the transmit aperture
N_rx_elements=128;          %  Number of physical elements in the receive aperture
%  Generate aperture for reception
rx=xdc_linear_array(N_rx_elements,width,element_height,kerf,1,1,focus);

% Set the relevent simulation parameters
set_sampling(fs);                   %  Sets sampling frequency
set_field('use_triangles',0);       %  Tells whether to use triangles (1) or not (0)
set_field('use_rectangles',1);      %  Tells whether to use rectangles (1) or not (0)
set_field('use_att',0);             %  Tells whether to use attenuation (1) or not (0)
set_field('c',c);                   %  Sets the speed of sound

% Generate aperture for transmission
tx=xdc_linear_array(N_tx_elements,width,element_height,kerf,1,1,focus);

% Set the impulse response of the transmit aperture
impulse_response=sin(2*pi*f0*(0:1/fs:2/f0));
impulse_response=impulse_response.*hanning(max(size(impulse_response)))';
xdc_impulse(tx,impulse_response);

% Set the excitation of the transmit aperture
excitation=sin(2*pi*f0*(0:1/fs:2/f0));
xdc_excitation(tx,excitation);
%  Set the impulse response for the receive aperture
xdc_impulse(rx,impulse_response);

% creating phantom 
% Define the parameters
num_points = 100;  % Number of random points
range_x = [-5/1000 5/1000]; % Range of x coordinates
range_z = [30/1000 50/1000]; % Range of z coordinates
x_axis = linspace(range_x(1), range_x(2), 100);
z_axis = linspace(range_z(1), range_z(2), 100);
scatter_points = {[3/1000 0 40/1000; -3/1000 0 40/1000], [10; 10]}; % Scatter points and amplitudes
num_lines = 50; % Number of scanning lines
x_scan = linspace(range_x(1), range_x(2), num_lines); % X coordinates of scanning lines
focus_point = [x_scan; zeros(1, num_lines); repmat(40/1000, 1, num_lines)]; % Focal points for each scanning line

% Generate random points within the specified range
random_x = (range_x(2) - range_x(1)) * rand(num_points, 1) + range_x(1);
random_z = (range_z(2) - range_z(1)) * rand(num_points, 1) + range_z(1);
random_y = zeros(num_points, 1); % Set y coordinate to 0 for all points
random_amplitudes = zeros(num_points, 1); % Set amplitudes to 0 for all points

% Concatenate scatter points with random points
points_x = [random_x; scatter_points{1}(:, 1)];
points_y = [random_y; scatter_points{1}(:, 2)];
points_z = [random_z; scatter_points{1}(:, 3)];
points = [points_x points_y points_z];
amplitudes = [random_amplitudes; scatter_points{2}];

imagesc(range_x*1000, range_z*1000, points);
colormap("gray");
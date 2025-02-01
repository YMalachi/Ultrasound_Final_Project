field_init(0);

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

%% 3.a
X_A_1 = zeros(N_tx_elements,1);
for i = 1:N_tx_elements
    X_A_1(i) = (i-N_tx_elements/2)*(kerf+width);
end
delay = Delay(N_tx_elements,c,[8/1000,0,40/1000],X_A_1);
delay = transpose(repmat(delay,1,N_tx_elements));
time = zeros(size(delay,1),1);
xdc_focus_times(tx,time,delay);

x_A_1 = linspace(-15/1000, 15/1000, 200);
z_A_1 = linspace(5/1000, 70/1000, 200);
[x_A_2, y_A_2, z_A_2] = meshgrid(x_A_1, 0, z_A_1);
points_A = [x_A_2(:) y_A_2(:) z_A_2(:)];
[hp_A, ~] = calc_hp(tx, points_A);
norm_hp_A = zeros(1, 40000);
for i = 1:size(hp_A,2)
    norm_hp_A(i) = norm(hp_A(:,i));
end 
norm_hp_A = reshape(norm_hp_A,[200,200])';

norm_hp_A = (norm_hp_A-min(norm_hp_A,[],'all'))/(max(norm_hp_A,[],'all')-min(norm_hp_A,[],'all'));
dB = 30;
b = 10.^(-dB/20);
a = 1-b;
hp_A_dB = 20*log10(a.*norm_hp_A+b);

figure(1);
imagesc(x_A_1*1000, z_A_1*1000, hp_A_dB); 
colormap('hot'); 
colorbar;
xlabel('Lateral axis X [mm]'); ylabel('Axial axis Z [mm]');
title('The Acoustic Field, XZ axis, 30dB, Focus at (8,0,40)mm')
%export_fig(gcf, 'Q3_a_Field_XZ_new_focus.png', '-png', '-m2')


%% 3.b

%  Generate aperture for reception
rx=xdc_linear_array(N_rx_elements,width,element_height,kerf,1,1,focus);

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

image_data=zeros(1,num_lines);
times=zeros(1,num_lines);
% Loop through each scanning line
for i = 1:num_lines
    % Calculate delays for focusing
    delays_B = Delay(N_tx_elements, c,focus_point(:, i),X_A_1);
    delays_B = transpose(repmat(delays_B,1,N_tx_elements));
    time_B = zeros(size(delays_B,1),1);
    xdc_focus_times(tx,time_B,delays_B);
    xdc_focus_times(rx,time_B,delays_B);

    % Calculate the received response
    [v, t1]=calc_scat(tx, rx, points, amplitudes);
    % Store the result
    image_data(1:max(size(v)),i)=v;
    times(i) = t1;
end

% Apply Hilbert transform to calculate envelope
envelope_data = abs(hilbert(image_data));

% Normalize the envelope data
envelope_data_norm=(envelope_data-min(envelope_data,[],'all'))/(max(envelope_data,[],'all')-min(envelope_data,[],'all'));

% Convert to logarithmic scale of 30dB
envelope_data_log = 20*log10(a.*envelope_data_norm+b);


% Display the image
figure(2);
imagesc(envelope_data_log, 'XData', x_axis*1000, 'YData', z_axis*1000); colormap(gray);
xlabel('X [mm]'); ylabel('Z [mm]'); title('Phantom image - 50 scanning lines');colorbar;
%export_fig(gcf, 'Q3_b_Field_XZ_phantom.png', '-png', '-m2')

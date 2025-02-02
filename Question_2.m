clear; clc; 
%% flags
plot_flag = 1; 

field_init(0);

%% params
% initiating the transducer and params
linear_128_trans = xdc_linear_array(128,0.1/1000,1/1000,0.05/1000,1,1,[0 0 80/1000]);
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
xdc_excitation(linear_128_trans, excitation);
xdc_impulse(linear_128_trans, impulse);

%% 2a

focus = [0 0 0];
focus_lst = [80/1000, 40/1000];
for i = 1:2
    focus(1,3) = focus_lst(i);

    xdc_focus(linear_128_trans,0,focus);
    % initiating vectors for field matrix
    x_field_vector = linspace(-15/1000, 15/1000, 100);
    z_field_vector = linspace(10/1000,100/1000, 100);
    [x, y, z] = meshgrid(x_field_vector, 0, z_field_vector);
    points = [x(:) y(:) z(:)];
    
    % calculating the pressure field
    pressure_field_2a = calc_hp(linear_128_trans, points);
    log_scale_pressure_field_2a = log_scale_field(pressure_field_2a);

    if plot_flag
        figure;
        imagesc(x_field_vector*1000, z_field_vector*1000, log_scale_pressure_field_2a);
        colormap('hot');
        colorbar;
        title('Pressure Field Image, 128 Element Transducer, Log Scaled 30dB', FontSize=15);
        xlabel('X [mm]');
        ylabel('Z [mm]');
    end
    
    %% 2b
    
    z_lateral = 40/1000; % mm
    z_lateral_idx = find(z_field_vector == z_lateral); % row number
    
    % plot
    z_lateral_vector = log_scale_pressure_field_2a(z_lateral_idx,:);
    if plot_flag
        figure;
        plot(x_field_vector,z_lateral_vector+30); % +30 because it was negative
        title('Lateral Cut of Pressure Field at 40mm Depth, 128 Elements Linear Transducer', FontSize=15);
        xlabel('X [mm]');
        ylabel('Amplitude [dB]');
    end
end

%% 2c
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% it was implemented with the for loop at the beginning

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2d

clear;
plot_flag = 1;
% initiating the transducer and params
linear_128_trans = xdc_linear_array(128,0.1/1000,1/1000,0.05/1000,1,1,[0 0 40/1000]);
f0 = 3e6; % main frequency
t0 = 1/f0; % period time
fs = 100e6; % sampling rate
ts = 1/fs; % sampling time
set_sampling(fs);
c = 1540; % m/sec, speed of sound in soft tissue

% initiating apodizations
rect_apod = ones(1,128);
hamming_apod = transpose(hamming(128,"symmetric"));
for i = 1:2
    if i == 1
        xdc_apodization(linear_128_trans, 0, rect_apod);
        time_vec = 0:ts:t0; % for one cycle
        excitation = sin(2*pi*time_vec*f0);
        impulse = sin(2*pi*time_vec*f0);
        xdc_excitation(linear_128_trans, excitation);
        xdc_impulse(linear_128_trans, impulse);
        
        % initiating vectors for field matrix
        x_field_vector = linspace(-15/1000, 15/1000, 100);
        z_field_vector = linspace(10/1000,100/1000, 100);
        [x, y, z] = meshgrid(x_field_vector, 0, z_field_vector);
        points = [x(:) y(:) z(:)];
        
        % calculating the pressure field
        pressure_field_2d_rect = calc_hp(linear_128_trans, points);
        log_scale_pressure_field_2d_rect = log_scale_field(pressure_field_2d_rect);
        if plot_flag
            figure(5);
            subplot(1,2,1);
            imagesc(x_field_vector*1000, z_field_vector*1000, log_scale_pressure_field_2d_rect);
            colormap('hot');
            colorbar;
            title('Pressure Field Image, Rect Apodization', FontSize=15);
            xlabel('X [mm]');
            ylabel('Z [mm]');
            z_lateral = 40/1000; % mm
            z_lateral_idx = find(z_field_vector == z_lateral); % row number
        end
        % plot
        if plot_flag
            z_lateral_vector = log_scale_pressure_field_2d_rect(z_lateral_idx,:);
            figure(6);
            subplot(1,2,1);
            plot(x_field_vector,z_lateral_vector+30); % +30 because it was negative
            title('Lateral Cut of Pressure Field at 40mm Depth, Rect Apodization', FontSize=15);
            xlabel('X [mm]');
            ylabel('Amplitude [dB]');
        end
    else
        % hamming apodization
        xdc_apodization(linear_128_trans, 0, hamming_apod);
        time_vec = 0:ts:t0; % for one cycle
        excitation = sin(2*pi*time_vec*f0);
        impulse = sin(2*pi*time_vec*f0);
        xdc_excitation(linear_128_trans, excitation);
        xdc_impulse(linear_128_trans, impulse);
        
        % initiating vectors for field matrix
        x_field_vector = linspace(-15/1000, 15/1000, 100);
        z_field_vector = linspace(10/1000,100/1000, 100);
        [x, y, z] = meshgrid(x_field_vector, 0, z_field_vector);
        points = [x(:) y(:) z(:)];
        
        % calculating the pressure field
        pressure_field_2d_hamm = calc_hp(linear_128_trans, points);
        log_scale_pressure_field_2d_hamm = log_scale_field(pressure_field_2d_hamm);
        if plot_flag
            figure(5);
            subplot(1,2,2);
            imagesc(x_field_vector*1000, z_field_vector*1000, log_scale_pressure_field_2d_hamm);
            colormap('hot');
            colorbar;
            title('Pressure Field Image, Hamming Apodization', FontSize=15);
            xlabel('X [mm]');
            ylabel('Z [mm]');
            z_lateral = 40/1000; % mm
            z_lateral_idx = find(z_field_vector == z_lateral); % row number
        end
        % plot
        
        if plot_flag
            z_lateral_vector = log_scale_pressure_field_2d_hamm(z_lateral_idx,:);
            figure(6);
            subplot(1,2,2);
            plot(x_field_vector,z_lateral_vector+30); % +30 because it was negative
            title('Lateral Cut of Pressure Field at 40mm Depth, Hamming Apodization', FontSize=15);
            xlabel('X [mm]');
            ylabel('Amplitude [dB]');
        end
    end
end

%% 2e

clear;
plot_flag = 1;

% test
% initiating the transducer and params
linear_128_trans = xdc_linear_array(128,0.1/1000,1/1000,0.05/1000,1,1,[0 0 80/1000]);
linear_128_trans_return = xdc_linear_array(128,0.1/1000,1/1000,0.05/1000,1,1,[0 0 80/1000]);
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
xdc_excitation(linear_128_trans, excitation);
xdc_impulse(linear_128_trans, impulse);
xdc_excitation(linear_128_trans_return, excitation);
xdc_impulse(linear_128_trans_return, impulse);

% values for element distance vector calculation
ele_num = 128;
pitch = 0.15/1000;
ele_pos_vec = -64:1:64;
ele_pos_vec(65) = [];
ele_dist(1:64) = ele_pos_vec(1:64)*pitch + 0.5*pitch;
ele_dist(65:128) = ele_pos_vec(65:128)*pitch - 0.5*pitch;
% calculating using the Delay function
time_vec_delay = Delay(ele_num, ele_dist, c, [10/1000 0 40/1000]); % see Delay function in the Delay.m function file
xdc_focus_times(linear_128_trans, 0, time_vec_delay);
xdc_focus_times(linear_128_trans_return, 0, time_vec_delay);

% initiating vectors for field matrix
x_field_vector = linspace(-15/1000, 15/1000, 100);
z_field_vector = linspace(10/1000,100/1000, 100);
[x, y, z] = meshgrid(x_field_vector, 0, z_field_vector);
points = [x(:) y(:) z(:)];

% calculating the pressure field
pressure_field_2e = calc_hhp(linear_128_trans, linear_128_trans, points);
log_scale_pressure_field_2e = log_scale_field(pressure_field_2e);

if plot_flag
    figure;
    imagesc(x_field_vector*1000, z_field_vector*1000, log_scale_pressure_field_2e);
    colormap('hot');
    colorbar;
    title('Pressure Field Image, 128 Element Transducer, Log Scaled 30dB', FontSize=15);
    xlabel('X [mm]');
    ylabel('Z [mm]');
end


z_lateral = 40/1000; % mm
z_lateral_idx = find(z_field_vector == z_lateral); % row number

% plot
z_lateral_vector = log_scale_pressure_field_2e(z_lateral_idx,:) + 30; % so it wont be negative

if plot_flag
    figure;
    plot(x_field_vector, z_lateral_vector); % +30 because it was negative
    title('Lateral Cut of Pressure Field at 40mm Depth, 128 Elements Linear Transducer', FontSize=15);
    xlabel('X [mm]');
    ylabel('Amplitude [dB]');
end
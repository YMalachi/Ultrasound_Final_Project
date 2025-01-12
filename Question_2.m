clear;
%flags
plot_flag = 0; 

%Th = xdc linear array (no elements, width, height, kerf, no sub x, no sub y, focus
field_init(0);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
focus = [0 0 0];
focus_lst = [40/1000, 80/1000];
for i = 1:2
    focus(1,3) = focus_lst(i);

    xdc_focus(linear_128_trans,0,focus);
    % initiating vectors for field matrix
    x_field_vector = linspace(-15/1000, 15/1000, 100);
    z_field_vector = linspace(10/1000,100/1000, 100);
    [x, y, z] = meshgrid(x_field_vector, 0, z_field_vector);
    points = [x(:) y(:) z(:)];
    
    % calculating the pressure field
    [hp, start_time] = calc_hp(linear_128_trans, points);
    
    % normalizing
    normalized_hp = [];
    for idx = 1:length(hp(1,:))
        normalized_hp(idx) = norm(hp(:,idx));
    end
    scaled_hp = reshape(normalized_hp,[100 100]); % scaling
    min_scaled_hp = min(scaled_hp,[], "all");
    max_scaled_hp = max(scaled_hp,[], "all");
    delta_max_min = max_scaled_hp - min_scaled_hp;
    linear_norm_scaled_hp = (scaled_hp - min_scaled_hp)/(delta_max_min);
    
    % log scale
    a = (10^-1.5);
    b = 1-a;
    log_scale_hp = transpose(20*log10(b*linear_norm_scaled_hp+a));
    if plot_flag
        figure;
        imagesc(x_field_vector*1000, z_field_vector*1000, log_scale_hp);
        colormap('hot');
        title('Pressure Field Image, 128 Element Transducer, Log Scaled 30dB', FontSize=15);
        xlabel('X [mm]');
        ylabel('Z [mm]');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 2b
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    z_lateral = 40/1000; % mm
    z_lateral_idx = find(z_field_vector == z_lateral); % row number
    
    % plot
    z_lateral_vector = log_scale_hp(z_lateral_idx,:);
    if plot_flag
        figure;
        plot(x_field_vector,z_lateral_vector+30); % +30 because it was negative
        title('Lateral Cut of Pressure Field at 40mm Depth, 128 Elements Linear Transducer', FontSize=15);
        xlabel('X [mm]');
        ylabel('Amplitude [dB]');
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2c
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% it was implemented with the for loop at the beginning

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
plot_flag = 0;
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
        [hp, start_time] = calc_hp(linear_128_trans, points);
        % normalizing
        normalized_hp = [];
        for idx = 1:length(hp(1,:))
            normalized_hp(idx) = norm(hp(:,idx));
        end
        scaled_hp = reshape(normalized_hp,[100 100]); % scaling
        min_scaled_hp = min(scaled_hp,[], "all");
        max_scaled_hp = max(scaled_hp,[], "all");
        delta_max_min = max_scaled_hp - min_scaled_hp;
        linear_norm_scaled_hp = (scaled_hp - min_scaled_hp)/(delta_max_min);
        
        % log scale
        a = (10^-1.5);
        b = 1-a;
        log_scale_hp = transpose(20*log10(b*linear_norm_scaled_hp+a));
        %if plot_flag
        figure;
        imagesc(x_field_vector*1000, z_field_vector*1000, log_scale_hp);
        colormap('hot');
        title('Pressure Field Image, Rect Apodization', FontSize=15);
        xlabel('X [mm]');
        ylabel('Z [mm]');
        z_lateral = 40/1000; % mm
        z_lateral_idx = find(z_field_vector == z_lateral); % row number

        % plot
        z_lateral_vector = log_scale_hp(z_lateral_idx,:);
        if plot_flag
            figure;
            plot(x_field_vector,z_lateral_vector+30); % +30 because it was negative
            title('Lateral Cut of Pressure Field at 40mm Depth, Rect Apodization', FontSize=15);
            xlabel('X [mm]');
            ylabel('Amplitude [dB]');
        end
    else
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
        [hp, start_time] = calc_hp(linear_128_trans, points);
        % normalizing
        normalized_hp = [];
        for idx = 1:length(hp(1,:))
            normalized_hp(idx) = norm(hp(:,idx));
        end
        scaled_hp = reshape(normalized_hp,[100 100]); % scaling
        min_scaled_hp = min(scaled_hp,[], "all");
        max_scaled_hp = max(scaled_hp,[], "all");
        delta_max_min = max_scaled_hp - min_scaled_hp;
        linear_norm_scaled_hp = (scaled_hp - min_scaled_hp)/(delta_max_min);
        
        % log scale
        a = (10^-1.5);
        b = 1-a;
        log_scale_hp = transpose(20*log10(b*linear_norm_scaled_hp+a));
        if plot_flag
            figure;
            imagesc(x_field_vector*1000, z_field_vector*1000, log_scale_hp);
            colormap('hot');
            title('Pressure Field Image, Hamming Apodization', FontSize=15);
            xlabel('X [mm]');
            ylabel('Z [mm]');
            z_lateral = 40/1000; % mm
            z_lateral_idx = find(z_field_vector == z_lateral); % row number
        end
        % plot
        z_lateral_vector = log_scale_hp(z_lateral_idx,:);
        if plot_flag
            figure;
            plot(x_field_vector,z_lateral_vector+30); % +30 because it was negative
            title('Lateral Cut of Pressure Field at 40mm Depth, Hamming Apodization', FontSize=15);
            xlabel('X [mm]');
            ylabel('Amplitude [dB]');
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2e
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
plot_flag = 1;

function [time_delay_vec] = Delay(element_number, element_distances, c, focus_point)
    % focus point is a vector that looks like that [x,0,z]
    time_delay_vec = [];
    x_co = focus_point(1);
    y_co = focus_point(2);
    z_co = focus_point(3);
    dist_to_focus_from_center = sqrt(x_co^2 + z_co^2);
    for curr_ele = 1:element_number % assuming there are always an even number of elements, for symmetry
        dist_to_focus_from_curr_ele = sqrt((element_distances(curr_ele)-x_co)^2+z_co^2);
        delay_ele = (dist_to_focus_from_center - dist_to_focus_from_curr_ele)/c; % c is [m/sec]
        time_delay_vec(curr_ele) = delay_ele;
    end
end

% test
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

ele_num = 128;
pitch = 0.15/1000;
ele_pos_vec = -64:1:64;
ele_pos_vec(65) = [];
ele_dist = ele_pos_vec * pitch;

time_vec_delay = Delay(ele_num, ele_dist, 1540, [10 0 40]);
xdc_focus_times(linear_128_trans,0,time_vec_delay);

% initiating vectors for field matrix
x_field_vector = linspace(-15/1000, 15/1000, 100);
z_field_vector = linspace(10/1000,100/1000, 100);
[x, y, z] = meshgrid(x_field_vector, 0, z_field_vector);
points = [x(:) y(:) z(:)];

% calculating the pressure field
[hp, start_time] = calc_hp(linear_128_trans, points);

% normalizing
normalized_hp = [];
for idx = 1:length(hp(1,:))
    normalized_hp(idx) = norm(hp(:,idx));
end
scaled_hp = reshape(normalized_hp,[100 100]); % scaling
min_scaled_hp = min(scaled_hp,[], "all");
max_scaled_hp = max(scaled_hp,[], "all");
delta_max_min = max_scaled_hp - min_scaled_hp;
linear_norm_scaled_hp = (scaled_hp - min_scaled_hp)/(delta_max_min);

% log scale
a = (10^-1.5);
b = 1-a;
log_scale_hp = transpose(20*log10(b*linear_norm_scaled_hp+a));
if plot_flag
    figure;
    imagesc(x_field_vector*1000, z_field_vector*1000, log_scale_hp);
    colormap('hot');
    title('Pressure Field Image, 128 Element Transducer, Log Scaled 30dB', FontSize=15);
    xlabel('X [mm]');
    ylabel('Z [mm]');
end


z_lateral = 40/1000; % mm
z_lateral_idx = find(z_field_vector == z_lateral); % row number

% plot
z_lateral_vector = log_scale_hp(z_lateral_idx,:);
if plot_flag
    figure;
    plot(x_field_vector,z_lateral_vector+30); % +30 because it was negative
    title('Lateral Cut of Pressure Field at 40mm Depth, 128 Elements Linear Transducer', FontSize=15);
    xlabel('X [mm]');
    ylabel('Amplitude [dB]');
end
clear;
%addpath('C:\Users\yotam\Downloads\Field_II_Toolbox') % we don't need this
%because we added it permanently but add yours so it works

% functions we've written and are used in this script:
function h = PlotSignal(data, fs, titleText, xLabel, yLabel, TimeOrFreq)
    % immportant: you need to initiate your own figures!!!!!
    % data is the signal data (1D array)
    % fs is the sampling frequency, number
    % titleText is string
    % xLabel is string
    % yLabel is string
    % TimeOrFreq = 'time' -> time domain analysis
    % TimeOrFreq = 'frequency' -> frequency domain analysis

    if TimeOrFreq == "time"
        L = length(data);
        time_vec = (0:L-1)/fs;
        h = plot(time_vec, data);
        title(titleText, FontSize=15);
        xlabel(xLabel);
        ylabel(yLabel);
        axis tight;
        linkaxes(findall(gcf, 'Type', 'axes'), 'xy');
    end

    if TimeOrFreq == "frequency"
        L = length(data);
        freq_vec = (0:L-1) * (fs/2*(L-1)); % only positive frequencies
        h = plot(freq_vec, data);
        title(titleText, FontSize=15);
        xlabel(xLabel);
        ylabel(yLabel);
        axis tight;
        linkaxes(findall(gcf, 'Type', 'axes'), 'xy');
    end
    
end

% initiate the toolbox
field_init(0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define transducer params
radius = 10/1000; % meters
focal_radius = 40/1000; % meters
element_size = 1/1000; % meters

% initiate transducer
round_trans = xdc_concave(radius, focal_radius, element_size);
show_xdc_geir(round_trans, 1);
axis equal; view(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define further params for simulation
f0 = 3e6; % main frequency
t0 = 1/f0; % period time
fs = 50e6; % sampling rate
ts = 1/fs; % sampling time
set_sampling(fs); 
c = 1482; % m/sec, speed of sound in water
t_ten_cycles = 0:1/fs:10*t0; % time vector for ten cycles
t_two_cycles = 0:1/fs:2*t0; % time vector for two cycles
% define points
point1 = [0 0 20]/1000; % meters
point2 = [0 0 35]/1000; % meters
% initiating excitation and impulse response
excitation = sin(2*pi*f0*t_ten_cycles); % transducer excitation 
impulse = sin(2*pi*f0*t_two_cycles); % transducer impulse response
xdc_excitation(round_trans, excitation); % setting this excitation to the transducer
xdc_impulse(round_trans, impulse); % setting this impulse response to the transducer

% calculate the pressure field at point 1 and point 2
[pressure_field_point1, start_time1] = calc_hp(round_trans, point1);
[pressure_field_point2, start_time2] = calc_hp(round_trans, point2);
% prepare for plotting in time domain
titles = struct('pressure_field_point1', 'Pressure Field 2cm Deep, Time Domain', ...
    'pressure_field_point2', 'Pressure Field 3.5cm Deep, Time Domain');
fields = {'pressure_field_point1', 'pressure_field_point2'};

% plots time domain
figure;
for i = 1:length(fields)
    field_name = fields{i}; 
    field_data = eval(field_name); 
    header = titles.(field_name);
    x_label = 'Time [sec]';
    y_label = 'Amplitude [MPa]';
    subplot(1,length(fields),i);
    PlotSignal(field_data, fs, header, x_label, y_label, "time");
end

% prepare for plotting in frequency domain
L1 = length(pressure_field_point1);
L2 = length(pressure_field_point2);
magnitude1 = abs(fft(pressure_field_point1))/L1; % looking solely at magnitude, normalized
magnitude2 = abs(fft(pressure_field_point2))/L2;
fft_result1 = magnitude1(1:(L1/2)+1); % keeping only positive frequencies
fft_result2 = magnitude2(1:(L2/2)+1);
titles = struct('fft_result1', 'Pressure Field 2cm Deep, Frequency Domain', ...
    'fft_result2', 'Pressure Field 3.5cm Deep, Frequency Domain');
fields_freq = {'fft_result1', 'fft_result2'};

% plots, frequency domain
figure;
for i = 1:length(fields_freq)
    field_name = fields_freq{i}; 
    field_data = eval(field_name); 
    header = titles.(field_name);
    x_label = 'Frequency [Hz]';
    y_label = 'Amplitude [MPa]';
    subplot(1,length(fields_freq),i);
    PlotSignal(field_data, fs, header, x_label, y_label, "frequency");
    xlim([0 100000000000]); % specific for our data, this makes it look better. fs is not high enough tho, looks kinda spiky.
end

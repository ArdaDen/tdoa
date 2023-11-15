close all;
clear all;
clc;

target_pos = [7500,1500];
fs = 8e6;
snr = 10;
min_range = 100;
max_range = 1e6;
num_of_range_bins = 100;
min_degree = 5;
num_of_elements = 4;
min_num_of_elements = 3;
max_num_of_elements = 6;
num_of_degree_bins = 40;
max_degree = 360/num_of_elements-360/num_of_elements/num_of_degree_bins;
min_snr = -10;
max_snr = 20;
snr_int = 5;
offset_range = 0;
offset_angle = 0;
shape = "Circle";
frame_time_length = 0.1;
signal_time_length = 1e-3;
starting_time = 1e-3;
theta = pi/4;
radius = 50000;

rms_vec = zeros(1,length(min_num_of_elements:1:max_num_of_elements));
count1 = 1;
for num_of_elements = min_num_of_elements:1:max_num_of_elements
        y = sensor_pos(radius,theta,offset_range,offset_angle*pi/180,shape,num_of_elements);
        frame_signal = randn(1,round(frame_time_length*fs));
        signal = frame_signal(ceil(starting_time*fs):ceil(starting_time*fs + signal_time_length*fs));
        match_filter = signal(end:-1:1).';
        match_filter = match_filter';
        filtered_signals = zeros(length(y(1,:)),length(signal));
        for i = 1:length(y(:,1))
            t_i = norm(y(i,:)-target_pos)/physconst("LightSpeed");
            sample_delay = ceil(t_i*fs);
            signal_i = frame_signal((ceil(starting_time*fs) + sample_delay):(ceil(starting_time*fs + signal_time_length*fs) + sample_delay));
            signal_i_w_noise = awgn(signal_i,snr,"measured");
            signal_i_w_noise_filtered = ifft(fft(match_filter).*fft(signal_i_w_noise));
            filtered_signals(i,:) = signal_i_w_noise_filtered(end:-1:1);
        end
        sig = filtered_signals;

        light_speed = physconst("LightSpeed");
        time_values = zeros(1,length(sig(:,1)));
        for i = 1:length(sig(:,1))
            [pk,idx] = max(sig(i,:));
            t_m = (idx-1)/fs;
            time_values(i) = t_m;
        end

        A_matrix = zeros(length(y(:,1))-1,3);
        b_matrix = zeros(length(y(:,1))-1,1);

        for i = 2:length(y(:,1))
            A_matrix((i-1),:) = [y(1,:)-y(i,:),-light_speed*(time_values(i)-time_values(1))];
            b_matrix((i-1),:) = light_speed^2*(time_values(i)-time_values(1))^2 + y(1,1)^2 + y(1,2)^2 - y(i,1)^2 - y(i,2)^2;
        end

        A_matrix = 2*A_matrix;
        pseudo_inv = (((A_matrix.')*A_matrix))\(A_matrix.');
        sol = pseudo_inv*b_matrix;
        x_target = sol(1);
        y_target = sol(2);
        result = [x_target,y_target];
        rms = sqrt(((x_target-target_pos(1))^2 + (y_target-target_pos(2))^2)/2);
        rms_vec(count1) = rms;
        count1 = count1 + 1;
end

num_of_elements_vec = min_num_of_elements:1:max_num_of_elements;
figure;
plot(num_of_elements_vec,rms_vec)
xlabel("Number of sensors")
ylabel("RMS error")
title("RMS values according to number of elements")
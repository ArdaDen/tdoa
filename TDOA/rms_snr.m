
clear all;
clc;

target_pos = [7500,1500];
fs = 8e6;
%snr = 10;
min_range = 100;
max_range = 1e6;
num_of_range_bins = 100;
min_degree = 5;
num_of_elements = 4;
num_of_degree_bins = 40;
max_degree = 360/num_of_elements-360/num_of_elements/num_of_degree_bins;
min_snr = -50;
max_snr = 20;
snr_int = 5;
offset_range = 0;
offset_angle = 0;
shape = "Circle";
frame_time_length = 0.1;
signal_time_length = 1e-3;
starting_time = 1e-3;
theta = pi/6;
radius = 50000;
num_of_monte_carlo_runs = 1000;

target_range_error = zeros(num_of_monte_carlo_runs,length(min_snr:snr_int:max_snr));
count1 = 1;
count2 = 1;
for snr = min_snr:snr_int:max_snr
    for x = 1:num_of_monte_carlo_runs
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
        target_range_error(count1,count2) = sqrt((x_target-target_pos(1))^2 + (y_target-target_pos(2))^2);
        count1 = count1 + 1;
    end
    count1 = 1;
    count2 = count2 + 1;
end
rms_vec = zeros(1,length(min_snr:snr_int:max_snr));
for m = 1:length(rms_vec(1,:))
    rms_vec(1,m) = sqrt(sum((target_range_error(:,m).^2))/num_of_monte_carlo_runs);
end
snr_vec = min_snr:snr_int:max_snr;
figure;
semilogy(snr_vec,rms_vec)
xlabel("SNR")
ylabel("RMS error")
title("RMS values according to SNR values")
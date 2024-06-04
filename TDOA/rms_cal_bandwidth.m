close all;
clear all;
clc;

%% Input section
fs = 1e5; % Sampling rate
light_speed = physconst("LightSpeed"); % Speed of light
theta = 2*pi/3;
r = 20000;
emitter = [0,0]; % Emmitter location
number_of_sensors = 5; % Number of sensors
ref = [1000,0]; % Signal model collection location
targets = [0,0]; % Target location
resampling_rate_max = 10; % Maximum resampling rate
rms_errors = zeros(1,resampling_rate_max);

%% Constructing noise
noise = randn(1,15000);

%% Locations
receiver_ref = [-r/3,0]; % Reference receiver location
receiver1 = [r*cos(theta),r*sin(theta)]; % Receiver 1 location
receiver2 = [r,0]; % Receiver 2 location
receiver3 = [r*cos(theta),-r*sin(theta)]; % Receiver 3 location

r_all = [];
zero_cross_l = [];
bw_l = [];

for resampling_rate = 1:resampling_rate_max

    noise_resampled = resample(noise,resampling_rate,1); % Resampled noise
    t = 0:1/fs:(0.1-1/fs); % Time vector
    r1_t = [];
    for target_data = 1:1
        tic
        targets = [targets(1)-r/5*target_data,targets(2)-r/5*target_data];
        %% Delay and range calculations
        t_ref = (norm(emitter-ref))/light_speed; % Delay in reference receiver
        t_r = (norm(targets-receiver_ref) + norm(emitter-targets))/light_speed;
        t_1 = (norm(targets-receiver1) + norm(emitter-targets))/light_speed;
        t_2 = (norm(targets-receiver2) + norm(emitter-targets))/light_speed;
        t_3 = (norm(targets-receiver3) + norm(emitter-targets))/light_speed;

        ref_delay_sample = ceil(t_ref*fs);
        r1_delay_sample = ceil(t_1*fs);
        r2_delay_sample = ceil(t_2*fs);
        r3_delay_sample = ceil(t_3*fs);
        rref_delay_sample = ceil(t_r*fs);

        ref_signal = noise_resampled(ref_delay_sample:ref_delay_sample+10000-1);
        r1_signal = noise_resampled(r1_delay_sample+3000:r1_delay_sample+13000-1);
        r2_signal = noise_resampled(r2_delay_sample+3000:r2_delay_sample+13000-1);
        r3_signal = noise_resampled(r3_delay_sample+3000:r3_delay_sample+13000-1);
        rref_signal = noise_resampled(rref_delay_sample+3000:rref_delay_sample+13000-1);

        ref_signal_w_noise = awgn(ref_signal,30,'measured');
        r1_signal_w_noise = awgn(r1_signal,10,'measured');
        r2_signal_w_noise = awgn(r2_signal,10,'measured');
        r3_signal_w_noise = awgn(r3_signal,10,'measured');
        rref_signal_w_noise = awgn(rref_signal,10,'measured');

        match_filter = ref_signal_w_noise(end:-1:1).';
        match_filter = match_filter';

        r1_signal_w_noise_filtered = ifft(fft(r1_signal_w_noise).*fft(match_filter));
        r2_signal_w_noise_filtered = ifft(fft(r2_signal_w_noise).*fft(match_filter));
        r3_signal_w_noise_filtered = ifft(fft(r3_signal_w_noise).*fft(match_filter));
        rref_signal_w_noise_filtered = ifft(fft(rref_signal_w_noise).*fft(match_filter));

        r1_signal_filtered = conv(match_filter,r1_signal);
        r1_signal_filtered = r1_signal_filtered(round(length(r1_signal_filtered)/2):end);
        fft_r1 = mag2db(abs(fftshift(fft(r1_signal)))/length(t));
        fft_r1_noise = mag2db(abs(fftshift(fft(r1_signal_w_noise)))/length(t));
        freq = linspace(-fs/2,fs/2-fs/length(t),length(t));
%         figure;
%         plot(freq,fft_r1);
        [xcor,k] = xcorr(r1_signal,ref_signal);
        xcor = abs(xcor);
      %  figure;
        plot(k/fs,xcor);
        title("Cross correlations")
        hold all;
        [pk,id] = max(xcor);
        id_up = id;
        id_down = id;
        for z = 1:200
            if xcor(id_up) > xcor(id_up+1)
                id_up = id_up + 1;
            else 
                break
            end
        end
        for z = 1:200
            if xcor(id_down) > xcor(id_down-1)
                id_down = id_down - 1;
            else 
                break
            end
        end
        zero_c = id_up - id_down;
        zero_cross_l = [zero_cross_l;zero_c];
        [pk_fft,id_fft] = max(fft_r1);
        db3 = pk_fft-3;
        [point1,idd1] = min(abs(fft_r1-db3));
        if idd1 > id_fft
            fft_r1(id_fft:end) = [];
            [point2,idd2] = min(abs(fft_r1-db3));
            bw = abs(idd1 - idd2);
        else
            fft_r1(1:id_fft) = [];
            [point2,idd2] = min(abs(fft_r1-db3));
            bw = abs(idd1 - (idd2 + id_fft));
        end

        bw_l = [bw_l;bw];

        [pk1,id1] = max(abs(r1_signal_w_noise_filtered));
        [pk2,id2] = max(abs(r2_signal_w_noise_filtered));
        [pk3,id3] = max(abs(r3_signal_w_noise_filtered));
        [pkref,idref] = max(abs(rref_signal_w_noise_filtered));
        t_values = [id1,id2,id3,idref]/fs;

        r1_t_values_tdoa_ref = (t_values(1));
        r2_t_values_tdoa_ref = (t_values(2));
        r3_t_values_tdoa_ref = (t_values(3));
        rref_t_values_tdoa_ref = (t_values(4));

        r1_t_values_delay = r1_t_values_tdoa_ref + t_ref;
        r2_t_values_delay = r2_t_values_tdoa_ref + t_ref;
        r3_t_values_delay = r3_t_values_tdoa_ref + t_ref;
        rref_t_values_delay = rref_t_values_tdoa_ref + t_ref;

        r1_t = [r1_t;[r1_t_values_delay,t_1 + 3000/fs]];

        A_matrix = 2*[receiver_ref-receiver1,-light_speed*(r1_t_values_delay-rref_t_values_delay);...
            receiver_ref-receiver2,-light_speed*(r2_t_values_delay-rref_t_values_delay);...
            receiver_ref-receiver3,-light_speed*(r3_t_values_delay-rref_t_values_delay)];
        b_matrix = [light_speed^2*(r1_t_values_delay-rref_t_values_delay)^2 + receiver_ref(1)^2 + receiver_ref(2)^2 - receiver1(1)^2 - receiver1(2)^2;...
            light_speed^2*(r2_t_values_delay-rref_t_values_delay)^2 + receiver_ref(1)^2 + receiver_ref(2)^2 - receiver2(1)^2 - receiver2(2)^2;...
            light_speed^2*(r3_t_values_delay-rref_t_values_delay)^2 + receiver_ref(1)^2 + receiver_ref(2)^2 - receiver3(1)^2 - receiver3(2)^2];
        sol = pinv(A_matrix)*b_matrix;
        x = sol(1);
        y = sol(2);

        rms = sqrt(((x-targets(1))^2 + (y-targets(2))^2)/2);
      %  rms_errors(resampling_rate) = rms;
        %   close all;

    end
    r_all = [r_all;r1_t];
end

% re = 1:resampling_rate_max;
% figure;
% plot(re,rms_errors);
zero_cross_l = zero_cross_l';
bw_l = bw_l';
figure;
plot(bw_l,zero_cross_l);
set(gca,"XDir","reverse");
% for l = 1:10
%     plot(r_all(((11*l-10):(12*l-10)),1),r_all(((11*l-10):(12*l-10)),2));
%     hold on;
% end
% hold off;

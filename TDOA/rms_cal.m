close all;
clear all;
clc;

%% Input section
fs = 1e5; % Sampling rate
t = 0:1/fs:(0.1-1/fs); % Time vector
light_speed = physconst("LightSpeed"); % Speed of light
theta_div = 12;
r_start = 100;
r_div = 100;
r_finish = 20000;
rms_errors = zeros(length(r_start:r_div:r_finish),length(0:(pi/theta_div):pi));
emitter = [0,0]; % Emmitter location
number_of_sensors = 5; % Number of sensors
ref = [1000,0]; % Signal model collection location
targets = [0,0]; % Target location
i = 1;

%% Constructing the received signals
filename = "C:\Users\user\Desktop\staj\Recorded sounds\SDRuno_20230815_110704Z_102389kHz.wav"; % FM signal
[y,fs] = audioread(filename); % Reading the signal
y = y(:,1) + 1i*y(:,2); % IQ data
y_selec = y(50000:55000)';

for r = r_start:r_div:r_finish
    p = 1;
    for theta = 0:pi/theta_div:pi
        tic
        receiver_ref = [-r/3,0]; % Reference receiver location
        receiver1 = [r*cos(theta),r*sin(theta)]; % Receiver 1 location
        receiver2 = [r,0]; % Receiver 2 location
        receiver3 = [r*cos(theta),-r*sin(theta)]; % Receiver 3 location

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

        ref_signal = [zeros(1,ref_delay_sample),y_selec,zeros(1,length(t)-ref_delay_sample-length(y_selec))];
        r1_signal = [zeros(1,r1_delay_sample),y_selec,zeros(1,length(t)-r1_delay_sample-length(y_selec))];
        r2_signal = [zeros(1,r2_delay_sample),y_selec,zeros(1,length(t)-r2_delay_sample-length(y_selec))];
        r3_signal = [zeros(1,r3_delay_sample),y_selec,zeros(1,length(t)-r3_delay_sample-length(y_selec))];
        rref_signal = [zeros(1,rref_delay_sample),y_selec,zeros(1,length(t)-rref_delay_sample-length(y_selec))];

        ref_signal_w_noise = awgn(ref_signal,30,'measured');
        r1_signal_w_noise = awgn(r1_signal,-10,'measured');
        r2_signal_w_noise = awgn(r2_signal,-10,'measured');
        r3_signal_w_noise = awgn(r3_signal,-10,'measured');
        rref_signal_w_noise = awgn(rref_signal,-10,'measured');

        match_filter = ref_signal_w_noise(end:-1:1).';
        match_filter = match_filter';

        r1_signal_w_noise_filtered = ifft(fft(r1_signal_w_noise).*fft(match_filter));
        r2_signal_w_noise_filtered = ifft(fft(r2_signal_w_noise).*fft(match_filter));
        r3_signal_w_noise_filtered = ifft(fft(r3_signal_w_noise).*fft(match_filter));
        rref_signal_w_noise_filtered = ifft(fft(rref_signal_w_noise).*fft(match_filter));


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
        rms_errors(i,p) = rms;
        p = p + 1;
    end
    i = i + 1;
end

r = r_start:r_div:r_finish;
theta = 0:pi/theta_div:pi;
figure;
imagesc(theta,r,rms_errors);
colorbar
caxis([0 200])
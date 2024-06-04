close all;
clear all;
 clc;

frame_time_length = 0.1;
starting_time = 1e-3;
signal_time_length = 1e-3;
fs = 8e6;
snr = 10;
y = sensor_pos(2000,pi/6,0,0,"Hexagon");
target_pos = [4678.97,8799.78];

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
figure;
plot(sig(2,:))

signals = filtered_signals;


light_speed = physconst("LightSpeed");
time_values = zeros(1,length(signals(:,1)));
for i = 1:length(signals(:,1))
    [pk,idx] = max(signals(i,:));
    t_m = idx/fs;
    time_values(i) = t_m;
end

A_matrix = zeros(length(y(:,1))-1,3);
b_matrix = zeros(length(y(:,1))-1,1);

for i = 2:length(y(:,1))
A_matrix((i-1),:) = [y(1,:)-y(i,:),-light_speed*(time_values(i)-time_values(1))];
b_matrix((i-1),:) = light_speed^2*(time_values(i)-time_values(1))^2 + y(1,1)^2 + y(1,2)^2 - y(i,1)^2 - y(i,2)^2;
end
A_matrix = 2*A_matrix;
sol = pinv(A_matrix)*b_matrix;
x = sol(1);
y = sol(2);
result = [x,y]
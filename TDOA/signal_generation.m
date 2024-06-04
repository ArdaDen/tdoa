function y = signal_generation(frame_time_length,signal_time_length,starting_time,fs,sensor_pos,target_pos,snr)
frame_signal = randn(1,round(frame_time_length*fs));
signal = frame_signal(ceil(starting_time*fs):ceil(starting_time*fs + signal_time_length*fs));
match_filter = signal(end:-1:1).';
match_filter = match_filter';
filtered_signals = zeros(length(sensor_pos(1,:)),length(signal));
for i = 1:length(sensor_pos(:,1))
    t_i = norm(sensor_pos(i,:)-target_pos)/physconst("LightSpeed");
    sample_delay = ceil(t_i*fs);
    signal_i = frame_signal((ceil(starting_time*fs) + sample_delay):(ceil(starting_time*fs + signal_time_length*fs) + sample_delay));
    signal_i_w_noise = awgn(signal_i,snr,"measured");
    signal_i_w_noise_filtered = ifft(fft(match_filter).*fft(signal_i_w_noise));
    filtered_signals(i,:) = signal_i_w_noise_filtered(end:-1:1);
end
y = filtered_signals;
end
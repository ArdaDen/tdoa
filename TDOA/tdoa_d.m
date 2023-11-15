close all;
clear;
clc;

%% Input section
fs = 1e5; % Sampling rate
t = 0:1/fs:(0.1-1/fs); % Time vector
pfa = 1e-6; % Probability of false alarm
light_speed = physconst("LightSpeed"); % Speed of light
r = 60000; % Radius of sensor elements
theta = pi/3; % Angle of two sided sensor elements
emitter = [0,0,0]; % Emmitter location
number_of_sensors = 6; % Number of sensors 
ref = [r/30,0,0]; % Signal model collection location
receiver_ref = [-r/3,0,0]; % Reference receiver location
receiver1 = [r*cos(theta),r*sin(theta),0]; % Receiver 1 location
receiver2 = [r,0,0]; % Receiver 2 location
receiver3 = [r*cos(theta),-r*sin(theta),0]; % Receiver 3 location
receiver4 = [-r,0,0]; % Receiver 4 location
targets = [3086.9,-5508.74,506.89];...
        %    20489.98,-3097.6,678.3]; % Target locations

%% Delay and range calculations
t_ref = sqrt((emitter(1)-ref(1))^2+(emitter(2)-ref(2))^2+(emitter(3)-ref(3))^2)/light_speed; % Delay in reference receiver
target_ranges = zeros(1,length(targets(:,1))); % Target ranges
rref_delays = zeros(1,length(targets(:,1))); % Delays in Receiver 1
r1_delays = zeros(1,length(targets(:,1))); % Delays in Receiver 1
r2_delays = zeros(1,length(targets(:,1))); % Delays in Receiver 2
r3_delays = zeros(1,length(targets(:,1))); % Delays in Receiver 3
r4_delays = zeros(1,length(targets(:,1))); % Delays in Receiver 3

for i = 1:length(targets(:,1))
    
target_ranges_i = sqrt(targets(i,1)^2 + targets(i,2)^2 + targets(i,3)^2);

t_r_i = (sqrt((targets(i,1)-receiver_ref(1))^2+(targets(i,2)-receiver_ref(2))^2+(targets(i,3)-receiver_ref(3))^2)...
    + sqrt((emitter(1)-targets(i,1))^2+(emitter(2)-targets(i,2))^2+(emitter(3)-targets(i,3))^2))/light_speed;
t_1_i = (sqrt((targets(i,1)-receiver1(1))^2+(targets(i,2)-receiver1(2))^2+(targets(i,3)-receiver1(3))^2)...
    + sqrt((emitter(1)-targets(i,1))^2+(emitter(2)-targets(i,2))^2+(emitter(3)-targets(i,3))^2))/light_speed;
t_2_i = (sqrt((targets(i,1)-receiver2(1))^2+(targets(i,2)-receiver2(2))^2+(targets(i,3)-receiver2(3))^2)...
    + sqrt((emitter(1)-targets(i,1))^2+(emitter(2)-targets(i,2))^2+(emitter(3)-targets(i,3))^2))/light_speed;
t_3_i = (sqrt((targets(i,1)-receiver3(1))^2+(targets(i,2)-receiver3(2))^2+(targets(i,3)-receiver3(3))^2)...
    + sqrt((emitter(1)-targets(i,1))^2+(emitter(2)-targets(i,2))^2+(emitter(3)-targets(i,3))^2))/light_speed;
t_4_i = (sqrt((targets(i,1)-receiver4(1))^2+(targets(i,2)-receiver4(2))^2+(targets(i,3)-receiver4(3))^2)...
    + sqrt((emitter(1)-targets(i,1))^2+(emitter(2)-targets(i,2))^2+(emitter(3)-targets(i,3))^2))/light_speed;

target_ranges(i) = target_ranges_i;
rref_delays(i) = t_r_i;
r1_delays(i) = t_1_i;
r2_delays(i) = t_2_i;
r3_delays(i) = t_3_i;
r4_delays(i) = t_4_i;
end

%% Constructing the received signals
filename = "C:\Users\user\Desktop\staj\Recorded sounds\SDRuno_20230815_110704Z_102389kHz.wav"; % FM signal
[y,fs] = audioread(filename); % Reading the signal
y = y(:,1) + 1i*y(:,2); % IQ data

y_selec = y(50000:55000)';

ref_delay_sample = ceil(t_ref*fs);
ref_delay_sample_error = (ref_delay_sample-t_ref*fs)/fs;

r1_delay_samples = zeros(1,length(targets(:,1)));
r1_delay_samples_errors = zeros(1,length(targets(:,1)));
r2_delay_samples = zeros(1,length(targets(:,1)));
r2_delay_samples_errors = zeros(1,length(targets(:,1)));
r3_delay_samples = zeros(1,length(targets(:,1)));
r3_delay_samples_errors = zeros(1,length(targets(:,1)));
r4_delay_samples = zeros(1,length(targets(:,1)));
r4_delay_samples_errors = zeros(1,length(targets(:,1)));
rref_delay_samples = zeros(1,length(targets(:,1)));
rref_delay_samples_errors = zeros(1,length(targets(:,1)));

for k = 1:length(targets(:,1))

    r1_delay_sample_k = ceil(r1_delays(k)*fs);
    r1_delay_sample_error_k = (r1_delay_sample_k-r1_delays(k)*fs)/fs;
    r2_delay_sample_k = ceil(r2_delays(k)*fs);
    r2_delay_sample_error_k = (r2_delay_sample_k-r2_delays(k)*fs)/fs;
    r3_delay_sample_k = ceil(r3_delays(k)*fs);
    r3_delay_sample_error_k = (r3_delay_sample_k-r3_delays(k)*fs)/fs;
    r4_delay_sample_k = ceil(r4_delays(k)*fs);
    r4_delay_sample_error_k = (r4_delay_sample_k-r4_delays(k)*fs)/fs;
    rref_delay_sample_k = ceil(rref_delays(k)*fs);
    rref_delay_sample_error_k = (rref_delay_sample_k-rref_delays(k)*fs)/fs;
    r1_delay_samples(k) = r1_delay_sample_k;
    r1_delay_samples_errors(i) = r1_delay_sample_error_k;
    r2_delay_samples(k) = r2_delay_sample_k;
    r2_delay_samples_errors(i) = r2_delay_sample_error_k;
    r3_delay_samples(k) = r3_delay_sample_k;
    r3_delay_samples_errors(i) = r3_delay_sample_error_k;
    r4_delay_samples(k) = r4_delay_sample_k;
    r4_delay_samples_errors(i) = r4_delay_sample_error_k;
    rref_delay_samples(k) = rref_delay_sample_k;
    rref_delay_samples_errors(i) = rref_delay_sample_error_k;

end

ref_signal = [zeros(1,ref_delay_sample),y_selec,zeros(1,length(t)-ref_delay_sample-length(y_selec))];
r1_signal = zeros(1,length(t));
r2_signal = zeros(1,length(t));
r3_signal = zeros(1,length(t));
r4_signal = zeros(1,length(t));
rref_signal = zeros(1,length(t));

for a = 1:length(targets(:,1))
    r1_signal_a = [zeros(1,r1_delay_samples(a)),y_selec,zeros(1,length(t)-r1_delay_samples(a)-length(y_selec))];
    r2_signal_a = [zeros(1,r2_delay_samples(a)),y_selec,zeros(1,length(t)-r2_delay_samples(a)-length(y_selec))];
    r3_signal_a = [zeros(1,r3_delay_samples(a)),y_selec,zeros(1,length(t)-r3_delay_samples(a)-length(y_selec))];
    r4_signal_a = [zeros(1,r4_delay_samples(a)),y_selec,zeros(1,length(t)-r4_delay_samples(a)-length(y_selec))];
    rref_signal_a = [zeros(1,rref_delay_samples(a)),y_selec,zeros(1,length(t)-rref_delay_samples(a)-length(y_selec))];
    r1_signal = r1_signal + r1_signal_a;
    r2_signal = r2_signal + r2_signal_a;
    r3_signal = r3_signal + r3_signal_a;
    r4_signal = r4_signal + r4_signal_a;
    rref_signal = rref_signal + rref_signal_a;
end

ref_signal_w_noise = awgn(ref_signal,30,'measured');
r1_signal_w_noise = awgn(r1_signal,0,'measured');
r2_signal_w_noise = awgn(r2_signal,0,'measured');
r3_signal_w_noise = awgn(r3_signal,0,'measured');
r4_signal_w_noise = awgn(r4_signal,0,'measured');
rref_signal_w_noise = awgn(rref_signal,0,'measured');

figure;
tiledlayout(6,2);
ax1 = nexttile;
plot(t,abs(ref_signal));
ax2 = nexttile;
plot(t,abs(ref_signal_w_noise));
ax3 = nexttile;
plot(t,abs(r1_signal));
ax4 = nexttile;
plot(t,abs(r1_signal_w_noise));
ax5 = nexttile;
plot(t,abs(r2_signal));
ax6 = nexttile;
plot(t,abs(r2_signal_w_noise));
ax7 = nexttile;
plot(t,abs(r3_signal));
ax8 = nexttile;
plot(t,abs(r3_signal_w_noise));
ax9 = nexttile;
plot(t,abs(r4_signal));
ax10 = nexttile;
plot(t,abs(r4_signal_w_noise));
ax11 = nexttile;
plot(t,abs(rref_signal));
ax12 = nexttile;
plot(t,abs(rref_signal_w_noise));

match_filter = ref_signal_w_noise(end:-1:1).';
match_filter = match_filter';
r1_signal_w_noise_filtered = conv(match_filter,r1_signal_w_noise);
r1_signal_w_noise_filtered = r1_signal_w_noise_filtered(round(length(r1_signal_w_noise_filtered)/2):end);
r2_signal_w_noise_filtered = conv(match_filter,r2_signal_w_noise);
r2_signal_w_noise_filtered = r2_signal_w_noise_filtered(round(length(r2_signal_w_noise_filtered)/2):end);
r3_signal_w_noise_filtered = conv(match_filter,r3_signal_w_noise);
r3_signal_w_noise_filtered = r3_signal_w_noise_filtered(round(length(r3_signal_w_noise_filtered)/2):end);
r4_signal_w_noise_filtered = conv(match_filter,r4_signal_w_noise);
r4_signal_w_noise_filtered = r4_signal_w_noise_filtered(round(length(r4_signal_w_noise_filtered)/2):end);
rref_signal_w_noise_filtered = conv(match_filter,rref_signal_w_noise);
rref_signal_w_noise_filtered = rref_signal_w_noise_filtered(round(length(rref_signal_w_noise_filtered)/2):end);

figure;
tiledlayout(5,1);
ax13 = nexttile;
plot(t,abs(r1_signal_w_noise_filtered));
ax14 = nexttile;
plot(t,abs(r2_signal_w_noise_filtered));
ax15 = nexttile;
plot(t,abs(r3_signal_w_noise_filtered));
ax16 = nexttile;
plot(t,abs(r4_signal_w_noise_filtered));
ax17 = nexttile;
plot(t,abs(rref_signal_w_noise_filtered));

guard_cell_number = 3000;
reference_cell_number = 50;
cfar_vec = [ones(1,reference_cell_number),zeros(1,guard_cell_number),ones(1,reference_cell_number)];
cfar_vec = cfar_vec/sum(cfar_vec);

N = 2*reference_cell_number;
alpha = N*(pfa^(-1/N)-1);

cfar_t_r1 = alpha*conv(abs(r1_signal_w_noise_filtered),cfar_vec,'same');
cfar_t_r2 = alpha*conv(abs(r2_signal_w_noise_filtered),cfar_vec,'same');
cfar_t_r3 = alpha*conv(abs(r3_signal_w_noise_filtered),cfar_vec,'same');
cfar_t_r4 = alpha*conv(abs(r4_signal_w_noise_filtered),cfar_vec,'same');
cfar_t_rref = alpha*conv(abs(rref_signal_w_noise_filtered),cfar_vec,'same');

figure;
tiledlayout(5,2);
ax18 = nexttile;
plot(t,cfar_t_r1);
hold on;
plot(t,abs(r1_signal_w_noise_filtered));
ax19 = nexttile;
plot(t,20*log10(cfar_t_r1));
hold on;
plot(t,20*log10(abs(r1_signal_w_noise_filtered)));
ax20 = nexttile;
plot(t,cfar_t_r2);
hold on;
plot(t,abs(r2_signal_w_noise_filtered));
ax21 = nexttile;
plot(t,20*log10(cfar_t_r2));
hold on;
plot(t,20*log10(abs(r2_signal_w_noise_filtered)));
ax22 = nexttile;
plot(t,cfar_t_r3);
hold on;
plot(t,abs(r3_signal_w_noise_filtered));
ax23 = nexttile;
plot(t,20*log10(cfar_t_r3));
hold on;
plot(t,20*log10(abs(r3_signal_w_noise_filtered)));
ax24 = nexttile;
plot(t,cfar_t_r4);
hold on;
plot(t,abs(r4_signal_w_noise_filtered));
ax25 = nexttile;
plot(t,20*log10(cfar_t_r4));
hold on;
plot(t,20*log10(abs(r4_signal_w_noise_filtered)));
ax26 = nexttile;
plot(t,cfar_t_rref);
hold on;
plot(t,abs(rref_signal_w_noise_filtered));
ax27 = nexttile;
plot(t,20*log10(cfar_t_rref));
hold on;
plot(t,20*log10(abs(rref_signal_w_noise_filtered)));

r1_threshold = r1_signal_w_noise_filtered.*(abs(r1_signal_w_noise_filtered)>cfar_t_r1);
r2_threshold = r2_signal_w_noise_filtered.*(abs(r2_signal_w_noise_filtered)>cfar_t_r2);
r3_threshold = r3_signal_w_noise_filtered.*(abs(r3_signal_w_noise_filtered)>cfar_t_r3);
r4_threshold = r4_signal_w_noise_filtered.*(abs(r4_signal_w_noise_filtered)>cfar_t_r4);
rref_threshold = rref_signal_w_noise_filtered.*(abs(rref_signal_w_noise_filtered)>cfar_t_rref);
ri_threshold = [r1_threshold;r2_threshold;r3_threshold;r4_threshold;rref_threshold];


t_values = {};
for i = 1:number_of_sensors-1
    ind = 1;
    ri_t_indexes = {};
    ri_detected = abs(ri_threshold(i,:))>0;
    ri_detected_indexes = strfind(ri_detected,1);
    for l = 1:length(ri_detected_indexes)-1
        if ri_detected_indexes(l+1)-ri_detected_indexes(l) > 50
            ri_t_indexes{end+1} = ri_detected_indexes(ind:l);
            ind = l + 1;
        end
        if l == length(ri_detected_indexes)-1
            sum = 0;
            for z = 1:length(ri_t_indexes)
                sum = sum + length(cell2mat(ri_t_indexes(z)));
            end
            if sum ~= length(ri_detected_indexes)
                ri_t_indexes{end+1} = ri_detected_indexes(ind:end);
            end
            if length(ri_t_indexes) == 0
                fprintf("Sensor ",i," could not detect anything.")
                continue
            end
            t_values_ri = zeros(1,length(ri_t_indexes));
            for m = 1:length(ri_t_indexes)
                array = cell2mat(ri_t_indexes(m));
                target_indexes = ri_threshold(i,[array]);
                [pk,idx] = max(abs(target_indexes));
                ind_of_max = array(idx);
                t_value_m = ind_of_max/fs;
                t_values_ri(m) = t_value_m;
            end
            t_values{1,i} = t_values_ri;
            break
        else
            continue
        end
    end
end


if length(cell2mat(t_values(1))) ~= length(cell2mat(t_values(2))) || length(cell2mat(t_values(1))) ~= length(cell2mat(t_values(3))) || length(cell2mat(t_values(1))) ~= length(cell2mat(t_values(4))) || length(cell2mat(t_values(1))) ~= length(cell2mat(t_values(5)))
    disp("False Alarm or Undetected Target")
end

r1_t_values_tdoa_ref = cell2mat(t_values(1));
r2_t_values_tdoa_ref = cell2mat(t_values(2));
r3_t_values_tdoa_ref = cell2mat(t_values(3));
r4_t_values_tdoa_ref = cell2mat(t_values(4));
rref_t_values_tdoa_ref = cell2mat(t_values(5));

r1_t_values_delay = r1_t_values_tdoa_ref + t_ref;
r2_t_values_delay = r2_t_values_tdoa_ref + t_ref;
r3_t_values_delay = r3_t_values_tdoa_ref + t_ref;
r4_t_values_delay = r4_t_values_tdoa_ref + t_ref;
rref_t_values_delay = rref_t_values_tdoa_ref + t_ref;

A_matrix = 2*[receiver_ref-receiver1,-light_speed*(r1_t_values_delay(1)-rref_t_values_delay(1));...
    receiver_ref-receiver2,-light_speed*(r2_t_values_delay(1)-rref_t_values_delay(1));...
    receiver_ref-receiver3,-light_speed*(r3_t_values_delay(1)-rref_t_values_delay(1));...
    receiver_ref-receiver4,-light_speed*(r4_t_values_delay(1)-rref_t_values_delay(1))];
b_matrix = [light_speed^2*(r1_t_values_delay(1)-rref_t_values_delay(1))^2 + receiver_ref(1)^2 + receiver_ref(2)^2 + receiver_ref(3)^2 - receiver1(1)^2 - receiver1(2)^2 - receiver1(3)^2;...
    light_speed^2*(r2_t_values_delay(1)-rref_t_values_delay(1))^2 + receiver_ref(1)^2 + receiver_ref(2)^2 + receiver_ref(3)^2 - receiver2(1)^2 - receiver2(2)^2 - receiver2(3)^2;...
    light_speed^2*(r3_t_values_delay(1)-rref_t_values_delay(1))^2 + receiver_ref(1)^2 + receiver_ref(2)^2 + receiver_ref(3)^2 - receiver3(1)^2 - receiver3(2)^2 - receiver3(3)^2;...
    light_speed^2*(r4_t_values_delay(1)-rref_t_values_delay(1))^2 + receiver_ref(1)^2 + receiver_ref(2)^2 + receiver_ref(3)^2 - receiver4(1)^2 - receiver4(2)^2 - receiver4(3)^2];
sol = pinv(A_matrix)*b_matrix;
x = sol(1);
y = sol(2);
z = sol(3);

error = sqrt((x-targets(1,1))^2 + (y-targets(1,2))^2 + (z-targets(1,3))^2);

cal_range = sqrt(x^2 + y^2 + z^2);
fprintf("X coordinate of the target = %d metres.\n" , x)
fprintf("Y coordinate of the target = %d metres.\n" , y)
fprintf("Z coordinate of the target = %d metres.\n" , z)
fprintf("Error in range = %d metres.\n" , error)

fprintf("Calculated range of target from center receiver = %d metres.\n" , cal_range);
fprintf("True range of target from center receiver = %d metres.\n" , target_ranges(1));

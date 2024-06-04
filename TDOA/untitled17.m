close all;
clear all;
clc;
target_pos = [5000,15000];
fs = 8e6;
snr = 10;
min_range = 10;
max_range = 1e6;
num_of_range_bins = 50;
min_degree = 5;
num_of_elements = 4;
max_degree = 360/num_of_elements;
num_of_degree_bins = 40;
offset_range = 0;
offset_angle = 0;
shape = "Circle";
frame_time_length = 0.1;
signal_time_length = 1e-3;
starting_time = 1e-3;

rms_vec = zeros(length(logspace(log10(min_range),log10(max_range),num_of_range_bins)),length(linspace(min_degree,max_degree,num_of_degree_bins)));
count1 = 1;
count2 = 1;
for radius = logspace(log10(min_range),log10(max_range),num_of_range_bins)
    for theta = (linspace(min_degree,max_degree,num_of_degree_bins)*pi/180)
        y = sensor_pos(radius,theta,offset_range,offset_angle*pi/180,shape,num_of_elements);
        signals = signal_generation(frame_time_length,signal_time_length,starting_time,fs,y,target_pos,snr);
        result = time_dif_of_arrival(signals,fs,y);
        rms = sqrt(((result(1)-target_pos(1))^2 + (result(2)-target_pos(2))^2)/2);
        rms_vec(count2,count1) = rms;
        count1 = count1 + 1;
    end
    count1 = 1;
    count2 = count2 + 1;
end
radius = logspace(log10(min_range),log10(max_range),num_of_range_bins);
theta = linspace(min_degree,max_degree,num_of_degree_bins);
figure;
imagesc(theta,radius,rms_vec);
colorbar
set(gca,"ColorScale","log")
xlabel("Degree")
ylabel("Range(m)")
title("RMS values according to ranges and degrees of the sensors")
caxis([1,10000])
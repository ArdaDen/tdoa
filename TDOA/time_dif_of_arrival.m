function result = time_dif_of_arrival(signals,fs,sensor_pos)
light_speed = physconst("LightSpeed");
time_values = zeros(1,length(signals(:,1)));
for i = 1:length(signals(:,1))
    [pk,idx] = max(signals(i,:));
    t_m = (idx-1)/fs;
    time_values(i) = t_m;
end

A_matrix = zeros(length(sensor_pos(:,1))-1,3);
b_matrix = zeros(length(sensor_pos(:,1))-1,1);

for i = 2:length(sensor_pos(:,1))
A_matrix((i-1),:) = [sensor_pos(1,:)-sensor_pos(i,:),-light_speed*(time_values(i)-time_values(1))];
b_matrix((i-1),:) = light_speed^2*(time_values(i)-time_values(1))^2 + sensor_pos(1,1)^2 + sensor_pos(1,2)^2 - sensor_pos(i,1)^2 - sensor_pos(i,2)^2;
end

A_matrix = 2*A_matrix;
pseudo_inv = (((A_matrix.')*A_matrix))\(A_matrix.');
sol = pseudo_inv*b_matrix;
x_target = sol(1);
y_target = sol(2);
result = [x_target,y_target];
end
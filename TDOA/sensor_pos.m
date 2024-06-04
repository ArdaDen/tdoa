function y = sensor_pos(r,yaw_angle,offset_range,offset_angle,name,num_of_elements)
if name == "Triangle"
    y = triangle(r,yaw_angle,offset_range,offset_angle);
elseif name == "Rectangle"
    y = rectangle(r,yaw_angle,offset_range,offset_angle);
elseif name == "Pentagon"
    y = pentagon(r,yaw_angle,offset_range,offset_angle);
elseif name == "Hexagon"
    y = hexagon(r,yaw_angle,offset_range,offset_angle);
elseif name == "Circle"
    y = circle(r,yaw_angle,num_of_elements,offset_range,offset_angle);
else 
    error("Not defined")
end

end

function y = triangle(r,yaw_angle,offset_range,offset_angle)
sensor_1_pos = [r*cos(yaw_angle + pi/2) + offset_range*cos(offset_angle),r*sin(yaw_angle + pi/2) + offset_range*sin(offset_angle)];
sensor_2_pos = [r*cos(yaw_angle + pi + pi/6) + offset_range*cos(offset_angle),r*sin(yaw_angle + pi + pi/6) + offset_range*sin(offset_angle)];
sensor_3_pos = [r*cos(yaw_angle + 2*pi - pi/6) + offset_range*cos(offset_angle),r*sin(yaw_angle + 2*pi - pi/6) + offset_range*sin(offset_angle)];
y = [sensor_1_pos;sensor_2_pos;sensor_3_pos];
end

function y = rectangle(r,yaw_angle,offset_range,offset_angle)
sensor_1_pos = [r*cos(yaw_angle + pi/4) + offset_range*cos(offset_angle),r*sin(yaw_angle + pi/4) + offset_range*sin(offset_angle)];
sensor_2_pos = [r*cos(yaw_angle + 3*pi/4) + offset_range*cos(offset_angle),r*sin(yaw_angle + 3*pi/4) + offset_range*sin(offset_angle)];
sensor_3_pos = [r*cos(yaw_angle + 5*pi/4) + offset_range*cos(offset_angle),r*sin(yaw_angle + 5*pi/4) + offset_range*sin(offset_angle)];
sensor_4_pos = [r*cos(yaw_angle + 7*pi/4) + offset_range*cos(offset_angle),r*sin(yaw_angle + 7*pi/4) + offset_range*sin(offset_angle)];
y = [sensor_1_pos;sensor_2_pos;sensor_3_pos;sensor_4_pos];
end

function y = pentagon(r,yaw_angle,offset_range,offset_angle)
sensor_1_pos = [r*cos(yaw_angle + deg2rad(18)) + offset_range*cos(offset_angle),r*sin(yaw_angle + deg2rad(18)) + offset_range*sin(offset_angle)];
sensor_2_pos = [r*cos(yaw_angle + pi/2) + offset_range*cos(offset_angle),r*sin(yaw_angle + pi/2) + offset_range*sin(offset_angle)];
sensor_3_pos = [r*cos(yaw_angle + deg2rad(162)) + offset_range*cos(offset_angle),r*sin(yaw_angle + deg2rad(162)) + offset_range*sin(offset_angle)];
sensor_4_pos = [r*cos(yaw_angle + deg2rad(234)) + offset_range*cos(offset_angle),r*sin(yaw_angle + deg2rad(234)) + offset_range*sin(offset_angle)];
sensor_5_pos = [r*cos(yaw_angle + deg2rad(306)) + offset_range*cos(offset_angle),r*sin(yaw_angle + deg2rad(306)) + offset_range*sin(offset_angle)];
y = [sensor_1_pos;sensor_2_pos;sensor_3_pos;sensor_4_pos;sensor_5_pos];
end

function y = hexagon(r,yaw_angle,offset_range,offset_angle)
sensor_1_pos = [r*cos(yaw_angle) + offset_range*cos(offset_angle),r*sin(yaw_angle) + offset_range*sin(offset_angle)];
sensor_2_pos = [r*cos(yaw_angle + pi/3) + offset_range*cos(offset_angle),r*sin(yaw_angle + pi/3) + offset_range*sin(offset_angle)];
sensor_3_pos = [r*cos(yaw_angle + 2*pi/3) + offset_range*cos(offset_angle),r*sin(yaw_angle + 2*pi/3) + offset_range*sin(offset_angle)];
sensor_4_pos = [r*cos(yaw_angle + 3*pi/3) + offset_range*cos(offset_angle),r*sin(yaw_angle + 3*pi/3) + offset_range*sin(offset_angle)];
sensor_5_pos = [r*cos(yaw_angle + 4*pi/3) + offset_range*cos(offset_angle),r*sin(yaw_angle + 4*pi/3) + offset_range*sin(offset_angle)];
sensor_6_pos = [r*cos(yaw_angle + 5*pi/3) + offset_range*cos(offset_angle),r*sin(yaw_angle + 5*pi/3) + offset_range*sin(offset_angle)];
y = [sensor_1_pos;sensor_2_pos;sensor_3_pos;sensor_4_pos;sensor_5_pos;sensor_6_pos];
end

function y = circle(r,theta,num_of_elements,offset_range,offset_angle)
sensor_poses = zeros(num_of_elements,2);
c = 1;
for i = 1:fix(num_of_elements/2)
    sensor_pos_i = [r*cos(theta*i) + offset_range*cos(offset_angle),r*sin(theta*i) + offset_range*sin(offset_angle);...
        r*cos((-theta)*i) + offset_range*cos(offset_angle),r*sin((-theta)*i) + offset_range*sin(offset_angle)];
    sensor_poses(c:c+1,:) = sensor_pos_i;
    c = c + 2;
end
if rem(num_of_elements,2) == 1
    sensor_poses(num_of_elements,:) = [r*cos(theta*(i+1)) + offset_range*cos(offset_angle),r*sin(theta*(i+1)) + offset_range*sin(offset_angle)];
end
y = sensor_poses;
end
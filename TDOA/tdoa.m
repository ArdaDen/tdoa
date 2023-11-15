close all;
clear;
clc;

Range_map = 5000;
dt = 1;
receiver1 = [-1000,-1000*sqrt(3),0];
receiver2 = [1000,-1000*sqrt(3),0];
receiver3 = [0,0,0];
receiver4 = [0,1000,0];
target = [4868.88,-4787.36,3745.3];
target_range = sqrt(target(1)^2 + target(2)^2 + target(3)^2);
t_1 = sqrt((target(1)-receiver1(1))^2+(target(2)-receiver1(2))^2+(target(3)-receiver1(3))^2)/physconst("LightSpeed");
t_2 = sqrt((target(1)-receiver2(1))^2+(target(2)-receiver2(2))^2+(target(3)-receiver2(3))^2)/physconst("LightSpeed");
t_3 = sqrt((target(1)-receiver3(1))^2+(target(2)-receiver3(2))^2+(target(3)-receiver3(3))^2)/physconst("LightSpeed");
t_4 = sqrt((target(1)-receiver4(1))^2+(target(2)-receiver4(2))^2+(target(3)-receiver4(3))^2)/physconst("LightSpeed");
x_t = -Range_map:dt:Range_map;
y_t = -Range_map:dt:Range_map;
[xx,yy] = meshgrid(x_t,y_t);
func1 = sqrt(abs((physconst("LightSpeed")*t_1)^2-(xx-receiver1(1)).^2-(yy-receiver1(2)).^2)) + receiver1(3);
func2 = sqrt(abs((physconst("LightSpeed")*t_2)^2-(xx-receiver2(1)).^2-(yy-receiver2(2)).^2)) + receiver2(3);
func3 = sqrt(abs((physconst("LightSpeed")*t_3)^2-(xx-receiver3(1)).^2-(yy-receiver3(2)).^2)) + receiver3(3);
func4 = sqrt(abs((physconst("LightSpeed")*t_4)^2-(xx-receiver4(1)).^2-(yy-receiver4(2)).^2)) + receiver4(3);

% figure;
% surf(func1-func2,"EdgeColor","none");
% hold on;
% surf(func2-func3,"EdgeColor","none");
% surf(abs(func1-func2),"EdgeColor","none")
% surf(func2,"EdgeColor","none");
% surf(func3,"EdgeColor","none");
% plot3(target(1),target(2),target(3),".");
% plot3(receiver1(1),receiver1(2),receiver1(3),"*");
% plot3(receiver2(1),receiver2(2),receiver2(3),"*");
% plot3(receiver3(1),receiver3(2),receiver3(3),"*");

% figure;
% surf(func1-func2,"EdgeColor","none");
% hold on;
% surf(func2-func3,"EdgeColor","none");
% surf(abs(func2-func3),"EdgeColor","none")
% surf(func2,"EdgeColor","none");
% surf(func3,"EdgeColor","none");
% plot3(target(1),target(2),target(3),".");
% plot3(receiver1(1),receiver1(2),receiver1(3),"*");
% plot3(receiver2(1),receiver2(2),receiver2(3),"*");
% plot3(receiver3(1),receiver3(2),receiver3(3),"*");

% figure;
% surf(func1-func2,"EdgeColor","none");
% hold on;
% surf(func2-func3,"EdgeColor","none");
% surf(abs(func1-func3),"EdgeColor","none")
% surf(func2,"EdgeColor","none");
% surf(func3,"EdgeColor","none");
% plot3(target(1),target(2),target(3),".");
% plot3(receiver1(1),receiver1(2),receiver1(3),"*");
% plot3(receiver2(1),receiver2(2),receiver2(3),"*");
% plot3(receiver3(1),receiver3(2),receiver3(3),"*");

% figure;
% surf(func1-func2,"EdgeColor","none");
% hold on;
% surf(func2-func3,"EdgeColor","none");
% surf(abs(func1-func3) + abs(func2-func3) + abs(func1-func2) + abs(func4-func2) + abs(func4-func3) + abs(func4-func1),"EdgeColor","none");
% hold on;
% surf(func2,"EdgeColor","none");
% surf(func3,"EdgeColor","none");
% plot3(target(1),target(2),target(3),"*");
% plot3(receiver1(1),receiver1(2),receiver1(3),"*");
% plot3(receiver2(1),receiver2(2),receiver2(3),"*");
% plot3(receiver3(1),receiver3(2),receiver3(3),"*");

f = abs(func1-func3) + abs(func2-func3) + abs(func1-func2) + abs(func4-func2) + abs(func4-func3) + abs(func4-func1);
x = min(min(f));
[y,x] = find(f == x);
y = y*dt-Range_map;
x = x*dt-Range_map;
z = sqrt((physconst("LightSpeed")*t_1)^2-(x-receiver1(1)).^2-(y-receiver1(2)).^2) + receiver1(3);
error = sqrt((x-target(1))^2 + (y-target(2))^2 + (z-target(3))^2);
cal_range = sqrt(x^2 + y^2 + z^2);
fprintf("X coordinate of the target = %d metres.\n" , x)
fprintf("Y coordinate of the target = %d metres.\n" , y)
fprintf("Z coordinate of the target = %d metres.\n" , z)
fprintf("Error in range = %d metres.\n" , error)
fprintf("Calculated range of target from center receiver = %d metres.\n" , cal_range);
fprintf("True range of target from center receiver = %d metres.\n" , target_range);

% figure;
% plot3(target(1),target(2),target(3),"*");
% hold on;
% plot3(receiver1(1),receiver1(2),receiver1(3),"*");
% plot3(receiver2(1),receiver2(2),receiver2(3),"*");
% plot3(receiver3(1),receiver3(2),receiver3(3),"*");
% plot3(receiver4(1),receiver4(2),receiver4(3),"*");

% figure;
% contour(xx,yy,func1);
% contour(xx,yy,func1-func2);
% hold on;
% contour(xx,yy,func2);
% contour(xx,yy,func3);
% contour(xx,yy,func2-func3);
% contour(xx,yy,func1-func3);
% plot3(target(1),target(2),target(3),".");
% plot3(receiver1(1),receiver1(2),receiver1(3),"*");
% plot3(receiver2(1),receiver2(2),receiver2(3),"*");
% plot3(receiver3(1),receiver3(2),receiver3(3),"*");

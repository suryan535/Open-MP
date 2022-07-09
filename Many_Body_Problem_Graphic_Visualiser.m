close all;
clear all;
clc;

filepath = 'Coordinates.txt';   %Coordinates of body path
fileId = fopen(filepath, 'r');  %Opening the file

line = fgetl(fileId);
w = str2double(regexp(line, '[\d.]+', 'match'));    %Width of the file

line = fgetl(fileId);
l = str2double(regexp(line, '[\d.]+', 'match'));    %Length of the file

line = fgetl(fileId);
d = str2double(regexp(line, '[\d.]+', 'match'));    %Depth of the file

line = fgetl(fileId);
N = str2double(regexp(line, '[\d.]+', 'match'));    %Number of bodies

line = fgetl(fileId);
r = str2double(regexp(line, '[\d.]+', 'match'));    %Radius of a body

line = fgetl(fileId);
t_step = str2double(regexp(line, '[\d.]+', 'match'));    %Time Step

line = fgetl(fileId);
num_step = str2double(regexp(line, '[\d.]+', 'match'));    %Number of time steps or iterations

data = zeros(3, N, num_step);

C = rand(N, 3); %Random Color for each ball

figure("Name", "Many Body Graphics", 'units','normalized','outerposition',[0 0 1 1]);
for i = 1:num_step
    clf;
    line = fgetl(fileId);
    iter_num = str2double(regexp(line, '[\d.]+', 'match'));    %Iteration number  
    
    coord = fscanf(fileId, '%f %f %f', [3, N*num_step]);    %Coordinates of the ith iteration
%     plot3(coord(1, :), coord(2, :), coord(3, :), 'o', 'LineWidth', r*8, 'MarkerSize', r*8, 'MarkerFaceColor', 'auto');  %plotting all the n bodies
    scatter3(coord(1, :), coord(2, :), coord(3, :), r*70, C, 'filled');
    hold on;
    grid on;
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title(["Many Bodies Positions at iteration = ",num2str(i)], "Color", 'r');
    axis([0, w, 0, l, 0, d]);
    hold off;
    
    currframe(i) = getframe;
    drawnow();
%     pause(0.25);
    data(:, :, i) = coord; %Storing the data for every time step in one 3d matrix
    disp("Iteration "+int2str(i)+" Data Loaded");
end

myVideo = VideoWriter('ManyBodyGraphical');
myVideo.FrameRate = 10;

open(myVideo);
writeVideo(myVideo, currframe);
close(myVideo);

fclose(fileId);
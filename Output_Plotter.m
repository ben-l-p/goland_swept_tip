clear all; close all;

load("case_output_forces.mat")

%Create vectors for X and Y directions
ang_Y = rad2deg(ang_flutter(:, 1));
pos_X = pos_flutter(1, :);

%Set divergence and flutter velocities to be infinite when they did not
%occur (v=0)
v_min_thresh = 0;                     
v_div(v_div<=v_min_thresh) = inf;
v_flutter(v_flutter<=v_min_thresh) = inf;
v_flutter_filt(v_flutter_filt<=v_min_thresh) = inf;
v_max = min(v_div, v_flutter_filt);

%Create boolean array for where flutter and divergence occur
is_flutter = zeros([length(ang_Y) length(pos_X)]);
for i = 1:length(ang_Y)
    for j = 1:length(pos_X)
        is_flutter(i, j) = (v_flutter_filt(i, j) < v_div(i, j));
    end
end

%Plot flutter, divergence and maximum speeds
subplot(1, 3, 1)
surf(pos_X, ang_Y, v_div)
xlabel("Fraction of Span to Kink");
ylabel("Sweep Angle (deg)");
zlabel("Divergence Speed (m/s)");
title("Divergence Speed");
xlim([min(pos_X) max(pos_X)]);
ylim([min(ang_Y) max(ang_Y)]);
zlim([100 500]);
clim([150 350])
view([0 90]);

subplot(1, 3, 2)
surf(pos_X, ang_Y, v_flutter_filt)
xlabel("Fraction of Span to Kink");
ylabel("Sweep Angle (deg)");
zlabel("Flutter Speed (m/s)");
title("Flutter Speed");
xlim([min(pos_X) max(pos_X)]);
ylim([min(ang_Y) max(ang_Y)]);
zlim([100 500]);
clim([150 350])
view([0 90]);

subplot(1, 3, 3)
surf(pos_X, ang_Y, v_max)
xlabel("Fraction of Span to Kink");
ylabel("Sweep Angle (deg)");
zlabel("Maximum Speed (m/s)");
title("Maximum Speed");
xlim([min(pos_X) max(pos_X)]);
ylim([min(ang_Y) max(ang_Y)]);
zlim([100 500]);
clim([150 350])
view([0 90]);

%Plot moment
figure();
for i = 1:3
    subplot(1, 3, i)
    surf(pos_X, ang_Y, moments_a(:, :, i))
    xlabel("Fraction of Span to Kink");
    ylabel("Sweep Angle (deg)");
    zlabel("Total Moment (Nm)");
    title("Total Moments");
    view([0 90]);
    colorbar
end

%Plot force
figure();
for i = 1:3
    subplot(1, 3, i)
    surf(pos_X, ang_Y, forces_a(:, :, i))
    xlabel("Fraction of Span to Kink");
    ylabel("Sweep Angle (deg)");
    zlabel("Total Force (N)");
    title("Total Forces");
    view([0 90]);
    colorbar
end

% c = colorbar;
% c.Label.String = 'Maximum Velocity (m/s)';
% clim([150 350])

% figure()
% surf(pos_X, ang_Y, is_flutter)
% view([0 90]);

% figure()
% hold on;
% surf(pos_X, ang_Y, v_div)
% surf(pos_X, ang_Y, v_flutter)
% hold off;


%shading interp
%h = gca;
%set(h,'zscale','log')



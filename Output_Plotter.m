clear all; close all; clc;

load("case_output.mat");

%Create vectors for X and Y directions
ang_Y = rad2deg(ang_flutter(:, 1));
pos_X = pos_flutter(1, :);
[~, ~, n_tstep, ~] = size(moments_a_s);

%% Forces and Moments

moments_tot = moments_g_u + moments_g_s + moments_a_u + moments_a_s;

max_moms = squeeze(max(moments_tot, [], 3));
min_moms = squeeze(min(moments_tot, [], 3));


hold on;
for i = 1:3
    subplot(1, 3, i);
    hold on;
    for j = 1:length(ang_Y)
        for k = 1:length(pos_X)
            plot((1:n_tstep).*dt, squeeze(moments_tot(j, k, :, i)));
        end
    end
    hold off;
end

%% Stability
%Set divergence and flutter velocities to be infinite when they did not
%occur (v=0)
% v_min_thresh = 0;                     
% v_div(v_div<=v_min_thresh) = inf;
% v_flutter(v_flutter<=v_min_thresh) = inf;
% v_flutter_filt(v_flutter_filt<=v_min_thresh) = inf;
% v_max = min(v_div, v_flutter_filt);
% 
% %Create boolean array for where flutter and divergence occur
% is_flutter = zeros([length(ang_Y) length(pos_X)]);
% for i = 1:length(ang_Y)
%     for j = 1:length(pos_X)
%         is_flutter(i, j) = (v_flutter_filt(i, j) < v_div(i, j));
%     end
% end
% 
% %Plot stability
% figure();
% v_vect = cat(3, v_div, v_flutter_filt, v_max);
% velocity_titles = ["Divergence", "Flutter", "Maximum"];
% for i = 1:3
%     subplot(1, 3, i)
%     surf(pos_X, ang_Y, v_vect(:, :, i))
%     xlabel("Fraction of Span to Kink");
%     ylabel("Sweep Angle (deg)");
%     zlabel("Velocity (m/s)");
%     title(velocity_titles(i));
%     view([0 90]);
%     c = colorbar;
%     c.Label.String = 'Maximum Stable Velocity (m/s)';
% end

%Plot forces and moments
% figure();
% sgtitle('Total Applied Aerodynamic Forces and Moments');
% 
% force_titles = ["X Force", "Y Force", "Z Force"];
% for i = 1:3
%     subplot(2, 3, i)
%     surf(pos_X, ang_Y, forces_a(:, :, i))
%     xlabel("Fraction of Span to Kink");
%     ylabel("Sweep Angle (deg)");
%     zlabel("Total Force (N)");
%     title(force_titles(i));    
%     view([0 90]);
%     c = colorbar;
%     c.Label.String = 'Total Applied Force (N)';
%     shading interp
% end
% 
% moment_titles = ["X Moment", "Y Moment", "Z Moment"];
% for i = 1:3
%     subplot(2, 3, i+3)
%     surf(pos_X, ang_Y, moments_a(:, :, i))
%     xlabel("Fraction of Span to Kink");
%     ylabel("Sweep Angle (deg)");
%     zlabel("Total Moment (Nm)");
%     title(moment_titles(i));
%     view([0 90]);
%     c = colorbar;
%     c.Label.String = 'Total Applied Moment (Nm)';
%     shading interp
% end


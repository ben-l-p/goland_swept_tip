clear all; close all; clc;

load("case_output.mat");

%Create vectors for X and Y directions
ang_Y = rad2deg(ang_flutter(:, 1));
pos_X = pos_flutter(1, :);
[~, ~, n_tstep, ~] = size(moments_a_s);

straight_index = find(ang_Y==0);

%% Forces and Moments

moments_tot = moments_g_u + moments_g_s + moments_a_u + moments_a_s;                                                                

max_moms = squeeze(max(moments_tot, [], 3));
min_moms = squeeze(min(moments_tot, [], 3));
abs_max_moms = squeeze(max(abs(moments_tot), [], 3));
signed_max_moms = (abs(min_moms)>abs(max_moms))...
    .*min_moms+(abs(min_moms)<abs(max_moms)).*max_moms;


norm_abs_max_moms = zeros(size(signed_max_moms));
for i = 1:3
    norm_abs_max_moms(:, :, i) = signed_max_moms(:, :, i)...
        ./signed_max_moms(straight_index, end, i);
end


figure();
sgtitle("Applied Moments during Gust Response");
mom_titles = ["Out-of-Plane", "Twisting", "Yawing"];

for i = 1:3
    subplot(2, 3, i);
    hold on
    title(mom_titles(i));
    for j = 1:length(ang_Y)
        for k = 1:length(pos_X)
            if j ~= straight_index && k ~= length(pos_X)
                plot((1:n_tstep).*dt, squeeze(moments_tot(j, k, :, i)), 'k-');
            end
        end
    end
    
    plot((1:n_tstep).*dt, squeeze(moments_tot(straight_index, length(pos_X), :, i)), 'r-', LineWidth=5);

    xlabel("Time (s)");
    ylabel("Applied Moment (Nm)")
    hold off;
end

for i = 1:3
    subplot(2, 3, i+3);
    hold on;
    surf(pos_X, ang_Y, norm_abs_max_moms(:, :, i));
    title(mom_titles(i));
    xlabel("Fraction of Span to Kink");
    ylabel("Sweep Angle (deg)");
    zlabel("Normalised Maximum Moment");
    c = colorbar;
    c.Label.String = 'Normalised Maximum Moment';
    shading interp;
    axis([min(pos_X), max(pos_X), min(ang_Y), max(ang_Y)]);
    hold off;
end

forces_a_u_padded = cat(4, forces_a_u, zeros(length(ang_Y), length(pos_X), n_tstep, 1));
forces_tot = forces_a_s + forces_g_s + forces_g_u + forces_a_u_padded;
forces_steady = squeeze(forces_tot(:, :, 1, :));
norm_forces_steady = zeros(size(forces_steady));

max_forces = squeeze(max(forces_tot, [], 3));
min_forces = squeeze(min(forces_tot, [], 3));
abs_max_forces = squeeze(max(abs(forces_tot), [], 3));
signed_max_forces = (abs(min_forces)>abs(max_forces))...
    .*min_forces+(abs(min_forces)<abs(max_forces)).*max_forces;

max_over_steady_forces = signed_max_forces./forces_steady;

for i = 1:3
    norm_forces_steady(:, :, i) = forces_steady(:, :, i)./forces_steady(straight_index, end, i);
end

figure();
sgtitle("Applied Forces during Gust Response");
force_titles = ["Chordwise", "Spanwise", "Out-of-Plane"];

for i = 1:3
    subplot(2, 3, i);
    hold on;
    title(force_titles(i));
    for j = 1:length(ang_Y)
        for k = 1:length(pos_X)
            if j ~= straight_index && k ~= length(pos_X)
                plot((1:n_tstep).*dt, squeeze(forces_tot(j, k, :, i)), 'k-');
            end
        end
    end

    plot((1:n_tstep).*dt, squeeze(forces_tot(straight_index, length(pos_X), :, i)), 'r-', LineWidth=5);

    xlabel("Time (s)");
    ylabel("Applied Force (N)")
    hold off;
end

for i = 1:3
    subplot(2, 3, i+3);
    hold on;
    surf(pos_X, ang_Y, norm_forces_steady(:, :, i));
    title(force_titles(i));
    xlabel("Fraction of Span to Kink");
    ylabel("Sweep Angle (deg)");
    zlabel("Normalised Steady Force");
    c = colorbar;
    c.Label.String = 'Normalised Steady Force';
    shading interp;
    axis([min(pos_X), max(pos_X), min(ang_Y), max(ang_Y)]);
    hold off;
end

m_o_l = abs_max_moms./cat(3, forces_steady(:, :, 3),...
    forces_steady(:, :, 3), forces_steady(:, :, 3));

norm_m_o_l = m_o_l;
for i = 1:3
    norm_m_o_l(:, :, i) = m_o_l(:, :, i)./m_o_l(straight_index, end, i);
end

% figure()
% for i = 1:3
%     subplot(1, 3, i);
%     hold on;
%     surf(pos_X, ang_Y, max_over_steady_forces(:, :, i));
%     title(force_titles(i));
%     xlabel("Fraction of Span to Kink");
%     ylabel("Sweep Angle (deg)");
%     zlabel("Maximum to Steady Force Ratio");
%     c = colorbar;
%     c.Label.String = 'Maximum to Steady Force Ratio';
%     shading interp;
%     axis([min(pos_X), max(pos_X), min(ang_Y), max(ang_Y)]);
%     hold off;
% end

% figure();
% sgtitle("Normalised Moment to Lift Ratio");
% for i = 1:3
%     subplot(1, 3, i);
%     hold on;
%     surf(pos_X, ang_Y, norm_m_o_l(:, :, i));
%     title(mom_titles(i));
%     xlabel("Fraction of Span to Kink");
%     ylabel("Sweep Angle (deg)");
%     zlabel("Moment to Lift Ratio");
%     c = colorbar;
%     c.Label.String = 'Moment to Lift Ratio';
%     shading interp;
%     axis([min(pos_X), max(pos_X), min(ang_Y), max(ang_Y)]);
%     hold off;
% end

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


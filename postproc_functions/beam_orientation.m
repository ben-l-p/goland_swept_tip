close all; clear all;   clc;
load("beam_orientation.mat");

[n_tstep, n_node] = size(pos, [1, 2]);
n_elem = (n_node - 1)/2;

%% Convert psi values to rotation matrices
R_bs = zeros(n_tstep, n_elem, 3, 3, 3);

for i_ts = 1:n_tstep
    for i_elem = 1:n_elem
        for i_en = 1:3
            R_bs(i_ts, i_elem, i_en, :, :) = psi_to_rot(squeeze(psi(i_ts, i_elem, i_en, :)));
        end
    end
end

%% Extract one rotation matrix per node
R_bs_node = zeros(n_tstep, n_node, 3, 3);

i_node = 1;
for i_elem = 1:n_elem
    for i_en = 1:3
        if i_en == 1 && i_elem ~= 1
            continue;
        end
        
        R_bs_node(:, i_node, :, :) = R_bs(:, i_elem, i_en, :, :);
        i_node = i_node + 1;
    end
end

%% Plot 3D beam with orientation triads
% i_ts = 60;
c_ref = 1.8288;
ea = 0.33;
triad_scale = 0.5;
triad_colour = 'rgb';

frames = struct('cdata',[],'colormap',[]);
v = VideoWriter("beam_orientation.avi");
open(v);

fig = figure();

for i_ts = 1:n_tstep
    clf(fig);

    hold on;

    axis equal;
    axis padded;
    xlabel("X (m)");
    ylabel("Y (m)");
    zlabel("Z (m)");
    view(0, 90);

    xlim([-2, 11]);
    ylim([-1, 17]);

    % Plot beam line
    plot3(pos(i_ts, :, 1), pos(i_ts, :, 2), pos(i_ts, :, 3), 'k-');
    
    % Plot B triad
    for i_dir = 1:3
        vect_b = u_vect(i_dir)*3;
        quiver3(0, 0, 0, vect_b(1), vect_b(2), vect_b(3), ...
            "Color", triad_colour(i_dir), LineWidth=2);
    end
    
    for i_node = 1:n_node
        % Plot beam nodes
        pos_node = squeeze(pos(i_ts, i_node, :));
        scatter3(pos_node(1), pos_node(2), pos_node(3), 100, 'k.');
    
        % Plot S triads
        for i_dir = 1:3
            u_vect_s = squeeze(R_bs_node(i_ts, i_node, :, :))*u_vect(i_dir)*triad_scale;
            quiver3(pos_node(1), pos_node(2), pos_node(3), ...
                u_vect_s(1), u_vect_s(2), u_vect_s(3),...
                "Color", triad_colour(i_dir), LineWidth=2);
        end
    
        % Plot chord
        c_part = [-c_ref*ea, c_ref*(1-ea)];
        for i_dir = 1:2
            u_vect_s = squeeze(R_bs_node(i_ts, i_node, :, :))*u_vect(2)*c_part(i_dir);
            quiver3(pos_node(1), pos_node(2), pos_node(3), ...
                u_vect_s(1), u_vect_s(2), u_vect_s(3),'k');
        end
    end
    
    hold off;
    frames(i_ts) = getframe(fig);
    writeVideo(v, frames(i_ts));
end

close(v);
    

%% Skew-symmetric operator
function a_skew = skew(a)
    a_skew = [0 -a(3) a(2); a(3), 0, -a(1); -a(2), a(1), 0];
end

%% Split phi into its angle and unit vector
function [psi_ang, psi_n] = split_psi(psi)
    psi_ang = norm(psi);
    psi_n = psi/psi_ang;
end

%% Convert psi to rotation matrix
function R_ab = psi_to_rot(psi)
    psi_s = skew(psi);
    [psi_ang, ~] = split_psi(psi);

    R_ab = eye(3) + sin(psi_ang)/psi_ang*psi_s + ...
        (1-cos(psi_ang))/psi_ang^2*psi_s*psi_s;
end

%% Create unit vector in given direction
function e = u_vect(dim)
    e = zeros(3, 1);
    e(dim) = 1;
end

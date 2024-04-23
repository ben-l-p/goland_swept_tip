t_step = 2;

ang_i = 1;
pos_i = 1;

psi = data{ang_i, pos_i}.psi;

n_elem = data{1, 1}.n_elem_surf;
R_bs = zeros(n_elem, 3, 3, 3);
rot_angs = zeros(n_elem, 3, 3);
for i_elem = 1:n_elem
    for i_en = 1:3
        R_bs(i_elem, i_en, :, :) = psi_to_rot(squeeze(psi(t_step, i_elem, i_en, :)));
        rot_angs(i_elem, i_en, :) = rot_to_angs(squeeze(R_bs(i_elem, i_en, :, :)));
    end
end

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

%% Convert rotation matrix to angles
function [theta, phi, psi] = rot_to_angs(R_ab)
    theta = asin(-R_ab(3, 1));
    psi = atan(R_ab(3, 2)/R_ab(3, 3));
    phi = atan(R_ab(2, 1)/R_ab(1, 1));
end
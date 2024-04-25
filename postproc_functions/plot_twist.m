close all;  clc;
% load("case_data.mat");

[n_ang, n_pos] = size(data);
h_ang = zeros(n_ang, n_pos);
pos_frac_h = zeros(n_pos, 1);
ang_h = zeros(n_ang, 1);


for i_ang = 1:n_ang
    for i_pos = 1:n_pos
        pos_frac_h(i_pos) = data{i_ang, i_pos}.pos_frac_h;
        ang_h(i_ang) = data{i_ang, i_pos}.ang_h;

        psi_j_elem = data{i_ang, i_pos}.psi_init;
        psi_t_elem = squeeze(data{i_ang, i_pos}.psi);
        
        n_elem_tot = data{1, 1}.n_elem_tot;
        n_node_tot = data{1, 1}.n_node_tot;
        
        node_h = data{i_ang, i_pos}.node_h;
        
        psi_j_nodal = zeros(n_elem_tot*3, 3);
        psi_t_nodal = zeros(n_elem_tot*3, 3);
        
        
        for i_elem = 1:n_elem_tot
            psi_j_nodal(3*(i_elem-1)+1, :) = psi_j_elem(i_elem, 1, :);
            psi_j_nodal(3*(i_elem-1)+2, :) = psi_j_elem(i_elem, 3, :);
            psi_j_nodal(3*(i_elem-1)+3, :) = psi_j_elem(i_elem, 2, :);
            
            psi_t_nodal(3*(i_elem-1)+1, :) = psi_t_elem(i_elem, 1, :);
            psi_t_nodal(3*(i_elem-1)+2, :) = psi_t_elem(i_elem, 3, :);
            psi_t_nodal(3*(i_elem-1)+3, :) = psi_t_elem(i_elem, 2, :);
        end
        
        R_j_nodal = zeros(n_elem_tot*3, 3, 3);
        R_t_nodal = zeros(n_elem_tot*3, 3, 3);
        R_tj_nodal = zeros(n_elem_tot*3, 3, 3);
        theta_tj = zeros(n_elem_tot*3, 1);
        psi_tj = zeros(n_elem_tot*3, 1);
        phi_tj = zeros(n_elem_tot*3, 1);
        theta_t = zeros(n_elem_tot*3, 1);
        psi_t = zeros(n_elem_tot*3, 1);
        phi_t = zeros(n_elem_tot*3, 1);
        
        for i_node = 1:n_elem_tot*3
            R_j_nodal(i_node, :, :) = psi_to_R(squeeze(psi_j_nodal(i_node, :)));
            R_t_nodal(i_node, :, :) = psi_to_R(squeeze(psi_t_nodal(i_node, :)));
            R_tj_nodal(i_node, :, :) = squeeze(R_t_nodal(i_node, :, :))*squeeze(R_j_nodal(i_node, :, :))';
        
            theta_tj(i_node) = -asind(R_tj_nodal(i_node, 3, 1));
            psi_tj(i_node) = atand(R_tj_nodal(i_node, 3, 2)/R_tj_nodal(i_node, 3, 3));
            phi_tj(i_node) = atand(R_tj_nodal(i_node, 2, 1)/R_tj_nodal(i_node, 1, 1));
        
            theta_t(i_node) = -asind(R_t_nodal(i_node, 3, 1));
            psi_t(i_node) = atand(R_t_nodal(i_node, 3, 2)/R_t_nodal(i_node, 3, 3));
            phi_t(i_node) = atand(R_t_nodal(i_node, 2, 1)/R_t_nodal(i_node, 1, 1));
        end
        
        eta = zeros(n_elem_tot*3, 1);
        spc = 1/double(n_node_tot-1);
        for i_elem = 2:n_elem_tot*3
            if mod(i_elem, 3) == 1
                eta(i_elem) = eta(i_elem-1);
            else
                eta(i_elem) = eta(i_elem-1)+spc;
            end
        end
        
        if node_h == -1
            node_h_elem = n_elem_tot*3;
        else
            node_h_elem = node_h*3/2-2;
        end
        h_ang(i_ang, i_pos) = -theta_t(node_h_elem);
    end
end

figure();
hold on;
for i_ang = 2:n_ang-1
    plot(pos_frac_h, h_ang(i_ang, :));
end
xlim([0.5 0.95]);
title("Wing Twist at Hinge")
xlabel("Hinge Semispan Position \eta_h");
ylabel("Twist Angle (deg)");
[leg,att] = legend(["-30 deg", "-15 deg", "0 deg", "15 deg", "30 deg"]);

title(leg,'Tip Sweep')
leg.Title.Visible = 'on';
leg.Title.NodeChildren.Position = [0.5 1.15 0];


hold off;

%% Functions
function a_tilde = skew(a)
    a_tilde = [0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0];
end

function R = psi_to_R(psi)
    ang = norm(psi);
    R = eye(3) + sin(ang)/ang*skew(psi) + (1-cos(ang))/ang^2*skew(psi)^2;
end
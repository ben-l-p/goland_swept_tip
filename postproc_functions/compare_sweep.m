close all;
load("convergence_out.mat")

t_aligned = aligned.dt:aligned.dt:aligned.dt*double(aligned.n_tstep);
t_misaligned = misaligned.dt:misaligned.dt:misaligned.dt*double(misaligned.n_tstep);

%Moment of area
ts_mom = 1;

a_mom_a = 0;
a_tot_a = 0;
[~, M_a, N_a] = size(aligned.gamma);
for i_M = 1:M_a
    for i_N = 1:N_a
        pnts = squeeze(aligned.zeta(ts_mom, :, i_M:i_M+1, i_N:i_N+1));
        midpnt = (pnts(:, 1, 1) + pnts(:, 1, 2) + pnts(:, 2, 2) + pnts(:, 2, 1))/4;

        pnts = reshape(pnts, [3, 4]);
        pnts = [pnts(:, 1), pnts(:, 2), pnts(:, 4), pnts(:, 3)];
        area = area3D(pnts(1, :), pnts(2, :), pnts(3, :));
        a_mom_a = a_mom_a + midpnt*area;
        a_tot_a = a_tot_a + area;
    end
end

a_mom_m = 0;
a_tot_m = 0;
[~, M_m, N_m] = size(misaligned.gamma);
for i_M = 1:M_m
    for i_N = 1:N_m
        pnts = squeeze(misaligned.zeta(ts_mom, :, i_M:i_M+1, i_N:i_N+1));
        midpnt = (pnts(:, 1, 1) + pnts(:, 1, 2) + pnts(:, 2, 2) + pnts(:, 2, 1))/4;

        pnts = reshape(pnts, [3, 4]);
        pnts = [pnts(:, 1), pnts(:, 2), pnts(:, 4), pnts(:, 3)];
        area = area3D(pnts(1, :), pnts(2, :), pnts(3, :));
        a_mom_m = a_mom_m + midpnt*area;
        a_tot_m = a_tot_m + area;
    end
end

a_mom_err = abs((a_mom_a - a_mom_m)./a_mom_a);

%Tip displacement
figure();
hold on;
plot(t_aligned, aligned.beam_pos(:, end, 3), '-b');
plot(t_misaligned, misaligned.beam_pos(:, end, 3), '-r');
legend(["Aligned", "Misaligned"]);
xlabel("Time (s)");
ylabel("Vertical Displacement (m)");
title("Tip Displacement");
hold off;

%Moments
mom_aligned = aligned.AeroForcesCalculator.moments_a_s + ...
    aligned.AeroForcesCalculator.moments_a_u + ...
    aligned.AeroForcesCalculator.moments_g_s + ...
    aligned.AeroForcesCalculator.moments_a_u;

mom_misaligned = misaligned.AeroForcesCalculator.moments_a_s + ...
    misaligned.AeroForcesCalculator.moments_a_u + ...
    misaligned.AeroForcesCalculator.moments_g_s + ...
    misaligned.AeroForcesCalculator.moments_a_u;

figure();
sgtitle("Total Moments");
mom_dir = ["Bending" "Torsion" "In-Plane"];

for j = 1:3
    subplot(1, 3, j);
    hold on;
    plot(t_aligned, mom_aligned(:, j), "-b");
    plot(t_misaligned, mom_misaligned(:, j), "-r");

    ylabel(mom_dir(j) + " Moment (Nm)");
    xlabel("Time (s)");
    
    if j == 3
        legend(["Aligned", "Misaligned"]);
    end
end


%Induced velocity
figure();
sgtitle("UVLM Induced Velocity");
vel_dir = ['X', 'Y', 'Z'];
% for i = 1:size(induced_v_pnts, 2)
for i = 3
    for j = 1:3
        subplot(1, 3, j);
        hold on;
        plot(t_aligned, aligned.induced_v(:, i, j), "-b");
        plot(t_misaligned, misaligned.induced_v(:, i, j), "-r");

        ylabel("Induced " + vel_dir(j) + " Velocity (m/s)");
        xlabel("Time (s)");
        
        if j == 3
            legend(["Aligned", "Misaligned"]);
        end
    end
end

%Wake gamma
figure();
title("Spanwise Gamma Distribution in Wake")
aligned_t = 50;
aligned_M = 5;
fact = misaligned.dt/aligned.dt;

hold on;
plot(squeeze(aligned.gamma_star(aligned_t, aligned_M, :)), '-b');
plot(squeeze(misaligned.gamma_star(round(aligned_t/fact), round(aligned_M/fact), :)), '-r');
xlabel("Spanwise Panel Index");
ylabel("Gamma");
legend(["Aligned", "Misaligned"]);
hold off;

function plarea = area3D(x,y,z)

%AREA3D  Area of a 3D planar polygon which does not lie in the x-y plane.
%	AREA3D(X,Y,Z) calculates the area of a polygon in space
%	formed by vertices with coordinate vectors X,Y and Z.
%	If the coordinates of vertex v_i are x_i, y_i and z_i
%	twice the area of a polygon is given by:
%	2 A(P) = abs(N . (sum_{i=0}^{n-1} (v_i x v_{i+1})))
%	where N is a unit vector normal to the plane. The `.' represents the
%	dot product operator, the `x' represents the cross product operator,
%	and abs() is the absolute value function.	

%	----------------------------------------------------------------
%  Copyright (c) 2000 by Ioan M. Buciu
%	nelu@zeus.csd.auth.gr
%	10/04/2000
%	----------------------------------------------------------------

%	Length of vectors X,Y and Z
lx = length(x);
ly = length(y);
lz = length(z);

%	Auxilliars needed for normals length
edge0 = [x(2) - x(1),y(2) - y(1),z(2) - z(1)];
edge1 = [x(3) - x(1),y(3) - y(1),z(3) - z(1)]; 

%	Cross products
nor3 = [edge0(2)*edge1(3) - edge0(3)*edge1(2),...
      edge0(3)*edge1(1) - edge0(1)*edge1(3),...
      edge0(1)*edge1(2) - edge0(2)*edge1(1)];

%	Length of normal vectors
inveln = 1/(sqrt(nor3(1)*nor3(1) + nor3(2)*nor3(2) + nor3(3)*nor3(3)));

%	Make normals unit length
csumx = zeros([1, lx-1]);
csumy = zeros([1, lx-1]);
csumz = zeros([1, lx-1]);
nor3 = inveln*nor3;
for	i = 1:lx-1
    csumx(i) = y(i)*z(i+1) - z(i)*y(i+1);
    csumy(i) = z(i)*x(i+1) - x(i)*z(i+1);
    csumz(i) = x(i)*y(i+1) - y(i)*x(i+1);
end
csumx_ = sum([csumx y(ly)*z(1)-z(lz)*y(1)]);
csumy_ = sum([csumy z(lz)*x(1)-x(lx)*z(1)]);
csumz_ = sum([csumz x(lx)*y(1)-y(ly)*x(1)]);

%	Calculate area
plarea = (abs(nor3(1)*csumx_ + nor3(2)*csumy_ + nor3(3)*csumz_))/2;
end

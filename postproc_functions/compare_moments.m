close all;
% load("case_data.mat");
[n_ang, n_pos] = size(data);

angs = zeros([n_ang, 1]);
for i_ang = 1:n_ang
    angs(i_ang) = rad2deg(data{i_ang, 1}.ang_h);
end

poss = zeros([n_pos, 1]);
for i_pos = 1:n_pos
    poss(i_pos) = data{1, i_pos}.pos_frac_h;
end

max_mom = zeros(n_ang, n_pos, 3);
for i_ang = 1:n_ang
    for i_pos = 1:n_pos
        time_moms = data{i_ang, i_pos}.AeroForcesCalculator.moments_a_s + ...
            data{i_ang, i_pos}.AeroForcesCalculator.moments_a_u + ...
            data{i_ang, i_pos}.AeroForcesCalculator.moments_g_s + ...
            data{i_ang, i_pos}.AeroForcesCalculator.moments_g_u;


        max_mom(i_ang, i_pos, :) = max(abs(time_moms), [], 1);
    end
end

[X, Y] = meshgrid(poss, angs);
titles_mom = ["Bending", "Torsion", "In-Plane"];
striaght_max_mom = zeros([3, 1]);
for i_ang = 1:n_ang
    for i_pos = 1:n_pos
        if abs(angs(i_ang)) < 1e-5
            striaght_max_mom = max_mom(i_ang, i_pos, :);
        end
    end
end
figure();
for dir = 1:2
    sgtitle("Maximum Wing Root Moments");
    subplot(1, 2, dir);
    contourf(X(2:end-1, 1:end-1), Y(2:end-1, 1:end-1), log10(max_mom(2:end-1, 1:end-1, dir)), 16);
    % shading interp;
    % view(0, 90);
    title(titles_mom(dir));
    xlabel("Hinge Curvilinear Coordinate on Beam");
    if dir == 1
        ylabel("Angle of Hinge (deg)");
    end
    c = colorbar;
    c.Label.String = "log_{10}(Moment)";
    % set(gca,'ColorScale','log')
end

figure();
t = data{1, 1}.dt:data{1, 1}.dt:data{1, 1}.dt*double(data{1, 1}.n_tstep);
hold on;
for i_ang = 1:n_ang
    for i_pos = 1:n_pos
        plot(t, data{i_ang, i_pos}.beam_pos(:, end, 3), 'k-'); 
    end
end
for i_ang = 1:n_ang
    for i_pos = 1:n_pos
        if abs(angs(i_ang)) < 1e-5
            plot(t, data{i_ang, i_pos}.beam_pos(:, end, 3), 'r-', LineWidth=5); 
        end
    end
end
title("Tip Deflection");
xlabel("Time (s)");
ylabel("Tip Vertical Deflection (m)");
xlim([min(t), max(t)]);
hold off;

%Moment of area
ts_mom = 1;

mom_a = zeros(n_ang, n_pos, 3);
tot_a = zeros(n_ang, n_pos);
[~, M_a, N_a] = size(data{1, 1}.gamma);
for i_ang = 1:n_ang
    for i_pos = 1:n_pos
        for i_M = 1:M_a
            for i_N = 1:N_a
                pnts = squeeze(data{i_ang, i_pos}.zeta(ts_mom, :, i_M:i_M+1, i_N:i_N+1));
                midpnt = (pnts(:, 1, 1) + pnts(:, 1, 2) + pnts(:, 2, 2) + pnts(:, 2, 1))/4;
                midpnt = midpnt([2 1 3]);
        
                pnts = reshape(pnts, [3, 4]);
                pnts = [pnts(:, 1), pnts(:, 2), pnts(:, 4), pnts(:, 3)];
                area = area3D(pnts(1, :), pnts(2, :), pnts(3, :));
                mom_a(i_ang, i_pos, :) = squeeze(mom_a(i_ang, i_pos, :)) + midpnt*area;
                tot_a(i_ang, i_pos) = tot_a(i_ang, i_pos) + area;
            end
        end
    end
end

figure();
subplot(2, 2, 1);
contourf(X, Y, tot_a);
view(0, 90);
shading interp;
c = colorbar;
c.Label.String = "Total Area (m^2)";
xlabel("Fraction of Span to Hinge");
ylabel("Angle of Hinge (deg)");
title("Total Area");

titles_area = ["x Moment", "y Moment", "z Moment"];
for dir = 1:3
    subplot(2, 2, dir+1);
    contourf(X, Y, mom_a(:, :, dir));
    view(0, 90);
    shading interp;
    c = colorbar;
    c.Label.String = "First Moment of Area (m^3)";
    xlabel("Fraction of Span to Hinge");
    ylabel("Angle of Hinge (deg)");
    title(titles_area(dir));
end

% figure();
% hold on;
% for dir = 1:3
%     plot(angs, max_mom(:, 4, dir));
% end
% xlabel("Hinge Angle (deg)");
% ylabel("Maximum Moment (Nm)");
% legend(["Out-of-Plane", "Torsion", "In-Plane"]);
% title("Maximum Moments with varying Hinge Angle");
% yscale log;
% hold off;
% 
% figure();
% hold on;
% for dir = 1:3
%     plot(poss, max_mom(7, :, dir));
% end
% xlabel("Hinge Placement (Fraction of Semispan)");
% ylabel("Maximum Moment (Nm)");
% legend(["Out-of-Plane", "Torsion", "In-Plane"]);
% title("Maximum Moments with varying Hinge Placement");
% yscale log;
% hold off;

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

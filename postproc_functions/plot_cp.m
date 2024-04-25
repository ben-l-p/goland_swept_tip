[n_ang, n_pos] = size(data);
n_tstep = data{1, 1}.n_tstep;

m_t = zeros(n_ang, n_pos, n_tstep, 3);
f_t = zeros(n_ang, n_pos, n_tstep, 3);
r_cp = zeros(n_ang, n_pos, n_tstep, 3);
m_err = zeros(n_ang, n_pos, n_tstep);

for i_ang = 1:n_ang
    for i_pos = 1:n_pos
        m_t(i_ang, i_pos, :, :) = data{i_ang, i_pos}.AeroForcesCalculator.moments_a_s + ...
            data{i_ang, i_pos}.AeroForcesCalculator.moments_a_u + ...
            data{i_ang, i_pos}.AeroForcesCalculator.moments_g_s + ...
            data{i_ang, i_pos}.AeroForcesCalculator.moments_g_u;

        f_t(i_ang, i_pos, :, :) = data{i_ang, i_pos}.AeroForcesCalculator.forces_a_s + ...
            [data{i_ang, i_pos}.AeroForcesCalculator.forces_a_u 0] + ...
            data{i_ang, i_pos}.AeroForcesCalculator.forces_g_s + ...
            data{i_ang, i_pos}.AeroForcesCalculator.forces_g_u;
        
        for i_t = 1:n_tstep
            f = squeeze(-f_t(i_ang, i_pos, i_t, :));
            m = squeeze(m_t(i_ang, i_pos, i_t, :));
            
            r_cp(i_ang, i_pos, i_t, 3) = 0;
            r_cp(i_ang, i_pos, i_t, 1) = ...
                -(r_cp(i_ang, i_pos, i_t, 3)*f(1)-m(2))/f(3);
            r_cp(i_ang, i_pos, i_t, 2) = ...
                -(r_cp(i_ang, i_pos, i_t, 3)*f(2)+m(1))/f(3);

            m_err(i_ang, i_pos, i_t) = norm(skew(r)*f+m)/norm(m);
        end
    end
end

function a_skew = skew(a)
    a_skew = [0 -a(3) a(2); a(3), 0, -a(1); -a(2), a(1), 0];
end

figure();
hold on;
plot(aligned.dt:aligned.dt:aligned.dt*double(aligned.n_tstep), ...
    aligned.beam_pos(:, end, 3), LineWidth=1.5);
plot(misaligned.dt:misaligned.dt:misaligned.dt*double(misaligned.n_tstep), ...
    misaligned.beam_pos(:, end, 3), LineWidth=1.5);
legend(["Aligned", "Misaligned"]);
% title("Wing Response to '1-cos' Gust");
xlabel("Time (s)");
ylabel("Tip Vertical Displacement (m)");
hold off;
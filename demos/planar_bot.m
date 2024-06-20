bot = PlanarBot;

kin = bot.kin;

%% Define a few paths
xyz_start = [-1; 1;0];
xyz_outside = [0; 2.5;0];
xyz_internal = [0;0;0];
lambda = linspace(0,1,100);

XYZ_1 = lambda.*xyz_outside + (1-lambda).*xyz_start;
XYZ_2 = lambda.*xyz_internal + (1-lambda).*xyz_outside;
XYZ_3 = lambda.*xyz_start + (1-lambda).*xyz_internal;
XYZ_123 = [XYZ_1 XYZ_2 XYZ_3];

diagrams.setup();
view(2);
hold on;
bot.plot([0;0]);

diagrams.dot(xyz_start, color=diagrams.colors.red);
diagrams.dot(xyz_outside, color=diagrams.colors.red);
diagrams.dot(xyz_internal, color=diagrams.colors.red);

hold off
diagrams.redraw();


%%
Q_1 = NaN(2, length(XYZ_123));
Q_2 = NaN(2, length(XYZ_123));
Q_single = NaN(2, length(XYZ_123));
norm_x_q2_vec = NaN(1,length(XYZ_123));
internal_1_mat = NaN(2, length(XYZ_123));
internal_2_mat = NaN(2, length(XYZ_123));

for i = 1:length(XYZ_123)
    T = XYZ_123(:,i);
    [Q_i, is_LS_vec, diag_sp_3, internal_1_vec, internal_2_vec]  = bot.IK(T);
    norm_x_q2_vec(i) = diag_sp_3.norm_x;
    internal_1_mat(1:length(internal_1_vec), i) = internal_1_vec;
    internal_2_mat(1:length(internal_2_vec), i) = internal_2_vec;

    if width(Q_i) == 2
        Q_1(:,i) = Q_i(:,1);
        Q_2(:,i) = Q_i(:,2);
    else
        Q_single(:,i) = Q_i;
    end
end

%%
for i = 1:1:length(XYZ_123)
    diagrams.setup();
    view(2);
    hold on;
    if ~isnan(Q_1(:,i))
        bot.plot(Q_1(:,i));
        bot.plot(Q_2(:,i));
    end
    if ~isnan(Q_single(:,i))
        bot.plot(Q_single(:,i));
    end
    diagrams.utils.plot3_mat(XYZ_123, lineStyle=":");
    diagrams.utils.plot3_mat(XYZ_123(:,i), color=diagrams.colors.red, marker="x");
    diagrams.utils.plot3_mat([-1; 2.6; 0], marker=".");
    hold off
    diagrams.redraw();
    % exportgraphics(gca,"planar_demo.gif","Append",i ~= 1)
end
%%
Q_1_disp = Q_1;
Q_1_disp(1,200:end) = unwrap(Q_1_disp(1,200:end)) + 2*pi;

plot(Q_1_disp(1,:), 'r'); hold on
plot(Q_1_disp(2,:), 'r-.');
plot(Q_2(1,:), 'g');
plot(Q_2(2,:), 'g-.');
plot(Q_single(1,:), 'k'); 
plot(Q_single(2,:), 'k-.'); hold off
xline([200, 201]);
legend(["q_1^1", "q_1^2", "q_2^1", "q_2^2", "q_2^{1,2}", "q_2^{1,2}"], Orientation='Horizontal', Location='south')
ylim([-1.5*pi, 1.5*pi]);
yticks([-1.5*pi, -pi, -pi/2, 0, pi/2, pi, 1.5*pi]);
yticklabels(["-3\pi/2", "-\pi", "-\pi/2", "0", "\pi/2", "\pi", "3\pi/2"])
xlabel("Time (sec)")
ylabel("Joint Angle (rad)")
title("Joint Angles")

%%
plot(norm_x_q2_vec, 'k');
yline(1);
xlabel("Time (sec)")
ylabel("$\left\Vert  x_{min} \right\Vert$", Interpreter="latex")
title("Subproblem 3 Singularity Condition")

%%
plot(internal_1_mat(1,:), 'r'); hold on
plot(internal_2_mat(1,:)', 'g');
hold off

xlabel("Time (sec)")
ylabel("$\left\Vert  p ^\times k_i \right\Vert$", Interpreter="latex")
title("Subproblem 1 Singularity Conditions")
legend(["$i=1$", "$i=2$"], Interpreter="latex")
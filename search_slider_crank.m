kin = define_slider_crank();
%%
q1_vec = linspace(-pi,pi, 1e3);
ex = [1;0;0];
ey = [0;1;0];
ez = [0;0;1];
zv = [0;0;0];

p_0T = [1; 0.75; 0];

e_vec = NaN([3 length(q1_vec)]);
for i = 1:length(q1_vec)
    q1 = q1_vec(i);
    R_01 = rot(kin.H(:,1), q1);
    p_02 = kin.P(:,1)+R_01*kin.P(:,2);
    [t02, t02_is_ls] = subproblem.sp_4(ey, kin.P(:,3), kin.H(:,2), ey'*(p_0T - p_02));
    for i_t02 = 1:length(t02)
        R_02 = rot(kin.H(:,1), t02(i_t02));
        p_0T_i = p_02 + R_02 * kin.P(:,2);
        e_i = ex'*(p_0T_i-p_0T);
        if t02_is_ls
            e_vec(3, i) = e_i;
        else
            e_vec(i_t02, i) = e_i;
        end
    end
end
plot(q1_vec, e_vec);
yline(0);
xline(1.065)

%%
q1 = 1.065;
R_01 = rot(kin.H(:,1), q1);
p_02 = kin.P(:,1)+R_01*kin.P(:,2);
[t02, t02_is_ls] = subproblem.sp_4(ey, kin.P(:,3), kin.H(:,2), ey'*(p_0T - p_02));
q_A = [q1 t02(1)-q1];
q_B = [q1 t02(2)-q1];

diagrams.setup;
view(2)
hold on
options = {"show_arrows", true,...
            "unit_size", 0.25,...
           "show_arrow_labels", false,...
           "show_joint_labels", false,...
           "show_base_label", false,...
           "show_base_frame", false,...
           "show_task_frame", false,...
           "cyl_half_length", 0.125,...
           "cyl_radius", 0.05};

diagrams.robot_plot(kin, q_A, options{:});
% diagrams.robot_plot(kin, q_B, options{:});

hold off
yline(p_0T(2));
diagrams.redraw;

function kin = define_slider_crank()
ex = [1;0;0];
ey = [0;1;0];
ez = [0;0;1];
zv = [0;0;0];

kin.H = [ez ez];
kin.P = [zv 2*ex ex];
kin.joint_type = [0 0];
end
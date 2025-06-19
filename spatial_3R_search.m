% kin = define_3R_bot();
% kin = define_3R_bot_rand();
kin = define_3R_bot_angled();

%%
q_vec = linspace(-pi/10,pi/10, 10);
p_mat = [];
figure; hold on; view(3);
for q1 = Q(1,1) + q_vec
    for q2 = Q(2,1) + q_vec
        for q3 = Q(3,1) + q_vec
            [~, p_0T_i] = fwdkin(kin, [q1; q2; q3]);
            p_mat = [p_mat p_0T_i];
        end
    end
end
plot3(p_mat(1,:), p_mat(2,:), p_mat(3,:), 'kx')
plot3(p_0T(1), p_0T(2), p_0T(3), 'ro')

%%
diagrams.setup;
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

diagrams.robot_plot(kin, Q(:,1), options{:});

hold off
diagrams.redraw;

%%
% q = rand_angle([3 1]);
% q = [pi/4;pi/8;0];

%[~, p_0T] = fwdkin(kin, q);
 p_0T = [0.8;0;-1];
 % p_0T = [1.6;-1;-1.5];

Q = IK_3R_bot(kin, p_0T)


plot_error_q1(kin, p_0T)
xline(Q(1,:))

J = robotjacobian(kin, Q(:,1));
J_p = J(4:end,:);
svd_1 = svd(J_p)

J = robotjacobian(kin, Q(:,2));
J_p = J(4:end,:);
svd_2 = svd(J_p)

function plot_error_q1(kin, p_0T)
q1_vec = linspace(-pi,pi, 1e3);
e_vec = NaN([3 length(q1_vec)]);

for i = 1:length(q1_vec)
    q1 = q1_vec(i);
    R_01 = rot(kin.H(:,1), q1);
    p_02 = kin.P(:,1) + R_01*kin.P(:,2);
    [t3, t3_is_LS] = subproblem.sp_3(kin.P(:,4), -kin.P(:,3), kin.H(:,3), norm(p_0T-p_02));

    for i_t3 = 1:length(t3)
        R_23 = rot(kin.H(:,3), t3(i_t3));
        e_i = dot(kin.H(:,2), R_01'*(p_0T - p_02) - kin.P(:,3) - R_23*kin.P(:,4));
        if t3_is_LS
            e_vec(3,i) = e_i;
        else
            e_vec(i_t3,i) = e_i;
        end
    end
end

plot(q1_vec,e_vec');
yline(0);
xlabel("q_1")
ylabel("e(q_1)")

end

function Q = IK_3R_bot(kin, p_0T)
    [Q1, Q2, Q3] = subproblem.sp_5(-kin.P(:,2), p_0T - kin.P(:,1), kin.P(:,3), kin.P(:,4), -kin.H(:,1), kin.H(:,2), kin.H(:,3));
    Q = [Q1;Q2;Q3];
end

function kin = define_3R_bot_rand()
zv = [0;0;0];
kin.H = rand_normal_vec(3);
kin.P = [zv rand_vec(3)];
kin.joint_type = [0 0 0];
end

function kin = define_3R_bot_angled()
ex = [1;0;0];
ey = [0;1;0];
ez = [0;0;1];
zv = [0;0;0];

kin.H = [ez ey ez];
kin.P = [zv ex ex+ey ex+ez];
kin.joint_type = [0 0 0];
end

function kin = define_3R_bot()
ex = [1;0;0];
ey = [0;1;0];
ez = [0;0;1];
zv = [0;0;0];

kin.H = [ez ey ez];
kin.P = [zv ex ex ex];
kin.joint_type = [0 0 0];
end
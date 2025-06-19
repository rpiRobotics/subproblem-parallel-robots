kin = hardcoded_IK_setups.IRB_6640.get_kin();

q_home = zeros(6,1);

[R, T] = fwdkin(kin, q_home)


J = robotjacobian(kin, q_home);
[U, S, V] = svd(J)

[Q, is_LS_vec, is_LS_mat, diag_mat] = robot.IK_comp(R,T)
%% Test comprehensive IK
robot = IRB;

q_r = rand_angle([6,1]);
[R, T] = fwdkin(robot.kin, q_r);
[Q, is_LS_vec, is_LS_mat, diag_mat] = robot.IK_comp(R,T)

%% Make a path that enters a singularity

p_A = [1; 0; 2.0505];
p_B = [1.6625; 0 ; 2.0550];
R = eye(3);

lambda = linspace(0,1);
q_path_enter = NaN(6, length(lambda));
diag_path_enter = NaN([1,length(lambda)]);
for i = 1:length(lambda)
    p_i = lambda(i)*p_B + (1-lambda(i))*p_A;
    [Q_i, ~, ~, diag_mat_i] =  robot.IK_comp(R, p_i);
    q_path_enter(:,i) = Q_i(:,3);
    diag_path_enter(i) = diag_mat_i{4,3}.norm_k_X_p1;
end

plot(q_path_enter');
%%
semilogy(diag_path_enter);

%% Make the second half of the path to leave the singularity
p_C = [1; 0; 2.0505];

q_path_exit = NaN(6, length(lambda));
diag_path_exit = NaN([1,length(lambda)]);
N_soln = 4;
for i = 1:length(lambda)
    p_i = lambda(i)*p_C + (1-lambda(i))*p_B;
    [Q_i, ~, ~, diag_mat_i] =  robot.IK_comp(R, p_i);
    q_path_exit(:,i) = Q_i(:,N_soln);
    diag_path_exit(i) = diag_mat_i{4,N_soln}.norm_k_X_p1;
end
%%
plot([q_path_enter q_path_exit]');
legend(string(1:6));
%%
plot([q_path_enter(:,1:end-3) q_path_exit(:,4:end)]');
legend(string(1:6));
%% Create the in-singularity path

% Project the enter and exit poses onto the feasible subspace
% Fesible subspace is [x;y;z; q4; q5; theta_46]

%i_enter = find(diag_path_enter < 1e-4, 1);
%q_enter = q_path_enter(:,i_enter)
q_enter = q_path_enter(:,end-3);

%i_exit = find(diag_path_exit < 1e-4, 1, "last");
%q_exit = q_path_exit(:,i_exit)
q_exit = q_path_exit(:,4);

% Linearly interpolate
q4_singularity = lambda*q_exit(4) + (1-lambda)*q_enter(4);

% Do IK from the in-singularity task space
N_soln = 3;
q_path_singularity = NaN([6,length(lambda)]);
for i = 1:length(lambda)
    Q_i =  robot.IK_comp(R, p_B, q4_singularity(i));
    q_path_singularity(:,i) = Q_i(:,N_soln);
end
plot(q_path_singularity')

%%
q_path_total = [q_path_enter(:,1:end-5) q_path_singularity q_path_exit(:,5:end)];
plot(q_path_total');
legend(string(1:6));
%%
diagrams.setup;
    campos([-7 -10 8]);
    camva(9);
    camtarget([0.65 0 1]);
hold on
options = {"show_arrows", true,...
            "unit_size", 0.25,...
           "show_arrow_labels", false,...
           "show_joint_labels", false,...
           "show_base_label", false,...
           "show_base_frame", false,...
           "show_task_frame", true,...
           "cyl_half_length", 0.125,...
           "cyl_radius", 0.05};

diagrams.robot_plot(kin, q_home, options{:});

hold off
diagrams.redraw;

%%
for i = 1:5:length(q_path_total)
    diagrams.setup();
    annotation('rectangle',[0 0 1 1]);
    campos([-7 -10 8]);
    camva(9);
    camtarget([0.65 0 1]);
    hold on;
    diagrams.robot_plot(kin, q_path_total(:,i), options{:});
    hold off
    diagrams.redraw();
    drawnow
    exportgraphics(gca,"irb_wrist_singularity_resolution.gif","Append",i ~= 1)
end
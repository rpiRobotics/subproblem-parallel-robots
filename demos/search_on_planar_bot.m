bot = PlanarBot;
kin = bot.kin;

p_0T = [0.25; 0.25; 0];
t = "Case 1: Nonsingular";

% p_0T = [0; 2; 0];
% t = "Case 2: Boundary Singularity"

% p_0T = [0; 0; 0];
% t = "Case 3: Internal Singularity";

Q = bot.IK(p_0T);
% [R, T] = fwdkin(kin, Q(:,1))

l = tiledlayout(1,2, TileSpacing="tight", Padding="tight");
title(l,t);

nexttile
plot_error_q1(kin, p_0T)
xline(Q(1,:))
nexttile

plot_error_q2(kin, p_0T)
xline(Q(2,:))   
fig=gcf;
fig.Position(3:4)=[800,400];

function plot_error_q1(kin, p_0T)
q1_vec = linspace(-pi,pi);
e_vec = NaN(size(q1_vec));

for i = 1:length(q1_vec)
    q1 = q1_vec(i);
    p_01 = rot(kin.H(:,1), q1)*kin.P(:,2);
    e_i = norm(p_01 - p_0T) - norm(kin.P(:,2));
    e_vec(i) = e_i;
end

plot(q1_vec,e_vec);
yline(0);
xlabel("q_1")
ylabel("e(q_1)")
end

function plot_error_q2(kin, p_0T)
q2_vec = linspace(-pi,pi);
e_vec = NaN(size(q2_vec));

for i = 1:length(q2_vec)
    q2 = q2_vec(i);
    e_i = norm(kin.P(:,2) + rot(kin.H(:,2), q2)*kin.P(:,3))^2 - norm(p_0T)^2;
    e_vec(i) = e_i;
end

plot(q2_vec,e_vec);
yline(0);
xlabel("q_2")
ylabel("e(q_2)")
end
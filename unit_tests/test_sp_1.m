zv = [0;0;0];

% Define a circle
k = [0 0 1]';
p = [1.2 0 0]';

% Make a path through the circle
xyz_start = [0.5 -2 0]';
xyz_end   = [0.5  2 0]';
lambda = linspace(0,1,100);
xyz_path = (1-lambda).*xyz_start + lambda.*xyz_end;

% Plot the circle and path
diagrams.setup();
hold on
diagrams.cone(zv, p, k);
diagrams.utils.plot3_mat(xyz_path);
hold off
diagrams.redraw();

%%
N = length(xyz_path);
theta_vec = NaN(1,N);
is_LS_vec = NaN(1,N);
norm_x_vec = NaN(1,N);

for i = 1:N
    [theta_vec(i), is_LS_vec(i), norm_x_vec(i)] = subproblem_comp.sp_1(p, xyz_path(:,i), k);
end
%%
plot(theta_vec); hold on
plot(is_LS_vec);
plot(norm_x_vec);
hold off

%% Test on exact soln
[theta, is_LS, norm_x] = subproblem_comp.sp_1(p, rot(k,pi/10)*p, k)
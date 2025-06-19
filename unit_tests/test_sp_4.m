zv = [0;0;0];

% Define a circle
k = [0 0 1]';
p = [1.2 0 0]';

% Make a path through the circle
xyz_start = [0.5 -2 0]';
xyz_end   = [0.5  2 0]';
lambda = linspace(0,1,100);
xyz_path = (1-lambda).*xyz_start + lambda.*xyz_end;

h = [0;1;0];

% Plot the circle and path
diagrams.setup();
hold on
diagrams.cone(zv, p, k);
diagrams.utils.plot3_mat(xyz_path, linestyle=":");
diagrams.arrow(xyz_start, xyz_start+h);
diagrams.plane(xyz_start, h, cross(xyz_start, h), 1)
hold off
diagrams.redraw();
%%
N = length(xyz_path);
theta_vec = NaN(2,N);
is_LS_vec = NaN(1,N);
norm_x_vec = NaN(1,N);

for i = 1:N
    [theta_vec(:,i), is_LS_vec(i), norm_x_vec(i)] = subproblem_comp.sp_4(h, p, k, h'*xyz_path(:,i));
end
%%
plot(theta_vec'); hold on
plot(is_LS_vec);
plot(norm_x_vec);
hold off
%% Plot intermediate result
i = 40;

diagrams.setup();
hold on
diagrams.cone(zv, p, k);
diagrams.arrow(xyz_path(:,i), xyz_path(:,i)+h);
diagrams.plane(xyz_path(:,i), h, cross(xyz_path(:,i), h), 2)
diagrams.dot(rot(k,theta_vec(1,i))*p);
diagrams.dot(rot(k,theta_vec(2,i))*p);
hold off
diagrams.redraw();
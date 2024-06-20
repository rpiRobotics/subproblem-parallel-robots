function plot(obj, q)
kin = obj.kin;
ex = [1;0;0];
ey = [0;1;0];
ez = [0;0;1];
zv = [0;0;0];

% Workspace boundary circle
diagrams.circle(zv, ez, sum(vecnorm(kin.P)), LineStyle='--');

[~, T, P_inter] = fwdkin_inter(kin, q, 2);
diagrams.line(zv, P_inter(:,1));
diagrams.line(P_inter, T);
diagrams.dot(zv);
diagrams.dot(P_inter);
diagrams.dot(T);
end
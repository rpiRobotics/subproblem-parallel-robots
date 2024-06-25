rob = Eclipse;

%%
codegen -report mex_funcs/eclipse_FK.m -args {rob, q_A, q_B, q_C}
%%
q_A = [0; 0]; % Fixed w.l.o.g.
q_B = [deg2rad(120); 0.2]; % Fixed for this graph

N1 = 100;
N2 = 100;
q1_C_vec = linspace(-pi, pi, N1);
q2_C_vec = linspace(-1.5, 1.5, N2);
min_mat = NaN(N1, N2, 16);

for i_1 = 1:N1
    i_1
    for i_2 = 1:N2
    q_C = [q1_C_vec(i_1); q2_C_vec(i_2)];
    [R_0T_vec, p_0T_vec, q_A_vec, q_B_vec, q_C_vec, slopes] = eclipse_FK_mex(rob,q_A, q_B, q_C);
    % if ~isempty(slopes)
    %     min_i = slopes(1);
    %     %min_i = min(abs(slopes));
    %     min_mat(i_1, i_2) = min_i;    
    % end
    min_mat(i_1, i_2, 1:length(slopes)) = slopes;
    end
end
%%
[X, Y] = meshgrid(q1_C_vec, q2_C_vec);
surf(X, Y, abs(min_mat(:,:,1)), EdgeColor="none")
view(2)
colorbar
%%
[X, Y] = meshgrid(q1_C_vec, q2_C_vec);
surf(X, Y, min_mat)
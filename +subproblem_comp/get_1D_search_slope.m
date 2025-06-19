function slopes = get_1D_search_slope(fun, x_vec, soln_num_vec)

DELTA = 1e-5;
slopes = NaN(size(x_vec));

for i = 1:length(x_vec)
    x_i = x_vec(i);
    soln_num_i = soln_num_vec(i);
    v1 = fun(x_i - DELTA);
    v2 = fun(x_i + DELTA);
    y1 = v1(soln_num_i);
    y2 = v2(soln_num_i);
    slopes(i) = (y2 - y1)/ (2*DELTA);
end
end
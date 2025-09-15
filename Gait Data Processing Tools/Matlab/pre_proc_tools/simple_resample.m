% Resample tabular data to different frequency
function data = simple_resample(data, f_orig, f_new)
    dt_orig = 1/f_orig;
    n_orig = size(data, 1);
    T_orig = dt_orig*(n_orig-1);
    t_orig = linspace(0, T_orig, n_orig);

    n_down = int32((f_new / f_orig) * n_orig);
    t_down = zeros(n_down, 1);
    for i = 1:n_down
        t_down(i) = t_orig(1+int32((i-1)*(f_orig/f_new)));
    end

    n_col = size(data, 2);
    [X1, Y1] = meshgrid(1:n_col, t_orig);
    [X2, Y2] = meshgrid(1:n_col, t_down);
    data = interp2(X1, Y1, data, X2, Y2, "nearest");
end
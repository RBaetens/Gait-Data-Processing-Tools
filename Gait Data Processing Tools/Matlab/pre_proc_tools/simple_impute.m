% Impute isolated points of missing data through linear interpolation
function data = simple_impute(data)
    n_pts = size(data, 1);
    n_col = size(data, 2);
    for i = 2:n_pts-1
        for j = 1:n_col
            if isnan(data(i, j)) && not(isnan(data(i-1, j)) || isnan(data(i+1, j)))
                data(i, j) = (data(i-1, j) + data(i+1, j)) / 2;
            end
        end
    end
end
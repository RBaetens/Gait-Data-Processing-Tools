function data = shift_impute(data)
    % Impute isolated points of missing data through linear interpolation
    n_pts = size(data, 1);
    for i = 2:n_pts-1
        if isnan(data(i)) && not(isnan(data(i-1)) || isnan(data(i+1)))
            data(i) = round((data(i-1) + data(i+1)) / 2);
        end
    end

    % Fist step doesn't catch evt --> change NaN values by one of the neighbouring shifts
    % Forward sweep
    for i = 2:n_pts
        if isnan(data(i)) && ~isnan(data(i-1))
            data(i) = data(i-1);
        end
    end
    % Backward sweep
    for i = 1:n_pts-1
        if isnan(data(n_pts-i)) && ~isnan(data(n_pts-i+1))
            data(n_pts-i) = data(n_pts-i+1);
        end
    end
end
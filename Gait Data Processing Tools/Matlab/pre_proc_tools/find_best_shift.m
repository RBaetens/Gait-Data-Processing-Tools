function best_shift = find_best_shift(Fv_FP1, Fv_FP2, Fv_FP3, Fv_FP4, ...
    Fv_INL, Fv_INR, ...
    matFP, matAllTraj, ThLOW, max_shift)

    % Foot on force plate
    foot_on_fp_cols = foot_on_fp(-[matFP(:, 3), matFP(:, 10), matFP(:, 17), matFP(:, 24)], ... % Fz data
        [matFP(:, 4:6), matFP(:, 11:13), matFP(:, 18:20), matFP(:, 25:27)], ... % CoP data
        matAllTraj(:, 25:27), matAllTraj(:, 28:30), matAllTraj(:, 43:45), matAllTraj(:, 46:48), ... % LHEE, LTOE, RHEE, RTOE
        ThLOW);
    
    % Make a wider version of 'foot_on_fp_cols': this expands the region where
    % we expect to be able to meaningfully look for correlation
    foot_on_fp_cols_wide = widen_blocks(foot_on_fp_cols, max_shift);
    
    % Find globally best shift
    n_shifts = 2*max_shift+1;
    shifts = linspace(-max_shift, max_shift, n_shifts);
    measures = zeros(n_shifts, 1);
    
    if (sum(Fv_FP1 > ThLOW) == 0) && (sum(Fv_FP2 > ThLOW) == 0) && (sum(Fv_FP3 > ThLOW) == 0) && (sum(Fv_FP4 > ThLOW) == 0)
        best_shift = NaN;
    else
        best_shift = NaN;
        best_measure = 0;
        for i=1:n_shifts
            measures(i) = nansum(Fv_INL(1+max_shift+shifts(i):end-max_shift+shifts(i)).*Fv_FP1(max_shift+1:end-max_shift).*foot_on_fp_cols_wide(max_shift+1:end-max_shift, 1));
            measures(i) = measures(i) + nansum(Fv_INR(1+max_shift+shifts(i):end-max_shift+shifts(i)).*Fv_FP1(max_shift+1:end-max_shift).*foot_on_fp_cols_wide(max_shift+1:end-max_shift, 2));
            measures(i) = measures(i) + nansum(Fv_INL(1+max_shift+shifts(i):end-max_shift+shifts(i)).*Fv_FP2(max_shift+1:end-max_shift).*foot_on_fp_cols_wide(max_shift+1:end-max_shift, 3));
            measures(i) = measures(i) + nansum(Fv_INR(1+max_shift+shifts(i):end-max_shift+shifts(i)).*Fv_FP2(max_shift+1:end-max_shift).*foot_on_fp_cols_wide(max_shift+1:end-max_shift, 4));
            measures(i) = measures(i) + nansum(Fv_INL(1+max_shift+shifts(i):end-max_shift+shifts(i)).*Fv_FP3(max_shift+1:end-max_shift).*foot_on_fp_cols_wide(max_shift+1:end-max_shift, 5));
            measures(i) = measures(i) + nansum(Fv_INR(1+max_shift+shifts(i):end-max_shift+shifts(i)).*Fv_FP3(max_shift+1:end-max_shift).*foot_on_fp_cols_wide(max_shift+1:end-max_shift, 6));
            measures(i) = measures(i) + nansum(Fv_INL(1+max_shift+shifts(i):end-max_shift+shifts(i)).*Fv_FP4(max_shift+1:end-max_shift).*foot_on_fp_cols_wide(max_shift+1:end-max_shift, 7));
            measures(i) = measures(i) + nansum(Fv_INR(1+max_shift+shifts(i):end-max_shift+shifts(i)).*Fv_FP4(max_shift+1:end-max_shift).*foot_on_fp_cols_wide(max_shift+1:end-max_shift, 8));
        
            if measures(i) > best_measure
                best_shift = shifts(i);
                best_measure = measures(i);
            end
        end
    end
end
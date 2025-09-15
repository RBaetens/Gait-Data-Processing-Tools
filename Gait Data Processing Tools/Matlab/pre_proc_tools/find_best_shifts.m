function [best_shift_L, best_shift_R] = find_best_shifts(Fv_FP1, Fv_FP2, Fv_FP3, Fv_FP4, ...
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
    n_shifts = 2*max_shift+1;
    shifts = linspace(-max_shift, max_shift, n_shifts);
    
    % Find best shift for the left insole
    measures_L = zeros(n_shifts, 1);
    
    if sum(foot_on_fp_cols(:, 1:2:end)) == 0
        best_shift_L = NaN;

    else
        best_shift_L = NaN;
        best_measure_L = 0;

        for i=1:n_shifts
            measures_L(i) = nansum(Fv_INL(1+max_shift+shifts(i):end-max_shift+shifts(i)).*Fv_FP1(max_shift+1:end-max_shift).*foot_on_fp_cols_wide(max_shift+1:end-max_shift, 1));
            measures_L(i) = measures_L(i) + nansum(Fv_INL(1+max_shift+shifts(i):end-max_shift+shifts(i)).*Fv_FP2(max_shift+1:end-max_shift).*foot_on_fp_cols_wide(max_shift+1:end-max_shift, 3));
            measures_L(i) = measures_L(i) + nansum(Fv_INL(1+max_shift+shifts(i):end-max_shift+shifts(i)).*Fv_FP3(max_shift+1:end-max_shift).*foot_on_fp_cols_wide(max_shift+1:end-max_shift, 5));
            measures_L(i) = measures_L(i) + nansum(Fv_INL(1+max_shift+shifts(i):end-max_shift+shifts(i)).*Fv_FP4(max_shift+1:end-max_shift).*foot_on_fp_cols_wide(max_shift+1:end-max_shift, 7));
            
            if measures_L(i) > best_measure_L
                best_shift_L = shifts(i);
                best_measure_L = measures_L(i);
            end
        end
    end

    % Find best shift for the right insole
    measures_R = zeros(n_shifts, 1);
    
    if sum(foot_on_fp_cols(:, 2:2:end)) == 0
        best_shift_R = NaN;
        
    else
        best_shift_R = NaN;
        best_measure_R = 0;

        for i=1:n_shifts
            measures_R(i) = nansum(Fv_INR(1+max_shift+shifts(i):end-max_shift+shifts(i)).*Fv_FP1(max_shift+1:end-max_shift).*foot_on_fp_cols_wide(max_shift+1:end-max_shift, 2));
            measures_R(i) = measures_R(i) + nansum(Fv_INR(1+max_shift+shifts(i):end-max_shift+shifts(i)).*Fv_FP2(max_shift+1:end-max_shift).*foot_on_fp_cols_wide(max_shift+1:end-max_shift, 4));
            measures_R(i) = measures_R(i) + nansum(Fv_INR(1+max_shift+shifts(i):end-max_shift+shifts(i)).*Fv_FP3(max_shift+1:end-max_shift).*foot_on_fp_cols_wide(max_shift+1:end-max_shift, 6));
            measures_R(i) = measures_R(i) + nansum(Fv_INR(1+max_shift+shifts(i):end-max_shift+shifts(i)).*Fv_FP4(max_shift+1:end-max_shift).*foot_on_fp_cols_wide(max_shift+1:end-max_shift, 8));

            if measures_R(i) > best_measure_R
                best_shift_R = shifts(i);
                best_measure_R = measures_R(i);
            end
        end
    end

    % See if feet were on the same force plate and put best_shift_L =
    % best_shift_R = NaN in that case

    % Step 1: find number of cycles on force plate
    fp_active = zeros(size(foot_on_fp_cols_wide, 1), 4);
    fp_active(:, 1) = foot_on_fp_cols_wide(:, 1) | foot_on_fp_cols_wide(:, 2);
    fp_active(:, 2) = foot_on_fp_cols_wide(:, 3) | foot_on_fp_cols_wide(:, 4);
    fp_active(:, 3) = foot_on_fp_cols_wide(:, 5) | foot_on_fp_cols_wide(:, 6);
    fp_active(:, 4) = foot_on_fp_cols_wide(:, 7) | foot_on_fp_cols_wide(:, 8);

    n_cycles = size(block_wave_edges(fp_active(:, 1)), 1);
    n_cycles = n_cycles + size(block_wave_edges(fp_active(:, 2)), 1);
    n_cycles = n_cycles + size(block_wave_edges(fp_active(:, 3)), 1);
    n_cycles = n_cycles + size(block_wave_edges(fp_active(:, 4)), 1);

    % Step 2: find number of times two feet were on the same force plate
    % Not exact, n_two_feet could be greater than actual number
    n_two_feet = 0;
    oog = eye(2);
    ander_oog = [[0, 1]; [1, 0]];
    for i=1:4
        for j=1:(size(foot_on_fp_cols, 1)-1)
            if all(oog == foot_on_fp_cols(j:j+1, 2*i-1:2*i), 'all') || all(ander_oog == foot_on_fp_cols(j:j+1, 2*i-1:2*i), 'all')
                n_two_feet = n_two_feet + 1;
            end
        end
    end

    % Check and correct
    if n_two_feet >= n_cycles
        best_shift_L = NaN;
        best_shift_R = NaN;
    end
end
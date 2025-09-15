% check validity: force plate does/does not give bull shit values
function bad_fp_mask = find_bad_fp(matFP, foot_on_fp_cols, pos_FP1, pos_FP2, pos_FP3, pos_FP4)
    n_down = size(matFP, 1);

    % Define corners of the force plates with extra safety margin
    corners_FP1 = [[pos_FP1(1)-190, pos_FP1(2)-240]; [pos_FP1(1)+190, pos_FP1(2)-240]; [pos_FP1(1)-190, pos_FP1(2)+240]; [pos_FP1(1)+190, pos_FP1(2)+240]];
    corners_FP2 = [[pos_FP2(1)-190, pos_FP2(2)-240]; [pos_FP2(1)+190, pos_FP2(2)-240]; [pos_FP2(1)-190, pos_FP2(2)+240]; [pos_FP2(1)+190, pos_FP2(2)+240]];
    corners_FP3 = [[pos_FP3(1)-190, pos_FP3(2)-240]; [pos_FP3(1)+190, pos_FP3(2)-240]; [pos_FP3(1)-190, pos_FP3(2)+240]; [pos_FP3(1)+190, pos_FP3(2)+240]];
    corners_FP4 = [[pos_FP4(1)-190, pos_FP4(2)-240]; [pos_FP4(1)+190, pos_FP4(2)-240]; [pos_FP4(1)-190, pos_FP4(2)+240]; [pos_FP4(1)+190, pos_FP4(2)+240]];

    % Rename force plate COP variables for clarity
    COP_FP1 = matFP(:, 4:6);
    COP_FP2 = matFP(:, 11:13);
    COP_FP3 = matFP(:, 18:20);
    COP_FP4 = matFP(:, 25:27);
    
    bad_fp_mask= ones(n_down, 1);
    
    foot_FP1_bool = (foot_on_fp_cols(:, 1) | foot_on_fp_cols(:, 2));
    foot_FP2_bool = (foot_on_fp_cols(:, 3) | foot_on_fp_cols(:, 4));
    foot_FP3_bool = (foot_on_fp_cols(:, 5) | foot_on_fp_cols(:, 6));
    foot_FP4_bool = (foot_on_fp_cols(:, 7) | foot_on_fp_cols(:, 8));
    
    % Check if CoPs measured are actually on the force plate
    for j=1:n_down
        check_FP1 = (pnt_in_rect(COP_FP1(j, 1:2), corners_FP1)) || (not(foot_FP1_bool(j)));
        check_FP2 = (pnt_in_rect(COP_FP2(j, 1:2), corners_FP2)) || (not(foot_FP2_bool(j)));
        check_FP3 = (pnt_in_rect(COP_FP3(j, 1:2), corners_FP3)) || (not(foot_FP3_bool(j)));
        check_FP4 = (pnt_in_rect(COP_FP4(j, 1:2), corners_FP4)) || (not(foot_FP4_bool(j)));
        
        if ~((check_FP1 && check_FP2) && (check_FP3 && check_FP4))
            bad_fp_mask(j, 1) = 0;
        end
    end
end
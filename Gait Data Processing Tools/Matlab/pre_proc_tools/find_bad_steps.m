function bad_step_mask = find_bad_steps(matMot, foot_on_fp_cols, insole_length_x, insole_length_y, transfos, Tifl, Tifr, pos_FP1, pos_FP2, pos_FP3, pos_FP4, matAllTraj, BW)
    n_down = size(matMot, 1);

    % Define corners of the force plates with extra safety margin
    corners_FP1 = [[pos_FP1(1)-190, pos_FP1(2)-240]; [pos_FP1(1)+190, pos_FP1(2)-240]; [pos_FP1(1)-190, pos_FP1(2)+240]; [pos_FP1(1)+190, pos_FP1(2)+240]];
    corners_FP2 = [[pos_FP2(1)-190, pos_FP2(2)-240]; [pos_FP2(1)+190, pos_FP2(2)-240]; [pos_FP2(1)-190, pos_FP2(2)+240]; [pos_FP2(1)+190, pos_FP2(2)+240]];
    corners_FP3 = [[pos_FP3(1)-190, pos_FP3(2)-240]; [pos_FP3(1)+190, pos_FP3(2)-240]; [pos_FP3(1)-190, pos_FP3(2)+240]; [pos_FP3(1)+190, pos_FP3(2)+240]];
    corners_FP4 = [[pos_FP4(1)-190, pos_FP4(2)-240]; [pos_FP4(1)+190, pos_FP4(2)-240]; [pos_FP4(1)-190, pos_FP4(2)+240]; [pos_FP4(1)+190, pos_FP4(2)+240]];
    corners_FPs = {corners_FP1, corners_FP1, corners_FP2, corners_FP2, corners_FP3, corners_FP3, corners_FP4, corners_FP4};
    corners_FP1_wide = [[pos_FP1(1)-250, pos_FP1(2)-300]; [pos_FP1(1)+250, pos_FP1(2)-300]; [pos_FP1(1)-250, pos_FP1(2)+300]; [pos_FP1(1)+250, pos_FP1(2)+300]];
    corners_FP2_wide = [[pos_FP2(1)-250, pos_FP2(2)-300]; [pos_FP2(1)+250, pos_FP2(2)-300]; [pos_FP2(1)-250, pos_FP2(2)+300]; [pos_FP2(1)+250, pos_FP2(2)+300]];
    corners_FP3_wide = [[pos_FP3(1)-250, pos_FP3(2)-300]; [pos_FP3(1)+250, pos_FP3(2)-300]; [pos_FP3(1)-250, pos_FP3(2)+300]; [pos_FP3(1)+250, pos_FP3(2)+300]];
    corners_FP4_wide = [[pos_FP4(1)-250, pos_FP4(2)-300]; [pos_FP4(1)+250, pos_FP4(2)-300]; [pos_FP4(1)-250, pos_FP4(2)+300]; [pos_FP4(1)+250, pos_FP4(2)+300]];
    corners_FPs_wide = {corners_FP1_wide, corners_FP2_wide, corners_FP3_wide, corners_FP4_wide};
    FP_Z_s = [pos_FP1(3), pos_FP2(3), pos_FP3(3), pos_FP4(3)];

    % Construct CoP points in local ref in homogeneous coords
    % With y-axis in direction of foot, z-axis upwards and x-axis to the right
    COP_INL = cat(2, matMot(:, 3)*insole_length_y, matMot(:, 2)*insole_length_x);
    COP_INL = cat(2, COP_INL, zeros(n_down, 1));
    COP_INR = cat(2, -matMot(:, 6)*insole_length_y, matMot(:, 5)*insole_length_x);
    COP_INR = cat(2, COP_INR, zeros(n_down, 1));
    
    COP_INL_HOM = cat(2, COP_INL, ones(n_down, 1));
    COP_INR_HOM = cat(2, COP_INR, ones(n_down, 1));
    
    % Mask and transformed CoPs
    bad_step_mask = ones(n_down, 1);
    cop_inl_gl = zeros(n_down, 2);
    cop_inr_gl = zeros(n_down, 2);
    
    for j=1:n_down
        for i = 1:2:8 % left
            if foot_on_fp_cols(j, i)
                cop_in = COP_INL_HOM(j, :);
                cop_gl = transfos.Tfgls(4*(j-1)+1:4*j, :)*Tifl*transpose(cop_in);
                cop_inl_gl(j, 1) = cop_gl(1);
                cop_inl_gl(j, 2) = cop_gl(2);
                if not(pnt_in_rect(cop_gl(1:2), corners_FPs{i}))
                    bad_step_mask(j) = 0;
                end
            end
        end
        for i = 2:2:8 % right
            if foot_on_fp_cols(j, i)
                cop_in = COP_INR_HOM(j, :);
                cop_gl = transfos.Tfgrs(4*(j-1)+1:4*j, :)*Tifr*transpose(cop_in);
                cop_inr_gl(j, 1) = cop_gl(1);
                cop_inr_gl(j, 2) = cop_gl(2);
                if not(pnt_in_rect(cop_gl(1:2), corners_FPs{i}))
                    bad_step_mask(j) = 0;
                end
            end
        end
    end

    % Look where two feet are on one force plate
    bad_step_mask_both = ones(n_down, 1);
    left_foot_planted = (matMot(:, 1) > 0.1*BW);
    right_foot_planted = (matMot(:, 4) > 0.1*BW);
    both_feet_planted = (left_foot_planted & right_foot_planted);

    LHEE = matAllTraj(:, 25:27);
    LTOE = matAllTraj(:, 28:30);
    RHEE = matAllTraj(:, 43:45);
    RTOE = matAllTraj(:, 46:48);
    LMFO = ((LTOE-LHEE)/2) + LHEE; % left middle of foot
    RMFO = ((RTOE-RHEE)/2) + RHEE; % right middle of foot

    for j=1:n_down
        for i = 1:4
            if (foot_on_fp_cols(j, 2*i-1) || foot_on_fp_cols(j, 2*i)) && both_feet_planted(j)
                if (pnt_in_rect(LMFO(j, 1:2), corners_FPs_wide{i})) && (pnt_in_rect(RMFO(j, 1:2), corners_FPs_wide{i})) && ((abs(LMFO(j, 3)-FP_Z_s(i)) < 100) && (abs(RMFO(j, 3)-FP_Z_s(i)) < 100))
                    bad_step_mask_both(j) = 0;
                end
            end
        end
    end

    % Join results
    bad_step_mask = bad_step_mask & bad_step_mask_both;
end
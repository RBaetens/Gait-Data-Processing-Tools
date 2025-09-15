function cop = get_cop_from_fm(fm, corners_force_plates, FP_num)
    az0 = -0.041-0.002; % Correction for protective plate
    
    if FP_num == 4
        % Put everything back in force plate frame of reference
        tilt_neg = -10*2*pi/360;
        Fy = cos(tilt_neg)*fm(2, :) - sin(tilt_neg)*fm(3, :);
        Fz = sin(tilt_neg)*fm(2, :)+ cos(tilt_neg)*fm(3, :);
        My = cos(tilt_neg)*fm(5, :) - sin(tilt_neg)*fm(6, :);
        
        % Calculate cop
        Mx_acc = fm(4, :) + Fy*az0;
        My_acc = My - fm(1, :)*az0;
        
        ax = My_acc./Fz; % Don't know why there shouldn't be a minus here, figure this out!!!
        ay_prim = Mx_acc./Fz;

        % Apply rotation
        tilt = 10*2*pi/360;
        az = sin(tilt)*ay_prim;
        ay = cos(tilt)*ay_prim;
    
        middle_FP = repmat(mean(corners_force_plates(:, :, FP_num), 2)', length(ax), 1);
        cop = [middle_FP(:, 1) + ax'*1000, middle_FP(:, 2) + ay'*1000, middle_FP(:, 3) + az'*1000];

    else
        Mx_acc = fm(4, :) + fm(2, :)*az0;
        My_acc = fm(5, :) - fm(1, :)*az0;
        
        ax = My_acc./fm(3, :); % Don't know why there shouldn't be a minus here, figure this out!!!
        ay = Mx_acc./fm(3, :);
    
        middle_FP = repmat(mean(corners_force_plates(:, :, FP_num), 2)', length(ax), 1);
    
        cop = [middle_FP(:, 1) + ax'*1000, middle_FP(:, 2) + ay'*1000, middle_FP(:, 3)];
    end
end
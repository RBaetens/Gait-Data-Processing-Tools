function cop_plate = get_fp_cop_meas(analog_signals, corners_force_plates, FP_num)
    % Returns the cops measured by the force plates
    fc = 25;
    fs = 1000;
    [but_b, but_a] = butter(2, fc/(fs/2), 'low');
    
    if FP_num == 1
        My_acc = filtfilt(but_b, but_a, analog_signals(:, 5));
        Mx_acc = filtfilt(but_b, but_a, analog_signals(:, 4));
        My_acc = My_acc(1:10:size(My_acc, 1));
        Mx_acc = Mx_acc(1:10:size(Mx_acc, 1));

        Fz = analog_signals(:, 29) + analog_signals(:, 30) + ...
            analog_signals(:, 31) + analog_signals(:, 32);
        Fz = filtfilt(but_b, but_a, Fz);
        Fz = Fz(1:10:size(Fz, 1));
        
        middle_FP = repmat(mean(corners_force_plates(:, :, 1), 2)', size(Fz, 1), 1);
        ax = My_acc./Fz;
        ay = Mx_acc./Fz;
        cop_plate = [middle_FP(:, 1) + ax, middle_FP(:, 2) + ay, middle_FP(:, 3)];
    
    elseif FP_num == 2
        My_acc = filtfilt(but_b, but_a, analog_signals(:, 11));
        Mx_acc = filtfilt(but_b, but_a, analog_signals(:, 10));
        My_acc = My_acc(1:10:size(My_acc, 1));
        Mx_acc = Mx_acc(1:10:size(Mx_acc, 1));

        Fz = analog_signals(:, 37) + analog_signals(:, 38) + ...
            analog_signals(:, 39) + analog_signals(:, 40);
        Fz = filtfilt(but_b, but_a, Fz);
        Fz = Fz(1:10:size(Fz, 1));
        
        middle_FP = repmat(mean(corners_force_plates(:, :, 2), 2)', size(Fz, 1), 1);
        ax = My_acc./Fz;
        ay = Mx_acc./Fz;
        cop_plate = [middle_FP(:, 1) + ax, middle_FP(:, 2) + ay, middle_FP(:, 3)];

    elseif FP_num == 3
        My_acc = filtfilt(but_b, but_a, analog_signals(:, 17));
        Mx_acc = filtfilt(but_b, but_a, analog_signals(:, 16));
        My_acc = My_acc(1:10:size(My_acc, 1));
        Mx_acc = Mx_acc(1:10:size(Mx_acc, 1));

        Fz = analog_signals(:, 45) + analog_signals(:, 46) + ...
            analog_signals(:, 47) + analog_signals(:, 48);
        Fz = filtfilt(but_b, but_a, Fz);
        Fz = Fz(1:10:size(Fz, 1));
        
        middle_FP = repmat(mean(corners_force_plates(:, :, 3), 2)', size(Fz, 1), 1);
        ax = My_acc./Fz;
        ay = Mx_acc./Fz;
        cop_plate = [middle_FP(:, 1) + ax, middle_FP(:, 2) + ay, middle_FP(:, 3)];

    elseif FP_num == 4
        % Parameters
        a = 0.21;
        b = 0.26;
        az0 = -0.041-0.002; % Correction for protective plate

        % Calculate everything in force plate frame of reference
        Mx = b * (analog_signals(:, 53) + analog_signals(:, 54) - ...
            analog_signals(:, 55) - analog_signals(:, 56));
        Mx = filtfilt(but_b, but_a, Mx);
        Mx = Mx(1:10:size(Mx, 1));
        My_prim = a * (-analog_signals(:, 53) + analog_signals(:, 54) + ...
            analog_signals(:, 55) - analog_signals(:, 56));
        My_prim = filtfilt(but_b, but_a, My_prim);
        My_prim = My_prim(1:10:size(My_prim, 1));

        Fx = (analog_signals(:, 49) + analog_signals(:, 50));
        Fx = filtfilt(but_b, but_a, Fx);
        Fx = Fx(1:10:size(Fx, 1));
        Fy_prim = (analog_signals(:, 51) + analog_signals(:, 52));
        Fy_prim = filtfilt(but_b, but_a, Fy_prim);
        Fy_prim = Fy_prim(1:10:size(Fy_prim, 1));
        Fz_prim = analog_signals(:, 53) + analog_signals(:, 54) + ...
            analog_signals(:, 55) + analog_signals(:, 56);
        Fz_prim = filtfilt(but_b, but_a, Fz_prim);
        Fz_prim = Fz_prim(1:10:size(Fz_prim, 1));

        Mx_acc = Mx + Fy_prim*az0;
        My_acc_prim = My_prim - Fx*az0;

        ax = My_acc_prim./Fz_prim;
        ay_prim = Mx_acc./Fz_prim;

        % Apply rotation
        tilt = 10*2*pi/360;
        az = sin(tilt)*ay_prim;
        ay = cos(tilt)*ay_prim;

        middle_FP = repmat(mean(corners_force_plates(:, :, 4), 2)', size(Fz_prim, 1), 1);
        cop_plate = [middle_FP(:, 1) + ax*1000, middle_FP(:, 2) + ay*1000, middle_FP(:, 3) + az*1000];
    end
end
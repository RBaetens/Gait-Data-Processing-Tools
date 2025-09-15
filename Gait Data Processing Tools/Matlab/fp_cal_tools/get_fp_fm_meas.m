function PILS_meas_mat = get_fp_fm_meas(analog_signals, FP_num)
    % Returns the forces and moments measured by the force plates
    a = 0.21;
    b = 0.26;

    fc = 25;
    fs = 100;
    [but_b, but_a] = butter(2, fc/(fs/2), 'low');

    if FP_num == 1
        Fx = -(analog_signals(:, 25) + analog_signals(:, 26)); % - for 180 deg rotation
        Fy = -(analog_signals(:, 27) + analog_signals(:, 28)); % - for 180 deg rotation
        Fz = analog_signals(:, 29) + analog_signals(:, 30) + ...
            analog_signals(:, 31) + analog_signals(:, 32);
        
        Mx = -b * (analog_signals(:, 29) + analog_signals(:, 30) - ...
            analog_signals(:, 31) - analog_signals(:, 32)); % - for 180 deg rotation
        My = -a * (-analog_signals(:, 29) + analog_signals(:, 30) + ...
            analog_signals(:, 31) - analog_signals(:, 32)); % - for 180 deg rotation
        Mz = b * (-analog_signals(:, 25) + analog_signals(:, 26)) + ...
            a * (analog_signals(:, 27) - analog_signals(:, 28));

    elseif FP_num == 2
        Fx = -(analog_signals(:, 33) + analog_signals(:, 34)); % - for 180 deg rotation
        Fy = -(analog_signals(:, 35) + analog_signals(:, 36)); % - for 180 deg rotation
        Fz = analog_signals(:, 37) + analog_signals(:, 38) + ...
            analog_signals(:, 39) + analog_signals(:, 40);
        
        Mx = -b * (analog_signals(:, 37) + analog_signals(:, 38) - ...
            analog_signals(:, 39) - analog_signals(:, 40)); % - for 180 deg rotation
        My = -a * (-analog_signals(:, 37) + analog_signals(:, 38) + ...
            analog_signals(:, 39) - analog_signals(:, 40)); % - for 180 deg rotation
        Mz = b * (-analog_signals(:, 33) + analog_signals(:, 34)) + ...
            a * (analog_signals(:, 35) - analog_signals(:, 36));
    
    elseif FP_num == 3
        Fx = -(analog_signals(:, 41) + analog_signals(:, 42)); % - for 180 deg rotation
        Fy = -(analog_signals(:, 43) + analog_signals(:, 44)); % - for 180 deg rotation
        Fz = analog_signals(:, 45) + analog_signals(:, 46) + ...
            analog_signals(:, 47) + analog_signals(:, 48);
        
        Mx = -b * (analog_signals(:, 45) + analog_signals(:, 46) - ...
            analog_signals(:, 47) - analog_signals(:, 48)); % - for 180 deg rotation
        My = -a * (-analog_signals(:, 45) + analog_signals(:, 46) + ...
            analog_signals(:, 47) - analog_signals(:, 48)); % - for 180 deg rotation
        Mz = b * (-analog_signals(:, 41) + analog_signals(:, 42)) + ...
            a * (analog_signals(:, 43) - analog_signals(:, 44));

    elseif FP_num == 4
        Fx = (analog_signals(:, 49) + analog_signals(:, 50));
        Fy_prim = (analog_signals(:, 51) + analog_signals(:, 52));
        Fz_prim = analog_signals(:, 53) + analog_signals(:, 54) + ...
            analog_signals(:, 55) + analog_signals(:, 56);
        
        Mx = b * (analog_signals(:, 53) + analog_signals(:, 54) - ...
            analog_signals(:, 55) - analog_signals(:, 56));
        My_prim = a * (-analog_signals(:, 53) + analog_signals(:, 54) + ...
            analog_signals(:, 55) - analog_signals(:, 56));
        Mz_prim = b * (-analog_signals(:, 49) + analog_signals(:, 50)) + ...
            a * (analog_signals(:, 51) - analog_signals(:, 52));
        
        tilt = 10*2*pi/360;
        Fy = cos(tilt)*Fy_prim - sin(tilt)*Fz_prim;
        Fz = sin(tilt)*Fy_prim + cos(tilt)*Fz_prim;
        My = cos(tilt)*My_prim - sin(tilt)*Mz_prim;
        Mz = sin(tilt)*My_prim + cos(tilt)*Mz_prim;
    end

    PILS_meas_mat = [Fx, Fy, Fz, Mx, My, Mz];
    PILS_meas_mat = PILS_meas_mat(1:10:size(PILS_meas_mat, 1), :);
    
    for i=1:6
        PILS_meas_mat(:, i) = filtfilt(but_b, but_a, PILS_meas_mat(:, i));
    end
end
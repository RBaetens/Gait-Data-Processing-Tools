function PILS_ref_mat = get_pole_fm_meas(analog_signals, markers, corners_force_plates, FP_num)
    % Returns the reference wrench

    % Rod: force
    pot_rod = analog_signals(:, 69);
    
    fc = 5;
    fs = 1000;
    [but_b, but_a] = butter(2, fc/(fs/2), 'low');
    pot_rod_filt = filtfilt(but_b, but_a, pot_rod);
    
    Ftot_rod = -9.81 * (11.65*pot_rod_filt + 0.0266);
    Ftot_rod = Ftot_rod(1:10:size(Ftot_rod));
    
    % Rod: markers
    side_marker = squeeze(markers(:, 1, :))';
    up_marker = squeeze(markers(:, 2, :))';
    center_marker = squeeze(markers(:, 3, :))';
    
    x_axis_rod = side_marker - center_marker;
    for i=1:size(x_axis_rod, 1)
        x_axis_rod(i, :) = x_axis_rod(i, :) / norm(x_axis_rod(i, :));
    end
    
    z_axis_rod = (up_marker - 4.1*x_axis_rod) - center_marker; % correction to make axis parallel with bar
    for i=1:size(z_axis_rod, 1)
        z_axis_rod(i, :) = z_axis_rod(i, :) / norm(z_axis_rod(i, :));
    end
    
    y_axis_rod = zeros(size(x_axis_rod, 1), 3);
    for i=1:size(y_axis_rod)
        y_axis_rod(i, :) = cross(z_axis_rod(i, :), x_axis_rod(i, :));
    end
    
    tip = center_marker + 30.1*x_axis_rod + 19.2*y_axis_rod - (361-0)*z_axis_rod; % correction for nut and protective plate
    
    % Get reference forces and moments
    Fx_rod = Ftot_rod .* -z_axis_rod(:, 1);
    Fy_rod = Ftot_rod .* -z_axis_rod(:, 2);
    Fz_rod = Ftot_rod .* -z_axis_rod(:, 3);
    Fref = [Fx_rod, Fy_rod, Fz_rod];

    g_dir = [zeros(size(Fz_rod, 1), 2), -ones(size(Fz_rod, 1), 1)];
    Non_axial_weight_dir = zeros(size(Fz_rod, 1), 3);
    for i=1:size(Fz_rod, 1)
        Non_axial_weight_dir(i, :) = g_dir(i, :) - dot(g_dir(i, :), -z_axis_rod(i, :))*(-z_axis_rod(i, :));
    end

    Fref = Fref + 1.91*9.81*(0.61/1.195)*Non_axial_weight_dir;
    Fref = -Fref;
    
    middle_FP = repmat(mean(corners_force_plates(:, :, FP_num), 2)', size(Fz_rod, 1), 1);
    moment_arm = (tip - middle_FP)*0.001; % from mm to m
    Mref = zeros(size(Fref, 1), 3);
    for i=1:size(Mref, 1)
        Mref(i, :) = cross(moment_arm(i, :), Fref(i, :));
    end

    PILS_ref_mat = [Fref, Mref];
    PILS_ref_mat(:, 2) = -PILS_ref_mat(:, 2); % I don't know why it reverses the sign for the y_axis only
    PILS_ref_mat(:, 5) = -PILS_ref_mat(:, 5); % Figure this out!!!
end
function tip = get_pole_tip_meas(markers)   
    % Return coordinates for the tip of the pole
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
end
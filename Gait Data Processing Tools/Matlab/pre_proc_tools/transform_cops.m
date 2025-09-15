% Transform C.o.P.s from global frame of reference to insole frame of
% reference and normalize them with respect to insole length

function COPS = transform_cops(matFP, transfos, Tfil, Tfir, foot_on_fp_cols, insole_length_x, insole_length_y)
    % Initialize result variables
    n_pts = size(matFP, 1);
    COPS_L = zeros(n_pts, 2);
    COPS_R = zeros(n_pts, 2);

    % Express COPS in homogeneous coordinates
    COPS_ORIG = [matFP(:, 4:6), ones(n_pts, 1), matFP(:, 11:13), ones(n_pts, 1), ...
        matFP(:, 18:20), ones(n_pts, 1), matFP(:, 25:27), ones(n_pts, 1)];
    
    for j=1:n_pts
        for i = 1:4
            % Left foot
            if foot_on_fp_cols(j, 2*i-1)
                % Determine transformation
                Tgil = Tfil*transfos.Tgfls(4*(j-1)+1:4*j, :);
                
                % Determine COP in original insole reference frame
                COP_L = Tgil*COPS_ORIG(j, 4*i-3:4*i)';
                COPS_L(j, :) = COP_L(1:2)';
            end

            % Right foot
            if foot_on_fp_cols(j, 2*i)
                % Determine transformation
                Tgir = Tfir*transfos.Tgfrs(4*(j-1)+1:4*j, :);
                
                % Determine COP in original insole reference frame
                COP_R = Tgir*COPS_ORIG(j, 4*i-3:4*i)';
                COPS_R(j, :) = COP_R(1:2)';
            end
        end
    end

    % Scale, switch some axes and combine
    COPS_L(:, 1) = COPS_L(:, 1)/insole_length_y;
    COPS_L(:, 2) = COPS_L(:, 2)/insole_length_x;
    COPS_R(:, 1) = COPS_R(:, 1)/insole_length_y;
    COPS_R(:, 2) = COPS_R(:, 2)/insole_length_x;
    COPS = [COPS_L(:, 2), COPS_L(:, 1), COPS_R(:, 2), -COPS_R(:, 1)];
end
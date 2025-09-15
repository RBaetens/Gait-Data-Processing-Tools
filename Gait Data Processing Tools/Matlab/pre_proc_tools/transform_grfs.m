function GRF_LR = transform_grfs(matFP, transfos, Tifl, Tifr, foot_on_fp_cols)
    n_down = size(matFP, 1);
    GRF_left = zeros(n_down, 4);
    GRF_right = zeros(n_down, 4);
    
    yax_glo_2d = [0; 1];
    
    for j=1:n_down
        for i = 1:4
            % Left foot
            if foot_on_fp_cols(j, 2*i-1)
                % Determine local y-axis
                Tigl = transfos.Tfgls(4*(j-1)+1:4*j, :)*Tifl;
                yax_loc_2d = Tigl(1:2, 2);
                % Determine foot progression angle
                alpha = acos( (dot(yax_loc_2d, yax_glo_2d)) / (norm(yax_loc_2d)*norm(yax_glo_2d)) );
                if yax_loc_2d(1) < 0
                    alpha = -alpha;
                end
                % Determine GRF in original insole reference frame
                rot_alpha = [[cos(alpha), -sin(alpha), 0]; [sin(alpha), cos(alpha), 0]; [0, 0 1]];
                grf_left = rot_alpha*transpose(matFP(j, 7*(i-1)+1:7*(i-1)+3));
                GRF_left(j, 1) = grf_left(2);
                GRF_left(j, 2) = grf_left(1);
                GRF_left(j, 3) = grf_left(3);
                GRF_left(j, 4) = matFP(j, 7*i); % Frictional torque along z-axis
            end
            % Right foot
            if foot_on_fp_cols(j, 2*i)
                % Determine local y-axis
                Tigr = transfos.Tfgrs(4*(j-1)+1:4*j, :)*Tifr;
                yax_loc_2d = Tigr(1:2, 2);
                % Determine foot progression angle
                alpha = acos( (dot(yax_loc_2d, yax_glo_2d)) / (norm(yax_loc_2d)*norm(yax_glo_2d)) );
                if yax_loc_2d(1) < 0
                    alpha = -alpha;
                end
                % Determine GRF in original insole reference frame
                rot_alpha = [[cos(alpha), -sin(alpha), 0]; [sin(alpha), cos(alpha), 0]; [0, 0 1]];
                grf_right = rot_alpha*transpose(matFP(j, 7*(i-1)+1:7*(i-1)+3));
                GRF_right(j, 1) = grf_right(2);
                GRF_right(j, 2) = -grf_right(1);
                GRF_right(j, 3) = grf_right(3);
                GRF_right(j, 4) = matFP(j, 7*i); % Frictional torque along z-axis
            end
        end
    end

    % Right now GRF_left and GRF_right express the force that the person exerts
    % on the ground for some reason, so reverse this
    GRF_left = -GRF_left;
    GRF_right = -GRF_right;
    GRF_LR = [GRF_left, GRF_right];
end
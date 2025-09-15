% Calculate transformations from ankle to global f.o.r.
function transfos = heel_to_global(matAllTraj, valid_mask_feet)
    %
    yax_glo_2d = [0; 1];

    % Rename some variables for clarity
    LHEE = matAllTraj(:, 25:27);
    LTOE = matAllTraj(:, 28:30);
    RHEE = matAllTraj(:, 43:45);
    RTOE = matAllTraj(:, 46:48);

    % Create variables to store transformations
    n_down = size(matAllTraj, 1);
    Tfgls = zeros(4*n_down, 4); % Transformations from left ankle to global f.o.r.
    Tfgrs = zeros(4*n_down, 4); % Transformations from right ankle to global f.o.r.
    Tgfls = zeros(4*n_down, 4); % Transformations from global to left ankle f.o.r.
    Tgfrs = zeros(4*n_down, 4); % Transformations from global to right ankle f.o.r.

    for i = 1:n_down
        if valid_mask_feet(i)
            try
                % Define axes of ankle f.o.r. for each foot with
                % y-axis from heel marker to toe marker
                % z-axis vertical
                % All vectors normalised
                % Origin in HEE
    
                % Calculate y-axes
                lht = LTOE(i, :)-LHEE(i, :);
                ly = lht/norm(lht);
                rht = RTOE(i, :)-RHEE(i, :);
                ry = rht/norm(rht);

                % Determine left foot progression angle
                lalpha = acos( (dot(ly(1:2), yax_glo_2d)) / (norm(ly(1:2))*norm(yax_glo_2d)) );
                if ly(1) < 0
                    lalpha = -lalpha;
                end
                rot_lalpha = [[cos(lalpha), -sin(lalpha), 0]; [sin(lalpha), cos(lalpha), 0]; [0, 0 1]];

                % Determine right foot progression angle
                ralpha = acos( (dot(ry(1:2), yax_glo_2d)) / (norm(ry(1:2))*norm(yax_glo_2d)) );
                if ry(1) < 0
                    ralpha = -ralpha;
                end
                rot_ralpha = [[cos(ralpha), -sin(ralpha), 0]; [sin(ralpha), cos(ralpha), 0]; [0, 0 1]];
    
                % Assemble transformation matrix T (ankle f.o.r. --> global f.o.r.)
                Tfgl = zeros(4, 4);
                Tfgr = zeros(4, 4);
                Tfgl(1:3, 1:3) = rot_lalpha;
                Tfgl(1:3, 4) = LHEE(i, :)';
                Tfgl(4, 4) = 1;
                Tfgr(1:3, 1:3) = rot_ralpha;
                Tfgr(1:3, 4) = RHEE(i, :)';
                Tfgr(4, 4) = 1;
            catch
                Tfgl = eye(4);
                Tfgr = eye(4);
            end
    
            % Find inverse transformation matrices T (global f.o.r. --> ankle f.o.r.)
            Tgfl = trans_mat_inv(Tfgl);
            Tgfr = trans_mat_inv(Tfgr);
    
            % Store transformation matrices
            Tfgls(4*(i-1)+1:4*i, :) = Tfgl;
            Tfgrs(4*(i-1)+1:4*i, :) = Tfgr;
            Tgfls(4*(i-1)+1:4*i, :) = Tgfl;
            Tgfrs(4*(i-1)+1:4*i, :) = Tgfr;
        end
    end
    
    % Assign transformations to structure
    transfos.Tfgls = Tfgls;
    transfos.Tfgrs = Tfgrs;
    transfos.Tgfls = Tgfls;
    transfos.Tgfrs = Tgfrs;
end
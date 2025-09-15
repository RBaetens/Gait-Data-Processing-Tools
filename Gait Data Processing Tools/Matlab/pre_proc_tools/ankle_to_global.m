% Calculate transformations from ankle to global f.o.r.
function transfos = ankle_to_global(matAllTraj, valid_mask_feet, AW, mark_diam, mark_base)
    % Rename some variables for clarity
    LANK = matAllTraj(:, 22:24);
    LHEE = matAllTraj(:, 25:27);
    LTOE = matAllTraj(:, 28:30);
    RANK = matAllTraj(:, 40:42);
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
                % z-axis perpendicular on y-axis and upwards through AJC
                % x-axis as cross product of y- and z-axis
                % All vectors normalised
                % Origin in AJC
    
                % Calculate y-axis
                lht = LTOE(i, :)-LHEE(i, :);
                ly = lht/norm(lht);
                rht = RTOE(i, :)-RHEE(i, :);
                ry = rht/norm(rht);
                
                % Calculate point on y-axis closest to ankle
                sl = (dot(ly, LANK(i, :)) - dot(ly, LHEE(i, :))) / dot(ly, lht);
                sr = (dot(ry, RANK(i, :)) - dot(ry, RHEE(i, :))) / dot(ry, rht);
                pl = [LHEE(i, 1)+sl*(LTOE(i, 1)-LHEE(i, 1)), LHEE(i, 2)+sl*(LTOE(i, 2)-LHEE(i, 2)), LHEE(i, 3)+sl*(LTOE(i, 3)-LHEE(i, 3))];
                pr = [RHEE(i, 1)+sr*(RTOE(i, 1)-RHEE(i, 1)), RHEE(i, 2)+sr*(RTOE(i, 2)-RHEE(i, 2)), RHEE(i, 3)+sr*(RTOE(i, 3)-RHEE(i, 3))];
                
                % Alternatively:
                % lah = LHEE(i, :)-LANK(i, :);
                % rah = RHEE(i, :)-RANK(i, :);
                % sl = -dot(ly, lah)/dot(ly, lht);
                % sr = -dot(ry, rah)/dot(ry, rht);
                % pl = [LHEE(i, 1)+sl*(LTOE(i, 1)-LHEE(i, 1)), LHEE(i, 2)+sl*(LTOE(i, 2)-LHEE(i, 2)), LHEE(i, 3)+sl*(LTOE(i, 3)-LHEE(i, 3))];
                % pr = [RHEE(i, 1)+sr*(RTOE(i, 1)-RHEE(i, 1)), RHEE(i, 2)+sr*(RTOE(i, 2)-RHEE(i, 2)), RHEE(i, 3)+sr*(RTOE(i, 3)-RHEE(i, 3))];
        
                % Calculate (estimate) AJC
                cl = (((AW+mark_diam)/2 + mark_base)^2) / (1+(ly(1)/ly(2))^2);
                cr = (((AW+mark_diam)/2 + mark_base)^2) / (1+(ry(1)/ry(2))^2);
                ajcl1 = [LANK(i, 1) + sqrt(cl), LANK(i, 2) - (ly(1)/ly(2))*sqrt(cl), LANK(i, 3)];
                ajcl2 = [LANK(i, 1) - sqrt(cl), LANK(i, 2) + (ly(1)/ly(2))*sqrt(cl), LANK(i, 3)];
                ajcr1 = [RANK(i, 1) + sqrt(cr), RANK(i, 2) - (ry(1)/ry(2))*sqrt(cr), RANK(i, 3)];
                ajcr2 = [RANK(i, 1) - sqrt(cr), RANK(i, 2) + (ry(1)/ry(2))*sqrt(cr), RANK(i, 3)];
                
                % Alternatively:
                % cl = dot(ly, LANK(i, :));
                % aal = (1 + (ly(1)/ly(2))^2);
                % bbl = ((2*ly(1)*ly(3)*LANK(i, 3) - 2*cl*ly(1)) / ly(2)^2) - 2*LANK(i, 1) + ((2*LANK(i, 2)*ly(1)) / ly(2));
                % ccl = LANK(i, 1)^2 + LANK(i, 2)^2 + ((cl^2 + (ly(3)^2)*(LANK(i, 3)^2) - 2*cl*ly(3)*LANK(i, 3)) / (ly(2)^2)) - (2*LANK(i, 2)/ly(2))*(cl-ly(3)*LANK(i, 3)) - ((AW+mark_diam)/2 + mark_base)^2;
                % 
                % Dl = bbl^2-4*aal*ccl;
                % ajcxl1 = (-bbl + sqrt(Dl))/(2*aal);
                % ajcxl2 = (-bbl - sqrt(Dl))/(2*aal);
                % ajcl1 = [ajcxl1, (cl - (ly(3)*LANK(i, 3) + ly(1)*ajcxl1)) / ly(2), LANK(i, 3)];
                % ajcl2 = [ajcxl2, (cl - (ly(3)*LANK(i, 3) + ly(1)*ajcxl2)) / ly(2), LANK(i, 3)];
                % 
                % cr = dot(ry, RANK(i, :));
                % aar = (1 + (ry(1)/ry(2))^2);
                % bbr = ((2*ry(1)*ry(3)*RANK(i, 3) - 2*cr*ry(1)) / ry(2)^2) - 2*RANK(i, 1) + ((2*RANK(i, 2)*ry(1)) / ry(2));
                % ccr = RANK(i, 1)^2 + RANK(i, 2)^2 + ((cr^2 + (ry(3)^2)*(RANK(i, 3)^2) - 2*cr*ry(3)*RANK(i, 3)) / (ry(2)^2)) - (2*RANK(i, 2)/ry(2))*(cr-ry(3)*RANK(i, 3)) - ((AW+mark_diam)/2 + mark_base)^2;
                % 
                % Dr = bbr^2-4*aar*ccr;
                % ajcxr1 = (-bbr + sqrt(Dr))/(2*aar);
                % ajcxr2 = (-bbr - sqrt(Dr))/(2*aar);
                % ajcr1 = [ajcxr1, (cr - (ry(3)*RANK(i, 3) + ry(1)*ajcxr1)) / ry(2), RANK(i, 3)];
                % ajcr2 = [ajcxr2, (cr - (ry(3)*RANK(i, 3) + ry(1)*ajcxr2)) / ry(2), RANK(i, 3)];

                if norm(ajcl1-pl) < norm(ajcl2-pl)
                    ajcl = ajcl1;
                else
                    ajcl = ajcl2;
                end

                if norm(ajcr1-pr) < norm(ajcr2-pr)
                    ajcr = ajcr1;
                else
                    ajcr = ajcr2;
                end
    
                % Calculate z-axis
                lz = ajcl-pl;
                lz = lz/norm(lz);
                rz = ajcr-pr;
                rz = rz/norm(rz);
    
                % Calculate x-axis
                lx = cross(ly, lz);
                rx = cross(ry, rz);
    
                % Assemble transformation matrix T (ankle f.o.r. --> global f.o.r.)
                Tfgl = zeros(4, 4);
                Tfgr = zeros(4, 4);
                Tfgl(1:3, 1) = lx;
                Tfgl(1:3, 2) = ly;
                Tfgl(1:3, 3) = lz;
                Tfgl(1:3, 4) = ajcl';
                Tfgl(4, 4) = 1;
                Tfgr(1:3, 1) = rx;
                Tfgr(1:3, 2) = ry;
                Tfgr(1:3, 3) = rz;
                Tfgr(1:3, 4) = ajcr';
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
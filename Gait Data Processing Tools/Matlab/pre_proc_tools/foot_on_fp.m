function foot_on_fp_cols = foot_on_fp(FZs, COPXYZs, LHEEXYZ, LTOEXYZ, RHEEXYZ, RTOEXYZ, thresh)
    % This function can be used to determine at what times, which foot is on which force plate.
    
    % Input arguments
    % FZs: vertical force component measured by each force plate (n_time_steps, n_force_plates)
    % COPXYZs: CoP coordinates for each force plate (n_time_steps, 3*n_force_plates)
    % LHEEXYZ: left heel marker coordinates (n_time_steps, 3)
    % LTOEXYZ: left toe marker coordinates (n_time_steps, 3)
    % RHEEXYZ: right heel marker coordinates (n_time_steps, 3)
    % RTOEXYZ: right toe marker coordinates (n_time_steps, 3)

    % Output
    % foot_on_fp_cols: mask indicating when which foot is on which force
    % plate, (left foot on FP1, right foot on FP1, left foot on FP2, ...),
    % (n_time_steps, 2*n_force_plates)

    LMOFXYZ = (LHEEXYZ + LTOEXYZ)/2; % coordinate of left 'middle of foot'
    RMOFXYZ = (RHEEXYZ + RTOEXYZ)/2; % coordinate of right 'middle of foot'
    
    n_time_steps = size(FZs, 1);
    n_fp = size(FZs, 2);
    foot_on_fp_cols = zeros(n_time_steps, 2*n_fp);
    
    for i = 1:n_fp
        % Distance of each foot to each FP CoP
        COPXYZ = COPXYZs(:, 3*(i-1)+1:3*i);

        left_dist2 = sum(((LMOFXYZ - COPXYZ).^2), 2);
        right_dist2 = sum(((RMOFXYZ - COPXYZ).^2), 2);
        
        left_closest = left_dist2 < right_dist2; % 1 if left foot is closest
        right_closest = not(left_closest); % 1 if right foot is closest

        % See if the FP is active
        FZ = FZs(:, i);
        fp_active = FZ > thresh; % 1 if FP is active
        
        % Add to output
        foot_on_fp_cols(:, 2*(i-1)+1:2*i) = [(fp_active & left_closest), (fp_active & right_closest)];
    end
end
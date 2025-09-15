function [tabAllTraj, tabBut, tabEMG, tabFP] = check_measurement_lengths(tabAllTraj, tabBut, tabEMG, tabFP, f_AllTraj, f_But, f_EMG, f_FP)
    n_check_AllTraj = size(tabAllTraj, 1)-3;
    n_check_But = size(tabBut, 1)-3;
    n_check_EMG = size(tabEMG, 1)-3;
    n_check_FP = size(tabFP, 1)-3;
    
    T_check_AllTraj = n_check_AllTraj/f_AllTraj;
    T_check_But = n_check_But/f_But;
    T_check_EMG = n_check_EMG/f_EMG;
    T_check_FP = n_check_FP/f_FP;
    T_list = [T_check_AllTraj, T_check_But, T_check_EMG, T_check_FP];
    [T_mode, f_mode] = mode(T_list);
    
    % More than two different lengths of measuring time
    if f_mode < 3
        error("Vicon really messed up, fix manually.")
    % Two different lengths of measuring time
    elseif f_mode == 3
        % Check if one is shorter (what you can expect)
        if T_check_AllTraj < T_mode
            error("Problem with the measurement time of markers, best to fix manually.")
        elseif T_check_But < T_mode
            % Get frame column
            frames_short = tabBut(4:end, 1);
            frames_long = tabEMG(4:end, 1);
            n_diff = n_check_EMG-n_check_But;
    
            % Check if the difference lies at the start or end
            check_start = sum((frames_long(1:2*n_diff) ~= frames_short(1:2*n_diff)));
            check_end = sum((frames_long(end-2*n_diff+1:end) ~= frames_short(end-2*n_diff+1:end)));
            
            % Append NaNs to match lengths
            if check_start ~= 0
                tabBut = [nan(n_diff, size(tabBut, 2)); tabBut];
            elseif check_end ~= 0
                tabBut = [tabBut; nan(n_diff, size(tabBut, 2))];
            end
    
        elseif T_check_EMG < T_mode
            % Get frame column
            frames_short = tabEMG(4:end, 1);
            frames_long = tabBut(4:end, 1);
            n_diff = n_check_But-n_check_EMG;
    
            % Check if the difference lies at the start or end
            check_start = sum((frames_long(1:2*n_diff) ~= frames_short(1:2*n_diff)));
            check_end = sum((frames_long(end-2*n_diff+1:end) ~= frames_short(end-2*n_diff+1:end)));
            
            % Append NaNs to match lengths
            if check_start ~= 0
                tabEMG = [nan(n_diff, size(tabEMG, 2)); tabEMG];
            elseif check_end ~= 0
                tabEMG = [tabEMG; nan(n_diff, size(tabEMG, 2))];
            end
    
        elseif T_check_FP < T_mode
            % Get frame column
            frames_short = tabFP(4:end, 1);
            if f_FP == f_AllTraj
                frames_long = tabAllTraj(4:end, 1);
                n_diff = n_check_AllTraj-n_check_FP;
            elseif f_FP == f_EMG
                frames_long = tabEMG(4:end, 1);
                n_diff = n_check_EMG-n_check_FP;
            else
                error("Frequencies do not match, fix manually.")
            end
    
            % Check if the difference lies at the start or end
            check_start = sum((frames_long(1:2*n_diff) ~= frames_short(1:2*n_diff)));
            check_end = sum((frames_long(end-2*n_diff+1:end) ~= frames_short(end-2*n_diff+1:end)));
            
            % Append NaNs to match lengths
            if check_start ~= 0
                tabFP = [nan(n_diff, size(tabFP, 2)); tabFP];
            elseif check_end ~= 0
                tabFP = [tabFP; nan(n_diff, size(tabFP, 2))];
            end
    
        % One is longer
        else
            error("Vicon really messed up, fix manually.")
        end
    end
end
function start_stop_reps = update_start_stop_reps(start_stop_reps, ON_OFF_ERROR)
    % Deletes reps where there was an error
    
    % See if any error happened
    error_not_pressed_ind = [];
    for i = 1:size(start_stop_reps, 1)
        error_pressed = false;
        for j = 1:size(ON_OFF_ERROR, 1)
            if ON_OFF_ERROR(j, 1) >= start_stop_reps(i, 1) && ON_OFF_ERROR(j, 1) <= start_stop_reps(i, 2)
                error_pressed = true;
            end
        end
        if not(error_pressed)
            error_not_pressed_ind = [error_not_pressed_ind, i];
        end
    end

    % Update start_stop_reps
    start_stop_reps = start_stop_reps(error_not_pressed_ind, :);
end
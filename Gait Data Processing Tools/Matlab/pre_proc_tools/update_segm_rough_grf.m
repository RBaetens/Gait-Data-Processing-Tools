function ON_OFF_SEGM_ROUGH = update_segm_rough_grf(ON_OFF_SEGM_ROUGH, ON_OFF_ERROR)
    % Deletes parts where error button was pressed

    % See if any fine segmentation happened
    ON_OFF_SEGM_ROUGH_NEW = ON_OFF_SEGM_ROUGH(1, :);
    for i = 2:size(ON_OFF_SEGM_ROUGH, 1)
        error_pressed = false;
        for j = 1:size(ON_OFF_ERROR, 1)
            if ON_OFF_SEGM_ROUGH(i-1, 2) <= ON_OFF_ERROR(j, 1) && ON_OFF_ERROR(j, 1) <= ON_OFF_SEGM_ROUGH(i, 1)
                error_pressed = true;
            end
        end
        if ~error_pressed
            ON_OFF_SEGM_ROUGH_NEW = [ON_OFF_SEGM_ROUGH_NEW; ON_OFF_SEGM_ROUGH(i, :)];
        else
            ON_OFF_SEGM_ROUGH_NEW(end, 2) = ON_OFF_SEGM_ROUGH(i, 2);
        end
    end

    % Update ON_OFF_SEGM_ROUGH
    ON_OFF_SEGM_ROUGH = ON_OFF_SEGM_ROUGH_NEW;
end
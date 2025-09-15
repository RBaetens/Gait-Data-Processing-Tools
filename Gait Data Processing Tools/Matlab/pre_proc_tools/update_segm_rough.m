function ON_OFF_SEGM_ROUGH = update_segm_rough(ON_OFF_SEGM_ROUGH, ON_OFF_SEGM_FINE)
    % Deletes presses where no fine segmentation happened (i.e. remove
    % accidental double presses)
    
    % See if any fine segmentation happened
    segm_fine_pressed_ind = [];
    for i = 1:size(ON_OFF_SEGM_ROUGH, 1)-1
        segm_fine_pressed = false;
        for j = 1:size(ON_OFF_SEGM_FINE, 1)
            if ON_OFF_SEGM_FINE(j, 1) >= ON_OFF_SEGM_ROUGH(i, 2) && ON_OFF_SEGM_FINE(j, 1) <= ON_OFF_SEGM_ROUGH(i+1, 1)
                segm_fine_pressed = true;
            end
        end
        if segm_fine_pressed
            segm_fine_pressed_ind = [segm_fine_pressed_ind, i];
        end
    end
    segm_fine_pressed_ind = [segm_fine_pressed_ind, i+1];

    % Update ON_OFF_SEGM_ROUGH
    ON_OFF_SEGM_ROUGH = ON_OFF_SEGM_ROUGH(segm_fine_pressed_ind, :);
end
% used to determine wether part of a force plate signal belongs to a task execution
function active = check_active(FP_Fv, ind, stop_ind, n_buf, fwbw, Th) % fwbw, 1 = forward check, 0 = backward check
    if fwbw
        if 1+ind+n_buf <= stop_ind
            active = sum((FP_Fv(1+ind:1+ind+n_buf, :) > Th), "all") > 1;
        else
            active = false;
        end
    else
        if stop_ind-ind-n_buf >= 1
            active = sum((FP_Fv(stop_ind-ind-n_buf:stop_ind-ind, :) > Th), "all") > 1;
        else
            active = false;
        end
    end
end
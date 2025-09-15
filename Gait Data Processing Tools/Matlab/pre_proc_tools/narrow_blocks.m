% This function makes blocks narrower in a block wave, per column
function narrow_waves = narrow_blocks(waves, n_narrower)
    n_pts = size(waves, 1);
    n_col = size(waves, 2);
    narrow_waves = zeros(n_pts, n_col);

    for i = 1:n_col
        j = 1;
        while (j <= n_pts)
            if waves(j, i)
                start_ind_block = j;
                j = j+1;
                while waves(j, i)
                    j = j+1;
                end
                j = j-1;
                stop_ind_block = j;
                narrow_waves(start_ind_block+n_narrower:stop_ind_block-n_narrower, i) = ones(stop_ind_block-start_ind_block+1 - 2*n_narrower, 1);
            end
            j = j+1;
        end
    end
end
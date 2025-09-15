% This function widens blocks in a block wave, per column
function wide_waves = widen_blocks(waves, n_wider)
    n_pts = size(waves, 1);
    n_col = size(waves, 2);
    wide_waves = zeros(n_pts, n_col);

    for i = 1:n_col
        j = 1;
        while (j <= n_pts)
            % Catch where the block starts
            if waves(j, i)
                % Look how long the block is
                start_ind_block = j;
                while waves(j, i) && (j < n_pts)
                    j = j+1;
                end

                % Treat last point separately
                if not((j == n_pts) && waves(j, i))
                    j = j-1;
                end
                stop_ind_block = j;
                
                % Actually widen the block
                if ((start_ind_block-n_wider) >= 1) && ((stop_ind_block+n_wider) <= n_pts)
                    wide_waves(start_ind_block-n_wider:stop_ind_block+n_wider, i) = ones(stop_ind_block-start_ind_block+1 + 2*n_wider, 1);
                elseif ((start_ind_block-n_wider) >= 1)
                    n_ones = size(wide_waves(start_ind_block-n_wider:end, i), 1);
                    wide_waves(start_ind_block-n_wider:end, i) = ones(n_ones, 1);
                elseif ((stop_ind_block+n_wider) <= n_pts)
                    n_ones = size(wide_waves(1:stop_ind_block+n_wider, i), 1);
                    wide_waves(1:stop_ind_block+n_wider, i) = ones(n_ones, 1);
                else
                    wide_waves(:, i) = ones(n_pts, 1);
                end
            end
            j = j+1;
        end
    end
end